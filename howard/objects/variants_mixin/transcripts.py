import logging as log
import random
import string
import os
import vcf  # type: ignore

from howard.functions.commons import (transcripts_file_to_df, clean_annotation_field, get_file_format, code_type_map, full_path)

class variants_transcripts:

    ###############
    # Transcripts #
    ###############

    # Transcripts view creation

    def create_transcript_view(
        self,
        transcripts_table: str = None,
        transcripts_table_drop: bool = False,
        param: dict = {},
    ) -> str:
        """
        Generates a transcript view by processing data from a specified table based on provided parameters and structural information.

        Args:
            transcripts_table (str, optional): The name of the table that will store the final transcript view data.
                If not provided, the function will create a new table to store the transcript view data. Defaults to "transcripts".
            transcripts_table_drop (bool, optional): Determines whether to drop the existing transcripts table before creating a new one.
                If set to True, the function will drop the existing transcripts table if it exists. Defaults to False.
            param (dict, optional): A dictionary that contains information needed to create a transcript view.
                It includes details such as the structure of the transcripts, columns mapping, column formats, and other necessary information
                for generating the view. This parameter allows for flexibility and customization.

        Returns:
            str: The name of the transcripts table that was created or modified during the execution of the function.
        """

        log.info("Transcripts view creation")

        # Default
        transcripts_table_default = "transcripts"

        # Param
        if not param:
            param = self.get_param()

        # Struct
        struct = param.get("transcripts", {}).get("struct", None)

        # Transcript veresion
        transcript_id_remove_version = param.get("transcripts", {}).get(
            "transcript_id_remove_version", False
        )

        # Transcripts mapping
        transcript_id_mapping_file = param.get("transcripts", {}).get(
            "transcript_id_mapping_file", None
        )

        # Transcripts mapping
        transcript_id_mapping_force = param.get("transcripts", {}).get(
            "transcript_id_mapping_force", None
        )

        # Transcripts table
        if transcripts_table is None:
            transcripts_table = param.get("transcripts", {}).get(
                "table", transcripts_table_default
            )

        # Check transcripts table exists
        if transcripts_table:

            # Query to check if transcripts table exists
            query_check_table = f"""
                SELECT * 
                FROM information_schema.tables 
                WHERE table_name = '{transcripts_table}'
            """
            df_check_table = self.get_query_to_df(query=query_check_table)

            # Check if transcripts table exists
            if len(df_check_table) > 0 and not transcripts_table_drop:
                log.debug(f"Table {transcripts_table} exists and not drop option")
                log.info("Transcripts view creation - already exists")
                return transcripts_table

        # Variants table
        variants_table = self.get_table_variants()

        if struct:

            # added_columns
            added_columns = []

            # Temporary tables
            temporary_tables = []
            temporary_intermediate_tables = []

            # Annotation fields
            annotation_fields = []

            # Annotation fields
            annotation_fields_type = {}

            # from columns map
            # temporary_tables and annotation_fields are appended within the function
            log.info("Transcripts view creation - Annotations mapping...")
            columns_maps = struct.get("from_columns_map", [])
            (
                added_columns_tmp,
                temporary_tables_tmp,
                temporary_intermediate_tables_tmp,
                annotation_fields_tmp,
                annotation_fields_type_tmp,
            ) = self.create_transcript_view_from_columns_map(
                transcripts_table=transcripts_table,
                columns_maps=columns_maps,
                added_columns=added_columns,
                temporary_tables=temporary_tables,
                annotation_fields=annotation_fields,
            )

            # Append temporary tables infos
            added_columns += added_columns_tmp
            temporary_intermediate_tables += temporary_intermediate_tables_tmp
            for field in annotation_fields_type_tmp:
                field_type = annotation_fields_type_tmp.get(field, "VARCHAR")
                annotation_fields_type[field] = field_type

            # from column format
            # temporary_tables and annotation_fields are appended within the function
            log.info("Transcripts view creation - Annotations in format field...")
            column_formats = struct.get("from_column_format", [])
            (
                added_columns,
                temporary_tables_tmp,
                annotation_fields_tmp,
                added_columns_type_list,
            ) = self.create_transcript_view_from_column_format(
                transcripts_table=transcripts_table,
                column_formats=column_formats,
                temporary_tables=temporary_tables,
                view_type="table",
                annotation_fields=annotation_fields,
            )

            # Append temporary tables infos
            added_columns += added_columns_tmp
            for field in added_columns_type_list:
                annotation_fields_type[field] = added_columns_type_list.get(
                    field, "VARCHAR"
                )

            # Remove some specific fields/column
            annotation_fields = list(set(annotation_fields))
            for field in ["#CHROM", "POS", "REF", "ALT", "INFO", "transcript"]:
                if field in annotation_fields:
                    annotation_fields.remove(field)

            # Merge temporary tables query
            query_merge = ""
            for temporary_table in list(set(temporary_tables)):

                # First temporary table
                if not query_merge:
                    query_merge = f"""
                        SELECT * FROM {temporary_table}
                    """
                # other temporary table (using UNION)
                else:
                    query_merge += f"""
                        UNION BY NAME SELECT * FROM {temporary_table}
                    """

            # Create final merge query with transcript handling
            # Field transcript as None or '' to 'UNKNOWN' prevent issues with group by and association of variants with all avaialble transcripts
            # A field 'transcript_1' will be created within this table
            query_merge = f"""
                SELECT CASE WHEN "transcript" IS NULL THEN 'UNKNOWN' ELSE "transcript" END AS "transcript", *
                FROM ({query_merge}) AS transcripts_merged
            """

            # transcript table tmp
            transcript_table_tmp = "transcripts_tmp"
            transcript_table_tmp2 = "transcripts_tmp2"
            transcript_table_tmp3 = "transcripts_tmp3"

            # Merge on transcript
            query_merge_on_transcripts_annotation_fields = []

            # Add transcript list
            query_merge_on_transcripts_annotation_fields.append(
                f""" list_aggregate(list_distinct(array_agg({transcript_table_tmp}.transcript)), 'string_agg', ',') AS transcript_list """
            )

            # Aggregate all annotations fields
            for annotation_field in set(annotation_fields):

                # Annotation field type
                annotation_field_type = "VARCHAR"

                # Aggregate field
                query_merge_on_transcripts_annotation_fields.append(
                    f""" list_aggregate(list_distinct(array_agg({transcript_table_tmp}.{annotation_field})), 'string_agg', ',')::{annotation_field_type} AS {annotation_field} """
                )

            # Transcripts mapping
            if transcript_id_mapping_file:

                # Transcript dataframe
                transcript_id_mapping_dataframe_name = "transcript_id_mapping_dataframe"
                transcript_id_mapping_dataframe = transcripts_file_to_df(
                    transcript_id_mapping_file, column_names=["transcript", "alias"]
                )

                # Transcript version remove
                if transcript_id_remove_version:
                    query_transcript_column_select = f"split_part({transcript_table_tmp}.transcript, '.', 1) AS transcript_original, split_part({transcript_id_mapping_dataframe_name}.transcript, '.', 1) AS transcript_mapped"
                    query_transcript_column_group_by = f"split_part({transcript_table_tmp}.transcript, '.', 1), split_part({transcript_id_mapping_dataframe_name}.transcript, '.', 1)"
                    query_left_join = f"""
                        LEFT JOIN {transcript_id_mapping_dataframe_name} ON (split_part({transcript_id_mapping_dataframe_name}.alias, '.', 1)=split_part({transcript_table_tmp}.transcript, '.', 1))
                    """
                else:
                    query_transcript_column_select = f"{transcript_table_tmp}.transcript AS transcript_original, {transcript_id_mapping_dataframe_name}.transcript AS transcript_mapped"
                    query_transcript_column_group_by = f"{transcript_table_tmp}.transcript, {transcript_id_mapping_dataframe_name}.transcript"
                    query_left_join = f"""
                        LEFT JOIN {transcript_id_mapping_dataframe_name} ON (split_part({transcript_id_mapping_dataframe_name}.alias, '.', 1)=split_part({transcript_table_tmp}.transcript, '.', 1))
                    """

                # Transcript column for group by merge
                query_transcript_merge_group_by = """
                        CASE
                            WHEN transcript_mapped NOT IN ('')
                            THEN split_part(transcript_mapped, '.', 1)
                            ELSE split_part(transcript_original, '.', 1)
                        END
                    """

                # Merge query
                transcripts_tmp2_query = f"""
                    SELECT "#CHROM", POS, REF, ALT, {query_transcript_column_select}, {", ".join(query_merge_on_transcripts_annotation_fields)}
                    FROM ({query_merge}) AS {transcript_table_tmp}
                    {query_left_join}
                    GROUP BY "#CHROM", POS, REF, ALT, {query_transcript_column_group_by}
                """

                # Retrive columns after mege
                transcripts_tmp2_describe_query = f"""
                    DESCRIBE {transcripts_tmp2_query}
                """
                transcripts_tmp2_describe_list = list(
                    self.get_query_to_df(query=transcripts_tmp2_describe_query)[
                        "column_name"
                    ]
                )

                # Create list of columns for select clause
                transcripts_tmp2_describe_select_clause = []
                for field in transcripts_tmp2_describe_list:
                    if field not in [
                        "#CHROM",
                        "POS",
                        "REF",
                        "ALT",
                        "INFO",
                        "transcript_mapped",
                    ]:
                        as_field = field
                        if field in ["transcript_original"]:
                            as_field = "transcripts_mapped"
                        transcripts_tmp2_describe_select_clause.append(
                            f""" list_aggregate(list_distinct(array_agg({transcript_table_tmp2}.{field})), 'string_agg', ',') AS {as_field} """
                        )

                # Merge with mapping
                query_merge_on_transcripts = f"""
                    SELECT "#CHROM", POS, REF, ALT, '' AS INFO,
                        CASE
                            WHEN ANY_VALUE(transcript_mapped) NOT IN ('')
                            THEN ANY_VALUE(transcript_mapped)
                            ELSE ANY_VALUE(transcript_original)
                        END AS transcript,
                        {", ".join(transcripts_tmp2_describe_select_clause)}
                    FROM ({transcripts_tmp2_query}) AS {transcript_table_tmp2}
                    GROUP BY "#CHROM", POS, REF, ALT,
                        {query_transcript_merge_group_by}
                """

                # Add transcript filter from mapping file
                if transcript_id_mapping_force:
                    query_merge_on_transcripts = f"""
                        SELECT *
                        FROM ({query_merge_on_transcripts}) AS {transcript_table_tmp3}
                        WHERE split_part({transcript_table_tmp3}.transcript, '.', 1) in (SELECT split_part(transcript, '.', 1) FROM transcript_id_mapping_dataframe)
                    """

            # No transcript mapping
            else:

                # Remove transcript version
                if transcript_id_remove_version:
                    query_transcript_column = f"""
                        split_part({transcript_table_tmp}.transcript, '.', 1)
                    """
                else:
                    query_transcript_column = """
                        transcript
                    """

                # Query sections
                query_transcript_column_select = (
                    f"{query_transcript_column} AS transcript"
                )
                query_transcript_column_group_by = query_transcript_column

                # Query for transcripts view
                query_merge_on_transcripts = f"""
                    SELECT "#CHROM", POS, REF, ALT, '' AS INFO, {query_transcript_column} AS transcript, NULL AS transcript_mapped, {", ".join(query_merge_on_transcripts_annotation_fields)}
                    FROM ({query_merge}) AS {transcript_table_tmp}
                    GROUP BY "#CHROM", POS, REF, ALT, {query_transcript_column}
                """

            # Drop transcript view is necessary
            if transcripts_table_drop:
                query_drop = f"""
                    DROP TABLE IF EXISTS {transcripts_table};
                """
                self.execute_query(query=query_drop)

            # Log
            log.info(f"Transcripts view creation - Create view...")

            # Create table with structure but without data, if not exists
            query_create_table = f"""
                CREATE TABLE IF NOT EXISTS {transcripts_table} AS
                SELECT * FROM ({query_merge_on_transcripts}) LIMIT 0
            """
            self.execute_query(query=query_create_table)

            # Evaluate block size
            batch_split = self.get_batch_split()

            # Insert by batch
            for batch_index in range(batch_split):
                # where clause
                if batch_split > 1:
                    where_clause = f" WHERE (POS % {batch_split}) = {batch_index} "
                else:
                    where_clause = ""
                # Insert data
                query_insert_chunk = f"""
                    INSERT INTO {transcripts_table}
                    SELECT * FROM ({query_merge_on_transcripts})
                    {where_clause}
                """
                # Log
                log.debug(
                    f"Transcripts view creation - Insert batch [{batch_index+1}/{batch_split}]..."
                )
                # Execute
                self.execute_query(query=query_insert_chunk)

            # Extract annotations from variants

            # Columns from variants parameters
            columns_from_variants = struct.get("from_variants", {})
            columns_from_variants_prefix = columns_from_variants.get("prefix", "")
            columns_from_variants_fields = columns_from_variants.get("fields", [])
            columns_from_variants_info = columns_from_variants.get("INFO", False)

            # Columns from variants processing
            if len(columns_from_variants):
                log.info(
                    "Transcripts view creation - Extract annotations from variants"
                )

                # Add INFO column from variants table
                if columns_from_variants_info:
                    query_update_info_column = f"""
                        UPDATE {transcripts_table}
                        SET "INFO" = {variants_table}."INFO"
                        FROM {variants_table}
                        WHERE {transcripts_table}."#CHROM" = {variants_table}."#CHROM"
                        AND {transcripts_table}."POS" = {variants_table}."POS"
                        AND {transcripts_table}."REF" = {variants_table}."REF"
                        AND {transcripts_table}."ALT" = {variants_table}."ALT"
                    """
                    # log.debug(f"query_update_info_column={query_update_info_column}")
                    log.info(
                        "Transcripts view creation - Extract annotations from variants - All INFO column..."
                    )
                    self.execute_query(query=query_update_info_column)

                # Add columns from variants table as exploded from a list of fields
                if len(columns_from_variants_fields) > 0:
                    log.info(
                        f"Transcripts view creation - Extract annotations from variants - Extract {len(columns_from_variants_fields)} fields..."
                    )
                    fields_exploded = self.explode_infos(
                        fields=columns_from_variants_fields,
                        prefix=columns_from_variants_prefix,
                        table_source=variants_table,
                        table_dest=transcripts_table,
                        table_key=["#CHROM", "POS", "REF", "ALT"],
                        proccess_all_fields_together=True,
                        fields_not_exists=False,
                        fields_forced_as_varchar=False,
                    )
                    log.debug(
                        f"Transcripts view creation - Extract annotations from variants - Extract {len(columns_from_variants_fields)} fields: {fields_exploded}"
                    )

            # Remove temporary tables
            self.remove_tables_or_views(
                tables=temporary_tables + temporary_intermediate_tables
            )

            # Remove added columns
            for added_column in added_columns:
                self.drop_column(column=added_column)

        else:

            transcripts_table = None

        return transcripts_table

    def create_transcript_view_from_columns_map(
        self,
        transcripts_table: str = "transcripts",
        columns_maps: dict = {},
        added_columns: list = [],
        temporary_tables: list = None,
        annotation_fields: list = None,
        column_rename: dict = {},
        column_clean: bool = False,
        column_case: str = None,
    ) -> tuple[list, list, list]:
        """
        Generates a temporary table view based on specified columns mapping for transcripts data.

        Args:
            transcripts_table (str, optional): The name of the table where the transcripts data is stored or will be stored in the database.
                This table typically contains information about transcripts such as Ensembl transcript IDs, gene names, scores, predictions, etc.
                Defaults to "transcripts".
            columns_maps (dict): A dictionary that contains information about how to map columns from a transcripts table to create a view.
                Each entry in the dictionary represents a mapping configuration for a specific set of columns.
            added_columns (list): A list that stores the additional columns that will be added to the view being created based on the columns map provided.
                These columns are generated by exploding the transcript information columns along with the main transcript column.
            temporary_tables (list, optional): A list that stores the names of temporary tables created during the process of creating a transcript view from a columns map.
                These temporary tables are used to store intermediate results or transformations before the final view is generated.
            annotation_fields (list, optional): A list that stores the fields that are used for annotation in the query view creation process.
                These fields are extracted from the `transcripts_column` and `transcripts_infos_columns` specified in the `columns_maps`.
            column_rename (dict, optional): A dictionary that allows you to specify custom renaming for columns during the creation of the temporary table view.
                This parameter provides a mapping of original column names to the desired renamed column names.
            column_clean (bool, optional): A boolean flag that determines whether the column values should be cleaned or not.
                If set to `True`, the column values will be cleaned by removing any non-alphanumeric characters from them. Defaults to False.
            column_case (str, optional): Specifies the case transformation to be applied to the columns during the view creation process.
                It allows you to control whether the column values should be converted to lowercase, uppercase, or remain unchanged.

        Returns:
            tuple[list, list, list]: The function returns a tuple containing three lists: `added_columns`, `temporary_tables`, and `annotation_fields`.
        """

        log.debug("Start transcrpts view creation from columns map...")

        # "from_columns_map": [
        #     {
        #         "transcripts_column": "Ensembl_transcriptid",
        #         "transcripts_infos_columns": [
        #             "genename",
        #             "Ensembl_geneid",
        #             "LIST_S2_score",
        #             "LIST_S2_pred",
        #         ],
        #     },
        #     {
        #         "transcripts_column": "Ensembl_transcriptid",
        #         "transcripts_infos_columns": [
        #             "genename",
        #             "VARITY_R_score",
        #             "Aloft_pred",
        #         ],
        #     },
        # ],

        # Init
        if temporary_tables is None:
            temporary_tables = []
        if annotation_fields is None:
            annotation_fields = []

        # Init
        annotation_fields_type = {}
        temporary_intermediate_tables = []

        # Variants table
        table_variants = self.get_table_variants()

        for columns_map in columns_maps:

            # Log
            log.debug(f"columns_map={columns_map}")

            # Transcript column
            transcripts_column = columns_map.get("transcripts_column", None)

            # Transcripts infos columns
            transcripts_infos_columns = columns_map.get("transcripts_infos_columns", [])

            # Transcripts infos columns rename
            column_rename = columns_map.get("column_rename", column_rename)

            # Transcripts infos columns clean
            column_clean = columns_map.get("column_clean", column_clean)

            # Transcripts infos columns case
            column_case = columns_map.get("column_case", column_case)

            if transcripts_column is not None:

                table_for_view = table_variants

                annotation_view_name_for_type = None

                if True:

                    # Create annotations view
                    annotation_view_name = (
                        table_variants
                        + "_view_"
                        + "".join(
                            random.choices(string.ascii_uppercase + string.digits, k=10)
                        )
                    )
                    annotation_view_fields = [
                        transcripts_column
                    ] + transcripts_infos_columns
                    annotation_view_name = self.create_annotations_view(
                        table=table_variants,
                        view=annotation_view_name,
                        view_type="table",
                        view_mode="full",
                        info_prefix_column="",
                        detect_type_list=False,
                        fields=annotation_view_fields,
                        fields_not_exists=True,
                        fields_forced_as_varchar=True,
                        fields_needed_all=False,
                    )
                    temporary_intermediate_tables.append(annotation_view_name)
                    table_for_view = annotation_view_name

                    # Create annotation view for field type
                    annotation_view_name_for_type = self.create_annotations_view(
                        table=table_variants,
                        view=annotation_view_name + "for_type",
                        view_type="view",
                        view_mode="full",
                        info_prefix_column="",
                        detect_type_list=True,
                        fields=annotation_view_fields,
                        fields_not_exists=True,
                        fields_needed_all=False,
                    )
                    temporary_intermediate_tables.append(annotation_view_name_for_type)

                # View clauses
                clause_select_variants = []
                clause_select_tanscripts = []
                for field in [transcripts_column] + transcripts_infos_columns:

                    # AS field
                    as_field = field

                    # Rename
                    if column_rename:
                        as_field = column_rename.get(as_field, as_field)

                    # Clean
                    if column_clean:
                        as_field = clean_annotation_field(as_field)

                    # Case
                    if column_case:
                        if column_case.lower() in ["lower"]:
                            as_field = as_field.lower()
                        elif column_case.lower() in ["upper"]:
                            as_field = as_field.upper()

                    # Field Type
                    if annotation_view_name_for_type:
                        field_type = self.get_query_to_df(
                            f""" 
                                    SELECT column_type
                                    FROM (
                                        DESCRIBE {annotation_view_name_for_type}
                                    )
                                    WHERE column_name == '{field}'
                                """
                        )["column_type"][0].replace("[]", "")

                        # If field type is "NULL" due to no data
                        if field_type == '"NULL"':
                            field_type = "VARCHAR"
                    else:
                        field_type = "VARCHAR"

                    # Clause select Variants
                    clause_select_variants.append(
                        f""" TRY_CAST(regexp_split_to_table(CAST("{field}" AS VARCHAR), ',') AS VARCHAR) AS '{field}' """
                    )

                    # Clause select Transcripts
                    if field in [transcripts_column]:
                        clause_select_tanscripts.append(
                            f""" TRY_CAST(regexp_split_to_table("{field}", ',') AS {field_type}) AS '{field}' """
                        )
                    else:
                        clause_select_tanscripts.append(
                            f""" TRY_CAST(regexp_split_to_table("{field}", ',') AS {field_type}) AS '{as_field}' """
                        )
                        annotation_fields.append(as_field)
                        annotation_fields_type[as_field] = field_type

                # Query View
                query = f""" 
                    SELECT
                        "#CHROM", POS, REF, ALT,
                        CASE -- prevent issues with group by and association of variants with all available transcripts when no transcript information is available (NULL or '.' or '')
                            WHEN {transcripts_column} IN ('.', '') OR {transcripts_column} IS NULL
                            THEN 'UN_' || CAST(CAST(random() * 1000000 AS INTEGER) AS VARCHAR)
                            ELSE {transcripts_column}
                        END AS 'transcript',
                        {", ".join(clause_select_tanscripts)}
                    FROM (
                        SELECT 
                            "#CHROM", POS, REF, ALT,
                            {", ".join(clause_select_variants)}
                        FROM {table_for_view}
                        )
                    WHERE "{transcripts_column}" IS NOT NULL
                """

                # Create temporary table
                temporary_table = transcripts_table + "".join(
                    random.choices(string.ascii_uppercase + string.digits, k=10)
                )

                # Temporary view
                temporary_tables.append(temporary_table)
                query_view = f"""
                    CREATE view {temporary_table}
                    AS ({query})
                """
                # log.debug(f"Create view:{query_view}")
                self.execute_query(query=query_view)

        return (
            added_columns,
            temporary_tables,
            temporary_intermediate_tables,
            annotation_fields,
            annotation_fields_type,
        )

    def create_transcript_view_from_column_format(
        self,
        transcripts_table: str = "transcripts",
        column_formats: dict = {},
        temporary_tables: list = None,
        annotation_fields: list = None,
        column_rename: dict = {},
        view_type: str = "view",
        column_clean: bool = False,
        column_case: str = None,
    ) -> tuple[list, list, list]:
        """
        Generates a transcript view based on specified column formats, adds additional columns and annotation fields,
        and returns the list of temporary tables and annotation fields.

        Args:
            transcripts_table (str, optional): The name of the table containing the transcripts data. This table will be used
                as the base table for creating the transcript view. Defaults to "transcripts".
            column_formats (dict): A dictionary that contains information about the columns to be used for creating the transcript view.
                Each entry in the dictionary specifies the mapping between a transcripts column and a transcripts infos column.
            temporary_tables (list, optional): A list that stores the names of temporary views created during the process of creating
                a transcript view from a column format. These temporary views are used to manipulate and extract data before generating
                the final transcript view.
            annotation_fields (list, optional): A list that stores the annotation fields that are extracted from the temporary views
                created during the process. These annotation fields are obtained by querying the temporary views and extracting the column
                names excluding specific columns.
            column_rename (dict, optional): A dictionary that allows you to specify custom renaming of columns in the transcripts infos table.
                By providing a mapping of original column names to new column names in this dictionary, you can rename specific columns during
                the process.
            view_type (str, optional): The type of the view to be created. Defaults to "view".
            column_clean (bool, optional): A flag that determines whether the transcripts infos columns should undergo a cleaning process.
                If set to True, the columns will be cleaned during the creation of the transcript view based on the specified column format.
                Defaults to False.
            column_case (str, optional): Specifies the case transformation to be applied to the columns in the transcript view.
                It can be set to either "upper" or "lower" to convert the column names to uppercase or lowercase, respectively.

        Returns:
            tuple[list, list, list]: The function returns two lists: `temporary_tables` and `annotation_fields`.
        """

        log.debug("Start transcrpts view creation from column format...")

        #  "from_column_format": [
        #     {
        #         "transcripts_column": "ANN",
        #         "transcripts_infos_column": "Feature_ID",
        #     }
        # ],

        # Init
        if temporary_tables is None:
            temporary_tables = []
        if annotation_fields is None:
            annotation_fields = []

        added_columns = []
        added_columns_type_list = {}

        for column_format in column_formats:

            # annotation field and transcript annotation field
            annotation_field = column_format.get("transcripts_column", "ANN")
            transcript_annotation = column_format.get(
                "transcripts_infos_column", "Feature_ID"
            )

            # Transcripts infos columns rename
            column_rename = column_format.get("column_rename", column_rename)

            # Transcripts infos columns clean
            column_clean = column_format.get("column_clean", column_clean)

            # Transcripts infos columns case
            column_case = column_format.get("column_case", column_case)

            # Temporary View name
            temporary_view_name = transcripts_table + "".join(
                random.choices(string.ascii_uppercase + string.digits, k=10)
            )

            # Create temporary view name
            temporary_view_name, added_columns, added_columns_type = (
                self.annotation_format_to_table(
                    annotation_field=annotation_field,
                    view_name=temporary_view_name,
                    view_type=view_type,
                    annotation_id=transcript_annotation,
                    column_rename=column_rename,
                    column_clean=column_clean,
                    column_case=column_case,
                )
            )

            # columns_types
            for column_type in added_columns_type:
                added_columns_type_list[column_type] = added_columns_type.get(
                    column_type, "VARCHAR"
                )

            # Annotation fields
            if temporary_view_name:
                query_annotation_fields = f"""
                    SELECT *
                    FROM (
                        DESCRIBE SELECT *
                        FROM {temporary_view_name}
                        )
                        WHERE column_name not in ('#CHROM', 'POS', 'REF', 'ALT')
                """
                df_annotation_fields = self.get_query_to_df(
                    query=query_annotation_fields
                )

                # Add temporary view and annotation fields
                temporary_tables.append(temporary_view_name)
                annotation_fields += list(set(df_annotation_fields["column_name"]))

        return (
            added_columns,
            temporary_tables,
            annotation_fields,
            added_columns_type_list,
        )


    # Transcripts operations
    #######################

    def transcripts_export(
        self,
        transcripts_table: str = None,
        param_export: dict = {},
        param_explode: dict = {},
    ) -> bool:
        """
        Exports transcript data from a table to a specified file, with options for formatting and additional information.

        :param transcripts_table: The name of the transcripts table to export data from. If None, it defaults to "transcripts".
        :type transcripts_table: str, optional
        :param param_export: A dictionary of parameters to customize the export process, such as output file path, header options, etc.
        :type param_export: dict, optional
        :param param_explode: A dictionary of parameters for exploding fields in the transcripts table, such as prefix and fields to explode.
        :type param_explode: dict, optional
        :return: Returns True if the export is successful, False otherwise.
        :rtype: bool
        """

        log.debug("Start transcripts export...")

        # Param
        param = self.get_param()

        # Transcripts table
        if transcripts_table is None:
            transcripts_table = param.get("transcripts", {}).get("table", "transcripts")

        # Param export
        if not param_export:
            param_export = self.get_param().get("transcripts", {}).get("export", {})
        transcripts_export_output = param_export.get("output", None)
        transcripts_export_header = param_export.get("export_header", False)
        transcripts_export_header_in_output = param_export.get(
            "header_in_output", False
        )
        transcripts_export_add_info = param_export.get("add_info", False)

        if not param_export or not transcripts_export_output:
            log.warning(f"No transcriipts export parameters defined!")
            return False

        # Param explode
        if not param_explode:
            param_explode = self.get_param().get("transcripts", {}).get("explode", {})

        # Explode fields
        if param_explode.get("explode_infos_fields", None) and param_explode.get(
            "explode_infos", True
        ):
            self.explode_infos(
                table=transcripts_table,
                prefix=param_explode.get("explode_infos_prefix", None),
                fields=param_explode.get("explode_infos_fields", None),
                force=False,
                fields_forced_as_varchar=False,
                table_key=["#CHROM", "POS", "REF", "ALT", "transcript"],
            )

        # Create transcripts table description
        query_describe = f"""
            SELECT *
            FROM (
                    DESCRIBE SELECT * FROM {transcripts_table}
                )
            WHERE column_name NOT IN ('#CHROM', 'POS', 'REF', 'ALT', 'INFO')
        """
        result_describe = self.execute_query(query=query_describe)
        description_dict = {
            row[0]: {"type": row[1]} for row in result_describe.fetchall()
        }
        transcripts_annotations_list = list(description_dict.keys())

        transcripts_annotations_list_columns = [
            f'"{field}"' for field in transcripts_annotations_list
        ]

        # Output file format
        transcripts_export_output_format = get_file_format(
            filename=transcripts_export_output
        )

        # Format VCF - construct INFO
        if transcripts_export_output_format in ["vcf"]:

            # Construct query update INFO and header
            query_update_info = []
            for field in transcripts_annotations_list:

                # If field not in header
                if field not in self.get_header_infos_list():

                    # Find previous desc
                    if self.get_header().infos.get(field, None) is not None:
                        field_description = self.get_header().infos.get(field).desc
                        field_number = self.get_header().infos.get(field).num
                        field_type = self.get_header().infos.get(field).type
                    else:
                        field_description = "Unknown annotation"
                        field_number = "."
                        field_type = "String"

                    # Add description about transription prioritization
                    field_description += f". Annotation '{field}' from transcript view"

                    # Add PZ Transcript in header
                    self.get_header().infos[field] = vcf.parser._Info(
                        field,
                        field_number,
                        field_type,
                        field_description,
                        "unknown",
                        "unknown",
                        code_type_map.get(field_type, 0),
                    )

                # Add field as INFO/tag
                column_type = description_dict.get(field, {}).get("type", "VARCHAR")
                if column_type.endswith("[]"):
                    column_type = "VARCHAR"
                    field_value = f""" list_aggregate("{field}", 'string_agg', ',') """
                else:
                    field_value = f""" "{field}" """

                # Add INFO field to query
                query_update_info.append(
                    f"""
                        CASE
                            WHEN "{field}" IS NOT NULL
                            THEN concat(
                                '{field}=',
                                {field_value},
                                ';'
                            )    
                            ELSE ''     
                        END
                        """
                )

            # Query param
            query_update_info_value = f""" regexp_replace(concat('',  {", ".join(query_update_info)}), ';$', '') """
            query_export_columns = f""" "#CHROM", "POS", '.' AS 'ID', "REF", "ALT", '.' AS 'QUAL', '.' AS 'FILTER', "INFO" """

        else:

            # Query param

            if transcripts_export_add_info:
                query_update_info_value = f""" INFO """
                query_export_columns = f""" "#CHROM", "POS", "REF", "ALT", "INFO", {', '.join(transcripts_annotations_list)} """
            else:
                query_update_info_value = f""" NULL """
                query_export_columns = f""" "#CHROM", "POS", "REF", "ALT", {', '.join(transcripts_annotations_list)} """

        # Query export
        query_export = f"""
                SELECT
                {query_export_columns}
                FROM (
                    SELECT "#CHROM", "POS", "REF", "ALT",
                    {query_update_info_value} 
                    AS 'INFO',
                    {', '.join(transcripts_annotations_list_columns)}
                    FROM {transcripts_table}
                    ORDER BY "#CHROM", "POS", "REF", "ALT"
                )
            """

        # Export
        self.export_output(
            output_file=transcripts_export_output,
            query=query_export,
            export_header=transcripts_export_header,
            header_in_output=transcripts_export_header_in_output,
        )

    def transcripts_prioritization(
        self, transcripts_table: str = None, param: dict = {}, strict: bool = False
    ) -> bool:
        """
        Prioritizes transcripts based on specified parameters and updates the variants table with the prioritized information.

        Args:
            transcripts_table (str, optional): The name of the table containing transcripts data. If not provided, it defaults to "transcripts".
                This parameter is used to identify the table where the transcripts data is stored for the prioritization process.
            param (dict, optional): A dictionary containing various configuration settings for the prioritization process of transcripts.
                It is used to customize the behavior of the prioritization algorithm and includes settings such as the prefix for prioritization fields,
                default profiles, and other relevant configurations.
            strict (bool, optional): A flag indicating whether to enforce strict prioritization criteria. Defaults to False.

        Returns:
            bool: True if the transcripts prioritization process is successfully completed, and False if there are any issues or if no profile is defined
                for transcripts prioritization.
        """

        log.debug("Start transcripts prioritization...")

        # Param
        if not param:
            param = self.get_param()

        # Variants table
        table_variants = self.get_table_variants()

        # Transcripts table
        if transcripts_table is None:
            transcripts_table = self.create_transcript_view(
                transcripts_table="transcripts", param=param
            )
        if transcripts_table is None:
            msg_err = "No Transcripts table availalble"
            log.error(msg_err)
            raise ValueError(msg_err)

        # Get transcripts columns
        columns_as_list_query = f"""
            DESCRIBE {transcripts_table}
        """
        columns_as_list = list(
            self.get_query_to_df(columns_as_list_query)["column_name"]
        )

        # Create INFO if not exists
        if "INFO" not in columns_as_list:
            query_add_info = f"""
                ALTER TABLE {transcripts_table} ADD COLUMN INFO STRING DEFAULT '';
            """
            self.execute_query(query_add_info)

        # Prioritization param and Force only PZ Score and Flag
        pz_param = param.get("transcripts", {}).get("prioritization", {})

        # PZ profile by default
        pz_profile_default = (
            param.get("transcripts", {}).get("prioritization", {}).get("profiles", None)
        )

        # Exit if no profile
        if pz_profile_default is None:
            log.warning("No profile defined for transcripts prioritization")
            return False

        # PZ fields
        pz_param_pzfields = {}

        # Order by
        pz_orders = (
            param.get("transcripts", {})
            .get("prioritization", {})
            .get("prioritization_transcripts_order", {})
        )
        if not pz_orders:
            pz_orders = {
                pz_param.get("pzprefix", "PTZ") + "Flag": "DESC",
                pz_param.get("pzprefix", "PTZ") + "Score": "DESC",
            }

        # PZ field transcripts
        pz_fields_transcripts = pz_param.get("pzprefix", "PTZ") + "Transcript"

        # Add description about transription prioritization
        pz_field_transcripts_description = f"Transcript selected from prioritization process, profile {pz_profile_default}"

        # Add PZ Transcript in header
        self.get_header().infos[pz_fields_transcripts] = vcf.parser._Info(
            pz_fields_transcripts,
            1,
            "String",
            pz_field_transcripts_description,
            "HOWARD transcript prioritization",
            "unknown",
            code_type_map.get("String", 0),
        )

        # Mandatory fields if asked in param
        pz_mandatory_fields_list = [
            "Score",
            "Flag",
            "Tags",
            "Comment",
            "Infos",
            "Class",
        ]
        pz_mandatory_fields = []
        for pz_mandatory_field in pz_mandatory_fields_list:
            pz_mandatory_fields.append(
                pz_param.get("pzprefix", "PTZ") + pz_mandatory_field
            )

        # PZ fields in param
        pz_param_mandatory_fields = []
        for pz_field in pz_param.get("pzfields", []):
            if pz_field in pz_mandatory_fields_list:
                pz_param_pzfields[pz_param.get("pzprefix", "PTZ") + pz_field] = (
                    pz_param.get("pzprefix", "PTZ") + pz_field
                )
                pz_param_mandatory_fields.append(
                    pz_param.get("pzprefix", "PTZ") + pz_field
                )
            else:
                pz_field_new = pz_param.get("pzprefix", "PTZ") + pz_field
                pz_param_pzfields[pz_field] = pz_field_new

                # Find previous desc and type
                if self.get_header().infos.get(pz_field, None) is not None:
                    pz_field_new_description = (
                        self.get_header().infos.get(pz_field).desc
                    )
                    pz_field_new_type = self.get_header().infos.get(pz_field).type
                else:
                    pz_field_new_description = "Unknown annotation"
                    pz_field_new_type = "String"

                # Add description about transription prioritization
                pz_field_new_description += f". Annotation '{pz_field}' from transcript selected from prioritization process, profile {pz_profile_default}"

                # Add PZ Transcript in header
                self.get_header().infos[pz_field_new] = vcf.parser._Info(
                    pz_field_new,
                    1,
                    pz_field_new_type,
                    pz_field_new_description,
                    "unknown",
                    "unknown",
                    code_type_map.get(pz_field_new_type, 0),
                )
        # Add order by fields in mandatory fields
        for pz_order in pz_orders:
            if pz_order not in pz_param_mandatory_fields:
                pz_param_mandatory_fields.append(pz_order)

        # PZ fields param
        pz_mandatory_fields = pz_param_mandatory_fields
        pz_param["pzfields"] = pz_mandatory_fields

        # Prioritization
        prioritization_result = self.prioritization(
            table=transcripts_table,
            pz_param=param.get("transcripts", {}).get("prioritization", {}),
            pz_keys=["#CHROM", "POS", "REF", "ALT", "transcript"],
            strict=strict,
        )
        if not prioritization_result:
            log.warning("Transcripts prioritization not processed")
            return False

        log.info(f"Update {table_variants} table with transcripts prioritization...")

        # PZ fields sql query
        query_update_select_list = []
        query_update_concat_list = []
        query_update_order_list = []
        for pz_param_pzfield in set(
            list(pz_param_pzfields.keys()) + pz_mandatory_fields
        ):
            query_update_select_list.append(f" {pz_param_pzfield}, ")

        for pz_param_pzfield in pz_param_pzfields:
            query_update_concat_list.append(
                f"""
                    , CASE 
                        WHEN {pz_param_pzfield} IS NOT NULL
                        THEN concat(';{pz_param_pzfields.get(pz_param_pzfield)}=', {pz_param_pzfield})
                        ELSE ''
                    END
                """
            )

        for pz_order in pz_orders:
            query_update_order_list.append(
                f""" {pz_order} {pz_orders.get(pz_order, "DESC")} """
            )

        # Fields to explode
        fields_to_explode = (
            list(pz_param_pzfields.keys())
            + pz_mandatory_fields
            + list(pz_orders.keys())
        )

        # Remove transcript column as a specific transcript column
        if "transcript" in fields_to_explode:
            fields_to_explode.remove("transcript")

        # Fields in transcripts table
        query_transcripts_table = f"""
            DESCRIBE SELECT * FROM {transcripts_table}
        """
        query_transcripts_table_description = self.get_query_to_df(
            query=query_transcripts_table
        )

        # Check fields to explode
        for field_to_explode in fields_to_explode:
            if field_to_explode not in self.get_header_infos_list() + list(
                query_transcripts_table_description.column_name
            ):
                msg_err = f"INFO/{field_to_explode} NOT IN header"
                log.error(msg_err)
                raise ValueError(msg_err)

        # Transcript preference file
        # First look for transcripts_preference, then for transcripts_preference_file for backward compatibility, and get full path
        transcripts_preference_file = (
            param.get("transcripts", {})
            .get("prioritization", {})
            .get("prioritization_transcripts", {})
        )
        transcripts_preference_file = (
            param.get("transcripts", {})
            .get("prioritization", {})
            .get("prioritization_transcripts_file", transcripts_preference_file)
        )
        transcripts_preference_file = full_path(transcripts_preference_file)

        # Transcript preference columns
        transcript_preference_columns = (
            param.get("transcripts", {})
            .get("prioritization", {})
            .get("prioritization_transcripts_columns", [])
        )
        if transcript_preference_columns is None:
            transcript_preference_columns = []

        # Transcript preference forced
        transcript_preference_force = (
            param.get("transcripts", {})
            .get("prioritization", {})
            .get("prioritization_transcripts_force", False)
        )
        # Transcript version forced
        transcript_version_force = (
            param.get("transcripts", {})
            .get("prioritization", {})
            .get("prioritization_transcripts_version_force", False)
        )

        # Create view as table
        annotation_view_name = "annotation_view_for_transcripts_prioritization_" + str(
            random.randrange(1000000)
        )
        annotation_view_name = self.create_annotations_view(
            table=transcripts_table,
            view=annotation_view_name,
            view_type="table",
            view_mode="explore",
            info_prefix_column="",
            detect_type_list=True,
            fields=fields_to_explode + ["transcript"] + transcript_preference_columns,
            fields_not_exists=False,
            fields_forced_as_varchar=False,
            fields_needed_all=False,
        )
        transcripts_table = annotation_view_name

        # Transcripts Ranking

        # Default values

        # Order by
        order_by = " , ".join(query_update_order_list)

        # Left join
        left_join = ""

        # where clause
        where_clause = ""

        if transcripts_preference_file:

            # Transcripts file to dataframe
            if os.path.exists(transcripts_preference_file):
                transcripts_preference_dataframe = transcripts_file_to_df(
                    transcripts_preference_file
                )
            else:
                log.error(
                    f"Transcript file '{transcripts_preference_file}' does NOT exist"
                )
                raise ValueError(
                    f"Transcript file '{transcripts_preference_file}' does NOT exist"
                )

            # Order by depending to transcript preference forcing
            if transcript_preference_force:
                order_by = f""" transcripts_preference.transcripts_preference_order ASC, {order_by} """
            else:
                order_by = f""" {order_by}, transcripts_preference.transcripts_preference_order ASC """

            # Transcript columns joined depend on version consideration
            if transcript_version_force:
                transcripts_version_join = f""" {transcripts_table}.transcript = transcripts_preference.transcripts_preference """
            else:
                transcripts_version_join = f""" split_part({transcripts_table}.transcript, '.', 1) = split_part(transcripts_preference.transcripts_preference, '.', 1) """

            # Left join
            left_join = f"""
                LEFT JOIN 
                    (
                        SELECT transcript AS 'transcripts_preference', row_number() OVER () AS transcripts_preference_order
                        FROM transcripts_preference_dataframe
                    ) AS transcripts_preference
                ON {transcripts_version_join}
            """

        if transcript_preference_columns:
            # Add transcript preference columns in select to filter only if exists in this column (can be varchar (list of values separated by ',') or array

            # List of transcript filter on columns depending on forcing or not of transcript preference and version
            transcripts_filter_array = []

            for transcript_preference_column in transcript_preference_columns:

                # Check if transcript preference column exists in transcripts table description and get its type
                try:
                    transcript_preference_column_type = (
                        query_transcripts_table_description[
                            query_transcripts_table_description.get("column_name", None)
                            == transcript_preference_column
                        ]
                        .get("column_type", None)
                        .iloc[0]
                    )
                except:
                    msg_err = f"Transcript preference column '{transcript_preference_column}' not found in transcripts table description"
                    log.error(msg_err)
                    raise ValueError(msg_err)

                # Transcript filter as list if is an array or not
                if transcript_preference_column_type.endswith("[]"):
                    transcript_filter_column_as_list = f""" {transcripts_table}."{transcript_preference_column}"::VARCHAR[] """
                else:
                    transcript_filter_column_as_list = f""" string_split({transcripts_table}."{transcript_preference_column}"::VARCHAR, ',')"""

                # Create transcript filter for current column depending on forcing or not of transcript preference and version
                if transcript_version_force:
                    transcripts_filter_array.append(f"""
                        array_contains({transcript_filter_column_as_list}, {transcripts_table}."transcript")
                    """)
                else:
                    transcripts_filter_array.append(f"""
                         array_contains(list_transform({transcript_filter_column_as_list}, x -> split_part(x, '.', 1)), split_part({transcripts_table}.transcript, '.', 1))
                    """)

            # Create transcript filter on columns
            if where_clause:
                where_clause += " AND "
            else:
                where_clause = " WHERE"
            where_clause += (
                "(" + " OR ".join(transcripts_filter_array) + ")"
                if transcripts_filter_array
                else ""
            )

        # Query ranking for update
        query_update_ranking = f"""
            SELECT
                "#CHROM", POS, REF, ALT,
                {transcripts_table}.transcript AS transcript,
                {" ".join(query_update_select_list)}
                ROW_NUMBER() OVER (
                    PARTITION BY "#CHROM", POS, REF, ALT
                    ORDER BY {order_by}
                ) AS rn
            FROM {transcripts_table}
            {left_join}
            {where_clause}
        """

        # Export Transcripts prioritization infos to variants table
        query_update = f"""
            WITH RankedTranscripts AS (
                {query_update_ranking}
            )
            UPDATE {table_variants}
                SET
                INFO = CONCAT(CASE
                            WHEN "INFO" IS NULL OR "INFO" IN ('','.')
                            THEN ''
                            ELSE concat("INFO", ';')
                        END,
                        concat('{pz_fields_transcripts}=', transcript {" ".join(query_update_concat_list)})
                        )
            FROM
                RankedTranscripts
            WHERE
                rn = 1
                AND variants."#CHROM" = RankedTranscripts."#CHROM"
                AND variants."POS" = RankedTranscripts."POS"
                AND variants."REF" = RankedTranscripts."REF"
                AND variants."ALT" = RankedTranscripts."ALT"     
        """

        # Query update
        self.execute_query(query=query_update)

        # Return
        return True

    def transcript_view_to_variants(
        self,
        transcripts_table: str = None,
        transcripts_column_id: str = None,
        transcripts_info_json: str = None,
        transcripts_info_field_json: str = None,
        transcripts_info_format: str = None,
        transcripts_info_field_format: str = None,
        param: dict = {},
    ) -> bool:
        """
        The `transcript_view_to_variants` function updates a variants table with information from
        transcripts in JSON format.

        :param transcripts_table: The `transcripts_table` parameter is used to specify the name of the
        table containing the transcripts data. If this parameter is not provided, the function will
        attempt to retrieve it from the `param` dictionary or use a default value of "transcripts"
        :type transcripts_table: str
        :param transcripts_column_id: The `transcripts_column_id` parameter is used to specify the
        column in the `transcripts_table` that contains the unique identifier for each transcript. This
        identifier is used to match transcripts with variants in the database
        :type transcripts_column_id: str
        :param transcripts_info_json: The `transcripts_info_json` parameter is used to specify the name
        of the column in the variants table where the transcripts information will be stored in JSON
        format. This parameter allows you to define the column in the variants table that will hold the
        JSON-formatted information about transcripts
        :type transcripts_info_json: str
        :param transcripts_info_field_json: The `transcripts_info_field_json` parameter is used to
        specify the field in the VCF header that will contain information about transcripts in JSON
        format. This field will be added to the VCF header as an INFO field with the specified name
        :type transcripts_info_field_json: str
        :param transcripts_info_format: The `transcripts_info_format` parameter is used to specify the
        format of the information about transcripts that will be stored in the variants table. This
        format can be used to define how the transcript information will be structured or displayed
        within the variants table
        :type transcripts_info_format: str
        :param transcripts_info_field_format: The `transcripts_info_field_format` parameter is used to
        specify the field in the VCF header that will contain information about transcripts in a
        specific format. This field will be added to the VCF header as an INFO field with the specified
        name
        :type transcripts_info_field_format: str
        :param param: The `param` parameter in the `transcript_view_to_variants` method is a dictionary
        that contains various configuration settings related to transcripts. It is used to provide
        default values for certain parameters if they are not explicitly provided when calling the
        method. The `param` dictionary can be passed as an argument
        :type param: dict
        :return: The function `transcript_view_to_variants` returns a boolean value. It returns `True`
        if the operation is successful and `False` if certain conditions are not met.
        """

        msg_info_prefix = "Start transcripts view to variants annotations"

        log.debug(f"{msg_info_prefix}...")

        # Default
        transcripts_table_default = "transcripts"
        transcripts_column_id_default = "transcript"
        transcripts_info_json_default = None
        transcripts_info_format_default = None
        transcripts_info_field_json_default = None
        transcripts_info_field_format_default = None

        # Param
        if not param:
            param = self.get_param()

        # Transcripts table
        if transcripts_table is None:
            transcripts_table = param.get("transcripts", {}).get(
                "table", transcripts_table_default
            )

        # Transcripts column ID
        if transcripts_column_id is None:
            transcripts_column_id = param.get("transcripts", {}).get(
                "column_id", transcripts_column_id_default
            )

        # Transcripts info json
        if transcripts_info_json is None:
            transcripts_info_json = param.get("transcripts", {}).get(
                "transcripts_info_json", transcripts_info_json_default
            )

        # Transcripts info field JSON
        if transcripts_info_field_json is None:
            transcripts_info_field_json = param.get("transcripts", {}).get(
                "transcripts_info_field_json", transcripts_info_field_json_default
            )
        # if transcripts_info_field_json is not None and transcripts_info_json is None:
        #     transcripts_info_json = transcripts_info_field_json

        # Transcripts info format
        if transcripts_info_format is None:
            transcripts_info_format = param.get("transcripts", {}).get(
                "transcripts_info_format", transcripts_info_format_default
            )

        # Transcripts info field FORMAT
        if transcripts_info_field_format is None:
            transcripts_info_field_format = param.get("transcripts", {}).get(
                "transcripts_info_field_format", transcripts_info_field_format_default
            )
        # if (
        #     transcripts_info_field_format is not None
        #     and transcripts_info_format is None
        # ):
        #     transcripts_info_format = transcripts_info_field_format

        # Variants table
        table_variants = self.get_table_variants()

        # Check info columns param
        if (
            transcripts_info_json is None
            and transcripts_info_field_json is None
            and transcripts_info_format is None
            and transcripts_info_field_format is None
        ):
            return False

        # Transcripts infos columns
        query_transcripts_infos_columns = f"""
            SELECT *
            FROM (
                DESCRIBE SELECT * FROM {transcripts_table}
                )
            WHERE "column_name" NOT IN ('#CHROM', 'POS', 'REF', 'ALT', '{transcripts_column_id}')
        """
        transcripts_infos_columns = list(
            self.get_query_to_df(query=query_transcripts_infos_columns)["column_name"]
        )

        # View results
        clause_select = []
        clause_to_json = []
        clause_to_format = []
        for field in transcripts_infos_columns:
            # Do not consider INFO field for export into fields
            if field not in ["INFO"]:
                clause_select.append(
                    f""" regexp_split_to_table(CAST("{field}" AS STRING), ',') AS '{field}' """
                )
                clause_to_json.append(f""" '{field}': "{field}" """)
                clause_to_format.append(f""" "{field}" """)

        # Update
        update_set_json = []
        update_set_format = []

        # VCF header
        vcf_reader = self.get_header()

        # Transcripts to info column in JSON
        if transcripts_info_json:

            # Create column on variants table
            self.add_column(
                table_name=table_variants,
                column_name=transcripts_info_json,
                column_type="JSON",
                default_value=None,
                drop=False,
            )

            # Add header
            vcf_reader.infos[transcripts_info_json] = vcf.parser._Info(
                transcripts_info_json,
                ".",
                "String",
                "Transcripts in JSON format",
                "unknwon",
                "unknwon",
                self.code_type_map["String"],
            )

            # Add to update
            update_set_json.append(
                f""" {transcripts_info_json}=t.{transcripts_info_json} """
            )

        # Transcripts to info field in JSON
        if transcripts_info_field_json:

            log.debug(f"{msg_info_prefix} - Annotation in JSON format...")

            # Add to update
            update_set_json.append(
                f""" 
                    INFO = concat(
                            CASE
                                WHEN INFO NOT IN ('', '.')
                                THEN INFO
                                ELSE ''
                            END,
                            CASE
                                WHEN CAST(t.{transcripts_info_json} AS VARCHAR) NOT IN ('', '.')
                                THEN concat(
                                    ';{transcripts_info_field_json}=',
                                    t.{transcripts_info_json}
                                )
                                ELSE ''
                            END
                            )
                """
            )

            # Add header
            vcf_reader.infos[transcripts_info_field_json] = vcf.parser._Info(
                transcripts_info_field_json,
                ".",
                "String",
                "Transcripts in JSON format",
                "unknwon",
                "unknwon",
                self.code_type_map["String"],
            )

        if update_set_json:

            # Update query
            query_update = f"""
                UPDATE {table_variants}
                    SET {", ".join(update_set_json)}
                FROM
                (
                    SELECT
                        "#CHROM", POS, REF, ALT,
                            concat(
                            '{{',
                            string_agg(
                                '"' || "{transcripts_column_id}" || '":' ||
                                to_json(json_output)
                            ),
                            '}}'
                            )::JSON AS {transcripts_info_json}
                    FROM
                        (
                        SELECT
                            "#CHROM", POS, REF, ALT,
                            "{transcripts_column_id}",
                            to_json(
                                {{{",".join(clause_to_json)}}}
                            )::JSON AS json_output
                        FROM
                            (SELECT "#CHROM", POS, REF, ALT, "{transcripts_column_id}", {", ".join(clause_select)} FROM {transcripts_table})
                        WHERE "{transcripts_column_id}" IS NOT NULL
                        )
                    GROUP BY "#CHROM", POS, REF, ALT
                ) AS t
                WHERE {table_variants}."#CHROM" = t."#CHROM"
                    AND {table_variants}."POS" = t."POS"
                    AND {table_variants}."REF" = t."REF"
                    AND {table_variants}."ALT" = t."ALT"
            """

            self.execute_query(query=query_update)

        # Transcripts to info column in FORMAT
        if transcripts_info_format:

            # Create column on variants table
            self.add_column(
                table_name=table_variants,
                column_name=transcripts_info_format,
                column_type="VARCHAR",
                default_value=None,
                drop=False,
            )

            # Add header
            vcf_reader.infos[transcripts_info_format] = vcf.parser._Info(
                transcripts_info_format,
                ".",
                "String",
                f"Transcripts annotations: 'transcript | {' | '.join(transcripts_infos_columns)}'",
                "unknwon",
                "unknwon",
                self.code_type_map["String"],
            )

            # Add to update
            update_set_format.append(
                f""" {transcripts_info_format}=t.{transcripts_info_format} """
            )

        else:

            # Set variable for internal queries
            transcripts_info_format = "transcripts_info_format"

        # Transcripts to info field in JSON
        if transcripts_info_field_format:

            log.debug(f"{msg_info_prefix} - Annotation in structured format...")

            # Add to update
            update_set_format.append(
                f""" 
                    INFO = concat(
                            CASE
                                WHEN INFO NOT IN ('', '.')
                                THEN INFO
                                ELSE ''
                            END,
                            CASE
                                WHEN CAST(t.{transcripts_info_format} AS VARCHAR) NOT IN ('', '.')
                                THEN concat(
                                    ';{transcripts_info_field_format}=',
                                    t.{transcripts_info_format}
                                )
                                ELSE ''
                            END
                            )
                """
            )

            # Add header
            vcf_reader.infos[transcripts_info_field_format] = vcf.parser._Info(
                transcripts_info_field_format,
                ".",
                "String",
                f"Transcripts annotations: 'transcript | {' | '.join(transcripts_infos_columns)}'",
                "unknwon",
                "unknwon",
                self.code_type_map["String"],
            )

        if update_set_format:

            # Update query
            query_update = f"""
                UPDATE {table_variants}
                    SET {", ".join(update_set_format)}
                FROM
                (
                    SELECT
                        "#CHROM", POS, REF, ALT,
                            string_agg({transcripts_info_format}) AS {transcripts_info_format}
                    FROM 
                        (
                        SELECT
                            "#CHROM", POS, REF, ALT,
                            "{transcripts_column_id}",
                            concat(
                                "{transcripts_column_id}",
                                '|',
                                {", '|', ".join(clause_to_format)}
                            ) AS {transcripts_info_format}
                        FROM
                            (SELECT "#CHROM", POS, REF, ALT, "{transcripts_column_id}", {", ".join(clause_select)} FROM {transcripts_table})
                        )
                    GROUP BY "#CHROM", POS, REF, ALT
                ) AS t
                WHERE {table_variants}."#CHROM" = t."#CHROM"
                    AND {table_variants}."POS" = t."POS"
                    AND {table_variants}."REF" = t."REF"
                    AND {table_variants}."ALT" = t."ALT"
            """

            self.execute_query(query=query_update)

        return True
