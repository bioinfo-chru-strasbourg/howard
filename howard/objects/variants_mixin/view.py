import logging as log

from howard.functions.commons import (code_type_map_to_sql)

class variants_view:

    ####################
    # Annotations view #
    ####################

    def create_annotations_view(
        self,
        table: str = None,
        view: str = None,
        view_type: str = None,
        view_mode: str = None,
        fields: list = None,
        fields_needed: list = None,
        fields_needed_all: bool = False,
        detect_type_list: bool = True,
        fields_not_exists: bool = True,
        only_in_columns: bool = False,
        formats: list = None,
        strict: bool = False,
        info_prefix_column: str = None,
        info_struct_column: str = None,
        sample_struct_column: str = None,
        drop_view: bool = False,
        fields_to_rename: dict = None,
        fields_forced_as_varchar: bool = False,
        limit: int = None,
    ) -> str:
        """
        Creates a SQL view from fields in a VCF INFO column, or already in a column.

        :param table: The name of the table containing variant data. If not provided, the default table is used.
        :type table: str, optional
        :param view: The name of the view that will be created based on the fields in the VCF INFO column. Defaults to None.
        :type view: str, optional
        :param view_type: The type of view to be created. It can be either a `VIEW` or a `TABLE`. Defaults to `VIEW`.
        :type view_type: str, optional
        :param view_mode: The mode of view to be created. It can be either `full` or `explore`. Defaults to `full`.
        :type view_mode: str, optional
        :param fields: A list of field names to be extracted from the INFO column in the VCF file. Defaults to None.
        :type fields: list, optional
        :param fields_needed: A list of fields that are required for the view. Defaults to None.
        :type fields_needed: list, optional
        :param fields_needed_all: A flag that determines whether to include all fields in the table in the view. Defaults to False.
        :type fields_needed_all: bool, optional
        :param detect_type_list: A flag that determines whether to detect the type of the fields extracted from the INFO column. Defaults to True.
        :type detect_type_list: bool, optional
        :param fields_not_exists: A flag that determines whether to include fields that do not exist in the table in the view. Defaults to True.
        :type fields_not_exists: bool, optional
        :param only_in_columns: A flag that determines whether to include only the fields that exist in the columns of the table. Defaults to False.
        :type only_in_columns: bool, optional
        :param formats: A list of field names to be extracted from the FORMAT column in the VCF file. Defaults to None.
        :type formats: list, optional
        :param strict: A flag that determines whether to enforce strict criteria for the fields in the view. Defaults to False.
        :type strict: bool, optional
        :param info_prefix_column: A prefix to be added to the field names in the view. Defaults to None.
        :type info_prefix_column: str, optional
        :param info_struct_column: The name of the column that will contain the extracted fields from the INFO column in the view. Defaults to None.
        :type info_struct_column: str, optional
        :param sample_struct_column: The name of the column that will contain the extracted formats from the samples columns in the view. Defaults to None.
        :type sample_struct_column: str, optional
        :param drop_view: A flag that determines whether to drop the existing view with the same name before creating a new view. Defaults to False.
        :type drop_view: bool, optional
        :param fields_to_rename: A dictionary that contains the mapping of fields to be renamed in the VCF file. Defaults to None.
        :type fields_to_rename: dict, optional
        :param fields_forced_as_varchar: A flag that forces fields to be treated as type VARCHAR. Defaults to False.
        :type fields_forced_as_varchar: bool, optional
        :param limit: The maximum number of rows to be included in the view. Defaults to None.
        :type limit: int, optional

        :return: The name of the view that is created based on the fields extracted from the INFO column in the VCF file.
        :rtype: str
        """

        # Create a sql view from fields in VCF INFO column, with each column is a field present in the VCF header (with a specific type from VCF header) and extracted from INFO column (with a regexp like in rename_info_fields), and each row is a variant.

        # Get table
        if table is None:
            table = self.get_table_variants()

        # Get view
        if view is None:
            view = f"{table}_view"

        # Get view type
        if view_type is None:
            view_type = "VIEW"

        # Get mode
        view_mode_allowed = ["full", "explore"]
        if view_mode is None:
            view_mode = "explore"

        # Mode lower
        view_mode = view_mode.lower()

        # Mode check
        if view_mode not in view_mode_allowed:
            msg_err = f"Invalid view mode: '{view_mode}' (either {view_mode_allowed})"
            log.error(msg_err)
            raise ValueError(msg_err)

        # Prefix
        if info_prefix_column is not None:
            prefix = info_prefix_column
        else:
            prefix = ""

        # Check view type value
        view_type_allowed = ["view", "table"]
        if view_type.lower() not in view_type_allowed:
            msg_err = f"Invalid view type: {view_type} (either {view_type_allowed})"
            log.error(msg_err)
            raise ValueError(msg_err)

        # Get header
        header = self.get_header()

        # Get fields
        if fields is None:
            fields = list(header.infos.keys())

        # # Get format fields
        # if formats is None:
        #     formats = list(header.formats.keys())
        #     # fields = list(header.infos.keys())

        # Get fields to rename
        if fields_to_rename is None:
            fields_to_rename = {}

        # If Samples structured columns
        if sample_struct_column:

            # # Get format
            # formats = list(header.formats.keys())
            if formats is None:
                formats = list(header.formats.keys())

            # Get samples
            samples = list(header.samples)

        else:

            # Empty format and samples
            formats = []
            samples = []

        log.debug(
            f"Create '{view}' view (as '{view_type}' mode '{view_mode}') from table '{table}' with {len(fields)} fields"
        )

        connexion_type = self.get_connexion_type()

        # Describe table
        if connexion_type in ["duckdb"]:
            table_describe_query = f"""
                DESCRIBE {table}
            """
            table_describe = self.get_query_to_df(query=table_describe_query)
        else:
            table_describe_query = f"""
                PRAGMA table_info({table})
            """
            table_describe = self.get_query_to_df(query=table_describe_query)
            table_describe["column_name"] = table_describe.get("name")

        # fields needed
        if fields_needed is None:
            if fields_needed_all:
                fields_needed = list(table_describe.get("column_name"))
            else:
                fields_needed = ["#CHROM", "POS", "REF", "ALT"]

        # Add samples in view mode 'full'
        # log.debug(f"samples={samples}")
        # log.debug(f"table_describe={table_describe}")
        if view_mode in ["full"] and sample_struct_column and len(samples):
            if "FORMAT" not in fields_needed and "FORMAT" in list(
                table_describe.get("column_name")
            ):  # in table_describe:
                fields_needed += ["FORMAT"]
            for field_needed in fields_needed:
                if field_needed not in fields_needed:
                    fields_needed += [field_needed]

        # Check needed fieds
        for field in fields_needed:
            if field not in list(table_describe.get("column_name")):
                msg_err = f"Field '{field}' is needed, but not in file"
                raise ValueError(msg_err)

        # Create fields for annotation view extracted from INFO column in table variants (with regexp_replace like in rename_info_fields), with column type from VCF header
        fields_columns = []
        fields_columns_annotations_struct = []
        samples_format_struct = []
        field_sql_type_list = False

        # Find "INFO" column
        if "INFO" in list(table_describe.get("column_name")):
            info_column = '"INFO"'
        else:
            info_column = "''"

        # Each field
        for field in set(fields):

            # Rename field
            field_to_rename = fields_to_rename.get(field, field)
            if field_to_rename is None:
                field_to_rename = field

            # Check field type

            # Field info
            field_infos = header.infos.get(field, None)

            # Field SQL type
            if field_infos is not None:

                # Field SQL type
                field_sql_type = code_type_map_to_sql.get(field_infos.type, "VARCHAR")

                # Column is a list
                if detect_type_list and str(field_infos.num) != "1":
                    field_sql_type_list = True
                else:
                    field_sql_type_list = False

            else:

                # Field SQL type
                field_sql_type = "VARCHAR"

                # Column is a list
                field_sql_type_list = False

            # fields_forced_as_varchar
            if fields_forced_as_varchar:
                field_sql_type = "VARCHAR"
                field_sql_type_list = False

            # Needed fields, not in other annotation fields (useful for DB with fields in column)
            if field in fields_needed and not field in list(
                table_describe.get("column_name")
            ):
                continue

            # Fields in table
            elif field in list(table_describe.get("column_name")):

                # Add field in needed if 'full' view mode
                if view_mode in ["full"]:
                    if field not in fields_needed:
                        fields_needed += [field]

                # Only if not needes (already in a column)
                if not field in fields_needed:
                    # log.debug(f"Filed '{field}' not in needed")
                    fields_columns.append(f"""
                            "{field}" AS '{prefix}{field_to_rename}' -- field in column but not in needed
                        """)

                # Flag
                if field_infos is not None and field_infos.type == "Flag":
                    fields_columns_annotations_struct.append(f"""
                            "{field_to_rename}":= TRY_CAST("{field}" AS BOOLEAN)
                        """)
                else:
                    if field_sql_type_list:
                        fields_columns_annotations_struct.append(f"""
                                "{field_to_rename}":= CAST(list_transform(string_split(CAST("{field}" AS VARCHAR), ','), x -> CASE WHEN x = '.' OR x = '' THEN NULL ELSE x END) AS {field_sql_type}[]) -- field in column
                            """)
                    else:
                        fields_columns_annotations_struct.append(f"""
                                "{field_to_rename}":= COALESCE(NULLIF(regexp_replace(CAST("{field}" AS VARCHAR), '^\\.$', ''), '')::{field_sql_type}, NULL)  -- field in column
                            """)

            # Fields in header
            elif (
                field in header.infos
                and not only_in_columns
                and "INFO" in list(table_describe.get("column_name"))
            ):

                # Colonne is a flag
                if field_infos.type == "Flag":

                    # Field pattern
                    field_pattern = rf"(^|;)({field})([^;]*)?"

                    if view_mode in ["explore"]:
                        fields_columns.append(f"""
                                regexp_matches({info_column}, '{field_pattern}')::BOOLEAN AS '{prefix}{field_to_rename}'
                            """)
                        fields_columns_annotations_struct.append(f"""
                                "{field_to_rename}":= regexp_matches({info_column}, '{field_pattern}')::BOOLEAN
                            """)
                    elif view_mode in ["full"]:
                        fields_columns.append(f"""
                                string_agg(CASE WHEN k = '{field}' THEN true END, ',')::BOOLEAN AS '{prefix}{field_to_rename}'
                            """)
                        fields_columns_annotations_struct.append(f"""
                                "{field_to_rename}":= string_agg(CASE WHEN k = '{field}' THEN true END, ',')::BOOLEAN
                            """)

                # Colonne with a type
                else:

                    # Field pattern
                    field_pattern = rf"(^|;)({field})=([^;]*)?"

                    if view_mode in ["explore"]:
                        field_source = (
                            f""" regexp_extract({info_column}, '{field_pattern}', 3) """
                        )
                    elif view_mode in ["full"]:
                        field_source = (
                            f""" string_agg(CASE WHEN k = '{field}' THEN v END, ',') """
                        )

                    # Field is a list
                    if field_sql_type_list:

                        fields_columns.append(f"""
                                CAST(list_transform(string_split({field_source}, ','), x -> CASE WHEN x = '.' OR x = '' THEN NULL ELSE x END) AS {field_sql_type}[]) AS '{prefix}{field_to_rename}'
                            """)
                        fields_columns_annotations_struct.append(f"""
                                "{field_to_rename}":= CAST(list_transform(string_split({field_source}, ','), x -> CASE WHEN x = '.' OR x = '' THEN NULL ELSE x END) AS {field_sql_type}[])
                            """)

                    # Field is a unique value
                    else:

                        fields_columns.append(f"""
                                NULLIF(regexp_replace({field_source}, '^\\.$', ''), '')::{field_sql_type} AS '{prefix}{field_to_rename}'
                            """)
                        fields_columns_annotations_struct.append(f"""
                                "{field_to_rename}":= COALESCE(NULLIF(regexp_replace({field_source}, '^\\.$', ''), '')::{field_sql_type}, NULL)
                            """)

            # Add field even if not exists
            elif fields_not_exists:

                fields_columns.append(f"""
                            null AS '{prefix}{field_to_rename}'
                        """)
                fields_columns_annotations_struct.append(f"""
                            "{field_to_rename}":= NULL
                        """)
                msg_err = f"Field '{field}' is not found (in table or header): '{field}' will be set to NULL"
                log.warning(msg=msg_err)

            else:

                # Field not found
                msg_err = f"Field '{field}' is not found (in table or header or column)"

                if strict:
                    log.error(msg=msg_err)
                    raise ValueError(msg_err)
                else:
                    log.warning(msg=msg_err)

        # If samples and struct as option
        if sample_struct_column and len(samples):

            # Format info
            format_infos = header.formats

            # For each sample
            for sample in samples:

                # Struct by format
                sample_format_struct = []

                # For each format
                for format in formats:
                    # for format in ["GT"]:

                    # Format cast and list
                    format_cast = ""
                    format_list = False
                    format_cast = code_type_map_to_sql.get(
                        format_infos.get(format).type, "VARCHAR"
                    )
                    if format_infos.get(format).num != 1:
                        format_list = True

                    # If format is a list
                    if format_list:
                        sample_format_struct.append(f""" 
                                "{format}":= 
                                    list_transform(
                                        string_split(
                                            NULLIF(
                                                string_split(CAST("{sample}" AS VARCHAR), ':')[list_position(string_split("FORMAT", ':'), '{format}')]
                                                , ''
                                            )
                                            , ',')
                                        , x -> CASE WHEN x = '.' OR x = '' THEN NULL ELSE x END
                                    )::{format_cast}[]
                            """)
                    # If format is NOT a list
                    else:
                        sample_format_struct.append(f""" 
                                "{format}":= 
                                    COALESCE(
                                        NULLIF(
                                            regexp_replace(
                                                string_split(CAST("{sample}" AS VARCHAR), ':')[list_position(string_split("FORMAT", ':'), '{format}')]
                                                , '^\\.$', ''
                                            )
                                        , ''
                                        )
                                    )::{format_cast}
                            """)

                # Add struct of the sample
                if len(sample_format_struct):
                    samples_format_struct.append(f"""
                        "{sample}":= STRUCT_PACK({", ".join(sample_format_struct)})
                    """)

        # Combine fields into columns
        if info_prefix_column is not None and len(fields_columns):
            annotations_column_annotations_columns = (
                f""", {", ".join(fields_columns)}"""
            )
        else:
            annotations_column_annotations_columns = ""

        # Combine fields into a STRUCT
        if info_struct_column and len(fields_columns_annotations_struct):
            annotations_column_annotations_struct = f""" 
                , STRUCT_PACK({", ".join(fields_columns_annotations_struct)}) AS {info_struct_column}
                """
        else:
            annotations_column_annotations_struct = ""

        # Combine samples into a STRUCT
        if sample_struct_column and len(samples_format_struct):
            samples_format_struct_clause = f""", STRUCT_PACK({", ".join(samples_format_struct)}) AS {sample_struct_column} """
        else:
            samples_format_struct_clause = ""

        # Limit
        limit_clause = ""
        if limit is not None:
            limit_clause = f" LIMIT {limit} "

        # Query select

        if view_mode in ["explore"]:
            query_select = f"""
                SELECT
                    {', '.join([f'"{field}"' for field in fields_needed])}  -- variant id
                    {annotations_column_annotations_columns}                -- annotations_column_annotations_columns
                    {annotations_column_annotations_struct}                 -- annotations_column_annotations_struct
                    {samples_format_struct_clause}                          -- samples_format_struct_clause
                FROM
                    {table}
                {limit_clause}
            """

        elif view_mode in ["full"]:
            query_select = f"""
                    SELECT
                        {', '.join([f'"{field}"' for field in fields_needed])}          -- variant id
                        {annotations_column_annotations_columns}                        -- annotations_column_annotations_columns
                        {annotations_column_annotations_struct}                         -- annotations_column_annotations_struct
                        {samples_format_struct_clause}                                  -- samples_format_struct_clause
                    FROM (
                        SELECT
                            {', '.join([f'"{field}"' for field in fields_needed])},     -- variant id
                            k,      -- key
                            v       -- value
                        FROM (
                            SELECT
                                {', '.join([f'"{field}"' for field in fields_needed])},     -- variant id
                                INFO,                                                       -- INFO
                                string_split(kv, '=')[1] AS k,  -- key
                                string_split(kv, '=')[2] AS v   -- value
                            FROM (
                                SELECT {', '.join([f'"{field}"' for field in fields_needed])},  -- variant id
                                {info_column} AS INFO,                                          -- INFO
                                -- unnest(string_split({info_column}, ';')) AS kv
                                unnest(string_split(concat({info_column}, ''), ';')) AS kv
                                FROM {table}
                                )
                            WHERE k in ('{"', '".join(fields)}') OR TRIM(INFO) in ('', '.')  OR INFO IS NULL
                            )
                        )
                    GROUP BY {', '.join([f'"{field}"' for field in fields_needed])}
                    {limit_clause}

            """

        # Drop if any
        if drop_view:
            log.debug(f"Drop view: {view}")
            try:
                query_create_view = f"""
                    DROP view IF EXISTS {view}
                """
                self.execute_query(query=query_create_view)
                log.debug(f"View dropped: {view}")
            except:
                try:
                    query_create_view = f"""
                        DROP table IF EXISTS {view}
                    """
                    self.execute_query(query=query_create_view)
                    log.debug(f"View dropped: {view}")
                except:
                    msg_err = f"View '{view}' can NOT be dropped"
                    log.error(msg_err)
                    raise ValueError(msg_err)

        # Create view
        log.debug(f"Create view: {view}")
        query_create_view = f"""
            CREATE {view_type} IF NOT EXISTS {view} AS {query_select}
        """
        # log.debug(f"Create view:{query_create_view}")
        self.execute_query(query=query_create_view)
        log.debug(f"View created: {view}")

        return view

    def remove_tables_or_views(self, tables: list = None, views: list = None) -> list:
        """
        Remove specified tables and views from the database.

        Args:
            tables (list): A list of table names to be removed. Default is None.
            views (list): A list of view names to be removed. Default is None.

        Returns:
            list: A list of tables and views that were successfully removed.

        This function attempts to remove the specified tables and views from the database.
        It first tries to drop each item as a table, and if that fails, it tries to drop it as a view.
        If an item is neither a table nor a view, an error is logged.
        """

        temporary_tables = (tables or []) + (views or [])
        removed_items = []

        # Remove temporary tables and views
        if temporary_tables:
            for temporary_table in set(temporary_tables):
                try:
                    query_drop_tmp_table = f"""
                        DROP TABLE IF EXISTS {temporary_table}
                    """
                    self.execute_query(query=query_drop_tmp_table)
                    log.debug(f"DROP TABLE '{temporary_table}' done.")
                    removed_items.append(temporary_table)
                except Exception as e:
                    log.debug(
                        f"DROP TABLE '{temporary_table}': Failed (not a table)! Try as a view."
                    )

                    try:
                        query_drop_tmp_view = f"""
                            DROP VIEW IF EXISTS {temporary_table}
                        """
                        self.execute_query(query=query_drop_tmp_view)
                        log.debug(f"DROP VIEW '{temporary_table}' done.")
                        removed_items.append(temporary_table)
                    except Exception as e:
                        log.debug(f"DROP VIEW '{temporary_table}': Failed (not a view)")
                        log.error(
                            f"DROP '{temporary_table}': Failed! Neither a table nor a view"
                        )

        return removed_items
