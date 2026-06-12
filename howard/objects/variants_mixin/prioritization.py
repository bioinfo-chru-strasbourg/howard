import random
import re
import numpy as np  # type: ignore
import vcf  # type: ignore
import logging as log

from howard.functions.commons import (
    full_path,
    code_type_map,
    comparison_map
)


class variants_prioritization:

    ###############
    # Calculation #
    ###############

    def get_config_prioritizations_default(self) -> dict:
        """
        The function `get_config_prioritizations_default` returns a dictionary containing default configurations for
        various prioritizations.

        :return: The function `get_config_prioritizations_default` returns a dictionary containing default configuration
        settings for different calculations.
        """

        config_default = {
            "default": {
                "ANN": [
                    {
                        "type": "contains",
                        "value": "HIGH",
                        "score": 5,
                        "flag": "PASS",
                        "comment": [
                            "The variant is assumed to have high (disruptive) impact in the protein, probably causing protein truncation, loss of function or triggering nonsense mediated decay"
                        ],
                    },
                    {
                        "type": "contains",
                        "value": "MODERATE",
                        "score": 3,
                        "flag": "PASS",
                        "comment": [
                            "A non-disruptive variant that might change protein effectiveness"
                        ],
                    },
                    {
                        "type": "contains",
                        "value": "LOW",
                        "score": 0,
                        "flag": "FILTERED",
                        "comment": [
                            "Assumed to be mostly harmless or unlikely to change protein behavior"
                        ],
                    },
                    {
                        "type": "contains",
                        "value": "MODIFIER",
                        "score": 0,
                        "flag": "FILTERED",
                        "comment": [
                            "Usually non-coding variants or variants affecting non-coding genes, where predictions are difficult or there is no evidence of impact"
                        ],
                    },
                ],
            }
        }

        return config_default

    def prioritization(
        self,
        table: str = None,
        pz_prefix: str = None,
        pz_param: dict = None,
        pz_keys: list = None,
        strict: bool = False,
    ) -> bool:
        """
        Processes VCF files, adds new INFO fields, and prioritizes variants based on configured profiles and criteria.

        Args:
            table (str, optional): The name of the table (presumably a VCF file) on which the prioritization operation will be performed.
                If not provided, the default variants table will be used.
            pz_prefix (str, optional): A prefix to be added to certain INFO fields in the VCF file during the prioritization process.
                Defaults to "PZ" if not provided.
            pz_param (dict, optional): Additional parameters specific to the prioritization process. These parameters can include settings
                related to prioritization profiles, fields, scoring modes, flags, comments, and other configurations needed for the prioritization
                of variants.
            pz_keys (list, optional): The keys used to join the prioritization table with the variant table. Defaults to ["#CHROM", "POS", "REF", "ALT"]
                if not provided.
            strict (bool, optional): Whether to enforce strict prioritization criteria availability in view (need to be in header and in column). Defaults to False.

        Returns:
            bool: True if the prioritization operation is successful, False otherwise.
        """

        # Config
        config = self.get_config()

        # Param
        param = self.get_param()

        # Prioritization param
        if pz_param is not None:
            prioritization_param = pz_param
        else:
            prioritization_param = param.get("prioritization", {})

        # Configuration profiles
        prioritization_config_file = prioritization_param.get(
            "prioritization_config", None
        )
        prioritization_config_file = full_path(prioritization_config_file)
        prioritizations_config = self.get_config_json(
            name="prioritizations", config_file=prioritization_config_file
        )

        # Prioritization prefix
        pz_prefix_default = "PZ"
        if pz_prefix is None:
            pz_prefix = prioritization_param.get("pzprefix", pz_prefix_default)

        # Prioritization options
        profiles = prioritization_param.get("profiles", [])
        if isinstance(profiles, str):
            profiles = profiles.split(",")
        pzfields = prioritization_param.get(
            "pzfields", [f"{pz_prefix}Flag", f"{pz_prefix}Score"]
        )
        if isinstance(pzfields, str):
            pzfields = pzfields.split(",")
        default_profile = prioritization_param.get("default_profile", None)
        pzfields_sep = prioritization_param.get("pzfields_sep", "_")
        prioritization_score_mode = prioritization_param.get(
            "prioritization_score_mode", "HOWARD"
        )

        # Quick Prioritizations
        prioritizations = param.get("prioritizations", None)
        if prioritizations:
            log.info("Quick Prioritization:")
            for profile in prioritizations.split(","):
                if profile not in profiles:
                    profiles.append(profile)
                    log.info(f"   {profile}")

        # Keys for prioritization join
        if pz_keys is None:
            pz_keys = ["#CHROM", "POS", "REF", "ALT"]

        # If profile "ALL" provided, all profiles in the config profiles
        if "ALL" in profiles:
            profiles = list(prioritizations_config.keys())

        for profile in profiles:
            if prioritizations_config.get(profile, None):
                log.debug(f"Profile '{profile}' configured")
            else:
                msg_error = f"Profile '{profile}' NOT configured"
                log.error(msg_error)
                raise ValueError(msg_error)

        if profiles:
            log.info(f"Prioritization... ")
        else:
            log.debug(f"No profile defined")
            return False

        if not default_profile and len(profiles):
            default_profile = profiles[0]

        log.debug("Profiles availables: " + str(list(prioritizations_config.keys())))
        log.debug("Profiles to check: " + str(list(profiles)))

        # Variables
        if table is not None:
            table_variants = table
        else:
            table_variants = self.get_table_variants(clause="update")
        log.debug(f"Table to prioritize: {table_variants}")

        # Added columns
        added_columns = []

        # Create list of PZfields
        # List of PZFields
        list_of_pzfields_original = pzfields + [
            pzfield + pzfields_sep + profile
            for pzfield in pzfields
            for profile in profiles
        ]
        list_of_pzfields = []
        log.debug(f"{list_of_pzfields_original}")

        # Remove existing PZfields to use if exists
        for pzfield in list_of_pzfields_original:
            if self.get_header().infos.get(pzfield, None) is None:
                list_of_pzfields.append(pzfield)
                log.debug(f"VCF Input - Header - PZfield '{pzfield}' not in VCF")
            else:
                log.debug(f"VCF Input - Header - PZfield '{pzfield}' already in VCF")

        if list_of_pzfields:

            # PZfields tags description
            PZfields_INFOS = {
                f"{pz_prefix}Tags": {
                    "ID": f"{pz_prefix}Tags",
                    "Number": ".",
                    "Type": "String",
                    "Description": "Variant tags based on annotation criteria",
                },
                f"{pz_prefix}Score": {
                    "ID": f"{pz_prefix}Score",
                    "Number": 1,
                    "Type": "Integer",
                    "Description": "Variant score based on annotation criteria",
                },
                f"{pz_prefix}Flag": {
                    "ID": f"{pz_prefix}Flag",
                    "Number": 1,
                    "Type": "String",
                    "Description": "Variant flag based on annotation criteria",
                },
                f"{pz_prefix}Comment": {
                    "ID": f"{pz_prefix}Comment",
                    "Number": ".",
                    "Type": "String",
                    "Description": "Variant comment based on annotation criteria",
                },
                f"{pz_prefix}Infos": {
                    "ID": f"{pz_prefix}Infos",
                    "Number": ".",
                    "Type": "String",
                    "Description": "Variant infos based on annotation criteria",
                },
                f"{pz_prefix}Class": {
                    "ID": f"{pz_prefix}Class",
                    "Number": ".",
                    "Type": "String",
                    "Description": "Variant class based on annotation criteria",
                },
            }

            # Create INFO fields if not exist
            for field in PZfields_INFOS:
                field_ID = PZfields_INFOS[field]["ID"]
                field_description = PZfields_INFOS[field]["Description"]
                if field_ID not in self.get_header().infos and field_ID in pzfields:
                    field_description = (
                        PZfields_INFOS[field]["Description"]
                        + f", profile {default_profile}"
                    )
                    self.get_header().infos[field_ID] = vcf.parser._Info(
                        field_ID,
                        PZfields_INFOS[field]["Number"],
                        PZfields_INFOS[field]["Type"],
                        field_description,
                        "unknown",
                        "unknown",
                        code_type_map[PZfields_INFOS[field]["Type"]],
                    )

            # Create INFO fields if not exist for each profile
            for profile in prioritizations_config:
                if profile in profiles or profiles == []:
                    for field in PZfields_INFOS:
                        field_ID = PZfields_INFOS[field]["ID"] + pzfields_sep + profile
                        field_description = (
                            PZfields_INFOS[field]["Description"]
                            + f", profile {profile}"
                        )
                        if (
                            field_ID not in self.get_header().infos
                            and field in pzfields
                        ):
                            self.get_header().infos[field_ID] = vcf.parser._Info(
                                field_ID,
                                PZfields_INFOS[field]["Number"],
                                PZfields_INFOS[field]["Type"],
                                field_description,
                                "unknown",
                                "unknown",
                                code_type_map[PZfields_INFOS[field]["Type"]],
                            )

            # Header
            for pzfield in list_of_pzfields:
                if re.match(f"{pz_prefix}Score.*", pzfield):
                    added_column = self.add_column(
                        table_name=table_variants,
                        column_name=pzfield,
                        column_type="INTEGER",
                        default_value="0",
                    )
                elif re.match(f"{pz_prefix}Flag.*", pzfield):
                    added_column = self.add_column(
                        table_name=table_variants,
                        column_name=pzfield,
                        column_type="BOOLEAN",
                        default_value="1",
                    )
                elif re.match(f"{pz_prefix}Class.*", pzfield):
                    added_column = self.add_column(
                        table_name=table_variants,
                        column_name=pzfield,
                        column_type="VARCHAR[]",
                        default_value="null",
                    )
                else:
                    added_column = self.add_column(
                        table_name=table_variants,
                        column_name=pzfield,
                        column_type="STRING",
                        default_value="''",
                    )
                added_columns.append(added_column)

            # Profiles
            if profiles:

                # foreach profile in configuration file
                for profile in prioritizations_config:

                    # If profile is asked in param, or ALL are asked (empty profile [])
                    if profile in profiles or profiles == []:
                        log.info(f"Profile '{profile}'")

                        sql_set_info_option = ""

                        sql_set_info = []

                        # PZ fields set

                        # PZScore
                        if (
                            f"{pz_prefix}Score{pzfields_sep}{profile}"
                            in list_of_pzfields
                        ):
                            sql_set_info.append(f"""
                                    concat(
                                        '{pz_prefix}Score{pzfields_sep}{profile}=',
                                        {pz_prefix}Score{pzfields_sep}{profile}
                                    ) 
                                """)
                            if (
                                profile == default_profile
                                and f"{pz_prefix}Score" in list_of_pzfields
                            ):
                                sql_set_info.append(f"""
                                        concat(
                                            '{pz_prefix}Score=',
                                            {pz_prefix}Score{pzfields_sep}{profile}
                                        )
                                    """)

                        # PZFlag
                        if (
                            f"{pz_prefix}Flag{pzfields_sep}{profile}"
                            in list_of_pzfields
                        ):
                            sql_set_info.append(f"""
                                    concat(
                                        '{pz_prefix}Flag{pzfields_sep}{profile}=',
                                        CASE 
                                            WHEN {pz_prefix}Flag{pzfields_sep}{profile}==1
                                            THEN 'PASS'
                                            WHEN {pz_prefix}Flag{pzfields_sep}{profile}==0
                                            THEN 'FILTERED'
                                        END
                                    ) 
                                """)
                            if (
                                profile == default_profile
                                and f"{pz_prefix}Flag" in list_of_pzfields
                            ):
                                sql_set_info.append(f"""
                                        concat(
                                            '{pz_prefix}Flag=',
                                            CASE 
                                                WHEN {pz_prefix}Flag{pzfields_sep}{profile}==1
                                                THEN 'PASS'
                                                WHEN {pz_prefix}Flag{pzfields_sep}{profile}==0
                                                THEN 'FILTERED'
                                            END
                                        )
                                    """)

                        # PZClass
                        if (
                            f"{pz_prefix}Class{pzfields_sep}{profile}"
                            in list_of_pzfields
                        ):
                            sql_set_info.append(f"""
                                    concat(
                                        '{pz_prefix}Class{pzfields_sep}{profile}=',
                                        CASE
                                            WHEN len({pz_prefix}Class{pzfields_sep}{profile}) > 0
                                            THEN list_aggregate(list_distinct({pz_prefix}Class{pzfields_sep}{profile}), 'string_agg', ',')
                                            ELSE '.'
                                        END 
                                    )
                                    
                                """)
                            if (
                                profile == default_profile
                                and f"{pz_prefix}Class" in list_of_pzfields
                            ):
                                sql_set_info.append(f"""
                                        concat(
                                            '{pz_prefix}Class=',
                                            CASE
                                                WHEN len({pz_prefix}Class{pzfields_sep}{profile}) > 0
                                                THEN list_aggregate(list_distinct({pz_prefix}Class{pzfields_sep}{profile}), 'string_agg', ',')
                                                ELSE '.'
                                            END 
                                        )
                                    """)

                        # PZComment
                        if (
                            f"{pz_prefix}Comment{pzfields_sep}{profile}"
                            in list_of_pzfields
                        ):
                            sql_set_info.append(f"""
                                    CASE
                                        WHEN {pz_prefix}Comment{pzfields_sep}{profile} NOT IN ('')
                                        THEN concat('{pz_prefix}Comment{pzfields_sep}{profile}=', {pz_prefix}Comment{pzfields_sep}{profile})
                                        ELSE ''
                                    END
                                """)
                            if (
                                profile == default_profile
                                and f"{pz_prefix}Comment" in list_of_pzfields
                            ):
                                sql_set_info.append(f"""
                                        CASE
                                            WHEN {pz_prefix}Comment{pzfields_sep}{profile} NOT IN ('')
                                            THEN concat('{pz_prefix}Comment=', {pz_prefix}Comment{pzfields_sep}{profile})
                                            ELSE ''
                                        END
                                    """)

                        # PZInfos
                        if (
                            f"{pz_prefix}Infos{pzfields_sep}{profile}"
                            in list_of_pzfields
                        ):
                            sql_set_info.append(f"""
                                    CASE
                                        WHEN {pz_prefix}Infos{pzfields_sep}{profile} NOT IN ('')
                                        THEN concat('{pz_prefix}Infos{pzfields_sep}{profile}=', {pz_prefix}Infos{pzfields_sep}{profile})
                                        ELSE ''
                                    END
                                """)
                            if (
                                profile == default_profile
                                and f"{pz_prefix}Infos" in list_of_pzfields
                            ):
                                sql_set_info.append(f"""
                                        CASE
                                            WHEN {pz_prefix}Infos{pzfields_sep}{profile} NOT IN ('')
                                            THEN concat('{pz_prefix}Infos=', {pz_prefix}Infos{pzfields_sep}{profile})
                                            ELSE ''
                                        END
                                    """)

                        # Merge PZfields
                        sql_set_info_option = ""
                        sql_set_sep = ""
                        for sql_set in sql_set_info:
                            if sql_set_sep:
                                sql_set_info_option += f"""
                                    , concat('{sql_set_sep}', {sql_set})
                                """
                            else:
                                sql_set_info_option += f"""
                                    , {sql_set}
                                """
                            sql_set_sep = ";"

                        sql_queries = []
                        criterion_fields_profile = []
                        annotation_view_name = (
                            "annotation_view_for_prioritization_"
                            + str(random.randrange(1000000))
                        )
                        annotations_view_prefix = ""
                        annotations_view_struct = "INFOS"
                        for annotation in prioritizations_config[profile]:

                            # skip special sections
                            if annotation.startswith("_"):
                                continue

                            # Log
                            log.info(f"Profile '{profile}' - Filter '{annotation}'")

                            # For each criterions
                            for criterion in prioritizations_config[profile][
                                annotation
                            ]:

                                # Criterion mode
                                criterion_mode = None
                                if np.any(
                                    np.isin(list(criterion.keys()), ["type", "value"])
                                ):
                                    criterion_mode = "operation"
                                elif np.any(
                                    np.isin(list(criterion.keys()), ["sql", "fields"])
                                ):
                                    criterion_mode = "sql"
                                log.debug(f"Criterion Mode: {criterion_mode}")

                                if criterion_mode in ["operation"]:
                                    log.warning(
                                        f"Prioritization criterion mode '{criterion_mode}' is deprecated. Please use 'sql' mode instead."
                                    )
                                    log.debug(f"Criterion: {criterion}")

                                # Criterion parameters
                                criterion_type = criterion.get("type", None)
                                criterion_value = criterion.get("value", None)
                                criterion_sql = criterion.get("sql", None)
                                criterion_fields = criterion.get("fields", None)
                                criterion_score = criterion.get("score", 0)
                                criterion_flag = criterion.get("flag", "PASS")
                                criterion_class = criterion.get("class", None)
                                criterion_flag_bool = criterion_flag == "PASS"
                                criterion_comment = (
                                    ", ".join(criterion.get("comment", []))
                                    .replace("'", "''")
                                    .replace(";", ",")
                                    .replace("\t", " ")
                                )
                                criterion_infos = (
                                    str(criterion)
                                    .replace("'", "''")
                                    .replace(";", ",")
                                    .replace("\t", " ")
                                )

                                # SQL
                                if criterion_sql is not None and isinstance(
                                    criterion_sql, list
                                ):
                                    criterion_sql = " ".join(criterion_sql)

                                # Fields and explode
                                if criterion_fields is None:
                                    criterion_fields = [annotation]
                                if not isinstance(criterion_fields, list):
                                    criterion_fields = str(criterion_fields).split(",")

                                # Class
                                if criterion_class is not None and not isinstance(
                                    criterion_class, list
                                ):
                                    criterion_class = str(criterion_class).split(",")

                                # Add criterion fields to the list of profile's criteria
                                criterion_fields_profile = list(
                                    set(criterion_fields_profile + criterion_fields)
                                )

                                # Create annotations view for prioritization
                                log.debug(
                                    f"""Profile '{profile}' - Prioritization - Create '{annotation_view_name}' view with '{criterion_fields_profile}'... """
                                )
                                annotation_view_name = self.create_annotations_view(
                                    view=annotation_view_name,
                                    table=table_variants,
                                    view_type="view",
                                    view_mode="explore",
                                    info_prefix_column=annotations_view_prefix,
                                    info_struct_column=annotations_view_struct,
                                    fields=criterion_fields_profile + pz_keys,
                                    fields_not_exists=(not strict),
                                    only_in_columns=strict,
                                    strict=strict,
                                    drop_view=True,
                                    detect_type_list=True,
                                )

                                # Describe annotation view and dict
                                annotation_view_describe = self.get_query_to_df(
                                    f"DESCRIBE {annotation_view_name}"
                                )
                                annotation_view_describe_dict = (
                                    annotation_view_describe.set_index("column_name")[
                                        "column_type"
                                    ].to_dict()
                                )

                                # Keys for join
                                clause_join = []
                                for key in pz_keys:
                                    if key in annotation_view_describe_dict:
                                        clause_join.append(
                                            f""" "{table_variants}"."{key}" == "{annotation_view_name}"."{key}" """
                                        )

                                sql_set = []
                                sql_set_info = []

                                # PZ fields set

                                # PZScore
                                if (
                                    f"{pz_prefix}Score{pzfields_sep}{profile}"
                                    in list_of_pzfields
                                ):
                                    # VaRank prioritization score mode
                                    if prioritization_score_mode.upper().strip() in [
                                        "VARANK",
                                        "MAX",
                                        "MAXIMUM",
                                        "TOP",
                                    ]:
                                        sql_set.append(
                                            f"{pz_prefix}Score{pzfields_sep}{profile} = CASE WHEN {criterion_score}>{pz_prefix}Score{pzfields_sep}{profile} THEN {criterion_score} ELSE {pz_prefix}Score{pzfields_sep}{profile} END "
                                        )
                                    # default HOWARD prioritization score mode
                                    else:
                                        sql_set.append(
                                            f"{pz_prefix}Score{pzfields_sep}{profile} = {pz_prefix}Score{pzfields_sep}{profile} + {criterion_score}"
                                        )

                                # PZFlag
                                if (
                                    f"{pz_prefix}Flag{pzfields_sep}{profile}"
                                    in list_of_pzfields
                                ):
                                    sql_set.append(
                                        f"{pz_prefix}Flag{pzfields_sep}{profile} = {pz_prefix}Flag{pzfields_sep}{profile} AND {criterion_flag_bool}"
                                    )

                                # PZClass
                                if (
                                    f"{pz_prefix}Class{pzfields_sep}{profile}"
                                    in list_of_pzfields
                                    and criterion_class is not None
                                ):
                                    sql_set.append(
                                        f" {pz_prefix}Class{pzfields_sep}{profile} = list_concat(list_distinct({pz_prefix}Class{pzfields_sep}{profile}), {criterion_class}) "
                                    )

                                # PZComment
                                if (
                                    f"{pz_prefix}Comment{pzfields_sep}{profile}"
                                    in list_of_pzfields
                                ):
                                    sql_set.append(f"""
                                            {pz_prefix}Comment{pzfields_sep}{profile} = 
                                                concat(
                                                    {pz_prefix}Comment{pzfields_sep}{profile},
                                                    CASE 
                                                        WHEN {pz_prefix}Comment{pzfields_sep}{profile}!=''
                                                        THEN ', '
                                                        ELSE ''
                                                    END,
                                                    '{criterion_comment}'
                                                )
                                        """)

                                # PZInfos
                                if (
                                    f"{pz_prefix}Infos{pzfields_sep}{profile}"
                                    in list_of_pzfields
                                ):
                                    sql_set.append(f"""
                                            {pz_prefix}Infos{pzfields_sep}{profile} = 
                                                concat(
                                                    {pz_prefix}Infos{pzfields_sep}{profile},
                                                    '{criterion_infos}'
                                                )
                                        """)
                                sql_set_option = ",".join(sql_set)

                                # Criterion and comparison
                                if sql_set_option:

                                    # Operation mode
                                    if criterion_mode in ["operation"]:

                                        # Check if value is a float
                                        try:

                                            # Test if criterion is a float
                                            float(criterion_value)

                                            # Query test cast as float
                                            query_test_cast = f"""
                                                SELECT "{annotation_view_name}"."{annotations_view_prefix}{annotation}"
                                                    FROM "{annotation_view_name}"
                                                    WHERE CAST("{annotation_view_name}"."{annotations_view_prefix}{annotation}" AS FLOAT) > 0
                                                LIMIT 1
                                            """
                                            self.execute_query(query_test_cast)

                                            sql_update = f"""
                                                UPDATE "{table_variants}"
                                                SET {sql_set_option}
                                                FROM (
                                                    SELECT *
                                                    FROM "{annotation_view_name}"
                                                    WHERE (
                                                        CAST("{annotation_view_name}"."{annotations_view_prefix}{annotation}" AS VARCHAR) NOT IN ('','.')
                                                        AND   CAST("{annotation_view_name}"."{annotations_view_prefix}{annotation}" AS FLOAT){comparison_map[criterion_type]}{criterion_value}
                                                        )
                                                    ) AS "{annotation_view_name}"
                                                WHERE ({" AND ".join(clause_join)})
                                                
                                            """
                                        # If not a float
                                        except:
                                            contains_option = ""
                                            if criterion_type == "contains":
                                                contains_option = ".*"
                                            sql_update = f"""
                                                UPDATE "{table_variants}"
                                                SET {sql_set_option}
                                                FROM (
                                                    SELECT *
                                                    FROM "{annotation_view_name}"
                                                    WHERE (
                                                    CAST("{annotation_view_name}"."{annotations_view_prefix}{annotation}" AS STRING) SIMILAR TO '{contains_option}{criterion_value}{contains_option}'
                                                        )
                                                    ) AS "{annotation_view_name}"
                                                WHERE ({" AND ".join(clause_join)})
                                                  
                                            """
                                        sql_queries.append(sql_update)

                                    # SQL mode
                                    elif criterion_mode in ["sql"]:

                                        sql_update = f"""
                                            UPDATE {table_variants}
                                            SET {sql_set_option}
                                            FROM (
                                                SELECT *
                                                FROM "{annotation_view_name}"
                                                WHERE ({criterion_sql})
                                                ) AS "{annotation_view_name}"
                                            WHERE ({" AND ".join(clause_join)})
                                        """
                                        sql_queries.append(sql_update)

                                    else:
                                        msg_err = f"Prioritization criterion mode failed (either 'operation' or 'sql')"
                                        log.error(msg_err)
                                        raise ValueError(msg_err)

                                else:
                                    log.warning(
                                        f"NO SQL SET option for '{annotation}' - '{criterion}'"
                                    )

                        # PZTags
                        if (
                            f"{pz_prefix}Tags{pzfields_sep}{profile}"
                            in list_of_pzfields
                        ):

                            # Create PZFalgs value
                            pztags_value = ""
                            pztags_sep_default = ","
                            pztags_sep = ""
                            for pzfield in pzfields:
                                if pzfield not in [f"{pz_prefix}Tags"]:
                                    if (
                                        f"{pzfield}{pzfields_sep}{profile}"
                                        in list_of_pzfields
                                    ):
                                        if pzfield in [f"{pz_prefix}Flag"]:
                                            pztags_value += f"""{pztags_sep}{pzfield}#', 
                                                CASE WHEN {pz_prefix}Flag{pzfields_sep}{profile}
                                                    THEN 'PASS'
                                                    ELSE 'FILTERED'
                                                END, '"""
                                        elif pzfield in [f"{pz_prefix}Class"]:
                                            pztags_value += f"""{pztags_sep}{pzfield}#', 
                                                CASE WHEN len({pz_prefix}Class{pzfields_sep}{profile}) > 0
                                                    THEN list_aggregate(list_distinct({pz_prefix}Class{pzfields_sep}{profile}), 'string_agg', ',')
                                                    ELSE '.'
                                                END, '"""
                                        else:
                                            pztags_value += f"{pztags_sep}{pzfield}#', {pzfield}{pzfields_sep}{profile}, '"
                                        pztags_sep = pztags_sep_default

                            # Add Query update for PZFlags
                            sql_update_pztags = f"""
                                UPDATE {table_variants}
                                SET INFO = concat(
                                        INFO,
                                        CASE WHEN INFO NOT in ('','.')
                                                THEN ';'
                                                ELSE ''
                                        END,
                                        '{pz_prefix}Tags{pzfields_sep}{profile}={pztags_value}'
                                    )
                                WHERE 1=1
                                """
                            sql_queries.append(sql_update_pztags)

                            # Add Query update for PZFlags for default
                            if profile == default_profile:
                                sql_update_pztags_default = f"""
                                UPDATE {table_variants}
                                SET INFO = concat(
                                        INFO,
                                        ';',
                                        '{pz_prefix}Tags={pztags_value}'
                                    )
                                    WHERE 1=1
                                """
                                sql_queries.append(sql_update_pztags_default)

                        log.info(f"""Profile '{profile}' - Prioritization... """)

                        # Chromosomes list
                        sql_uniq_chrom = f"""
                            SELECT DISTINCT "#CHROM"
                            FROM {table_variants}
                        """
                        chroms = self.get_query_to_df(sql_uniq_chrom)["#CHROM"].tolist()

                        for chrom in chroms:

                            log.debug(
                                f"""Profile '{profile}' - Prioritization query - Chromosome '{chrom}'... """
                            )

                            if sql_queries:

                                # Query num
                                num_query = 0

                                # For each query
                                for sql_query in sql_queries:

                                    # Query num
                                    num_query += 1

                                    sql_query_chrom = f"""
                                        {sql_query}
                                        AND {table_variants}."#CHROM" LIKE '{chrom}' 
                                    """
                                    log.debug(
                                        f"""Profile '{profile}' - Prioritization query - Chromosome '{chrom}' [{num_query}/{len(sql_queries)}]"""
                                    )
                                    # log.debug(
                                    #     f"""sql_query_chrom:\n{sql_query_chrom}"""
                                    # )
                                    self.execute_query(query=sql_query_chrom)

                        # Update INFO field
                        log.info(f"""Profile '{profile}' - Update... """)
                        sql_query_update = f"""
                            UPDATE {table_variants}
                            SET INFO =  
                                concat(
                                    CASE
                                        WHEN INFO NOT IN ('','.')
                                        THEN concat(INFO, ';')
                                        ELSE ''
                                    END
                                    {sql_set_info_option}
                                )
                        """
                        # log.debug(f"sql_query_update={sql_query_update}")
                        self.execute_query(query=sql_query_update)

                        # Remove annotations view for prioritization
                        self.remove_tables_or_views(tables=[annotation_view_name])

        else:

            log.warning(f"No profiles in parameters")

        # Remove added columns
        for added_column in added_columns:
            self.drop_column(column=added_column)

        return True
