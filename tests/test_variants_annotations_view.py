# -*- coding: utf-8 -*-
"""
Tests

Usage:
pytest tests/

Coverage:
coverage run -m pytest tests/test_variants_annotations_view.py -x -vv --log-cli-level=DEBUG --capture=tee-sys
coverage report --include=howard/* -m
"""

import logging as log
from tempfile import TemporaryDirectory
import pytest  # type: ignore


from howard.objects.variants import Variants
from test_needed import tests_folder, tests_config, tests_data_folder
from howard.functions.commons import set_log_level


def test_create_annotations_view_empty_info():
    """ """

    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # Init files
        input_vcf = tests_data_folder + "/example.empty_info.vcf"
        output_vcf = f"{tmp_dir}/output.vcf.gz"

        # config dict
        config = tests_config

        # Construct param dict
        param = {}

        # Create object
        variants = Variants(
            conn=None,
            input=input_vcf,
            output=output_vcf,
            config=config,
            param=param,
            load=True,
        )

        annotations_view_name = "annotations_view_test"

        # TEST 0
        ##########

        # Create annotations view
        annotations_view_name_result = variants.create_annotations_view(
            view=annotations_view_name,
            table="variants",
            view_type="view",
            view_mode="full",
            fields=None,
            info_prefix_column="",
            info_struct_column="INFOS",
            sample_struct_column="SAMPLES",
        )

        # Check annotations view name
        assert annotations_view_name == annotations_view_name_result

        # Check annotations_view content
        annotations_view_select = variants.get_query_to_df(
            query=f"""
            SELECT *
            FROM {annotations_view_name}
            LIMIT 100
            """
        )
        log.debug(f"annotations_view_select={annotations_view_select}")
        # Check shape
        assert annotations_view_select.shape == (10, 11)
        assert sorted(set(annotations_view_select.columns.to_list())) == sorted(
            set(
                [
                    "#CHROM",
                    "POS",
                    "REF",
                    "ALT",
                    "genome",
                    "uniprot_id",
                    "protein_variant",
                    "am_pathogenicity",
                    "transcript_id",
                    "am_class",
                    "INFOS",
                ]
            )
        )

        # TEST 1
        ##########

        # Create annotations view
        annotations_view_name_result = variants.create_annotations_view(
            view=annotations_view_name,
            table="variants",
            view_type="view",
            view_mode="explore",
            fields=None,
            info_prefix_column="",
            info_struct_column="INFOS",
            sample_struct_column="SAMPLES",
        )

        # Check annotations view name
        assert annotations_view_name == annotations_view_name_result

        # Check annotations_view content
        annotations_view_select = variants.get_query_to_df(
            query=f"""
            SELECT *
            FROM {annotations_view_name}
            LIMIT 100
            """
        )
        log.debug(f"annotations_view_select={annotations_view_select}")
        # Check shape
        assert annotations_view_select.shape == (10, 11)
        assert sorted(set(annotations_view_select.columns.to_list())) == sorted(
            set(
                [
                    "#CHROM",
                    "POS",
                    "REF",
                    "ALT",
                    "genome",
                    "uniprot_id",
                    "protein_variant",
                    "am_pathogenicity",
                    "transcript_id",
                    "am_class",
                    "INFOS",
                ]
            )
        )


def test_create_annotations_view_no_info():
    """ """

    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # Init files
        input_vcf = tests_data_folder + "/example.no_info.vcf"
        output_vcf = f"{tmp_dir}/output.vcf.gz"

        # config dict
        config = tests_config

        # Construct param dict
        param = {}

        # Create object
        variants = Variants(
            conn=None,
            input=input_vcf,
            output=output_vcf,
            config=config,
            param=param,
            load=True,
        )

        annotations_view_name = "annotations_view_test"

        # TEST 0
        ##########

        # Create annotations view
        annotations_view_name_result = variants.create_annotations_view(
            view=annotations_view_name,
            table="variants",
            view_type="view",
            view_mode="full",
            fields=None,
            info_prefix_column="",
            info_struct_column="INFOS",
            sample_struct_column="SAMPLES",
        )

        # Check annotations view name
        assert annotations_view_name == annotations_view_name_result

        # Check annotations_view content
        annotations_view_select = variants.get_query_to_df(
            query=f"""
            SELECT *
            FROM {annotations_view_name}
            LIMIT 100
            """
        )
        log.debug(f"annotations_view_select={annotations_view_select}")
        # Check shape
        assert annotations_view_select.shape == (10, 11)
        assert sorted(set(annotations_view_select.columns.to_list())) == sorted(
            set(
                [
                    "#CHROM",
                    "POS",
                    "REF",
                    "ALT",
                    "genome",
                    "uniprot_id",
                    "protein_variant",
                    "am_pathogenicity",
                    "transcript_id",
                    "am_class",
                    "INFOS",
                ]
            )
        )

        # TEST 1
        ##########

        # Create annotations view
        annotations_view_name_result = variants.create_annotations_view(
            view=annotations_view_name,
            table="variants",
            view_type="view",
            view_mode="explore",
            fields=None,
            info_prefix_column="",
            info_struct_column="INFOS",
            sample_struct_column="SAMPLES",
        )

        # Check annotations view name
        assert annotations_view_name == annotations_view_name_result

        # Check annotations_view content
        annotations_view_select = variants.get_query_to_df(
            query=f"""
            SELECT *
            FROM {annotations_view_name}
            LIMIT 100
            """
        )
        log.debug(f"annotations_view_select={annotations_view_select}")
        # Check shape
        assert annotations_view_select.shape == (10, 11)
        assert sorted(set(annotations_view_select.columns.to_list())) == sorted(
            set(
                [
                    "#CHROM",
                    "POS",
                    "REF",
                    "ALT",
                    "genome",
                    "uniprot_id",
                    "protein_variant",
                    "am_pathogenicity",
                    "transcript_id",
                    "am_class",
                    "INFOS",
                ]
            )
        )


def test_create_annotations_view_chrom_pos_ref_alt():
    """ """

    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # Init files
        input_vcf = tests_data_folder + "/example.chrom.pos.ref.alt.vcf"
        output_vcf = f"{tmp_dir}/output.vcf.gz"

        # config dict
        config = tests_config
        # config["access"] = "RO"

        # Construct param dict
        param = {}

        # Create object
        variants = Variants(
            conn=None,
            input=input_vcf,
            output=output_vcf,
            config=config,
            param=param,
            load=True,
        )

        annotations_view_name = "annotations_view_test"

        # TEST 0
        ##########

        # Create annotations view
        annotations_view_name_result = variants.create_annotations_view(
            view=annotations_view_name, table="variants", fields=None
        )

        # Check annotations view name
        assert annotations_view_name == annotations_view_name_result

        # Check annotations_view content
        annotations_view_select = variants.get_query_to_df(
            query=f"""
            SELECT *
            FROM {annotations_view_name}
            LIMIT 100
            """
        )
        # Check shape
        assert annotations_view_select.shape == (7, 4)
        assert sorted(set(annotations_view_select.columns.to_list())) == sorted(
            set(
                [
                    "#CHROM",
                    "POS",
                    "REF",
                    "ALT",
                ]
            )
        )

        # TEST 1
        ##########
        # Generates columns from fields
        # Not dropped! Same than before

        # Create annotations view
        annotations_view_name_result = variants.create_annotations_view(
            view=annotations_view_name,
            table="variants",
            fields=None,
            info_prefix_column="",
        )

        # Check annotations view name
        assert annotations_view_name == annotations_view_name_result

        # Check annotations_view content
        annotations_view_select = variants.get_query_to_df(
            query=f"""
            SELECT *
            FROM {annotations_view_name}
            LIMIT 100
            """
        )
        # Check shape
        assert annotations_view_select.shape == (7, 4)
        assert sorted(set(annotations_view_select.columns.to_list())) == sorted(
            set(
                [
                    "#CHROM",
                    "POS",
                    "REF",
                    "ALT",
                ]
            )
        )

        # TEST 2
        ##########
        # Add specific fields
        # Without drop

        # Create annotations view
        fields = [
            "CLNSIG",
            "SIFT",
        ]
        annotations_view_name_result = variants.create_annotations_view(
            view=annotations_view_name,
            table="variants",
            fields=fields,
            info_prefix_column="",
        )

        # Check annotations view name
        assert annotations_view_name == annotations_view_name_result

        # Check annotations_view content
        annotations_view_select = variants.get_query_to_df(
            query=f"""
            SELECT *
            FROM {annotations_view_name}
            LIMIT 100
            """
        )
        # Check shape
        assert annotations_view_select.shape == (7, 4)
        assert sorted(set(annotations_view_select.columns.to_list())) == sorted(
            set(
                [
                    "#CHROM",
                    "POS",
                    "REF",
                    "ALT",
                ]
            )
        )

        # TEST 3
        ##########
        # Add specific fields
        # With drop

        # Create annotations view
        fields = [
            "CLNSIG",
            "SIFT",
        ]
        annotations_view_name_result = variants.create_annotations_view(
            view=annotations_view_name,
            table="variants",
            fields=fields,
            info_prefix_column="",
            drop_view=True,
        )

        # Check annotations view name
        assert annotations_view_name == annotations_view_name_result

        # Check annotations_view content
        annotations_view_select = variants.get_query_to_df(
            query=f"""
            SELECT *
            FROM {annotations_view_name}
            LIMIT 100
            """
        )
        # Check shape
        assert annotations_view_select.shape == (7, 6)
        assert sorted(set(annotations_view_select.columns.to_list())) == sorted(
            set(
                [
                    "#CHROM",
                    "POS",
                    "REF",
                    "ALT",
                    "CLNSIG",
                    "SIFT",
                ]
            )
        )

        # TEST 4
        ##########
        # Add specific fields
        # Add specific fields needed
        # With drop

        # Create annotations view
        fields = [
            "CLNSIG",
            "SIFT",
        ]
        fields_needed = ["#CHROM", "POS", "ID", "REF", "ALT", "FILTER"]

        with pytest.raises(ValueError) as e:
            annotations_view_name_result = variants.create_annotations_view(
                view=annotations_view_name,
                table="variants",
                fields=fields,
                info_prefix_column="",
                fields_needed=fields_needed,
                drop_view=True,
            )
        assert str(e.value) == f"Field 'ID' is needed, but not in file"

        # TEST 5
        ##########
        # Add INFO struct column

        # Create annotations view
        fields = ["CLNSIG", "SIFT", "FIELD_THAT_NOT_EXISTS"]
        fields_needed = ["#CHROM", "POS", "ID", "REF", "ALT", "FILTER"]
        info_struct_column = "INFOS"
        sample_struct_column = "SAMPLES"
        annotations_view_name_result = variants.create_annotations_view(
            view=annotations_view_name,
            table="variants",
            fields=fields,
            info_prefix_column="INFOS_",
            info_struct_column=info_struct_column,
            sample_struct_column=sample_struct_column,
            # fields_needed=fields_needed,
            # fields_needed=None,
            fields_needed_all=True,
            fields_not_exists=False,
            detect_type_list=True,
            drop_view=True,
        )

        # Check annotations view name
        assert annotations_view_name == annotations_view_name_result

        # Check annotations_view content
        annotations_view_select = variants.get_query_to_df(
            query=f"""
            SELECT *
            FROM {annotations_view_name}
            LIMIT 100
            """
        )
        # Check shape
        assert annotations_view_select.shape == (7, 8)
        assert sorted(set(annotations_view_select.columns.to_list())) == sorted(
            set(
                [
                    "#CHROM",
                    "POS",
                    "REF",
                    "ALT",
                    "INFO",
                    "INFOS_CLNSIG",
                    "INFOS_SIFT",
                    "INFOS",
                ]
            )
        )


def test_create_annotations_view():
    """ """

    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # Init files
        input_vcf = tests_data_folder + "/example.annotation_names.vcf.gz"
        output_vcf = f"{tmp_dir}/output.vcf.gz"

        # config dict
        config = tests_config
        # config["access"] = "RO"

        # Construct param dict
        param = {}

        # Create object
        variants = Variants(
            conn=None,
            input=input_vcf,
            output=output_vcf,
            config=config,
            param=param,
            load=True,
        )

        annotations_view_name = "annotations_view_test"

        # TEST 0
        ##########

        # Create annotations view
        annotations_view_name_result = variants.create_annotations_view(
            view=annotations_view_name, table="variants", fields=None
        )

        # Check annotations view name
        assert annotations_view_name == annotations_view_name_result

        # Check annotations_view content
        annotations_view_select = variants.get_query_to_df(
            query=f"""
            SELECT *
            FROM {annotations_view_name}
            LIMIT 100
            """
        )
        # Check shape
        assert annotations_view_select.shape == (7, 4)
        assert sorted(set(annotations_view_select.columns.to_list())) == sorted(
            set(
                [
                    "#CHROM",
                    "POS",
                    "REF",
                    "ALT",
                ]
            )
        )

        # TEST 1
        ##########
        # Generates columns from fields
        # Not dropped! Same than before

        # Create annotations view
        annotations_view_name_result = variants.create_annotations_view(
            view=annotations_view_name,
            table="variants",
            fields=None,
            info_prefix_column="",
        )

        # Check annotations view name
        assert annotations_view_name == annotations_view_name_result

        # Check annotations_view content
        annotations_view_select = variants.get_query_to_df(
            query=f"""
            SELECT *
            FROM {annotations_view_name}
            LIMIT 100
            """
        )
        # Check shape
        assert annotations_view_select.shape == (7, 4)
        assert sorted(set(annotations_view_select.columns.to_list())) == sorted(
            set(
                [
                    "#CHROM",
                    "POS",
                    "REF",
                    "ALT",
                ]
            )
        )

        # TEST 2
        ##########
        # Add specific fields
        # Without drop

        # Create annotations view
        fields = [
            "CLNSIG",
            "SIFT",
        ]
        annotations_view_name_result = variants.create_annotations_view(
            view=annotations_view_name,
            table="variants",
            fields=fields,
            info_prefix_column="",
        )

        # Check annotations view name
        assert annotations_view_name == annotations_view_name_result

        # Check annotations_view content
        annotations_view_select = variants.get_query_to_df(
            query=f"""
            SELECT *
            FROM {annotations_view_name}
            LIMIT 100
            """
        )
        # Check shape
        assert annotations_view_select.shape == (7, 4)
        assert sorted(set(annotations_view_select.columns.to_list())) == sorted(
            set(
                [
                    "#CHROM",
                    "POS",
                    "REF",
                    "ALT",
                ]
            )
        )

        # TEST 3
        ##########
        # Add specific fields
        # With drop

        # Create annotations view
        fields = [
            "CLNSIG",
            "SIFT",
        ]
        annotations_view_name_result = variants.create_annotations_view(
            view=annotations_view_name,
            table="variants",
            fields=fields,
            info_prefix_column="",
            drop_view=True,
        )

        # Check annotations view name
        assert annotations_view_name == annotations_view_name_result

        # Check annotations_view content
        annotations_view_select = variants.get_query_to_df(
            query=f"""
            SELECT *
            FROM {annotations_view_name}
            LIMIT 100
            """
        )
        # Check shape
        assert annotations_view_select.shape == (7, 6)
        assert sorted(set(annotations_view_select.columns.to_list())) == sorted(
            set(
                [
                    "#CHROM",
                    "POS",
                    "REF",
                    "ALT",
                    "CLNSIG",
                    "SIFT",
                ]
            )
        )

        # TEST 4
        ##########
        # Add specific fields
        # Add specific fields needed
        # With drop

        # Create annotations view
        fields = [
            "CLNSIG",
            "SIFT",
        ]
        fields_needed = ["#CHROM", "POS", "ID", "REF", "ALT", "FILTER"]
        annotations_view_name_result = variants.create_annotations_view(
            view=annotations_view_name,
            table="variants",
            fields=fields,
            info_prefix_column="",
            fields_needed=fields_needed,
            drop_view=True,
        )

        # Check annotations view name
        assert annotations_view_name == annotations_view_name_result

        # Check annotations_view content
        annotations_view_select = variants.get_query_to_df(
            query=f"""
            SELECT *
            FROM {annotations_view_name}
            LIMIT 100
            """
        )
        # Check shape
        assert annotations_view_select.shape == (7, 8)
        assert sorted(set(annotations_view_select.columns.to_list())) == sorted(
            set(
                [
                    "#CHROM",
                    "POS",
                    "ID",
                    "REF",
                    "ALT",
                    "FILTER",
                    "CLNSIG",
                    "SIFT",
                ]
            )
        )

        # TEST 5
        ##########
        # Add specific fields
        # Add specific fields needed as all fields in table
        # With drop

        # Create annotations view
        fields = [
            "CLNSIG",
            "SIFT",
        ]
        fields_needed = None
        fields_needed_all = True
        annotations_view_name_result = variants.create_annotations_view(
            view=annotations_view_name,
            table="variants",
            fields=fields,
            info_prefix_column="",
            fields_needed=fields_needed,
            fields_needed_all=fields_needed_all,
            drop_view=True,
        )

        # Check annotations view name
        assert annotations_view_name == annotations_view_name_result

        # Check annotations_view content
        annotations_view_select = variants.get_query_to_df(
            query=f"""
            SELECT *
            FROM {annotations_view_name}
            LIMIT 100
            """
        )
        # Check shape
        assert annotations_view_select.shape == (7, 15)
        assert sorted(set(annotations_view_select.columns.to_list())) == sorted(
            set(
                [
                    "#CHROM",
                    "POS",
                    "ID",
                    "REF",
                    "ALT",
                    "QUAL",
                    "FILTER",
                    "INFO",
                    "FORMAT",
                    "sample1",
                    "sample2",
                    "sample3",
                    "sample4",
                    "CLNSIG",
                    "SIFT",
                ]
            )
        )

        # TEST 6
        ##########
        # Add specific fields
        # Add specific fields needed
        # Detect field type as list
        # With drop

        # Create annotations view
        fields = [
            "CLNSIG",
            "SIFT",
        ]
        fields_needed = ["#CHROM", "POS", "ID", "REF", "ALT", "FILTER"]
        annotations_view_name_result = variants.create_annotations_view(
            view=annotations_view_name,
            table="variants",
            fields=fields,
            info_prefix_column="",
            fields_needed=fields_needed,
            detect_type_list=True,
            drop_view=True,
        )

        # Check annotations view name
        assert annotations_view_name == annotations_view_name_result

        # Check annotations_view content
        annotations_view_select = variants.get_query_to_df(
            query=f"""
            SELECT *
            FROM {annotations_view_name}
            LIMIT 100
            """
        )
        # Check shape
        assert annotations_view_select.shape == (7, 8)
        assert sorted(set(annotations_view_select.columns.to_list())) == sorted(
            set(
                [
                    "#CHROM",
                    "POS",
                    "ID",
                    "REF",
                    "ALT",
                    "FILTER",
                    "CLNSIG",
                    "SIFT",
                ]
            )
        )

        # Check annotations_view content
        # CHeck row with #CHROM = chr1 and position 69101, if column SIFT is an array of 2 value [D, P]
        annotations_view_select = variants.get_query_to_df(
            query=f"""
            SELECT *
            FROM {annotations_view_name}
            WHERE "#CHROM" = 'chr1' AND POS = 69101
            """
        )
        # Check shape
        assert annotations_view_select.shape == (1, 8)
        # Compare list length
        assert len(annotations_view_select["SIFT"].values[0]) == 2
        # Compare list content
        assert all(
            item in annotations_view_select["SIFT"].values[0] for item in ["D", "P"]
        )

        # TEST 7
        ##########
        # Add specific fields with one does NOT exists
        # Add specific fields needed
        # Detect field type as list
        # With drop

        # Create annotations view
        fields = ["CLNSIG", "SIFT", "FIELD_THAT_NOT_EXISTS"]
        fields_needed = ["#CHROM", "POS", "ID", "REF", "ALT", "FILTER"]
        annotations_view_name_result = variants.create_annotations_view(
            view=annotations_view_name,
            table="variants",
            fields=fields,
            info_prefix_column="",
            fields_needed=fields_needed,
            detect_type_list=True,
            drop_view=True,
        )

        # Check annotations view name
        assert annotations_view_name == annotations_view_name_result

        # Check annotations_view content
        annotations_view_select = variants.get_query_to_df(
            query=f"""
            SELECT *
            FROM {annotations_view_name}
            LIMIT 100
            """
        )
        # Check shape
        assert annotations_view_select.shape == (7, 9)
        assert sorted(set(annotations_view_select.columns.to_list())) == sorted(
            set(
                [
                    "#CHROM",
                    "POS",
                    "ID",
                    "REF",
                    "ALT",
                    "FILTER",
                    "CLNSIG",
                    "SIFT",
                    "FIELD_THAT_NOT_EXISTS",
                ]
            )
        )

        # TEST 8
        ##########
        # Add specific fields without one does NOT exists
        # Add specific fields needed
        # Detect field type as list
        # With drop

        # Create annotations view
        fields = ["CLNSIG", "SIFT", "FIELD_THAT_NOT_EXISTS"]
        fields_needed = ["#CHROM", "POS", "ID", "REF", "ALT", "FILTER"]
        annotations_view_name_result = variants.create_annotations_view(
            view=annotations_view_name,
            table="variants",
            fields=fields,
            info_prefix_column="",
            fields_needed=fields_needed,
            fields_not_exists=False,
            detect_type_list=True,
            drop_view=True,
        )

        # Check annotations view name
        assert annotations_view_name == annotations_view_name_result

        # Check annotations_view content
        annotations_view_select = variants.get_query_to_df(
            query=f"""
            SELECT *
            FROM {annotations_view_name}
            LIMIT 100
            """
        )
        # Check shape
        assert annotations_view_select.shape == (7, 8)
        assert sorted(set(annotations_view_select.columns.to_list())) == sorted(
            set(
                [
                    "#CHROM",
                    "POS",
                    "ID",
                    "REF",
                    "ALT",
                    "FILTER",
                    "CLNSIG",
                    "SIFT",
                ]
            )
        )

        # TEST 9
        ##########
        # Add specific fields without one does NOT exists
        # Add specific fields needed
        # Detect field type as list
        # With drop
        # With prefix

        # Create annotations view
        fields = ["CLNSIG", "SIFT", "FIELD_THAT_NOT_EXISTS"]
        fields_needed = ["#CHROM", "POS", "ID", "REF", "ALT", "FILTER"]
        prefix = "PREFIX_"
        annotations_view_name_result = variants.create_annotations_view(
            view=annotations_view_name,
            table="variants",
            fields=fields,
            info_prefix_column=prefix,
            fields_needed=fields_needed,
            fields_not_exists=False,
            detect_type_list=True,
            drop_view=True,
        )

        # Check annotations view name
        assert annotations_view_name == annotations_view_name_result

        # Check annotations_view content
        annotations_view_select = variants.get_query_to_df(
            query=f"""
            SELECT *
            FROM {annotations_view_name}
            LIMIT 100
            """
        )
        # Check shape
        assert annotations_view_select.shape == (7, 8)
        assert sorted(set(annotations_view_select.columns.to_list())) == sorted(
            set(
                [
                    "#CHROM",
                    "POS",
                    "ID",
                    "REF",
                    "ALT",
                    "FILTER",
                    "PREFIX_CLNSIG",
                    "PREFIX_SIFT",
                ]
            )
        )

        # TEST 10
        ##########
        # Add specific fields without one does NOT exists
        # Add specific fields needed
        # Detect field type as list
        # With drop
        # With prefix
        # Limit 2

        # Create annotations view
        fields = ["CLNSIG", "SIFT", "FIELD_THAT_NOT_EXISTS"]
        fields_needed = ["#CHROM", "POS", "ID", "REF", "ALT", "FILTER"]
        prefix = "PREFIX_"
        limit = 2
        annotations_view_name_result = variants.create_annotations_view(
            view=annotations_view_name,
            table="variants",
            fields=fields,
            info_prefix_column=prefix,
            fields_needed=fields_needed,
            fields_not_exists=False,
            detect_type_list=True,
            drop_view=True,
            limit=limit,
        )

        # Check annotations view name
        assert annotations_view_name == annotations_view_name_result

        # Check annotations_view content
        annotations_view_select = variants.get_query_to_df(
            query=f"""
            SELECT *
            FROM {annotations_view_name}
            LIMIT 100
            """
        )
        # Check shape
        assert annotations_view_select.shape == (2, 8)
        assert sorted(set(annotations_view_select.columns.to_list())) == sorted(
            set(
                [
                    "#CHROM",
                    "POS",
                    "ID",
                    "REF",
                    "ALT",
                    "FILTER",
                    "PREFIX_CLNSIG",
                    "PREFIX_SIFT",
                ]
            )
        )

        # TEST 11
        ##########
        # Add INFO struct column

        # Create annotations view
        fields = ["CLNSIG", "SIFT", "FIELD_THAT_NOT_EXISTS"]
        fields_needed = ["#CHROM", "POS", "ID", "REF", "ALT", "FILTER"]
        info_struct_column = "INFOS"
        annotations_view_name_result = variants.create_annotations_view(
            view=annotations_view_name,
            table="variants",
            fields=fields,
            info_struct_column=info_struct_column,
            fields_needed=fields_needed,
            fields_not_exists=False,
            detect_type_list=True,
            drop_view=True,
        )

        # Check annotations view name
        assert annotations_view_name == annotations_view_name_result

        # Check annotations_view content
        annotations_view_select = variants.get_query_to_df(
            query=f"""
            SELECT *
            FROM {annotations_view_name}
            LIMIT 100
            """
        )
        # log.debug(annotations_view_select)
        # Check shape
        assert annotations_view_select.shape == (7, 7)
        assert sorted(set(annotations_view_select.columns.to_list())) == sorted(
            set(
                [
                    "#CHROM",
                    "POS",
                    "ID",
                    "REF",
                    "ALT",
                    "FILTER",
                    "INFOS",
                ]
            )
        )

        # TEST 12
        ##########
        # Add INFO struct column

        # Create annotations view
        fields = ["CLNSIG", "SIFT", "FIELD_THAT_NOT_EXISTS"]
        fields_needed = ["#CHROM", "POS", "ID", "REF", "ALT", "FILTER"]
        info_struct_column = "INFOS"
        sample_struct_column = "SAMPLES"
        annotations_view_name_result = variants.create_annotations_view(
            view=annotations_view_name,
            table="variants",
            fields=fields,
            info_struct_column=info_struct_column,
            sample_struct_column=sample_struct_column,
            # fields_needed=fields_needed,
            fields_needed=None,
            fields_needed_all=True,
            fields_not_exists=False,
            detect_type_list=True,
            drop_view=True,
        )

        # Check annotations view name
        assert annotations_view_name == annotations_view_name_result

        # Check annotations_view content
        annotations_view_select = variants.get_query_to_df(
            query=f"""
            SELECT *
            FROM {annotations_view_name}
            LIMIT 100
            """
        )
        # log.debug(annotations_view_select)
        # Check shape
        assert annotations_view_select.shape == (7, 15)
        assert sorted(set(annotations_view_select.columns.to_list())) == sorted(
            set(
                [
                    "#CHROM",
                    "POS",
                    "ID",
                    "REF",
                    "ALT",
                    "QUAL",
                    "FILTER",
                    "INFO",
                    "FORMAT",
                    "sample1",
                    "sample2",
                    "sample3",
                    "sample4",
                    "INFOS",
                    "SAMPLES",
                ]
            )
        )

        # Check struct
        annotations_view_select = variants.get_query_to_df(
            query=f"""
            SELECT "#CHROM", POS, REF, ALT, FORMAT, SAMPLES.sample1, SAMPLES.sample1.AD[2]/(SAMPLES.sample1.AD[1]+SAMPLES.sample1.AD[2]) AS 'sample1_VAF'
            FROM {annotations_view_name}
            WHERE SAMPLES.sample1.GQ > 90
              AND SAMPLES.sample1.DP > 300
              AND sample1_VAF >= 0
            LIMIT 100
            """
        )
        # log.debug(annotations_view_select.to_string())
        # Check shape
        assert annotations_view_select.shape == (6, 7)
        assert sorted(set(annotations_view_select.columns.to_list())) == sorted(
            set(
                [
                    "#CHROM",
                    "POS",
                    "REF",
                    "ALT",
                    "FORMAT",
                    "sample1",
                    "sample1_VAF",
                ]
            )
        )
