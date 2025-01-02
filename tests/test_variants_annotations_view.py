# -*- coding: utf-8 -*-
"""
Tests

Usage:
pytest tests/

Coverage:
coverage run -m pytest tests/test_variants_transcripts.py -x -vv --log-cli-level=DEBUG --capture=tee-sys
coverage report --include=howard/* -m
"""

import logging as log
from tempfile import TemporaryDirectory
import pytest
import vcf
import os

from howard.functions.commons import remove_if_exists, get_file_format
from howard.objects.variants import Variants
from test_needed import tests_folder, tests_config, tests_data_folder


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

        # TEST 1
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
        assert annotations_view_select.shape == (7, 12)
        assert annotations_view_select.columns.to_list() == [
            "#CHROM",
            "POS",
            "REF",
            "ALT",
            "NS",
            "DP",
            "AA",
            "CLNSIG",
            "PREFIXCLNSIG",
            "CLNSIGSUFFIX",
            "SPiP_Alt",
            "SIFT",
        ]

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
            view=annotations_view_name, table="variants", fields=fields
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
        assert annotations_view_select.shape == (7, 12)
        assert annotations_view_select.columns.to_list() == [
            "#CHROM",
            "POS",
            "REF",
            "ALT",
            "NS",
            "DP",
            "AA",
            "CLNSIG",
            "PREFIXCLNSIG",
            "CLNSIGSUFFIX",
            "SPiP_Alt",
            "SIFT",
        ]

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
            view=annotations_view_name, table="variants", fields=fields, drop_view=True
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
        assert annotations_view_select.columns.to_list() == [
            "#CHROM",
            "POS",
            "REF",
            "ALT",
            "CLNSIG",
            "SIFT",
        ]

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
        assert annotations_view_select.columns.to_list() == [
            "#CHROM",
            "POS",
            "ID",
            "REF",
            "ALT",
            "FILTER",
            "CLNSIG",
            "SIFT",
        ]

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
        assert annotations_view_select.columns.to_list() == [
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
        assert annotations_view_select.columns.to_list() == [
            "#CHROM",
            "POS",
            "ID",
            "REF",
            "ALT",
            "FILTER",
            "CLNSIG",
            "SIFT",
        ]

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
        assert annotations_view_select.columns.to_list() == [
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
        assert annotations_view_select.columns.to_list() == [
            "#CHROM",
            "POS",
            "ID",
            "REF",
            "ALT",
            "FILTER",
            "CLNSIG",
            "SIFT",
        ]

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
            fields_needed=fields_needed,
            fields_not_exists=False,
            detect_type_list=True,
            drop_view=True,
            prefix=prefix,
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
        log.debug(annotations_view_select.columns.to_list())
        assert annotations_view_select.columns.to_list() == [
            "#CHROM",
            "POS",
            "ID",
            "REF",
            "ALT",
            "FILTER",
            "PREFIX_CLNSIG",
            "PREFIX_SIFT",
        ]

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
            fields_needed=fields_needed,
            fields_not_exists=False,
            detect_type_list=True,
            drop_view=True,
            prefix=prefix,
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
        assert annotations_view_select.columns.to_list() == [
            "#CHROM",
            "POS",
            "ID",
            "REF",
            "ALT",
            "FILTER",
            "PREFIX_CLNSIG",
            "PREFIX_SIFT",
        ]
