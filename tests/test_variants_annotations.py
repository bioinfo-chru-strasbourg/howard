# -*- coding: utf-8 -*-
"""
Tests

Usage:
pytest tests/

Coverage:
coverage run -m pytest tests/test_objects_variants.py -x -v --log-cli-level=INFO --capture=tee-sys
coverage report --include=howard/* -m 
"""

import logging as log
import os
import sys
from tempfile import TemporaryDirectory
import duckdb
import re
import Bio.bgzf as bgzf
import gzip
import pytest

from howard.functions.commons import *
from howard.objects.variants import Variants
from howard.functions.databases import *
from test_needed import *


def test_annotations():
    """
    This function tests the annotation process of a VCF file with multiple annotations.

    The function initializes input VCF and annotation files, constructs a parameter dictionary with different annotation types,
    creates a Variants object with the input file, parameter dictionary, and output file, and tests the annotation process.

    The function then checks the parameter dictionary of the Variants object, and tests the output VCF file for the presence of
    annotated variants using SQL queries.

    Finally, the function exports the output VCF file and checks if it is in the correct format with pyVCF.

    :raises AssertionError: If any of the tests fail.
    """

    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # Init files
        input_vcf = tests_data_folder + "/example.vcf.gz"

        annotation1 = database_files.get("parquet")
        annotation2 = database_files.get("example_vcf_gz")
        annotation3 = database_files.get("refgene_gz")

        output_vcf = f"{tmp_dir}/output.vcf.gz"

        # Config
        config = tests_config.copy()

        # Construct param dict
        param_annotation = {
            "parquet": {
                "annotations": {
                    annotation1: {"INFO": None},
                    annotation2: {"CLNSIG": "CLNSIG_new"},
                },
            },
            "bcftools": {
                "annotations": {
                    annotation2: {"CLNSIG": "CLNSIG_new_bcftools"},
                    annotation3: {"symbol": "gene"},
                },
            },
        }
        param = {"annotation": param_annotation, "assembly": "hg19"}

        # Create object
        variants = Variants(
            conn=None,
            input=input_vcf,
            output=output_vcf,
            config=config,
            param=param,
            load=True,
        )

        # Remove if output file exists
        remove_if_exists([output_vcf])

        # Annotation
        variants.annotation()

        # check param
        param_input = variants.get_param()
        expected_param = param

        assert param_input == expected_param

        # Check annotation1
        result = variants.get_query_to_df(
            "SELECT 1 AS count FROM variants WHERE \"#CHROM\" = 'chr1' AND POS = 28736 AND REF = 'A' AND ALT = 'C' AND INFO LIKE '%CLNSIG_new=%'"
        )
        assert len(result) == 1

        # Check annotation2
        result = variants.get_query_to_df(
            "SELECT 1 AS count FROM variants WHERE \"#CHROM\" = 'chr7' AND POS = 55249063 AND REF = 'G' AND ALT = 'A' AND INFO = 'DP=125;nci60=0.66;gene=EGFR,EGFR-AS1'"
        )
        assert len(result) == 1

        # Check if VCF is in correct format with pyVCF
        variants.export_output()
        try:
            vcf.Reader(filename=output_vcf)
        except:
            assert False
