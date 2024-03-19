# -*- coding: utf-8 -*-
"""
Tests

Usage:
pytest tests/

Coverage:
coverage run -m pytest . -x -v
coverage report --include=howard/* -m
"""

import logging as log
import os
import sys
import duckdb
import re
import Bio.bgzf as bgzf
import gzip
import pytest
import pandas as pd
from pandas.testing import assert_frame_equal
from unittest.mock import patch

from howard.objects.variants import Variants
from howard.functions.commons import *
from howard.tools.tools import *
from test_needed import *


def test_calculation_tsv():

    # Init files
    input_vcf = tests_data_folder + "/example.vcf.gz"
    output_vcf = "/tmp/output_file.tsv"
    config = {}
    calculations = "VARTYPE,NOMEN,TRIO"

    # prepare arguments for the query function
    args = argparse.Namespace(
        input=input_vcf,
        output=output_vcf,
        config=config,
        calculations=calculations,
        hgvs_field="hgvs",
        transcripts=None,
        show_calculations=False,
        trio_pedigree='{"father":"sample1", "mother":"sample2", "child":"sample3"}',
        calculation_config=None,
        arguments_dict=arguments_dict,
    )

    # Remove if output file exists
    remove_if_exists([output_vcf])

    # Query
    calculation(args)

    # Check output file exists
    assert os.path.exists(output_vcf)

    # read the contents of the actual output file
    with open(output_vcf, "r") as f:
        result_output_nb_lines = 0
        result_output_nb_variants = 0
        for line in f:
            result_output_nb_lines += 1
            if not line.startswith("#"):
                result_output_nb_variants += 1

    # Expected result
    expected_result_nb_lines = 8
    expected_result_nb_variants = 7

    # Compare
    assert result_output_nb_lines == expected_result_nb_lines
    assert result_output_nb_variants == expected_result_nb_variants

    # Create object
    variants = Variants(conn=None, input=output_vcf, config=config, load=True)

    # Check annotation
    result = variants.get_query_to_df(
        "SELECT INFO FROM variants WHERE INFO LIKE '%VARTYPE=%'"
    )
    assert len(result) == 7

    # Check annotation
    result = variants.get_query_to_df(
        "SELECT INFO FROM variants WHERE INFO LIKE '%trio=%'"
    )
    assert len(result) == 7

    # Check annotation
    result = variants.get_query_to_df(
        "SELECT INFO FROM variants WHERE INFO LIKE '%hgvs=%'"
    )
    assert len(result) == 0


def test_calculation_vcf():

    # Init files
    input_vcf = tests_data_folder + "/example.vcf.gz"
    output_vcf = "/tmp/output_file.vcf"
    config = {}
    calculations = "VARTYPE,NOMEN,TRIO"

    # prepare arguments for the query function
    args = argparse.Namespace(
        input=input_vcf,
        output=output_vcf,
        config=config,
        calculations=calculations,
        hgvs_field="hgvs",
        transcripts=None,
        show_calculations=False,
        trio_pedigree='{"father":"sample1", "mother":"sample2", "child":"sample3"}',
        calculation_config=None,
        arguments_dict=arguments_dict,
    )

    # Remove if output file exists
    remove_if_exists([output_vcf])

    # Query
    calculation(args)

    # Check output file exists
    assert os.path.exists(output_vcf)

    # read the contents of the actual output file
    with open(output_vcf, "r") as f:
        result_output_nb_lines = 0
        result_output_nb_variants = 0
        for line in f:
            result_output_nb_lines += 1
            if not line.startswith("#"):
                result_output_nb_variants += 1

    # Expected result
    expected_result_nb_lines = 72
    expected_result_nb_variants = 7

    # Compare
    assert result_output_nb_lines == expected_result_nb_lines
    assert result_output_nb_variants == expected_result_nb_variants

    # Create object
    variants = Variants(conn=None, input=output_vcf, config=config, load=True)

    # Check annotation
    result = variants.get_query_to_df(
        "SELECT INFO FROM variants WHERE INFO LIKE '%VARTYPE=%'"
    )
    assert len(result) == 7

    # Check annotation
    result = variants.get_query_to_df(
        "SELECT INFO FROM variants WHERE INFO LIKE '%trio=%'"
    )
    assert len(result) == 7

    # Check annotation
    result = variants.get_query_to_df(
        "SELECT INFO FROM variants WHERE INFO LIKE '%hgvs=%'"
    )
    assert len(result) == 0
