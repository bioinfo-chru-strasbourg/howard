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



def test_prioritization_tsv():

    # Init files
    input_vcf = tests_data_folder + "/example.vcf.gz"
    output_vcf = "/tmp/output_file.tsv"
    config = {}
    prioritizations = "default,GERMLINE"
    prioritization_config = tests_data_folder + "/prioritization_profiles.json"

    # prepare arguments for the query function
    args = argparse.Namespace(
        input = input_vcf,
        output = output_vcf,
        config = config,
        prioritizations = prioritizations,
        prioritization_config = prioritization_config,
        pzfields = "PZScore,PZFlag",
        default_profile = "default",
        prioritization_score_mode = "HOWARD",
        arguments_dict = arguments_dict
    )

    # Remove if output file exists
    remove_if_exists([output_vcf])

    # Query
    prioritization(args)

    # read the contents of the actual output file
    with open(output_vcf, 'r') as f:
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
    result = variants.get_query_to_df("SELECT INFO FROM variants WHERE INFO LIKE '%PZScore=%'")
    assert len(result) == 7


def test_prioritization_vcf():

    # Init files
    input_vcf = tests_data_folder + "/example.vcf.gz"
    output_vcf = "/tmp/output_file.vcf"
    config = {}
    prioritizations = "default,GERMLINE"
    prioritization_config = tests_data_folder + "/prioritization_profiles.json"

    # prepare arguments for the query function
    args = argparse.Namespace(
        input = input_vcf,
        output = output_vcf,
        config = config,
        prioritizations = prioritizations,
        prioritization_config = prioritization_config,
        pzfields = "PZScore,PZFlag",
        default_profile = "default",
        prioritization_score_mode = "HOWARD",
        arguments_dict = arguments_dict
    )

    # Remove if output file exists
    remove_if_exists([output_vcf])

    # Query
    prioritization(args)

    # read the contents of the actual output file
    with open(output_vcf, 'r') as f:
        result_output_nb_lines = 0
        result_output_nb_variants = 0
        for line in f:
            result_output_nb_lines += 1
            if not line.startswith("#"):
                result_output_nb_variants += 1

    # Expected result
    expected_result_nb_lines = 66
    expected_result_nb_variants = 7

    # Compare
    assert result_output_nb_lines == expected_result_nb_lines
    assert result_output_nb_variants == expected_result_nb_variants

    # Create object
    variants = Variants(conn=None, input=output_vcf, config=config, load=True)
    
    # Check annotation
    result = variants.get_query_to_df("SELECT INFO FROM variants WHERE INFO LIKE '%PZScore=%'")
    assert len(result) == 7