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


def test_process():

    # Init files
    input_vcf = tests_data_folder + "/example.vcf.gz"
    output_vcf = "/tmp/output_file.vcf"
    config = {}
    param = "{}"
    annotations = database_files.get("parquet")
    calculations = "VARTYPE"
    prioritizations = "default"
    input_query = None

    # prepare arguments for the query function
    args = argparse.Namespace(
        input=input_vcf,
        output=output_vcf,
        config=config,
        param=param,
        annotations=annotations,
        calculations=calculations,
        prioritizations=prioritizations,
        query=input_query,
        explode_infos=False,
        explode_infos_prefix="",
        explode_infos_fields="*",
        include_header=False,
        arguments_dict=arguments_dict,
    )

    # Remove if output file exists
    remove_if_exists([output_vcf])

    # Query
    process(args)

    # Create object
    variants = Variants(conn=None, input=output_vcf, config=config, load=True)

    # Check annotation
    result = variants.get_query_to_df(
        "SELECT INFO FROM variants WHERE INFO LIKE '%VARTYPE=%'"
    )
    assert len(result) == 7

    # Check annotation
    result = variants.get_query_to_df(
        "SELECT INFO FROM variants WHERE INFO LIKE '%nci60=%'"
    )
    assert len(result) == 1


def test_process_with_param_file():

    # Init files
    input_vcf = tests_data_folder + "/example.vcf.gz"
    output_vcf = "/tmp/output_file.vcf"
    config = {}
    param = tests_folder + "/data/param.snpeff_hgvs.json"
    annotations = database_files.get("parquet")
    calculations = "VARTYPE"
    prioritizations = "default"
    input_query = None

    # prepare arguments for the query function
    args = argparse.Namespace(
        input=input_vcf,
        output=output_vcf,
        config=config,
        param=param,
        annotations=annotations,
        calculations=calculations,
        prioritizations=prioritizations,
        query=input_query,
        explode_infos=False,
        explode_infos_prefix="",
        explode_infos_fields="*",
        include_header=False,
        arguments_dict=arguments_dict,
    )

    # Remove if output file exists
    remove_if_exists([output_vcf])

    # Query
    process(args)

    # Create object
    variants = Variants(conn=None, input=output_vcf, config=config, load=True)

    # Check annotation
    result = variants.get_query_to_df(
        "SELECT INFO FROM variants WHERE INFO LIKE '%VARTYPE=%'"
    )
    assert len(result) == 7

    # Check annotation
    result = variants.get_query_to_df(
        "SELECT INFO FROM variants WHERE INFO LIKE '%nci60=%'"
    )
    assert len(result) == 1


def test_process_with_query():

    # Init files
    input_vcf = tests_data_folder + "/example.vcf.gz"
    output_vcf = "/tmp/output_file.tsv"
    config = {}
    param = "{}"
    annotations = database_files.get("parquet")
    calculations = "VARTYPE"
    prioritizations = "default"
    input_query = "SELECT count(*) AS 'count' FROM variants WHERE INFO LIKE '%VARTYPE%' AND INFO LIKE '%PZScore%'"

    # prepare arguments for the query function
    args = argparse.Namespace(
        input=input_vcf,
        output=output_vcf,
        config=config,
        param=param,
        annotations=annotations,
        calculations=calculations,
        prioritizations=prioritizations,
        query=input_query,
        explode_infos=False,
        explode_infos_prefix="",
        explode_infos_fields="*",
        include_header=False,
        arguments_dict=arguments_dict,
    )

    # Remove if output file exists
    remove_if_exists([output_vcf])

    # Query
    process(args)

    # Create object
    variants = Variants(conn=None, input=output_vcf, config=config, load=True)

    # Check annotation
    result = variants.get_query_to_df("SELECT count FROM variants")
    log.debug(f"result={result}")
    assert result["count"][0] == 7
