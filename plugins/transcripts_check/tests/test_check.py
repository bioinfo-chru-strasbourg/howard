# -*- coding: utf-8 -*-
"""
Tests

Usage:
pytest tests/

Coverage:
coverage run -m pytest tests/test_commons.py -x -vv --log-cli-level=INFO --capture=tee-sys
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
import pandas as pd
from pandas.testing import assert_frame_equal
from unittest.mock import patch

from howard.objects.variants import Variants
from howard.functions.commons import *


from plugins.transcripts_check.__main__ import *

from howard.tools.tools import *


def test_check():
    """
    The function `test_detect_column_type` contains test cases for detecting the type of data in a
    column using the `detect_column_type` function.
    """

    tests_folder = os.path.dirname(__file__)
    tests_data_folder = os.path.join(tests_folder, "data")

    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # Init files
        input_vcf = tests_data_folder + "/example.ann.transcripts.vcf.gz"
        transcripts_expected = tests_data_folder + "/transcripts.tsv"
        transcripts_missing = None
        stats_json = os.path.join(tmp_dir, "stats.json")
        config = {}
        param = tests_data_folder + "/param.transcripts.json"

        # Expected results
        stats_json_expected = {
            "available": 24,
            "list": 6,
            "intersection": 4,
            "union": 26,
            "percent": 0.6666666666666666,
            "missing": 2,
            "missing_list": ["NM_005228.5", "NM_0123456"],
        }

        # prepare arguments for the query function
        args = argparse.Namespace(
            input=input_vcf,
            config=config,
            param=param,
            arguments_dict=arguments_dict,
            transcripts_expected=transcripts_expected,
            transcripts_missing=transcripts_missing,
            stats_json=stats_json,
        )

        # Query
        try:
            main(args)
            assert True
        except:
            assert False

        if os.path.isfile(stats_json):
            with open(stats_json, "r") as file:
                stats_results = yaml.safe_load(file)
            log.debug(f"stats_results={stats_results}")
            log.debug(f"stats_json_expected={stats_json_expected}")
            assert stats_results == stats_json_expected
        else:
            assert False


def test_check_with_version_included():
    """
    The function `test_detect_column_type` contains test cases for detecting the type of data in a
    column using the `detect_column_type` function.
    """

    tests_folder = os.path.dirname(__file__)
    tests_data_folder = os.path.join(tests_folder, "data")

    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # Init files
        input_vcf = tests_data_folder + "/example.ann.transcripts.vcf.gz"
        transcripts_expected = tests_data_folder + "/transcripts.tsv"
        transcripts_missing = None
        stats_json = os.path.join(tmp_dir, "stats.json")
        config = {}
        param = tests_data_folder + "/param.transcripts_version_not_removed.json"

        # Expected results
        stats_json_expected = {
            "available": 24,
            "list": 6,
            "intersection": 1,
            "union": 29,
            "percent": 0.16666666666666666,
            "missing": 5,
            "missing_list": [
                "NR_024540",
                "NR_036266",
                "NM_001346897",
                "NM_005228",
                "NM_0123456",
            ],
        }

        # prepare arguments for the query function
        args = argparse.Namespace(
            input=input_vcf,
            config=config,
            param=param,
            arguments_dict=arguments_dict,
            transcripts_expected=transcripts_expected,
            transcripts_missing=transcripts_missing,
            stats_json=stats_json,
        )

        # Query
        try:
            main(args)
            assert True
        except:
            assert False

        if os.path.isfile(stats_json):
            with open(stats_json, "r") as file:
                stats_results = yaml.safe_load(file)
            log.debug(f"stats_results={stats_results}")
            log.debug(f"stats_json_expected={stats_json_expected}")
            assert stats_results == stats_json_expected
        else:
            assert False
