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
from howard.commons import *
from howard.tools.tools import *


tests_folder = os.path.dirname(__file__)


def test_query():

    # Init files
    input_vcf = tests_data_folder + "/example.vcf.gz"
    output_vcf = "/tmp/output_file.tsv"
    config = {'threads': 4}
    input_query = "SELECT count(*) AS '#count' FROM variants"

    for explode_infos in [True, False]:

        # prepare arguments for the query function
        args = argparse.Namespace(
            input = input_vcf,
            output = output_vcf,
            config = config,
            query = input_query,
            explode_infos = explode_infos
        )

        # Remove if output file exists
        remove_if_exists([output_vcf])

        # Query
        query(args)

        # read the contents of the actual output file
        with open(output_vcf, 'r') as f:
            result_output_nb_lines = 0
            result_output_nb_variants = 0
            result_lines = []
            for line in f:
                result_output_nb_lines += 1
                if not line.startswith("#"):
                    result_output_nb_variants += 1
                    result_lines.append(line.strip())

        # Expected result
        expected_result_nb_lines = 54
        expected_result_nb_variants = 1
        expected_result_lines = ["7"]

        # Compare
        assert result_output_nb_lines == expected_result_nb_lines
        assert result_output_nb_variants == expected_result_nb_variants
        assert result_lines == expected_result_lines