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


def test_annotation():

    # Init files
    input_vcf = tests_folder + "/data/example.vcf.gz"
    output_vcf = "/tmp/output_file.tsv"
    config = {}
    annotations = tests_folder + "/data/annotations/nci60.parquet"

    # prepare arguments for the query function
    args = argparse.Namespace(
        input = input_vcf,
        output = output_vcf,
        config = config,
        annotations = annotations,
    )

    # Remove if output file exists
    remove_if_exists([output_vcf])

    # Query
    annotation(args)

    # Check output file exists
    assert os.path.exists(output_vcf)

    # read the contents of the actual output file
    with open(output_vcf, 'r') as f:
        result_output_nb_lines = 0
        result_output_nb_variants = 0
        for line in f:
            result_output_nb_lines += 1
            if not line.startswith("#"):
                result_output_nb_variants += 1

    # Expected result
    expected_result_nb_lines = 61
    expected_result_nb_variants = 7

    # Compare
    assert result_output_nb_lines == expected_result_nb_lines
    assert result_output_nb_variants == expected_result_nb_variants

