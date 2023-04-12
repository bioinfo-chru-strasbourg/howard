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


def test_convert():

    # Init files
    input_vcf = tests_folder + "/data/example.vcf.gz"
    output_vcf = "/tmp/output_file.tsv"
    config = {'threads': 4}

    for export_infos in [True, False]:

        for export_infos_prefix in [None, "INFO/", "CUSTOM_"]:

            # prepare arguments for the query function
            args = argparse.Namespace(
                input = input_vcf,
                output = output_vcf,
                config = config,
                export_infos = export_infos,
                export_infos_prefix = export_infos_prefix
            )

            # Remove if output file exists
            remove_if_exists([output_vcf])

            # Query
            convert(args)

            # read the contents of the actual output file
            with open(output_vcf, 'r') as f:
                result_output_nb_lines = len(f.readlines())

            # Expected result
            expected_result_nb_lines = 8

            # Compare
            assert result_output_nb_lines == expected_result_nb_lines
