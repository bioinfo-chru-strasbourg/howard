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
from test_needed import *



def test_convert():

    # Init files
    input_vcf = tests_data_folder + "/example.vcf.gz"
    output_vcf = "/tmp/output_file.tsv"
    config = {'threads': 4}

    #for explode_infos in [True, False]:
    for explode_infos in [True]:

        #for explode_infos_prefix in ["", "INFO/", "CUSTOM_"]:
        for explode_infos_prefix in [""]:

            #for explode_infos_fields in ['ALL', 'DP,SIFT,AA']:
            for explode_infos_fields in ['DP,SIFT,AA']:

                # prepare arguments for the query function
                args = argparse.Namespace(
                    input = input_vcf,
                    output = output_vcf,
                    config = config,
                    explode_infos = explode_infos,
                    explode_infos_prefix = explode_infos_prefix,
                    explode_infos_fields = explode_infos_fields,
                    include_header = True,
                    arguments_dict = arguments_dict
                )

                # Remove if output file exists
                remove_if_exists([output_vcf])

                # Query
                convert(args)

                # Check output file exists
                assert os.path.exists(output_vcf)

                # DEVEL
                variants = Variants(input=output_vcf, load=True)
                res = variants.get_query_to_df("SELECT * FROM variants")
                log.debug(res)

                # read the contents of the actual output file
                with open(output_vcf, 'r') as f:
                    result_output_nb_lines = 0
                    result_output_nb_variants = 0
                    for line in f:
                        result_output_nb_lines += 1
                        if not line.startswith("#"):
                            result_output_nb_variants += 1

                # Expected result
                expected_result_nb_lines = 60
                expected_result_nb_variants = 7

                # Compare
                assert result_output_nb_lines == expected_result_nb_lines
                assert result_output_nb_variants == expected_result_nb_variants
