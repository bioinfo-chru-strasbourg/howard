# -*- coding: utf-8 -*-
"""
Tests

Usage:
pytest tests/

Coverage:
coverage run -m pytest . -x -v
coverage report --include=howard/* -m
"""

# from howard.functions.commons import *
# from howard.tools.tools import *

import argparse

from howard.tools.tools import arguments_dict
from howard.tools.stats import stats

from test_needed import tests_data_folder


def test_stats():

    # Init files
    input_vcf = tests_data_folder + "/example.vcf.gz"
    config = {}

    # prepare arguments for the query function
    args = argparse.Namespace(
        input=input_vcf, config=config, arguments_dict=arguments_dict
    )

    # Query
    try:
        stats(args)
        assert True
    except:
        assert False
