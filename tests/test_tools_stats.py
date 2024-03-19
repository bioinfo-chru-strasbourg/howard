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
