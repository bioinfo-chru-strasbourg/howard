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


def test_get_argument():
    # test case 1: argument exists in the dictionary
    arguments = {'arg1': {'type': 'string', 'help': 'argument 1'}}
    assert get_argument(arguments, 'arg1', required=True) == {'type': 'string', 'help': 'argument 1', 'required': True}

    # test case 2: argument does not exist in the dictionary
    arguments = {'arg1': {'type': 'string', 'help': 'argument 1'}}
    assert get_argument(arguments, 'arg2') == {}

    # test case 3: required is set to None
    arguments = {'arg1': {'type': 'string', 'help': 'argument 1'}}
    assert get_argument(arguments, 'arg1', required=None) == {'type': 'string', 'help': 'argument 1'}

    # test case 4: required is set to False
    arguments = {'arg1': {'type': 'string', 'help': 'argument 1'}}
    assert get_argument(arguments, 'arg1', required=False) == {'type': 'string', 'help': 'argument 1', 'required': False}
