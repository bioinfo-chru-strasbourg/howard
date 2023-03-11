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

from howard.objects.variants import Variants
from howard.commons import *


tests_folder = os.path.dirname(__file__)


def test_set_log_level():

    # define verbosity
    verbosity = "info"

    # set verbosity
    result_verbosity = set_log_level(verbosity)

    assert verbosity == result_verbosity

def test_set_log_level_error():

    # define verbosity
    verbosity = "not_a_level"

    # check verbosity
    with pytest.raises(ValueError) as e:
        set_log_level(verbosity)
    assert str(e.value) == "Unknown verbosity level:" + verbosity


def test_split_interval_either():

    start = 0
    end = 1000
    step = None
    ncuts = None
    with pytest.raises(ValueError) as e:
        split_interval(start, end, step=step, ncuts=ncuts)
    assert str(e.value) == "Either step or ncuts must be provided"


def test_split_interval_only():

    start = 0
    end = 1000
    step = 100
    ncuts = 4
    with pytest.raises(ValueError) as e:
        split_interval(start, end, step=step, ncuts=ncuts)
    assert str(e.value) == "Only one of step or ncuts must be provided"


def test_split_interval_step():

    start = 0
    end = 1000
    step = 100
    ncuts = None
    split = split_interval(start, end, step=step, ncuts=ncuts)
    split_expectd = [0, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000]
    
    assert split == split_expectd


def test_split_interval_ncuts():

    start = 0
    end = 1000
    step = None
    ncuts = 4
    split = split_interval(start, end, step=step, ncuts=ncuts)
    split_expectd = [0, 250, 500, 750, 1000]
    
    assert split == split_expectd


def test_merged_regions():

    # Define a list of genomic regions
    regions = [
        ('chr1', 100, 200),
        ('chr1', 150, 300),
        ('chr2', 500, 600),
        ('chr2', 550, 650),
        ('chr3', 800, 900)
    ]
    # Call the merge_regions function
    merged_regions = merge_regions(regions)
    # Define merged regions expected
    merged_regions_expected = [
        ('chr1', 100, 300),
        ('chr2', 500, 650),
        ('chr3', 800, 900)
    ]
    
    assert merged_regions == merged_regions_expected


def test_create_where_clause():

    # Define a list of merged regions
    merged_regions = [
        ('chr1', 100, 300),
        ('chr2', 500, 650),
        ('chr3', 800, 900),
        ('chr3', 1000, 1200)
    ]

    # define table
    table = "variants"

    # Call the create_where_clause function
    where_clause = create_where_clause(merged_regions, table=table)

    # Create expected where clause
    where_clause_expected = """ ( variants."#CHROM" = 'chr1' AND (   (variants.POS >= 100 AND variants.POS <= 300)  ) )   OR  ( variants."#CHROM" = 'chr2' AND (   (variants.POS >= 500 AND variants.POS <= 650)  ) )   OR  ( variants."#CHROM" = 'chr3' AND (   (variants.POS >= 800 AND variants.POS <= 900)   OR  (variants.POS >= 1000 AND variants.POS <= 1200)  ) ) """
    
    assert where_clause.strip() == where_clause_expected.strip()


def test_command():

    # Creation of a command
    cmd = "echo 'Command'"

    # Execute command
    command_output = command(cmd)

    # Expected result
    results_expected = "Command"

    assert command_output == results_expected


def test_run_parallel_commands():

    # Creation of a list of command
    commands = [
        "echo 'Command 1'",
        "echo 'Command 2'",
        "echo 'Command 3'"
    ]

    # Execute commands in parallele
    results = run_parallel_commands(commands, 2)

    # Expected results
    results_expected = ['Command 1', 'Command 2', 'Command 3']

    assert results == results_expected


def test_run_parallel_functions():

    # List of functions (examples)
    functions = [example_function(1, "hello"), example_function(2, "world")]
    threads = 2

    # Launch functions
    results = run_parallel_functions(functions, threads)

    # Number of result expected
    expected_output_length = 2

    assert len(results) == expected_output_length


def test_example_function():

    # Launch functions
    result = example_function(1, "hello")
                               
    # result expected
    expected_result = [1, "hello"]

    assert result == expected_result

