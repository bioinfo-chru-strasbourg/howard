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



def test_remove_if_exists():

    # filename
    filename = "/tmp/output.test"

    # create filename
    fhandle = open(filename, 'a')
    try:
        os.utime(filename, None)
    finally:
        fhandle.close()
    created = os.path.exists(filename)

    # remove filename
    remove_if_exists(filename)

    # check delete
    deleted = not os.path.exists(filename)

    assert created and deleted



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


def test_get_index():

    # Test with a list of values
    values = ['a', 'b', 'c', 'd']
    assert get_index('a', values) == 0
    assert get_index('d', values) == 3
    assert get_index('e', values) == -1

    # Test with an empty list
    assert get_index('a', []) == -1
    assert get_index('', []) == -1

    # Test with a None element
    assert get_index(None, values) == -1


# def test_find_nomen():
#     hgvs = "NM_001234.5:c.123A>T,NM_001234.6:g.456C>A"
#     pattern = "GNOMEN:TNOMEN:ENOMEN:CNOMEN:RNOMEN:NNOMEN:PNOMEN"
#     transcripts = ["NM_001234.5", "NM_001234.6"]
#     expected_output = {
#         "NOMEN": "A:123T",
#         "CNOMEN": "c.123A>T",
#         "RNOMEN": None,
#         "NNOMEN": None,
#         "PNOMEN": None,
#         "TVNOMEN": "NM_001234.6",
#         "TNOMEN": "NM_001234",
#         "VNOMEN": "5",
#         "ENOMEN": None,
#         "GNOMEN": "g.456C>A"
#     }
#     assert find_nomen(hgvs, pattern, transcripts) == expected_output


def test_find_nomen_full():
    # Test case 1
    print("test1")
    hgvs = "NM_001637.3:c.1582G>T"
    transcripts = ["NM_001637.3"]
    expected_output = {
        "NOMEN": "NM_001637:c.1582G>T",
        "CNOMEN": "c.1582G>T",
        "RNOMEN": None,
        "NNOMEN": None,
        "PNOMEN": None,
        "TVNOMEN": "NM_001637.3",
        "TNOMEN": "NM_001637",
        "VNOMEN": "3",
        "ENOMEN": None,
        "GNOMEN": None,
    }
    assert find_nomen(hgvs, transcripts=transcripts) == expected_output

    # Test case 2
    print("test2")
    hgvs = "NM_001637.3:c.1582G>T,NM_001637.3:c.1583G>T"
    transcripts = ["NM_001637.3"]
    expected_output = {
        "NOMEN": "NM_001637:c.1582G>T",
        "CNOMEN": "c.1582G>T",
        "RNOMEN": None,
        "NNOMEN": None,
        "PNOMEN": None,
        "TVNOMEN": "NM_001637.3",
        "TNOMEN": "NM_001637",
        "VNOMEN": "3",
        "ENOMEN": None,
        "GNOMEN": None,
    }
    assert find_nomen(hgvs, transcripts=transcripts) == expected_output

    # Test case 3
    print("test3")
    hgvs = "NM_001637.3:c.1582G>T,NM_001637.3:c.1583G>T,NM_001637.2:c.1582G>T:p.G12D"
    transcripts = ["NM_001637.2", "NM_001637.3"]
    expected_output = {
        "NOMEN": "NM_001637:c.1582G>T:p.G12D",
        "CNOMEN": "c.1582G>T",
        "RNOMEN": None,
        "NNOMEN": None,
        "PNOMEN": "p.G12D",
        "TVNOMEN": "NM_001637.2",
        "TNOMEN": "NM_001637",
        "VNOMEN": "2",
        "ENOMEN": None,
        "GNOMEN": None,
    }
    assert find_nomen(hgvs, transcripts=transcripts) == expected_output

    # Test case 4
    print("test3")
    hgvs = "Gene1:exon12:n.1582G>T:NR_001637.3"
    transcripts = []
    expected_output = {
        "NOMEN": "Gene1:NR_001637:exon12:n.1582G>T",
        "CNOMEN": None,
        "RNOMEN": None,
        "NNOMEN": "n.1582G>T",
        "PNOMEN": None,
        "TVNOMEN": "NR_001637.3",
        "TNOMEN": "NR_001637",
        "VNOMEN": "3",
        "ENOMEN": "exon12",
        "GNOMEN": "Gene1",
    }
    assert find_nomen(hgvs, transcripts=transcripts) == expected_output

    # Test case 5
    print("test3")
    hgvs = "Gene1:exon12:r.1582G>T"
    transcripts = []
    expected_output = {
        "NOMEN": "Gene1:exon12:r.1582G>T",
        "CNOMEN": None,
        "RNOMEN": "r.1582G>T",
        "NNOMEN": None,
        "PNOMEN": None,
        "TVNOMEN": None,
        "TNOMEN": None,
        "VNOMEN": None,
        "ENOMEN": "exon12",
        "GNOMEN": "Gene1",
    }
    assert find_nomen(hgvs, transcripts=transcripts) == expected_output


def test_get_gzip():

    command_gzip_expected = "bgzip -c  --threads=2 --compress-level=5"

    command_gzip = get_bgzip(threads=2, level=5)
    print(command_gzip)

    assert command_gzip.strip() == command_gzip_expected.strip()
          


def test_find_all():
    # Create a temporary directory structure
    with tempfile.TemporaryDirectory() as tmpdir:
        # Create some files with the name 'test_file' in different directories
        open(os.path.join(tmpdir, 'test_file'), 'a').close()
        os.makedirs(os.path.join(tmpdir, 'subdir'))
        open(os.path.join(tmpdir, 'subdir', 'test_file'), 'a').close()

        # Test that find_all returns the correct paths
        assert find_all('test_file', tmpdir) == [
            os.path.join(tmpdir, 'test_file'),
            os.path.join(tmpdir, 'subdir', 'test_file')
        ]


def test_find_genome():
    # create a temporary directory
    with tempfile.TemporaryDirectory() as tmpdir:
        # specify a non-existent path for the genome file
        genome_path = os.path.join(tmpdir, 'nonexistent_genome.fa')
        # call the function to find the genome file
        genome_path = find_genome(genome_path)
        # check if the genome file was found
        assert os.path.exists(genome_path)

