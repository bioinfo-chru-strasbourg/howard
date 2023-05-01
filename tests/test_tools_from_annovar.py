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


def test_from_annovar():

    # Init files
    input_vcf = tests_folder + "/data/annotations/hg19_nci60.txt"
    output_vcf = "/tmp/hg19_nci60.vcf.gz"
    output_parquet = "/tmp/hg19_nci60.parquet"
    genome = "/databases/genomes/current/hg19.fa"
    annovar_code = "nci60"
    config = {}

    # remove
    remove_if_exists([output_vcf])
    remove_if_exists([output_vcf+".hdr"])
    remove_if_exists([output_parquet])
    remove_if_exists([output_parquet+".hdr"])

    # Find genome
    genome = find_genome(genome)

    # prepare arguments for the query function
    args = argparse.Namespace(
        input = input_vcf,
        output = output_vcf,
        genome = genome,
        annovar_code = annovar_code,
        to_parquet = output_parquet,
        threads = "2",
        reduce_memory = "disable",
        multi_variant = "disable",
        config = config
    )

    # Process
    try:
        from_annovar(args)
        assert True
    except:
        assert False

    # Files exists
    assert os.path.exists(output_vcf)
    assert os.path.exists(output_vcf+".hdr")
    assert os.path.exists(output_parquet)
    assert os.path.exists(output_parquet+".hdr")


def test_from_annovar_reduce_memory():

    # Init files
    input_vcf = tests_folder + "/data/annotations/hg19_nci60.txt"
    output_vcf = "/tmp/hg19_nci60.vcf.gz"
    output_parquet = "/tmp/hg19_nci60.parquet"
    genome = "/databases/genomes/current/hg19.fa"
    annovar_code = "nci60"
    config = {}

    # remove
    remove_if_exists([output_vcf])
    remove_if_exists([output_vcf+".hdr"])
    remove_if_exists([output_parquet])
    remove_if_exists([output_parquet+".hdr"])

    # Find genome
    genome = find_genome(genome)

    # prepare arguments for the query function
    args = argparse.Namespace(
        input = input_vcf,
        output = output_vcf,
        genome = genome,
        annovar_code = annovar_code,
        to_parquet = output_parquet,
        threads = "2",
        reduce_memory = "enable",
        multi_variant = "disable",
        config = config
    )

    # Process
    try:
        from_annovar(args)
        assert True
    except:
        assert False

    # Files exists
    assert os.path.exists(output_vcf)
    assert os.path.exists(output_vcf+".hdr")
    assert os.path.exists(output_parquet)
    assert os.path.exists(output_parquet+".hdr")


def test_from_annovar_multi_variant():

    # Init files
    input_vcf = tests_folder + "/data/annotations/hg19_nci60.txt"
    output_vcf = "/tmp/hg19_nci60.vcf.gz"
    output_parquet = "/tmp/hg19_nci60.parquet"
    genome = "/databases/genomes/current/hg19.fa"
    annovar_code = "nci60"
    config = {}

    # remove
    remove_if_exists([output_vcf])
    remove_if_exists([output_vcf+".hdr"])
    remove_if_exists([output_parquet])
    remove_if_exists([output_parquet+".hdr"])

    # Find genome
    genome = find_genome(genome)

    # prepare arguments for the query function
    args = argparse.Namespace(
        input = input_vcf,
        output = output_vcf,
        genome = genome,
        annovar_code = annovar_code,
        to_parquet = output_parquet,
        threads = "2",
        reduce_memory = "disable",
        multi_variant = "enable",
        config = config
    )

    # Process
    try:
        from_annovar(args)
        assert True
    except:
        assert False

    # Files exists
    assert os.path.exists(output_vcf)
    assert os.path.exists(output_vcf+".hdr")
    assert os.path.exists(output_parquet)
    assert os.path.exists(output_parquet+".hdr")


def test_from_annovar_reduce_memory_multi_variant():

    # Init files
    input_vcf = tests_folder + "/data/annotations/hg19_nci60.txt"
    output_vcf = "/tmp/hg19_nci60.vcf.gz"
    output_parquet = "/tmp/hg19_nci60.parquet"
    genome = "/databases/genomes/current/hg19.fa"
    annovar_code = "nci60"
    config = {}

    # remove
    remove_if_exists([output_vcf])
    remove_if_exists([output_vcf+".hdr"])
    remove_if_exists([output_parquet])
    remove_if_exists([output_parquet+".hdr"])

    # Find genome
    genome = find_genome(genome)

    # prepare arguments for the query function
    args = argparse.Namespace(
        input = input_vcf,
        output = output_vcf,
        genome = genome,
        annovar_code = annovar_code,
        to_parquet = output_parquet,
        threads = "2",
        reduce_memory = "enable",
        multi_variant = "enable",
        config = config
    )

    # Process
    try:
        from_annovar(args)
        assert True
    except:
        assert False

    # Files exists
    assert os.path.exists(output_vcf)
    assert os.path.exists(output_vcf+".hdr")
    assert os.path.exists(output_parquet)
    assert os.path.exists(output_parquet+".hdr")
