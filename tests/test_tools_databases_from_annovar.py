"""
Tests

Usage:
pytest tests/

Coverage:
coverage run -m pytest tests/test_tools_databases_from_annovar.py -x -v --log-cli-level=DEBUG --capture=tee-sys
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
from howard.functions.from_annovar import *
from howard.tools.tools import *
from test_needed import *


tests_folder = os.path.dirname(__file__)


def test_from_annovar():

    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # Init files
        input_vcf = tests_databases_folder + "/others/hg19_nci60.txt"
        output_vcf = os.path.join(tmp_dir,"hg19_nci60.vcf.gz")
        output_parquet = os.path.join(tmp_dir,"hg19_nci60.parquet")
        genomes_folder = tests_config["folders"]["databases"]["genomes"]
        genome = genomes_folder + "/hg19//hg19.fa"
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
            input_annovar = input_vcf,
            output_annovar = output_vcf,
            genome = genome,
            annovar_code = annovar_code,
            annovar_to_parquet = output_parquet,
            threads = "2",
            annovar_reduce_memory = "disable",
            annovar_multi_variant = "disable",
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

    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # Init files
        input_vcf = tests_databases_folder + "/others/hg19_nci60.txt"
        output_vcf = os.path.join(tmp_dir,"hg19_nci60.vcf.gz")
        output_parquet = os.path.join(tmp_dir,"hg19_nci60.parquet")
        genomes_folder = tests_config["folders"]["databases"]["genomes"]
        genome = genomes_folder + "/hg19//hg19.fa"
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
            input_annovar = input_vcf,
            output_annovar = output_vcf,
            genome = genome,
            annovar_code = annovar_code,
            annovar_to_parquet = output_parquet,
            threads = "2",
            annovar_reduce_memory = "enable",
            annovar_multi_variant = "disable",
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

    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # Init files
        input_vcf = tests_databases_folder + "/others/hg19_nci60.txt"
        output_vcf = os.path.join(tmp_dir,"hg19_nci60.vcf.gz")
        output_parquet = os.path.join(tmp_dir,"hg19_nci60.parquet")
        genomes_folder = tests_config["folders"]["databases"]["genomes"]
        genome = genomes_folder + "/hg19//hg19.fa"
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
            input_annovar = input_vcf,
            output_annovar = output_vcf,
            genome = genome,
            annovar_code = annovar_code,
            annovar_to_parquet = output_parquet,
            threads = "2",
            annovar_reduce_memory = "disable",
            annovar_multi_variant = "enable",
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


    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # Init files
        input_vcf = tests_databases_folder + "/others/hg19_nci60.txt"
        output_vcf = os.path.join(tmp_dir,"hg19_nci60.vcf.gz")
        output_parquet = os.path.join(tmp_dir,"hg19_nci60.parquet")
        genomes_folder = tests_config["folders"]["databases"]["genomes"]
        genome = genomes_folder + "/hg19//hg19.fa"
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
            input_annovar = input_vcf,
            output_annovar = output_vcf,
            genome = genome,
            annovar_code = annovar_code,
            annovar_to_parquet = output_parquet,
            threads = "2",
            annovar_reduce_memory = "enable",
            annovar_multi_variant = "enable",
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

