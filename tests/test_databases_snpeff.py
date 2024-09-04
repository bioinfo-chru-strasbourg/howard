# -*- coding: utf-8 -*-
"""
Tests

Usage:
pytest tests/

Coverage:
coverage run -m pytest tests/test_databases.py -x -vv --log-cli-level=DEBUG --capture=tee-sys
coverage report --include=howard/* -m
"""

import logging as log
import os
import sys
from tempfile import TemporaryDirectory
import duckdb
import re
import Bio.bgzf as bgzf
import gzip
import pytest
import pandas as pd
from pandas.testing import assert_frame_equal
from unittest.mock import patch

from howard.objects.variants import Variants
from howard.objects.database import Database
from howard.functions.commons import *
from howard.tools.databases import *
from howard.tools.tools import arguments_dict

from test_needed import *


def test_databases_download_snpeff():
    """
    The function `test_databases_download_snpeff` downloads and prepares the snpEff database for specified
    assemblies.
    """

    # Init
    threads = 2

    # Full database generation, hg19 only (due to lack of hg38 assembly in tests)
    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # Assembly
        assemblies = "hg19"
        assemblies_list = [value for value in assemblies.split(",")]

        # snpEff folder
        snpeff_folder = tmp_dir

        # Create empty file
        open(os.path.join(snpeff_folder, "snpeff_databases.list"), "w").close()

        # Download and prepare database
        databases_download_snpeff(
            folder=snpeff_folder,
            assemblies=assemblies_list,
            config=tests_config,
            threads=threads,
        )

        # Check
        downloaded_files = os.listdir(snpeff_folder)
        for assembly in assemblies_list:
            assert assembly in downloaded_files
            downloaded_assembly_files = os.listdir(f"{snpeff_folder}/{assembly}")
            expected_files = ["sequence.bin"]
            for expected_file in expected_files:
                if expected_file not in downloaded_assembly_files:
                    assert False
            assert True


def test_databases_download_snpeff_mouse():
    """
    The function `test_databases_download_snpeff` downloads and prepares the snpEff database for specified
    assemblies.
    """

    # Init
    threads = 2

    # Full database generation, hg19 only (due to lack of hg38 assembly in tests)
    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # Assembly
        assemblies = "GRCm39.105"
        assemblies_list = [value for value in assemblies.split(",")]

        # snpEff folder
        snpeff_folder = tmp_dir

        # Create empty file
        open(os.path.join(snpeff_folder, "snpeff_databases.list"), "w").close()

        # Download and prepare database
        databases_download_snpeff(
            folder=snpeff_folder,
            assemblies=assemblies_list,
            config=tests_config,
            threads=threads,
        )

        # Check
        downloaded_files = os.listdir(snpeff_folder)
        for assembly in assemblies_list:
            assert assembly in downloaded_files
            downloaded_assembly_files = os.listdir(f"{snpeff_folder}/{assembly}")
            expected_files = ["sequence.bin"]
            for expected_file in expected_files:
                if expected_file not in downloaded_assembly_files:
                    assert False
            assert True
