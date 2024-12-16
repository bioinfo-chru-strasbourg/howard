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


def test_databases_download_exomiser():
    """
    The function `test_databases_download_exomiser` tests the `databases_download_exomiser` function by
    checking if the downloaded files match the expected number and if the specified assemblies are
    present.
    """

    # Init
    exomiser_data_hg19_source = os.path.join(
        tests_databases_folder, "exomiser", "test_hg19.zip"
    )
    exomiser_data_hg38_source = os.path.join(
        tests_databases_folder, "exomiser", "test_hg38.zip"
    )
    exomiser_phenotype_source = os.path.join(
        tests_databases_folder, "exomiser", "test_phenotype.zip"
    )
    threads = 2

    # Tmp folder
    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # Assembly
        assemblies = "hg19,hg38"
        assemblies_list = [value for value in assemblies.split(",")]

        # Exomiser folder
        exomiser_folder = tmp_dir

        # Exomiser release
        exomiser_release = "test"
        exomiser_phenotype_release = "test"

        # Download exomiser simulation
        exomiser_data_hg19_target = os.path.join(tmp_dir, "test_hg19.zip")
        shutil.copy(exomiser_data_hg19_source, exomiser_data_hg19_target)
        exomiser_data_hg38_target = os.path.join(tmp_dir, "test_hg38.zip")
        shutil.copy(exomiser_data_hg38_source, exomiser_data_hg38_target)
        exomiser_phenotype_target = os.path.join(tmp_dir, "test_phenotype.zip")
        shutil.copy(exomiser_phenotype_source, exomiser_phenotype_target)

        # Download and prepare database
        databases_download_exomiser(
            assemblies=assemblies_list,
            exomiser_folder=exomiser_folder,
            exomiser_release=exomiser_release,
            exomiser_phenotype_release=exomiser_phenotype_release,
            threads=threads,
        )

        # Check
        downloaded_files = os.listdir(exomiser_folder)
        for assembly in assemblies_list:
            assert assembly in downloaded_files
            downloaded_assembly_files = os.listdir(f"{exomiser_folder}/{assembly}")
            expected_files = [
                f"test_{assembly}",
                "application.properties",
                "test_phenotype",
            ]
            for expected_file in expected_files:
                if expected_file not in downloaded_assembly_files:
                    assert False
            assert True
