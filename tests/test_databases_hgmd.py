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


def test_download_hgmd():
    """
    The function `test_download_hgmd` downloads HGMD files for specified assemblies and checks if the
    downloaded files match the expected files.
    """

    # Init
    threads = 2

    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # Assembly
        assemblies = "hg19"
        assemblies_list = [value for value in assemblies.split(",")]

        # Exomiser folder
        hgmd_folder = tmp_dir

        # HGMD conversion
        hgmd_file_hg19 = os.path.join(
            tests_databases_folder, "hgmd", "HGMD_TEST_hg19.vcf.gz"
        )
        databases_download_hgmd(
            assemblies=assemblies_list,
            hgmd_file=hgmd_file_hg19,
            hgmd_folder=hgmd_folder,
            threads=threads,
        )

        # Check
        downloaded_files = os.listdir(hgmd_folder)
        for assembly in assemblies_list:
            assert assembly in downloaded_files
            downloaded_assembly_files = os.listdir(f"{hgmd_folder}/{assembly}")
            expected_files = [
                "HGMD_TEST.vcf.gz.tbi",
                "HGMD_TEST.tsv.hdr",
                "HGMD_TEST.parquet.hdr",
                "HGMD_TEST.vcf.gz",
                "HGMD_TEST.parquet",
                "HGMD_TEST.tsv",
            ]
            for expected_file in expected_files:
                if expected_file not in downloaded_assembly_files:
                    assert False
            assert True
