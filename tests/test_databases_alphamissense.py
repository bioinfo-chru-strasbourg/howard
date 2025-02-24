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
from tempfile import TemporaryDirectory

from howard.functions.databases import databases_download_alphamissense

from test_needed import tests_folder


def test_download_alphamissense():
    """
    The function `test_download_alphamissense` tests the `databases_download_alphamissense` function by
    downloading AlphaMissense databases for specified assemblies and checking if the files are
    downloaded correctly.
    """

    # Tmp folder
    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # Assembly
        assemblies = "hg19"
        assemblies_list = [value for value in assemblies.split(",")]

        # AlphaMissense folder
        alphamissense_folder = os.path.join(tmp_dir, "alphamissense")

        # Threads
        threads = 2

        # Download AlphaMissense
        databases_download_alphamissense(
            assemblies=assemblies_list,
            alphamissense_folder=alphamissense_folder,
            threads=threads,
        )

        # Check
        downloaded_files = os.listdir(alphamissense_folder)
        for assembly in assemblies_list:
            assert assembly in downloaded_files
            downloaded_assembly_files = os.listdir(f"{alphamissense_folder}/{assembly}")
            log.debug(downloaded_assembly_files)
            nb_files = 3
            assert len(downloaded_assembly_files) == nb_files

        # Download AlphaMissense again
        databases_download_alphamissense(
            assemblies=assemblies_list,
            alphamissense_folder=alphamissense_folder,
            threads=threads,
        )

        # Check
        downloaded_files = os.listdir(alphamissense_folder)
        for assembly in assemblies_list:
            assert assembly in downloaded_files
            downloaded_assembly_files = os.listdir(f"{alphamissense_folder}/{assembly}")
            log.debug(downloaded_assembly_files)
            nb_files = 3
            assert len(downloaded_assembly_files) == nb_files
