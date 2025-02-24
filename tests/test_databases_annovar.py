# -*- coding: utf-8 -*-
"""
Tests

Usage:
pytest tests/

Coverage:
coverage run -m pytest tests/test_databases.py -x -vv --log-cli-level=DEBUG --capture=tee-sys
coverage report --include=howard/* -m
"""

import fnmatch
import os
from tempfile import TemporaryDirectory

from howard.functions.databases import databases_download_annovar
from test_needed import tests_folder


def test_databases_download_annovar_multiple_assembly():
    """
    The function `test_databases_download_annovar_multiple_assembly` tests the functionality of
    downloading multiple files with different assemblies using the `databases_download_annovar`
    function.
    """

    # Test downloading an existing file
    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # assembly
        assemblies = ["hg19", "hg38"]

        # files
        file_list = ["nci60"]

        # Threads
        threads = 2

        # Download
        databases_download_annovar(
            folder=tmp_dir, files=file_list, assemblies=assemblies, threads=threads
        )

        # Dowloaded files
        downloaded_files = os.listdir(tmp_dir)

        # Check
        for assembly in assemblies:
            assert assembly in downloaded_files
            downloaded_files_assembly = os.listdir(f"{tmp_dir}/{assembly}")
            assert f"{assembly}_refGene.txt" in downloaded_files_assembly
            for file in file_list:
                downloaded_file = f"{assembly}_{file}.txt"
                assert downloaded_file in downloaded_files_assembly

        # Download
        databases_download_annovar(
            folder=tmp_dir, files=file_list, assemblies=assemblies, threads=threads
        )

        # Dowloaded files
        downloaded_files_bis = os.listdir(tmp_dir)

        # Check
        for assembly in assemblies:
            assert assembly in downloaded_files_bis
            downloaded_files_bis_assembly = os.listdir(f"{tmp_dir}/{assembly}")
            assert f"{assembly}_refGene.txt" in downloaded_files_bis_assembly
            for file in file_list:
                downloaded_file = f"{assembly}_{file}.txt"
                assert downloaded_file in downloaded_files_bis_assembly


def test_databases_download_annovar_mandatory_refgene():
    """
    The function `test_databases_download_annovar_mandatory_refgene` tests the downloading of the
    mandatory file `refGene` from the ANNOVAR databases.
    """

    # Test downloading mandatory file refGene (no file list in input)
    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # assembly
        assemblies = ["hg19"]

        # files
        file_list = None

        # Threads
        threads = 2

        # Download
        databases_download_annovar(
            folder=tmp_dir, files=file_list, assemblies=assemblies, threads=threads
        )

        # Dowloaded files
        downloaded_files = os.listdir(tmp_dir)

        # Check
        for assembly in assemblies:
            assert assembly in downloaded_files
            downloaded_files_assembly = os.listdir(f"{tmp_dir}/{assembly}")
            assert f"{assembly}_refGene.txt" in downloaded_files_assembly


def test_databases_download_annovar_pattern_files():
    """
    The function `test_databases_download_annovar_pattern_files` tests the functionality of
    downloading multiple files with a pattern using the `databases_download_annovar` function.
    """

    # Test downloading multiple files with pattern
    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # assembly
        assemblies = ["hg19"]

        # files
        file_list = ["cosmic68*"]

        # Threads
        threads = 2

        # Download
        databases_download_annovar(
            folder=tmp_dir, files=file_list, assemblies=assemblies, threads=threads
        )

        # Dowloaded files
        downloaded_files = os.listdir(tmp_dir)

        # Check
        for assembly in assemblies:
            assert assembly in downloaded_files
            downloaded_files_assembly = os.listdir(f"{tmp_dir}/{assembly}")
            for file in file_list:
                downloaded_file = f"{assembly}_{file}.txt"
                filtered_files = fnmatch.filter(
                    downloaded_files_assembly, downloaded_file
                )
                assert len(filtered_files) > 1
