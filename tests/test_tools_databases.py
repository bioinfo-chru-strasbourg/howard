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
from howard.commons import *
from howard.tools.databases import *


tests_folder = os.path.dirname(__file__)


def test_databases_download_annovar():
    """
    This function tests the functionality of the databases_download_annovar function by downloading and
    checking various files.
    """

    # Test downloading an existing file
    with TemporaryDirectory() as tmp_dir:
       
        # assembly
        assemblies = ["hg19","hg38"]
        
        # files
        file_list = ['nci60']
        
        # Download
        databases_download_annovar(folder=tmp_dir, files=file_list, assemblies=assemblies)
        
        # Dowloaded files
        downloaded_files = os.listdir(tmp_dir)
        
        # Check
        for assembly in assemblies:
            assert  f"{assembly}_refGene.txt" in downloaded_files
            for file in file_list:
                downloaded_file = f"{assembly}_{file}.txt"
                assert downloaded_file in downloaded_files

        # Download
        databases_download_annovar(folder=tmp_dir, files=file_list, assemblies=assemblies)
        
        # Dowloaded files
        downloaded_files_bis = os.listdir(tmp_dir)

        assert len(downloaded_files) == len(downloaded_files_bis)
        
            
    # Test downloading mandatory file refGene (no file list in input)
    with TemporaryDirectory() as tmp_dir:
        
        # assembly
        assemblies = ["hg19"]
        
        # files
        file_list = None
        
        # Download
        databases_download_annovar(folder=tmp_dir, files=file_list, assemblies=assemblies)
        
        # Dowloaded files
        downloaded_files = os.listdir(tmp_dir)
        
        # Check
        for assembly in assemblies:
            downloaded_file = f"{assembly}_refGene.txt"
            assert downloaded_file in downloaded_files

    # Test downloading multiple files with pattern
    with TemporaryDirectory() as tmp_dir:
        
        # assembly
        assemblies = ["hg19"]
        
        # files
        file_list = ['cosmic68*']
        
        # Download
        databases_download_annovar(folder=tmp_dir, files=file_list, assemblies=assemblies)
        
        # Dowloaded files
        downloaded_files = os.listdir(tmp_dir)
        
        # Check
        for assembly in assemblies:
            for file in file_list:
                downloaded_file = f"{assembly}_{file}.txt"
                filtered_files = fnmatch.filter(downloaded_files, downloaded_file)
                assert len(filtered_files) > 1

    # Test downloading an existing file in multiple assemblies
    # with TemporaryDirectory() as tmp_dir:
        
    #     # assembly
    #     assemblies = ["hg19","hg38"]
        
    #     # files
    #     file_list = ['nci60']
        
    #     # Download
    #     databases_download_annovar(folder=tmp_dir, files=file_list, assemblies=assemblies)
        
    #     # Dowloaded files
    #     downloaded_files = os.listdir(tmp_dir)
        
    #     # Check
    #     for assembly in assemblies:
    #         for file in file_list:
    #             downloaded_file = f"{assembly}_{file}.txt"
    #             assert downloaded_file in downloaded_files


def test_databases_download():
    """
    This function tests the download of databases for Annovar and snpEff tools.
    """

    # Tmp folder
    with TemporaryDirectory() as tmp_dir:

        # assembly
        assemblies = ["hg19"]

        # files
        annovar_file_list = ['nci60']
        annovar_file_list_list = [value for val in annovar_file_list for value in val.split(',')]

        # config
        config = {"tools": {"snpeff": {"jar": default_snpeff_bin}}}

        # annovar URL
        download_annovar_url = default_annovar_url

        # Arguments
        args = argparse.Namespace(
            assembly=assemblies,
            download_annovar=tmp_dir,
            download_annovar_files=annovar_file_list,
            download_annovar_url=download_annovar_url,
            download_snpeff=tmp_dir,
            config=config
        )

        # Download
        databases_download(args)

        # Dowloaded files
        downloaded_files = os.listdir(tmp_dir)

        # Check Annovar
        # Check
        for assembly in assemblies:
            for file in annovar_file_list_list:
                downloaded_file = f"{assembly}_{file}.txt"
                assert downloaded_file in downloaded_files

        # Check snpEff databases list file
        snpeff_databases_list_file = "snpeff_databases.list"
        assert snpeff_databases_list_file in downloaded_files
        
        # Check assembly folders
        for assembly in assemblies:
            assert assembly in downloaded_files


# def test_databases_download_snpeff():
#     """
#     This function tests the download of snpEff databases for specified assemblies.
#     """

#     with TemporaryDirectory() as tmp_dir:
        
#         # assembly
#         assemblies = ["hg19","hg38"]
        
#         # config
#         config = {"tools": {"snpeff": {"jar": default_snpeff_bin}}}
        
#         # Dowloaded
#         databases_download_snpeff(folder=tmp_dir, assemblies=assemblies, config=config)
        
#         # Dowloaded files
#         downloaded_files = os.listdir(tmp_dir)
        
#         # Check
        
#         # Check snpEff databases list file
#         snpeff_databases_list_file = "snpeff_databases.list"
#         assert snpeff_databases_list_file in downloaded_files
        
#         # Check assembly folders
#         for assembly in assemblies:
#             assert assembly in downloaded_files

