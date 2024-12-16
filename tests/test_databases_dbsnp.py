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


def test_databases_download_dbsnp_full():
    """
    The function `test_databases_download_dbsnp_full` downloads and prepares the dbsnp database for specified
    assemblies.
    """

    # Genomes
    genomes_folder = tests_config["folders"]["databases"]["genomes"]
    download_needed_databases()
    threads = 2

    # Full database generation, hg19 only (due to lack of hg38 assembly in tests)
    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # Assembly
        assemblies = "hg19"
        assemblies_list = [value for value in assemblies.split(",")]

        # Releases
        dbsnp_releases = ["b156"]

        # Exomiser folder
        dbsnp_folder = tmp_dir

        # Download exomiser simulation
        dnsnp_assemblies_map: dict = {"hg19": "25", "hg38": "40"}
        for assembly in assemblies_list:
            for dbsnp_release in dbsnp_releases:
                dbsnp_data_source = os.path.join(
                    tests_databases_folder,
                    "dbsnp",
                    f"GCF_000001405.{dnsnp_assemblies_map.get(assembly)}.gz",
                )
                dbsnp_data_target = os.path.join(
                    dbsnp_folder,
                    assembly,
                    dbsnp_release,
                    f"GCF_000001405.{dnsnp_assemblies_map.get(assembly)}.gz",
                )
                if not os.path.exists(
                    os.path.join(dbsnp_folder, assembly, dbsnp_release)
                ):
                    Path(os.path.join(dbsnp_folder, assembly, dbsnp_release)).mkdir(
                        parents=True, exist_ok=True
                    )
                if not os.path.exists(
                    os.path.join(dbsnp_folder, assembly, dbsnp_data_target)
                ):
                    shutil.copy(dbsnp_data_source, dbsnp_data_target)

        # Download and prepare database
        databases_download_dbsnp(
            assemblies=assemblies_list,
            dbsnp_folder=dbsnp_folder,
            dbsnp_vcf=True,
            dbsnp_parquet=True,
            genomes_folder=genomes_folder,
            threads=threads,
        )

        # Check
        downloaded_files = os.listdir(dbsnp_folder)
        for assembly in assemblies_list:
            assert assembly in downloaded_files
            downloaded_assembly_files = os.listdir(f"{dbsnp_folder}/{assembly}")
            assert "default" in downloaded_assembly_files
            for dbsnp_release in dbsnp_releases:
                assert dbsnp_release in downloaded_assembly_files
                downloaded_assembly_release_files = os.listdir(
                    f"{dbsnp_folder}/{assembly}/{dbsnp_release}"
                )
                expected_files = [
                    f"GCF_000001405.{dnsnp_assemblies_map.get(assembly)}.gz",
                    "dbsnp.vcf.gz",
                    "dbsnp.parquet.hdr",
                    "dbsnp.parquet",
                ]
                for expected_file in expected_files:
                    if expected_file not in downloaded_assembly_release_files:
                        assert False
                assert True

        # Download and prepare database again
        databases_download_dbsnp(
            assemblies=assemblies_list,
            dbsnp_folder=dbsnp_folder,
            dbsnp_vcf=True,
            dbsnp_parquet=True,
            genomes_folder=genomes_folder,
            threads=threads,
        )

        # Check
        downloaded_files = os.listdir(dbsnp_folder)
        for assembly in assemblies_list:
            assert assembly in downloaded_files
            downloaded_assembly_files = os.listdir(f"{dbsnp_folder}/{assembly}")
            assert "default" in downloaded_assembly_files
            for dbsnp_release in dbsnp_releases:
                assert dbsnp_release in downloaded_assembly_files
                downloaded_assembly_release_files = os.listdir(
                    f"{dbsnp_folder}/{assembly}/{dbsnp_release}"
                )
                expected_files = [
                    f"GCF_000001405.{dnsnp_assemblies_map.get(assembly)}.gz",
                    "dbsnp.vcf.gz",
                    "dbsnp.parquet.hdr",
                    "dbsnp.parquet",
                ]
                for expected_file in expected_files:
                    if expected_file not in downloaded_assembly_release_files:
                        assert False
                assert True


def test_databases_download_dbsnp_multi_assemblies():
    """
    The function `test_databases_download_dbsnp_multi_assemblies` downloads and prepares the dbsnp database for specified
    assemblies.
    """

    # Genomes
    genomes_folder = tests_config["folders"]["databases"]["genomes"]
    download_needed_databases()
    threads = 2

    # Multi assembly (without VCF and PArquet generation, due to lack of assembly hg38 in tests)
    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # Assembly
        assemblies = "hg19,hg38"
        assemblies_list = [value for value in assemblies.split(",")]

        # Releases
        dbsnp_releases = ["b156"]

        # Exomiser folder
        dbsnp_folder = tmp_dir

        # Download exomiser simulation
        dnsnp_assemblies_map: dict = {"hg19": "25", "hg38": "40"}
        for assembly in assemblies_list:
            for dbsnp_release in dbsnp_releases:
                dbsnp_data_source = os.path.join(
                    tests_databases_folder,
                    "dbsnp",
                    f"GCF_000001405.{dnsnp_assemblies_map.get(assembly)}.gz",
                )
                dbsnp_data_target = os.path.join(
                    dbsnp_folder,
                    assembly,
                    dbsnp_release,
                    f"GCF_000001405.{dnsnp_assemblies_map.get(assembly)}.gz",
                )
                if not os.path.exists(
                    os.path.join(dbsnp_folder, assembly, dbsnp_release)
                ):
                    Path(os.path.join(dbsnp_folder, assembly, dbsnp_release)).mkdir(
                        parents=True, exist_ok=True
                    )
                if not os.path.exists(
                    os.path.join(dbsnp_folder, assembly, dbsnp_data_target)
                ):
                    shutil.copy(dbsnp_data_source, dbsnp_data_target)

        # Download and prepare database
        databases_download_dbsnp(
            assemblies=assemblies_list,
            dbsnp_folder=dbsnp_folder,
            dbsnp_vcf=False,
            dbsnp_parquet=False,
            genomes_folder=genomes_folder,
            threads=threads,
        )

        # Check
        downloaded_files = os.listdir(dbsnp_folder)
        for assembly in assemblies_list:
            assert assembly in downloaded_files
            downloaded_assembly_files = os.listdir(f"{dbsnp_folder}/{assembly}")
            assert "default" in downloaded_assembly_files
            for dbsnp_release in dbsnp_releases:
                assert dbsnp_release in downloaded_assembly_files
                downloaded_assembly_release_files = os.listdir(
                    f"{dbsnp_folder}/{assembly}/{dbsnp_release}"
                )
                expected_files = [
                    f"GCF_000001405.{dnsnp_assemblies_map.get(assembly)}.gz"
                ]
                for expected_file in expected_files:
                    if expected_file not in downloaded_assembly_release_files:
                        assert False
                assert True
