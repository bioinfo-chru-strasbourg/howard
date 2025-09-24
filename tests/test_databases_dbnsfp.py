# -*- coding: utf-8 -*-
"""
Tests

Usage:
pytest tests/

Coverage:
coverage run -m pytest tests/test_databases.py -x -vv --log-cli-level=DEBUG --capture=tee-sys
coverage report --include=howard/* -m
"""

import os
import shutil
from tempfile import TemporaryDirectory
from howard.functions.databases import databases_download_dbnsfp
import pytest  # type: ignore

from howard.objects.variants import Variants

from test_needed import (
    download_needed_databases,
    tests_folder,
    tests_databases_folder,
    tests_config,
)


def test_database_dbnsfp_step_by_step_from_source():
    """
    This function tests the "databases" function with a set of arguments.
    """

    # Init
    dbnsfp_source = os.path.join(tests_databases_folder, "dbnsfp", "dbNSFP4.4a.zip")
    genomes_folder = tests_config["folders"]["databases"]["genomes"]
    download_needed_databases()
    uniquify = False
    threads = 8

    # Tmp folder
    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # Assembly
        assemblies = "hg19,hg38"
        assemblies_list = [value for value in assemblies.split(",")]

        dbnsfp_folder = tmp_dir

        # Try to convert
        try:
            databases_download_dbnsfp(
                assemblies=assemblies_list,
                dbnsfp_folder=dbnsfp_folder,
                dbnsfp_source=dbnsfp_source,
                generate_parquet_file=False,
                generate_sub_databases=False,
                generate_vcf_file=False,
                genomes_folder=genomes_folder,
                uniquify=uniquify,
                threads=threads,
            )
        except:
            assert False

        downloaded_files = os.listdir(dbnsfp_folder)
        for assembly in assemblies_list:
            assert assembly in downloaded_files
            downloaded_assembly_files = os.listdir(f"{dbnsfp_folder}/{assembly}")
            nb_files = 2
            assert len(downloaded_assembly_files) == nb_files


@pytest.mark.parametrize(
    "generate_parquet_file, generate_sub_databases, generate_vcf_file, add_info, only_info, uniquify",
    [
        (
            generate_parquet_file,
            generate_sub_databases,
            generate_vcf_file,
            add_info,
            only_info,
            uniquify,
        )
        for generate_parquet_file in [True, False]
        for generate_sub_databases in [True, False]
        for generate_vcf_file in [True, False]
        for add_info in [True, False]
        for only_info in [True, False]
        for uniquify in [True, False]
    ],
)
def test_database_dbnsfp_options_full(
    generate_parquet_file,
    generate_sub_databases,
    generate_vcf_file,
    add_info,
    only_info,
    uniquify,
):
    """
    This function tests the "databases" function with a set of arguments.
    """

    # Init
    dbnsfp_source = os.path.join(tests_databases_folder, "dbnsfp", "dbNSFP4.4a.zip")
    genomes_folder = tests_config["folders"]["databases"]["genomes"]
    download_needed_databases()

    # Tmp folder
    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # Assembly
        assemblies = "hg19"
        assemblies_list = [value for value in assemblies.split(",")]

        # Download dbnsfp simulation
        dbnsfp_target = os.path.join(tmp_dir, "dbNSFP4.4a.zip")
        shutil.copy(dbnsfp_source, dbnsfp_target)

        dbnsfp_folder = tmp_dir

        # Try to convert
        try:
            databases_download_dbnsfp(
                assemblies=assemblies_list,
                dbnsfp_folder=dbnsfp_folder,
                genomes_folder=genomes_folder,
                generate_parquet_file=generate_parquet_file,
                generate_sub_databases=generate_sub_databases,
                generate_vcf_file=generate_vcf_file,
                add_info=add_info,
                only_info=only_info,
                uniquify=uniquify,
                threads=2,
            )
        except:
            assert False


def test_database_dbnsfp_step_by_step():
    """
    This function tests the "databases" function with a set of arguments.
    """

    # Init
    dbnsfp_source = os.path.join(tests_databases_folder, "dbnsfp", "dbNSFP4.4a.zip")
    genomes_folder = tests_config["folders"]["databases"]["genomes"]
    download_needed_databases()
    uniquify = False
    threads = 2

    # Tmp folder
    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # Assembly
        assemblies = "hg19,hg38"
        assemblies_list = [value for value in assemblies.split(",")]

        # Download dbnsfp simulation
        dbnsfp_target = os.path.join(tmp_dir, "dbNSFP4.4a.zip")
        shutil.copy(dbnsfp_source, dbnsfp_target)

        dbnsfp_folder = tmp_dir

        # Try to convert
        try:
            databases_download_dbnsfp(
                assemblies=assemblies_list,
                dbnsfp_folder=dbnsfp_folder,
                generate_parquet_file=False,
                generate_sub_databases=False,
                generate_vcf_file=False,
                genomes_folder=genomes_folder,
                uniquify=uniquify,
                threads=threads,
            )
        except:
            assert False

        downloaded_files = os.listdir(dbnsfp_folder)
        for assembly in assemblies_list:
            assert assembly in downloaded_files
            downloaded_assembly_files = os.listdir(f"{dbnsfp_folder}/{assembly}")
            nb_files = 2
            assert len(downloaded_assembly_files) == nb_files

        # Try again to generate parquet
        try:
            databases_download_dbnsfp(
                assemblies=assemblies_list,
                dbnsfp_folder=dbnsfp_folder,
                generate_parquet_file=True,
                generate_sub_databases=False,
                generate_vcf_file=False,
                genomes_folder=genomes_folder,
                uniquify=uniquify,
                threads=threads,
            )
        except:
            assert False

        downloaded_files = os.listdir(dbnsfp_folder)
        for assembly in assemblies_list:
            assert assembly in downloaded_files
            downloaded_assembly_files = os.listdir(f"{dbnsfp_folder}/{assembly}")
            nb_files = 4
            assert len(downloaded_assembly_files) == nb_files

        # Try again to generate sub databases parquet
        try:
            databases_download_dbnsfp(
                assemblies=assemblies_list,
                dbnsfp_folder=dbnsfp_folder,
                generate_parquet_file=True,
                generate_sub_databases=True,
                generate_vcf_file=False,
                genomes_folder=genomes_folder,
                uniquify=uniquify,
                threads=threads,
            )
        except:
            assert False

        downloaded_files = os.listdir(dbnsfp_folder)
        for assembly in assemblies_list:
            assert assembly in downloaded_files
            downloaded_assembly_files = os.listdir(f"{dbnsfp_folder}/{assembly}")
            nb_files = 108
            assert len(downloaded_assembly_files) == nb_files

        # Try again to generate VCF
        try:
            databases_download_dbnsfp(
                assemblies=assemblies_list,
                dbnsfp_folder=dbnsfp_folder,
                generate_parquet_file=True,
                generate_sub_databases=True,
                generate_vcf_file=True,
                genomes_folder=genomes_folder,
                uniquify=uniquify,
                threads=threads,
            )
        except:
            assert False

        downloaded_files = os.listdir(dbnsfp_folder)
        for assembly in assemblies_list:
            assert assembly in downloaded_files
            downloaded_assembly_files = os.listdir(f"{dbnsfp_folder}/{assembly}")
            nb_files = 189
            assert len(downloaded_assembly_files) == nb_files

        # Try again to generate nothing more
        try:
            databases_download_dbnsfp(
                assemblies=assemblies_list,
                dbnsfp_folder=dbnsfp_folder,
                generate_parquet_file=True,
                generate_sub_databases=True,
                generate_vcf_file=True,
                genomes_folder=genomes_folder,
                uniquify=uniquify,
                threads=threads,
            )
        except:
            assert False

        downloaded_files = os.listdir(dbnsfp_folder)
        for assembly in assemblies_list:
            assert assembly in downloaded_files
            downloaded_assembly_files = os.listdir(f"{dbnsfp_folder}/{assembly}")
            nb_files = 189
            assert len(downloaded_assembly_files) == nb_files


def test_database_dbnsfp_without_parquet_size_10k():
    """
    This function tests the "databases" function with a set of arguments.
    """

    # Init
    dbnsfp_source = os.path.join(tests_databases_folder, "dbnsfp", "dbNSFP4.4a.zip")
    genomes_folder = tests_config["folders"]["databases"]["genomes"]
    download_needed_databases()
    threads = 2

    # Tmp folder
    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # Assembly
        assemblies = "hg19"
        assemblies_list = [value for value in assemblies.split(",")]

        # Download dbnsfp simulation
        dbnsfp_target = os.path.join(tmp_dir, "dbNSFP4.4a.zip")
        shutil.copy(dbnsfp_source, dbnsfp_target)

        dbnsfp_folder = tmp_dir

        # Try to generate all files in one time with parquet size of 1Mb
        try:
            databases_download_dbnsfp(
                assemblies=assemblies_list,
                dbnsfp_folder=dbnsfp_folder,
                generate_parquet_file=True,
                generate_sub_databases=True,
                generate_vcf_file=True,
                parquet_size=0.01,
                genomes_folder=genomes_folder,
                threads=threads,
            )
        except:
            assert False

        downloaded_files = os.listdir(dbnsfp_folder)
        for assembly in assemblies_list:
            assert assembly in downloaded_files
            downloaded_assembly_files = os.listdir(f"{dbnsfp_folder}/{assembly}")
            nb_files = 189
            assert len(downloaded_assembly_files) == nb_files


def test_database_dbnsfp_parquet_with_info():
    """
    This function tests the "databases" function with a set of arguments.
    """

    # Init
    dbnsfp_source = os.path.join(tests_databases_folder, "dbnsfp", "dbNSFP4.4a.zip")
    genomes_folder = tests_config["folders"]["databases"]["genomes"]
    download_needed_databases()
    threads = 2

    # Tmp folder
    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # Assembly
        assemblies = "hg19"
        assemblies_list = [value for value in assemblies.split(",")]

        # Download dbnsfp simulation
        dbnsfp_target = os.path.join(tmp_dir, "dbNSFP4.4a.zip")
        shutil.copy(dbnsfp_source, dbnsfp_target)

        dbnsfp_folder = tmp_dir

        # Try to generate ALL and sub-database parquet folders with INFO column
        try:
            databases_download_dbnsfp(
                assemblies=assemblies_list,
                dbnsfp_folder=dbnsfp_folder,
                generate_parquet_file=True,
                generate_sub_databases=False,
                generate_vcf_file=False,
                add_info=True,
                genomes_folder=genomes_folder,
                threads=threads,
            )
        except:
            assert False

        db_file = f"{dbnsfp_folder}/hg19/dbNSFP4.4a.ALL.partition.parquet"
        db = Variants(input=db_file, load=True)
        query = " SELECT * FROM (DESCRIBE SELECT * FROM variants) WHERE column_name == 'INFO'"
        res = db.get_query_to_df(query)
        assert len(res) == 1
