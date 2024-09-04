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


def test_databases_infos():
    """
    The function `test_databases_infos` retrieves information about databases based on specified
    configurations and parameters.
    """

    databases_infos_dict = databases_infos(
        config=tests_config,
        database_folder_releases=["current"],
        database_formats=["parquet", "vcf"],
        assembly="hg19",
    )

    assert len(databases_infos_dict)


def test_databases_param():
    """
    The `test_databases_param` function retrieves information about databases based on specified
    configurations and parameters and asserts that the output files exist and contain the expected
    information.
    """

    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        output = os.path.join(tmp_dir, "param.json")
        output_description = os.path.join(tmp_dir, "description.json")

        databases_infos_dict = databases_infos(
            config=tests_config,
            database_folder_releases=["current"],
            database_formats=["parquet", "vcf"],
            assembly="hg19",
        )
        databases_param_dict = databases_param(
            databases_infos_dict=databases_infos_dict,
            output=output,
            output_description=output_description,
        )

        assert len(databases_infos_dict)
        assert len(databases_param_dict)
        assert len(databases_infos_dict) == databases_param_dict.get(
            "Number of databases"
        )
        assert os.path.exists(output)
        assert os.path.exists(output_description)


def test_databases():
    """
    This function tests the "databases" function with a set of arguments.
    """

    # Prepare arguments for the query function
    args = argparse.Namespace(
        assembly="hg19",
        download_genomes=None,
        download_annovar_files=None,
        download_annovar=None,
        download_annovar_url=None,
        download_snpeff=None,
        download_refseq=None,
        download_dbnsfp=None,
        download_alphamissense=None,
        download_exomiser=None,
        download_dbsnp=None,
        convert_hgmd=None,
        input_annovar=None,
        output_annovar=None,
        input_extann=None,
        output_extann=None,
        genome=None,
        generate_param=None,
        config=None,
        arguments_dict=arguments_dict,
    )

    try:
        databases(args)
        assert True
    except:
        assert False


def test_databases_download_genomes_only():
    """
    This function tests the download of genomes.
    """

    # Tmp folder
    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # assembly
        assemblies = "sacCer3"
        assemblies_list = [value for value in assemblies.split(",")]

        # Genome
        genome_provider = None
        genome_contig_regex = None

        # threads
        threads = 2

        # Arguments
        args = argparse.Namespace(
            assembly=assemblies,
            download_genomes=tmp_dir,
            download_genomes_provider=genome_provider,
            download_genomes_contig_regex=genome_contig_regex,
            download_annovar=None,
            download_snpeff=None,
            download_refseq=None,
            download_dbnsfp=None,
            download_alphamissense=None,
            download_exomiser=None,
            download_dbsnp=None,
            convert_hgmd=None,
            input_annovar=None,
            output_annovar=None,
            input_extann=None,
            output_extann=None,
            genome=None,
            generate_param=None,
            threads=threads,
            arguments_dict=arguments_dict,
        )

        # Download
        databases(args)

        # Dowloaded files
        downloaded_files = os.listdir(tmp_dir)

        # Check Genome
        for assembly in assemblies_list:
            assert assembly in downloaded_files
            downloaded_assembly_files = os.listdir(f"{tmp_dir}/{assembly}")
            genome_file = f"{assembly}.fa"
            assert genome_file in downloaded_assembly_files


def test_databases_download():
    """
    This function tests the download of databases for Annovar and snpEff tools.
    """

    # Tmp folder
    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # assembly
        assemblies = "hg19"
        assemblies_list = [value for value in assemblies.split(",")]

        # Genome
        genome_provider = None
        genome_contig_regex = None

        # tmp dir
        genomes_tmp_dir = os.path.join(tmp_dir, "genomes")
        os.mkdir(genomes_tmp_dir)
        annovar_tmp_dir = os.path.join(tmp_dir, "annovar")
        os.mkdir(annovar_tmp_dir)
        snpeff_tmp_dir = os.path.join(tmp_dir, "snpeff")
        os.mkdir(snpeff_tmp_dir)
        refseq_tmp_dir = os.path.join(tmp_dir, "refseq")
        os.mkdir(refseq_tmp_dir)

        # files
        annovar_file_list = "nci60"
        annovar_file_list_list = [value for value in annovar_file_list.split(",")]

        # config
        config = {"tools": {"snpeff": {"jar": tests_config.get("tools").get("snpeff")}}}

        # threads
        threads = 2

        # annovar URL
        download_annovar_url = DEFAULT_ANNOVAR_URL

        # Arguments
        args = argparse.Namespace(
            assembly=assemblies,
            download_genomes=genomes_tmp_dir,
            download_genomes_provider=genome_provider,
            download_genomes_contig_regex=genome_contig_regex,
            download_annovar=annovar_tmp_dir,
            download_annovar_files=annovar_file_list,
            download_annovar_url=download_annovar_url,
            download_snpeff=snpeff_tmp_dir,
            download_refseq=refseq_tmp_dir,
            download_refseq_files="ncbiRefSeq.txt,ncbiRefSeqLink.txt",
            download_refseq_url=None,
            download_refseq_prefix="ncbiRefSeq",
            download_refseq_format_file="ncbiRefSeq.txt",
            download_refseq_include_utr5=True,
            download_refseq_include_utr3=True,
            download_refseq_include_chrM=True,
            download_refseq_include_non_canonical_chr=True,
            download_refseq_include_non_coding_transcripts=True,
            download_refseq_include_transcript_version=True,
            download_dbnsfp=None,  # Too long...
            download_alphamissense=None,
            download_exomiser=None,
            download_dbsnp=None,
            convert_hgmd=None,
            input_annovar=None,
            output_annovar=None,
            input_extann=None,
            output_extann=None,
            genome=None,
            generate_param=None,
            config=config,
            threads=threads,
            arguments_dict=arguments_dict,
        )

        # Download
        try:
            databases(args)
            assert True
        except:
            assert False

        # Check Genome
        downloaded_files = os.listdir(genomes_tmp_dir)
        for assembly in assemblies_list:
            assert assembly in downloaded_files
            downloaded_assembly_files = os.listdir(f"{genomes_tmp_dir}/{assembly}")
            genome_file = f"{assembly}.fa"
            assert genome_file in downloaded_assembly_files

        # Check Annovar
        downloaded_files = os.listdir(annovar_tmp_dir)
        for assembly in assemblies_list:
            assert assembly in downloaded_files
            downloaded_files_assembly = os.listdir(f"{annovar_tmp_dir}/{assembly}")
            assert f"{assembly}_refGene.txt" in downloaded_files_assembly
            for file in annovar_file_list_list:
                downloaded_file = f"{assembly}_{file}.txt"
                assert downloaded_file in downloaded_files_assembly

        # Check snpEff databases list file
        downloaded_files = os.listdir(snpeff_tmp_dir)
        snpeff_databases_list_file = "snpeff_databases.list"
        assert snpeff_databases_list_file in downloaded_files
        # Check assembly folders of snpEff
        for assembly in assemblies_list:
            assert assembly in downloaded_files

        # Check refSeq
        downloaded_files = os.listdir(refseq_tmp_dir)
        for assembly in assemblies_list:
            assert assembly in downloaded_files
            downloaded_refseq_files = os.listdir(f"{refseq_tmp_dir}/{assembly}")
            refseq_file = "ncbiRefSeq.txt"
            assert refseq_file in downloaded_refseq_files
