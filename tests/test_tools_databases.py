# -*- coding: utf-8 -*-
"""
Tests

Usage:
pytest tests/

Coverage:
coverage run -m pytest . -x -v --log-cli-level=INFO --capture=tee-sys
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
from howard.commons import *
from howard.tools.databases import *
from test_needed import *





def test_download_hgmd():
    """
    The function `test_download_hgmd` downloads HGMD files for specified assemblies and checks if the
    downloaded files match the expected files.
    """

    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # Assembly
        assemblies = 'hg19'
        assemblies_list = [value for value in assemblies.split(',')]

        # Exomiser folder
        hgmd_folder = tmp_dir

        # HGMD conversion
        hgmd_file_hg19 = os.path.join(tests_databases_folder, "hgmd", "HGMD_TEST_hg19.vcf.gz")
        databases_download_hgmd(assemblies=assemblies_list, hgmd_file=hgmd_file_hg19, hgmd_folder=hgmd_folder, threads=8)

        # Check
        downloaded_files = os.listdir(hgmd_folder)
        for assembly in assemblies_list:
            assert assembly in downloaded_files
            downloaded_assembly_files = os.listdir(f"{hgmd_folder}/{assembly}")
            expected_files = ['HGMD_TEST.vcf.gz.tbi', 'HGMD_TEST.tsv.hdr', 'HGMD_TEST.parquet.hdr', 'HGMD_TEST.vcf.gz', 'HGMD_TEST.parquet', 'HGMD_TEST.tsv']
            for expected_file in expected_files:
                if expected_file not in downloaded_assembly_files:
                    assert False
            assert True


def test_databases_download_snpeff():
    """
    The function `test_databases_download_snpeff` downloads and prepares the snpEff database for specified
    assemblies.
    """

    # Full database generation, hg19 only (due to lack of hg38 assembly in tests)
    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # Assembly
        assemblies = 'hg19'
        assemblies_list = [value for value in assemblies.split(',')]

        # snpEff folder
        snpeff_folder = tmp_dir

        # Download and prepare database
        databases_download_snpeff(folder=snpeff_folder, assemblies=assemblies_list, config=tests_config)

        # Check
        downloaded_files = os.listdir(snpeff_folder)
        for assembly in assemblies_list:
            assert assembly in downloaded_files
            downloaded_assembly_files = os.listdir(f"{snpeff_folder}/{assembly}")
            expected_files = ['sequence.bin']
            for expected_file in expected_files:
                if expected_file not in downloaded_assembly_files:
                    assert False
            assert True


def test_databases_download_dbsnp():
    """
    The function `test_databases_download_dbsnp` downloads and prepares the dbsnp database for specified
    assemblies.
    """

    # Genomes
    genomes_folder = tests_config["folders"]["databases"]["genomes"]
    download_needed_databases()

    # Full database generation, hg19 only (due to lack of hg38 assembly in tests)
    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # Assembly
        assemblies = 'hg19'
        assemblies_list = [value for value in assemblies.split(',')]

        # Releases
        dbsnp_releases = ["b156"]

        # Threads
        threads = 2

        # Exomiser folder
        dbsnp_folder = tmp_dir
        
        # Download exomiser simulation
        dnsnp_assemblies_map:dict = {"hg19": "25", "hg38": "40"}
        for assembly in assemblies_list:
            for dbsnp_release in dbsnp_releases:
                dbsnp_data_source = os.path.join(tests_databases_folder, "dbsnp", f"GCF_000001405.{dnsnp_assemblies_map.get(assembly)}.gz")
                dbsnp_data_target = os.path.join(dbsnp_folder, assembly, dbsnp_release, f"GCF_000001405.{dnsnp_assemblies_map.get(assembly)}.gz")
                if not os.path.exists(os.path.join(dbsnp_folder, assembly, dbsnp_release)):
                    Path(os.path.join(dbsnp_folder, assembly, dbsnp_release)).mkdir(parents=True, exist_ok=True)
                if not os.path.exists(os.path.join(dbsnp_folder, assembly, dbsnp_data_target)):
                    shutil.copy(dbsnp_data_source, dbsnp_data_target)

        # Download and prepare database
        databases_download_dbsnp(assemblies=assemblies_list, dbsnp_folder=dbsnp_folder, threads=threads, dbsnp_vcf=True, dbsnp_parquet=True, genomes_folder=genomes_folder)

        # Check
        downloaded_files = os.listdir(dbsnp_folder)
        for assembly in assemblies_list:
            assert assembly in downloaded_files
            downloaded_assembly_files = os.listdir(f"{dbsnp_folder}/{assembly}")
            assert "default" in downloaded_assembly_files
            for dbsnp_release in dbsnp_releases:
                assert dbsnp_release in downloaded_assembly_files
                downloaded_assembly_release_files = os.listdir(f"{dbsnp_folder}/{assembly}/{dbsnp_release}")
                expected_files = [f'GCF_000001405.{dnsnp_assemblies_map.get(assembly)}.gz', 'dbsnp.vcf.gz', 'dbsnp.parquet.hdr', 'dbsnp.parquet']
                for expected_file in expected_files:
                    if expected_file not in downloaded_assembly_release_files:
                        assert False
                assert True

        # Download and prepare database again
        databases_download_dbsnp(assemblies=assemblies_list, dbsnp_folder=dbsnp_folder, threads=threads, dbsnp_vcf=True, dbsnp_parquet=True, genomes_folder=genomes_folder)

        # Check
        downloaded_files = os.listdir(dbsnp_folder)
        for assembly in assemblies_list:
            assert assembly in downloaded_files
            downloaded_assembly_files = os.listdir(f"{dbsnp_folder}/{assembly}")
            assert "default" in downloaded_assembly_files
            for dbsnp_release in dbsnp_releases:
                assert dbsnp_release in downloaded_assembly_files
                downloaded_assembly_release_files = os.listdir(f"{dbsnp_folder}/{assembly}/{dbsnp_release}")
                expected_files = [f'GCF_000001405.{dnsnp_assemblies_map.get(assembly)}.gz', 'dbsnp.vcf.gz', 'dbsnp.parquet.hdr', 'dbsnp.parquet']
                for expected_file in expected_files:
                    if expected_file not in downloaded_assembly_release_files:
                        assert False
                assert True


    # Multi assembly (without VCF and PArquet generation, due to lack of assembly hg38 in tests)
    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # Assembly
        assemblies = 'hg19,hg38'
        assemblies_list = [value for value in assemblies.split(',')]

        # Releases
        dbsnp_releases = ["b156"]

        # Threads
        threads = 2

        # Exomiser folder
        dbsnp_folder = tmp_dir
        
        # Download exomiser simulation
        dnsnp_assemblies_map:dict = {"hg19": "25", "hg38": "40"}
        for assembly in assemblies_list:
            for dbsnp_release in dbsnp_releases:
                dbsnp_data_source = os.path.join(tests_databases_folder, "dbsnp", f"GCF_000001405.{dnsnp_assemblies_map.get(assembly)}.gz")
                dbsnp_data_target = os.path.join(dbsnp_folder, assembly, dbsnp_release, f"GCF_000001405.{dnsnp_assemblies_map.get(assembly)}.gz")
                if not os.path.exists(os.path.join(dbsnp_folder, assembly, dbsnp_release)):
                    Path(os.path.join(dbsnp_folder, assembly, dbsnp_release)).mkdir(parents=True, exist_ok=True)
                if not os.path.exists(os.path.join(dbsnp_folder, assembly, dbsnp_data_target)):
                    shutil.copy(dbsnp_data_source, dbsnp_data_target)

        # Download and prepare database
        databases_download_dbsnp(assemblies=assemblies_list, dbsnp_folder=dbsnp_folder, threads=threads, dbsnp_vcf=False, dbsnp_parquet=False, genomes_folder=genomes_folder)

        # Check
        downloaded_files = os.listdir(dbsnp_folder)
        for assembly in assemblies_list:
            assert assembly in downloaded_files
            downloaded_assembly_files = os.listdir(f"{dbsnp_folder}/{assembly}")
            assert "default" in downloaded_assembly_files
            for dbsnp_release in dbsnp_releases:
                assert dbsnp_release in downloaded_assembly_files
                downloaded_assembly_release_files = os.listdir(f"{dbsnp_folder}/{assembly}/{dbsnp_release}")
                expected_files = [f'GCF_000001405.{dnsnp_assemblies_map.get(assembly)}.gz']
                for expected_file in expected_files:
                    if expected_file not in downloaded_assembly_release_files:
                        assert False
                assert True


def test_databases_download_exomiser():
    """
    The function `test_databases_download_exomiser` tests the `databases_download_exomiser` function by
    checking if the downloaded files match the expected number and if the specified assemblies are
    present.
    """

    # Init
    exomiser_data_hg19_source = os.path.join(tests_databases_folder, "exomiser", "test_hg19.zip")
    exomiser_data_hg38_source = os.path.join(tests_databases_folder, "exomiser", "test_hg38.zip")
    exomiser_phenotype_source = os.path.join(tests_databases_folder, "exomiser", "test_phenotype.zip")

    # Tmp folder
    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # Assembly
        assemblies = 'hg19,hg38'
        assemblies_list = [value for value in assemblies.split(',')]

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
        databases_download_exomiser(assemblies=assemblies_list, exomiser_folder=exomiser_folder, exomiser_release=exomiser_release, exomiser_phenotype_release=exomiser_phenotype_release)

        # Check
        downloaded_files = os.listdir(exomiser_folder)
        for assembly in assemblies_list:
            assert assembly in downloaded_files
            downloaded_assembly_files = os.listdir(f"{exomiser_folder}/{assembly}")
            expected_files = [f'test_{assembly}', 'application.properties', 'test_phenotype']
            for expected_file in expected_files:
                if expected_file not in downloaded_assembly_files:
                    assert False
            assert True


def test_download_alphamissense():
    """
    The function `test_download_alphamissense` tests the `databases_download_alphamissense` function by
    downloading AlphaMissense databases for specified assemblies and checking if the files are
    downloaded correctly.
    """

    # Tmp folder
    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # Assembly
        assemblies = 'hg19'
        assemblies_list = [value for value in assemblies.split(',')]

        # AlphaMissense folder
        alphamissense_folder = os.path.join(tmp_dir,"alphamissense")

        # Threads
        threads=4

        # Download AlphaMissense
        databases_download_alphamissense(assemblies=assemblies_list, alphamissense_folder=alphamissense_folder, threads=threads)

        # Check
        downloaded_files = os.listdir(alphamissense_folder)
        for assembly in assemblies_list:
            assert assembly in downloaded_files
            downloaded_assembly_files = os.listdir(f"{alphamissense_folder}/{assembly}")
            log.debug(downloaded_assembly_files)
            nb_files = 3
            assert len(downloaded_assembly_files) == nb_files

        # Download AlphaMissense again
        databases_download_alphamissense(assemblies=assemblies_list, alphamissense_folder=alphamissense_folder, threads=threads)

        # Check
        downloaded_files = os.listdir(alphamissense_folder)
        for assembly in assemblies_list:
            assert assembly in downloaded_files
            downloaded_assembly_files = os.listdir(f"{alphamissense_folder}/{assembly}")
            log.debug(downloaded_assembly_files)
            nb_files = 3
            assert len(downloaded_assembly_files) == nb_files


def test_database_dbnsfp():
    """
    This function tests the "databases" function with a set of arguments.
    """

    # Init
    dbnsfp_source = os.path.join(tests_databases_folder, "dbnsfp", "dbNSFP4.4a.zip")
    genomes_folder = tests_config["folders"]["databases"]["genomes"]

    # Tmp folder
    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # Assembly
        assemblies = 'hg19,hg38'
        assemblies_list = [value for value in assemblies.split(',')]

        # Download dbnsfp simulation
        dbnsfp_target = os.path.join(tmp_dir, "dbNSFP4.4a.zip")
        shutil.copy(dbnsfp_source, dbnsfp_target)

        dbnsfp_folder = tmp_dir

        # Try to convert 
        try:
            databases_download_dbnsfp(assemblies=assemblies_list, dbnsfp_folder=dbnsfp_folder, generate_parquet_file = False, generate_sub_databases = False, generate_vcf_file = False, genomes_folder=genomes_folder)
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
            databases_download_dbnsfp(assemblies=assemblies_list, dbnsfp_folder=dbnsfp_folder, generate_parquet_file = True, generate_sub_databases = False, generate_vcf_file = False, genomes_folder=genomes_folder)
        except:
            assert False

        downloaded_files = os.listdir(dbnsfp_folder)
        for assembly in assemblies_list:
            assert assembly in downloaded_files
            downloaded_assembly_files = os.listdir(f"{dbnsfp_folder}/{assembly}")
            nb_files = 4
            assert len(downloaded_assembly_files) == nb_files


        # Try again to generate parquet
        try:
            databases_download_dbnsfp(assemblies=assemblies_list, dbnsfp_folder=dbnsfp_folder, generate_parquet_file = True, generate_sub_databases = True, generate_vcf_file = False, genomes_folder=genomes_folder)
        except:
            assert False

        downloaded_files = os.listdir(dbnsfp_folder)
        for assembly in assemblies_list:
            assert assembly in downloaded_files
            downloaded_assembly_files = os.listdir(f"{dbnsfp_folder}/{assembly}")
            nb_files = 316
            assert len(downloaded_assembly_files) == nb_files

        # Try again to generate VCF
        try:
            databases_download_dbnsfp(assemblies=assemblies_list, dbnsfp_folder=dbnsfp_folder, generate_parquet_file = True, generate_sub_databases = True, generate_vcf_file = True, genomes_folder=genomes_folder)
        except:
            assert False

        downloaded_files = os.listdir(dbnsfp_folder)
        for assembly in assemblies_list:
            assert assembly in downloaded_files
            downloaded_assembly_files = os.listdir(f"{dbnsfp_folder}/{assembly}")
            nb_files = 553
            assert len(downloaded_assembly_files) == nb_files

        # Try again to generate nothing more
        try:
            databases_download_dbnsfp(assemblies=assemblies_list, dbnsfp_folder=dbnsfp_folder, generate_parquet_file = True, generate_sub_databases = True, generate_vcf_file = True, genomes_folder=genomes_folder)
        except:
            assert False

        downloaded_files = os.listdir(dbnsfp_folder)
        for assembly in assemblies_list:
            assert assembly in downloaded_files
            downloaded_assembly_files = os.listdir(f"{dbnsfp_folder}/{assembly}")
            nb_files = 553
            assert len(downloaded_assembly_files) == nb_files


    # Tmp folder
    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # Assembly
        assemblies = 'hg19'
        assemblies_list = [value for value in assemblies.split(',')]

        # Download dbnsfp simulation
        dbnsfp_target = os.path.join(tmp_dir, "dbNSFP4.4a.zip")
        shutil.copy(dbnsfp_source, dbnsfp_target)

        dbnsfp_folder = tmp_dir

        # Try to generate vcf file without parquet
        try:
            databases_download_dbnsfp(assemblies=assemblies_list, dbnsfp_folder=dbnsfp_folder, generate_parquet_file = False, generate_sub_databases = True, generate_vcf_file = True, genomes_folder=genomes_folder)
        except:
            assert False

        downloaded_files = os.listdir(dbnsfp_folder)
        for assembly in assemblies_list:
            assert assembly in downloaded_files
            downloaded_assembly_files = os.listdir(f"{dbnsfp_folder}/{assembly}")
            nb_files = 395
            assert len(downloaded_assembly_files) == nb_files


    # Tmp folder
    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # Assembly
        assemblies = 'hg19'
        assemblies_list = [value for value in assemblies.split(',')]

        # Download dbnsfp simulation
        dbnsfp_target = os.path.join(tmp_dir, "dbNSFP4.4a.zip")
        shutil.copy(dbnsfp_source, dbnsfp_target)

        dbnsfp_folder = tmp_dir

        # Try to generate all files in one time with parquet size of 1Mb
        try:
            databases_download_dbnsfp(assemblies=assemblies_list, dbnsfp_folder=dbnsfp_folder, generate_parquet_file = True, generate_sub_databases = True, generate_vcf_file = True, parquet_size = 1, genomes_folder=genomes_folder)
        except:
            assert False

        downloaded_files = os.listdir(dbnsfp_folder)
        for assembly in assemblies_list:
            assert assembly in downloaded_files
            downloaded_assembly_files = os.listdir(f"{dbnsfp_folder}/{assembly}")
            nb_files = 553
            assert len(downloaded_assembly_files) == nb_files


    # Tmp folder
    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # Assembly
        assemblies = 'hg19'
        assemblies_list = [value for value in assemblies.split(',')]

        # Download dbnsfp simulation
        dbnsfp_target = os.path.join(tmp_dir, "dbNSFP4.4a.zip")
        shutil.copy(dbnsfp_source, dbnsfp_target)

        dbnsfp_folder = tmp_dir

        # Try to generate ALL and sub-database parquet folders but only sub-database parquet files
        try:
            databases_download_dbnsfp(assemblies=assemblies_list, dbnsfp_folder=dbnsfp_folder, generate_parquet_file = True, generate_sub_databases = True, generate_vcf_file = False, not_generate_files_all = True, genomes_folder=genomes_folder)
        except:
            assert False

        downloaded_files = os.listdir(dbnsfp_folder)
        for assembly in assemblies_list:
            assert assembly in downloaded_files
            downloaded_assembly_files = os.listdir(f"{dbnsfp_folder}/{assembly}")
            nb_files = 314
            assert len(downloaded_assembly_files) == nb_files


    # Tmp folder
    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # Assembly
        assemblies = 'hg19'
        assemblies_list = [value for value in assemblies.split(',')]

        # Download dbnsfp simulation
        dbnsfp_target = os.path.join(tmp_dir, "dbNSFP4.4a.zip")
        shutil.copy(dbnsfp_source, dbnsfp_target)

        dbnsfp_folder = tmp_dir

        # Try to generate ALL and sub-database parquet folders with INFO column
        try:
            databases_download_dbnsfp(assemblies=assemblies_list, dbnsfp_folder=dbnsfp_folder, generate_parquet_file = True, generate_sub_databases = True, generate_vcf_file = True, add_info=True, genomes_folder=genomes_folder)
        except:
            assert False

        downloaded_files = os.listdir(dbnsfp_folder)
        for assembly in assemblies_list:
            assert assembly in downloaded_files
            downloaded_assembly_files = os.listdir(f"{dbnsfp_folder}/{assembly}")
            nb_files = 553
            assert len(downloaded_assembly_files) == nb_files


def test_database():
    """
    This function tests the "databases" function with a set of arguments.
    """

    # Prepare arguments for the query function
    args = argparse.Namespace(
        assembly = 'hg19',
        download_genomes = None,
        download_annovar_files = None,
        download_annovar = None,
        download_annovar_url = None,
        download_snpeff = None,
        download_refseq = None,
        download_dbnsfp = None,
        download_alphamissense = None,
        download_exomiser = None,
        download_dbsnp = None,
        convert_hgmd = None,
        config = None,
    )

    try:
        databases(args)
        assert True
    except:
        assert False
    

def test_databases_download():
    """
    This function tests the download of databases for Annovar and snpEff tools.
    """

    # Tmp folder
    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # assembly
        assemblies = 'hg19'
        assemblies_list = [value for value in assemblies.split(',')]

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
        annovar_file_list = 'nci60'
        annovar_file_list_list = [value for value in annovar_file_list.split(',')]

        # config
        config = {"tools": {"snpeff": {"jar": DEFAULT_SNPEFF_BIN}}}

        # threads
        threads = 1

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
            download_dbnsfp = None, # Too long...
            download_alphamissense = None,
            download_exomiser = None,
            download_dbsnp = None,
            convert_hgmd = None,
            config=config,
            threads=threads
        )

        # Download
        try:
            databases_download(args)
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


def test_databases_download_genomes_only():
    """
    This function tests the download of genomes.
    """

    # Tmp folder
    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # assembly
        assemblies = 'sacCer3'
        assemblies_list = [value for value in assemblies.split(',')]

        # Genome
        genome_provider = None
        genome_contig_regex = None

        # threads
        threads = 1

        # Arguments
        args = argparse.Namespace(
            assembly=assemblies,
            download_genomes=tmp_dir,
            download_genomes_provider=genome_provider,
            download_genomes_contig_regex=genome_contig_regex,
            download_annovar=None,
            download_snpeff=None,
            download_refseq=None,
            download_dbnsfp = None,
            download_alphamissense = None,
            download_exomiser = None,
            download_dbsnp = None,
            convert_hgmd = None,
            threads=threads
        )

        # Download
        databases_download(args)

        # Dowloaded files
        downloaded_files = os.listdir(tmp_dir)
        
        # Check Genome
        for assembly in assemblies_list:
            assert assembly in downloaded_files
            downloaded_assembly_files = os.listdir(f"{tmp_dir}/{assembly}")
            genome_file = f"{assembly}.fa"
            assert genome_file in downloaded_assembly_files


def test_databases_format_refseq():
    """
    The function `test_databases_format_refseq` tests the `databases_format_refseq` function by
    formatting a refSeq file in different ways and checking the output.
    """

    # Test downloading one file one assembly, and format in BED
    with TemporaryDirectory(dir=tests_folder) as tmp_dir:
       
        # refSeq files
        refseq_file_to_format = f"{tests_databases_folder}/others/ncbiRefSeq.test.txt"

        # Format refSeq by default
        
        # Param
        refseq_file_formatted = f"{tmp_dir}/test1.bed"
        include_utr_5 = True
        include_utr_3 = True
        include_chrM = True
        include_non_canonical_chr = True
        include_non_coding_transcripts = True
        include_transcript_ver = True

        # Format
        databases_format_refseq(refseq_file=refseq_file_to_format, output_file=refseq_file_formatted, include_utr_3=include_utr_3, include_chrM=include_chrM, include_non_canonical_chr=include_non_canonical_chr, include_non_coding_transcripts=include_non_coding_transcripts, include_transcript_ver=include_transcript_ver)

        # Check
        assert os.path.exists(refseq_file_formatted) and os.stat(refseq_file_formatted).st_size > 0
        # Dataframe
        df = pd.read_csv(refseq_file_formatted, sep="\t", header=None)
        # Count number of transcript with version
        count_transcript_with_ver = df[df[4].str.contains(r'\.')].shape[0]
        # Count non canonical chromosome
        count_non_canonical_chr = sum("_" in key for key in set(df[0]))
        # Check
        assert len(df) > 0
        assert count_non_canonical_chr > 0
        assert count_transcript_with_ver > 0

        # Format refSeq all parameters False

        # Param
        refseq_file_formatted = f"{tmp_dir}/test2.bed"
        include_utr_5 = False
        include_utr_3 = False
        include_chrM = False
        include_non_canonical_chr = False
        include_non_coding_transcripts = False
        include_transcript_ver = False

        # Format
        databases_format_refseq(refseq_file=refseq_file_to_format, output_file=refseq_file_formatted, include_utr_5=include_utr_5, include_utr_3=include_utr_3, include_chrM=include_chrM, include_non_canonical_chr=include_non_canonical_chr, include_non_coding_transcripts=include_non_coding_transcripts, include_transcript_ver=include_transcript_ver)

        # Check
        assert os.path.exists(refseq_file_formatted) and os.stat(refseq_file_formatted).st_size > 0
        # Dataframe
        df = pd.read_csv(refseq_file_formatted, sep="\t", header=None)
        # Count number of transcript with version
        count_transcript_with_ver = df[df[4].str.contains(r'\.')].shape[0]
        # Start end position of the first gene
        first_gene_start, first_gene_end = min(df[df[3]==df[3][0]][1]), max(df[df[3]==df[3][0]][2])
        # Count non canonical chromosome
        count_non_canonical_chr = sum("_" in key for key in set(df[0]))
        # Check
        assert len(df) > 0
        assert count_non_canonical_chr == 0
        assert count_transcript_with_ver == 0
        
        # Format refSeq only UTR 5

        # Param
        refseq_file_formatted = f"{tmp_dir}/test3.bed"
        include_utr_5 = True
        include_utr_3 = False
        include_chrM = False
        include_non_canonical_chr = False
        include_non_coding_transcripts = False
        include_transcript_ver = False
        
        # Format
        databases_format_refseq(refseq_file=refseq_file_to_format, output_file=refseq_file_formatted, include_utr_5=include_utr_5, include_utr_3=include_utr_3, include_chrM=include_chrM, include_non_canonical_chr=include_non_canonical_chr, include_non_coding_transcripts=include_non_coding_transcripts, include_transcript_ver=include_transcript_ver)
        
        # Check
        assert os.path.exists(refseq_file_formatted) and os.stat(refseq_file_formatted).st_size > 0
        # Dataframe
        df = pd.read_csv(refseq_file_formatted, sep="\t", header=None)
        # Count number of transcript with version
        count_transcript_with_ver = df[df[4].str.contains(r'\.')].shape[0]
        # Start end position of the first gene
        first_gene_start_utr5, first_gene_end_utr5 = min(df[df[3]==df[3][0]][1]), max(df[df[3]==df[3][0]][2])
        # Count non canonical chromosome
        count_non_canonical_chr = sum("_" in key for key in set(df[0]))
        # Check
        assert first_gene_start_utr5 < first_gene_start and first_gene_end_utr5 == first_gene_end
        assert len(df) > 0
        assert count_non_canonical_chr == 0
        assert count_transcript_with_ver == 0

        # Format refSeq only UTR 3

        # Param
        refseq_file_formatted = f"{tmp_dir}/test4.bed"
        include_utr_5 = False
        include_utr_3 = True
        include_chrM = False
        include_non_canonical_chr = False
        include_non_coding_transcripts = False
        include_transcript_ver = False

        # Format
        databases_format_refseq(refseq_file=refseq_file_to_format, output_file=refseq_file_formatted, include_utr_5=include_utr_5, include_utr_3=include_utr_3, include_chrM=include_chrM, include_non_canonical_chr=include_non_canonical_chr, include_non_coding_transcripts=include_non_coding_transcripts, include_transcript_ver=include_transcript_ver)
        
        # Check
        assert os.path.exists(refseq_file_formatted) and os.stat(refseq_file_formatted).st_size > 0
        # Dataframe
        df = pd.read_csv(refseq_file_formatted, sep="\t", header=None)
        # Count number of transcript with version
        count_transcript_with_ver = df[df[4].str.contains(r'\.')].shape[0]
        # Start end position of the first gene
        first_gene_start_utr3, first_gene_end_utr3 = min(df[df[3]==df[3][0]][1]), max(df[df[3]==df[3][0]][2])
        # Count non canonical chromosome
        count_non_canonical_chr = sum("_" in key for key in set(df[0]))
        # Check
        assert first_gene_start_utr3 == first_gene_start and first_gene_end_utr3 > first_gene_end
        assert len(df) > 0
        assert count_non_canonical_chr == 0
        assert count_transcript_with_ver == 0


def test_databases_download_refseq():
    """
    The function `test_databases_download_refseq` tests the `databases_download_refseq` function by
    downloading different files for different assemblies and checking if the downloaded files match the
    expected files.
    """

    with TemporaryDirectory(dir=tests_folder) as tmp_dir:
       
        # Test downloading all files for 2 assemblies

        # assembly
        assemblies = ["hg19","hg38"]

        # refSeq Folder
        refseq_folder = f"{tmp_dir}/test1"

        # downloaded_files expected 
        downloaded_files_expected = {'hg19': ['ncbiRefSeq.txt', 'ncbiRefSeqLink.txt'], 'hg38': ['ncbiRefSeq.txt', 'ncbiRefSeqLink.txt']}

        # download
        downloaded_files = databases_download_refseq(assemblies=assemblies, refseq_folder=refseq_folder)
        assert downloaded_files == downloaded_files_expected

        # check if each files exists
        for assembly in downloaded_files:
            for file in downloaded_files[assembly]:
                refseq_file = f"{refseq_folder}/{assembly}/{file}"
                assert os.path.exists(refseq_file) and os.stat(refseq_file).st_size > 0
                
        # Test downloading all files no assembly
       
        # assembly
        assemblies = []

        # refSeq Folder
        refseq_folder = f"{tmp_dir}/test2"

        # downloaded_files expected 
        downloaded_files_expected = {}

        # download
        downloaded_files = databases_download_refseq(assemblies=assemblies, refseq_folder=refseq_folder)
        assert downloaded_files == downloaded_files_expected

        # check if each files exists
        for assembly in downloaded_files:
            for file in downloaded_files[assembly]:
                refseq_file = f"{refseq_folder}/{assembly}/{file}"
                assert os.path.exists(refseq_file) and os.stat(refseq_file).st_size > 0

        # Test downloading one file one assembly, and re-download test
       
        # assembly
        assemblies = ["hg19"]

        # refSeq Folder
        refseq_folder = f"{tmp_dir}/test3"

        # refSeq files
        refseq_files = ["ncbiRefSeq.txt"]

        # downloaded_files expected 
        downloaded_files_expected = {"hg19" : refseq_files}

        # download
        downloaded_files = databases_download_refseq(assemblies=assemblies, refseq_folder=refseq_folder, refseq_files=refseq_files)
        assert downloaded_files == downloaded_files_expected

        # check if each files exists
        for assembly in downloaded_files:
            for file in downloaded_files[assembly]:
                refseq_file = f"{refseq_folder}/{assembly}/{file}"
                assert os.path.exists(refseq_file) and os.stat(refseq_file).st_size > 0

        # Re-download
        downloaded_files = databases_download_refseq(assemblies=assemblies, refseq_folder=refseq_folder, refseq_files=refseq_files)
        assert downloaded_files == downloaded_files_expected

        # check if each files exists
        for assembly in downloaded_files:
            for file in downloaded_files[assembly]:
                refseq_file = f"{refseq_folder}/{assembly}/{file}"
                assert os.path.exists(refseq_file) and os.stat(refseq_file).st_size > 0

        # And format

        # Param
        refseq_format_file = "ncbiRefSeq.txt"
        refseq_format_file_output = f"{tmp_dir}/test.bed"
        include_utr_5 = True
        include_utr_3 = True
        include_chrM = True
        include_non_canonical_chr = True
        include_non_coding_transcripts = True
        include_transcript_ver = True

        # Format
        databases_download_refseq(assemblies=assemblies, refseq_folder=refseq_folder, refseq_files=refseq_files, refseq_format_file=refseq_format_file, refseq_format_file_output=refseq_format_file_output, include_utr_5=include_utr_5, include_utr_3=include_utr_3, include_chrM=include_chrM, include_non_canonical_chr=include_non_canonical_chr, include_non_coding_transcripts=include_non_coding_transcripts, include_transcript_ver=include_transcript_ver)
        
        # Check
        assert os.path.exists(refseq_format_file_output) and os.stat(refseq_format_file_output).st_size > 0


def test_databases_download_genomes():
    """
    The function tests the databases_download_genomes function by checking if genomes are downloaded correctly for
    different assemblies and contig filters.
    """

    import genomepy

    # Init
    assemblies_config = {
            "sacCer3": {
                "assembly": "sacCer3",
                "contigs": ['chrM', 'chrXI', 'chrII', 'chrXVI', 'chrIII', 'chrVI', 'chrV', 'chrXII', 'chrVIII', 'chrXV', 'chrIV', 'chrI', 'chrXIII', 'chrX', 'chrIX', 'chrVII', 'chrXIV']
            },
            "sacCer2": {
                "assembly": "sacCer2",
                "contigs": ['chrM', '2micron', 'chrXI', 'chrII', 'chrXVI', 'chrIII', 'chrVI', 'chrV', 'chrXII', 'chrVIII', 'chrXV', 'chrIV', 'chrI', 'chrXIII', 'chrX', 'chrIX', 'chrVII', 'chrXIV']
            }
        }

    # Uniq assembly not folder provided
    with tempfile.TemporaryDirectory() as tmpdir:

        assemblies = ["sacCer3"]
        
        genomes_folder = None
        provider = "UCSC"
        contig_regex = None
        threads = 1
        try:
            genome = databases_download_genomes(assemblies=assemblies, genomes_folder=genomes_folder, provider=provider, contig_regex=contig_regex, threads=threads)
            for assembly in assemblies:
                genome = genomepy.Genome(assembly, genomes_dir=DEFAULT_GENOME_FOLDER)
                assert os.path.exists(genome.genome_file)
                assert list(genome.keys()).sort() == assemblies_config.get(assembly).get("contigs",[]).sort()
        except:
            assert False

    # Uniq assembly
    with tempfile.TemporaryDirectory() as tmpdir:

        assemblies = ["sacCer3"]
        
        genomes_folder = tmpdir
        provider = "UCSC"
        contig_regex = None
        threads = 1
        try:
            genome = databases_download_genomes(assemblies=assemblies, genomes_folder=genomes_folder, provider=provider, contig_regex=contig_regex, threads=threads)
            for assembly in assemblies:
                genome = genomepy.Genome(assembly, genomes_dir=genomes_folder)
                assert os.path.exists(genome.genome_file)
                assert list(genome.keys()).sort() == assemblies_config.get(assembly).get("contigs",[]).sort()
        except:
            assert False

    # Multiple assemblies
    with tempfile.TemporaryDirectory() as tmpdir:

        assemblies = ["sacCer2", "sacCer3"]
        
        genomes_folder = tmpdir
        provider = "UCSC"
        contig_regex = None
        threads = 1
        try:
            genome = databases_download_genomes(assemblies=assemblies, genomes_folder=genomes_folder, provider=provider, contig_regex=contig_regex, threads=threads)
            for assembly in assemblies:
                genome = genomepy.Genome(assembly, genomes_dir=genomes_folder)
                assert os.path.exists(genome.genome_file)
                assert list(genome.keys()).sort() == assemblies_config.get(assembly).get("contigs",[]).sort()
        except:
            assert False

    # Filtered assembl
    with tempfile.TemporaryDirectory() as tmpdir:

        assemblies = ["sacCer3"]
        
        genomes_folder = tmpdir
        provider = "UCSC"
        contig_regex = "^>chrX.*$"
        threads = 1
        try:
            genome = databases_download_genomes(assemblies=assemblies, genomes_folder=genomes_folder, provider=provider, contig_regex=contig_regex, threads=threads)
            for assembly in assemblies:
                genome = genomepy.Genome(assembly, genomes_dir=genomes_folder)
                assert os.path.exists(genome.genome_file)
                assert list(genome.keys()).sort() == ['chrXI', 'chrXVI', 'chrXII', 'chrXV', 'chrXIII', 'chrX', 'chrXIV'].sort()
        except:
            assert False 


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
        file_list = ['nci60']
        
        # Download
        databases_download_annovar(folder=tmp_dir, files=file_list, assemblies=assemblies)
        
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
        databases_download_annovar(folder=tmp_dir, files=file_list, assemblies=assemblies)
        
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
        
        # Download
        databases_download_annovar(folder=tmp_dir, files=file_list, assemblies=assemblies)
        
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
        file_list = ['cosmic68*']
        
        # Download
        databases_download_annovar(folder=tmp_dir, files=file_list, assemblies=assemblies)
        
        # Dowloaded files
        downloaded_files = os.listdir(tmp_dir)
        
        # Check
        for assembly in assemblies:
            assert assembly in downloaded_files
            downloaded_files_assembly = os.listdir(f"{tmp_dir}/{assembly}")
            for file in file_list:
                downloaded_file = f"{assembly}_{file}.txt"
                filtered_files = fnmatch.filter(downloaded_files_assembly, downloaded_file)
                assert len(filtered_files) > 1
