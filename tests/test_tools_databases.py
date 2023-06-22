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
        config = None,
    )

    try:
        databases(args)
    except:
        assert False

    assert True


def test_databases_download_genomes():
    """
    The function tests the databases_download_genomes function by checking if genomes are downloaded correctly for
    different assemblies and contig filters.
    """

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
        
        genome_folder = None
        provider = "UCSC"
        contig_regex = None
        threads = 1
        try:
            genome = databases_download_genomes(assemblies=assemblies, genome_folder=genome_folder, provider=provider, contig_regex=contig_regex, threads=threads)
            for assembly in assemblies:
                genome = genomepy.Genome(assembly, genomes_dir=DEFAULT_GENOME_FOLDER)
                assert os.path.exists(genome.genome_file)
                assert list(genome.keys()).sort() == assemblies_config.get(assembly).get("contigs",[]).sort()
        except:
            assert False

    # Uniq assembly
    with tempfile.TemporaryDirectory() as tmpdir:

        assemblies = ["sacCer3"]
        
        genome_folder = tmpdir
        provider = "UCSC"
        contig_regex = None
        threads = 1
        try:
            genome = databases_download_genomes(assemblies=assemblies, genome_folder=genome_folder, provider=provider, contig_regex=contig_regex, threads=threads)
            for assembly in assemblies:
                genome = genomepy.Genome(assembly, genomes_dir=genome_folder)
                assert os.path.exists(genome.genome_file)
                assert list(genome.keys()).sort() == assemblies_config.get(assembly).get("contigs",[]).sort()
        except:
            assert False

    # Multiple assemblies
    with tempfile.TemporaryDirectory() as tmpdir:

        assemblies = ["sacCer2", "sacCer3"]
        
        genome_folder = tmpdir
        provider = "UCSC"
        contig_regex = None
        threads = 1
        try:
            genome = databases_download_genomes(assemblies=assemblies, genome_folder=genome_folder, provider=provider, contig_regex=contig_regex, threads=threads)
            for assembly in assemblies:
                genome = genomepy.Genome(assembly, genomes_dir=genome_folder)
                assert os.path.exists(genome.genome_file)
                assert list(genome.keys()).sort() == assemblies_config.get(assembly).get("contigs",[]).sort()
        except:
            assert False

    # Filtered assembl
    with tempfile.TemporaryDirectory() as tmpdir:

        assemblies = ["sacCer3"]
        
        genome_folder = tmpdir
        provider = "UCSC"
        contig_regex = "^>chrX.*$"
        threads = 1
        try:
            genome = databases_download_genomes(assemblies=assemblies, genome_folder=genome_folder, provider=provider, contig_regex=contig_regex, threads=threads)
            for assembly in assemblies:
                genome = genomepy.Genome(assembly, genomes_dir=genome_folder)
                assert os.path.exists(genome.genome_file)
                assert list(genome.keys()).sort() == ['chrXI', 'chrXVI', 'chrXII', 'chrXV', 'chrXIII', 'chrX', 'chrXIV'].sort()
        except:
            assert False 


def test_databases_download_annovar():
    """
    This function tests the functionality of the databases_download_annovar function by downloading and
    checking various files.
    """

    # Test downloading an existing file
    with TemporaryDirectory(dir=".") as tmp_dir:
       
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
    with TemporaryDirectory(dir=".") as tmp_dir:
        
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
    with TemporaryDirectory(dir=".") as tmp_dir:
        
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


def test_databases_download_genomes_only():
    """
    This function tests the download of genomes.
    """

    # Tmp folder
    with TemporaryDirectory(dir=".") as tmp_dir:

        # assembly
        assemblies = 'sacCer3'
        assemblies_list = [value for value in assemblies.split(',')]

        # Genome
        genome_provider = None
        genome_contig_regex = None

        # Arguments
        args = argparse.Namespace(
            assembly=assemblies,
            download_genomes=tmp_dir,
            download_genomes_provider=genome_provider,
            download_genomes_contig_regex=genome_contig_regex,
            download_annovar=None,
            download_snpeff=None
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


def test_databases_download():
    """
    This function tests the download of databases for Annovar and snpEff tools.
    """

    # Tmp folder
    with TemporaryDirectory(dir=".") as tmp_dir:

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

        # files
        annovar_file_list = 'nci60'
        annovar_file_list_list = [value for value in annovar_file_list.split(',')]

        # config
        config = {"tools": {"snpeff": {"jar": DEFAULT_SNPEFF_BIN}}}

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
            config=config
        )

        # Download
        databases_download(args)

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
            for file in annovar_file_list_list:
                downloaded_file = f"{assembly}_{file}.txt"
                assert downloaded_file in downloaded_files

        # Check snpEff databases list file
        downloaded_files = os.listdir(snpeff_tmp_dir)
        snpeff_databases_list_file = "snpeff_databases.list"
        assert snpeff_databases_list_file in downloaded_files
        # Check assembly folders of snpEff
        for assembly in assemblies_list:
            assert assembly in downloaded_files
