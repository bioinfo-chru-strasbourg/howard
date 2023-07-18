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




def test_databases_format_refseq():
    """
    The function `test_databases_format_refseq` tests the `databases_format_refseq` function by
    formatting a refSeq file in different ways and checking the output.
    """

    # Test downloading one file one assembly, and format in BED
    with TemporaryDirectory(dir=".") as tmp_dir:
       
        # refSeq files
        refseq_file_to_format = f"{tests_folder}/data/ncbiRefSeq.test.txt"

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

    with TemporaryDirectory(dir=".") as tmp_dir:
       
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
