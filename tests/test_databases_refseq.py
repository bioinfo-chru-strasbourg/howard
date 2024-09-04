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
        sort = False
        header = False
        header_first_line = True

        # Format
        databases_format_refseq(
            refseq_file=refseq_file_to_format,
            output_file=refseq_file_formatted,
            include_utr_3=include_utr_3,
            include_chrM=include_chrM,
            include_non_canonical_chr=include_non_canonical_chr,
            include_non_coding_transcripts=include_non_coding_transcripts,
            include_transcript_ver=include_transcript_ver,
            sort=sort,
            header=header,
            header_first_line=header_first_line,
        )

        # Check
        assert (
            os.path.exists(refseq_file_formatted)
            and os.stat(refseq_file_formatted).st_size > 0
        )
        # Dataframe
        df = pd.read_csv(refseq_file_formatted, sep="\t", header=None, skiprows=1)
        # Count number of transcript with version
        count_transcript_with_ver = df[df[4].str.contains(r"\.")].shape[0]
        # Count non canonical chromosome
        count_non_canonical_chr = sum("_" in key for key in set(df[0]))
        # Check
        assert len(df) > 0
        assert count_non_canonical_chr > 0
        assert count_transcript_with_ver > 0

        # Format refSeq by default with sort and header

        # Param
        refseq_file_formatted = f"{tmp_dir}/test2.bed"
        include_utr_5 = True
        include_utr_3 = True
        include_chrM = True
        include_non_canonical_chr = True
        include_non_coding_transcripts = True
        include_transcript_ver = True
        sort = True
        header = True
        header_first_line = True

        # Format
        databases_format_refseq(
            refseq_file=refseq_file_to_format,
            output_file=refseq_file_formatted,
            include_utr_3=include_utr_3,
            include_chrM=include_chrM,
            include_non_canonical_chr=include_non_canonical_chr,
            include_non_coding_transcripts=include_non_coding_transcripts,
            include_transcript_ver=include_transcript_ver,
            sort=sort,
            header=header,
            header_first_line=header_first_line,
        )

        # Check
        assert (
            os.path.exists(refseq_file_formatted)
            and os.stat(refseq_file_formatted).st_size > 0
        )
        assert (
            os.path.exists(f"{refseq_file_formatted}.hdr")
            and os.stat(f"{refseq_file_formatted}.hdr").st_size > 0
        )
        # Dataframe
        df = pd.read_csv(refseq_file_formatted, sep="\t", header=None, skiprows=1)
        # Count number of transcript with version
        count_transcript_with_ver = df[df[4].str.contains(r"\.")].shape[0]
        # Count non canonical chromosome
        count_non_canonical_chr = sum("_" in key for key in set(df[0]))
        # Check
        assert len(df) > 0
        assert count_non_canonical_chr > 0
        assert count_transcript_with_ver > 0

        # Format refSeq all parameters False

        # Param
        refseq_file_formatted = f"{tmp_dir}/test3.bed"
        include_utr_5 = False
        include_utr_3 = False
        include_chrM = False
        include_non_canonical_chr = False
        include_non_coding_transcripts = False
        include_transcript_ver = False
        sort = False
        header = False
        header_first_line = False

        # Format
        databases_format_refseq(
            refseq_file=refseq_file_to_format,
            output_file=refseq_file_formatted,
            include_utr_5=include_utr_5,
            include_utr_3=include_utr_3,
            include_chrM=include_chrM,
            include_non_canonical_chr=include_non_canonical_chr,
            include_non_coding_transcripts=include_non_coding_transcripts,
            include_transcript_ver=include_transcript_ver,
            sort=sort,
            header=header,
            header_first_line=header_first_line,
        )

        # Check
        assert (
            os.path.exists(refseq_file_formatted)
            and os.stat(refseq_file_formatted).st_size > 0
        )
        # Dataframe
        df = pd.read_csv(refseq_file_formatted, sep="\t", header=None, skiprows=0)
        # Count number of transcript with version
        count_transcript_with_ver = df[df[4].str.contains(r"\.")].shape[0]
        # Start end position of the first gene
        first_gene_start, first_gene_end = min(df[df[3] == df[3][0]][1]), max(
            df[df[3] == df[3][0]][2]
        )
        # Count non canonical chromosome
        count_non_canonical_chr = sum("_" in key for key in set(df[0]))
        # Check
        assert len(df) > 0
        assert count_non_canonical_chr == 0
        assert count_transcript_with_ver == 0

        # Format refSeq only UTR 5

        # Param
        refseq_file_formatted = f"{tmp_dir}/test4.bed"
        include_utr_5 = True
        include_utr_3 = False
        include_chrM = False
        include_non_canonical_chr = False
        include_non_coding_transcripts = False
        include_transcript_ver = False
        sort = False
        header = False
        header_first_line = True

        # Format
        databases_format_refseq(
            refseq_file=refseq_file_to_format,
            output_file=refseq_file_formatted,
            include_utr_5=include_utr_5,
            include_utr_3=include_utr_3,
            include_chrM=include_chrM,
            include_non_canonical_chr=include_non_canonical_chr,
            include_non_coding_transcripts=include_non_coding_transcripts,
            include_transcript_ver=include_transcript_ver,
            sort=sort,
            header=header,
            header_first_line=header_first_line,
        )

        # Check
        assert (
            os.path.exists(refseq_file_formatted)
            and os.stat(refseq_file_formatted).st_size > 0
        )
        # Dataframe
        df = pd.read_csv(refseq_file_formatted, sep="\t", header=None, skiprows=1)
        # Count number of transcript with version
        count_transcript_with_ver = df[df[4].str.contains(r"\.")].shape[0]
        # Start end position of the first gene
        first_gene_start_utr5, first_gene_end_utr5 = min(df[df[3] == df[3][0]][1]), max(
            df[df[3] == df[3][0]][2]
        )
        # Count non canonical chromosome
        count_non_canonical_chr = sum("_" in key for key in set(df[0]))
        # Check
        assert (
            first_gene_start_utr5 < first_gene_start
            and first_gene_end_utr5 == first_gene_end
        )
        assert len(df) > 0
        assert count_non_canonical_chr == 0
        assert count_transcript_with_ver == 0

        # Format refSeq only UTR 3

        # Param
        refseq_file_formatted = f"{tmp_dir}/test5.bed"
        include_utr_5 = False
        include_utr_3 = True
        include_chrM = False
        include_non_canonical_chr = False
        include_non_coding_transcripts = False
        include_transcript_ver = False
        sort = False
        header = False
        header_first_line = True

        # Format
        databases_format_refseq(
            refseq_file=refseq_file_to_format,
            output_file=refseq_file_formatted,
            include_utr_5=include_utr_5,
            include_utr_3=include_utr_3,
            include_chrM=include_chrM,
            include_non_canonical_chr=include_non_canonical_chr,
            include_non_coding_transcripts=include_non_coding_transcripts,
            include_transcript_ver=include_transcript_ver,
            sort=sort,
            header=header,
            header_first_line=header_first_line,
        )

        # Check
        assert (
            os.path.exists(refseq_file_formatted)
            and os.stat(refseq_file_formatted).st_size > 0
        )
        # Dataframe
        df = pd.read_csv(refseq_file_formatted, sep="\t", header=None, skiprows=1)
        # Count number of transcript with version
        count_transcript_with_ver = df[df[4].str.contains(r"\.")].shape[0]
        # Start end position of the first gene
        first_gene_start_utr3, first_gene_end_utr3 = min(df[df[3] == df[3][0]][1]), max(
            df[df[3] == df[3][0]][2]
        )
        # Count non canonical chromosome
        count_non_canonical_chr = sum("_" in key for key in set(df[0]))
        # Check
        assert (
            first_gene_start_utr3 == first_gene_start
            and first_gene_end_utr3 > first_gene_end
        )
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
        assemblies = ["hg19", "hg38"]

        # refSeq Folder
        refseq_folder = f"{tmp_dir}/test1"

        # downloaded_files expected
        downloaded_files_expected = {
            "hg19": ["ncbiRefSeq.txt", "ncbiRefSeqLink.txt"],
            "hg38": ["ncbiRefSeq.txt", "ncbiRefSeqLink.txt"],
        }

        # download
        downloaded_files = databases_download_refseq(
            assemblies=assemblies, refseq_folder=refseq_folder
        )
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
        downloaded_files = databases_download_refseq(
            assemblies=assemblies, refseq_folder=refseq_folder
        )
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
        downloaded_files_expected = {"hg19": refseq_files}

        # download
        downloaded_files = databases_download_refseq(
            assemblies=assemblies,
            refseq_folder=refseq_folder,
            refseq_files=refseq_files,
        )
        assert downloaded_files == downloaded_files_expected

        # check if each files exists
        for assembly in downloaded_files:
            for file in downloaded_files[assembly]:
                refseq_file = f"{refseq_folder}/{assembly}/{file}"
                assert os.path.exists(refseq_file) and os.stat(refseq_file).st_size > 0

        # Re-download
        downloaded_files = databases_download_refseq(
            assemblies=assemblies,
            refseq_folder=refseq_folder,
            refseq_files=refseq_files,
        )
        assert downloaded_files == downloaded_files_expected

        # check if each files exists
        for assembly in downloaded_files:
            for file in downloaded_files[assembly]:
                refseq_file = f"{refseq_folder}/{assembly}/{file}"
                assert os.path.exists(refseq_file) and os.stat(refseq_file).st_size > 0

        # And format

        # Param
        refseq_format_file = "ncbiRefSeq.txt"
        refseq_format_file_output = f"{tmp_dir}/test.bed.gz"
        include_utr_5 = True
        include_utr_3 = True
        include_chrM = True
        include_non_canonical_chr = True
        include_non_coding_transcripts = True
        include_transcript_ver = True

        # Format
        databases_download_refseq(
            assemblies=assemblies,
            refseq_folder=refseq_folder,
            refseq_files=refseq_files,
            refseq_format_file=refseq_format_file,
            refseq_format_file_output=refseq_format_file_output,
            include_utr_5=include_utr_5,
            include_utr_3=include_utr_3,
            include_chrM=include_chrM,
            include_non_canonical_chr=include_non_canonical_chr,
            include_non_coding_transcripts=include_non_coding_transcripts,
            include_transcript_ver=include_transcript_ver,
        )

        # Check
        assert (
            os.path.exists(refseq_format_file_output)
            and os.stat(refseq_format_file_output).st_size > 0
        )
