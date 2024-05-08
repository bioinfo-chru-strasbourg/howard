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
import duckdb
import re
import Bio.bgzf as bgzf
import gzip
import pytest
import pandas as pd
from pandas.testing import assert_frame_equal
from unittest.mock import patch

from howard.objects.variants import Variants
from howard.functions.utils import *
from test_needed import *


def test_json_perfect_exons_to_cdna_match():
    # Test without single
    ordered_exons = [(1, 10), (15, 20), (25, 30)]
    expected_result = [
        [1, 10, 1, 9, None],
        [15, 20, 10, 14, None],
        [25, 30, 15, 19, None],
    ]
    result = json_perfect_exons_to_cdna_match(ordered_exons)
    assert result == expected_result

    # Test with single
    ordered_exons_single = [(1, 10), (15, 20), (25, 30)]
    expected_result_single = [[1, 30, 1, 19, "M9 D5 M5 D5 M5"]]
    result_single = json_perfect_exons_to_cdna_match(ordered_exons_single, single=True)
    assert result_single == expected_result_single


def test_get_genomic_sequence():
    # Définissez les données d'entrée pour le test
    genome = {"chr1": "ATCGATCGATCG", "chr2": "GCTAGCTAGCTA"}
    chrom = "chr1"
    start = 2
    end = 5

    # Définissez le résultat attendu pour le test
    expected_result = "TCGA"

    # Appelez la fonction que vous souhaitez tester avec les données d'entrée
    result = get_genomic_sequence(genome, chrom, start, end)

    # Vérifiez si le résultat correspond au résultat attendu
    assert result == expected_result


def test_hgvs_justify_dup():

    genome = {"chr1": "ATCGATCGATCG", "chr2": "GCTAGCTAGCTA", "chr3": "TTTTGCTATTTT"}

    # Cas où ref et alt sont vides
    chrom = "chr1"
    offset = 10
    ref = ""
    alt = ""

    expected_result = ("chr1", 10, "", "", ">")

    result = hgvs_justify_dup(chrom, offset, ref, alt, genome)

    assert result == expected_result

    # Cas où ref et alt ne sont pas vides
    chrom = "chr2"
    offset = 5
    ref = "CTA"
    alt = "TCT"

    expected_result = ("chr2", 5, "CTA", "TCT", "delins")

    result = hgvs_justify_dup(chrom, offset, ref, alt, genome)

    assert result == expected_result

    # Cas où len(ref) > len(alt)
    chrom = "chr1"
    offset = 8
    ref = "ATCG"
    alt = ""

    expected_result = ("chr1", 8, "ATCG", "", "del")

    result = hgvs_justify_dup(chrom, offset, ref, alt, genome)

    assert result == expected_result

    # Cas dup previous
    chrom = "chr1"
    offset = 5
    ref = ""
    alt = "TTTT"

    expected_result = ("chr1", 5, "", "TTTT", "ins")

    result = hgvs_justify_dup(chrom, offset, ref, alt, genome)

    assert result == expected_result

    # Cas dup previous
    chrom = "chr1"
    offset = 8
    ref = ""
    alt = "GATC"

    expected_result = ("chr1", 4, "GATC", "GATCGATC", "dup")

    result = hgvs_justify_dup(chrom, offset, ref, alt, genome)

    assert result == expected_result

    # Cas dup next
    chrom = "chr3"
    offset = 5
    ref = ""
    alt = "GCTA"

    expected_result = ("chr3", 5, "GCTA", "GCTAGCTA", "dup")

    result = hgvs_justify_dup(chrom, offset, ref, alt, genome)

    assert result == expected_result


def test_hgvs_justify_indel():

    genome = {"chr1": "ATCGATCGATCG", "chr2": "GCTAGCTAGCTA"}

    # ref and alt empty
    chrom = "chr1"
    offset = 10
    ref = ""
    alt = ""
    strand = "+"

    expected_result = ("chr1", 10, "", "")

    result = hgvs_justify_indel(chrom, offset, ref, alt, strand, genome)

    assert result == expected_result

    # ref and alt not empty
    chrom = "chr2"
    offset = 5
    ref = "A"
    alt = "T"
    strand = "+"

    expected_result = ("chr2", 5, "A", "T")

    result = hgvs_justify_indel(chrom, offset, ref, alt, strand, genome)

    assert result == expected_result

    # ref and alt not empty
    chrom = "chr2"
    offset = 5
    ref = ""
    alt = "T"
    strand = "+"

    expected_result = ("chr2", 5, "", "T")

    result = hgvs_justify_indel(chrom, offset, ref, alt, strand, genome)

    assert result == expected_result

    # ref and alt not empty
    chrom = "chr2"
    offset = 5
    ref = "TTT"
    alt = ""
    strand = "+"

    expected_result = ("chr2", 5, "TTT", "")

    result = hgvs_justify_indel(chrom, offset, ref, alt, strand, genome)

    assert result == expected_result
