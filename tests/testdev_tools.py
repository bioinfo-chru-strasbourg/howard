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
from tempfile import TemporaryDirectory, tempdir
import duckdb
import re
import Bio.bgzf as bgzf
import gzip
import pytest

from howard.commons import *
from howard.tools.tsv_to_vcf import *


tests_folder = os.path.dirname(__file__)




def test_replace_dot():
    assert replace_dot('.') is np.nan
    assert replace_dot('2.3') == '2.3'
    assert replace_dot(4) == 4


def test_replace_char():
    # Test when input is a string
    assert replace_char('Hello ; World =') == 'Hello___World__'

    # Test when input contains no semicolon, space, or equal sign
    assert replace_char('HelloWorld') == 'HelloWorld'

    # Test when input is an integer
    assert replace_char(123) == 123

    # Test when input is a list of strings
    assert replace_char('Hello ; World =') == 'Hello___World__'


# def test_process_chunk():
#     # Create a test dataframe with some sample data
#     test_df = pd.DataFrame({'A': [1, 2, 3],
#                             'B': ['Hello ; World =', 'Python ; Programming =', 'Goodbye ; World ='],
#                             'C': ['Hello ; Python ; World =', 'Python ; Programming ; Language =', 'Goodbye ; World ; Python ='],
#                             'D': [4.2, 3.1, 2.7]})
    
#     # Define the indexes of the columns that contain the information we want to extract
#     info_idxs = [1, 2]
    
#     # Define the names of the columns that we want to extract from the file
#     annotations = ['B', 'C']
    
#     # Call the function and store the output
#     processed_df = process_chunk(test_df, info_idxs, annotations)
    
#     # Define the expected output
#     expected_output = pd.DataFrame({'B': ['Hello___World___', 'Python___Programming___', 'Goodbye___World___'],
#                                     'C': ['Hello___Python___World___', 'Python___Programming___Language___', 'Goodbye___World___Python___']})
    
#     # Check if the output matches the expected output
#     assert processed_df.equals(expected_output)


def test_read_line():
    chunk_data = pd.DataFrame({
        '#CHROM': ['1', '2', '3', '4'],
        'POS': [100, 200, 300, 400],
        'REF': ['A', 'C', 'G', 'T'],
        'ALT': ['C', 'T', 'A', 'G'],
        'ANN_1': ['foo', 'bar', '.', 'baz'],
        'ANN_2': ['.', 'qux', 'quux', 'corge'],
        'ANN_3': ['grault', '', 'garply', '.']
    })
    info_idxs = [4, 5, 6]
    annotations = ['ANN_1', 'ANN_2', 'ANN_3']
    result = read_line(chunk_data, info_idxs, annotations)
    assert result == [
        'chr1\t100\t.\tA\tC\t.\t.\tANN_1=foo;ANN_3=grault',
        'chr2\t200\t.\tC\tT\t.\t.\tANN_1=bar;ANN_2=qux',
        'chr3\t300\t.\tG\tA\t.\t.\tANN_2=quux;ANN_3=garply',
        'chr4\t400\t.\tT\tG\t.\t.\tANN_1=baz;ANN_2=corge'
    ]


def test_read_line_chr_indel():
    chunk_data = pd.DataFrame({
        '#CHROM': ['1', '2', 'chr3', '4'],
        'POS': [100, 200, 300, 400],
        'REF': ['A', '-', 'G', 'T'],
        'ALT': ['C', 'T', '-', 'G'],
        'ANN_1': ['foo', 'bar', '.', 'baz'],
        'ANN_2': ['.', 'qux', 'quux', 'corge'],
        'ANN_3': ['grault', '', 'garply', '.']
    })
    info_idxs = [4, 5, 6]
    annotations = ['ANN_1', 'ANN_2', 'ANN_3']
    result = read_line(chunk_data, info_idxs, annotations)
    assert result == [
        'chr1\t100\t.\tA\tC\t.\t.\tANN_1=foo;ANN_3=grault',
        'chr2\t199\t.\tN\tNT\t.\t.\tANN_1=bar;ANN_2=qux',
        'chr3\t299\t.\tNG\tN\t.\t.\tANN_2=quux;ANN_3=garply',
        'chr4\t400\t.\tT\tG\t.\t.\tANN_1=baz;ANN_2=corge'
    ]


def test_tsv_to_vcf():

    # input
    input_tsv = tests_folder + "/data/annotations/hg19_nci60.txt"
    output_vcf = tests_folder + "/tmp/hg19_nci60.txt.vcf.gz"
    genome = tests_folder + "/data/annotations/hg19.fa"
    database_name = "nci60"
    
    # check genome
    genome = find_genome(genome)

    #remove
    remove_if_exists([output_vcf])

    # transforme
    tsv_to_vcf(input_tsv, output_vcf, database_name=database_name, genome=genome)

    assert os.path.exists(output_vcf)


def test_tsv_to_vcf_with_correct_header():

    # input
    input_tsv = tests_folder + "/data/annotations/hg19_nci60.correct_header.txt"
    output_vcf = tests_folder + "/tmp/hg19_nci60.txt.vcf.gz"
    genome = tests_folder + "/data/annotations/hg19.fa"
    database_name = "nci60"

    # check genome
    genome = find_genome(genome)
    
    #remove
    remove_if_exists([output_vcf])

    # transforme
    tsv_to_vcf(input_tsv, output_vcf, database_name=database_name, genome=genome)

    assert os.path.exists(output_vcf)


def test_tsv_to_vcf_with_annovar_header():

    # input
    input_tsv = tests_folder + "/data/annotations/hg19_nci60.annovar_header.txt"
    output_vcf = tests_folder + "/tmp/hg19_nci60.txt.vcf.gz"
    genome = tests_folder + "/data/annotations/hg19.fa"
    database_name = "nci60"

    # check genome
    genome = find_genome(genome)
    
    #remove
    remove_if_exists([output_vcf])

    # transforme
    tsv_to_vcf(input_tsv, output_vcf, database_name=database_name, genome=genome)

    assert os.path.exists(output_vcf)



def test_tsv_to_vcf_with_strange_header():

    # input
    input_tsv = tests_folder + "/data/annotations/hg19_nci60.strange_header.txt"
    output_vcf = tests_folder + "/tmp/hg19_nci60.txt.vcf.gz"
    genome = tests_folder + "/data/annotations/hg19.fa"
    database_name = "nci60"

    # check genome
    genome = find_genome(genome)
    
    #remove
    remove_if_exists([output_vcf])

    # transforme
    tsv_to_vcf(input_tsv, output_vcf, database_name=database_name, genome=genome)

    assert os.path.exists(output_vcf)




def test_tsv_to_vcf_with_incorrect_header():

    # input
    input_tsv = tests_folder + "/data/annotations/hg19_nci60.incorrect_header.txt"
    output_vcf = tests_folder + "/tmp/hg19_nci60.txt.vcf.gz"
    genome = tests_folder + "/data/annotations/hg19.fa"
    database_name = "1nci60"

    # check genome
    genome = find_genome(genome)
    
    #remove
    remove_if_exists([output_vcf])

    # transforme
    #tsv_to_vcf(input_tsv, output_vcf, database_name=database_name, genome=genome)

    with pytest.raises(ValueError) as e:
        tsv_to_vcf(input_tsv, output_vcf, database_name=database_name, genome=genome)
    assert str(e.value) == f"Error in header"

    #assert not os.path.exists(output_vcf)


