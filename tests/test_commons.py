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

from howard.objects.variants import Variants
from howard.commons import *


tests_folder = os.path.dirname(__file__)


def test_get_file_compressed():
    # Test pour un fichier compressé .gz
    assert get_file_compressed("testfile.gz") == True
    
    # Test pour un fichier compressé .bcf
    assert get_file_compressed("testfile.bcf") == True
    
    # Test pour un fichier non compressé .vcf
    assert get_file_compressed("testfile.vcf") == False
    
    # Test pour un fichier non compressé .tsv
    assert get_file_compressed("testfile.tsv") == False
    
    # Test pour un fichier compressé .csv.gz
    assert get_file_compressed("testfile.csv.gz") == True
    
    # Test pour un fichier compressé .parquet
    assert get_file_compressed("testfile.parquet") == False
    
    # Test pour un fichier compressé .parquet.gz
    assert get_file_compressed("testfile.parquet.gz") == True


def test_get_file_format():
    # Test pour un fichier .vcf
    assert get_file_format("testfile.vcf") == "vcf"
    
    # Test pour un fichier .tsv.gz
    assert get_file_format("testfile.tsv.gz") == "tsv"
    
    # Test pour un fichier .csv
    assert get_file_format("testfile.csv") == "csv"
    
    # Test pour un fichier .bcf.gz
    assert get_file_format("testfile.bcf") == "bcf"
    
    # Test pour un fichier .bcf.gz
    assert get_file_format("testfile.bcf.gz") == "bcf"
    
    # Test pour un fichier .parquet
    assert get_file_format("testfile.parquet") == "parquet"
    
    # Test pour un fichier .txt.gz
    assert get_file_format("testfile.txt.gz") == "txt"


def test_remove_if_exists():

    # filename
    filename = "/tmp/output.test"

    # create filename
    fhandle = open(filename, 'a')
    try:
        os.utime(filename, None)
    finally:
        fhandle.close()
    created = os.path.exists(filename)

    # remove filename
    remove_if_exists([filename])

    # check delete
    deleted = not os.path.exists(filename)

    assert created and deleted


def test_remove_if_exists_complete():
    # Crée un fichier pour tester la fonction
    test_path = "test_dir"
    os.makedirs(test_path, exist_ok=True)
    test_file_path = os.path.join(test_path, "test_file.txt")
    with open(test_file_path, "w") as f:
        f.write("Test file content")
    
    # Teste la fonction avec un fichier qui existe
    remove_if_exists(test_file_path)
    assert not os.path.exists(test_file_path)
    
    # Teste la fonction avec un fichier qui n'existe pas
    remove_if_exists("nonexistent_file.txt")
    
    # Teste la fonction avec une liste de fichiers dont certains existent
    filepaths = [test_file_path, "nonexistent_file.txt"]
    remove_if_exists(filepaths)
    assert not os.path.exists(test_file_path)


def test_set_log_level():

    # define verbosity
    verbosity = "info"

    # set verbosity
    result_verbosity = set_log_level(verbosity)

    assert verbosity == result_verbosity


def test_set_log_level_error():

    # define verbosity
    verbosity = "not_a_level"

    # check verbosity
    with pytest.raises(ValueError) as e:
        set_log_level(verbosity)
    assert str(e.value) == "Unknown verbosity level:" + verbosity


def test_split_interval_either():

    start = 0
    end = 1000
    step = None
    ncuts = None
    with pytest.raises(ValueError) as e:
        split_interval(start, end, step=step, ncuts=ncuts)
    assert str(e.value) == "Either step or ncuts must be provided"


def test_split_interval_only():

    start = 0
    end = 1000
    step = 100
    ncuts = 4
    with pytest.raises(ValueError) as e:
        split_interval(start, end, step=step, ncuts=ncuts)
    assert str(e.value) == "Only one of step or ncuts must be provided"


def test_split_interval_step():

    start = 0
    end = 1000
    step = 100
    ncuts = None
    split = split_interval(start, end, step=step, ncuts=ncuts)
    split_expectd = [0, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000]
    
    assert split == split_expectd


def test_split_interval_ncuts():

    start = 0
    end = 1000
    step = None
    ncuts = 4
    split = split_interval(start, end, step=step, ncuts=ncuts)
    split_expectd = [0, 250, 500, 750, 1000]
    
    assert split == split_expectd


def test_merged_regions():

    # Define a list of genomic regions
    regions = [
        ('chr1', 100, 200),
        ('chr1', 150, 300),
        ('chr2', 500, 600),
        ('chr2', 550, 650),
        ('chr3', 800, 900)
    ]
    # Call the merge_regions function
    merged_regions = merge_regions(regions)
    # Define merged regions expected
    merged_regions_expected = [
        ('chr1', 100, 300),
        ('chr2', 500, 650),
        ('chr3', 800, 900)
    ]
    
    assert merged_regions == merged_regions_expected


def test_create_where_clause():

    # Define a list of merged regions
    merged_regions = [
        ('chr1', 100, 300),
        ('chr2', 500, 650),
        ('chr3', 800, 900),
        ('chr3', 1000, 1200)
    ]

    # define table
    table = "variants"

    # Call the create_where_clause function
    where_clause = create_where_clause(merged_regions, table=table)

    # Create expected where clause
    where_clause_expected = """ ( variants."#CHROM" = 'chr1' AND (   (variants.POS >= 100 AND variants.POS <= 300)  ) )   OR  ( variants."#CHROM" = 'chr2' AND (   (variants.POS >= 500 AND variants.POS <= 650)  ) )   OR  ( variants."#CHROM" = 'chr3' AND (   (variants.POS >= 800 AND variants.POS <= 900)   OR  (variants.POS >= 1000 AND variants.POS <= 1200)  ) ) """
    
    assert where_clause.strip() == where_clause_expected.strip()


def test_command():

    # Creation of a command
    cmd = "echo 'Command'"

    # Execute command
    command_output = command(cmd)

    # Expected result
    results_expected = "Command"

    assert command_output == results_expected


def test_run_parallel_commands():

    # Creation of a list of command
    commands = [
        "echo 'Command 1'",
        "echo 'Command 2'",
        "echo 'Command 3'"
    ]

    # Execute commands in parallele
    results = run_parallel_commands(commands, 2)

    # Expected results
    results_expected = ['Command 1', 'Command 2', 'Command 3']

    assert results == results_expected


def test_run_parallel_functions():

    # List of functions (examples)
    functions = [example_function(1, "hello"), example_function(2, "world")]
    threads = 2

    # Launch functions
    results = run_parallel_functions(functions, threads)

    # Number of result expected
    expected_output_length = 2

    assert len(results) == expected_output_length


def test_example_function():

    # Launch functions
    result = example_function(1, "hello")
                               
    # result expected
    expected_result = [1, "hello"]

    assert result == expected_result


def test_get_index():

    # Test with a list of values
    values = ['a', 'b', 'c', 'd']
    assert get_index('a', values) == 0
    assert get_index('d', values) == 3
    assert get_index('e', values) == -1

    # Test with an empty list
    assert get_index('a', []) == -1
    assert get_index('', []) == -1

    # Test with a None element
    assert get_index(None, values) == -1


def test_find_nomen_full():
    # Test case 1
    print("test1")
    hgvs = "NM_001637.3:c.1582G>T"
    transcripts = ["NM_001637.3"]
    expected_output = {
        "NOMEN": "NM_001637:c.1582G>T",
        "CNOMEN": "c.1582G>T",
        "RNOMEN": None,
        "NNOMEN": None,
        "PNOMEN": None,
        "TVNOMEN": "NM_001637.3",
        "TNOMEN": "NM_001637",
        "VNOMEN": "3",
        "ENOMEN": None,
        "GNOMEN": None,
    }
    assert find_nomen(hgvs, transcripts=transcripts) == expected_output

    # Test case 2
    print("test2")
    hgvs = "NM_001637.3:c.1582G>T,NM_001637.3:c.1583G>T"
    transcripts = ["NM_001637.3"]
    expected_output = {
        "NOMEN": "NM_001637:c.1582G>T",
        "CNOMEN": "c.1582G>T",
        "RNOMEN": None,
        "NNOMEN": None,
        "PNOMEN": None,
        "TVNOMEN": "NM_001637.3",
        "TNOMEN": "NM_001637",
        "VNOMEN": "3",
        "ENOMEN": None,
        "GNOMEN": None,
    }
    assert find_nomen(hgvs, transcripts=transcripts) == expected_output

    # Test case 3
    print("test3")
    hgvs = "NM_001637.3:c.1582G>T,NM_001637.3:c.1583G>T,NM_001637.2:c.1582G>T:p.G12D"
    transcripts = ["NM_001637.2", "NM_001637.3"]
    expected_output = {
        "NOMEN": "NM_001637:c.1582G>T:p.G12D",
        "CNOMEN": "c.1582G>T",
        "RNOMEN": None,
        "NNOMEN": None,
        "PNOMEN": "p.G12D",
        "TVNOMEN": "NM_001637.2",
        "TNOMEN": "NM_001637",
        "VNOMEN": "2",
        "ENOMEN": None,
        "GNOMEN": None,
    }
    assert find_nomen(hgvs, transcripts=transcripts) == expected_output

    # Test case 4
    print("test3")
    hgvs = "Gene1:exon12:n.1582G>T:NR_001637.3"
    transcripts = []
    expected_output = {
        "NOMEN": "Gene1:NR_001637:exon12:n.1582G>T",
        "CNOMEN": None,
        "RNOMEN": None,
        "NNOMEN": "n.1582G>T",
        "PNOMEN": None,
        "TVNOMEN": "NR_001637.3",
        "TNOMEN": "NR_001637",
        "VNOMEN": "3",
        "ENOMEN": "exon12",
        "GNOMEN": "Gene1",
    }
    assert find_nomen(hgvs, transcripts=transcripts) == expected_output

    # Test case 5
    print("test3")
    hgvs = "Gene1:exon12:r.1582G>T"
    transcripts = []
    expected_output = {
        "NOMEN": "Gene1:exon12:r.1582G>T",
        "CNOMEN": None,
        "RNOMEN": "r.1582G>T",
        "NNOMEN": None,
        "PNOMEN": None,
        "TVNOMEN": None,
        "TNOMEN": None,
        "VNOMEN": None,
        "ENOMEN": "exon12",
        "GNOMEN": "Gene1",
    }
    assert find_nomen(hgvs, transcripts=transcripts) == expected_output


def test_get_gzip():

    command_gzip_expected = "bgzip -c  --threads=2 --compress-level=5"

    command_gzip = get_bgzip(threads=2, level=5)
    print(command_gzip)

    assert command_gzip.strip() == command_gzip_expected.strip()
          

def test_find():
    # Crée un fichier pour tester la fonction
    test_path = "test_dir"
    os.makedirs(test_path, exist_ok=True)
    test_file_path = os.path.join(test_path, "test_file.txt")
    with open(test_file_path, "w") as f:
        f.write("Test file content")
    
    # Teste la fonction avec un nom de fichier qui n'existe pas
    assert find("nonexistent_file.txt", test_path) == ""
    
    # Teste la fonction avec un nom de fichier qui existe dans le dossier de base
    assert find("test_file.txt", test_path) == test_file_path
    
    # Teste la fonction avec un nom de fichier qui existe dans un sous-dossier
    sub_dir_path = os.path.join(test_path, "sub_dir")
    os.makedirs(sub_dir_path, exist_ok=True)
    sub_dir_file_path = os.path.join(sub_dir_path, "sub_dir_file.txt")
    with open(sub_dir_file_path, "w") as f:
        f.write("Sub-directory file content")
    assert find("sub_dir_file.txt", test_path) == sub_dir_file_path
    
    # Nettoie le dossier de test
    os.remove(test_file_path)
    os.remove(sub_dir_file_path)
    os.rmdir(sub_dir_path)
    os.rmdir(test_path)


def test_find_all():
    # Create a temporary directory structure
    with tempfile.TemporaryDirectory() as tmpdir:
        # Create some files with the name 'test_file' in different directories
        open(os.path.join(tmpdir, 'test_file'), 'a').close()
        os.makedirs(os.path.join(tmpdir, 'subdir'))
        open(os.path.join(tmpdir, 'subdir', 'test_file'), 'a').close()

        # Test that find_all returns the correct paths
        assert find_all('test_file', tmpdir) == [
            os.path.join(tmpdir, 'test_file'),
            os.path.join(tmpdir, 'subdir', 'test_file')
        ]


def test_find_genome():

    # Default genome filenaùe
    genome_filename = "hg19.fa"

    # Create existing genome
    tmp_genome = NamedTemporaryFile(prefix="genome.", suffix=".fa", dir="/tmp")
    tmp_genome_name = tmp_genome.name
    # call the function to find the genome file
    genome_path_found = find_genome(genome_path = tmp_genome_name, genome = genome_filename)
    # check if the genome file was found
    assert os.path.exists(genome_path_found) and genome_path_found == tmp_genome_name


    # Either genome in the system or not

    # create a temporary directory
    with tempfile.TemporaryDirectory() as tmpdir:
        # specify a non-existent path for the genome file
        genome_path_nonexistent = os.path.join(tmpdir, 'nonexistent_genome.fa')
        # call the function to find the genome file
        error = None
        try:
            genome_path_found = find_genome(genome_path = genome_path_nonexistent, genome = genome_filename)
        except:
            with pytest.raises(ValueError) as e:
                genome_path_found = find_genome(genome_path = genome_path_nonexistent, genome = genome_filename)
                error = str(e.value)
        
        # check if the genome file was found

        assert (os.path.exists(genome_path_found) and genome_path_found != genome_path_nonexistent) or error == f"Genome failed: no genome '{genome_filename}'"

    
