# -*- coding: utf-8 -*-
"""
Tests

Usage:
pytest tests/

Coverage:
coverage run -m pytest . -x -v --log-cli-level=INFO --capture=tee-sys
coverage report --include=howard/* -m
"""

import argparse
import os

from howard.functions.commons import remove_if_exists
from howard.objects.variants import Variants
from howard.tools.hgvs import hgvs
from howard.tools.tools import arguments_dict

from test_needed import tests_data_folder, tests_config


def test_hgvs_tsv():

    # Init files
    input_vcf = tests_data_folder + "/example.vcf.gz"
    output_vcf = "/tmp/output_file.tsv"
    config = {}

    genomes_folder = tests_config["folders"]["databases"]["genomes"]
    refseq_folder = tests_config["folders"]["databases"]["refseq"]

    # prepare arguments for the query function
    args = argparse.Namespace(
        input=input_vcf,
        output=output_vcf,
        hgvs_options="full_format",
        config=config,
        genomes_folder=genomes_folder,
        refseq_folder=refseq_folder,
        arguments_dict=arguments_dict,
    )

    # Remove if output file exists
    remove_if_exists([output_vcf])

    # Query
    hgvs(args)

    # Check output file exists
    assert os.path.exists(output_vcf)

    # read the contents of the actual output file
    with open(output_vcf, "r") as f:
        result_output_nb_lines = 0
        result_output_nb_variants = 0
        for line in f:
            result_output_nb_lines += 1
            if not line.startswith("#"):
                result_output_nb_variants += 1

    # Expected result
    expected_result_nb_lines = 8
    expected_result_nb_variants = 7

    # Compare
    assert result_output_nb_lines == expected_result_nb_lines
    assert result_output_nb_variants == expected_result_nb_variants

    # Create object
    variants = Variants(conn=None, input=output_vcf, config=config, load=True)

    # Check annotation
    result = variants.get_query_to_df(
        "SELECT INFO FROM variants WHERE INFO LIKE '%hgvs=%'"
    )
    assert len(result) == 7


def test_hgvs_vcf():

    # Init files
    input_vcf = tests_data_folder + "/example.vcf.gz"
    output_vcf = "/tmp/output_file.vcf"
    config = {}

    genomes_folder = tests_config["folders"]["databases"]["genomes"]
    refseq_folder = tests_config["folders"]["databases"]["refseq"]

    # prepare arguments for the query function
    args = argparse.Namespace(
        input=input_vcf,
        output=output_vcf,
        hgvs_options="full_format",
        config=config,
        genomes_folder=genomes_folder,
        refseq_folder=refseq_folder,
        arguments_dict=arguments_dict,
    )

    # Remove if output file exists
    remove_if_exists([output_vcf])

    # Query
    hgvs(args)

    # Check output file exists
    assert os.path.exists(output_vcf)

    # read the contents of the actual output file
    with open(output_vcf, "r") as f:
        result_output_nb_lines = 0
        result_output_nb_variants = 0
        for line in f:
            result_output_nb_lines += 1
            if not line.startswith("#"):
                result_output_nb_variants += 1

    # Expected result
    expected_result_nb_lines = 62
    expected_result_nb_variants = 7

    # Compare
    assert result_output_nb_lines == expected_result_nb_lines
    assert result_output_nb_variants == expected_result_nb_variants

    # Create object
    variants = Variants(conn=None, input=output_vcf, config=config, load=True)

    # Check annotation
    result = variants.get_query_to_df(
        "SELECT INFO FROM variants WHERE INFO LIKE '%hgvs=%'"
    )
    assert len(result) == 7
