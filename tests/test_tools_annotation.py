# -*- coding: utf-8 -*-
"""
Tests

Usage:
pytest tests/

Coverage:
coverage run -m pytest tests/test_tools_annotation.py -x -v --log-cli-level=INFO --capture=tee-sys
coverage report --include=howard/* -m
"""

import argparse
import os
from tempfile import TemporaryDirectory

from howard.functions.commons import remove_if_exists
from howard.tools.annotation import annotation
from howard.tools.tools import arguments_dict
from test_needed import tests_data_folder, database_files, tests_folder


def test_annotation_tsv_update():

    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # Init files
        input_vcf = tests_data_folder + "/example.vcf.gz"
        output_vcf = os.path.join(tmp_dir, "output_file.tsv")
        config = {}
        annotations = database_files.get("parquet")

        # prepare arguments for the query function
        args = argparse.Namespace(
            input=input_vcf,
            output=output_vcf,
            config=config,
            annotations=annotations,
            annotations_update=True,
            annotations_append=False,
            arguments_dict=arguments_dict,
        )

        # Remove if output file exists
        remove_if_exists([output_vcf])

        # Query
        annotation(args)

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


def test_annotation_tsv_append():

    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # Init files
        input_vcf = tests_data_folder + "/example.vcf.gz"
        output_vcf = os.path.join(tmp_dir, "output_file.tsv")
        config = {}
        annotations = database_files.get("parquet")

        # prepare arguments for the query function
        args = argparse.Namespace(
            input=input_vcf,
            output=output_vcf,
            config=config,
            annotations=annotations,
            annotations_update=False,
            annotations_append=True,
            arguments_dict=arguments_dict,
        )

        # Remove if output file exists
        remove_if_exists([output_vcf])

        # Query
        annotation(args)

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


def test_annotation_tsv_update_append():

    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # Init files
        input_vcf = tests_data_folder + "/example.vcf.gz"
        output_vcf = os.path.join(tmp_dir, "output_file.tsv")
        config = {}
        annotations = database_files.get("parquet")

        # prepare arguments for the query function
        args = argparse.Namespace(
            input=input_vcf,
            output=output_vcf,
            config=config,
            annotations=annotations,
            annotations_update=True,
            annotations_append=True,
            arguments_dict=arguments_dict,
        )

        # Remove if output file exists
        remove_if_exists([output_vcf])

        # Query
        annotation(args)

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


def test_annotation_tsv():

    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # Init files
        input_vcf = tests_data_folder + "/example.vcf.gz"
        output_vcf = os.path.join(tmp_dir, "output_file.tsv")
        config = {}
        annotations = database_files.get("parquet")

        # prepare arguments for the query function
        args = argparse.Namespace(
            input=input_vcf,
            output=output_vcf,
            config=config,
            annotations=annotations,
            annotations_update=False,
            annotations_append=False,
            arguments_dict=arguments_dict,
        )

        # Remove if output file exists
        remove_if_exists([output_vcf])

        # Query
        annotation(args)

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


def test_annotation_vcf():

    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # Init files
        input_vcf = tests_data_folder + "/example.vcf.gz"
        output_vcf = os.path.join(tmp_dir, "output_file.vcf")
        config = {}
        annotations = database_files.get("parquet")

        # prepare arguments for the query function
        args = argparse.Namespace(
            input=input_vcf,
            output=output_vcf,
            config=config,
            annotations=annotations,
            annotations_update=False,
            annotations_append=False,
            arguments_dict=arguments_dict,
        )

        # Remove if output file exists
        remove_if_exists([output_vcf])

        # Query
        annotation(args)

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


def test_annotation_vcf_fast():

    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # Init files
        input_vcf = tests_data_folder + "/example.vcf.gz"
        output_vcf = os.path.join(tmp_dir, "output_file.vcf")
        config = {}
        annotations = database_files.get("parquet")

        # prepare arguments for the query function
        args = argparse.Namespace(
            input=input_vcf,
            output=output_vcf,
            config=config,
            annotations=annotations,
            annotations_update=False,
            annotations_append=False,
            fast=True,
            arguments_dict=arguments_dict,
        )

        # Remove if output file exists
        remove_if_exists([output_vcf])

        # Query
        annotation(args)

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
