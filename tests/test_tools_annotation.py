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
import logging as log
import os
from tempfile import TemporaryDirectory

from howard.functions.commons import remove_if_exists
from howard.objects.variants import Variants
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


def test_annotation_pipeline():
    """
    Test processing with pipeline enabled and specific parameters.
    Including exporting transcripts to a TSV file.
    """

    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # Init files
        input_vcf = tests_data_folder + "/example.vcf.gz"
        output_vcf = os.path.join(tmp_dir, "output_file.vcf")
        config = {}
        param = tests_folder + "/data/param.pipeline.json"
        annotations = None
        calculations = None
        prioritizations = None
        input_query = None
        output_transcripts_tsv = os.path.join(tmp_dir, "output.transcripts.tsv")

        # prepare arguments for the query function
        args = argparse.Namespace(
            input=input_vcf,
            output=output_vcf,
            config=config,
            param=param,
            annotations=annotations,
            calculations=calculations,
            prioritizations=prioritizations,
            query=input_query,
            explode_infos=False,
            explode_infos_prefix="",
            explode_infos_fields="*",
            include_header=False,
            arguments_dict=arguments_dict,
        )

        # Remove if output file exists
        remove_if_exists([output_vcf, output_transcripts_tsv])

        # Query
        annotation(args)

        # Create object
        variants = Variants(conn=None, input=output_vcf, config=config, load=True)

        # Check annotation
        result = variants.get_query_to_df(
            "SELECT INFO FROM variants WHERE INFO LIKE '%hgvs=%'"
        )
        assert len(result) == 7, f"Expected 7 variants with hgvs, but got {len(result)}"

        # Check annotation
        result = variants.get_query_to_df(
            "SELECT INFO FROM variants WHERE INFO LIKE '%SIFT_score=%'"
        )
        assert len(result) == 1, f"Expected 1 variant with SIFT_score, but got {len(result)}"
        
        # Check annotation
        result = variants.get_query_to_df(
            "SELECT INFO FROM variants WHERE INFO LIKE '%avsnp150=%'"
        )
        assert len(result) == 1, f"Expected 1 variant with avsnp150, but got {len(result)}"

        # Check annotation
        result = variants.get_query_to_df(
            "SELECT INFO FROM variants WHERE INFO LIKE '%cosmic70=%'"
        )
        assert len(result) == 1, f"Expected 1 variant with cosmic70, but got {len(result)}"