# -*- coding: utf-8 -*-
"""
Tests

Usage:
pytest tests/

Coverage:
coverage run -m pytest tests/test_tools_calculation.py -x -v --log-cli-level=INFO --capture=tee-sys
coverage report --include=howard/* -m
"""

import argparse
import logging as log
import os
from tempfile import TemporaryDirectory

from howard.functions.commons import remove_if_exists
from howard.objects.variants import Variants
from howard.tools.calculation import calculation
from howard.tools.tools import arguments_dict

from test_needed import tests_data_folder, tests_folder


def test_calculation_tsv():

    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # Init files
        input_vcf = tests_data_folder + "/example.vcf.gz"
        output_vcf = os.path.join(tmp_dir, "output_file.tsv")
        config = {}
        calculations = "VARTYPE,NOMEN,TRIO"

        # prepare arguments for the query function
        args = argparse.Namespace(
            input=input_vcf,
            output=output_vcf,
            config=config,
            calculations=calculations,
            hgvs_field="hgvs",
            transcripts=None,
            show_calculations=False,
            trio_pedigree='{"father":"sample1", "mother":"sample2", "child":"sample3"}',
            calculation_config=None,
            arguments_dict=arguments_dict,
        )

        # Remove if output file exists
        remove_if_exists([output_vcf])

        # Query
        calculation(args)

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
            "SELECT INFO FROM variants WHERE INFO LIKE '%VARTYPE=%'"
        )
        assert len(result) == 7

        # Check annotation
        result = variants.get_query_to_df(
            "SELECT INFO FROM variants WHERE INFO LIKE '%trio=%'"
        )
        assert len(result) == 7

        # Check annotation
        result = variants.get_query_to_df(
            "SELECT INFO FROM variants WHERE INFO LIKE '%hgvs=%'"
        )
        assert len(result) == 0


def test_calculation_vcf():

    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # Init files
        input_vcf = tests_data_folder + "/example.vcf.gz"
        output_vcf = os.path.join(tmp_dir, "output_file.vcf")
        config = {}
        calculations = "VARTYPE,NOMEN,TRIO"

        # prepare arguments for the query function
        args = argparse.Namespace(
            input=input_vcf,
            output=output_vcf,
            config=config,
            calculations=calculations,
            hgvs_field="hgvs",
            transcripts=None,
            show_calculations=False,
            trio_pedigree='{"father":"sample1", "mother":"sample2", "child":"sample3"}',
            calculation_config=None,
            arguments_dict=arguments_dict,
        )

        # Remove if output file exists
        remove_if_exists([output_vcf])

        # Query
        calculation(args)

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
        expected_result_nb_lines = 77
        expected_result_nb_variants = 7

        # Compare
        assert result_output_nb_lines == expected_result_nb_lines
        assert result_output_nb_variants == expected_result_nb_variants

        # Create object
        variants = Variants(conn=None, input=output_vcf, config=config, load=True)

        # Check annotation
        result = variants.get_query_to_df(
            "SELECT INFO FROM variants WHERE INFO LIKE '%VARTYPE=%'"
        )
        assert len(result) == 7

        # Check annotation
        result = variants.get_query_to_df(
            "SELECT INFO FROM variants WHERE INFO LIKE '%trio=%'"
        )
        assert len(result) == 7

        # Check annotation
        result = variants.get_query_to_df(
            "SELECT INFO FROM variants WHERE INFO LIKE '%hgvs=%'"
        )
        assert len(result) == 0



def test_calculation_show_calculations():
    """
    Test processing with pipeline enabled and specific parameters.
    Including exporting transcripts to a TSV file.
    """

    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # Init files
        input_vcf = tests_data_folder + "/example.ann.vcf.gz"
        output_vcf = os.path.join(tmp_dir, "output_file.vcf")
        config = {}
        param = tests_folder + "/data/param.pipeline.json"
        annotations = None
        calculations = None
        prioritizations = None
        input_query = None
        show_calculations = True
        show_calculations_md = os.path.join(tmp_dir, "annotations.md")


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
            show_calculations=show_calculations,
            show_calculations_md=show_calculations_md
        )

        # Remove if output file exists
        remove_if_exists([output_vcf, show_calculations_md])

        # Query
        try:
            calculation(args)
        except SystemExit:
            log.info("SystemExit raised as expected due to show_calculations flag.")
            assert True

        assert os.path.exists(show_calculations_md), f"Expected markdown file {show_calculations_md} to be created."


def test_calculation_pipeline():
    """
    Test processing with pipeline enabled and specific parameters.
    Including exporting transcripts to a TSV file.
    """

    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # Init files
        input_vcf = tests_data_folder + "/example.ann.vcf.gz"
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
        calculation(args)

        # Create object
        variants = Variants(conn=None, input=output_vcf, config=config, load=True)

        # Check annotation
        result = variants.get_query_to_df(
            "SELECT INFO FROM variants WHERE INFO LIKE '%VARTYPE=%'"
        )
        assert len(result) == 7, f"Expected 7 variants with VARTYPE, but got {len(result)}"

        # Check annotation
        result = variants.get_query_to_df(
            "SELECT INFO FROM variants WHERE INFO LIKE '%snpeff_Annotation=%'"
        )
        assert len(result) == 7, f"Expected 7 variants with snpeff_Annotation, but got {len(result)}"

        # Check annotation
        result = variants.get_query_to_df(
            "SELECT INFO FROM variants WHERE INFO LIKE '%;NOMEN=%'"
        )
        assert len(result) == 7, f"Expected 7 variants with NOMEN, but got {len(result)}"
