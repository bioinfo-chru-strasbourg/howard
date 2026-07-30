# -*- coding: utf-8 -*-
"""
Tests

Usage:
pytest tests/

Coverage:
coverage run -m pytest tests/test_tools_sort.py -x -v --log-cli-level=INFO --capture=tee-sys
coverage report --include=howard/* -m
"""

import logging as log
import os
from tempfile import TemporaryDirectory
import argparse

from howard.objects.variants import Variants
from howard.functions.commons import remove_if_exists
from howard.tools.tools import arguments_dict
from howard.tools.sort import sort as vcf_sort

from test_needed import tests_folder, tests_data_folder


def test_filter():
    """
    The `test_filter` function tests filter of variants and exporting the output in correct
    format using pyVCF.
    """

    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # Init files
        input_vcf = tests_data_folder + "/example.vcf.gz"
        output_vcf = os.path.join(tmp_dir, "output_file.tsv")

        # prepare arguments for the query function
        args = argparse.Namespace(
            input=input_vcf,
            output=output_vcf,
            include_header=True,
            arguments_dict=arguments_dict,
        )

        # Remove if output file exists
        remove_if_exists([output_vcf])

        # Filter
        vcf_sort(args)

        # Variants
        variants = Variants(input=output_vcf, load=True)
        assert list(variants.get_header().contigs.keys()) == [
            "1",
            "chr1",
            "2",
            "3",
            "4",
            "5",
            "6",
            "7",
            "chr7",
            "8",
            "9",
            "10",
            "11",
            "12",
            "13",
            "14",
            "15",
            "16",
            "17",
            "18",
            "19",
            "20",
            "21",
            "22",
            "X",
            "Y",
            "M",
        ]


def test_sort_pipeline():
    """
    Test processing with pipeline enabled and specific parameters.
    Including exporting transcripts to a TSV file.
    """

    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # Init files
        input_vcf = tests_data_folder + "/example.ann.vcf.gz"
        output_vcf = os.path.join(tmp_dir, "output_file.tsv")
        config = {}
        param = tests_folder + "/data/param.pipeline.json"
        annotations = None
        calculations = None
        prioritizations = None
        input_query = "SELECT * FROM variants WHERE REF = 'A' AND POS < 100000"
        explode_infos_fields = "DP,SIFT,AA"

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
            explode_infos=True,
            explode_infos_prefix="",
            explode_infos_fields=explode_infos_fields,
            include_header=False,
            arguments_dict=arguments_dict,
        )

        # Remove if output file exists
        remove_if_exists([output_vcf])

        # Query
        vcf_sort(args)

        # Create object
        variants = Variants(conn=None, input=output_vcf, config=config, load=True)

        # Check annotation
        result = variants.get_query_to_df(
            "SELECT * FROM variants"
        )
        assert len(result) == 7, f"Expected 7 variants with sorting, but got {len(result)}"

        # Chromosomes list
        assert list(variants.get_header().contigs.keys()) == [
            "1",
            "chr1",
            "2",
            "3",
            "4",
            "5",
            "6",
            "7",
            "chr7",
            "8",
            "9",
            "10",
            "11",
            "12",
            "13",
            "14",
            "15",
            "16",
            "17",
            "18",
            "19",
            "20",
            "21",
            "22",
            "X",
            "Y",
            "M",
        ]

