# -*- coding: utf-8 -*-
"""
Tests

Usage:
pytest tests/

Coverage:
coverage run -m pytest tests/test_tools_query.py -x -v --log-cli-level=INFO --capture=tee-sys
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
