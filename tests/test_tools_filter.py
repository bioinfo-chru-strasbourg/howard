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
import pytest  # type: ignore

from howard.objects.variants import Variants
from howard.functions.commons import remove_if_exists
from howard.tools.tools import arguments_dict
from howard.tools.filter import filter as vcf_filter

from test_needed import tests_folder, tests_data_folder


@pytest.mark.parametrize(
    "filter, samples, expected_results",
    [
        # Filter 1
        (
            "POS < 100000",
            None,
            {
                "nb_lines": 57,
                "nb_variants": 3,
                "samples": [
                    "sample1",
                    "sample2",
                    "sample3",
                    "sample4",
                ],
            },
        ),
        # Filter with INFOS
        (
            "INFOS.CLNSIG LIKE 'pathogenic'",
            None,
            {
                "nb_lines": 55,
                "nb_variants": 1,
                "samples": [
                    "sample1",
                    "sample2",
                    "sample3",
                    "sample4",
                ],
            },
        ),
        # Filter with SAMPLES
        (
            "SAMPLES.sample2.GT != './.'",
            "sample1,sample2",
            {
                "nb_lines": 57,
                "nb_variants": 3,
                "samples": [
                    "sample1",
                    "sample2",
                ],
            },
        ),
    ],
)
def test_filter(filter, samples, expected_results):
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
            filter=filter,
            samples=samples,
            include_header=True,
            arguments_dict=arguments_dict,
        )

        # Remove if output file exists
        remove_if_exists([output_vcf])

        # Filter
        vcf_filter(args)

        # read the contents of the actual output file
        with open(output_vcf, "r") as f:
            result_output_nb_lines = 0
            result_output_nb_variants = 0
            result_lines = []
            for line in f:
                if not result_output_nb_lines:
                    log.debug(line)
                result_output_nb_lines += 1
                if not line.startswith("#"):
                    result_output_nb_variants += 1
                    result_lines.append(line.strip())

        # Expected result
        expected_result_nb_lines = expected_results.get("nb_lines", None)
        expected_result_nb_variants = expected_results.get("nb_variants", None)
        expected_result_samples = expected_results.get("samples", None)

        # Compare
        assert result_output_nb_lines == expected_result_nb_lines
        assert result_output_nb_variants == expected_result_nb_variants

        # Variants
        variants = Variants(input=output_vcf, load=True)
        assert variants.get_header_sample_list() == expected_result_samples
