# -*- coding: utf-8 -*-
"""
Tests

Usage:
pytest tests/

Coverage:
coverage run -m pytest tests/test_tools_filter.py -x -v --log-cli-level=INFO --capture=tee-sys
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



def test_filter_pipeline():
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
        #explode_infos_fields = "DP,SIFT,AA"
        filter="INFOS.CLNSIG LIKE 'pathogenic'"
        samples="sample1,sample2"

        # prepare arguments for the filter function
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
            #explode_infos_fields=explode_infos_fields,
            filter=filter,
            samples=samples,
            include_header=False,
            arguments_dict=arguments_dict,
        )

        # Remove if output file exists
        remove_if_exists([output_vcf])

        # Filter
        vcf_filter(args)

        # Create object
        variants = Variants(conn=None, input=output_vcf, config=config, load=True)

        # DEVEL
        result = variants.get_query_to_df("SELECT * FROM variants")
        log.debug(f"result:\n{result.to_string()}")

        # Check annotation
        result = variants.get_query_to_df(
            "SELECT * FROM variants"
        )
        assert len(result) == 1, f"Expected 1 variant with filter '{filter}', but got {len(result)}"

        # Check annotation
        result = variants.get_query_to_df(
            "DESCRIBE SELECT * FROM variants"
        )
        log.debug(f"""result:\n{list(result.get("column_name"))}""")
        assert "sample1" in list(result.get("column_name")), f"Expected 'sample1' in columns, but got {list(result.get('column_name'))}"
        assert "sample2" in list(result.get("column_name")), f"Expected 'sample2' in columns, but got {list(result.get('column_name'))}"
        assert "sample3" not in list(result.get("column_name")), f"Expected 'sample3' not in columns, but got {list(result.get('column_name'))}"
        assert "sample4" not in list(result.get("column_name")), f"Expected 'sample4' not in columns, but got {list(result.get('column_name'))}"
