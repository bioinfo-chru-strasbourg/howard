# -*- coding: utf-8 -*-
"""
Tests

Usage:
pytest tests/

Coverage:
coverage run -m pytest tests/test_tools_convert.py -x -v --log-cli-level=INFO --capture=tee-sys
coverage report --include=howard/* -m
"""

import argparse
import logging as log
import os
from tempfile import TemporaryDirectory

from howard.functions.commons import remove_if_exists
from howard.objects.variants import Variants
from howard.tools.convert import convert
from howard.tools.tools import arguments_dict

from test_needed import tests_data_folder, tests_folder


def test_convert():

    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # Init files
        input_vcf = tests_data_folder + "/example.vcf.gz"
        output_vcf = os.path.join(tmp_dir, "output_file.tsv")
        config = {"threads": 4}

        # for explode_infos in [True, False]:
        for explode_infos in [True]:

            # for explode_infos_prefix in ["", "INFO/", "CUSTOM_"]:
            for explode_infos_prefix in [""]:

                # for explode_infos_fields in ['ALL', 'DP,SIFT,AA']:
                for explode_infos_fields in ["DP,SIFT,AA"]:

                    # prepare arguments for the query function
                    args = argparse.Namespace(
                        input=input_vcf,
                        output=output_vcf,
                        config=config,
                        explode_infos=explode_infos,
                        explode_infos_prefix=explode_infos_prefix,
                        explode_infos_fields=explode_infos_fields,
                        include_header=True,
                        arguments_dict=arguments_dict,
                    )

                    # Remove if output file exists
                    remove_if_exists([output_vcf])

                    # Query
                    convert(args)

                    # Check output file exists
                    assert os.path.exists(output_vcf)

                    # DEVEL
                    variants = Variants(input=output_vcf, load=True)
                    res = variants.get_query_to_df("SELECT * FROM variants")
                    log.debug(res)

                    # read the contents of the actual output file
                    with open(output_vcf, "r") as f:
                        result_output_nb_lines = 0
                        result_output_nb_variants = 0
                        for line in f:
                            result_output_nb_lines += 1
                            if not line.startswith("#"):
                                result_output_nb_variants += 1

                    # Expected result
                    expected_result_nb_lines = 61
                    expected_result_nb_variants = 7

                    # Compare
                    assert result_output_nb_lines == expected_result_nb_lines
                    assert result_output_nb_variants == expected_result_nb_variants


def test_convert_pipeline():
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
        input_query = None
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
        convert(args)

        # Create object
        variants = Variants(conn=None, input=output_vcf, config=config, load=True)

        # Check annotation
        result = variants.get_query_to_df(
            "DESCRIBE SELECT * FROM variants"
        )
        log.debug(f"""result:\n{list(result.get("column_name"))}""")
        assert "DP" in list(result.get("column_name")), f"Expected 'DP' in columns, but got {list(result.get('column_name'))}"
        assert "SIFT" in list(result.get("column_name")), f"Expected 'SIFT' in columns, but got {list(result.get('column_name'))}"
        assert "AA" in list(result.get("column_name")), f"Expected 'AA' in columns, but got {list(result.get('column_name'))}"
        assert "CLNSIG" not in list(result.get("column_name")), f"Expected 'CLNSIG' not in columns, but got {list(result.get('column_name'))}"
        assert "ANN" not in list(result.get("column_name")), f"Expected 'ANN' not in columns, but got {list(result.get('column_name'))}"

        # Check annotation
        result = variants.get_query_to_df(
            "SELECT DP FROM variants WHERE DP IS NOT NULL"
        )
        assert len(result) == 2, f"Expected 2 variants with DP, but got {len(result)}"
        assert list(set(result.get("DP").tolist())) == list(set([50, 125])), f"Expected DP values [50, 125], but got {result.get('DP').tolist()}"
        log.debug(f"result:\n{result.to_string()}")

