# -*- coding: utf-8 -*-
"""
Tests

Usage:
pytest tests/

Coverage:
coverage run -m pytest tests/test_variants_annotations_bcftools.py -x -v --log-cli-level=INFO --capture=tee-sys
coverage report --include=howard/* -m
"""

import copy
import json
import os
from tempfile import TemporaryDirectory
import pytest  # type: ignore
import vcf  # type: ignore
import logging as log
import yaml  # type: ignore

from howard.objects.variants import Variants
from howard.functions.commons import DEFAULT_DATABASE_FOLDER, full_path, remove_if_exists

from test_needed import (
    tests_data_folder,
    tests_folder,
    tests_annotations_folder,
    tests_config
)


@pytest.mark.parametrize(
    "parameters, mode, tool, config_vep",
    [
        # Test with param mode (into param.json file/dict)
        # None
        (
            None,
            "param",
            "vep",
            None,
        ),
        # Empty string
        (
            "",
            "param",
            "vep",
            None,
        ),
        # Empty list
        (
            [],
            "param",
            "vep",
            None,
        ),
        # Symbol with dashes --, only for uniq option
        (
            "--symbol",
            "param",
            "vep",
            None,
        ),
        # Everything with dashes -- in list
        (
            ["--symbol", "--hgvs"],
            "param",
            "vep",
            None,
        ),
        # With config with vep_online
        (
            None,
            "param",
            "vep_online",
            os.path.join(tests_folder, "config", "config.vep.json"),
        ),
        # Test with parameters mode (with --annotations command line option, i.e. "annotations" section)
        # None
        (
            None,
            "annotations",
            "vep",
            None,
        ),
        # Empty string
        (
            "",
            "annotations",
            "vep",
            None,
        ),
        # symbol without dashes --
        (
            "symbol",
            "annotations",
            "vep",
            None,
        ),
        # symbol and hgvs without dashes --
        (
            "symbol:hgvs",
            "annotations",
            "vep",
            None,
        ),
        # Symbol with dashes --
        (
            "--symbol",
            "annotations",
            "vep",
            None,
        ),
        # Symbol and hgvs with dashes --
        (
            "--symbol --hgvs --refseq",
            "annotations",
            "vep",
            None,
        ),
    ],
)
def test_annotation_vep_annotations(parameters, mode, tool, config_vep):
    """
    Test VEP annotation with different parameter sets.
    """

    with TemporaryDirectory(dir=tests_folder) as tmp_dir:
    
        # Init files
        input_vcf = tests_data_folder + "/example.vcf.gz"
        output_vcf = f"{tmp_dir}/output.vcf.gz"

        # Copy config
        if config_vep is None:
            tests_config_vep = tests_config.copy()
        else:
            log.debug("Using VEP config: %s", config_vep)
            tests_config_vep = yaml.safe_load(open(config_vep))

        # Number of threads
        tests_config_vep["threads"] = 2

        # Memory
        tests_config_vep["memory"] = "4G"

        # Construct param dict
        if mode == "param":
            param = {
                "annotation": {
                    "vep": {
                        "tool": tool,
                        "parameters": parameters
                    }
                }
            }
        elif mode == "annotations":
            param = {
                "annotations": ":".join(["vep", f"{parameters}"]) if parameters else "vep"
            }
        else:
            pytest.fail("Invalid mode {}".format(mode))

        # Create object
        variants = Variants(
            conn=None,
            input=input_vcf,
            output=output_vcf,
            config=tests_config_vep,
            param=param,
            load=True,
        )

        # Remove if output file exists
        remove_if_exists([output_vcf])

        # Annotation
        variants.annotation()

        # query annotated variant
        result = variants.get_query_to_df(
            """ SELECT * FROM variants WHERE "INFO" LIKE '%CSQ%' """
        )
        assert len(result) == 7, f"Expected 7 annotated variants, got {len(result)}, with tests parameters {parameters}, mode {mode}, tool {tool}, config {config_vep}"

        # Check if VCF is in correct format with pyVCF
        variants.export_output()
        try:
            vcf.Reader(filename=output_vcf)
        except:
            assert False

