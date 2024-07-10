"""
Tests

Usage:
pytest tests/

Coverage:
coverage run -m pytest tests/test_variants_annotations_splice.py -x -v --log-cli-level=DEBUG --capture=tee-sys
coverage report --include=howard/* -m
"""

import logging as log
from tempfile import TemporaryDirectory
from howard.functions.commons import remove_if_exists, identical
from howard.objects.variants import Variants
from test_needed import (
    tests_folder,
    tests_data_folder,
    tests_config,
    tests_databases_folder,
    tests_databases_release,
)

"""
Aim: test splice annotation
"""


def test_annotation_splice():
    """
    The function `test_annotation_splice` performs annotation and re-annotation of VCF files using the
    Splice tool with specified parameters and checks the results against expected output.
    """

    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # Init files
        input_vcf = tests_data_folder + "/example.vcf"
        output_vcf = f"{tmp_dir}/output.splice.vcf"
        output_reannotated_vcf = f"{tmp_dir}/output.splice.reannotated.vcf"
        expected_vcf = tests_data_folder + "/example.splice.vcf"

        # Construct param dict
        param = {
            "annotation": {
                "splice": {
                    "split_mode": "one",
                    "spliceai_distance": 500,
                    "spliceai_mask": 1,
                },
            }
        }
        config = tests_config.copy()
        config["tools"] = {
            "splice": {
                "docker": {
                    "image": "bioinfochrustrasbourg/splice:0.2.1",
                    "entrypoint": "/bin/bash",
                },
                "mount": {
                    f"{tests_databases_folder}/genomes/{tests_databases_release}": "ro",
                    tests_data_folder: "ro",
                },
                "tmp": "/tmp",
            },
            "docker": "docker",
        }

        # Create object
        variants = Variants(
            conn=None,
            input=input_vcf,
            output=output_vcf,
            param=param,
            load=True,
            config=config,
        )

        # Remove if output file exists
        remove_if_exists([output_vcf])

        # Annotation
        variants.annotation_splice()

        # Export output
        variants.export_output()

        # Ensure results
        assert identical([output_vcf, expected_vcf])

        # Check re-annotation

        # Create object
        variants_reannotated = Variants(
            conn=None,
            input=output_vcf,
            output=output_reannotated_vcf,
            param=param,
            load=True,
            config=config,
        )

        # Remove if output file exists
        remove_if_exists([output_reannotated_vcf])

        # Re-Annotation
        variants_reannotated.annotation_splice()

        # Export output
        variants_reannotated.export_output()

        # Ensure results
        assert identical([output_reannotated_vcf, expected_vcf])
