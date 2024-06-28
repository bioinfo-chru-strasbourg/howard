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
    with TemporaryDirectory(dir=tests_folder) as tmp_dir:
        input_vcf = tests_data_folder + "/example.vcf"
        output_vcf = f"{tmp_dir}/output.splice.vcf"
        expected_vcf = tests_data_folder + "/example.splice.vcf"

        # Construct param dict
        param = {
            "annotation": {
                "splice": {
                    "options": {
                        "genome": "hg19",
                        "threads": 2,
                        "split_mode": "one",
                        "spliceai_distance": 500,
                        "spliceai_mask": 1,
                    }
                },
            }
        }
        config = tests_config.copy()
        config["tools"] = {
            "splice": {
                "docker": {
                    "image": "jblamouche/splice:0.2.1",
                    "entrypoint": "/bin/bash",
                },
                "mount": {
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
        identical([output_vcf, expected_vcf])
