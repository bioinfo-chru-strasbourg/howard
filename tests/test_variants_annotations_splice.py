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
from howard.functions.commons import identical
from howard.objects.variants import Variants
from test_needed import (
    tests_data_folder,
    tests_folder,
    tests_databases_folder,
    tests_databases_release,
)
import os

"""
Aim: test splice annotation
"""


PARAM = {
    "annotation": {
        "splice": {
            "options": {
                "split_mode": "one",
                "spliceai_distance": 500,
                "spliceai_mask": 1,
                "workdir": "/work",
                "condacachedir": "/work",
                "threads": 8,
            }
        }
    }
}
CONFIG = {
    "folders": {
        "databases": {
            "genomes": f"{tests_databases_folder}/genomes/{tests_databases_release}"
        }
    },
    "tools": {
        "splice": {
            "docker": {
                "image": "bioinfochrustrasbourg/splice:0.2.3",
                "entrypoint": "/bin/bash",
            }
        }
    },
}

CONFIG_DOCKER = {
    "folders": {
        "databases": {
            "genomes": f"{tests_databases_folder}/genomes/{tests_databases_release}",
        },
    },
    "tools": {
        "splice": {
            "docker": {
                "image": "bioinfochrustrasbourg/splice:0.2.3",
                "entrypoint": "/bin/bash",
                "config": {"automount": True, "notremove": True, "tmp": False},
            }
        }
    },
}


def test_annotation_splice():
    """
    The function `test_annotation_splice` performs annotation and re-annotation of VCF files using the
    Splice tool with specified parameters and checks the results against expected output.
    """
    if not os.path.exists("/.dockerenv"):
        config = CONFIG
    else:
        config = CONFIG_DOCKER
    param = PARAM
    with TemporaryDirectory(dir=tests_folder) as tmp_dir:
        param["annotation"]["splice"]["options"]["output_folder"] = tmp_dir
        # Init files
        input_vcf = tests_data_folder + "/example.vcf"
        output_vcf = f"{tmp_dir}/output.splice.vcf"
        expected_vcf = tests_data_folder + "/example.splice.vcf"
        # Construct param dict
        # Create object
        variants = Variants(
            conn=None,
            input=input_vcf,
            output=output_vcf,
            param=param,
            load=True,
            config=config,
        )
        # Annotation
        variants.annotation_splice()
        # Export output
        variants.export_output()
        # Ensure results
        assert identical([output_vcf, expected_vcf])
