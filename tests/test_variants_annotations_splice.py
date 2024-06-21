import logging as log
from tempfile import TemporaryDirectory
from howard.functions.commons import *
from howard.objects.variants import Variants
from howard.functions.databases import *
from test_needed import *


def test_annotation_splice():
    with TemporaryDirectory(dir=tests_folder) as tmp_dir:
        input_vcf = tests_data_folder + "/example.vcf"
        output_vcf = f"{tmp_dir}/output.splice.vcf"

        # Construct param dict
        param = {
            "annotation": {
                "splice": {"options": {"genome": "hg19", "threads": 2}},
            }
        }
        config = tests_config.copy()
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
