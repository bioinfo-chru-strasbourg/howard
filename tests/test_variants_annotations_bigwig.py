# -*- coding: utf-8 -*-
"""
Tests

Usage:
pytest tests/

Coverage:
coverage run -m pytest tests/test_variants_annotations_snpsift.py -x -v --log-cli-level=INFO --capture=tee-sys
coverage report --include=howard/* -m 
"""

import logging as log
import os
from tempfile import TemporaryDirectory
import vcf  # type: ignore

from howard.objects.variants import Variants
from howard.functions.commons import remove_if_exists

from test_needed import tests_folder, tests_data_folder, tests_annotations_folder


def test_annotation_bigwig():
    """
    This function tests the annotation of a VCF file using snpsift annotations.
    """

    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # Init files
        input_vcf = tests_data_folder + "/example.vcf.gz"
        input_vcf = tests_data_folder + "/example.full.vcf.gz"
        annotation_bigwig = os.path.join(tests_annotations_folder, "gerp.bw")
        annotation_bigwig2 = os.path.join(tests_annotations_folder, "gerp2.bw")
        annotation_bigbed_regions = os.path.join(
            tests_annotations_folder, "annotation_regions.bb"
        )
        output_vcf = f"{tmp_dir}/output.vcf.gz"

        # Construct param dict
        param = {
            "annotation": {
                "bigwig": {
                    "annotations": {
                        annotation_bigwig: {"INFO": None},
                        annotation_bigwig2: {"gerp": "gerp_renamed"},
                        annotation_bigbed_regions: {
                            "score1": "score_one",
                            "score2": "score2",
                            "strand": "strand",
                            "itemRgb": "itemRgb",
                        },
                    },
                    "param": {
                        "default": {
                            "method": {
                                "Integer": "mean",
                                "Float": "mean",
                                "String": "uniq",
                            }
                        },
                        annotation_bigwig2: {
                            "method": {
                                "gerp_renamed": "max",
                            }
                        },
                        annotation_bigbed_regions: {
                            "method": {
                                # "score_one": "max",
                                "score2": "min",
                                "strand": "join",
                                "itemRgb": "uniq",
                            }
                        },
                    },
                }
            }
        }

        # Create object
        variants = Variants(
            conn=None, input=input_vcf, output=output_vcf, param=param, load=True
        )

        # Remove if output file exists
        remove_if_exists([output_vcf])

        # Annotation
        variants.annotation()

        # DEVEL
        result = variants.get_query_to_df("""SELECT  * FROM variants""")
        # log.debug(f"result={result.to_string()}")

        # query annotated variant with BigWig - SNP
        result = variants.get_query_to_df(
            """SELECT 1 AS count FROM variants WHERE "#CHROM" = 'chr7' AND POS = 55249063 AND REF = 'G' AND ALT = 'A' AND INFO = 'DP=125;gerp=4.86;gerp_renamed=4.86'"""
        )
        assert len(result) == 1

        # query annotated variant with BigWig - MNP - gerp mean and gerp_renamed max
        result = variants.get_query_to_df(
            """SELECT 1 AS count FROM variants WHERE "#CHROM" = 'chr7' AND POS = 55249063 AND REF = 'GC' AND ALT = 'TT' AND INFO = 'gerp=5.39;gerp_renamed=5.92'"""
        )
        assert len(result) == 1

        # query annotated variant with BigBed
        result = variants.get_query_to_df(
            """SELECT 1 AS count FROM variants WHERE "#CHROM" = 'chr9' AND POS = 37020622 AND INFO LIKE '%score_one=878522;score2=176%'"""
        )
        assert len(result) == 1

        # Check if VCF is in correct format with pyVCF
        variants.export_output()
        try:
            vcf.Reader(filename=output_vcf)
        except:
            assert False
