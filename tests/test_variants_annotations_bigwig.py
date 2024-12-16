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
import sys
from tempfile import TemporaryDirectory
import duckdb
import re
import Bio.bgzf as bgzf
import gzip
import pytest

from howard.functions.commons import *
from howard.objects.variants import Variants
from howard.functions.databases import *
from test_needed import *


def test_annotation_bigwig():
    """
    This function tests the annotation of a VCF file using snpsift annotations.
    """

    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # Init files
        input_vcf = tests_data_folder + "/example.vcf.gz"
        annotation_bigwig = os.path.join(tests_annotations_folder, "gerp.bw")
        annotation_bigwig2 = os.path.join(tests_annotations_folder, "gerp2.bw")
        output_vcf = f"{tmp_dir}/output.vcf.gz"

        # Construct param dict
        param = {
            "annotation": {
                "bigwig": {"annotations": {annotation_bigwig: {"INFO": None}, annotation_bigwig2: {"gerp": "gerp_renamed"}}}
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
        log.debug(f"result={result}")

        # query annotated variant
        result = variants.get_query_to_df(
            """SELECT 1 AS count FROM variants WHERE "#CHROM" = 'chr7' AND POS = 55249063 AND REF = 'G' AND ALT = 'A' AND INFO = 'DP=125;gerp=4.86;gerp_renamed=4.86'"""
        )
        assert len(result) == 1

        # Check if VCF is in correct format with pyVCF
        variants.export_output()
        try:
            vcf.Reader(filename=output_vcf)
        except:
            assert False
