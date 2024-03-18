# -*- coding: utf-8 -*-
"""
Tests

Usage:
pytest tests/

Coverage:
coverage run -m pytest tests/test_objects_variants.py -x -v --log-cli-level=INFO --capture=tee-sys
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

from howard.commons import *
from howard.objects.variants import Variants
from howard.tools.databases import *
from test_needed import *



    
def test_annotation_hgvs():
    """
    The function `test_annotation_hgvs` tests the annotation of a VCF file using bcftools and SQLite.
    """

    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # Init files
        input_vcf = tests_data_folder + "/example.vcf.gz"
        output_vcf = f"{tmp_dir}/output.vcf.gz"

        # Construct config dict
        config = tests_config.copy()

        # Download database
        download_needed_databases()

        # Construct param dict
        param = {"hgvs": {"use_exon": True, "use_version": True}}

        # Create object
        variants = Variants(conn=None, input=input_vcf, output=output_vcf, param=param, config=config, load=True)

        # Remove if output file exists
        remove_if_exists([output_vcf])

        # Annotation
        variants.annotation_hgvs()

        # Check
        result = variants.get_query_to_df("""SELECT * FROM variants WHERE INFO LIKE '%hgvs%'""")
        assert len(result) == 7
        result = variants.get_query_to_df("""SELECT INFO FROM variants WHERE "#CHROM" = 'chr7' AND POS = 55249063 AND REF = 'G' AND ALT = 'A' AND INFO LIKE '%NM_001346897.2(exon19):c.2226G>A%'""")
        assert len(result) == 1

        # Gene Protein

        # Construct param dict
        param = {"hgvs": {"add_protein": True, "use_gene": True}}

        # Create object
        variants = Variants(conn=None, input=input_vcf, output=output_vcf, param=param, config=config, load=True)

        # Remove if output file exists
        remove_if_exists([output_vcf])

        # Annotation
        variants.annotation_hgvs()

        #  Check
        result = variants.get_query_to_df("""SELECT * FROM variants WHERE INFO LIKE '%hgvs%'""")
        assert len(result) == 7
        result = variants.get_query_to_df("""SELECT INFO FROM variants WHERE "#CHROM" = 'chr7' AND POS = 55249063 AND REF = 'G' AND ALT = 'A' AND INFO LIKE '%NM_001346897(EGFR):c.2226G>A%' AND INFO LIKE '%NP_001333826(EGFR):p.Gln742Gln%'""")
        assert len(result) == 1

