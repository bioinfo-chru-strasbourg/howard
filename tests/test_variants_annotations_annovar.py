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

from howard.functions.commons import *
from howard.objects.variants import Variants
from howard.functions.databases import *
from test_needed import *




def test_annotation_annovar():
    """
    This function tests the annotation of variants using Annovar and checks if the output VCF file is in
    the correct format.
    """

    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # Init files
        input_vcf = tests_data_folder + "/example.vcf.gz"
        annotation_annovar = "nci60"
        output_vcf = f"{tmp_dir}/output.vcf.gz"

        # Construct param dict
        param = {"annotation": {"annovar": {"annotations":  {annotation_annovar: {"INFO": None}}}}}

        # Create object
        variants = Variants(conn=None, input=input_vcf, output=output_vcf, config=tests_config, param=param, load=True)

        # Remove if output file exists
        remove_if_exists([output_vcf])

        # Annotation
        variants.annotation()

        # query annotated variant
        result = variants.get_query_to_df("""SELECT 1 AS count FROM variants WHERE "#CHROM" = 'chr7' AND POS = 55249063 AND REF = 'G' AND ALT = 'A' AND INFO = 'DP=125;nci60=0.66'""")
        assert len(result) == 1

        # Check if VCF is in correct format with pyVCF
        variants.export_output()
        try:
            vcf.Reader(filename=output_vcf)
        except:
            assert False


def test_annotation_annovar_full_unsorted():
    """
    The function tests the annotation of variants using Annovar and checks if the output VCF file is in
    the correct format.
    Test with a VCF full variants type: SNV, INDEL, MNV, SV
    This VCF is unsorted
    """

    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # Init files
        input_vcf = tests_data_folder + "/example.full.unsorted.vcf.gz"
        annotation_annovar = "nci60"
        output_vcf = f"{tmp_dir}/output.vcf.gz"

        # Construct param dict
        param = {"annotation": {"annovar": {"annotations":  {annotation_annovar: {"INFO": None}}}}}

        # Create object
        variants = Variants(conn=None, input=input_vcf, output=output_vcf, config=tests_config, param=param, load=True)

        # Remove if output file exists
        remove_if_exists([output_vcf])

        # Annotation
        variants.annotation()

        # query annotated variant
        result = variants.get_query_to_df("""SELECT 1 AS count FROM variants WHERE "#CHROM" = 'chr7' AND POS = 55249063 AND REF = 'G' AND ALT = 'A' AND INFO = 'DP=125;nci60=0.66'""")
        assert len(result) == 1

        # Check if VCF is in correct format with pyVCF
        variants.export_output()
        try:
            vcf.Reader(filename=output_vcf)
        except:
            assert False


def test_annotation_annovar_no_samples():
    """
    This function tests the annotation of a VCF file using Annovar when there are no samples present.
    """

    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # Init files
        input_vcf = tests_data_folder + "/example.no_samples.vcf.gz"
        annotation_annovar = "nci60"
        output_vcf = f"{tmp_dir}/output.vcf.gz"

        # Construct param dict
        param = {"annotation": {"annovar": {"annotations":  {annotation_annovar: {"INFO": None}}}}}

        # Create object
        variants = Variants(conn=None, input=input_vcf, output=output_vcf, config=tests_config, param=param, load=True)

        # Remove if output file exists
        remove_if_exists([output_vcf])

        # Annotation
        variants.annotation()

        # query annotated variant
        result = variants.get_query_to_df(""" SELECT 1 AS count FROM variants WHERE "#CHROM" = 'chr12' AND POS = 68724951 AND REF = 'G' AND ALT = 'T' AND INFO = 'nci60=0.77' """)
        assert len(result) == 1

        # Check if VCF is in correct format with pyVCF
        variants.export_output()
        try:
            vcf.Reader(filename=output_vcf)
        except:
            assert False


def test_annotation_annovar_sqlite():
    """
    This function tests the annotation of variants using Annovar and SQLite database.
    """

    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # Init files
        input_vcf = tests_data_folder + "/example.vcf.gz"
        annotation_annovar = "nci60"
        output_vcf = f"{tmp_dir}/output.vcf.gz"

        # Construct config dict
        config = tests_config.copy()
        config["connexion_format"] = "sqlite"

        # Construct param dict
        param = {"annotation": {"annovar": {"annotations":  {annotation_annovar: {"INFO": None}}}}}

        # Create object
        variants = Variants(conn=None, input=input_vcf, output=output_vcf, config=config, param=param, load=True)

        # Remove if output file exists
        remove_if_exists([output_vcf])

        # Annotation
        variants.annotation()

        # query annotated variant
        result = variants.get_query_to_df("""SELECT 1 AS count FROM variants WHERE "#CHROM" = 'chr7' AND POS = 55249063 AND REF = 'G' AND ALT = 'A' AND INFO = 'DP=125;nci60=0.66'""")
        assert len(result) == 1

        # Check if VCF is in correct format with pyVCF
        variants.export_output()
        try:
            vcf.Reader(filename=output_vcf)
        except:
            assert False


def test_annotation_quick_annovar():
    """
    This function tests the annotation of a VCF file using Annovar.
    """

    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # Init files
        input_vcf = tests_data_folder + "/example.vcf.gz"
        annotation_annovar = "nci60"
        output_vcf = f"{tmp_dir}/output.vcf.gz"

        # Construct param dict
        param = {"annotations": {
                    f"annovar:{annotation_annovar}": None
                    }
        }

        # Create object
        variants = Variants(conn=None, input=input_vcf, output=output_vcf, config=tests_config, param=param, load=True)

        # Remove if output file exists
        remove_if_exists([output_vcf])

        # Annotation
        variants.annotation()

        # query annotated variant
        result = variants.get_query_to_df("""SELECT 1 AS count FROM variants WHERE "#CHROM" = 'chr7' AND POS = 55249063 AND REF = 'G' AND ALT = 'A' AND INFO = 'DP=125;nci60=0.66'""")
        assert len(result) == 1
        
        # Check if VCF is in correct format with pyVCF
        variants.export_output()
        try:
            vcf.Reader(filename=output_vcf)
        except:
            assert False

