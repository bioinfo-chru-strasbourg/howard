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


def test_annotation_bcftools():
    """
    This function tests the annotation of a VCF file using bcftools annotations.
    """

    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # Init files
        input_vcf = tests_data_folder + "/example.vcf.gz"
        annotation_parquet = os.path.join(tests_annotations_folder, "nci60.vcf.gz")
        output_vcf = f"{tmp_dir}/output.vcf.gz"

        # Construct param dict
        param = {
            "annotation": {
                "bcftools": {"annotations": {annotation_parquet: {"INFO": None}}}
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

        # query annotated variant
        result = variants.get_query_to_df(
            """SELECT 1 AS count FROM variants WHERE "#CHROM" = 'chr7' AND POS = 55249063 AND REF = 'G' AND ALT = 'A' AND INFO = 'DP=125;nci60=0.66'"""
        )
        assert len(result) == 1

        # Check if VCF is in correct format with pyVCF
        variants.export_output()
        try:
            vcf.Reader(filename=output_vcf)
        except:
            assert False


def test_annotation_bcftools_rename_field():
    """
    This function tests the annotation of a VCF file using bcftools annotations.
    """

    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # Init files
        input_vcf = tests_data_folder + "/example.vcf.gz"
        annotation_parquet = os.path.join(tests_annotations_folder, "nci60.vcf.gz")
        output_vcf = f"{tmp_dir}/output.vcf.gz"

        # Construct param dict
        param = {
            "annotation": {
                "bcftools": {"annotations": {annotation_parquet: {"nci60": "nci61"}}}
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

        # query annotated variant
        result = variants.get_query_to_df(
            """SELECT 1 AS count FROM variants WHERE "#CHROM" = 'chr7' AND POS = 55249063 AND REF = 'G' AND ALT = 'A' AND INFO = 'DP=125;nci61=0.66'"""
        )
        assert len(result) == 1

        # Check if VCF is in correct format with pyVCF
        variants.export_output()
        try:
            vcf.Reader(filename=output_vcf)
        except:
            assert False


def test_annotation_bcftools_bed():
    """
    This function tests the annotation of a VCF file using bcftools and a bed file.
    """

    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # Init files
        input_vcf = tests_data_folder + "/example.vcf.gz"
        annotation_parquet = os.path.join(tests_annotations_folder, "refGene.bed.gz")
        output_vcf = f"{tmp_dir}/output.vcf.gz"

        # Construct param dict
        param = {
            "annotation": {
                "bcftools": {"annotations": {annotation_parquet: {"symbol": None}}}
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

        # query annotated variant
        result = variants.get_query_to_df(
            """SELECT 1 AS count FROM variants WHERE "#CHROM" = 'chr7' AND POS = 55249063 AND REF = 'G' AND ALT = 'A' AND INFO = 'DP=125;symbol=EGFR,EGFR-AS1'"""
        )
        assert len(result) == 1

        # Check if VCF is in correct format with pyVCF
        variants.export_output()
        try:
            vcf.Reader(filename=output_vcf)
        except:
            assert False


def test_annotation_bcftools_sqlite():
    """
    This function tests the annotation of a VCF file using bcftools and SQLite.
    """

    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # Init files
        input_vcf = tests_data_folder + "/example.vcf.gz"
        annotation_parquet = os.path.join(tests_annotations_folder, "nci60.vcf.gz")
        output_vcf = f"{tmp_dir}/output.vcf.gz"

        # Construct config dict
        config = tests_config.copy()
        config["connexion_format"] = "sqlite"

        # Construct param dict
        param = {
            "annotation": {
                "bcftools": {"annotations": {annotation_parquet: {"INFO": None}}}
            }
        }

        # Create object
        variants = Variants(
            conn=None,
            input=input_vcf,
            output=output_vcf,
            param=param,
            config=config,
            load=True,
        )

        # Remove if output file exists
        remove_if_exists([output_vcf])

        # Annotation
        variants.annotation()

        result = variants.get_query_to_df(
            """SELECT * FROM variants WHERE "#CHROM" = 'chr7' AND POS = 55249063 AND REF = 'G' AND ALT = 'A' AND INFO = 'DP=125;nci60=0.66'"""
        )
        assert len(result) == 1

        # Check if VCF is in correct format with pyVCF
        variants.export_output()
        try:
            vcf.Reader(filename=output_vcf)
        except:
            assert False
