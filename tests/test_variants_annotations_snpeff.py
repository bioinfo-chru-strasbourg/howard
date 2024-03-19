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


def test_annotation_snpeff():
    """
    This function tests the annotation of variants using the snpEff tool.
    """

    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # Init files
        input_vcf = tests_data_folder + "/example.vcf.gz"
        output_vcf = f"{tmp_dir}/output.vcf.gz"

        # Construct param dict
        param = {
            "annotation": {
                "snpeff": {
                    "options": "-lof -hgvs -oicr -noShiftHgvs -spliceSiteSize 3 "
                }
            }
        }

        # Create object
        variants = Variants(
            conn=None,
            input=input_vcf,
            output=output_vcf,
            config=tests_config,
            param=param,
            load=True,
        )

        # Remove if output file exists
        remove_if_exists([output_vcf])

        # Annotation
        variants.annotation()

        # query annotated variant
        result = variants.get_query_to_df(""" SELECT * FROM variants """)
        assert len(result) == 7

        # Check if VCF is in correct format with pyVCF
        variants.export_output()
        try:
            vcf.Reader(filename=output_vcf)
        except:
            assert False


def test_annotation_snpeff_full_unsorted():
    """
    This function tests the annotation of variants using the snpEff tool with specific options and
    checks if the output VCF file is in the correct format using pyVCF.
    Test with a VCF full variants type: SNV, INDEL, MNV, SV
    This VCF is unsorted
    """

    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # Init files
        input_vcf = tests_data_folder + "/example.full.unsorted.vcf.gz"
        output_vcf = f"{tmp_dir}/output.vcf.gz"

        # Construct param dict
        param = {
            "annotation": {
                "snpeff": {
                    "options": "-lof -hgvs -oicr -noShiftHgvs -spliceSiteSize 3 "
                }
            },
            "explode": {"explode_infos": True},
        }

        # Create object
        variants = Variants(
            conn=None,
            input=input_vcf,
            output=output_vcf,
            config=tests_config,
            param=param,
            load=True,
        )

        # Remove if output file exists
        remove_if_exists([output_vcf])

        # Annotation
        variants.annotation()

        # query annotated variant
        result = variants.get_query_to_df(""" SELECT INFO, "ANN" FROM variants """)
        assert len(result) == 36

        # query annotated variant as gene_fusion
        result = variants.get_query_to_df(
            """ SELECT INFO FROM variants WHERE INFO LIKE '%gene_fusion%'"""
        )
        assert len(result) == 7

        # Check if VCF is in correct format with pyVCF
        variants.export_output()
        try:
            vcf.Reader(filename=output_vcf)
        except:
            assert False


def test_annotation_snpeff_no_samples():
    """
    This function tests the annotation of variants using snpEff when there are no samples in the input
    VCF file.
    """

    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # Init files
        input_vcf = tests_data_folder + "/example.no_samples.vcf.gz"
        output_vcf = f"{tmp_dir}/output.vcf.gz"

        # Construct param dict
        param = {
            "annotation": {
                "snpeff": {
                    "options": "-lof -hgvs -oicr -noShiftHgvs -spliceSiteSize 3 "
                }
            }
        }

        # Create object
        variants = Variants(
            conn=None,
            input=input_vcf,
            output=output_vcf,
            config=tests_config,
            param=param,
            load=True,
        )

        # Remove if output file exists
        remove_if_exists([output_vcf])

        # Annotation
        variants.annotation()

        # query annotated variant
        result = variants.get_query_to_df(
            """ SELECT * FROM variants WHERE "#CHROM" = 'chr12' AND POS = 68724951 AND REF = 'G' AND ALT = 'T' AND INFO LIKE '%T|synonymous_variant|LOW|MDM1|MDM1|transcript|NM_001354969.1|protein_coding|2/15|c.69C>A|p.Ser23Ser|238/3032|69/2175|23/724||%' """
        )
        assert len(result) == 1

        # query annotated variant
        result = variants.get_query_to_df(""" SELECT * FROM variants """)
        assert len(result) == 10

        # Check if VCF is in correct format with pyVCF
        variants.export_output()
        try:
            vcf.Reader(filename=output_vcf)
        except:
            assert False


def test_annotation_quick_snpeff():
    """
    This function tests the annotation of a VCF file using the snpEff tool.
    """

    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # Init files
        input_vcf = tests_data_folder + "/example.vcf.gz"
        annotation_snpeff = "snpeff"
        output_vcf = f"{tmp_dir}/output.vcf.gz"

        # Construct param dict
        param = {"annotations": {f"{annotation_snpeff}": None}}

        # Create object
        variants = Variants(
            conn=None,
            input=input_vcf,
            output=output_vcf,
            config=tests_config,
            param=param,
            load=True,
        )

        # Remove if output file exists
        remove_if_exists([output_vcf])

        # Annotation
        variants.annotation()

        # query annotated variant
        result = variants.get_query_to_df(""" SELECT 1 AS count FROM variants """)
        assert len(result) == 7

        # Check if VCF is in correct format with pyVCF
        variants.export_output()
        try:
            vcf.Reader(filename=output_vcf)
        except:
            assert False


def test_annotation_snpeff_sqlite():
    """
    This function tests the annotation of variants using snpEff and SQLite database.
    """

    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # Init files
        input_vcf = tests_data_folder + "/example.vcf.gz"
        output_vcf = f"{tmp_dir}/output.vcf.gz"

        # Construct config dict
        config = tests_config.copy()
        config["connexion_format"] = "sqlite"

        # Construct param dict
        param = {
            "annotation": {
                "snpeff": {
                    "options": "-lof -hgvs -oicr -noShiftHgvs -spliceSiteSize 3 "
                }
            }
        }

        # Create object
        variants = Variants(
            conn=None,
            input=input_vcf,
            output=output_vcf,
            config=config,
            param=param,
            load=True,
        )

        # Remove if output file exists
        remove_if_exists([output_vcf])

        # Annotation
        variants.annotation()

        # query annotated variant
        result = variants.get_query_to_df(
            """ SELECT INFO FROM variants WHERE INFO LIKE '%ANN=%' """
        )
        assert len(result) == 7

        # Check if VCF is in correct format with pyVCF
        variants.export_output()
        try:
            vcf.Reader(filename=output_vcf)
        except:
            assert False
