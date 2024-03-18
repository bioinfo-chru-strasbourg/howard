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




def test_annotation_parquet_append():
    """
    The function `test_annotation_parquet_append` tests the annotation functionality for appending data
    to a Parquet file in a VCF file.
    """

    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # Init files
        input_vcf = tests_data_folder + "/example.nci60_1.vcf"
        annotation1 = os.path.join(tests_annotations_folder, "nci60.parquet")
        output_vcf = f"{tmp_dir}/output.vcf.gz"

        # Construct param dict
        param = {
                    'annotation': {
                        'parquet': {
                            'annotations': {
                                annotation1: {
                                    "nci60": "nci60"
                                }
                            },
                        },
                        'options': {
                            'annotations_append': False
                        }
                    }
                }
        param_update = {
                    'annotation': {
                        'parquet': {
                            'annotations': {
                                annotation1: {
                                    "nci60": "nci60"
                                }
                            },
                        },
                        'options': {
                            'annotations_append': True
                        }
                        
                    }
                }
        log.debug(f"param={param}")
        log.debug(f"param_update={param_update}")

        # Create object
        variants = Variants(conn=None, input=input_vcf, output=output_vcf, param=param, load=True)

        # Remove if output file exists
        remove_if_exists([output_vcf])

        # Annotation
        variants.annotation()

        # Check annotation not changed
        
        result = variants.get_query_to_df("SELECT INFO FROM variants")
        log.debug(result)
        result1 = variants.get_query_to_df("SELECT 1 AS count FROM variants WHERE \"#CHROM\" = 'chr1' AND POS = 768253 AND REF = 'A' AND ALT = 'G' AND INFO LIKE '%nci60=0.321%'")
        result2 = variants.get_query_to_df("SELECT 1 AS count FROM variants WHERE \"#CHROM\" = 'chr7' AND POS = 55249063 AND REF = 'G' AND ALT = 'A' AND INFO LIKE '%nci60=0.66%'")
        #log.debug(result1)
        assert len(result1) == 1
        assert len(result2) == 0

        variants.set_param(param=param_update)
        variants.annotation()

        # Check annotation changed (existing kept, one annotation added)
        result = variants.get_query_to_df("SELECT INFO FROM variants")
        log.debug(result)
        result1 = variants.get_query_to_df("SELECT 1 AS count FROM variants WHERE \"#CHROM\" = 'chr1' AND POS = 768253 AND REF = 'A' AND ALT = 'G' AND INFO LIKE '%nci60=0.321%'")
        result2 = variants.get_query_to_df("SELECT 1 AS count FROM variants WHERE \"#CHROM\" = 'chr7' AND POS = 55249063 AND REF = 'G' AND ALT = 'A' AND INFO LIKE '%nci60=0.66%'")
        #log.debug(result)
        assert len(result1) == 1
        assert len(result2) == 1

        # Check if VCF is in correct format with pyVCF
        variants.export_output()
        try:
            vcf.Reader(filename=output_vcf)
        except:
            assert False


def test_annotation_parquet_update():
    """
    The function `test_annotation_parquet_update` tests the updating functionality of annotations in a
    VCF file using Parquet format.
    """

    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # Init files
        input_vcf = tests_data_folder + "/example.nci60_2.vcf"
        annotation1 = os.path.join(tests_annotations_folder, "nci60.parquet")
        output_vcf = f"{tmp_dir}/output.vcf.gz"

        # Construct param dict
        param = {
                    'annotation': {
                        'parquet': {
                            'annotations': {
                                annotation1: {
                                    "nci60": "nci60"
                                }
                            },
                        },
                        'options': {
                            'annotations_update': False
                        }
                    }
                }
        param_update = {
                    'annotation': {
                        'parquet': {
                            'annotations': {
                                annotation1: {
                                    "nci60": "nci60"
                                }
                            },
                        },
                        'options': {
                            'annotations_update': True
                        }
                    }
                }

        # Create object
        variants = Variants(conn=None, input=input_vcf, output=output_vcf, param=param, load=True)

        # Remove if output file exists
        remove_if_exists([output_vcf])

        # Annotation
        variants.annotation()

        # Check annotation not changed
        result1 = variants.get_query_to_df("SELECT 1 AS count FROM variants WHERE \"#CHROM\" = 'chr1' AND POS = 768253 AND REF = 'A' AND ALT = 'G' AND INFO LIKE '%nci60=0.321%'")
        result2 = variants.get_query_to_df("SELECT 1 AS count FROM variants WHERE \"#CHROM\" = 'chr7' AND POS = 55249063 AND REF = 'G' AND ALT = 'A' AND INFO LIKE '%nci60=0.123%'")
        assert len(result1) == 1
        assert len(result2) == 1

        variants.set_param(param=param_update)
        variants.annotation()

        # Check annotation changed (all removed, but one added)
        result1 = variants.get_query_to_df("SELECT 1 AS count FROM variants WHERE \"#CHROM\" = 'chr1' AND POS = 768253 AND REF = 'A' AND ALT = 'G' AND INFO LIKE '%nci60=0.321%'")
        result2 = variants.get_query_to_df("SELECT 1 AS count FROM variants WHERE \"#CHROM\" = 'chr7' AND POS = 55249063 AND REF = 'G' AND ALT = 'A' AND INFO LIKE '%nci60=0.66%'")
        assert len(result1) == 0
        assert len(result2) == 1

        # Check if VCF is in correct format with pyVCF
        variants.export_output()
        try:
            vcf.Reader(filename=output_vcf)
        except:
            assert False


def test_annotations_parquet_all_available_annotations_databases():
    """
    The function `test_annotations_all_available_annotations_databases` tests the annotation
    functionality of a `Variants` object using all available databases.
    """

    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # Init files
        input_vcf = tests_data_folder + "/example.vcf.gz"
        output_vcf = f"{tmp_dir}/output.vcf.gz"

        # Config
        config = tests_config.copy()

        # Construct param annotations with all available parquet databases (with header)
        parquet_files = list(set(
            glob.glob(rf"{tests_annotations_folder}/*parquet", recursive=False)
            + glob.glob(rf"{tests_annotations_folder}/*duckdb", recursive=False)
            + glob.glob(rf"{tests_annotations_folder}/*vcf", recursive=False)
            + glob.glob(rf"{tests_annotations_folder}/*vcf.gz", recursive=False)
            + glob.glob(rf"{tests_annotations_folder}/*tsv", recursive=False)
            + glob.glob(rf"{tests_annotations_folder}/*tsv.gz", recursive=False)
            + glob.glob(rf"{tests_annotations_folder}/*csv", recursive=False)
            + glob.glob(rf"{tests_annotations_folder}/*csv.gz", recursive=False)
            + glob.glob(rf"{tests_annotations_folder}/*json", recursive=False)
            + glob.glob(rf"{tests_annotations_folder}/*json.gz", recursive=False)
        ))

        param_annotation = {
                                'parquet': {
                                    'annotations': {}
                                }
                            }
        for parquet_file in parquet_files:
            if os.path.exists(parquet_file) and os.path.exists(parquet_file+".hdr") and "fail" not in parquet_file:
                param_annotation["parquet"]["annotations"][parquet_file] = {'INFO': None}

        param = {
            "annotation": param_annotation,
            "assembly": "hg19"
            }

        # Create object
        variants = Variants(conn=None, input=input_vcf, output=output_vcf, config=config, param=param, load=True)

        # Remove if output file exists
        remove_if_exists([output_vcf])

        # Annotation
        try:
            variants.annotation()
            assert True
        except:
            assert False

        # check param
        param_input = variants.get_param()
        expected_param = param

        assert param_input == expected_param

        # Check if VCF is in correct format with pyVCF
        variants.export_output()
        try:
            vcf.Reader(filename=output_vcf)
        except:
            assert False


def test_annotations_parquet_no_samples():
    """
    This function tests the annotation of a VCF file without samples using various annotations
    and checks if the output VCF file is in the correct format.
    """

    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # Init files
        input_vcf = tests_data_folder + "/example.no_samples.vcf.gz"
        annotation1 = os.path.join(tests_annotations_folder, "nci60.parquet")
        annotation2 = tests_data_folder + "/example.vcf.gz"
        annotation3 = os.path.join(tests_annotations_folder, "refGene.bed.gz")
        output_vcf = f"{tmp_dir}/output.vcf.gz"

        # Construct param dict
        param_annotations = {
                annotation1: {"INFO": None},
                annotation2: {"CLNSIG": "CLNSIG_new"},
                annotation3: {"symbol": "gene"},
                }
        param = {"annotations": param_annotations }

        # Create object
        variants = Variants(conn=None, input=input_vcf, output=output_vcf, param=param, load=True)

        # Remove if output file exists
        remove_if_exists([output_vcf])

        # Annotation
        variants.annotation()

        # check param
        param_input = variants.get_param()
        expected_param = {'annotations': param_annotations,
                'annotation': {
                    'parquet': {
                        'annotations': {
                            annotation1: {'INFO': None},
                            annotation2: {'CLNSIG': 'CLNSIG_new'},
                            annotation3: {'symbol': 'gene'}
                        }
                    }
                }
            }

        assert param_input == expected_param

        result = variants.get_query_to_df("SELECT * FROM variants")
        assert len(result) == 10

        # Check if VCF is in correct format with pyVCF
        variants.export_output()
        try:
            vcf.Reader(filename=output_vcf)
        except:
            assert False


def test_annotation_parquet_with_all_formats():
    """
    This function tests the `annotation()` method of the `Variants` class using a Parquet file as
    annotation source with various formats.
    """
    
    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        for annotation_format in ["vcf", "vcf.gz", "tsv", "tsv.gz", "csv", "csv.gz", "json", "json.gz", "tbl", "tbl.gz", "parquet", "partition.parquet", "duckdb"]:

            # Init files
            input_vcf = tests_data_folder + "/example.vcf.gz"
            annotation_parquet = os.path.join(tests_annotations_folder, f"nci60.{annotation_format}")
            output_vcf = f"{tmp_dir}/output.vcf.gz"

            # Construct param dict
            param = {"annotation": {"parquet": {"annotations": {annotation_parquet: {"INFO": None}}}}}

            # Create object
            variants = Variants(conn=None, input=input_vcf, output=output_vcf, param=param, load=True)

            # Remove if output file exists
            remove_if_exists([output_vcf])

            # Annotation
            variants.annotation()

            # query annotated variant
            result = variants.get_query_to_df("SELECT 1 AS count FROM variants WHERE \"#CHROM\" = 'chr7' AND POS = 55249063 AND REF = 'G' AND ALT = 'A' AND INFO = 'DP=125;nci60=0.66'")
            length = len(result)
            
            assert length == 1

            # Check if VCF is in correct format with pyVCF
            variants.export_output()
            try:
                vcf.Reader(filename=output_vcf)
            except:
                assert False


def test_annotation_parquet_regions():
    """
    The `test_annotation_parquet_regions()` function tests the `annotation()` method of the `Variants`
    class using a Parquet file as an annotation source with various formats.
    """
    
    # Init files
    input_vcf = tests_data_folder + "/example.vcf.gz"

    # annotation regions
    with TemporaryDirectory(dir=tests_folder) as tmp_dir:
        
        # Init
        annotation_parquet = os.path.join(tests_annotations_folder, f"annotation_regions.bed.gz")
        output_vcf = f"{tmp_dir}/output.vcf.gz"

        # Construct param dict
        param = {"annotation": {"parquet": {"annotations": {annotation_parquet: {"INFO": None}}}}}

        # Create object
        variants = Variants(conn=None, input=input_vcf, output=output_vcf, param=param, load=True)

        # Annotation
        variants.annotation()

        # query annotated variant
        result = variants.get_query_to_df("SELECT 1 AS count FROM variants WHERE INFO LIKE '%blue%' OR INFO LIKE '%red%' OR INFO LIKE '%orange%' OR INFO LIKE '%cherry%'")
        length = len(result)
        
        assert length == 3

        # query annotated variant
        result = variants.get_query_to_df("SELECT 1 AS count FROM variants WHERE INFO LIKE '%yellow%' OR INFO LIKE '%banana%'")
        length = len(result)
        
        assert length == 3

        # query annotated variant
        result = variants.get_query_to_df("SELECT 1 AS count FROM variants WHERE INFO NOT LIKE '%annot1%' AND INFO NOT LIKE '%annot2%'")
        length = len(result)
        
        assert length == 1

        # Check if VCF is in correct format with pyVCF
        variants.export_output()
        try:
            vcf.Reader(filename=output_vcf)
        except:
            assert False

    # annotation regions with refgene and an associated header hdr
    with TemporaryDirectory(dir=tests_folder) as tmp_dir:
        
        # Init
        annotation_parquet = os.path.join(tests_annotations_folder, f"refGene.bed")
        output_vcf = f"{tmp_dir}/output.vcf.gz"

        # Construct param dict
        param = {"annotation": {"parquet": {"annotations": {annotation_parquet: {"INFO": None}}}}}

        # Create object
        variants = Variants(conn=None, input=input_vcf, output=output_vcf, param=param, load=True)

        # Annotation
        variants.annotation()

        # query annotated variant
        result = variants.get_query_to_df("SELECT 1 AS count FROM variants WHERE INFO LIKE '%symbol%' OR INFO LIKE '%transcripts%' OR INFO LIKE '%strand%'")
        length = len(result)
        
        assert length == 3

        # Check if VCF is in correct format with pyVCF
        variants.export_output()
        try:
            vcf.Reader(filename=output_vcf)
        except:
            assert False

    # annotation regions with refgene compressed gz and an associated header hdr
    with TemporaryDirectory(dir=tests_folder) as tmp_dir:
        
        # Init
        annotation_parquet = os.path.join(tests_annotations_folder, f"refGene.bed.gz")
        output_vcf = f"{tmp_dir}/output.vcf.gz"

        # Construct param dict
        param = {"annotation": {"parquet": {"annotations": {annotation_parquet: {"INFO": None}}}}}

        # Create object
        variants = Variants(conn=None, input=input_vcf, output=output_vcf, param=param, load=True)

        # Annotation
        variants.annotation()

        # query annotated variant
        result = variants.get_query_to_df("SELECT 1 AS count FROM variants WHERE INFO LIKE '%symbol%' OR INFO LIKE '%transcripts%' OR INFO LIKE '%strand%'")
        length = len(result)
        
        assert length == 3

        # Check if VCF is in correct format with pyVCF
        variants.export_output()
        try:
            vcf.Reader(filename=output_vcf)
        except:
            assert False

    # annotation regions with refgene without any associated header hdr and INFO annotation (no annotation)
    with TemporaryDirectory(dir=tests_folder) as tmp_dir:
        
        # Init
        annotation_parquet = os.path.join(tests_annotations_folder, f"refGene.without_header.bed")
        output_vcf = f"{tmp_dir}/output.vcf.gz"

        # Construct param dict
        param = {"annotation": {"parquet": {"annotations": {annotation_parquet: {"INFO": None}}}}}

        # Create object
        variants = Variants(conn=None, input=input_vcf, output=output_vcf, param=param, load=True)

        # Annotation
        variants.annotation()

        # query annotated variant
        result = variants.get_query_to_df("SELECT 1 AS count FROM variants WHERE INFO LIKE '%column3%' OR INFO LIKE '%column4%' OR INFO LIKE '%column5%'")
        length = len(result)
        
        assert length == 0

        # Check if VCF is in correct format with pyVCF
        variants.export_output()
        try:
            vcf.Reader(filename=output_vcf)
        except:
            assert False

    # annotation regions with refgene without any associated header hdr and ALL annotation (annotation with "column...)")
    with TemporaryDirectory(dir=tests_folder) as tmp_dir:
        
        # Init
        annotation_parquet = os.path.join(tests_annotations_folder, f"refGene.without_header.bed")
        output_vcf = f"{tmp_dir}/output.vcf.gz"

        # Construct param dict
        param = {"annotation": {"parquet": {"annotations": {annotation_parquet: {"ALL": None}}}}}

        # Create object
        variants = Variants(conn=None, input=input_vcf, output=output_vcf, param=param, load=True)

        # Annotation
        variants.annotation()

        # query annotated variant
        result = variants.get_query_to_df("SELECT 1 AS count FROM variants WHERE INFO LIKE '%column3%' OR INFO LIKE '%column4%' OR INFO LIKE '%column5%'")
        length = len(result)
        
        assert length == 3

        # Check if VCF is in correct format with pyVCF
        variants.export_output()
        try:
            vcf.Reader(filename=output_vcf)
        except:
            assert False


def test_annotation_parquet_field_already_in_vcf():
    """
    This function tests if a field already present in a VCF file is not changed during annotation with a
    Parquet file.
    """

    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # Init files
        input_vcf = tests_data_folder + "/example.vcf.gz"
        annotation1 = os.path.join(tests_annotations_folder, "nci60.parquet")
        output_vcf = f"{tmp_dir}/output.vcf.gz"

        # Construct param dict
        param_annotations = {
                annotation1: {"nci60": "DP"}
                }
        param = {"annotations": param_annotations }

        # Create object
        variants = Variants(conn=None, input=input_vcf, output=output_vcf, param=param, load=True)

        # Remove if output file exists
        remove_if_exists([output_vcf])

        # Annotation
        variants.annotation()

        # check param
        param_input = variants.get_param()
        expected_param = {'annotations': param_annotations,
                        'annotation': {
                            'parquet': {'annotations': {annotation1: {"nci60": "DP"}}}
                        }
                        }

        assert param_input and expected_param

        # Check annotation not changed
        result = variants.get_query_to_df("SELECT 1 AS count FROM variants WHERE \"#CHROM\" = 'chr7' AND POS = 55249063 AND REF = 'G' AND ALT = 'A' AND INFO = 'DP=125'")
        assert len(result) == 1

        # Check if VCF is in correct format with pyVCF
        variants.export_output()
        try:
            vcf.Reader(filename=output_vcf)
        except:
            assert False


def test_annotation_parquet_duckdb():
    """
    This function tests the annotation of variants using DuckDB.
    """

    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # Init files
        annotation_parquet = os.path.join(tests_annotations_folder, "nci60.parquet")
        annotation_duckdb = f"{tmp_dir}/nci60.duckdb"

        remove_if_exists([annotation_duckdb])

        annotation_database = Variants(input=annotation_parquet, output=annotation_duckdb, load=True)
        annotation_database.export_output()

        # Test annotation with duckdb database

        # Init files
        input_vcf = tests_data_folder + "/example.vcf.gz"
        output_vcf = f"{tmp_dir}/output.vcf.gz"

        # Construct param dict
        param = {"annotation": {"parquet": {"annotations": {annotation_duckdb: {"INFO": None}}}}}

        # Create object
        variants = Variants(conn=None, input=input_vcf, output=output_vcf, param=param, load=True)

        # Remove if output file exists
        remove_if_exists([output_vcf])

        # Annotation
        variants.annotation()

        # query annotated variant
        result = variants.get_query_to_df("SELECT 1 AS count FROM variants WHERE \"#CHROM\" = 'chr7' AND POS = 55249063 AND REF = 'G' AND ALT = 'A' AND INFO = 'DP=125;nci60=0.66'")
        assert len(result) == 1

        # Check if VCF is in correct format with pyVCF
        variants.export_output()
        try:
            vcf.Reader(filename=output_vcf)
        except:
            assert False

    