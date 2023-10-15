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




def test_annotations():
    """
    This function tests the annotation process of a VCF file with multiple annotations.

    The function initializes input VCF and annotation files, constructs a parameter dictionary with different annotation types,
    creates a Variants object with the input file, parameter dictionary, and output file, and tests the annotation process.

    The function then checks the parameter dictionary of the Variants object, and tests the output VCF file for the presence of
    annotated variants using SQL queries.

    Finally, the function exports the output VCF file and checks if it is in the correct format with pyVCF.

    :raises AssertionError: If any of the tests fail.
    """

    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # Init files
        input_vcf = tests_data_folder + "/example.vcf.gz"
        # annotation1 = "nci60.parquet"
        # annotation2 = tests_data_folder + "/example.vcf.gz"
        # annotation3 = "refGene.bed.gz"

        annotation1 = database_files.get("parquet")
        annotation2 = database_files.get("example_vcf_gz")
        annotation3 = database_files.get("refgene_gz")


        output_vcf = f"{tmp_dir}/output.vcf.gz"

        # Config
        config = tests_config.copy()

        # Construct param dict
        param_annotation = {
                                'parquet': {
                                    'annotations': {
                                        annotation1: {
                                            'INFO': None
                                        },
                                        annotation2: {
                                            'CLNSIG': 'CLNSIG_new'
                                        }
                                    },
                                },
                                'bcftools': {
                                    'annotations': {
                                        annotation2: {
                                            'CLNSIG': 'CLNSIG_new_bcftools'
                                        },
                                        annotation3: {
                                            'symbol': 'gene'
                                        }
                                    },
                                }
                            }
        param = {
            "annotation": param_annotation,
            "assembly": "hg19"
            }

        # Create object
        variants = Variants(conn=None, input=input_vcf, output=output_vcf, config=config, param=param, load=True)

        # Remove if output file exists
        remove_if_exists([output_vcf])

        # Annotation
        variants.annotation()

        # check param
        param_input = variants.get_param()
        expected_param = param

        assert param_input == expected_param

        # Check annotation1
        result = variants.get_query_to_df("SELECT 1 AS count FROM variants WHERE \"#CHROM\" = 'chr1' AND POS = 28736 AND REF = 'A' AND ALT = 'C' AND INFO LIKE '%CLNSIG_new=%'")
        assert len(result) == 1

        # Check annotation2
        result = variants.get_query_to_df("SELECT 1 AS count FROM variants WHERE \"#CHROM\" = 'chr7' AND POS = 55249063 AND REF = 'G' AND ALT = 'A' AND INFO = 'DP=125;nci60=0.66;gene=EGFR,EGFR-AS1'")
        assert len(result) == 1

        # Check if VCF is in correct format with pyVCF
        variants.export_output()
        try:
            vcf.Reader(filename=output_vcf)
        except:
            assert False



def test_annotations_all_available_annotations_databases():
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


def test_annotations_no_samples():
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
                            'parquet': {'annotations': {annotation1: {'INFO': None}}},
                            'bcftools': {'annotations': {annotation2: {'CLNSIG': 'CLNSIG_new'}, annotation3: {'symbol': 'gene'}}}}
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


def test_annotation_duckdb():
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
        param = {"annotation": {"bcftools": {"annotations":  {annotation_parquet: {"INFO": None}}}}}

        # Create object
        variants = Variants(conn=None, input=input_vcf, output=output_vcf, param=param, load=True)

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
        param = {"annotation": {"bcftools": {"annotations":  {annotation_parquet: {"symbol": None}}}}}

        # Create object
        variants = Variants(conn=None, input=input_vcf, output=output_vcf, param=param, load=True)

        # Remove if output file exists
        remove_if_exists([output_vcf])

        # Annotation
        variants.annotation()

        # query annotated variant
        result = variants.get_query_to_df("""SELECT 1 AS count FROM variants WHERE "#CHROM" = 'chr7' AND POS = 55249063 AND REF = 'G' AND ALT = 'A' AND INFO = 'DP=125;symbol=EGFR,EGFR-AS1'""")
        assert len(result) == 1

        # Check if VCF is in correct format with pyVCF
        variants.export_output()
        try:
            vcf.Reader(filename=output_vcf)
        except:
            assert False


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


def test_annotation_snpeff():
    """
    This function tests the annotation of variants using the snpEff tool.
    """

    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # Init files
        input_vcf = tests_data_folder + "/example.vcf.gz"
        output_vcf = f"{tmp_dir}/output.vcf.gz"

        # Construct param dict
        param = {"annotation": {"snpeff": {"options": "-lof -hgvs -oicr -noShiftHgvs -spliceSiteSize 3 "}}}

        # Create object
        variants = Variants(conn=None, input=input_vcf, output=output_vcf, config=tests_config, param=param, load=True)

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
            "annotation": {"snpeff": {"options": "-lof -hgvs -oicr -noShiftHgvs -spliceSiteSize 3 "}},
            "explode_infos": True     
                }

        # Create object
        variants = Variants(conn=None, input=input_vcf, output=output_vcf, config=tests_config, param=param, load=True)

        # Remove if output file exists
        remove_if_exists([output_vcf])

        # Annotation
        variants.annotation()

        # query annotated variant
        result = variants.get_query_to_df(""" SELECT INFO, "ANN" FROM variants """)
        assert len(result) == 36
        
        # query annotated variant as gene_fusion
        result = variants.get_query_to_df(""" SELECT INFO FROM variants WHERE INFO LIKE '%gene_fusion%'""")
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
        param = {"annotation": {"snpeff": {"options": "-lof -hgvs -oicr -noShiftHgvs -spliceSiteSize 3 "}}}

        # Create object
        variants = Variants(conn=None, input=input_vcf, output=output_vcf, config=tests_config, param=param, load=True)

        # Remove if output file exists
        remove_if_exists([output_vcf])

        # Annotation
        variants.annotation()

        # query annotated variant
        result = variants.get_query_to_df(""" SELECT * FROM variants WHERE "#CHROM" = 'chr12' AND POS = 68724951 AND REF = 'G' AND ALT = 'T' AND INFO LIKE '%T|synonymous_variant|LOW|MDM1|MDM1|transcript|NM_001354969.1|protein_coding|2/15|c.69C>A|p.Ser23Ser|238/3032|69/2175|23/724||%' """)
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
        param = {"annotations": {
                    f"{annotation_snpeff}": None
                    }
        }

        # Create object
        variants = Variants(conn=None, input=input_vcf, output=output_vcf, config=tests_config, param=param, load=True)

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
        param = {"annotation": {"snpeff": {"options": "-lof -hgvs -oicr -noShiftHgvs -spliceSiteSize 3 "}}}

        # Create object
        variants = Variants(conn=None, input=input_vcf, output=output_vcf, config=config, param=param, load=True)

        # Remove if output file exists
        remove_if_exists([output_vcf])

        # Annotation
        variants.annotation()

        # query annotated variant
        result = variants.get_query_to_df(""" SELECT INFO FROM variants WHERE INFO LIKE '%ANN=%' """)
        assert len(result) == 7

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
        param = {"annotation": {"bcftools": {"annotations":  {annotation_parquet: {"INFO": None}}}}}

        # Create object
        variants = Variants(conn=None, input=input_vcf, output=output_vcf, param=param, config=config, load=True)

        # Remove if output file exists
        remove_if_exists([output_vcf])

        # Annotation
        variants.annotation()

        result = variants.get_query_to_df("""SELECT * FROM variants WHERE "#CHROM" = 'chr7' AND POS = 55249063 AND REF = 'G' AND ALT = 'A' AND INFO = 'DP=125;nci60=0.66'""")
        assert len(result) == 1

        # Check if VCF is in correct format with pyVCF
        variants.export_output()
        try:
            vcf.Reader(filename=output_vcf)
        except:
            assert False
    
    
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

