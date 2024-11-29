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


def test_calculation_sql_on_table():
    """
    The function `test_calculation_sql_on_table` performs various calculations and checks on a table of
    genetic variants using SQL operations.
    """

    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        tmp_dir = "/tmp"

        # Init files
        input_vcf = tests_data_folder + "/example.vcf.gz"
        output_vcf = f"{tmp_dir}/output.vcf.gz"

        # Create object
        variants = Variants(input=input_vcf, output=output_vcf, load=True)

        # Operations config file
        operations_config_file = os.path.join(
            tests_data_folder, "operations_config_test.json"
        )

        # New table
        new_table = "transcripts"

        # Operations config dict
        operations_config_dict = {
            "variant_chr_pos_alt_ref_dict": {
                "type": "sql",
                "name": "variant_chr_pos_alt_ref_dict",
                "description": "Create a variant ID with chromosome, position, alt and ref",
                "available": False,
                "table": new_table,
                "output_column_name": "variant_chr_pos_alt_ref_dict",
                "output_column_type": "String",
                "output_column_description": "variant ID with chromosome, position, alt and ref",
                "operation_query": """ concat("#CHROM", '_', "POS", '_', "REF", '_', "ALT") """,
                "operation_info": True,
            }
        }

        # Operations
        operations = {
            "variant_chr_pos_alt_ref": None,
            "variant_chr_pos_alt_ref_file": None,
            "variant_chr_pos_alt_ref_dict": None,
        }

        # Create another table
        query_create_table = f"""
            CREATE TABLE {new_table} AS SELECT * FROM variants
        """
        variants.execute_query(query=query_create_table)

        # Calculation
        variants.calculation(
            operations=operations,
            operations_config_dict=operations_config_dict,
            operations_config_file=operations_config_file,
        )

        # Check number of variant_chr_pos_alt_ref
        result = variants.get_query_to_df(
            f"""SELECT INFO FROM variants WHERE INFO LIKE '%variant_chr_pos_alt_ref=%' """
        )
        assert len(result) == 7

        # Check number of variant_chr_pos_alt_ref_dict with new table
        result = variants.get_query_to_df(
            f"""SELECT INFO FROM {new_table} WHERE INFO LIKE '%variant_chr_pos_alt_ref_dict=%' """
        )
        assert len(result) == 7

        # Check number of variant_chr_pos_alt_ref_file
        result = variants.get_query_to_df(
            f"""SELECT INFO FROM variants WHERE INFO LIKE '%variant_chr_pos_alt_ref_file=%' """
        )
        assert len(result) == 7

        # Check number of middle (7)
        result = variants.get_query_to_df(
            f"""SELECT INFO FROM variants WHERE INFO LIKE '%variant_chr_pos_alt_ref=chr1_28736_A_C%' """
        )
        assert len(result) == 1

        # Check if VCF is in correct format with pyVCF
        remove_if_exists([output_vcf])
        variants.export_output()
        try:
            vcf.Reader(filename=output_vcf)
        except:
            assert False


def test_calculation_sql():
    """
    This function tests the calculation and annotation of genetic variants using input parameters and
    checks if the output VCF file is in the correct format.
    """

    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        tmp_dir = "/tmp"

        # Init files
        input_vcf = tests_data_folder + "/example.vcf.gz"
        output_vcf = f"{tmp_dir}/output.vcf.gz"

        # Create object
        variants = Variants(input=input_vcf, output=output_vcf, load=True)

        # Operations config file
        operations_config_file = os.path.join(
            tests_data_folder, "operations_config_test.json"
        )

        # Operations config dict
        operations_config_dict = {
            "variant_chr_pos_alt_ref_dict": {
                "type": "sql",
                "name": "variant_chr_pos_alt_ref_dict",
                "description": "Create a variant ID with chromosome, position, alt and ref",
                "available": False,
                "output_column_name": "variant_chr_pos_alt_ref_dict",
                "output_column_type": "String",
                "output_column_description": "variant ID with chromosome, position, alt and ref",
                "operation_query": """ concat("#CHROM", '_', "POS", '_', "REF", '_', "ALT") """,
                "operation_info": True,
            }
        }

        # Operations
        operations = {
            "variant_chr_pos_alt_ref": None,
            "variant_chr_pos_alt_ref_file": None,
            "variant_chr_pos_alt_ref_dict": None,
        }

        # Calculation
        variants.calculation(
            operations=operations,
            operations_config_dict=operations_config_dict,
            operations_config_file=operations_config_file,
        )

        # Check number of variant_chr_pos_alt_ref
        result = variants.get_query_to_df(
            """SELECT INFO FROM variants WHERE INFO LIKE '%variant_chr_pos_alt_ref=%' """
        )
        assert len(result) == 7

        # Check number of variant_chr_pos_alt_ref_dict
        result = variants.get_query_to_df(
            """SELECT INFO FROM variants WHERE INFO LIKE '%variant_chr_pos_alt_ref_dict=%' """
        )
        assert len(result) == 7

        # Check number of variant_chr_pos_alt_ref_file
        result = variants.get_query_to_df(
            """SELECT INFO FROM variants WHERE INFO LIKE '%variant_chr_pos_alt_ref_file=%' """
        )
        assert len(result) == 7

        # Check number of middle (7)
        result = variants.get_query_to_df(
            """SELECT INFO FROM variants WHERE INFO LIKE '%variant_chr_pos_alt_ref=chr1_28736_A_C%' """
        )
        assert len(result) == 1

        # Check if VCF is in correct format with pyVCF
        remove_if_exists([output_vcf])
        variants.export_output()
        try:
            vcf.Reader(filename=output_vcf)
        except:
            assert False


def test_calculation_sql_config_file_in_param():
    """
    This function tests the calculation and annotation of genetic variants using input parameters and
    checks if the output VCF file is in the correct format.
    """

    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        tmp_dir = "/tmp"

        # Init files
        input_vcf = tests_data_folder + "/example.vcf.gz"
        output_vcf = f"{tmp_dir}/output.vcf.gz"

        # Create object
        variants = Variants(input=input_vcf, output=output_vcf, load=True)

        # Operations config file
        operations_config_file = os.path.join(
            tests_data_folder, "operations_config_test.json"
        )

        # Load param
        param = {
            "calculation": {
                "calculation_config" : operations_config_file
            }
        }
        variants.set_param(param=param)

        # Operations config dict
        operations_config_dict = {
            "variant_chr_pos_alt_ref_dict": {
                "type": "sql",
                "name": "variant_chr_pos_alt_ref_dict",
                "description": "Create a variant ID with chromosome, position, alt and ref",
                "available": False,
                "output_column_name": "variant_chr_pos_alt_ref_dict",
                "output_column_type": "String",
                "output_column_description": "variant ID with chromosome, position, alt and ref",
                "operation_query": """ concat("#CHROM", '_', "POS", '_', "REF", '_', "ALT") """,
                "operation_info": True,
            }
        }

        # Operations
        operations = {
            "variant_chr_pos_alt_ref": None,
            "variant_chr_pos_alt_ref_file": None,
            "variant_chr_pos_alt_ref_dict": None,
        }

        # Calculation
        variants.calculation(
            operations=operations,
            operations_config_dict=operations_config_dict,
        )

        # Check number of variant_chr_pos_alt_ref
        result = variants.get_query_to_df(
            """SELECT INFO FROM variants WHERE INFO LIKE '%variant_chr_pos_alt_ref=%' """
        )
        assert len(result) == 7

        # Check number of variant_chr_pos_alt_ref_dict
        result = variants.get_query_to_df(
            """SELECT INFO FROM variants WHERE INFO LIKE '%variant_chr_pos_alt_ref_dict=%' """
        )
        assert len(result) == 7

        # Check number of variant_chr_pos_alt_ref_file
        result = variants.get_query_to_df(
            """SELECT INFO FROM variants WHERE INFO LIKE '%variant_chr_pos_alt_ref_file=%' """
        )
        assert len(result) == 7

        # Check number of middle (7)
        result = variants.get_query_to_df(
            """SELECT INFO FROM variants WHERE INFO LIKE '%variant_chr_pos_alt_ref=chr1_28736_A_C%' """
        )
        assert len(result) == 1

        # Check if VCF is in correct format with pyVCF
        remove_if_exists([output_vcf])
        variants.export_output()
        try:
            vcf.Reader(filename=output_vcf)
        except:
            assert False


def test_calculation_sql_fail():
    """
    This function tests the calculation and annotation of genetic variants using input parameters and
    checks if the output VCF file is in the correct format.
    """

    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        tmp_dir = "/tmp"

        # Init files
        input_vcf = tests_data_folder + "/example.vcf.gz"
        output_vcf = f"{tmp_dir}/output.vcf.gz"

        # Create object
        variants = Variants(input=input_vcf, output=output_vcf, load=True)

        # Operations config dict
        operations_config_dict = {
            "QUERY_FAILED": {
                "type": "sql",
                "name": "QUERY_FAILED",
                "description": "Variant type (e.g. SNV, INDEL, MNV, BND...)",
                "available": True,
                "output_column_name": "QUERY_FAILED",
                "output_column_type": "String",
                "output_column_description": "Variant type: SNV if X>Y, MOSAIC if X>Y,Z or X,Y>Z, INDEL if XY>Z or X>YZ",
                "operation_query": "blabla",
                "info_fields": ["FAILED"],
                "operation_info": True,
            }
        }

        # Operations
        operations = {"QUERY_FAILED": None}

        # Calculation
        with pytest.raises(ValueError) as e:
            variants.calculation(
                operations=operations, operations_config_dict=operations_config_dict
            )
        assert (
            str(e.value) == "Operations config: Calculation 'QUERY_FAILED' query failed"
        )


def test_calculation_sql_info_fields_check():
    """
    This function tests the calculation and annotation of genetic variants using input parameters and
    checks if the output VCF file is in the correct format.
    """

    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        tmp_dir = "/tmp"

        # Init files
        input_vcf = tests_data_folder + "/example.vcf.gz"
        output_vcf = f"{tmp_dir}/output.vcf.gz"

        # Create object
        variants = Variants(input=input_vcf, output=output_vcf, load=True)

        # Operations config dict
        operations_config_dict = {
            "variant_chr_pos_alt_ref_dict": {
                "type": "sql",
                "name": "variant_chr_pos_alt_ref_dict",
                "description": "Create a variant ID with chromosome, position, alt and ref",
                "available": False,
                "output_column_name": "variant_chr_pos_alt_ref_dict",
                "output_column_type": "String",
                "output_column_description": "variant ID with chromosome, position, alt and ref",
                "operation_query": """ concat("#CHROM", '_', "POS", '_', "REF", '_', "ALT") """,
                "operation_info": True,
                "info_fields": ["SVTYPE"],
                "info_fields_check": False,
            }
        }

        # Operations
        operations = {"variant_chr_pos_alt_ref_dict": None}

        # Calculation
        variants.calculation(
            operations=operations, operations_config_dict=operations_config_dict
        )

        # Check number of variant_chr_pos_alt_ref
        result = variants.get_query_to_df(
            """SELECT INFO FROM variants WHERE INFO LIKE '%variant_chr_pos_alt_ref_dict=%' """
        )
        assert len(result) == 7

        # Operations config dict
        operations_config_dict = {
            "variant_chr_pos_alt_ref_dict": {
                "type": "sql",
                "name": "variant_chr_pos_alt_ref_dict",
                "description": "Create a variant ID with chromosome, position, alt and ref",
                "available": False,
                "output_column_name": "variant_chr_pos_alt_ref_dict",
                "output_column_type": "String",
                "output_column_description": "variant ID with chromosome, position, alt and ref",
                "operation_query": """ concat("#CHROM", '_', "POS", '_', "REF", '_', "ALT") """,
                "operation_info": True,
                "info_fields": ["SVTYPE"],
                "info_fields_check": True,
            }
        }

        # Operations
        operations = {"variant_chr_pos_alt_ref_dict": None}

        # Calculation
        with pytest.raises(ValueError) as e:
            variants.calculation(
                operations=operations, operations_config_dict=operations_config_dict
            )
        assert (
            str(e.value)
            == "Operations config: Calculation 'variant_chr_pos_alt_ref_dict' DOES NOT contain all mandatory fields ['SVTYPE']"
        )


def test_calculation_nomen():
    """
    This function tests the calculation and annotation of genetic variants using input parameters and
    checks if the output VCF file is in the correct format.
    """

    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # Init files
        input_vcf = tests_data_folder + "/example.vcf.gz"
        output_vcf = f"{tmp_dir}/output.vcf.gz"
        input_param = {
            "annotation": {
                "annovar": {
                    "annotations": {
                        "refGene": {
                            "Func_refGene": "location",
                            "Gene_refGene": "gene",
                            "GeneDetail_refGene": "GeneDetail",
                            "ExonicFunc_refGene": "outcome",
                            "AAChange_refGene": "hgvs",
                        }
                    },
                    "options": {
                        "genebase": "-hgvs -splicing_threshold 3 ",
                        "intronhgvs": 10,
                    },
                }
            },
            "calculation": {
                "calculations": {"NOMEN": {"options": {"hgvs_field": "hgvs"}}}
            },
        }

        # Create object
        variants = Variants(
            input=input_vcf,
            output=output_vcf,
            config=tests_config,
            param=input_param,
            load=True,
        )

        # Annotation
        variants.annotation()

        # Calculation
        variants.calculation()

        # Check number of NOMEN (2)
        result = variants.get_query_to_df(
            """SELECT INFO FROM variants WHERE INFO LIKE '%NOMEN=%' """
        )
        assert len(result) == 2

        # Check if VCF is in correct format with pyVCF
        remove_if_exists([output_vcf])
        variants.export_output()
        try:
            vcf.Reader(filename=output_vcf)
        except:
            assert False


def test_calculation_vartype():
    """
    This function tests the calculation of variant types in a VCF file using the Variants class.
    """

    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # Init files
        input_vcf = tests_data_folder + "/example.snv.indel.mosaic.vcf"
        output_vcf = f"{tmp_dir}/output.vcf.gz"

        # Construct param dict
        param = {"calculation": {"calculations": {"VARTYPE": None}}}

        # Create object
        variants = Variants(
            conn=None, input=input_vcf, output=output_vcf, param=param, load=True
        )

        # Calculation
        variants.calculation()

        result = variants.get_query_to_df(
            """ SELECT * FROM variants WHERE INFO LIKE '%VARTYPE=SNV%' """
        )
        assert len(result) == 5

        result = variants.get_query_to_df(
            """ SELECT * FROM variants WHERE INFO LIKE '%VARTYPE=INDEL%' """
        )
        assert len(result) == 1

        result = variants.get_query_to_df(
            """ SELECT * FROM variants WHERE INFO LIKE '%VARTYPE=MOSAIC%' """
        )
        assert len(result) == 1

        # Check if VCF is in correct format with pyVCF
        remove_if_exists([output_vcf])
        variants.export_output()
        try:
            vcf.Reader(filename=output_vcf)
        except:
            assert False


def test_calculation_vartype_full():
    """
    This function tests the calculation of variant types in a VCF file using the Variants class.
    """

    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # Init files
        input_vcf = tests_data_folder + "/example.full.unsorted.vcf.gz"
        output_vcf = f"{tmp_dir}/output.vcf.gz"

        # Construct param dict
        param = {"calculation": {"calculations": {"VARTYPE": None}}}

        # Create object
        variants = Variants(
            conn=None, input=input_vcf, output=output_vcf, param=param, load=True
        )

        # Calculation
        variants.calculation()

        result = variants.get_query_to_df(
            """ SELECT * FROM variants WHERE INFO LIKE '%VARTYPE=SNV%' """
        )
        assert len(result) == 4

        result = variants.get_query_to_df(
            """ SELECT * FROM variants WHERE INFO LIKE '%VARTYPE=INDEL%' """
        )
        assert len(result) == 2

        result = variants.get_query_to_df(
            """ SELECT * FROM variants WHERE INFO LIKE '%VARTYPE=CNV%' """
        )
        assert len(result) == 3

        result = variants.get_query_to_df(
            """ SELECT * FROM variants WHERE INFO LIKE '%VARTYPE=INV%' """
        )
        assert len(result) == 3

        result = variants.get_query_to_df(
            """ SELECT * FROM variants WHERE INFO LIKE '%VARTYPE=DEL%' """
        )
        assert len(result) == 3

        result = variants.get_query_to_df(
            """ SELECT * FROM variants WHERE INFO LIKE '%VARTYPE=INS%' """
        )
        assert len(result) == 5

        result = variants.get_query_to_df(
            """ SELECT * FROM variants WHERE INFO LIKE '%VARTYPE=DUP%' """
        )
        assert len(result) == 6

        result = variants.get_query_to_df(
            """ SELECT * FROM variants WHERE INFO LIKE '%VARTYPE=BND%' """
        )
        assert len(result) == 7

        result = variants.get_query_to_df(
            """ SELECT * FROM variants WHERE INFO LIKE '%VARTYPE=MNV%' """
        )
        assert len(result) == 1

        # Check if VCF is in correct format with pyVCF
        remove_if_exists([output_vcf])
        variants.export_output()
        try:
            vcf.Reader(filename=output_vcf)
        except:
            assert False


@pytest.mark.parametrize(
    "calculation",
    ["snpeff_ann_explode", "snpeff_ann_explode_uniquify", "snpeff_ann_explode_json"],
)
def test_calculation_snpeff_ann_explode(calculation):
    """
    This function is a test for calculating snpeff_hgvs in a VCF file using the Variants class.

    :param calculation: It looks like the `calculation` parameter is used to specify a particular
    calculation method in the test function `test_calculation_snpeff_ann_explode`. This parameter is
    then used to construct a parameter dictionary `param` with the calculation method specified
    """

    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # Init files
        input_vcf = tests_data_folder + "/example.ann.vcf.gz"
        output_vcf = f"{tmp_dir}/output.{calculation}.vcf"

        # Construct param dict
        param = {"calculation": {"calculations": {calculation: None}}}

        # Create object
        variants = Variants(
            conn=None, input=input_vcf, output=output_vcf, param=param, load=True
        )

        # Check if no snpeff_hgvs
        result = variants.get_query_to_df(
            """ SELECT INFO FROM variants WHERE regexp_matches(INFO,'snpeff[^=]') """
        )
        assert len(result) == 0

        # Calculation
        variants.calculation()

        # query annotated variant
        result = variants.get_query_to_df(
            """ SELECT * FROM variants WHERE regexp_matches(INFO,'snpeff[^=]') """
        )
        assert len(result) == 7

        # Check if VCF is in correct format with pyVCF
        remove_if_exists([output_vcf])
        variants.export_output()
        try:
            vcf.Reader(filename=output_vcf)
        except:
            assert False


def test_calculation_snpeff_hgvs():
    """
    This function is a test for calculating snpeff_hgvs in a VCF file using the Variants class.
    """

    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # Init files
        input_vcf = tests_data_folder + "/example.ann.vcf.gz"
        output_vcf = f"{tmp_dir}/output.vcf.gz"

        # Construct param dict
        param = {"calculation": {"calculations": {"snpeff_hgvs": None}}}

        # Create object
        variants = Variants(
            conn=None, input=input_vcf, output=output_vcf, param=param, load=True
        )

        # Check if no snpeff_hgvs
        result = variants.get_query_to_df(
            """ SELECT INFO FROM variants WHERE INFO LIKE '%snpeff_hgvs=%' """
        )
        assert len(result) == 0

        # Calculation
        variants.calculation()

        # query annotated variant
        result = variants.get_query_to_df(
            """ SELECT * FROM variants WHERE INFO LIKE '%snpeff_hgvs=%' """
        )
        assert len(result) == 7

        # Check if VCF is in correct format with pyVCF
        remove_if_exists([output_vcf])
        variants.export_output()
        try:
            vcf.Reader(filename=output_vcf)
        except:
            assert False


def test_calculation_snpeff_hgvs_no_ann():
    """
    This function tests the calculation of SNPEff HGVS annotations on a VCF file with no annotations.
    """

    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # Init files
        input_vcf = tests_data_folder + "/example.vcf.gz"
        output_vcf = f"{tmp_dir}/output.vcf.gz"

        # Construct param dict
        param = {"calculation": {"calculations": {"snpeff_hgvs": None}}}

        # Create object
        variants = Variants(
            conn=None, input=input_vcf, output=output_vcf, param=param, load=True
        )

        # Calculation
        variants.calculation()

        # query annotated variant
        result = variants.get_query_to_df(
            """ SELECT * FROM variants WHERE INFO LIKE '%snpeff_hgvs=%' """
        )
        assert len(result) == 0

        # Check if VCF is in correct format with pyVCF
        remove_if_exists([output_vcf])
        variants.export_output()
        try:
            vcf.Reader(filename=output_vcf)
        except:
            assert False


def test_calculation_snpeff_hgvs_transcripts():
    """
    This function tests the calculation of SNPEff HGVS transcripts using a VCF file and a transcripts
    file.
    """

    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # Init files
        input_vcf = tests_data_folder + "/example.snpeff.vcf.gz"
        transcripts_file = tests_data_folder + "/transcripts.tsv"
        output_vcf = f"{tmp_dir}/output.vcf.gz"

        # Construct param dict
        param = {
            "calculation": {
                "calculations": {
                    "NOMEN": {
                        "options": {
                            "hgvs_field": "snpeff_hgvs",
                            "transcripts": transcripts_file,
                        }
                    }
                }
            }
        }

        # Create object
        variants = Variants(
            conn=None, input=input_vcf, output=output_vcf, param=param, load=True
        )

        # Calculation
        variants.calculation()

        # query annotated variant
        result = variants.get_query_to_df(
            """ SELECT * FROM variants WHERE INFO LIKE '%;NOMEN=%' """
        )
        assert len(result) == 7

        # query annotated variant
        result = variants.get_query_to_df(
            """ SELECT * FROM variants WHERE INFO LIKE '%;TNOMEN=%' """
        )
        assert len(result) == 7

        # Check transcript priority
        result = variants.get_query_to_df(
            """ SELECT * FROM variants WHERE INFO LIKE '%;NOMEN=EGFR:NM_001346897%' """
        )
        assert len(result) == 1

        # Check if VCF is in correct format with pyVCF
        remove_if_exists([output_vcf])
        variants.export_output()
        try:
            vcf.Reader(filename=output_vcf)
        except:
            assert False


def test_calculation_snpeff_hgvs_notranscripts():
    """
    This function tests the calculation of SNPEff HGVS notranscripts in a VCF file.
    """

    # Init files
    input_vcf = tests_data_folder + "/example.snpeff.vcf.gz"
    transcripts_file = tests_data_folder + "/notranscripts.tsv"

    # Construct param dict
    param = {
        "calculation": {
            "calculations": {
                "NOMEN": {
                    "options": {
                        "hgvs_field": "snpeff_hgvs",
                        "transcripts": transcripts_file,
                    }
                }
            }
        }
    }

    # Create object
    variants = Variants(conn=None, input=input_vcf, param=param, load=True)

    # Calculation
    with pytest.raises(ValueError) as e:
        variants.calculation()
    assert str(e.value) == f"Transcript file '{transcripts_file}' does NOT exist"


def test_calculation_findbypipeline():
    """
    This is a test function for the "FINDBYPIPELINE" calculation in the Variants class, which checks if
    the calculation is performed correctly and the output VCF file is in the correct format.
    """

    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # Init files
        input_vcf = tests_data_folder + "/example.vcf.gz"
        output_vcf = f"{tmp_dir}/output.vcf.gz"

        # Construct param dict
        param = {"calculation": {"calculations": {"FINDBYPIPELINE": None}}}

        # Create object
        variants = Variants(
            conn=None, input=input_vcf, output=output_vcf, param=param, load=True
        )

        # Calculation
        variants.calculation()

        result = variants.get_query_to_df(
            """ SELECT INFO FROM variants WHERE INFO LIKE '%findbypipeline%' """
        )
        assert len(result) == 7

        result = variants.get_query_to_df(
            """ SELECT * FROM variants WHERE INFO LIKE '%findbypipeline=4/4%' """
        )
        assert len(result) == 1

        result = variants.get_query_to_df(
            """ SELECT * FROM variants WHERE INFO LIKE '%findbypipeline=3/4%' """
        )
        assert len(result) == 6

        # Check if VCF is in correct format with pyVCF
        remove_if_exists([output_vcf])
        variants.export_output()
        try:
            vcf.Reader(filename=output_vcf)
        except:
            assert False


def test_calculation_findbysample():
    """
    This is a test function for the "FINDBYPIPELINE" calculation in the Variants class, which checks if
    the calculation is performed correctly and the output VCF file is in the correct format.
    """

    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # Init files
        input_vcf = tests_data_folder + "/example.vcf.gz"
        output_vcf = f"{tmp_dir}/output.vcf.gz"

        # Construct param dict
        param = {"calculation": {"calculations": {"FINDBYSAMPLE": None}}}

        # Create object
        variants = Variants(
            conn=None, input=input_vcf, output=output_vcf, param=param, load=True
        )

        # Calculation
        variants.calculation()

        result = variants.get_query_to_df(
            """ SELECT INFO FROM variants WHERE INFO LIKE '%findbysample%' """
        )
        assert len(result) == 7

        result = variants.get_query_to_df(
            """ SELECT * FROM variants WHERE INFO LIKE '%findbysample=4/4%' """
        )
        assert len(result) == 1

        result = variants.get_query_to_df(
            """ SELECT * FROM variants WHERE INFO LIKE '%findbysample=3/4%' """
        )
        assert len(result) == 6

        # Check if VCF is in correct format with pyVCF
        remove_if_exists([output_vcf])
        variants.export_output()
        try:
            vcf.Reader(filename=output_vcf)
        except:
            assert False


def test_calculation_genotype_concordance():
    """
    This is a test function for calculating genotype concordance in a VCF file using the Variants class.
    """

    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # Init files
        input_vcf = tests_data_folder + "/example.vcf.gz"
        output_vcf = f"{tmp_dir}/output.vcf.gz"

        # Construct param dict
        param = {"calculation": {"calculations": {"GENOTYPECONCORDANCE": None}}}

        # Create object
        variants = Variants(
            conn=None, input=input_vcf, output=output_vcf, param=param, load=True
        )

        # Calculation
        variants.calculation()

        result = variants.get_query_to_df(
            """ SELECT INFO FROM variants WHERE INFO LIKE '%genotypeconcordance%' """
        )
        assert len(result) == 7

        result = variants.get_query_to_df(
            """ SELECT * FROM variants WHERE INFO LIKE '%genotypeconcordance=FALSE%' """
        )
        assert len(result) == 1

        result = variants.get_query_to_df(
            """ SELECT * FROM variants WHERE INFO LIKE '%genotypeconcordance=TRUE%' """
        )
        assert len(result) == 6

        # Check if VCF is in correct format with pyVCF
        remove_if_exists([output_vcf])
        variants.export_output()
        try:
            vcf.Reader(filename=output_vcf)
        except:
            assert False


def test_calculation_barcode():
    """
    This is a test function for a Python script that calculates barcode information from a VCF file and
    checks if the output is correct.
    """

    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # Init files
        input_vcf = tests_data_folder + "/example.vcf.gz"
        output_vcf = f"{tmp_dir}/output.vcf.gz"

        # Construct param dict
        param = {"calculation": {"calculations": {"BARCODE": None}}}

        # Create object
        variants = Variants(
            conn=None, input=input_vcf, output=output_vcf, param=param, load=True
        )

        # Calculation
        variants.calculation()

        result = variants.get_query_to_df(
            """ SELECT INFO FROM variants WHERE INFO LIKE '%barcode%' """
        )
        assert len(result) == 7

        result = variants.get_query_to_df(
            """ SELECT * FROM variants WHERE INFO LIKE '%barcode=1122%' """
        )
        assert len(result) == 1

        result = variants.get_query_to_df(
            """ SELECT * FROM variants WHERE INFO LIKE '%barcode=0111%' """
        )
        assert len(result) == 1

        result = variants.get_query_to_df(
            """ SELECT * FROM variants WHERE INFO LIKE '%barcode=1011%' """
        )
        assert len(result) == 4

        result = variants.get_query_to_df(
            """ SELECT * FROM variants WHERE INFO LIKE '%barcode=1101%' """
        )
        assert len(result) == 1

        # Check if VCF is in correct format with pyVCF
        remove_if_exists([output_vcf])
        variants.export_output()
        try:
            vcf.Reader(filename=output_vcf)
        except:
            assert False


def test_calculation_barcode_genotype():
    """
    The function `test_calculation_barcode_genotype` is a test function in Python that calculates
    barcode information from a VCF file and checks if the output is correct.
    """

    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        tmp_dir = "/tmp"

        # Init files
        input_vcf = tests_data_folder + "/example.vcf.gz"
        output_vcf = f"{tmp_dir}/output.vcf.gz"

        # Construct param dict
        param = {
            "calculation": {"calculations": {"BARCODEFAMILY": {"family_pedigree": ""}}}
        }

        # Create object
        variants = Variants(
            conn=None, input=input_vcf, output=output_vcf, param=param, load=True
        )

        # Calculation
        variants.calculation()

        # DEBUG
        result = variants.get_query_to_df(""" SELECT * FROM variants """)
        log.debug(result)

        # Construct param dict
        params = {
            "param_str_list": {
                "calculation": {
                    "calculations": {
                        "BARCODEFAMILY": {"family_pedigree": "sample1,sample2,sample4"}
                    }
                }
            },
            "param_str_json": {
                "calculation": {
                    "calculations": {
                        "BARCODEFAMILY": {
                            "family_pedigree": """{
                                "father": "sample1", "mother": "sample2", "child": "sample4"}"""
                        }
                    }
                }
            },
            "param_dict": {
                "calculation": {
                    "calculations": {
                        "BARCODEFAMILY": {
                            "family_pedigree": {
                                "father": "sample1",
                                "mother": "sample2",
                                "child": "sample4",
                            }
                        }
                    }
                }
            },
            "param_file": {
                "calculation": {
                    "calculations": {
                        "BARCODEFAMILY": {
                            "family_pedigree": os.path.join(
                                tests_data_folder, "trio.json"
                            )
                        }
                    }
                }
            },
        }

        for param in params:
            param = params[param]

            # Create object
            variants = Variants(
                conn=None, input=input_vcf, output=output_vcf, param=param, load=True
            )

            # Calculation
            variants.calculation()

            # Check if BCF and BCFS are in FORMAT
            result = variants.get_query_to_df(
                """ SELECT FORMAT FROM variants WHERE FORMAT LIKE '%:BCF:BCFS' """
            )
            assert len(result) == 7

            # Check if VCF is in correct format with pyVCF
            remove_if_exists([output_vcf])
            variants.export_output()
            try:
                vcf.Reader(filename=output_vcf)
            except:
                assert False


def test_calculation_trio():
    """
    This is a test function for the calculation of trio variants in a VCF file using specific
    parameters.
    """

    params = {
        "param_str_list": {
            "calculation": {
                "calculations": {"TRIO": {"trio_pedigree": "sample1,sample2,sample4"}}
            }
        },
        "param_str_json": {
            "calculation": {
                "calculations": {
                    "TRIO": {
                        "trio_pedigree": """{
                            "father": "sample1", "mother": "sample2", "child": "sample4"}"""
                    }
                }
            }
        },
        "param_dict": {
            "calculation": {
                "calculations": {
                    "TRIO": {
                        "trio_pedigree": {
                            "father": "sample1",
                            "mother": "sample2",
                            "child": "sample4",
                        }
                    }
                }
            }
        },
        "param_file": {
            "calculation": {
                "calculations": {
                    "TRIO": {
                        "trio_pedigree": os.path.join(tests_data_folder, "trio.json")
                    }
                }
            }
        },
    }

    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # Init files
        input_vcf = tests_data_folder + "/example.vcf.gz"
        output_vcf = f"{tmp_dir}/output.vcf.gz"

        # Param NO
        param = {"calculation": {"calculations": {"TRIO": {"trio_pedigree": None}}}}

        # Create object
        variants = Variants(
            conn=None, input=input_vcf, output=output_vcf, param=param, load=True
        )

        # Calculation
        variants.calculation()

        result = variants.get_query_to_df(
            """ SELECT INFO FROM variants WHERE INFO LIKE '%trio=recessive%' """
        )
        assert len(result) == 1

        result = variants.get_query_to_df(
            """ SELECT * FROM variants WHERE INFO LIKE '%trio=dominant%' """
        )
        assert len(result) == 5

        result = variants.get_query_to_df(
            """ SELECT * FROM variants WHERE INFO LIKE '%trio=unknown%' """
        )
        assert len(result) == 1

        # Check if VCF is in correct format with pyVCF
        remove_if_exists([output_vcf])
        variants.export_output()
        try:
            vcf.Reader(filename=output_vcf)
        except:
            assert False

        # Construct param dict
        for param_id in params:

            # Param
            param = params.get(param_id)

            # Create object
            variants = Variants(
                conn=None, input=input_vcf, output=output_vcf, param=param, load=True
            )

            # Calculation
            variants.calculation()

            result = variants.get_query_to_df(
                """ SELECT INFO FROM variants WHERE INFO LIKE '%trio=recessive%' """
            )
            assert len(result) == 1

            result = variants.get_query_to_df(
                """ SELECT * FROM variants WHERE INFO LIKE '%trio=dominant%' """
            )
            assert len(result) == 6

            result = variants.get_query_to_df(
                """ SELECT * FROM variants WHERE INFO LIKE '%trio=unknown%' """
            )
            assert len(result) == 0

            # Check if VCF is in correct format with pyVCF
            remove_if_exists([output_vcf])
            variants.export_output()
            try:
                vcf.Reader(filename=output_vcf)
            except:
                assert False


def test_calculation_vaf_normalization():
    """
    This is a test function for the calculation of variant allele frequency normalization in a VCF file.
    """

    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # Init files
        input_vcf = tests_data_folder + "/example.vcf.gz"
        output_vcf = f"{tmp_dir}/output.vcf.gz"

        # Construct param dict
        param = {"calculation": {"calculations": {"vaf": None}}}

        # Create object
        variants = Variants(
            conn=None, input=input_vcf, output=output_vcf, param=param, load=True
        )

        # Calculation
        variants.calculation()

        result = variants.get_query_to_df(
            """ SELECT INFO FROM variants WHERE FORMAT LIKE '%:VAF' """
        )
        assert len(result) == 7

        result = variants.get_query_to_df(
            """ SELECT * FROM variants WHERE "#CHROM" = 'chr1' AND POS = 28736 AND sample1 LIKE '%:0.279835' """
        )
        assert len(result) == 1

        result = variants.get_query_to_df(
            """ SELECT * FROM variants WHERE "#CHROM" = 'chr1' AND POS = 28736 AND sample2 LIKE '%:0.282898' """
        )
        assert len(result) == 1

        result = variants.get_query_to_df(
            """ SELECT * FROM variants WHERE "#CHROM" = 'chr1' AND POS = 28736 AND sample3 LIKE '%:0.282955' """
        )
        assert len(result) == 1

        result = variants.get_query_to_df(
            """ SELECT * FROM variants WHERE "#CHROM" = 'chr1' AND POS = 28736 AND sample4 LIKE '%:0.303819' """
        )
        assert len(result) == 1

        # Check if VCF is in correct format with pyVCF
        remove_if_exists([output_vcf])
        variants.export_output()
        try:
            vcf.Reader(filename=output_vcf)
        except:
            assert False


def test_calculation_vaf_stats():
    """
    This is a test function for the calculation of variant allele frequency (VAF) statistics in a VCF
    file using the Variants class in Python.
    """

    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # Init files
        input_vcf = tests_data_folder + "/example.vcf.gz"
        output_vcf = f"{tmp_dir}/output.vcf.gz"

        # Construct param dict
        param = {"calculation": {"calculations": {"vaf": None, "vaf_stats": None}}}

        # Create object
        variants = Variants(
            conn=None, input=input_vcf, output=output_vcf, param=param, load=True
        )

        # Calculation
        variants.calculation()

        result = variants.get_query_to_df(
            """ SELECT INFO FROM variants WHERE INFO LIKE '%VAF_stats%' """
        )
        assert len(result) == 7

        result = variants.get_query_to_df(
            """ SELECT * FROM variants WHERE "#CHROM" = 'chr1' AND POS = 28736 AND INFO LIKE '%VAF_stats_nb=4%' """
        )
        assert len(result) == 1

        result = variants.get_query_to_df(
            """ SELECT * FROM variants WHERE "#CHROM" = 'chr1' AND POS = 28736 AND INFO LIKE '%VAF_stats_min=0.279835%' """
        )
        assert len(result) == 1

        result = variants.get_query_to_df(
            """ SELECT * FROM variants WHERE "#CHROM" = 'chr1' AND POS = 28736 AND INFO LIKE '%VAF_stats_max=0.303819%' """
        )
        assert len(result) == 1

        result = variants.get_query_to_df(
            """ SELECT * FROM variants WHERE "#CHROM" = 'chr1' AND POS = 28736 AND INFO LIKE '%VAF_stats_mean=0.28737675%' """
        )
        assert len(result) == 1

        # Check if VCF is in correct format with pyVCF
        remove_if_exists([output_vcf])
        variants.export_output()
        try:
            vcf.Reader(filename=output_vcf)
        except:
            assert False


def test_calculation_dp_stats():
    """
    This is a test function for the calculation of depth statistics in a VCF file using the Variants
    class in Python.
    """

    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # Init files
        input_vcf = tests_data_folder + "/example.vcf.gz"
        output_vcf = f"{tmp_dir}/output.vcf.gz"

        # Construct param dict
        param = {"calculation": {"calculations": {"dp_stats": None}}}

        # Create object
        variants = Variants(
            conn=None, input=input_vcf, output=output_vcf, param=param, load=True
        )

        # Remove if output file exists
        remove_if_exists([output_vcf])

        # Calculation
        variants.calculation()

        result = variants.get_query_to_df(
            """ SELECT INFO FROM variants WHERE INFO LIKE '%DP_stats%' """
        )
        assert len(result) == 7

        result = variants.get_query_to_df(
            """ SELECT * FROM variants WHERE "#CHROM" = 'chr1' AND POS = 28736 AND INFO LIKE '%DP_stats_nb=4%' """
        )
        assert len(result) == 1

        result = variants.get_query_to_df(
            """ SELECT * FROM variants WHERE "#CHROM" = 'chr1' AND POS = 28736 AND INFO LIKE '%DP_stats_min=576.0%' """
        )
        assert len(result) == 1

        result = variants.get_query_to_df(
            """ SELECT * FROM variants WHERE "#CHROM" = 'chr1' AND POS = 28736 AND INFO LIKE '%DP_stats_max=17664.0%' """
        )
        assert len(result) == 1

        result = variants.get_query_to_df(
            """ SELECT * FROM variants WHERE "#CHROM" = 'chr1' AND POS = 28736 AND INFO LIKE '%DP_stats_mean=9158.0%' """
        )
        assert len(result) == 1

        # Check if VCF is in correct format with pyVCF
        remove_if_exists([output_vcf])
        variants.export_output()
        try:
            vcf.Reader(filename=output_vcf)
        except:
            assert False


def test_calculation_variant_id():
    """
    This is a test function for the calculation of depth statistics in a VCF file using the Variants
    class in Python.
    """

    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # Init files
        input_vcf = tests_data_folder + "/example.vcf.gz"
        output_vcf = f"{tmp_dir}/output.vcf.gz"

        # Construct param dict
        param = {"calculation": {"calculations": {"variant_id": None}}}

        # Create object
        variants = Variants(
            conn=None, input=input_vcf, output=output_vcf, param=param, load=True
        )

        # Remove if output file exists
        remove_if_exists([output_vcf])

        # Calculation
        variants.calculation()

        # Check if all variant have variant_id
        result = variants.get_query_to_df(
            """ SELECT INFO FROM variants WHERE INFO LIKE '%variant_id%' """
        )
        assert len(result) == 7

        # Explode info
        variants.explode_infos(prefix="INFO/")

        # Check distinct variant_id
        result = variants.get_query_to_df(
            """ SELECT distinct "INFO/variant_id" FROM variants """
        )
        assert len(result) == 7

        # Check if VCF is in correct format with pyVCF
        remove_if_exists([output_vcf])
        variants.export_output()
        try:
            vcf.Reader(filename=output_vcf)
        except:
            assert False
