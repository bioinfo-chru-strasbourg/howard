# -*- coding: utf-8 -*-
"""
Tests

Usage:
pytest tests/

Coverage:
coverage run -m pytest tests/test_variants_calculations.py -x -v --log-cli-level=INFO --capture=tee-sys
coverage report --include=howard/* -m
"""

import logging as log
import os
from tempfile import TemporaryDirectory
import pytest  # type: ignore
import vcf  # type: ignore

from howard.objects.variants import Variants
from howard.functions.commons import remove_if_exists

from test_needed import tests_folder, tests_data_folder, tests_config, database_files


def test_calculation_transcript_annotations():
    """
    Test calculation `transcript_annotations` which performs transcripts annotations and generate a transcripts table/view, using dict parameters directly in calculation paramters.
    """

    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # Init files
        input_vcf = f"{tests_data_folder}/example.ann.transcripts.vcf.gz"
        output_vcf = f"{tmp_dir}/output.vcf.gz"

        # Construct param dict
        param = {
            "calculation": {
                "calculations": {
                    "transcripts_annotations": {
                        "table": "transcripts",
                        "column_id": "transcript",
                        "transcripts_info_json": "transcripts_json",
                        "transcripts_info_field": "transcripts_json",
                        "struct": {
                            "from_column_format": [
                                {
                                    "transcripts_column": "ANN",
                                    "transcripts_infos_column": "Feature_ID",
                                    "column_rename": None,
                                    "column_clean": True,
                                    "column_case": None,
                                },
                                {
                                    "transcripts_column": "ANN",
                                    "transcripts_infos_column": "Feature_ID",
                                    "column_rename": None,
                                    "column_clean": True,
                                    "column_case": None,
                                },
                            ],
                            "from_columns_map": [
                                {
                                    "transcripts_column": "Ensembl_transcriptid",
                                    "transcripts_infos_columns": [
                                        "genename",
                                        "Ensembl_geneid",
                                        "LIST_S2_score",
                                        "LIST_S2_pred",
                                    ],
                                    "column_rename": None,
                                    "column_clean": False,
                                    "column_case": None,
                                },
                                {
                                    "transcripts_column": "Ensembl_transcriptid",
                                    "transcripts_infos_columns": [
                                        "genename",
                                        "VARITY_R_score",
                                        "Aloft_pred",
                                    ],
                                    "column_rename": None,
                                    "column_clean": False,
                                    "column_case": None,
                                },
                            ],
                            "from_variants": {
                                "fields": ["CLNSIG", "DP"],
                                "prefix": "",
                            },
                        },
                    }
                }
            }
        }

        # Create object
        variants = Variants(
            input=input_vcf,
            output=output_vcf,
            config=tests_config,
            param=param,
            load=True,
        )

        # Calculation
        variants.calculation()

        # Check if transcript table exists and is not empty
        result = variants.get_query_to_df("""SELECT INFO FROM transcripts """)
        assert len(result) > 0

        # Check if VCF is in correct format with pyVCF
        remove_if_exists([output_vcf])
        variants.export_output()
        try:
            vcf.Reader(filename=output_vcf)
        except:
            assert False


def test_calculation_transcript_annotations_and_prioritization():
    """
    Test calculation `transcript_annotations` which performs transcripts annotations and generate a transcripts table/view, using dict parameters directly in calculation paramters. Then check if the generated transcripts table can be used for the `transcripts_prioritization` calculation and if the output VCF is in the correct format with pyVCF.
    """

    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # Init files
        input_vcf = f"{tests_data_folder}/example.ann.transcripts.vcf.gz"
        output_vcf = f"{tmp_dir}/output.vcf.gz"

        # Construct param dict
        param = {
            "calculation": {
                "calculations": {
                    "transcripts_annotations": {
                        "table": "transcripts",
                        "column_id": "transcript",
                        "transcripts_info_json": "transcripts_json",
                        "transcripts_info_field": "transcripts_json",
                        "struct": {
                            "from_column_format": [
                                {
                                    "transcripts_column": "ANN",
                                    "transcripts_infos_column": "Feature_ID",
                                    "column_rename": None,
                                    "column_clean": True,
                                    "column_case": None,
                                },
                                {
                                    "transcripts_column": "ANN",
                                    "transcripts_infos_column": "Feature_ID",
                                    "column_rename": None,
                                    "column_clean": True,
                                    "column_case": None,
                                },
                            ],
                            "from_columns_map": [
                                {
                                    "transcripts_column": "Ensembl_transcriptid",
                                    "transcripts_infos_columns": [
                                        "genename",
                                        "Ensembl_geneid",
                                        "LIST_S2_score",
                                        "LIST_S2_pred",
                                    ],
                                    "column_rename": None,
                                    "column_clean": False,
                                    "column_case": None,
                                },
                                {
                                    "transcripts_column": "Ensembl_transcriptid",
                                    "transcripts_infos_columns": [
                                        "genename",
                                        "VARITY_R_score",
                                        "Aloft_pred",
                                    ],
                                    "column_rename": None,
                                    "column_clean": False,
                                    "column_case": None,
                                },
                            ],
                            "from_variants": {
                                "fields": ["CLNSIG", "DP"],
                                "prefix": "",
                            },
                        },
                    },
                    "transcripts_prioritization": {
                        "prioritization": {
                            "profiles": ["transcripts"],
                            "prioritization_config": f"{tests_data_folder}/prioritization_transcripts_profiles.json",
                            "pzprefix": "PZT",
                            "prioritization_score_mode": "HOWARD",
                        },
                    },
                }
            }
        }

        # Create object
        variants = Variants(
            input=input_vcf,
            output=output_vcf,
            config=tests_config,
            param=param,
            load=True,
        )

        # Calculation
        variants.calculation()

        # Check if transcript table exists and is not empty
        result = variants.get_query_to_df("""SELECT INFO FROM transcripts """)
        assert len(result) > 0

        # Check number of variant with transcript prioritization
        result = variants.get_query_to_df(
            """SELECT INFO FROM variants WHERE INFO LIKE '%PZTTranscript=%' """
        )
        assert len(result) == 7

        # Check if VCF is in correct format with pyVCF
        remove_if_exists([output_vcf])
        variants.export_output()
        try:
            vcf.Reader(filename=output_vcf)
        except:
            assert False


def test_calculation_transcript_annotations_and_export():
    """
    Test calculation `transcript_annotations` which performs transcripts annotations and generate a transcripts table/view, using dict parameters directly in calculation paramters. Then check if the generated transcripts table can be used for the `transcripts_export` calculation and if the output VCF is in the correct format with pyVCF.
    """

    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # Init files
        input_vcf = f"{tests_data_folder}/example.ann.transcripts.vcf.gz"
        output_vcf = f"{tmp_dir}/output.vcf.gz"

        # Construct param dict
        param = {
            "calculation": {
                "calculations": {
                    "transcripts_annotations": {
                        "table": "transcripts",
                        "column_id": "transcript",
                        "transcripts_info_json": "transcripts_json",
                        "transcripts_info_field": "transcripts_json",
                        "struct": {
                            "from_column_format": [
                                {
                                    "transcripts_column": "ANN",
                                    "transcripts_infos_column": "Feature_ID",
                                    "column_rename": None,
                                    "column_clean": True,
                                    "column_case": None,
                                },
                                {
                                    "transcripts_column": "ANN",
                                    "transcripts_infos_column": "Feature_ID",
                                    "column_rename": None,
                                    "column_clean": True,
                                    "column_case": None,
                                },
                            ],
                            "from_columns_map": [
                                {
                                    "transcripts_column": "Ensembl_transcriptid",
                                    "transcripts_infos_columns": [
                                        "genename",
                                        "Ensembl_geneid",
                                        "LIST_S2_score",
                                        "LIST_S2_pred",
                                    ],
                                    "column_rename": None,
                                    "column_clean": False,
                                    "column_case": None,
                                },
                                {
                                    "transcripts_column": "Ensembl_transcriptid",
                                    "transcripts_infos_columns": [
                                        "genename",
                                        "VARITY_R_score",
                                        "Aloft_pred",
                                    ],
                                    "column_rename": None,
                                    "column_clean": False,
                                    "column_case": None,
                                },
                            ],
                            "from_variants": {
                                "fields": ["CLNSIG", "DP"],
                                "prefix": "",
                            },
                        },
                    },
                    "transcripts_export": {
                        "export": {
                            "output": f"{tmp_dir}/output.transcripts.tsv",
                            "export_header": False,
                            "add_info": True,
                        },
                        "explode": {
                            "explode_infos": True,
                            "explode_infos_fields": "PZTScore,variant_type_rank",
                        },
                    },
                }
            }
        }

        # Create object
        variants = Variants(
            input=input_vcf,
            output=output_vcf,
            config=tests_config,
            param=param,
            load=True,
        )

        # Calculation
        variants.calculation()

        # Check if transcript table exists and is not empty
        result = variants.get_query_to_df("""SELECT INFO FROM transcripts """)
        assert len(result) > 0

        # Check if output transcripts exists and as content (at least one line with a transcript, so 2 lines with header)
        assert os.path.exists(f"{tmp_dir}/output.transcripts.tsv")
        if os.path.exists(f"{tmp_dir}/output.transcripts.tsv"):
            with open(f"{tmp_dir}/output.transcripts.tsv", "r") as f:
                lines = f.readlines()
                content_lines = [line for line in lines if line.strip()]
                assert len(content_lines) > 1

        # Check if VCF is in correct format with pyVCF
        remove_if_exists([output_vcf])
        variants.export_output()
        try:
            vcf.Reader(filename=output_vcf)
        except:
            assert False


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
            "VARIANT_CHROM_POS_REF_ALT": None,
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
            f"""SELECT INFO FROM variants WHERE INFO LIKE '%VARIANT_CHROM_POS_REF_ALT=%' """
        )
        assert len(result) == 7, f"Expected 7 variants with 'VARIANT_CHROM_POS_REF_ALT' but got {len(result)}"

        # Check number of variant_chr_pos_alt_ref_dict with new table
        result = variants.get_query_to_df(
            f"""SELECT INFO FROM {new_table} WHERE INFO LIKE '%variant_chr_pos_alt_ref_dict=%' """
        )
        assert len(result) == 7, f"Expected 7 variants with 'variant_chr_pos_alt_ref_dict' but got {len(result)}"

        # Check number of variant_chr_pos_alt_ref_file
        result = variants.get_query_to_df(
            f"""SELECT INFO FROM variants WHERE INFO LIKE '%variant_chr_pos_alt_ref_file=%' """
        )
        assert len(result) == 7, f"Expected 7 variants with 'variant_chr_pos_alt_ref_file' but got {len(result)}"

        # Check number of middle (7)
        result = variants.get_query_to_df(
            f"""SELECT INFO FROM variants WHERE INFO LIKE '%VARIANT_CHROM_POS_REF_ALT=chr1_28736_A_C%' """
        )
        assert len(result) == 1, f"Expected 1 variant with 'VARIANT_CHROM_POS_REF_ALT=chr1_28736_A_C' but got {len(result)}"

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
            "VARIANT_CHROM_POS_REF_ALT": None,
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
            """SELECT INFO FROM variants WHERE INFO LIKE '%VARIANT_CHROM_POS_REF_ALT=%' """
        )
        assert len(result) == 7, f"Expected 7 variants with 'VARIANT_CHROM_POS_REF_ALT' but got {len(result)}"

        # Check number of variant_chr_pos_alt_ref_dict
        result = variants.get_query_to_df(
            """SELECT INFO FROM variants WHERE INFO LIKE '%variant_chr_pos_alt_ref_dict=%' """
        )
        assert len(result) == 7, f"Expected 7 variants with 'variant_chr_pos_alt_ref_dict' but got {len(result)}"

        # Check number of variant_chr_pos_alt_ref_file
        result = variants.get_query_to_df(
            """SELECT INFO FROM variants WHERE INFO LIKE '%variant_chr_pos_alt_ref_file=%' """
        )
        assert len(result) == 7, f"Expected 7 variants with 'variant_chr_pos_alt_ref_file' but got {len(result)}"

        # Check number of middle (7)
        result = variants.get_query_to_df(
            """SELECT INFO FROM variants WHERE INFO LIKE '%VARIANT_CHROM_POS_REF_ALT=chr1_28736_A_C%' """
        )
        assert len(result) == 1, f"Expected 1 variant with 'VARIANT_CHROM_POS_REF_ALT=chr1_28736_A_C' but got {len(result)}"

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
        param = {"calculation": {"calculation_config": operations_config_file}}
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
            "VARIANT_CHROM_POS_REF_ALT": None,
            "variant_chr_pos_alt_ref_file": None,
            "variant_chr_pos_alt_ref_dict": None,
        }

        # Calculation
        variants.calculation(
            operations=operations,
            operations_config_dict=operations_config_dict,
        )

        # Check number of VARIANT_CHROM_POS_REF_ALT
        result = variants.get_query_to_df(
            """SELECT INFO FROM variants WHERE INFO LIKE '%VARIANT_CHROM_POS_REF_ALT=%' """
        )
        assert len(result) == 7, f"Expected 7 variants with 'VARIANT_CHROM_POS_REF_ALT' but got {len(result)}"

        # Check number of variant_chr_pos_alt_ref_dict
        result = variants.get_query_to_df(
            """SELECT INFO FROM variants WHERE INFO LIKE '%variant_chr_pos_alt_ref_dict=%' """
        )
        assert len(result) == 7, f"Expected 7 variants with 'variant_chr_pos_alt_ref_dict' but got {len(result)}"

        # Check number of variant_chr_pos_alt_ref_file
        result = variants.get_query_to_df(
            """SELECT INFO FROM variants WHERE INFO LIKE '%variant_chr_pos_alt_ref_file=%' """
        )
        assert len(result) == 7, f"Expected 7 variants with 'variant_chr_pos_alt_ref_file' but got {len(result)}"

        # Check number of middle (7)
        result = variants.get_query_to_df(
            """SELECT INFO FROM variants WHERE INFO LIKE '%VARIANT_CHROM_POS_REF_ALT=chr1_28736_A_C%' """
        )
        assert len(result) == 1, f"Expected 1 variant with 'VARIANT_CHROM_POS_REF_ALT=chr1_28736_A_C' but got {len(result)}"

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
        input_vcf = tests_data_folder + "/example.annovar.vcf"
        output_vcf = f"{tmp_dir}/output.vcf"
        input_param = {
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


def test_calculation_nomen_one_hgvs_field():
    """
    This function tests the calculation and annotation of genetic variants using input parameters and
    checks if the output VCF file is in the correct format.
    """

    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # Init files
        input_vcf = tests_data_folder + "/example.annotated.hgvs.annovar.snpeff.vcf"
        output_vcf = f"{tmp_dir}/output.vcf.gz"
        input_param = {
            "calculation": {
                "calculations": {
                    "NOMEN": {
                        "options": {
                            "hgvs_fields": ["AAChange_refGene"],
                        }
                    }
                }
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


def test_calculation_nomen_two_hgvs_fields():
    """
    This function tests the calculation and annotation of genetic variants using input parameters and
    checks if the output VCF file is in the correct format.
    """

    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # Init files
        input_vcf = tests_data_folder + "/example.annotated.hgvs.annovar.snpeff.vcf"
        output_vcf = f"{tmp_dir}/output.vcf.gz"
        input_param = {
            "calculation": {
                "calculations": {
                    "NOMEN": {
                        "options": {
                            "hgvs_fields": ["AAChange_refGene", "snpeff_hgvs"],
                        }
                    }
                }
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
        assert len(result) == 7

        # Check if VCF is in correct format with pyVCF
        remove_if_exists([output_vcf])
        variants.export_output()
        try:
            vcf.Reader(filename=output_vcf)
        except:
            assert False


def test_calculation_nomen_two_hgvs_fields_uniquify():
    """
    This function tests the calculation and annotation of genetic variants using input parameters and
    checks if the output VCF file is in the correct format.
    """

    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # Init files
        input_vcf = tests_data_folder + "/example.annotated.hgvs.annovar.snpeff.vcf"
        output_vcf = f"{tmp_dir}/output.vcf.gz"
        input_param = {
            "calculation": {
                "calculations": {
                    "NOMEN": {
                        "options": {
                            "hgvs_fields": ["AAChange_refGene", "snpeff_hgvs"],
                            "uniquify_hgvs": False,
                        }
                    }
                }
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
        assert len(result) == 7

        # Check if VCF is in correct format with pyVCF
        remove_if_exists([output_vcf])
        variants.export_output()
        try:
            vcf.Reader(filename=output_vcf)
        except:
            assert False


# def test_calculation_merged_hgvs():
#     """
#     This function tests the calculation of merging multiple HGVS annotations.
#     """

#     with TemporaryDirectory(dir=tests_folder) as tmp_dir:

#         # Init files
#         # input_vcf = tests_data_folder + "/example.annotated.annovar.snpeff.vcf"
#         input_vcf = tests_data_folder + "/example.annotated.hgvs.annovar.snpeff.vcf"
#         output_vcf = f"{tmp_dir}/output.vcf.gz"
#         input_param = {
#             "calculation": {"calculations": {"MERGED_HGVS": None}},
#         }

#         # Create object
#         variants = Variants(
#             input=input_vcf,
#             output=output_vcf,
#             config=tests_config,
#             param=input_param,
#             load=True,
#         )

#         # Annotation
#         variants.annotation()

#         # Calculation
#         variants.calculation()

#         # Check number of NOMEN
#         annovar_hgvs = "AAChange_refGene=EGFR:NM_001346941:exon14:c.1560G>A:p.Q520Q,EGFR:NM_001346897:exon19:c.2226G>A:p.Q742Q,EGFR:NM_001346899:exon19:c.2226G>A:p.Q742Q,EGFR:NM_001346898:exon20:c.2361G>A:p.Q787Q,EGFR:NM_001346900:exon20:c.2202G>A:p.Q734Q,EGFR:NM_005228:exon20:c.2361G>A:p.Q787Q"
#         snpeff_hgvs = "snpeff_hgvs=EGFR:NM_005228.5:exon20:c.2361G>A:p.Gln787Gln,EGFR:NM_001346897.2:exon19:c.2226G>A:p.Gln742Gln,EGFR:NM_001346898.2:exon20:c.2361G>A:p.Gln787Gln,EGFR:NM_001346941.2:exon14:c.1560G>A:p.Gln520Gln,EGFR:NM_001346899.1:exon19:c.2226G>A:p.Gln742Gln,EGFR:NM_001346900.2:exon20:c.2202G>A:p.Gln734Gln,EGFR-AS1:NR_047551.1:exon2:n.1201C>T"
#         # merged_hgvs = "merged_hgvs=EGFR:NM_005228:exon20:c.2361G>A:p.Q787Q,EGFR:NM_001346898:exon20:c.2361G>A:p.Q787Q,EGFR:NM_001346898.2:exon20:c.2361G>A:p.Gln787Gln,EGFR:NM_001346899:exon19:c.2226G>A:p.Q742Q,EGFR:NM_001346900.2:exon20:c.2202G>A:p.Gln734Gln,EGFR:NM_001346941:exon14:c.1560G>A:p.Q520Q,EGFR:NM_001346900:exon20:c.2202G>A:p.Q734Q,EGFR:NM_005228.5:exon20:c.2361G>A:p.Gln787Gln,EGFR-AS1:NR_047551.1:exon2:n.1201C>T,EGFR:NM_001346899.1:exon19:c.2226G>A:p.Gln742Gln,EGFR:NM_001346941.2:exon14:c.1560G>A:p.Gln520Gln,EGFR:NM_001346897:exon19:c.2226G>A:p.Q742Q,EGFR:NM_001346897.2:exon19:c.2226G>A:p.Gln742Gln"
#         merged_hgvs = "merged_hgvs=EGFR-AS1:NR_047551.1:exon2:n.1201C>T,EGFR:NM_001346897.2:exon19:c.2226G>A:p.Gln742Gln,EGFR:NM_001346897:exon19:c.2226G>A:p.Q742Q,EGFR:NM_001346898.2:exon20:c.2361G>A:p.Gln787Gln,EGFR:NM_001346898:exon20:c.2361G>A:p.Q787Q,EGFR:NM_001346899.1:exon19:c.2226G>A:p.Gln742Gln,EGFR:NM_001346899:exon19:c.2226G>A:p.Q742Q,EGFR:NM_001346900.2:exon20:c.2202G>A:p.Gln734Gln,EGFR:NM_001346900:exon20:c.2202G>A:p.Q734Q,EGFR:NM_001346941.2:exon14:c.1560G>A:p.Gln520Gln,EGFR:NM_001346941:exon14:c.1560G>A:p.Q520Q,EGFR:NM_005228.5:exon20:c.2361G>A:p.Gln787Gln,EGFR:NM_005228:exon20:c.2361G>A:p.Q787Q"
#         result = variants.get_query_to_df(
#             f"""SELECT INFO FROM variants WHERE INFO LIKE '%{annovar_hgvs}%' AND INFO LIKE '%{snpeff_hgvs}%' AND INFO LIKE '%{merged_hgvs}%' """
#         )
#         assert len(result) == 1

#         # Check if VCF is in correct format with pyVCF
#         remove_if_exists([output_vcf])
#         variants.export_output()
#         try:
#             vcf.Reader(filename=output_vcf)
#         except:
#             assert False


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
        assert len(result) == 4  # 6 if duplication !!!

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
    "input_file, expected",
    [
        (
            "example.vcf.gz",
            0,
        ),
        (
            "example.ann.vcf.gz",
            7
        )
    ],
)
def test_calculation_snpeff_hgvs(input_file, expected):
    """
    This function is a test for calculating snpeff_hgvs in a VCF file using the Variants class.
    """

    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # Init files
        input_vcf = tests_data_folder + f"/{input_file}"
        output_vcf = f"{tmp_dir}/output.vcf.gz"

        # Construct param dict
        #param = {"calculation": {"calculations": {"snpeff_hgvs": None}}}
        param = {
            "calculation": {
                "calculations": {
                    "annotation_extract": {
                        "annotation_field": "ANN",
                        "annotation_hgvs": "snpeff_hgvs"
                    }
                }
            }
        }

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
        assert len(result) == expected

        # Check if VCF is in correct format with pyVCF
        remove_if_exists([output_vcf])
        variants.export_output()
        try:
            vcf.Reader(filename=output_vcf)
        except:
            assert False


@pytest.mark.parametrize(
    "input_file, expected, explode, json, uniquify",
    [
        (
            input_file[0],
            input_file[1],
            explode,
            json,
            uniquify,
        )
        for input_file in [("example.ann.vcf.gz", 7), ("example.chrM.bad_ann.vcf", 10), ("example.chrM.ann.vcf", 10)]
        for explode in [True, False]
        for json in [True, False]
        for uniquify in [True, False]
    ],
)
def test_calculation_snpeff_explode(input_file, expected, explode, json, uniquify):
    """
    This function is a test for calculating snpeff_hgvs in a VCF file using the Variants class.
    """

    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # Init files
        input_vcf = tests_data_folder + f"/{input_file}"
        output_vcf = f"{tmp_dir}/output.vcf.gz"

        # Construct param dict
        param = {
            "calculation": {
                "calculations": {
                    "annotation_extract": {
                        "annotation_field": "ANN",
                        "annotation_explode": "snpeff_" if explode else None,
                        "annotation_json": "snpeff_json" if json else None,
                        "uniquify": uniquify
                    }
                }
            }
        }

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
            """ SELECT * FROM variants WHERE INFO LIKE '%snpeff_%' """
        )
        assert len(result) == expected, f"For 'snpeff_', Expected {expected} but got {len(result)} for input_file={input_file}, explode={explode}, json={json}, uniquify={uniquify}"

        if explode:

            # query annotated variant
            result = variants.get_query_to_df(
                """ SELECT * FROM variants WHERE INFO LIKE '%snpeff_Allele=%' """
            )
            assert len(result) == expected, f"For 'snpeff_Allele', Expected {expected} but got {len(result)} for input_file={input_file}, explode={explode}, json={json}, uniquify={uniquify}"

            # query annotated variant
            result = variants.get_query_to_df(
                """ SELECT * FROM variants WHERE INFO LIKE '%snpeff_Annotation%' """
            )
            assert len(result) == expected, f"For 'snpeff_Annotation', Expected {expected} but got {len(result)} for input_file={input_file}, explode={explode}, json={json}, uniquify={uniquify}"

        if json:

            # query annotated variant
            result = variants.get_query_to_df(
                """ SELECT * FROM variants WHERE INFO LIKE '%snpeff_json=%' """
            )
            assert len(result) == expected, f"For 'snpeff_json', Expected {expected} but got {len(result)} for input_file={input_file}, explode={explode}, json={json}, uniquify={uniquify}"

        # Uniquify check: For the variant chr1:28736, there are 5 annotations with snpeff_Allele=C,C,C,C,C. If uniquify is True, we should only have one entry for this variant.
        check_variant_exists = variants.get_query_to_df(
            """ SELECT * FROM variants WHERE "#CHROM" = 'chr1' AND POS = 28736"""
        )
        if len(check_variant_exists) == 1:
            if uniquify:
                result = variants.get_query_to_df(
                    """ SELECT * FROM variants WHERE "#CHROM" = 'chr1' AND POS = 28736 AND INFO LIKE '%snpeff_Allele=C;%' """
                )
                assert len(result) == 1, f"For 'snpeff_Allele=C;', Expected 1 but got {len(result)} for input_file={input_file}, explode={explode}, json={json}, uniquify={uniquify}"
            else:
                result = variants.get_query_to_df(
                    """ SELECT * FROM variants WHERE "#CHROM" = 'chr1' AND POS = 28736 AND INFO LIKE '%snpeff_Allele=C,C,C,C,C;%' """
                )
                assert len(result) == 1, f"For 'snpeff_Allele=C,C,C,C,C;', Expected 1 but got {len(result)} for input_file={input_file}, explode={explode}, json={json}, uniquify={uniquify}"

        # Check if VCF is in correct format with pyVCF
        remove_if_exists([output_vcf])
        variants.export_output()
        try:
            vcf.Reader(filename=output_vcf)
        except:
            assert False


def test_calculation_nomen_snpeff_hgvs_transcripts():
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
        # result = variants.get_query_to_df(""" SELECT INFO FROM variants  """)
        # log.debug(f"result=\n{result.to_string()}")

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


def test_calculation_nomen_snpeff_hgvs_notranscripts():
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


def test_calculation_find_samples():
    """
    This is a test function for the "FIND_SAMPLES" calculation in the Variants class, which checks if
    the calculation is performed correctly and the output VCF file is in the correct format.
    """

    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # Init files
        input_vcf = tests_data_folder + "/example.vcf.gz"
        output_vcf = f"{tmp_dir}/output.vcf.gz"

        # Construct param dict
        param = {"calculation": {"calculations": {"FIND_SAMPLES": None}}}

        # Create object
        variants = Variants(
            conn=None, input=input_vcf, output=output_vcf, param=param, load=True
        )

        # Calculation
        variants.calculation()

        result = variants.get_query_to_df(
            """ SELECT INFO FROM variants WHERE INFO LIKE '%count_samples%' AND INFO LIKE '%list_samples%' """
        )
        assert len(result) == 7

        result = variants.get_query_to_df(
            """ SELECT * FROM variants WHERE INFO LIKE '%count_samples=4/4%'  """
        )
        assert len(result) == 1

        result = variants.get_query_to_df(
            """ SELECT * FROM variants WHERE INFO LIKE '%count_samples=3/4%' """
        )
        assert len(result) == 6

        # Check if VCF is in correct format with pyVCF
        remove_if_exists([output_vcf])
        variants.export_output()
        try:
            vcf.Reader(filename=output_vcf)
        except:
            assert False


def test_calculation_find_samples_options_in_calculation():
    """
    This is a test function for the "FIND_SAMPLES" calculation in the Variants class, which checks if
    the calculation is performed correctly and the output VCF file is in the correct format.
    """

    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # Init files
        input_vcf = tests_data_folder + "/example.vcf.gz"
        output_vcf = f"{tmp_dir}/output.vcf.gz"

        # Construct param dict
        param = {"calculation": {"calculations": {"FIND_SAMPLES": {"tags": {"count_samples_for_variant": "count", "list_samples_for_variant": "list"}}}}}

        # Create object
        variants = Variants(
            conn=None, input=input_vcf, output=output_vcf, param=param, load=True
        )

        # Calculation
        variants.calculation()

        result = variants.get_query_to_df(
            """ SELECT INFO FROM variants WHERE INFO LIKE '%count_samples_for_variant%' AND INFO LIKE '%list_samples_for_variant%' """
        )
        assert len(result) == 7

        result = variants.get_query_to_df(
            """ SELECT * FROM variants WHERE INFO LIKE '%count_samples_for_variant=4/4%' """
        )
        assert len(result) == 1

        result = variants.get_query_to_df(
            """ SELECT * FROM variants WHERE INFO LIKE '%count_samples_for_variant=3/4%' """
        )
        assert len(result) == 6

        # Check if VCF is in correct format with pyVCF
        remove_if_exists([output_vcf])
        variants.export_output()
        try:
            vcf.Reader(filename=output_vcf)
        except:
            assert False


def test_calculation_find_samples_options_in_calculation_only_count_samples():
    """
    This is a test function for the "FIND_SAMPLES" calculation in the Variants class, which checks if
    the calculation is performed correctly and the output VCF file is in the correct format.
    """

    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # Init files
        input_vcf = tests_data_folder + "/example.vcf.gz"
        output_vcf = f"{tmp_dir}/output.vcf.gz"

        # Construct param dict
        param = {"calculation": {"calculations": {"FIND_SAMPLES": {"tags": {"count_samples_for_variant": "count"}}}}}

        # Create object
        variants = Variants(
            conn=None, input=input_vcf, output=output_vcf, param=param, load=True
        )

        # Calculation
        variants.calculation()

        result = variants.get_query_to_df(
            """ SELECT INFO FROM variants WHERE INFO LIKE '%count_samples_for_variant%' AND INFO NOT LIKE '%list_samples_for_variant%' """
        )
        assert len(result) == 7

        result = variants.get_query_to_df(
            """ SELECT * FROM variants WHERE INFO LIKE '%count_samples_for_variant=4/4%' """
        )
        assert len(result) == 1

        result = variants.get_query_to_df(
            """ SELECT * FROM variants WHERE INFO LIKE '%count_samples_for_variant=3/4%' """
        )
        assert len(result) == 6

        # Check if VCF is in correct format with pyVCF
        remove_if_exists([output_vcf])
        variants.export_output()
        try:
            vcf.Reader(filename=output_vcf)
        except:
            assert False


def test_calculation_find_samples_options_in_calculation_only_list_samples():
    """
    This is a test function for the "FIND_SAMPLES" calculation in the Variants class, which checks if
    the calculation is performed correctly and the output VCF file is in the correct format.
    """

    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # Init files
        input_vcf = tests_data_folder + "/example.vcf.gz"
        output_vcf = f"{tmp_dir}/output.vcf.gz"

        # Construct param dict
        param = {"calculation": {"calculations": {"FIND_SAMPLES": {"tags": {"list_samples_for_variant": "list"}}}}}

        # Create object
        variants = Variants(
            conn=None, input=input_vcf, output=output_vcf, param=param, load=True
        )

        # Calculation
        variants.calculation()

        result = variants.get_query_to_df(
            """ SELECT INFO FROM variants WHERE INFO NOT LIKE '%count_samples%' AND INFO LIKE '%list_samples_for_variant%' """
        )
        assert len(result) == 7

        # Check if VCF is in correct format with pyVCF
        remove_if_exists([output_vcf])
        variants.export_output()
        try:
            vcf.Reader(filename=output_vcf)
        except:
            assert False


def test_calculation_find_samples_options_in_calculation_two_list_samples():
    """
    This is a test function for the "FIND_SAMPLES" calculation in the Variants class, which checks if
    the calculation is performed correctly and the output VCF file is in the correct format.
    """

    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # Init files
        input_vcf = tests_data_folder + "/example.vcf.gz"
        output_vcf = f"{tmp_dir}/output.vcf.gz"

        # Construct param dict
        param = {"calculation": {"calculations": {"FIND_SAMPLES": {"tags": {"list_samples_for_variant": "list", "list_samples_for_variant_bis": "list"}}}}}

        # Create object
        variants = Variants(
            conn=None, input=input_vcf, output=output_vcf, param=param, load=True
        )

        # Calculation
        variants.calculation()

        result = variants.get_query_to_df(
            """ SELECT INFO FROM variants WHERE INFO NOT LIKE '%count_samples%' AND INFO LIKE '%list_samples_for_variant%' AND INFO LIKE '%list_samples_for_variant_bis%' """
        )
        assert len(result) == 7

        # Check if VCF is in correct format with pyVCF
        remove_if_exists([output_vcf])
        variants.export_output()
        try:
            vcf.Reader(filename=output_vcf)
        except:
            assert False



def test_calculation_find_samples_count_samples():
    """
    This is a test function for the "FIND_SAMPLES" calculation in the Variants class, which checks if
    the calculation is performed correctly and the output VCF file is in the correct format.
    """

    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # Init files
        input_vcf = tests_data_folder + "/example.vcf.gz"
        output_vcf = f"{tmp_dir}/output.vcf.gz"

        # Construct param dict
        #param = {"calculation": {"calculations": {"COUNT_SAMPLES": None}}}
        param = {
            "calculation": {
                "calculations": {
                    "FIND_SAMPLES": {
                        "tags": {"count_samples": "count"}
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

        result = variants.get_query_to_df(
            """ SELECT INFO FROM variants WHERE INFO LIKE '%count_samples%' AND INFO NOT LIKE '%list_samples%' """
        )
        assert len(result) == 7

        result = variants.get_query_to_df(
            """ SELECT * FROM variants WHERE INFO LIKE '%count_samples=4/4%'  """
        )
        assert len(result) == 1

        result = variants.get_query_to_df(
            """ SELECT * FROM variants WHERE INFO LIKE '%count_samples=3/4%' """
        )
        assert len(result) == 6

        # Check if VCF is in correct format with pyVCF
        remove_if_exists([output_vcf])
        variants.export_output()
        try:
            vcf.Reader(filename=output_vcf)
        except:
            assert False

@pytest.mark.parametrize(
    "tag",
    [
            "count_samples_for_variant",
            "count_samples_for_variant_other_tag"
    ],
)
def test_calculation_find_samples_count_samples_options_in_calculation(tag):
    """
    This is a test function for the "FIND_SAMPLES" calculation in the Variants class, which checks if
    the calculation is performed correctly and the output VCF file is in the correct format.
    """

    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # Init files
        input_vcf = tests_data_folder + "/example.vcf.gz"
        output_vcf = f"{tmp_dir}/output.vcf.gz"

        # Construct param dict
        param = {"calculation": {"calculations": {"FIND_SAMPLES": {"tags": {tag: "count"}}}}}

        # Create object
        variants = Variants(
            conn=None, input=input_vcf, output=output_vcf, param=param, load=True
        )

        # Calculation
        variants.calculation()

        result = variants.get_query_to_df(
            f""" SELECT INFO FROM variants WHERE INFO LIKE '%{tag}%' AND INFO NOT LIKE '%list_samples%' """
        )
        assert len(result) == 7, f"With '{tag}', Expected 7 but got {len(result)}"

        result = variants.get_query_to_df(
            f""" SELECT * FROM variants WHERE INFO LIKE '%{tag}=4/4%'  """
        )
        assert len(result) == 1, f"With '{tag}', Expected 1 but got {len(result)}"

        result = variants.get_query_to_df(
            f""" SELECT * FROM variants WHERE INFO LIKE '%{tag}=3/4%' """
        )
        assert len(result) == 6, f"With '{tag}', Expected 6 but got {len(result)}"

        # Check if VCF is in correct format with pyVCF
        remove_if_exists([output_vcf])
        variants.export_output()
        try:
            vcf.Reader(filename=output_vcf)
        except:
            assert False


@pytest.mark.parametrize(
    "tag",
    [
            "list_samples_for_variant",
            "list_samples_for_variant_other_tag"
    ],
)
def test_calculation_find_samples_list_samples(tag):
    """
    This is a test function for the "FIND_SAMPLES" calculation in the Variants class, which checks if
    the calculation is performed correctly and the output VCF file is in the correct format.
    """

    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # Init files
        input_vcf = tests_data_folder + "/example.vcf.gz"
        output_vcf = f"{tmp_dir}/output.vcf.gz"

        # Construct param dict
        param = {"calculation": {"calculations": {"FIND_SAMPLES": {"tags": {tag: "list"}}}}}

        # Create object
        variants = Variants(
            conn=None, input=input_vcf, output=output_vcf, param=param, load=True
        )

        # Calculation
        variants.calculation()

        result = variants.get_query_to_df(
            f""" SELECT INFO FROM variants WHERE INFO NOT LIKE '%count_samples%' AND INFO LIKE '%{tag}=%' """
        )
        assert len(result) == 7, f"With '{tag}', Expected 7 but got {len(result)}"

        # Check if VCF is in correct format with pyVCF
        remove_if_exists([output_vcf])
        variants.export_output()
        try:
            vcf.Reader(filename=output_vcf)
        except:
            assert False


@pytest.mark.parametrize(
    "tag",
    [
            "genotype_concordance",
            "genotype_concordance_other_tag"
    ],
)
def test_calculation_genotype_concordance(tag):
    """
    This is a test function for calculating genotype concordance in a VCF file using the Variants class.
    """

    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # Init files
        input_vcf = tests_data_folder + "/example.vcf.gz"
        output_vcf = f"{tmp_dir}/output.vcf.gz"

        # Construct param dict
        param = {"calculation": {"calculations": {"GENOTYPE_CONCORDANCE": {"tag": tag}}}}

        # Create object
        variants = Variants(
            conn=None, input=input_vcf, output=output_vcf, param=param, load=True
        )

        # Calculation
        variants.calculation()

        result = variants.get_query_to_df(
            f""" SELECT INFO FROM variants WHERE INFO LIKE '%{tag}=%' """
        )
        assert len(result) == 7, f"With '{tag}', Expected 7 but got {len(result)}"

        result = variants.get_query_to_df(
            f""" SELECT * FROM variants WHERE INFO LIKE '%{tag}=FALSE%' """
        )
        assert len(result) == 1, f"With '{tag}', Expected 1 but got {len(result)}"

        result = variants.get_query_to_df(
            f""" SELECT * FROM variants WHERE INFO LIKE '%{tag}=TRUE%' """
        )
        assert len(result) == 6, f"With '{tag}', Expected 6 but got {len(result)}"

        # Check if VCF is in correct format with pyVCF
        remove_if_exists([output_vcf])
        variants.export_output()
        try:
            vcf.Reader(filename=output_vcf)
        except:
            assert False

@pytest.mark.parametrize(
    "input_file, tag",
    [
        *[
            (input_file, tag)
            for input_file in ["example.vcf.gz", "example.name_with_special_char.vcf"]
            for tag in ["barcode", "barcode_other_tag", "BC"]
        ]
    ],
)
def test_calculation_barcode(input_file, tag):
    """
    This is a test function for a Python script that calculates barcode information from a VCF file and
    checks if the output is correct.
    """

    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # Init files
        input_vcf = tests_data_folder + f"/{input_file}"
        output_vcf = f"{tmp_dir}/output.vcf.gz"

        # Construct param dict
        param = {"calculation": {"calculations": {"BARCODE": {"tag": tag}}}}

        # Create object
        variants = Variants(
            conn=None, input=input_vcf, output=output_vcf, param=param, load=True
        )

        # Calculation
        variants.calculation()

        result = variants.get_query_to_df(
            f""" SELECT INFO FROM variants WHERE INFO LIKE '%{tag}=%' """
        )
        assert len(result) == 7

        result = variants.get_query_to_df(
            f""" SELECT * FROM variants WHERE INFO LIKE '%{tag}=1122%' """
        )
        assert len(result) == 1

        result = variants.get_query_to_df(
            f""" SELECT * FROM variants WHERE INFO LIKE '%{tag}=0111%' """
        )
        assert len(result) == 1

        result = variants.get_query_to_df(
            f""" SELECT * FROM variants WHERE INFO LIKE '%{tag}=1011%' """
        )
        assert len(result) == 4

        result = variants.get_query_to_df(
            f""" SELECT * FROM variants WHERE INFO LIKE '%{tag}=1101%' """
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
    "tag, tag_samples, family_pedigree",
    [
        (tag, tag_samples, family_pedigree)
        for tag in [None, "BCF", "BarCodeFamily"]
        for tag_samples in [None, "BCFS", "BarCodeFamilySamples"]
        for family_pedigree in [None, "", "sample1,sample2,sample4", {"father": "sample1", "mother": "sample2", "child": "sample4"}, os.path.join(tests_data_folder, "trio.json")]
    ],
)
def test_calculation_barcode_family(tag, tag_samples, family_pedigree):
    """
    The function `test_calculation_barcode_family` is a test function in Python that calculates
    barcode information from a VCF file and checks if the output is correct.
    """

    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # Init files
        input_vcf = tests_data_folder + "/example.vcf.gz"
        output_vcf = f"{tmp_dir}/output.vcf.gz"

        # Construct param dict
        param = {
            "calculation": {"calculations": {"BARCODE_FAMILY": {"tag": tag, "tag_samples": tag_samples, "family_pedigree": family_pedigree}}}
        }

        # Expected for default valeurs
        if tag is None:
            tag = "BCF"
        if tag_samples is None:
            tag_samples = f"{tag}S"

        # Create object
        variants = Variants(
            conn=None, input=input_vcf, output=output_vcf, param=param, load=True
        )

        # Calculation
        variants.calculation()

        # Create object
        variants = Variants(
            conn=None, input=input_vcf, output=output_vcf, param=param, load=True
        )

        # Calculation
        variants.calculation()

        # Check if BCF and BCFS are in FORMAT
        result = variants.get_query_to_df(
            f""" SELECT FORMAT FROM variants WHERE FORMAT LIKE '%:{tag}:{tag_samples}' """
        )
        assert len(result) == 7, f"With '{tag}' and '{tag_samples}', Expected 7 but got {len(result)}"

        # Check if VCF is in correct format with pyVCF
        remove_if_exists([output_vcf])
        variants.export_output()
        try:
            vcf.Reader(filename=output_vcf)
        except:
            assert False


@pytest.mark.parametrize(
    "tag, trio_pedigree",
    [
        (tag, trio_pedigree)
        for tag in [None, "trio", "TrioCalculation"]
        for trio_pedigree in [None, "", "sample1,sample2,sample4", {"father": "sample1", "mother": "sample2", "child": "sample4"}, os.path.join(tests_data_folder, "trio.json"), {"child": "sample4", "mother": "sample2", "father": "sample1"}]
    ],
)
def test_calculation_trio(tag, trio_pedigree):
    """
    This is a test function for the calculation of trio variants in a VCF file using specific
    parameters.
    """

    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # Init files
        input_vcf = tests_data_folder + "/example.vcf.gz"
        output_vcf = f"{tmp_dir}/output.vcf.gz"

        # Param NO
        param = {
                    "calculation": {"calculations": {"TRIO": {"tag": tag, "trio_pedigree": trio_pedigree}}}
                }

        # Expected for default valeurs
        if tag is None:
            tag = "trio"

        # Expected results for trio calculation depending on the trio_pedigree parameter
        if trio_pedigree is None or trio_pedigree == "":
            expected = {    
                "recessive": 1,
                "dominant": 5,
                "unknown": 1
            }
        else:
            expected = {
                "recessive": 1,
                "dominant": 6,
                "unknown": 0
            }

        # Create object
        variants = Variants(
            conn=None, input=input_vcf, output=output_vcf, param=param, load=True
        )

        # Calculation
        variants.calculation()

        # Check results

        result = variants.get_query_to_df(
            f""" SELECT INFO FROM variants WHERE INFO LIKE '%{tag}=recessive%' """
        )
        assert len(result) == expected.get("recessive"), f"With '{tag}' and '{trio_pedigree}', Expected {expected.get('recessive')} but got {len(result)}"

        result = variants.get_query_to_df(
            f""" SELECT * FROM variants WHERE INFO LIKE '%{tag}=dominant%' """
        )
        assert len(result) == expected.get("dominant"), f"With '{tag}' and '{trio_pedigree}', Expected {expected.get('dominant')} but got {len(result)}"

        result = variants.get_query_to_df(
            f""" SELECT * FROM variants WHERE INFO LIKE '%{tag}=unknown%' """
        )
        assert len(result) == expected.get("unknown"), f"With '{tag}' and '{trio_pedigree}', Expected {expected.get('unknown')} but got {len(result)}"

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
        param = {"calculation": {"calculations": {"vaf_normalization": None}}}

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


def test_calculation_vaf_normalization_ad():
    """
    This is a test function for the calculation of variant allele frequency normalization in a VCF file.
    """

    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # Init files
        input_vcf = tests_data_folder + "/example.AD.vcf"
        output_vcf = f"{tmp_dir}/output.vcf.gz"

        # Construct param dict
        param = {"calculation": {"calculations": {"vaf_normalization": None}}}

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

        result = variants.get_query_to_df(
            """ SELECT * FROM variants WHERE "#CHROM" = 'chr1' AND POS = 35144 AND sample4 LIKE '%:.' """
        )
        assert len(result) == 1

        result = variants.get_query_to_df(
            """ SELECT * FROM variants WHERE sample4 LIKE '%:.' """
        )
        assert len(result) == 4

        result = variants.get_query_to_df(
            """ SELECT * FROM variants WHERE sample4 NOT LIKE '%:.' """
        )
        assert len(result) == 3

        result = variants.get_query_to_df(""" SELECT * FROM variants""")
        # log.debug(result.to_string())

        # Check if VCF is in correct format with pyVCF
        remove_if_exists([output_vcf])
        variants.export_output()
        try:
            vcf.Reader(filename=output_vcf)
        except:
            assert False


def test_calculation_vaf_normalization_empty():
    """
    This is a test function for the calculation of variant allele frequency normalization in a VCF file.
    """

    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # Init files
        input_vcf = tests_data_folder + "/example.empty.vcf"
        output_vcf = f"{tmp_dir}/output.vcf.gz"

        # Construct param dict
        param = {"calculation": {"calculations": {"vaf_normalization": None}}}

        # Create object
        variants = Variants(
            conn=None, input=input_vcf, output=output_vcf, param=param, load=True
        )

        # Calculation
        variants.calculation()

        result = variants.get_query_to_df(
            """ SELECT INFO FROM variants WHERE FORMAT LIKE '%:VAF' """
        )
        assert len(result) == 0

        # Check if VCF is in correct format with pyVCF
        remove_if_exists([output_vcf])
        variants.export_output()
        try:
            vcf.Reader(filename=output_vcf)
        except:
            assert False


@pytest.mark.parametrize(
    "infos, vaf_normalization",
    [
        (infos, vaf_normalization)
        for infos in [None, "VAF", "DP", ["GQ"], ["GQ", "DP"], ["GQ", "DP", "VAF"], "GQ,DP"]
        for vaf_normalization in [True, False]
    ],
)
def test_calculation_genotype_stats(infos, vaf_normalization):
    """
    This is a test function for the calculation of genotype statistics in a VCF
    file using the Variants class in Python.
    """

    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # Init files
        input_vcf = tests_data_folder + "/example.vcf.gz"
        output_vcf = f"{tmp_dir}/output.vcf.gz"

        # Construct param dict
        # Use VAF normalization to calculate statistics on VAF
        if vaf_normalization:
            param = {
                "calculation": {
                    "calculations": {
                        "vaf_normalization": None,
                        "genotype_stats": {"infos": infos}
                    }
                }
            }
        else:
            param = {
                "calculation": {
                    "calculations": {
                        "genotype_stats": {"infos": infos}
                    }
                }
            }

        # Create object
        variants = Variants(
            conn=None, input=input_vcf, output=output_vcf, param=param, load=True
        )

        # Expected info tags
        if infos is None:
            infos = ["VAF"]

        # Clean up infos parameter to ensure it's a list
        if isinstance(infos, str):
            infos = infos.split(",")

        # Expected results for genotype_stats calculation depending on the infos parameter
        expected = {
            "VAF": {
                "nb": 4,
                "min": 0.279835,
                "max": 0.303819,
                "mean": 0.28737675
            } if vaf_normalization else {
                "nb": 0
            },
            "GQ": {
                "nb": 4,
                "min": 99.0,
                "max": 99.0,
                "mean": 99.0
            },
            "DP": {
                "nb": 4,
                "min": 576.0,
                "max": 17664.0,
                "mean": 9158.0
            }
        }

        # Calculation
        variants.calculation()

        for info in infos:

            # Stats in INFO field
            result = variants.get_query_to_df(
                f""" SELECT INFO FROM variants WHERE INFO LIKE '%{info}_stats%' """
            )
            assert len(result) == 7, f"With '{info}', Expected 7 but got {len(result)}"

            # Stats values NB in INFO field at least
            result = variants.get_query_to_df(
                f""" SELECT * FROM variants WHERE "#CHROM" = 'chr1' AND POS = 28736 AND INFO LIKE '%{info}_stats_nb={expected.get(info, {}).get("nb")}%' """
            )
            assert len(result) == 1, f"With '{info}' and {vaf_normalization}, Expected {expected.get(info, {}).get('nb')} but got {len(result)}"

            # If NB > 0, check min, max and mean values in INFO field
            if expected.get(info, {}).get("nb") > 0:

                # Stats values MIN in INFO field
                result = variants.get_query_to_df(
                    f""" SELECT * FROM variants WHERE "#CHROM" = 'chr1' AND POS = 28736 AND INFO LIKE '%{info}_stats_min={expected.get(info, {}).get("min")}%' """
                )
                assert len(result) == 1, f"With '{info}' and {vaf_normalization}, Expected {expected.get(info, {}).get('min')} but got {len(result)}"

                # Stats values MAX in INFO field
                result = variants.get_query_to_df(
                    f""" SELECT * FROM variants WHERE "#CHROM" = 'chr1' AND POS = 28736 AND INFO LIKE '%{info}_stats_max={expected.get(info, {}).get("max")}%' """
                )
                assert len(result) == 1, f"With '{info}' and {vaf_normalization}, Expected {expected.get(info, {}).get('max')} but got {len(result)}"

                # Stats values MEAN in INFO field
                result = variants.get_query_to_df(
                    f""" SELECT * FROM variants WHERE "#CHROM" = 'chr1' AND POS = 28736 AND INFO LIKE '%{info}_stats_mean={expected.get(info, {}).get("mean")}%' """
                )
                assert len(result) == 1, f"With '{info}' and {vaf_normalization}, Expected {expected.get(info, {}).get('mean')} but got {len(result)}"

        # Check if VCF is in correct format with pyVCF
        remove_if_exists([output_vcf])
        variants.export_output()
        try:
            vcf.Reader(filename=output_vcf)
        except:
            assert False



@pytest.mark.parametrize(
    "options, expect_in_info_fields, expect_not_in_info_fields, expected_in_format_fields, expected_not_in_format_fields, pattern_format_column, pattern_samples_column_variant1, pattern_samples_column_variant2, pattern_samples_column_variant3, pattern_samples_column_variant4",
    [
        (
            None,
            ["CLNSIG", "DP"],
            [],
            ["DP"],
            ["CLNSIG"],
            "GT:AD:DP:GQ",
            {
                "sample1": "0/1:525,204:729:99",
                "sample2": "0/1:12659,4994:17664:99",
                "sample3": "1/1:12658,4995:17663:99",
                "sample4": "1/1:401,175:576:99"
            },
            {
                "sample1": "./.:.:.:.",
                "sample2": "0/1:12659,4994:17664:99",
                "sample3": "0/1:12658,4995:17663:99",
                "sample4": "0/1:401,175:576:99"
            },
            {
                "sample1": "0/1:525,204:729:99",
                "sample2": "./.:.:.:.",
                "sample3": "0/1:12658,4995:17663:99",
                "sample4": "0/1:401,175:576:99"
            },
            {
                "sample1": "0/1:525,204:729:99",
                "sample2": "0/1:12659,4994:17664:99",
                "sample3": "./.:.:.:.",
                "sample4": "0/1:401,175:576:99"
            }
        ),
        (
            {},
            ["CLNSIG", "DP"],
            [],
            ["DP"],
            ["CLNSIG"],
            "GT:AD:DP:GQ",
            {
                "sample1": "0/1:525,204:729:99",
                "sample2": "0/1:12659,4994:17664:99",
                "sample3": "1/1:12658,4995:17663:99",
                "sample4": "1/1:401,175:576:99"
            },
            {
                "sample1": "./.:.:.:.",
                "sample2": "0/1:12659,4994:17664:99",
                "sample3": "0/1:12658,4995:17663:99",
                "sample4": "0/1:401,175:576:99"
            },
            {
                "sample1": "0/1:525,204:729:99",
                "sample2": "./.:.:.:.",
                "sample3": "0/1:12658,4995:17663:99",
                "sample4": "0/1:401,175:576:99"
            },
            {
                "sample1": "0/1:525,204:729:99",
                "sample2": "0/1:12659,4994:17664:99",
                "sample3": "./.:.:.:.",
                "sample4": "0/1:401,175:576:99"
            }
        ),
        (
            {
                "annotation_fields": {},
                "samples": [],
                "remove_info_fields": False
            },
            ["CLNSIG", "DP"],
            [],
            ["DP"],
            ["CLNSIG"],
            "GT:AD:DP:GQ",
            {
                "sample1": "0/1:525,204:729:99",
                "sample2": "0/1:12659,4994:17664:99",
                "sample3": "1/1:12658,4995:17663:99",
                "sample4": "1/1:401,175:576:99"
            },
            {
                "sample1": "./.:.:.:.",
                "sample2": "0/1:12659,4994:17664:99",
                "sample3": "0/1:12658,4995:17663:99",
                "sample4": "0/1:401,175:576:99"
            },
            {
                "sample1": "0/1:525,204:729:99",
                "sample2": "./.:.:.:.",
                "sample3": "0/1:12658,4995:17663:99",
                "sample4": "0/1:401,175:576:99"
            },
            {
                "sample1": "0/1:525,204:729:99",
                "sample2": "0/1:12659,4994:17664:99",
                "sample3": "./.:.:.:.",
                "sample4": "0/1:401,175:576:99"
            }
        ),
        (
            {
                "samples": ["sample5"],
            },
            ["CLNSIG", "DP"],
            [],
            ["DP"],
            ["CLNSIG"],
            "GT:AD:DP:GQ",
            {
                "sample1": "0/1:525,204:729:99",
                "sample2": "0/1:12659,4994:17664:99",
                "sample3": "1/1:12658,4995:17663:99",
                "sample4": "1/1:401,175:576:99"
            },
            {
                "sample1": "./.:.:.:.",
                "sample2": "0/1:12659,4994:17664:99",
                "sample3": "0/1:12658,4995:17663:99",
                "sample4": "0/1:401,175:576:99"
            },
            {
                "sample1": "0/1:525,204:729:99",
                "sample2": "./.:.:.:.",
                "sample3": "0/1:12658,4995:17663:99",
                "sample4": "0/1:401,175:576:99"
            },
            {
                "sample1": "0/1:525,204:729:99",
                "sample2": "0/1:12659,4994:17664:99",
                "sample3": "./.:.:.:.",
                "sample4": "0/1:401,175:576:99"
            }
        ),
        (
            {
                "annotation_fields": {
                    "CLNSIG": "CLINVAR",
                },
                "samples": [],
                "remove_info_fields": False
            },
            ["CLNSIG", "DP"],
            ["CLINVAR"],
            ["CLINVAR", "DP"],
            ["CLNSIG"],
            "GT:AD:DP:GQ:CLINVAR",
            {
                "sample1": "0/1:525,204:729:99:pathogenic",
                "sample2": "0/1:12659,4994:17664:99:pathogenic",
                "sample3": "1/1:12658,4995:17663:99:pathogenic",
                "sample4": "1/1:401,175:576:99:pathogenic"
            },
            {
                "sample1": "./.:.:.:.:non-pathogenic",
                "sample2": "0/1:12659,4994:17664:99:non-pathogenic",
                "sample3": "0/1:12658,4995:17663:99:non-pathogenic",
                "sample4": "0/1:401,175:576:99:non-pathogenic"
            },
            {
                "sample1": "0/1:525,204:729:99:.",
                "sample2": "./.:.:.:.:.",
                "sample3": "0/1:12658,4995:17663:99:.",
                "sample4": "0/1:401,175:576:99:."
            },
            {
                "sample1": "0/1:525,204:729:99:.",
                "sample2": "0/1:12659,4994:17664:99:.",
                "sample3": "./.:.:.:.:.",
                "sample4": "0/1:401,175:576:99:."
            }
        ),
        (
            {
                "annotation_fields": {
                    "CLNSIG": "CLINVAR",
                },
                "samples": ["sample1", "sample2", "sample3", "sample4", "sample5"],
                "remove_info_fields": False
            },
            ["CLNSIG", "DP"],
            ["CLINVAR"],
            ["CLINVAR", "DP"],
            ["CLNSIG"],
            "GT:AD:DP:GQ:CLINVAR",
            {
                "sample1": "0/1:525,204:729:99:pathogenic",
                "sample2": "0/1:12659,4994:17664:99:pathogenic",
                "sample3": "1/1:12658,4995:17663:99:pathogenic",
                "sample4": "1/1:401,175:576:99:pathogenic"
            },
            {
                "sample1": "./.:.:.:.:non-pathogenic",
                "sample2": "0/1:12659,4994:17664:99:non-pathogenic",
                "sample3": "0/1:12658,4995:17663:99:non-pathogenic",
                "sample4": "0/1:401,175:576:99:non-pathogenic"
            },
            {
                "sample1": "0/1:525,204:729:99:.",
                "sample2": "./.:.:.:.:.",
                "sample3": "0/1:12658,4995:17663:99:.",
                "sample4": "0/1:401,175:576:99:."
            },
            {
                "sample1": "0/1:525,204:729:99:.",
                "sample2": "0/1:12659,4994:17664:99:.",
                "sample3": "./.:.:.:.:.",
                "sample4": "0/1:401,175:576:99:."
            }
        ),
        (
            {
                "annotation_fields": {},
                "samples": ["sample1", "sample2"],
                "remove_info_fields": False
            },
            ["CLNSIG", "DP"],
            ["CLINVAR"],
            ["DP"],
            ["CLNSIG"],
            "GT:AD:DP:GQ",
            {
                "sample1": "0/1:525,204:729:99",
                "sample2": "0/1:12659,4994:17664:99",
                "sample3": "1/1:12658,4995:17663:99",
                "sample4": "1/1:401,175:576:99"
            },
            {
                "sample1": "./.:.:.:.",
                "sample2": "0/1:12659,4994:17664:99",
                "sample3": "0/1:12658,4995:17663:99",
                "sample4": "0/1:401,175:576:99"
            },
            {
                "sample1": "0/1:525,204:729:99",
                "sample2": "./.:.:.:.",
                "sample3": "0/1:12658,4995:17663:99",
                "sample4": "0/1:401,175:576:99"
            },
            {
                "sample1": "0/1:525,204:729:99",
                "sample2": "0/1:12659,4994:17664:99",
                "sample3": "./.:.:.:.",
                "sample4": "0/1:401,175:576:99"
            }
        ),
        (
            {
                "annotation_fields": {},
                "samples": ["sample1", "sample2"],
                "remove_info_fields": True
            },
            ["CLNSIG", "DP"],
            ["CLINVAR"],
            ["DP"],
            ["CLNSIG"],
            "GT:AD:DP:GQ",
            {
                "sample1": "0/1:525,204:729:99",
                "sample2": "0/1:12659,4994:17664:99",
                "sample3": "1/1:12658,4995:17663:99",
                "sample4": "1/1:401,175:576:99"
            },
            {
                "sample1": "./.:.:.:.",
                "sample2": "0/1:12659,4994:17664:99",
                "sample3": "0/1:12658,4995:17663:99",
                "sample4": "0/1:401,175:576:99"
            },
            {
                "sample1": "0/1:525,204:729:99",
                "sample2": "./.:.:.:.",
                "sample3": "0/1:12658,4995:17663:99",
                "sample4": "0/1:401,175:576:99"
            },
            {
                "sample1": "0/1:525,204:729:99",
                "sample2": "0/1:12659,4994:17664:99",
                "sample3": "./.:.:.:.",
                "sample4": "0/1:401,175:576:99"
            }
        ),
        (
            {
                "annotation_fields": {
                    "CLNSIG": None,
                    "DP": "DP2"
                },
                "samples": ["sample1", "sample2"],
                "remove_info_fields": True
            },
            [],
            ["CLNSIG", "DP"],
            ["CLNSIG", "DP", "DP2"],
            [],
            "GT:AD:DP:GQ:CLNSIG:DP2",
            {
                "sample1": "0/1:525,204:729:99:pathogenic:.",
                "sample2": "0/1:12659,4994:17664:99:pathogenic:.",
                "sample3": "1/1:12658,4995:17663:99:.:.",
                "sample4": "1/1:401,175:576:99:.:."
            },
            {
                "sample1": "./.:.:.:.:non-pathogenic:.",
                "sample2": "0/1:12659,4994:17664:99:non-pathogenic:.",
                "sample3": "0/1:12658,4995:17663:99:.:.",
                "sample4": "0/1:401,175:576:99:.:."
            },
            {
                "sample1": "0/1:525,204:729:99:.:.",
                "sample2": "./.:.:.:.:.:.",
                "sample3": "0/1:12658,4995:17663:99:.:.",
                "sample4": "0/1:401,175:576:99:.:."
            },
            {
                "sample1": "0/1:525,204:729:99:.:125",
                "sample2": "0/1:12659,4994:17664:99:.:125",
                "sample3": "./.:.:.:.:.:.",
                "sample4": "0/1:401,175:576:99:.:."
            }
        )
    ],
)
def test_calculation_info_to_format(options, expect_in_info_fields, expect_not_in_info_fields, expected_in_format_fields, expected_not_in_format_fields, pattern_format_column, pattern_samples_column_variant1, pattern_samples_column_variant2, pattern_samples_column_variant3, pattern_samples_column_variant4):
    """
    This is a test function for the calculation of INFO to FORMAT conversion in a VCF file.
    """

    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # Init files
        input_vcf = tests_data_folder + "/example.vcf.gz"
        output_vcf = f"{tmp_dir}/output.vcf.gz"

        # Construct param dict
        param = {
            "calculation": {
                "calculations": {
                    "INFO_TO_FORMAT": options
                }
            }
        }

        # Create object
        variants = Variants(
            conn=None, input=input_vcf, output=output_vcf, param=param, load=True
        )

        # Calculation
        variants.calculation()

        # Test if INFO fields are present or not in header of variants
        for field in expect_in_info_fields:
            assert field in variants.get_header().infos, f"Field {field} should be present in INFO header"
        for field in expect_not_in_info_fields:
            assert field not in variants.get_header().infos, f"Field {field} should not be present in INFO header"
        
        # Test if FORMAT fields are present or not in header of variants
        for field in expected_in_format_fields:
            assert field in variants.get_header().formats, f"Field {field} should be present in FORMAT header"
        for field in expected_not_in_format_fields:
            assert field not in variants.get_header().formats, f"Field {field} should not be present in FORMAT header"

        # Test if FORMAT column is correct for all variants
        result = variants.get_query_to_df(
            """ SELECT FORMAT FROM variants """
        )
        assert len(result) == 7
        for _, row in result.iterrows():
            assert row["FORMAT"] == pattern_format_column, f"FORMAT column value {row['FORMAT']} does not match expected {pattern_format_column}"

        # Test if sample columns are correct for variant 1
        log.debug(f"Testing sample columns for variant 1 at chr1:28736")
        result = variants.get_query_to_df(
            """ SELECT * FROM variants WHERE "#CHROM" = 'chr1' AND POS = 28736 """
        )
        assert len(result) == 1
        for _, row in result.iterrows():
            for sample in pattern_samples_column_variant1:
                assert row[sample] == pattern_samples_column_variant1[sample], f"Sample {sample} value {row[sample]} does not match expected {pattern_samples_column_variant1[sample]}"

        # Test if sample columns are correct for variant 2
        log.debug(f"Testing sample columns for variant 2 at chr1:35144")
        result = variants.get_query_to_df(
            """ SELECT * FROM variants WHERE "#CHROM" = 'chr1' AND POS = 35144 """
        )
        assert len(result) == 1
        for _, row in result.iterrows():
            for sample in pattern_samples_column_variant2:
                assert row[sample] == pattern_samples_column_variant2[sample], f"Sample {sample} value {row[sample]} does not match expected {pattern_samples_column_variant2[sample]}"

        # Test if sample columns are correct for variant 3
        log.debug(f"Testing sample columns for variant 3 at chr1:768251")
        result = variants.get_query_to_df(
            """ SELECT * FROM variants WHERE "#CHROM" = 'chr1' AND POS = 768251 """
        )
        assert len(result) == 1
        for _, row in result.iterrows():
            for sample in pattern_samples_column_variant3:
                assert row[sample] == pattern_samples_column_variant3[sample], f"Sample {sample} value {row[sample]} does not match expected {pattern_samples_column_variant3[sample]}"

        # Test if sample columns are correct for variant 4
        log.debug(f"Testing sample columns for variant 4 at chr7:55249063")
        result = variants.get_query_to_df(
            """ SELECT * FROM variants WHERE "#CHROM" = 'chr7' AND POS = 55249063 """
        )
        assert len(result) == 1
        for _, row in result.iterrows():
            for sample in pattern_samples_column_variant4:
                assert row[sample] == pattern_samples_column_variant4[sample], f"Sample {sample} value {row[sample]} does not match expected {pattern_samples_column_variant4[sample]}"

        # Check if VCF is in correct format with pyVCF
        remove_if_exists([output_vcf])
        variants.export_output()
        try:
            vcf.Reader(filename=output_vcf)
        except:
            assert False

        # Print content of vcf file
        import gzip
        with gzip.open(output_vcf, "rt") as f:
            for line in f:
                log.debug(line.strip())


def test_calculation_info_to_format_flag_info_field():
    """
    This is a test function for the calculation of INFO to FORMAT conversion in a VCF file.
    """

    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # Init files
        input_vcf = tests_data_folder + "/example.flag.vcf"
        output_vcf = f"{tmp_dir}/output.vcf"

        # Construct param dict
        param = {
            "calculation": {
                "calculations": {
                    "INFO_TO_FORMAT": {
                        "annotation_fields": {
                            "FLAG": None
                        },
                        "samples": [],
                        "remove_info_fields": False
                    }
                }
            }
        }

        # Create object
        variants = Variants(
            conn=None, input=input_vcf, output=output_vcf, param=param, load=True
        )

        # Calculation
        try:
            variants.calculation()
            assert False, "Expected NotImplementedError for flag INFO field conversion to FORMAT"
        except NotImplementedError as e:
            log.debug(f"Calculation not implemented: {e}")
            assert True, "Expected NotImplementedError for flag INFO field conversion to FORMAT"


@pytest.mark.parametrize(
    "calculation_name, variant_id_tag, variant_id_tag_info, keep_variant_id_tag_column",
    [
        (
            calculation_name,
            variant_id_tag,
            variant_id_tag_info,
            keep_variant_id_tag_column,
        )
        for calculation_name in ["variant_id"]
        for variant_id_tag in [
            None,
            "variant_id",
            "var_id",
            "id",
            "variantID",
            "varID",
        ]
        for variant_id_tag_info in [None, "HOWARD variant ID generation"]
        for keep_variant_id_tag_column in [False, True]
    ],
)
def test_calculation_variant_id(
    calculation_name, variant_id_tag, variant_id_tag_info, keep_variant_id_tag_column
):
    """
    This is a test function for the calculation of variant IDs in a VCF file using the Variants class in Python.
    """

    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # Init files
        input_vcf = tests_data_folder + "/example.vcf.gz"
        output_vcf = f"{tmp_dir}/output.vcf.gz"

        # Check if variant_id_tag is None and set default tag name and params
        if variant_id_tag is None:
            variant_id_tag_param = {}
            if calculation_name == "variant_id":
                variant_id_tag = "variant_id"
        else:
            variant_id_tag_param = {"variant_id_tag": variant_id_tag}

        # Add info and keep column param
        variant_id_tag_param["variant_id_tag_info"] = variant_id_tag_info
        variant_id_tag_param["keep_variant_id_tag_column"] = keep_variant_id_tag_column

        # Construct param dict
        param = {
            "calculation": {"calculations": {calculation_name: variant_id_tag_param}}
        }  # test custom tag name

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
            f""" SELECT INFO FROM variants WHERE INFO LIKE '%{variant_id_tag}%' """
        )
        assert len(result) == 7

        # Explode info
        variants.explode_infos(prefix="INFO/", fields=[variant_id_tag])

        # Check distinct variant_id
        result = variants.get_query_to_df(
            f""" SELECT distinct "INFO/{variant_id_tag}" FROM variants """
        )
        assert len(result) == 7

        # Check if variant ID column exists or not
        if keep_variant_id_tag_column:
            result = variants.get_query_to_df(
                f""" SELECT "{variant_id_tag}" FROM variants """
            )
            assert len(result) == 7
        else:
            try:
                result = variants.get_query_to_df(
                    f""" SELECT "{variant_id_tag}" FROM variants """
                )
                assert False
            except:
                assert True

        # Check if VCF is in correct format with pyVCF
        remove_if_exists([output_vcf])
        variants.export_output()
        try:
            vcf.Reader(filename=output_vcf)
        except:
            assert False

@pytest.mark.parametrize(
    "input_vcf, filter_name, where_clause, expected_nb_variant, sample_list, genotype_filter",
    [
        (
            tests_data_folder + "/example.vcf.gz",
            "No filter on variants",
            None,
            7,
            [],
            False
        ),
        (
            tests_data_folder + "/example.vcf.gz",
            "No filter on variants with empty where clause",
            "",
            7,
            [],
            False
        ),
        (
            tests_data_folder + "/example.vcf.gz",
            "Filter on chromosome 1",
            "\"#CHROM\" = 'chr1'",
            6,
            [],
            False
        ),
        (
            tests_data_folder + "/example.vcf.gz",
            "Filter on CLNSIG with pathogenic value",
            "INFOS.CLNSIG = 'pathogenic'",
            1,
            [],
            False
        ),
        (
            tests_data_folder + "/example.vcf.gz",
            "Filter on sample1 DP with value greater than 100",
            "SAMPLES.sample1.DP >= 100 AND SAMPLES.sample2.DP >= 100",
            2,
            [],
            False
        ),
        (
            tests_data_folder + "/example.vcf.gz",
            "Filter on samples with only sample1 and sample2",
            "",
            7,
            ["sample1", "sample2"],
            False
        ),
        (
            tests_data_folder + "/example.vcf.gz",
            "Filter on samples with only sample1 and sample2, with sample list as string",
            "",
            7,
            "sample1,sample2",
            False
        ),
        (
            tests_data_folder + "/example.vcf.gz",
            "Filter on chromosome 1 and on samples with only sample1 and sample2, with sample list as string",
            "\"#CHROM\" = 'chr1'",
            6,
            "sample1,sample2",
            False
        ),
        (
            tests_data_folder + "/example.vcf.gz",
            "Filter on chromosome 1 and on samples with only sample1 and sample2, with sample list as string",
            "\"#CHROM\" = 'chr1'",
            6,
            "sample2",
            False
        ),
        (
            tests_data_folder + "/example.vcf.gz",
            "Filter on chromosome 1 and on samples with only sample1 and sample2, with sample list as string, and with genotype filter",
            "\"#CHROM\" = 'chr1'",
            2,
            "sample2",
            True
        ),
        (
            tests_data_folder + "/example.vcf.gz",
            "Filter on chromosome 1 and on samples with only sample1 and sample2, with sample list as string",
            "\"#CHROM\" = 'chr1'",
            6,
            "sample1,sample2",
            False
        ),
        (
            tests_data_folder + "/example.with_allowed_genotypes.vcf",
            "Filter on chromosome 1 and on samples with only sample1 and sample2, with sample list as string, and with genotype filter",
            "\"#CHROM\" = 'chr1'",
            2,
            "sample2,sample3",
            True
        ),
        (
            tests_data_folder + "/example.with_allowed_genotypes.vcf",
            "Filter on samples with only sample2 and sample3, with sample list as string, and with genotype filter",
            "",
            3,
            "sample2,sample3",
            True
        ),
        (
        tests_data_folder + "/example.with_allowed_genotypes.vcf",
        "Filter only with genotype filter",
        "",
        7,
        None,
        True
    ),
    ],
)
def test_calculation_variant_filter(
    input_vcf, filter_name, where_clause, expected_nb_variant, sample_list, genotype_filter
):
    """
    This is a test function for the calculation of variant IDs in a VCF file using the Variants class in Python.
    """

    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # Init files
        output_vcf = f"{tmp_dir}/output.vcf.gz"

        # Construct param dict
        param = {
            "calculation_test": {
                "calculations": {
                    "variant_filter": {
                        "filter_name": filter_name,
                        "where_clause": where_clause,
                        "sample_list": sample_list,
                        "genotype_filter": genotype_filter
                    }
                }
            }
        }  # test custom tag name

        # Create object
        variants = Variants(
            conn=None, input=input_vcf, output=output_vcf, param=param, load=True
        )

        # Remove if output file exists
        remove_if_exists([output_vcf])

        # Calculation
        variants.calculation(section="calculation_test")

        # Check if all variant have variant_id
        result = variants.get_query_to_df(
            f""" SELECT * FROM variants """
        )
        assert len(result) == expected_nb_variant, f"Expected {expected_nb_variant} variants with filter_name '{filter_name}', but got {len(result)}"

        # Check if VCF is in correct format with pyVCF
        remove_if_exists([output_vcf])
        variants.export_output()
        try:
            vcf.Reader(filename=output_vcf)
        except:
            assert False

def test_calculation_annotation():
    """
    This is a test function for the calculation of annotation in a VCF file using the Variants class in Python.
    """

    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # Init files
        input_vcf = tests_data_folder + "/example.vcf.gz"
        output_vcf = f"{tmp_dir}/output.vcf.gz"

        # Annotation databases
        annotation1 = database_files.get("parquet")
        annotation2 = database_files.get("example_vcf_gz")
        annotation3 = database_files.get("refgene_gz")

        # Param annotation dict
        # Construct param dict
        param_annotation = {
            "parquet": {
                "annotations": {
                    annotation1: {"INFO": None},
                    annotation2: {"CLNSIG": "CLNSIG_new"},
                },
            },
            "bcftools": {
                "annotations": {
                    annotation2: {"CLNSIG": "CLNSIG_new_bcftools"},
                    annotation3: {"symbol": "gene"},
                },
            },
        }

        # Construct param dict
        param = {
            "calculation": {
                "calculations": {
                    "annotation": param_annotation
                }
            }
        }

        # Create object
        variants = Variants(
            conn=None, input=input_vcf, output=output_vcf, param=param, load=True
        )

        # Remove if output file exists
        remove_if_exists([output_vcf])

        # Calculation
        variants.calculation()

        result = variants.get_query_to_df(
            """ SELECT INFO FROM variants """
        )
        assert len(result) == 7, f"Expected 7 variants, but got {len(result)}"

        result = variants.get_query_to_df(
            """ SELECT INFO FROM variants WHERE INFO LIKE '%CLNSIG_new=%' """
        )
        assert len(result) == 2, f"Expected 2 variants with INFO LIKE '%CLNSIG_new=%', but got {len(result)}"

        result = variants.get_query_to_df(
            """ SELECT INFO FROM variants WHERE INFO LIKE '%CLNSIG_new_bcftools=%' """
        )
        assert len(result) == 2, f"Expected 2 variants with INFO LIKE '%CLNSIG_new_bcftools=%', but got {len(result)}"

        result = variants.get_query_to_df(
            """ SELECT INFO FROM variants WHERE INFO LIKE '%gene=%' """
        )
        assert len(result) == 3, f"Expected 3 variants with INFO LIKE '%gene=%', but got {len(result)}"

        # Check if VCF is in correct format with pyVCF
        remove_if_exists([output_vcf])
        variants.export_output()
        try:
            vcf.Reader(filename=output_vcf)
        except:
            assert False


def test_calculation_prioritization():
    """
    This is a test function for the calculation of prioritization in a VCF file using the Variants class in Python.
    """

    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # Init files
        input_vcf = tests_data_folder + "/example.vcf.gz"
        output_vcf = f"{tmp_dir}/output.vcf.gz"

        # Construct config dict
        config = {}

        # Construct param prioritization dict
        param_prioritization = {
            "prioritization_config": tests_data_folder
            + "/prioritization_profiles.json",
            "profiles": ["default", "GERMLINE", "sql_class"],
            "pzfields": [
                "PZFlag",
                "PZScore",
                "PZClass",
                "PZComment",
                "PZInfos",
                "PZTags",
            ],
        }

        # Construct param dict
        param = {
            "calculation": {
                "calculations": {
                    "prioritization": param_prioritization
                }
            }
        }

        # Create object
        variants = Variants(
            input=input_vcf, output=output_vcf, load=True, config=config, param=param
        )

        # Remove if output file exists
        remove_if_exists([output_vcf])

        # Calculation
        variants.calculation()

        # Check if all variants have INFO field
        result = variants.get_query_to_df(
            """ SELECT INFO FROM variants"""
        )
        assert len(result) == 7, f"Expected 7 variants, but got {len(result)}"

        # Check all priorized default profile
        result = variants.get_query_to_df(
            """
            SELECT * FROM variants
            WHERE INFO LIKE '%PZFlag_default=%'
            AND INFO LIKE '%PZScore_default=%'
            AND INFO LIKE '%PZClass_default=%'
            AND INFO LIKE '%PZComment_default=%'
            AND INFO LIKE '%PZInfos_default=%'
            AND INFO LIKE '%PZTags_default=%'
            """
        )
        assert len(result) == 4, f"Expected 4 variants with default profile, but got {len(result)}"

        # Check all priorized GERMLINE profile
        result = variants.get_query_to_df(
            """
            SELECT * FROM variants
            WHERE INFO LIKE '%PZFlag_GERMLINE=%'
            AND INFO LIKE '%PZScore_GERMLINE=%'
            AND INFO LIKE '%PZClass_GERMLINE=%'
            AND INFO LIKE '%PZComment_GERMLINE=%'
            AND INFO LIKE '%PZInfos_GERMLINE=%'
            AND INFO LIKE '%PZTags_GERMLINE=%'
            """
        )
        assert len(result) == 2, f"Expected 2 variants with GERMLINE profile, but got {len(result)}"

        # Check all priorized sql_class profile
        result = variants.get_query_to_df(
            """
            SELECT * FROM variants
            WHERE INFO LIKE '%PZFlag_sql_class=%'
            AND INFO LIKE '%PZScore_sql_class=%'
            AND INFO LIKE '%PZClass_sql_class=%'
            AND INFO LIKE '%PZComment_sql_class=%'
            AND INFO LIKE '%PZInfos_sql_class=%'
            """
        )
        assert len(result) == 2, f"Expected 2 variants with sql_class profile, but got {len(result)}"

        # Check all priorized default profile (as default)
        result = variants.get_query_to_df(
            """
            SELECT * FROM variants
            WHERE INFO LIKE '%PZFlag=%'
            AND INFO LIKE '%PZScore=%'
            AND INFO LIKE '%PZClass=%'
            AND INFO LIKE '%PZComment_default=%'
            AND INFO LIKE '%PZInfos_default=%'
            AND INFO LIKE '%PZTags_default=%'
            """
        )
        assert len(result) == 4, f"Expected 4 variants with default profile (as default), but got {len(result)}"

        # Check all priorized default profile
        result = variants.get_query_to_df(
            """
            SELECT * FROM variants
            WHERE INFO LIKE '%PZFlag_default=%'
            AND INFO LIKE '%PZScore_default=%'
            """
        )
        assert len(result) == 7, f"Expected 7 variants with default profile, but got {len(result)}"

        # Check all priorized GERMLINE profile
        result = variants.get_query_to_df(
            """
            SELECT * FROM variants
            WHERE INFO LIKE '%PZFlag_GERMLINE=%'
            AND INFO LIKE '%PZScore_GERMLINE=%'
            """
        )
        assert len(result) == 7, f"Expected 7 variants with GERMLINE profile, but got {len(result)}"

        # Check all priorized default profile (as default)
        result = variants.get_query_to_df(
            """
            SELECT * FROM variants
            WHERE INFO LIKE '%PZFlag=%'
            AND INFO LIKE '%PZScore=%'
            """
        )
        assert len(result) == 7, f"Expected 7 variants with default profile (as default), but got {len(result)}"

        # Check annotation default
        result = variants.get_query_to_df(
            """ SELECT * FROM variants WHERE "#CHROM" = 'chr7' AND POS = 55249063 AND REF = 'G' AND ALT = 'A' AND INFO LIKE '%PZScore_default=105%' """
        )
        assert len(result) == 1, f"Expected 1 variant with PZScore_default=105, but got {len(result)}"

        # Check annotation GERMILNE
        result = variants.get_query_to_df(
            """ SELECT * FROM variants WHERE "#CHROM" = 'chr7' AND POS = 55249063 AND REF = 'G' AND ALT = 'A' AND INFO LIKE '%PZScore_GERMLINE=5%' """
        )
        assert len(result) == 1, f"Expected 1 variant with PZScore_GERMLINE=5, but got {len(result)}"

        # Check annotation sql_class
        result = variants.get_query_to_df(
            """ SELECT * FROM variants WHERE "#CHROM" = 'chr7' AND POS = 55249063 AND REF = 'G' AND ALT = 'A' AND INFO LIKE '%PZScore_sql_class=60%' """
        )
        assert len(result) == 1, f"Expected 1 variant with PZScore_sql_class=60, but got {len(result)}"

        # Check FILTERED
        result = variants.get_query_to_df(
            """ SELECT INFO FROM variants WHERE INFO LIKE '%FILTERED%' """
        )
        assert len(result) == 1, f"Expected 1 variant with FILTERED in INFO, but got {len(result)}"

        # Check if VCF is in correct format with pyVCF
        remove_if_exists([output_vcf])
        variants.export_output()
        try:
            vcf.Reader(filename=output_vcf)
        except:
            assert False


@pytest.mark.parametrize(
    "param",
    [
        (
            None
        ),
        (
            {}
        ),
        (
            {
                "stats_md": "/tmp/stats.md",
            }
        ),
        (
            {
                "stats_stdout": True,
            }
        ),
        (
            {
                "stats_md": "/tmp/stats.md",
                "stats_json": "/tmp/stats.json",
            }
        ),
        (
            {
                "stats_html": "/tmp/stats.html",
            }
        ),
        (
            {
                "stats_pdf": "/tmp/stats.pdf",
            }
        ),
        (
            {
            "stats_md": "/tmp/stats.md",
                "annotations_stats": True,
            }
        ),
        (
            {
                "stats_md": "/tmp/stats.md",
                "queries": {
                    "First variants": "SELECT \"#CHROM\", POS, REF, ALT FROM variants_view LIMIT 5",
                    "First INFO tags": "SELECT INFOS.* FROM variants_view LIMIT 5",
                }
            }
        ),
    ],
)
def test_calculation_export_stats(param):

    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # Init files
        input_vcf = tests_data_folder + "/example.vcf.gz"

        # Construct config dict
        config = {}

        # Construct param dict
        param = {
            "calculation": {
                "calculations": {
                    "EXPORT_STATS": param
                }
            }
        }

        # Create object
        variants = Variants(
            input=input_vcf, load=True, config=config, param=param
        )

        # Calculation
        variants.calculation()


@pytest.mark.parametrize(
    "input_file, output_file, export, explode, query",
    [
        (
            input_file,
            output_file,
            export,
            explode, query
        )
        for input_file in ["example.vcf.gz"]
        for output_file in [
            "variants.vcf",
            "variants.vcf.gz",
            "variants.tsv",
            "variants.parquet",
            "variants.partition.parquet",
            ]
        for export in [
                None,
                {},
                {
                    "fields_to_rename": {
                        "CLNSIG": "CLNSIG_renamed",
                        "SIFT": None
                    }
                },
                {
                    "order_by": "POS DESC"
                },
                {
                    "force_cast_as_flat": True
                },
                {
                    "include_header": True,
                    "order_by": "DP ASC",
                    "parquet_partitions": "#CHROM"
                }
            ]
        for explode in [
            None,
            {},
            {
                "explode_infos": True,
                "explode_infos_prefix": "test_",
                "explode_infos_fields": ["CLNSIG", "SIFT", "DP"],
            }
        ]
        for query in [
            None,
            "SELECT * FROM variants WHERE INFO LIKE '%CLNSIG%'",
        ]
    ]
)
def test_calculation_export_variants(input_file, output_file, export, explode, query):

    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        tmp_dir = "/tmp"

        # Init files
        input_vcf = tests_data_folder + f"/{input_file}"
        output_vcf = f"{tmp_dir}/{output_file}"

        # Construct config dict
        config = {}

         # Construct param dict
        param = {
            "calculation": {
                "calculations": {
                    "EXPORT_VARIANTS": {
                        "file": output_vcf,
                        "query": query,
                        "export": export,
                        "explode": explode
                    }
                }
            }
        }

        # Create object
        variants = Variants(
            input=input_vcf, load=True, config=config, param=param
        )

        # Remove if output file exists
        remove_if_exists([output_vcf])

        # Calculation
        variants.calculation()

        # Check output file exists
        assert os.path.exists(output_vcf), f"Output file {output_vcf} does not exist"

        # Check content of file with object 
        # Only if not partition due to VCF partition failed
        param_export = (((param.get("calculation") or {}).get("calculations") or {}).get("EXPORT_VARIANTS") or {}).get("export") or {}
        if param_export.get("parquet_partitions", None) is None:
            
            log.debug(f"Results from output file {output_vcf} with param {param}:")
            
            variants_output = Variants(
                input=output_vcf, load=True, config=config, param=param
            )

            results = variants_output.get_query_to_df(
                """ SELECT * FROM variants """
            )
            #log.debug(f"Results from output file {output_vcf}: {results.to_string()}")
            assert len(results) > 0, f"Expected variants in output file {output_vcf}, but got {len(results)}"