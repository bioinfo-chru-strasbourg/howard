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

from test_needed import tests_folder, tests_data_folder, tests_config


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

        # result = variants.get_query_to_df(f"""SELECT * FROM variants """)
        # log.debug(result.to_string())

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


def test_calculation_merged_hgvs():
    """
    This function tests the calculation of merging multiple HGVS annotations.
    """

    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # Init files
        # input_vcf = tests_data_folder + "/example.annotated.annovar.snpeff.vcf"
        input_vcf = tests_data_folder + "/example.annotated.hgvs.annovar.snpeff.vcf"
        output_vcf = f"{tmp_dir}/output.vcf.gz"
        input_param = {
            "calculation": {"calculations": {"MERGED_HGVS": None}},
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

        # Check number of NOMEN
        annovar_hgvs = "AAChange_refGene=EGFR:NM_001346941:exon14:c.1560G>A:p.Q520Q,EGFR:NM_001346897:exon19:c.2226G>A:p.Q742Q,EGFR:NM_001346899:exon19:c.2226G>A:p.Q742Q,EGFR:NM_001346898:exon20:c.2361G>A:p.Q787Q,EGFR:NM_001346900:exon20:c.2202G>A:p.Q734Q,EGFR:NM_005228:exon20:c.2361G>A:p.Q787Q"
        snpeff_hgvs = "snpeff_hgvs=EGFR:NM_005228.5:exon20:c.2361G>A:p.Gln787Gln,EGFR:NM_001346897.2:exon19:c.2226G>A:p.Gln742Gln,EGFR:NM_001346898.2:exon20:c.2361G>A:p.Gln787Gln,EGFR:NM_001346941.2:exon14:c.1560G>A:p.Gln520Gln,EGFR:NM_001346899.1:exon19:c.2226G>A:p.Gln742Gln,EGFR:NM_001346900.2:exon20:c.2202G>A:p.Gln734Gln,EGFR-AS1:NR_047551.1:exon2:n.1201C>T"
        merged_hgvs = "merged_hgvs=EGFR:NM_005228:exon20:c.2361G>A:p.Q787Q,EGFR:NM_001346898:exon20:c.2361G>A:p.Q787Q,EGFR:NM_001346898.2:exon20:c.2361G>A:p.Gln787Gln,EGFR:NM_001346899:exon19:c.2226G>A:p.Q742Q,EGFR:NM_001346900.2:exon20:c.2202G>A:p.Gln734Gln,EGFR:NM_001346941:exon14:c.1560G>A:p.Q520Q,EGFR:NM_001346900:exon20:c.2202G>A:p.Q734Q,EGFR:NM_005228.5:exon20:c.2361G>A:p.Gln787Gln,EGFR-AS1:NR_047551.1:exon2:n.1201C>T,EGFR:NM_001346899.1:exon19:c.2226G>A:p.Gln742Gln,EGFR:NM_001346941.2:exon14:c.1560G>A:p.Gln520Gln,EGFR:NM_001346897:exon19:c.2226G>A:p.Q742Q,EGFR:NM_001346897.2:exon19:c.2226G>A:p.Gln742Gln"
        result = variants.get_query_to_df(
            f"""SELECT INFO FROM variants WHERE INFO LIKE '%{annovar_hgvs}%' AND INFO LIKE '%{snpeff_hgvs}%' AND INFO LIKE '%{merged_hgvs}%' """
        )
        assert len(result) == 1

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
    "input_file, calculation",
    [
        (
            input_file,
            calculation,
        )
        for input_file in ["example.ann.vcf.gz", "example.chrM.ann.vcf.gz"]
        for calculation in [
            "snpeff_extract",
            "snpeff_hgvs",
            "snpeff_ann_explode",
            "snpeff_ann_explode_uniquify",
            "snpeff_ann_explode_json",
        ]
    ],
)
def test_calculation_snpeff_ann_explode(input_file, calculation):
    """
    This function is a test for calculating snpeff_hgvs in a VCF file using the Variants class.

    :param calculation: It looks like the `calculation` parameter is used to specify a particular
    calculation method in the test function `test_calculation_snpeff_ann_explode`. This parameter is
    then used to construct a parameter dictionary `param` with the calculation method specified
    """

    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # Init files
        input_vcf = tests_data_folder + "/" + input_file
        output_vcf = f"{tmp_dir}/output.{calculation}.vcf"

        # Construct param dict
        param = {"calculation": {"calculations": {calculation: None}}}

        # Create object
        variants = Variants(
            conn=None, input=input_vcf, output=output_vcf, param=param, load=True
        )

        # Check nb variants
        result = variants.get_query_to_df(""" SELECT INFO FROM variants""")
        nb_variants = len(result)

        # Check if no snpeff_hgvs
        result = variants.get_query_to_df(
            """ SELECT INFO FROM variants WHERE regexp_matches(INFO,'snpeff[^=]') """
        )
        assert len(result) == 0

        # Calculation
        variants.calculation()

        # query annotated variant (0 if no snpeff annotation like in chrM)
        result = variants.get_query_to_df(
            """ SELECT * FROM variants WHERE regexp_matches(INFO,'snpeff[^=]') """
        )
        assert len(result) == nb_variants or len(result) == 0

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
        param = {"calculation": {"calculations": {"COUNT_SAMPLES": None}}}

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


def test_calculation_find_samples_count_samples_options_in_calculation():
    """
    This is a test function for the "FIND_SAMPLES" calculation in the Variants class, which checks if
    the calculation is performed correctly and the output VCF file is in the correct format.
    """

    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # Init files
        input_vcf = tests_data_folder + "/example.vcf.gz"
        output_vcf = f"{tmp_dir}/output.vcf.gz"

        # Construct param dict
        param = {"calculation": {"calculations": {"COUNT_SAMPLES": {"tags": {"count_samples_for_variant": "count"}}}}}

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


def test_calculation_find_samples_list_samples():
    """
    This is a test function for the "FIND_SAMPLES" calculation in the Variants class, which checks if
    the calculation is performed correctly and the output VCF file is in the correct format.
    """

    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # Init files
        input_vcf = tests_data_folder + "/example.vcf.gz"
        output_vcf = f"{tmp_dir}/output.vcf.gz"

        # Construct param dict
        param = {"calculation": {"calculations": {"LIST_SAMPLES": None}}}

        # Create object
        variants = Variants(
            conn=None, input=input_vcf, output=output_vcf, param=param, load=True
        )

        # Calculation
        variants.calculation()

        result = variants.get_query_to_df(
            """ SELECT INFO FROM variants WHERE INFO NOT LIKE '%count_samples%' AND INFO LIKE '%list_samples%' """
        )
        assert len(result) == 7

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


def test_calculation_barcode_sample_name_special_char():
    """
    This is a test function for a Python script that calculates barcode information from a VCF file and
    checks if the output is correct.
    """

    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # Init files
        input_vcf = tests_data_folder + "/example.name_with_special_char.vcf"
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
        # log.debug(result)

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


def test_calculation_vaf_normalization_ad():
    """
    This is a test function for the calculation of variant allele frequency normalization in a VCF file.
    """

    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # Init files
        input_vcf = tests_data_folder + "/example.AD.vcf"
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
        assert len(result) == 0

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


def test_calculation_genotype_stats():
    """
    This is a test function for the calculation of genotype statistics in a VCF
    file using the Variants class in Python.
    """

    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # Init files
        input_vcf = tests_data_folder + "/example.vcf.gz"
        output_vcf = f"{tmp_dir}/output.vcf.gz"

        # Construct param dict
        param = {"calculation": {"calculations": {"genotype_stats": {"info": "GQ"}}}}

        # Create object
        variants = Variants(
            conn=None, input=input_vcf, output=output_vcf, param=param, load=True
        )

        # Calculation
        variants.calculation()

        result = variants.get_query_to_df(
            """ SELECT INFO FROM variants WHERE INFO LIKE '%GQ_stats%' """
        )
        assert len(result) == 7

        result = variants.get_query_to_df(
            """ SELECT * FROM variants WHERE "#CHROM" = 'chr1' AND POS = 28736 AND INFO LIKE '%GQ_stats_nb=4%' """
        )
        assert len(result) == 1

        result = variants.get_query_to_df(
            """ SELECT * FROM variants WHERE "#CHROM" = 'chr1' AND POS = 28736 AND INFO LIKE '%GQ_stats_min=99.0%' """
        )
        assert len(result) == 1

        result = variants.get_query_to_df(
            """ SELECT * FROM variants WHERE "#CHROM" = 'chr1' AND POS = 28736 AND INFO LIKE '%GQ_stats_max=99.0%' """
        )
        assert len(result) == 1

        result = variants.get_query_to_df(
            """ SELECT * FROM variants WHERE "#CHROM" = 'chr1' AND POS = 28736 AND INFO LIKE '%GQ_stats_mean=99.0%' """
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
        for calculation_name in ["variant_id", "variant_id_varid"]
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
            elif calculation_name == "variant_id_varid":
                variant_id_tag = "varid"
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
    "where_clause, expected_nb_variant",
    [
        (
            None, 7
        ),
        (
            "", 7
        ),
        (
            "\"#CHROM\" = 'chr1'", 6
        ),
        (
            "CLNSIG = 'pathogenic'", 1
        ),
        (
            "SAMPLES.sample1.DP >= 100 AND SAMPLES.sample2.DP >= 100", 2
        )
    ],
)
def test_calculation_variant_filter(
    where_clause, expected_nb_variant
):
    """
    This is a test function for the calculation of variant IDs in a VCF file using the Variants class in Python.
    """

    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # Init files
        input_vcf = tests_data_folder + "/example.vcf.gz"
        output_vcf = f"{tmp_dir}/output.vcf.gz"

        # Construct param dict
        param = {
            "calculation_test": {"calculations": {"variant_filter": {"where_clause": where_clause}}}
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
        assert len(result) == expected_nb_variant, f"Expected {expected_nb_variant} variants with where_clause '{where_clause}', but got {len(result)}"

        # Check if VCF is in correct format with pyVCF
        remove_if_exists([output_vcf])
        variants.export_output()
        try:
            vcf.Reader(filename=output_vcf)
        except:
            assert False