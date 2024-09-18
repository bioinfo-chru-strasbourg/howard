# -*- coding: utf-8 -*-
"""
Tests

Usage:
pytest tests/

Coverage:
coverage run -m pytest tests/test_variants_transcripts.py -x -vv --log-cli-level=DEBUG --capture=tee-sys
coverage report --include=howard/* -m
"""

import logging as log
from tempfile import TemporaryDirectory
import pytest
import vcf

from howard.functions.commons import remove_if_exists
from howard.objects.variants import Variants
from test_needed import tests_folder, tests_config


@pytest.mark.parametrize(
    "input_vcf",
    [
        "tests/data/example.ann.transcripts.vcf.gz",
        "tests/data/example.ann.vcf.gz",
        "tests/data/example.dbnsfp.transcripts.vcf.gz",
        "tests/data/example.dbnsfp.no_transcripts.vcf.gz",
    ],
)
def test_create_transcript_view(input_vcf):
    """
    The function `test_devel_create_transcript_view` creates a transcript view from a VCF file using
    specified parameters and checks the resulting table for data.

    :param input_vcf: It seems like the `input_vcf` parameter is missing in the provided code snippet.
    Could you please provide the value or path that should be assigned to the `input_vcf` variable in
    the `test_devel_create_transcript_view` function?
    """

    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # Init files
        output_vcf = f"{tmp_dir}/output.vcf"

        # Construct param dict
        param = {
            "transcripts": {
                "table": "transcripts",
                "struct": {
                    "from_column_format": [  # format List, e.g. snpEff
                        {
                            "transcripts_column": "ANN",
                            "transcripts_infos_column": "Feature_ID",
                        }
                    ],
                    "from_columns_map": [  # format List, e.g. dbNSFP columns
                        {
                            "transcripts_column": "Ensembl_transcriptid",
                            "transcripts_infos_columns": [
                                "genename",
                                "Ensembl_geneid",
                                "LIST_S2_score",
                                "LIST_S2_pred",
                            ],
                        },
                        {
                            "transcripts_column": "Ensembl_transcriptid",
                            "transcripts_infos_columns": [
                                "genename",
                                "VARITY_R_score",
                                "Aloft_pred",
                            ],
                        },
                    ],
                },
            }
        }

        # Create object
        variants = Variants(
            conn=None, input=input_vcf, output=output_vcf, param=param, load=True
        )

        # Create transcript view
        transcripts_table = variants.create_transcript_view()

        # Check table exists
        assert transcripts_table is not None

        # Check table content
        query_check = f"""
            SELECT * FROM {transcripts_table}
            ORDER BY "#CHROM", POS, REF, ALT, transcript
        """
        check = variants.get_query_to_df(query=query_check)
        log.debug(f"check={check}")
        assert len(check) > 0


@pytest.mark.parametrize(
    "input_vcf",
    [
        "tests/data/example.ann.transcripts.vcf.gz",
        "tests/data/example.ann.vcf.gz",
        "tests/data/example.dbnsfp.transcripts.vcf.gz",
        "tests/data/example.dbnsfp.no_transcripts.vcf.gz",
    ],
)
def test_create_transcript_view_to_variants(input_vcf):
    """ """

    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # Init files
        output_vcf = f"{tmp_dir}/output.vcf"

        # Construct param dict
        param = {
            "transcripts": {
                "table": "transcripts",
                # "transcripts_info_field_json": "transcripts_json",
                # "transcripts_info_field_format": "transcripts_ann",
                "struct": {
                    "from_column_format": [  # format List, e.g. snpEff
                        {
                            "transcripts_column": "ANN",
                            "transcripts_infos_column": "Feature_ID",
                        }
                    ],
                    "from_columns_map": [  # format List, e.g. dbNSFP columns
                        {
                            "transcripts_column": "Ensembl_transcriptid",
                            "transcripts_infos_columns": [
                                "genename",
                                "Ensembl_geneid",
                                "LIST_S2_score",
                                "LIST_S2_pred",
                            ],
                        },
                        {
                            "transcripts_column": "Ensembl_transcriptid",
                            "transcripts_infos_columns": [
                                "genename",
                                "VARITY_R_score",
                                "Aloft_pred",
                            ],
                        },
                    ],
                },
            }
        }

        # Create object
        variants = Variants(
            conn=None, input=input_vcf, output=output_vcf, param=param, load=True
        )

        # Create transcript view
        transcripts_table = variants.create_transcript_view()

        # Check table exists
        assert transcripts_table is not None

        # Check table content
        query_check = f"""
            SELECT * FROM {transcripts_table}
            ORDER BY "#CHROM", POS, REF, ALT, transcript
        """
        check = variants.get_query_to_df(query=query_check)
        assert len(check) > 0

        # Param in function param - generate column json
        assert variants.transcript_view_to_variants(
            transcripts_info_json="transcripts_in_json_as_column",
            transcripts_info_field_json=None,
            transcripts_info_format=None,
            transcripts_info_field_format=None,
        )
        query_check = """
            SELECT transcripts_in_json_as_column FROM variants
        """
        try:
            check = variants.get_query_to_df(query=query_check)
            assert len(check) > 0
        except:
            assert False

        # Param in function param - generate INFO field json
        assert variants.transcript_view_to_variants(
            transcripts_info_json=None,
            transcripts_info_field_json="transcripts_in_json",
            transcripts_info_format=None,
            transcripts_info_field_format=None,
        )
        query_check = """
            SELECT * FROM variants
            WHERE contains(INFO, 'transcripts_in_json=')
        """
        check = variants.get_query_to_df(query=query_check)
        assert len(check) > 0

        # Param in function param - generate column format
        assert variants.transcript_view_to_variants(
            transcripts_info_json=None,
            transcripts_info_field_json=None,
            transcripts_info_format="transcripts_in_format_as_column",
            transcripts_info_field_format=None,
        )
        query_check = """
            SELECT transcripts_in_format_as_column FROM variants
        """
        try:
            check = variants.get_query_to_df(query=query_check)
            assert len(check) > 0
        except:
            assert False

        # Param in function param - generate INFO field format
        assert variants.transcript_view_to_variants(
            transcripts_info_json=None,
            transcripts_info_field_json=None,
            transcripts_info_format=None,
            transcripts_info_field_format="transcripts_in_format",
        )
        query_check = """
            SELECT * FROM variants
            WHERE contains(INFO, 'transcripts_in_format=')
        """
        check = variants.get_query_to_df(query=query_check)
        assert len(check) > 0

        # Param not available
        assert not variants.transcript_view_to_variants()

        # Param in param dict for JSON
        param["transcripts"][
            "transcripts_info_field_json"
        ] = "transcripts_json_in_param_dict"
        assert variants.transcript_view_to_variants(param=param)
        query_check = """
            SELECT * FROM variants
            WHERE contains(INFO, 'transcripts_json_in_param_dict=')
        """
        check = variants.get_query_to_df(query=query_check)
        assert len(check) > 0

        # Param in param dict for JSON with another field name
        param["transcripts"][
            "transcripts_info_field_json"
        ] = "transcripts_json_in_param_dict_duplicate"
        variants.set_param(param=param)
        assert variants.transcript_view_to_variants()
        query_check = """
            SELECT * FROM variants
            WHERE contains(INFO, 'transcripts_json_in_param_dict_duplicate=')
        """
        check = variants.get_query_to_df(query=query_check)
        assert len(check) > 0

        # Param in param dict for JSON with another field name
        param["transcripts"][
            "transcripts_info_field_format"
        ] = "transcripts_format_in_param_dict"
        variants.set_param(param=param)
        assert variants.transcript_view_to_variants()
        query_check = """
            SELECT * FROM variants
            WHERE contains(INFO, 'transcripts_format_in_param_dict=')
        """
        check = variants.get_query_to_df(query=query_check)
        assert len(check) > 0

        # Check if VCF is in correct format with pyVCF
        remove_if_exists([output_vcf])
        variants.export_output()
        try:
            vcf.Reader(filename=output_vcf)
        except:
            assert False


@pytest.mark.parametrize(
    "input_vcf",
    [
        "tests/data/example.ann.transcripts.vcf.gz",
        "tests/data/example.ann.vcf.gz",
        "tests/data/example.dbnsfp.transcripts.vcf.gz",
        "tests/data/example.dbnsfp.no_transcripts.vcf.gz",
    ],
)
def test_transcripts_prioritization(input_vcf):
    """ """

    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # Init files
        output_vcf = f"{tmp_dir}/output.vcf"

        # Construct param dict
        param = {}
        param_struct = {
            "table": "transcripts",
            "column_id": "transcript",
            "transcripts_info_json": "transcripts_json",
            "transcripts_info_field": "transcripts_json",
            "struct": {
                "from_column_format": [
                    {
                        "transcripts_column": "ANN",
                        "transcripts_infos_column": "Feature_ID",
                    },
                    {
                        "transcripts_column": "ANN",
                        "transcripts_infos_column": "Feature_ID",
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
                    },
                    {
                        "transcripts_column": "Ensembl_transcriptid",
                        "transcripts_infos_columns": [
                            "genename",
                            "VARITY_R_score",
                            "Aloft_pred",
                        ],
                    },
                ],
            },
        }
        param_prioritization = {
            "profiles": ["transcripts"],
            "prioritization_config": "config/prioritization_transcripts_profiles.json",
            "pzprefix": "PZT",
            "prioritization_score_mode": "HOWARD",
        }

        # Param without prioritization
        param_without_prioritization = {"transcripts": dict(param_struct)}

        # Param with prioritization
        param_with_prioritization = {"transcripts": dict(param_struct)}
        param_with_prioritization["transcripts"]["prioritization"] = dict(
            param_prioritization
        )

        # Create object
        variants = Variants(
            conn=None, input=input_vcf, output=output_vcf, param=param, load=True
        )

        # Create transcript view
        transcripts_table = variants.create_transcript_view(
            param=param_without_prioritization
        )

        # Check table exists
        assert transcripts_table is not None

        # Transcripts without prioritization
        assert not variants.transcripts_prioritization(
            param=param_without_prioritization
        )

        # Transcripts with prioritization
        assert variants.transcripts_prioritization(param=param_with_prioritization)

        # Check transcript prioritization result
        # Check table content
        query_check = """
            SELECT * FROM variants
            WHERE "#CHROM" = 'chr1'
              AND POS = 69101
              AND contains(INFO, 'PZTTranscript')
        """
        check = variants.get_query_to_df(query=query_check)
        assert len(check) > 0

        # Check if VCF is in correct format with pyVCF
        remove_if_exists([output_vcf])
        variants.export_output(output_file=output_vcf)
        try:
            vcf.Reader(filename=output_vcf)
        except:
            assert False


@pytest.mark.parametrize(
    "input_vcf, param_prioritization, where_clause, raise_value_error",
    [
        (
            # No PZfields
            "tests/data/example.ann.transcripts.vcf.gz",
            {
                "profiles": ["transcripts"],
                "prioritization_config": "config/prioritization_transcripts_profiles.json",
                "pzprefix": "PZT",
                "prioritization_score_mode": "HOWARD",
            },
            """
                "#CHROM" = 'chr1'
                AND POS = 69101
                AND contains(INFO, 'PZTTranscript=ENST00000641515')
                AND NOT contains(INFO, 'PZTScore')
            """,
            None,
        ),
        (
            # Add PZfields
            "tests/data/example.ann.transcripts.vcf.gz",
            {
                "profiles": ["transcripts"],
                "prioritization_config": "config/prioritization_transcripts_profiles.json",
                "pzprefix": "PZT",
                "pzfields": ["Score", "Flag"],
                "prioritization_score_mode": "HOWARD",
            },
            """
                "#CHROM" = 'chr1'
                AND POS = 69101
                AND contains(INFO, 'PZTTranscript=ENST00000641515')
                AND contains(INFO, 'PZTScore')
                AND contains(INFO, 'PZTFlag')
            """,
            None,
        ),
        (
            # Add PZfields plus
            "tests/data/example.ann.transcripts.vcf.gz",
            {
                "profiles": ["transcripts"],
                "prioritization_config": "config/prioritization_transcripts_profiles.json",
                "pzprefix": "PZT",
                "pzfields": ["Score", "Flag", "LIST_S2_score", "LIST_S2_pred"],
                "prioritization_score_mode": "HOWARD",
            },
            """
                "#CHROM" = 'chr1'
                AND POS = 69101
                AND contains(INFO, 'PZTTranscript=ENST00000641515')
                AND contains(INFO, 'PZTScore')
                AND contains(INFO, 'PZTFlag')
                AND contains(INFO, 'LIST_S2_score')
                AND contains(INFO, 'LIST_S2_pred')
            """,
            None,
        ),
        (
            # Add PZfields plus with field not present
            "tests/data/example.ann.transcripts.vcf.gz",
            {
                "profiles": ["transcripts"],
                "prioritization_config": "config/prioritization_transcripts_profiles.json",
                "pzprefix": "PZT",
                "pzfields": [
                    "Score",
                    "Flag",
                    "LIST_S2_score",
                    "field_not_present_in_header",
                ],
                "prioritization_score_mode": "HOWARD",
            },
            None,
            "INFO/field_not_present_in_header NOT IN header",
        ),
        (
            # Order
            "tests/data/example.ann.transcripts.vcf.gz",
            {
                "profiles": ["transcripts"],
                "prioritization_config": "config/prioritization_transcripts_profiles.json",
                "pzprefix": "PZT",
                "pzorder": {
                    "LIST_S2_score": "ASC",
                    "transcript": "ASC",
                },
                "prioritization_score_mode": "HOWARD",
            },
            """
                "#CHROM" = 'chr1'
                AND POS = 69101
                AND contains(INFO, 'PZTTranscript=ENST00000335137')
            """,
            None,
        ),
        (
            # Order with explode field
            "tests/data/example.ann.transcripts.vcf.gz",
            {
                "profiles": ["transcripts"],
                "prioritization_config": "config/prioritization_transcripts_profiles.json",
                "pzprefix": "PZT",
                "pzorder": {
                    "CADD_raw": "ASC",
                    "transcript": "ASC",
                },
                "prioritization_score_mode": "HOWARD",
            },
            """
                "#CHROM" = 'chr1'
                AND POS = 69101
                AND contains(INFO, 'PZTTranscript=ENST00000335137')
            """,
            None,
        ),
        (
            # Order with explode field not present
            "tests/data/example.ann.transcripts.vcf.gz",
            {
                "profiles": ["transcripts"],
                "prioritization_config": "config/prioritization_transcripts_profiles.json",
                "pzprefix": "PZT",
                "pzorder": {
                    "field_not_present_in_header": "ASC",
                    "transcript": "ASC",
                },
                "prioritization_score_mode": "HOWARD",
            },
            None,
            "INFO/field_not_present_in_header NOT IN header",
        ),
    ],
)
def test_transcripts_prioritization_multiple_param(
    input_vcf, param_prioritization, where_clause, raise_value_error
):
    """
    The function `test_transcripts_prioritization_multiple_param` tests transcript prioritization functionality in a
    genetic variant analysis pipeline.

    :param input_vcf: It seems like you were about to provide some information about the `input_vcf`
    parameter in the `test_transcripts_prioritization_multiple_param` function. Could you please provide more details
    or let me know how I can assist you further with this function?
    :param param_prioritization: param_prioritization is a dictionary containing information about
    prioritization configuration for transcripts. It includes details such as profiles, prioritization
    configuration file path, prefix, and score mode
    :param where_clause: The `where_clause` parameter is a SQL WHERE clause that is used to filter the
    results of a query. It specifies a condition that must be met for a row to be included in the result
    set. For example, if you have a WHERE clause like `genename = 'BRCA1'
    :param raise_value_error: The `raise_value_error` parameter in the `test_transcripts_prioritization_multiple_param`
    function is a boolean flag that determines whether the test should raise a `ValueError` and check if
    the raised error message matches a specific value. If `raise_value_error` is `True`, the test will
    raise
    """

    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # Init files
        output_vcf = f"{tmp_dir}/output.vcf"

        # Construct param dict
        param = {}
        param_struct = {
            "table": "transcripts",
            "column_id": "transcript",
            "transcripts_info_json": "transcripts_json",
            "transcripts_info_field": "transcripts_json",
            "struct": {
                "from_column_format": [
                    {
                        "transcripts_column": "ANN",
                        "transcripts_infos_column": "Feature_ID",
                    },
                    {
                        "transcripts_column": "ANN",
                        "transcripts_infos_column": "Feature_ID",
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
                    },
                    {
                        "transcripts_column": "Ensembl_transcriptid",
                        "transcripts_infos_columns": [
                            "genename",
                            "VARITY_R_score",
                            "Aloft_pred",
                        ],
                    },
                ],
            },
        }

        # No pzfields
        ##############

        # Param without prioritization
        param_without_prioritization = {"transcripts": dict(param_struct)}

        # Param with prioritization
        param_with_prioritization = {"transcripts": dict(param_struct)}
        param_with_prioritization["transcripts"]["prioritization"] = dict(
            param_prioritization
        )

        # Create object
        variants = Variants(
            conn=None, input=input_vcf, output=output_vcf, param=param, load=True
        )

        # Create transcript view
        transcripts_table = variants.create_transcript_view(
            param=param_without_prioritization
        )

        # Check table exists
        assert transcripts_table is not None

        # If Raise with Value Error
        if raise_value_error:

            # Catch ValueError
            with pytest.raises(ValueError) as excinfo:

                # Prioritization
                variants.transcripts_prioritization(param=param_with_prioritization)

            assert str(excinfo.value) == raise_value_error

        # If expected results
        if where_clause:

            # Prioritization
            assert variants.transcripts_prioritization(param=param_with_prioritization)

            # Check transcript prioritization result
            # Check table content
            query_check = f"""
                SELECT * FROM variants
                WHERE {where_clause}
            """
            check = variants.get_query_to_df(query=query_check)
            assert len(check) > 0

            # Export
            ########

            # Check if VCF is in correct format with pyVCF
            remove_if_exists([output_vcf])
            variants.export_output(output_file=output_vcf)
            try:
                vcf.Reader(filename=output_vcf)
            except:
                assert False
