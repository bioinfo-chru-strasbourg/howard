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
from test_needed import tests_folder, tests_config, tests_data_folder


@pytest.mark.parametrize(
    "input_vcf",
    [
        f"{tests_data_folder}/example.ann.transcripts.vcf.gz",
        f"{tests_data_folder}/example.ann.vcf.gz",
        f"{tests_data_folder}/example.dbnsfp.transcripts.vcf.gz",
        f"{tests_data_folder}/example.dbnsfp.no_transcripts.vcf.gz",
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
                            "column_rename": None,
                            "column_clean": True,
                            "column_case": None,
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


@pytest.mark.parametrize(
    "input_vcf",
    [
        f"{tests_data_folder}/example.ann.transcripts.vcf.gz",
        f"{tests_data_folder}/example.ann.vcf.gz",
        f"{tests_data_folder}/example.dbnsfp.transcripts.vcf.gz",
        f"{tests_data_folder}/example.dbnsfp.no_transcripts.vcf.gz",
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
                            "column_rename": None,
                            "column_clean": True,
                            "column_case": None,
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
        f"{tests_data_folder}/example.ann.transcripts.vcf.gz",
        f"{tests_data_folder}/example.ann.vcf.gz",
        f"{tests_data_folder}/example.dbnsfp.transcripts.vcf.gz",
        f"{tests_data_folder}/example.dbnsfp.no_transcripts.vcf.gz",
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
            },
        }
        param_prioritization = {
            "profiles": ["transcripts"],
            "prioritization_config": f"{tests_data_folder}/prioritization_transcripts_profiles.json",
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
        (  # Without transcripts of preference (selected after ordering) and without forcing transcript preference order
            f"{tests_data_folder}/example.ann.transcripts.vcf.gz",
            {
                "profiles": ["transcripts"],
                "prioritization_config": f"{tests_data_folder}/prioritization_transcripts_profiles.json",
                "pzprefix": "PZT",
                "pzfields": ["Score", "Flag"],
                "prioritization_transcripts_order": {
                    "PZTFlag": "ASC",
                    "PZTScore": "DESC",
                },
                "prioritization_score_mode": "HOWARD",
                "prioritization_transcripts": None,
                "prioritization_transcripts_force": False,
                "prioritization_transcripts_version_force": False,
            },
            """
                "#CHROM" = 'chr7'
                AND POS = 55249063
                AND contains(INFO, 'PZTTranscript=NM_005228.5')
            """,
            None,
        ),
        (  # With transcripts of preference selected after ordering
            f"{tests_data_folder}/example.ann.transcripts.vcf.gz",
            {
                "profiles": ["transcripts"],
                "prioritization_config": f"{tests_data_folder}/prioritization_transcripts_profiles.json",
                "pzprefix": "PZT",
                "pzfields": ["Score", "Flag"],
                "prioritization_score_mode": "HOWARD",
                "prioritization_transcripts_order": {
                    "PZTFlag": "ASC",
                    "PZTScore": "DESC",
                },
                "prioritization_transcripts": f"{tests_data_folder}/transcripts.for_prioritization.tsv",
                "prioritization_transcripts_force": False,
                "prioritization_transcripts_version_force": False,
            },
            """
                "#CHROM" = 'chr7'
                AND POS = 55249063
                AND contains(INFO, 'PZTTranscript=NR_047551.1')
            """,
            None,
        ),
        (  # With transcripts of preference forced
            f"{tests_data_folder}/example.ann.transcripts.vcf.gz",
            {
                "profiles": ["transcripts"],
                "prioritization_config": f"{tests_data_folder}/prioritization_transcripts_profiles.json",
                "pzprefix": "PZT",
                "pzfields": ["Score", "Flag"],
                "prioritization_score_mode": "HOWARD",
                "prioritization_transcripts": f"{tests_data_folder}/transcripts.for_prioritization.tsv",
                "prioritization_transcripts_force": True,
                "prioritization_transcripts_version_force": False,
            },
            """
                "#CHROM" = 'chr7'
                AND POS = 55249063
                AND contains(INFO, 'PZTTranscript=NM_001346900.2')
            """,
            None,
        ),
        (  # With transcripts of preference and forcing version
            f"{tests_data_folder}/example.ann.transcripts.vcf.gz",
            {
                "profiles": ["transcripts"],
                "prioritization_config": f"{tests_data_folder}/prioritization_transcripts_profiles.json",
                "pzprefix": "PZT",
                "pzfields": ["Score", "Flag"],
                "prioritization_score_mode": "HOWARD",
                "prioritization_transcripts": f"{tests_data_folder}/transcripts.for_prioritization.tsv",
                "prioritization_transcripts_force": False,
                "prioritization_transcripts_version_force": True,
            },
            """
                "#CHROM" = 'chr7'
                AND POS = 55249063
                AND contains(INFO, 'PZTTranscript=NM_001346941.2')
            """,
            None,
        ),
        (  # No PZfields
            f"{tests_data_folder}/example.ann.transcripts.vcf.gz",
            {
                "profiles": ["transcripts"],
                "prioritization_config": f"{tests_data_folder}/prioritization_transcripts_profiles.json",
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
        (  # Add PZfields
            f"{tests_data_folder}/example.ann.transcripts.vcf.gz",
            {
                "profiles": ["transcripts"],
                "prioritization_config": f"{tests_data_folder}/prioritization_transcripts_profiles.json",
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
        (  # Add PZfields plus
            f"{tests_data_folder}/example.ann.transcripts.vcf.gz",
            {
                "profiles": ["transcripts"],
                "prioritization_config": f"{tests_data_folder}/prioritization_transcripts_profiles.json",
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
        (  # Add PZfields plus with field not present
            f"{tests_data_folder}/example.ann.transcripts.vcf.gz",
            {
                "profiles": ["transcripts"],
                "prioritization_config": f"{tests_data_folder}/prioritization_transcripts_profiles.json",
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
        (  # Order
            f"{tests_data_folder}/example.ann.transcripts.vcf.gz",
            {
                "profiles": ["transcripts"],
                "prioritization_config": f"{tests_data_folder}/prioritization_transcripts_profiles.json",
                "pzprefix": "PZT",
                "prioritization_transcripts_order": {
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
        (  # Order with explode field
            f"{tests_data_folder}/example.ann.transcripts.vcf.gz",
            {
                "profiles": ["transcripts"],
                "prioritization_config": f"{tests_data_folder}/prioritization_transcripts_profiles.json",
                "pzprefix": "PZT",
                "prioritization_transcripts_order": {
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
        (  # Order with explode field not present
            f"{tests_data_folder}/example.ann.transcripts.vcf.gz",
            {
                "profiles": ["transcripts"],
                "prioritization_config": f"{tests_data_folder}/prioritization_transcripts_profiles.json",
                "pzprefix": "PZT",
                "prioritization_transcripts_order": {
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
    The `test_transcripts_prioritization_multiple_param` function tests transcript prioritization
    functionality in a genetic variant analysis pipeline with configurable parameters.

    :param input_vcf: It seems like the `input_vcf` parameter is the path or reference to the VCF
    (Variant Call Format) file that contains genetic variant data. This file is likely used as input for
    the genetic variant analysis pipeline where the transcript prioritization functionality is being
    tested
    :param param_prioritization: The `param_prioritization` parameter is a dictionary that contains
    information about the prioritization configuration for transcripts in a genetic variant analysis
    pipeline. It includes details such as profiles, prioritization configuration file path, prefix, and
    score mode. This parameter is used to customize how transcripts are prioritized during the
    :param where_clause: The `where_clause` parameter in the `test_transcripts_prioritization` function
    is a SQL WHERE clause that is used to filter the results of a query. It specifies a condition that
    must be met for a row to be included in the result set
    :param raise_value_error: The `raise_value_error` parameter in the `test_transcripts_prioritization`
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
                        "column_rename": None,
                        "column_clean": True,
                        "column_case": None,
                    }
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
            },
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


@pytest.mark.parametrize(
    "struct, fields_list",
    [
        (  # By default, no rename, no clean, no case (except clean for snpEff because mandatory)
            {
                "from_column_format": [  # format List, e.g. snpEff
                    {
                        "transcripts_column": "ANN",
                        "transcripts_infos_column": "Feature_ID",
                        "column_clean": True,
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
            [
                "FeatureID",
                "Ensembl_geneid",
                "LIST_S2_score",
                "LIST_S2_pred",
                "VARITY_R_score",
                "Aloft_pred",
            ],
        ),
        (  # No rename, no clean, nor case (except clean for snpEff because mandatory)
            {
                "from_column_format": [  # format List, e.g. snpEff
                    {
                        "transcripts_column": "ANN",
                        "transcripts_infos_column": "Feature_ID",
                        "column_rename": None,
                        "column_clean": True,
                        "column_case": None,
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
            },
            [
                "FeatureID",
                "Ensembl_geneid",
                "LIST_S2_score",
                "LIST_S2_pred",
                "VARITY_R_score",
                "Aloft_pred",
            ],
        ),
        (  # No rename, clean all and nocase (except clean for snpEff because mandatory)
            {
                "from_column_format": [  # format List, e.g. snpEff
                    {
                        "transcripts_column": "ANN",
                        "transcripts_infos_column": "Feature_ID",
                        "column_rename": None,
                        "column_clean": True,
                        "column_case": None,
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
                        "column_rename": None,
                        "column_clean": True,
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
                        "column_clean": True,
                        "column_case": None,
                    },
                ],
            },
            [
                "FeatureID",
                "Ensemblgeneid",
                "LISTS2score",
                "LISTS2pred",
                "VARITYRscore",
                "Aloftpred",
            ],
        ),
        (  # No rename, no clean, and case all (except clean for snpEff because mandatory)
            {
                "from_column_format": [  # format List, e.g. snpEff
                    {
                        "transcripts_column": "ANN",
                        "transcripts_infos_column": "Feature_ID",
                        "column_rename": None,
                        "column_clean": True,
                        "column_case": "lower",
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
                        "column_rename": None,
                        "column_clean": False,
                        "column_case": "lower",
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
                        "column_case": "lower",
                    },
                ],
            },
            [
                "featureid",
                "ensembl_geneid",
                "list_s2_score",
                "list_s2_pred",
                "varity_r_score",
                "aloft_pred",
            ],
        ),
        (  # No rename, clean all and case all (except clean for snpEff because mandatory)
            {
                "from_column_format": [  # format List, e.g. snpEff
                    {
                        "transcripts_column": "ANN",
                        "transcripts_infos_column": "Feature_ID",
                        "column_rename": None,
                        "column_clean": True,
                        "column_case": "lower",
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
                        "column_rename": None,
                        "column_clean": True,
                        "column_case": "lower",
                    },
                    {
                        "transcripts_column": "Ensembl_transcriptid",
                        "transcripts_infos_columns": [
                            "genename",
                            "VARITY_R_score",
                            "Aloft_pred",
                        ],
                        "column_rename": None,
                        "column_clean": True,
                        "column_case": "lower",
                    },
                ],
            },
            [
                "featureid",
                "ensemblgeneid",
                "lists2score",
                "lists2pred",
                "varityrscore",
                "aloftpred",
            ],
        ),
        (  # Rename "genename" columns to merge, transcript ANN id, extra columns on struct_map
            {
                "from_column_format": [  # format List, e.g. snpEff
                    {
                        "transcripts_column": "ANN",
                        "transcripts_infos_column": "Feature_ID",
                        "column_rename": {
                            "Gene_Name": "genename",
                            "Feature_ID": "THETRANSCRIPTOFSNPEFF",
                        },
                        "column_clean": True,
                        "column_case": None,
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
                        "column_clean": False,
                        "column_case": None,
                        "column_rename": {
                            "LIST_S2_score": "LISTScore",
                            "LIST_S2_pred": "LISTPred",
                        },
                    },
                    {
                        "transcripts_column": "Ensembl_transcriptid",
                        "transcripts_infos_columns": [
                            "genename",
                            "VARITY_R_score",
                            "Aloft_pred",
                        ],
                        "column_clean": False,
                        "column_case": None,
                    },
                ],
            },
            [
                "genename",
                "THETRANSCRIPTOFSNPEFF",
                "LISTScore",
                "LISTPred",
                "VARITY_R_score",
                "Aloft_pred",
            ],
        ),
    ],
)
def test_create_transcript_view_rename_clean_case(struct, fields_list):
    """
    The function `test_devel_create_transcript_view` creates a transcript view from a VCF file using
    specified parameters and checks the resulting table for data.

    :param input_vcf: It seems like the `input_vcf` parameter is missing in the provided code snippet.
    Could you please provide the value or path that should be assigned to the `input_vcf` variable in
    the `test_devel_create_transcript_view` function?
    """

    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # Init files
        input_vcf = f"{tests_data_folder}/example.ann.transcripts.vcf.gz"
        output_vcf = f"{tmp_dir}/output.vcf"

        # Construct param dict
        param = {"transcripts": {"table": "transcripts", "struct": struct}}

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
            SELECT column_name
            FROM (
                DESCRIBE SELECT * FROM {transcripts_table}
            )
            WHERE column_name in ('{"', '".join(fields_list)}')
        """
        check = variants.get_query_to_df(query=query_check)

        assert len(check) == len(list(set(fields_list)))


@pytest.mark.parametrize(
    "input_vcf, param_prioritization, where_clause, raise_value_error",
    [
        (  # Add PZfields plus
            f"{tests_data_folder}/example.ann.transcripts.vcf.gz",
            {
                "profiles": ["transcripts"],
                "prioritization_config": f"{tests_data_folder}/prioritization_transcripts_profiles_fields_renamed.json",
                "pzprefix": "PZT",
                "pzfields": ["Score", "Flag", "LISTScore", "LISTPred"],
                "prioritization_score_mode": "HOWARD",
            },
            """
                "#CHROM" = 'chr1'
                AND POS = 69101
                AND contains(INFO, 'PZTTranscript=ENST00000641515')
                AND contains(INFO, 'PZTScore')
                AND contains(INFO, 'PZTFlag')
                AND contains(INFO, 'PZTLISTScore')
                AND contains(INFO, 'PZTLISTPred')
            """,
            None,
        ),
    ],
)
def test_transcripts_prioritization_multiple_param_fields_renamed(
    input_vcf, param_prioritization, where_clause, raise_value_error
):
    """
    The `test_transcripts_prioritization_multiple_param` function tests transcript prioritization
    functionality in a genetic variant analysis pipeline with configurable parameters.

    :param input_vcf: It seems like the `input_vcf` parameter is the path or reference to the VCF
    (Variant Call Format) file that contains genetic variant data. This file is likely used as input for
    the genetic variant analysis pipeline where the transcript prioritization functionality is being
    tested
    :param param_prioritization: The `param_prioritization` parameter is a dictionary that contains
    information about the prioritization configuration for transcripts in a genetic variant analysis
    pipeline. It includes details such as profiles, prioritization configuration file path, prefix, and
    score mode. This parameter is used to customize how transcripts are prioritized during the
    :param where_clause: The `where_clause` parameter in the `test_transcripts_prioritization` function
    is a SQL WHERE clause that is used to filter the results of a query. It specifies a condition that
    must be met for a row to be included in the result set
    :param raise_value_error: The `raise_value_error` parameter in the `test_transcripts_prioritization`
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
                        "column_rename": None,
                        "column_clean": True,
                        "column_case": None,
                    }
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
                        "column_rename": {
                            "LIST_S2_score": "LISTScore",
                            "LIST_S2_pred": "LISTPred",
                        },
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
                        "column_clean": False,
                        "column_case": None,
                    },
                ],
            },
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
