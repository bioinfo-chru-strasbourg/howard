# -*- coding: utf-8 -*-
"""
Tests

Usage:
pytest tests/

Coverage:
coverage run -m pytest tests/test_tools_process.py -x -v --log-cli-level=INFO --capture=tee-sys
coverage report --include=howard/* -m
"""

import argparse
import logging as log
import os
from tempfile import TemporaryDirectory
import vcf  # type: ignore

from howard.functions.commons import remove_if_exists
from howard.objects.variants import Variants
from howard.tools.process import process
from howard.tools.tools import arguments_dict

from test_needed import tests_data_folder, database_files, tests_folder


def test_process():

    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # Init files
        input_vcf = tests_data_folder + "/example.vcf.gz"
        output_vcf = os.path.join(tmp_dir, "output_file.vcf")
        config = {}
        param = "{}"
        annotations = database_files.get("parquet")
        calculations = "VARTYPE"
        prioritizations = "default"
        input_query = None

        # prepare arguments for the query function
        args = argparse.Namespace(
            input=input_vcf,
            output=output_vcf,
            config=config,
            param=param,
            annotations=annotations,
            calculations=calculations,
            prioritizations=prioritizations,
            query=input_query,
            explode_infos=False,
            explode_infos_prefix="",
            explode_infos_fields="*",
            include_header=False,
            arguments_dict=arguments_dict,
        )

        # Remove if output file exists
        remove_if_exists([output_vcf])

        # Query
        process(args)

        # Create object
        variants = Variants(conn=None, input=output_vcf, config=config, load=True)

        # Check annotation
        result = variants.get_query_to_df(
            "SELECT INFO FROM variants WHERE INFO LIKE '%VARTYPE=%'"
        )
        assert len(result) == 7

        # Check annotation
        result = variants.get_query_to_df(
            "SELECT INFO FROM variants WHERE INFO LIKE '%nci60=%'"
        )
        assert len(result) == 1


def test_process_with_param_file():

    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # Init files
        input_vcf = tests_data_folder + "/example.vcf.gz"
        output_vcf = os.path.join(tmp_dir, "output_file.vcf")
        config = {}
        param = tests_folder + "/data/param.snpeff_hgvs.json"
        annotations = database_files.get("parquet")
        calculations = "VARTYPE"
        prioritizations = "default"
        input_query = None

        # prepare arguments for the query function
        args = argparse.Namespace(
            input=input_vcf,
            output=output_vcf,
            config=config,
            param=param,
            annotations=annotations,
            calculations=calculations,
            prioritizations=prioritizations,
            query=input_query,
            explode_infos=False,
            explode_infos_prefix="",
            explode_infos_fields="*",
            include_header=False,
            arguments_dict=arguments_dict,
        )

        # Remove if output file exists
        remove_if_exists([output_vcf])

        # Query
        process(args)

        # Create object
        variants = Variants(conn=None, input=output_vcf, config=config, load=True)

        # Check annotation
        result = variants.get_query_to_df(
            "SELECT INFO FROM variants WHERE INFO LIKE '%VARTYPE=%'"
        )
        assert len(result) == 7

        # Check annotation
        result = variants.get_query_to_df(
            "SELECT INFO FROM variants WHERE INFO LIKE '%nci60=%'"
        )
        assert len(result) == 1


def test_process_with_query():

    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # Init files
        input_vcf = tests_data_folder + "/example.vcf.gz"
        output_vcf = os.path.join(tmp_dir, "output_file.tsv")
        config = {}
        param = "{}"
        annotations = database_files.get("parquet")
        calculations = "VARTYPE"
        prioritizations = "default"
        input_query = "SELECT count(*) AS 'count' FROM variants WHERE INFO LIKE '%VARTYPE%' AND INFO LIKE '%PZScore%'"

        # prepare arguments for the query function
        args = argparse.Namespace(
            input=input_vcf,
            output=output_vcf,
            config=config,
            param=param,
            annotations=annotations,
            calculations=calculations,
            prioritizations=prioritizations,
            query=input_query,
            explode_infos=False,
            explode_infos_prefix="",
            explode_infos_fields="*",
            include_header=False,
            arguments_dict=arguments_dict,
        )

        # Remove if output file exists
        remove_if_exists([output_vcf])

        # Query
        process(args)

        # Create object
        variants = Variants(conn=None, input=output_vcf, config=config, load=True)

        # Check annotation
        result = variants.get_query_to_df("SELECT count FROM variants")
        log.debug(f"result={result}")
        assert result["count"][0] == 7


def test_process_with_chunking():

    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # Init files
        input_vcf = tests_data_folder + "/example.vcf.gz"
        output_vcf = os.path.join(tmp_dir, "output_file.vcf")
        config = {}
        param = tests_folder + "/data/param.snpeff_hgvs.json"
        annotations = database_files.get("parquet")
        calculations = "VARTYPE"
        prioritizations = "default"
        input_query = None

        # prepare arguments for the query function
        args = argparse.Namespace(
            input=input_vcf,
            output=output_vcf,
            config=config,
            param=param,
            annotations=annotations,
            calculations=calculations,
            prioritizations=prioritizations,
            query=input_query,
            explode_infos=False,
            explode_infos_prefix="",
            explode_infos_fields="*",
            include_header=False,
            chunking_enable=True,
            arguments_dict=arguments_dict,
        )

        # Remove if output file exists
        remove_if_exists([output_vcf])

        # Query
        process(args)

        # Create object
        variants = Variants(conn=None, input=output_vcf, config=config, load=True)

        # Check annotation
        result = variants.get_query_to_df(
            "SELECT INFO FROM variants WHERE INFO LIKE '%VARTYPE=%'"
        )
        assert len(result) == 7

        # Check annotation
        result = variants.get_query_to_df(
            "SELECT INFO FROM variants WHERE INFO LIKE '%nci60=%'"
        )
        assert len(result) == 1

        # Check header
        assert "INFO" in variants.get_header_columns()
        assert set(variants.get_header_infos_list()) == set(
            [
                "NS",
                "DP",
                "AA",
                "CLNSIG",
                "SIFT",
                "nci60",
                "VARTYPE",
                "NOMEN",
                "CNOMEN",
                "RNOMEN",
                "NNOMEN",
                "PNOMEN",
                "UPNOMEN",
                "TVNOMEN",
                "TNOMEN",
                "VNOMEN",
                "TPVNOMEN",
                "TPNOMEN",
                "TPVVNOMEN",
                "ENOMEN",
                "GNOMEN",
                "PZScore",
                "PZFlag",
                "PZScore_default",
                "PZFlag_default",
            ]
        )


def test_process_with_chunking_and_chunk_size():

    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # Init files
        input_vcf = tests_data_folder + "/example.vcf.gz"
        output_vcf = os.path.join(tmp_dir, "output_file.vcf")
        config = {}
        param = tests_folder + "/data/param.snpeff_hgvs.json"
        annotations = database_files.get("parquet")
        calculations = "VARTYPE"
        prioritizations = "default"
        input_query = None

        # prepare arguments for the query function
        args = argparse.Namespace(
            input=input_vcf,
            output=output_vcf,
            config=config,
            param=param,
            annotations=annotations,
            calculations=calculations,
            prioritizations=prioritizations,
            query=input_query,
            explode_infos=False,
            explode_infos_prefix="",
            explode_infos_fields="*",
            include_header=False,
            chunking_enable=True,
            chunk_size=4,
            arguments_dict=arguments_dict,
        )

        # Remove if output file exists
        remove_if_exists([output_vcf])

        # Query
        process(args)

        # Create object
        variants = Variants(conn=None, input=output_vcf, config=config, load=True)

        # Check annotation
        result = variants.get_query_to_df(
            "SELECT INFO FROM variants WHERE INFO LIKE '%VARTYPE=%'"
        )
        assert len(result) == 7

        # Check annotation
        result = variants.get_query_to_df(
            "SELECT INFO FROM variants WHERE INFO LIKE '%nci60=%'"
        )
        assert len(result) == 1

        # Check header
        assert "INFO" in variants.get_header_columns()
        assert set(variants.get_header_infos_list()) == set(
            [
                "NS",
                "DP",
                "AA",
                "CLNSIG",
                "SIFT",
                "nci60",
                "VARTYPE",
                "NOMEN",
                "CNOMEN",
                "RNOMEN",
                "NNOMEN",
                "PNOMEN",
                "UPNOMEN",
                "TVNOMEN",
                "TNOMEN",
                "VNOMEN",
                "TPVNOMEN",
                "TPNOMEN",
                "TPVVNOMEN",
                "ENOMEN",
                "GNOMEN",
                "PZScore",
                "PZFlag",
                "PZScore_default",
                "PZFlag_default",
            ]
        )


def test_process_with_chunking_empty_input_file():

    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # Init files
        input_vcf = tests_data_folder + "/example.empty.vcf"
        output_vcf = os.path.join(tmp_dir, "output_file.vcf")
        config = {}
        param = tests_folder + "/data/param.snpeff_hgvs.json"
        annotations = database_files.get("parquet")
        calculations = "VARTYPE"
        prioritizations = "default"
        input_query = None

        # prepare arguments for the query function
        args = argparse.Namespace(
            input=input_vcf,
            output=output_vcf,
            config=config,
            param=param,
            annotations=annotations,
            calculations=calculations,
            prioritizations=prioritizations,
            query=input_query,
            explode_infos=False,
            explode_infos_prefix="",
            explode_infos_fields="*",
            include_header=False,
            chunking_enable=True,
            chunk_size=4,
            arguments_dict=arguments_dict,
        )

        # Remove if output file exists
        remove_if_exists([output_vcf])

        # Query
        process(args)

        # Create object
        variants = Variants(conn=None, input=output_vcf, config=config, load=True)

        # Check annotation
        result = variants.get_query_to_df(
            "SELECT INFO FROM variants WHERE INFO LIKE '%VARTYPE=%'"
        )
        assert len(result) == 0

        # Check header
        assert "INFO" in variants.get_header_columns()
        assert set(variants.get_header_infos_list()) == set(
            [
                "CLNSIG",
                "SIFT",
                "PZScore",
                "DP",
                "PZScore_default",
                "AA",
                "PZFlag_default",
                "PZFlag",
                "NS",
            ]
        )


def test_process_with_chunking_param():
    """
    Test processing with chunking enabled and specific parameters.
    Including exporting transcripts to a TSV file.
    Chunking size and partitioned process on #CHROM.
    """

    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # Init files
        input_vcf = tests_data_folder + "/example.ann.vcf.gz"
        output_vcf = os.path.join(tmp_dir, "output_file.vcf")
        config = {}
        param = tests_folder + "/data/param.transcripts.calculation.chunking.json"
        annotations = database_files.get("parquet")
        calculations = "VARTYPE"
        prioritizations = "default"
        input_query = None
        output_transcripts_tsv = os.path.join(tmp_dir, "output.transcripts.tsv")

        # prepare arguments for the query function
        args = argparse.Namespace(
            input=input_vcf,
            output=output_vcf,
            config=config,
            param=param,
            annotations=annotations,
            calculations=calculations,
            prioritizations=prioritizations,
            query=input_query,
            explode_infos=False,
            explode_infos_prefix="",
            explode_infos_fields="*",
            include_header=False,
            arguments_dict=arguments_dict,
        )

        # Remove if output file exists
        remove_if_exists([output_vcf, output_transcripts_tsv])

        # Query
        process(args)

        # Create object
        variants = Variants(conn=None, input=output_vcf, config=config, load=True)

        # Check annotation
        result = variants.get_query_to_df(
            "SELECT INFO FROM variants WHERE INFO LIKE '%VARTYPE=%'"
        )

        assert len(result) == 7

        # Check annotation
        result = variants.get_query_to_df(
            "SELECT INFO FROM variants WHERE INFO LIKE '%nci60=%'"
        )
        assert len(result) == 1

        if os.path.exists(output_transcripts_tsv):
            # Create object
            transcripts = Variants(
                conn=None, input=output_transcripts_tsv, config=config, load=True
            )
            # Check calculation
            result = transcripts.get_query_to_df("SELECT * FROM variants")
            assert len(result) == 38

        # Check header
        assert "INFO" in variants.get_header_columns()
        assert set(variants.get_header_infos_list()) == set(
            [
                "NS",
                "DP",
                "AA",
                "CLNSIG",
                "SIFT",
                "ANN",
                "LOF",
                "NMD",
                "nci60",
                "VARTYPE",
                "transcripts_json",
                "transcripts_ann",
                "PZTTranscript",
                "PZTScore",
                "PZTFlag",
                "PZTScore_default",
                "PZTFlag_default",
                "PZTScore_transcripts",
                "PZTFlag_transcripts",
                "PZScore",
                "PZFlag",
                "PZScore_default",
                "PZFlag_default",
            ]
        )


def test_process_pipeline():
    """
    Test processing with pipeline enabled and specific parameters.
    Including exporting transcripts to a TSV file.
    """

    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # Init files
        input_vcf = tests_data_folder + "/example.ann.vcf.gz"
        output_vcf = os.path.join(tmp_dir, "output_file.vcf")
        config = {}
        param = tests_folder + "/data/param.pipeline.json"
        annotations = None
        calculations = None
        prioritizations = None
        input_query = None
        output_transcripts_tsv = os.path.join(tmp_dir, "output.transcripts.tsv")

        # prepare arguments for the query function
        args = argparse.Namespace(
            input=input_vcf,
            output=output_vcf,
            config=config,
            param=param,
            annotations=annotations,
            calculations=calculations,
            prioritizations=prioritizations,
            query=input_query,
            explode_infos=False,
            explode_infos_prefix="",
            explode_infos_fields="*",
            include_header=False,
            arguments_dict=arguments_dict,
        )

        # Remove if output file exists
        remove_if_exists([output_vcf, output_transcripts_tsv])

        # Query
        process(args)

        # Create object
        variants = Variants(conn=None, input=output_vcf, config=config, load=True)

        # Check annotation
        result = variants.get_query_to_df(
            "SELECT INFO FROM variants WHERE INFO LIKE '%VARTYPE=%'"
        )
        assert len(result) == 7, f"Expected 7 variants with VARTYPE, but got {len(result)}"

        # Check annotation
        result = variants.get_query_to_df(
            "SELECT INFO FROM variants WHERE INFO LIKE '%snpeff_Annotation=%'"
        )
        assert len(result) == 7, f"Expected 7 variants with snpeff_Annotation, but got {len(result)}"

        # Check annotation
        result = variants.get_query_to_df(
            "SELECT INFO FROM variants WHERE INFO LIKE '%snpeff_hgvs=%'"
        )
        assert len(result) == 7, f"Expected 7 variants with snpeff_hgvs, but got {len(result)}"

        # Check annotation
        result = variants.get_query_to_df(
            "SELECT INFO FROM variants WHERE INFO LIKE '%;hgvs=%'"
        )
        assert len(result) == 7, f"Expected 7 variants with hgvs, but got {len(result)}"

        # Check annotation
        result = variants.get_query_to_df(
            "SELECT INFO FROM variants WHERE INFO LIKE '%;NOMEN=%'"
        )
        assert len(result) == 7, f"Expected 7 variants with NOMEN, but got {len(result)}"

        # Check annotation
        result = variants.get_query_to_df(
            "SELECT INFO FROM variants WHERE INFO LIKE '%PZScore=%'"
        )
        assert len(result) == 7, f"Expected 7 variants with PZScore, but got {len(result)}"

        # Check annotation
        result = variants.get_query_to_df(
            "SELECT INFO FROM variants WHERE INFO LIKE '%SIFT_score=%'"
        )
        assert len(result) == 1, f"Expected 1 variant with SIFT_score, but got {len(result)}"

        # Check annotation
        result = variants.get_query_to_df(
            "SELECT INFO FROM variants WHERE INFO LIKE '%avsnp150=%'"
        )
        assert len(result) == 1, f"Expected 1 variant with avsnp150, but got {len(result)}"

        # Check annotation
        result = variants.get_query_to_df(
            "SELECT INFO FROM variants WHERE INFO LIKE '%cosmic70=%'"
        )
        assert len(result) == 1, f"Expected 1 variant with cosmic70, but got {len(result)}"
         

def test_process_pipeline_full_description():
    """
    Test processing with pipeline enabled and specific parameters.
    Including exporting transcripts to a TSV file.
    """

    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # Init files
        input_vcf = tests_data_folder + "/example.ann.vcf.gz"
        output_vcf = os.path.join(tmp_dir, "output_file.vcf")
        config = {}
        param = tests_folder + "/data/param.pipeline_full_description.json"
        annotations = None
        calculations = None
        prioritizations = None
        input_query = None
        output_transcripts_tsv = os.path.join(tmp_dir, "output.transcripts.tsv")

        # prepare arguments for the query function
        args = argparse.Namespace(
            input=input_vcf,
            output=output_vcf,
            config=config,
            param=param,
            annotations=annotations,
            calculations=calculations,
            prioritizations=prioritizations,
            query=input_query,
            explode_infos=False,
            explode_infos_prefix="",
            explode_infos_fields="*",
            include_header=False,
            arguments_dict=arguments_dict,
        )

        # Remove if output file exists
        remove_if_exists([output_vcf, output_transcripts_tsv])

        # Query
        process(args)

        # Create object
        variants = Variants(conn=None, input=output_vcf, config=config, load=True)

        # Check annotation
        result = variants.get_query_to_df(
            "SELECT INFO FROM variants WHERE INFO LIKE '%VARTYPE=%'"
        )
        assert len(result) == 7, f"Expected 7 variants with VARTYPE, but got {len(result)}"

        # Check annotation
        result = variants.get_query_to_df(
            "SELECT INFO FROM variants WHERE INFO LIKE '%snpeff_Annotation=%'"
        )
        assert len(result) == 7, f"Expected 7 variants with snpeff_Annotation, but got {len(result)}"

        # Check annotation
        result = variants.get_query_to_df(
            "SELECT INFO FROM variants WHERE INFO LIKE '%snpeff_hgvs=%'"
        )
        assert len(result) == 7, f"Expected 7 variants with snpeff_hgvs, but got {len(result)}"

        # Check annotation
        result = variants.get_query_to_df(
            "SELECT INFO FROM variants WHERE INFO LIKE '%;hgvs=%'"
        )
        assert len(result) == 7, f"Expected 7 variants with hgvs, but got {len(result)}"

        # Check annotation
        result = variants.get_query_to_df(
            "SELECT INFO FROM variants WHERE INFO LIKE '%;NOMEN=%'"
        )
        assert len(result) == 7, f"Expected 7 variants with NOMEN, but got {len(result)}"

        # Check annotation
        result = variants.get_query_to_df(
            "SELECT INFO FROM variants WHERE INFO LIKE '%PZScore=%'"
        )
        assert len(result) == 7, f"Expected 7 variants with PZScore, but got {len(result)}"

        # Check annotation
        result = variants.get_query_to_df(
            "SELECT INFO FROM variants WHERE INFO LIKE '%SIFT_score=%'"
        )
        assert len(result) == 1, f"Expected 1 variant with SIFT_score, but got {len(result)}"

        # Check annotation
        result = variants.get_query_to_df(
            "SELECT INFO FROM variants WHERE INFO LIKE '%avsnp150=%'"
        )
        assert len(result) == 1, f"Expected 1 variant with avsnp150, but got {len(result)}"

        # Check annotation
        result = variants.get_query_to_df(
            "SELECT INFO FROM variants WHERE INFO LIKE '%cosmic70=%'"
        )
        assert len(result) == 1, f"Expected 1 variant with cosmic70, but got {len(result)}"


def test_process_pipeline_with_annotation_prioritization_within_calculation():
    """
    Test processing with pipeline enabled and specific parameters.
    Including exporting transcripts to a TSV file.
    """

    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # Init files
        input_vcf = tests_data_folder + "/example.ann.vcf.gz"
        output_vcf = os.path.join(tmp_dir, "output_file.vcf")
        config = {}
        param = tests_folder + "/data/param.pipeline_with_annotation_prioritization.json"
        annotations = None
        calculations = None
        prioritizations = None
        input_query = None
        output_transcripts_tsv = os.path.join(tmp_dir, "output.transcripts.tsv")

        # prepare arguments for the query function
        args = argparse.Namespace(
            input=input_vcf,
            output=output_vcf,
            config=config,
            param=param,
            annotations=annotations,
            calculations=calculations,
            prioritizations=prioritizations,
            query=input_query,
            explode_infos=False,
            explode_infos_prefix="",
            explode_infos_fields="*",
            include_header=False,
            arguments_dict=arguments_dict,
        )

        # Remove if output file exists
        remove_if_exists([output_vcf, output_transcripts_tsv])

        # Query
        process(args)

        # Create object
        variants = Variants(conn=None, input=output_vcf, config=config, load=True)

        # Check if all variants have INFO field
        result = variants.get_query_to_df(
            "SELECT INFO FROM variants"
        )
        assert len(result) == 7, f"Expected 7 variants with INFO, but got {len(result)}"

        # Check annotation
        result = variants.get_query_to_df(
            "SELECT INFO FROM variants WHERE INFO LIKE '%VARTYPE=%'"
        )
        assert len(result) == 7, f"Expected 7 variants with VARTYPE, but got {len(result)}"

        # Check annotation
        result = variants.get_query_to_df(
            "SELECT INFO FROM variants WHERE INFO LIKE '%snpeff_Annotation=%'"
        )
        assert len(result) == 7, f"Expected 7 variants with snpeff_Annotation, but got {len(result)}"

        # Check annotation
        result = variants.get_query_to_df(
            "SELECT INFO FROM variants WHERE INFO LIKE '%snpeff_hgvs=%'"
        )
        assert len(result) == 7, f"Expected 7 variants with snpeff_hgvs, but got {len(result)}"

        # Check annotation
        result = variants.get_query_to_df(
            "SELECT INFO FROM variants WHERE INFO LIKE '%;hgvs=%'"
        )
        assert len(result) == 7, f"Expected 7 variants with hgvs, but got {len(result)}"

        # Check annotation
        result = variants.get_query_to_df(
            "SELECT INFO FROM variants WHERE INFO LIKE '%;NOMEN=%'"
        )
        assert len(result) == 7, f"Expected 7 variants with NOMEN, but got {len(result)}"

        # Check annotation
        result = variants.get_query_to_df(
            "SELECT INFO FROM variants WHERE INFO LIKE '%PZScore=%'"
        )
        assert len(result) == 7, f"Expected 7 variants with PZScore, but got {len(result)}"

        # Check annotation
        result = variants.get_query_to_df(
            "SELECT INFO FROM variants WHERE INFO LIKE '%SIFT_score=%'"
        )
        assert len(result) == 1, f"Expected 1 variant with SIFT_score, but got {len(result)}"

        # Check annotation
        result = variants.get_query_to_df(
            "SELECT INFO FROM variants WHERE INFO LIKE '%avsnp150=%'"
        )
        assert len(result) == 1, f"Expected 1 variant with avsnp150, but got {len(result)}"

        # Check annotation
        result = variants.get_query_to_df(
            "SELECT INFO FROM variants WHERE INFO LIKE '%cosmic70=%'"
        )
        assert len(result) == 1, f"Expected 1 variant with cosmic70, but got {len(result)}"
