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


def test_genotype_format():
    """
    The function `test_genotype_format` tests if a VCF file is in the correct format with specified
    sample names using pyVCF library.
    """

    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # Init files
        input_vcf = os.path.join(
            tests_data_folder, "example.with_allowed_genotypes.vcf"
        )
        output_vcf = f"{tmp_dir}/output.vcf"
        expected_samples = [
            "sample1",
            "sample2",
            "sample3",
            "sample4",
        ]

        # Construct param dict
        param = {}

        # Create object
        variants = Variants(
            conn=None, input=input_vcf, output=output_vcf, param=param, load=True
        )

        # Check if VCF is in correct format with pyVCF
        remove_if_exists([output_vcf])
        variants.export_output(output_file=output_vcf)
        try:
            vcf_obj = vcf.Reader(filename=output_vcf)
            assert vcf_obj.samples == expected_samples
        except:
            assert False


@pytest.mark.parametrize(
    "check, samples, samples_force, samples_list_expected",
    [
        # No list of samples
        (
            True,
            None,
            False,
            [
                "sample1",
                "sample2",
                "sample3",
                "sample4",
            ],
        ),
        (
            True,
            None,
            True,
            [
                "sample1",
                "sample2",
                "sample3",
                "sample4",
            ],
        ),
        (
            False,
            None,
            False,
            [
                "sample1",
                "sample2",
                "sample3",
                "sample4",
                "ann1",
                "ann2",
            ],
        ),
        (
            False,
            None,
            True,
            [
                "sample1",
                "sample2",
                "sample3",
                "sample4",
                "ann1",
                "ann2",
            ],
        ),
        # Add a list of samples
        (
            True,
            ["sample1", "sample2"],
            False,
            [
                "sample1",
                "sample2",
            ],
        ),
        (
            True,
            ["sample1", "sample2"],
            True,
            [
                "sample1",
                "sample2",
            ],
        ),
        (
            False,
            ["sample1", "sample2"],
            False,
            [
                "sample1",
                "sample2",
            ],
        ),
        (
            False,
            ["sample1", "sample2"],
            True,
            [
                "sample1",
                "sample2",
            ],
        ),
        # List of samples with not sample 'ann1'
        (
            True,
            ["sample1", "sample2", "ann1"],
            False,
            [
                "sample1",
                "sample2",
            ],
        ),
        (
            True,
            ["sample1", "sample2", "ann1"],
            True,
            ["sample1", "sample2", "ann1"],
        ),
        (
            False,
            ["sample1", "sample2", "ann1"],
            False,
            ["sample1", "sample2", "ann1"],
        ),
        (
            False,
            ["sample1", "sample2", "ann1"],
            True,
            ["sample1", "sample2", "ann1"],
        ),
    ],
)
def test_get_header_sample_list(check, samples, samples_force, samples_list_expected):
    """
    The function `test_get_header_sample_list` tests the `get_header_sample_list` method of the
    `Variants` class.

    :param check: The `check` parameter is likely a boolean flag that determines whether a certain
    condition should be checked or not within the `get_header_sample_list` method of the `Variants`
    class. It is used to control the behavior of the method based on the value passed to it
    :param samples: The `samples` parameter in the `test_get_header_sample_list` function likely refers
    to a list of sample names or identifiers that are being passed as an argument to the
    `get_header_sample_list` method of the `Variants` class. This list is used to specify which samples
    should be included
    :param samples_force: The `samples_force` parameter in the `test_get_header_sample_list` function is
    likely a list of sample names that you want to force inclusion in the header sample list. This
    parameter is used in the `get_header_sample_list` method of the `Variants` class to ensure that
    specific samples
    :param samples_list_expected: The `samples_list_expected` parameter in the
    `test_get_header_sample_list` function is the expected list of samples that should be returned by
    the `get_header_sample_list` method of the `Variants` class. This list is compared with the actual
    result obtained from calling the method with the provided
    """

    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # Init files
        input_vcf = os.path.join(
            tests_data_folder, "example.with_annotation_columns.tsv"
        )
        output_vcf = f"{tmp_dir}/output.vcf"

        # Construct param dict
        param = {}

        # Create object
        variants = Variants(
            conn=None, input=input_vcf, output=output_vcf, param=param, load=True
        )

        # samples list
        sample_list_result = variants.get_header_sample_list(
            check=check,
            samples=samples,
            samples_force=samples_force,
        )
        assert sample_list_result == samples_list_expected


@pytest.mark.parametrize(
    "input_vcf, param_samples",
    [
        (
            os.path.join(tests_data_folder, input_vcf),
            param_samples,
        )
        for input_vcf in ["example.with_annotation_columns.tsv"]
        for param_samples in [
            {
                "samples_list_expected": [
                    "sample1",
                    "sample2",
                    "sample3",
                    "sample4",
                ]
            },
            {
                "check": False,
                "samples_list_expected": [
                    "sample1",
                    "sample2",
                    "sample3",
                    "sample4",
                    "ann1",
                    "ann2",
                ],
            },
            {
                "check": True,
                "samples_list_expected": [
                    "sample1",
                    "sample2",
                    "sample3",
                    "sample4",
                ],
            },
            {
                "list": ["sample1", "sample2", "ann1"],
                "check": True,
                "samples_list_expected": ["sample1", "sample2", "ann1"],
            },
            {
                "list": ["sample1", "sample2", "ann1"],
                "check": False,
                "samples_list_expected": ["sample1", "sample2", "ann1"],
            },
            {
                "list": None,
                "check": False,
                "samples_list_expected": [
                    "sample1",
                    "sample2",
                    "sample3",
                    "sample4",
                    "ann1",
                    "ann2",
                ],
            },
        ]
    ],
)
def test_export_samples(input_vcf, param_samples):
    """
    The function `test_export_samples` exports samples from a VCF file and checks if the output VCF is
    in the correct format using pyVCF.

    :param input_vcf: The `input_vcf` parameter in the `test_export_samples` function is likely a file
    path to a VCF (Variant Call Format) file. This file contains genetic variant information such as
    SNPs, insertions, deletions, and structural variations for a set of samples
    :param param_samples: It seems like you haven't provided the value for the `param_samples`
    parameter. Could you please provide the specific value or details for `param_samples` that you would
    like to use in the `test_export_samples` function?
    """

    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # Init files
        output_vcf = f"{tmp_dir}/output.vcf"

        # Construct param dict
        param = {"samples": param_samples}

        # Create object
        variants = Variants(
            conn=None, input=input_vcf, output=output_vcf, param=param, load=True
        )

        samples_list_expected = param.get("samples", {}).get(
            "samples_list_expected", None
        )

        # Check if VCF is in correct format with pyVCF
        remove_if_exists([output_vcf])
        variants.export_output(output_file=output_vcf)
        try:
            vcf_obj = vcf.Reader(filename=output_vcf)
            assert vcf_obj.samples == samples_list_expected
        except:
            assert False


@pytest.mark.parametrize(
    "input_vcf, remove_info, add_samples, where_clause",
    [
        (
            os.path.join(tests_data_folder, input_vcf),
            remove_info,
            add_samples,
            where_clause,
        )
        for input_vcf in ["example.vcf.gz", "example.without_sample.vcf"]
        for remove_info in [True, False]
        for add_samples in [True, False]
        for where_clause in [
            None,
            "",
            """WHERE "#CHROM" == 'chr1' """,
            """WHERE regexp_matches("INFO",'DP=') """,
        ]
    ],
)
def test_export_variant_vcf(input_vcf, remove_info, add_samples, where_clause):

    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # Init files
        output_vcf = f"{tmp_dir}/output.vcf.gz"

        # Config
        config = tests_config.copy()

        # Create object
        variants = Variants(
            conn=None,
            input=input_vcf,
            output=output_vcf,
            config=config,
            load=True,
        )
        try:
            # Export varaint VCF
            variants.export_variant_vcf(
                vcf_file=output_vcf,
                remove_info=remove_info,
                add_samples=add_samples,
                list_samples=[],
                where_clause=where_clause,
                index=False,
                threads=1,
            )

            # Check number of variants
            variants_output = Variants(
                conn=None,
                input=output_vcf,
                config=config,
                load=True,
            )

            # Input filtered
            if not where_clause:
                where_clause = ""
            query = f"SELECT * FROM variants {where_clause}"
            df_input_filtered = variants.get_query_to_df(query=query)

            # Output
            df_output = variants_output.get_query_to_df(query="SELECT * FROM variants")

            assert len(df_input_filtered) == len(df_output)

            assert True
        except:
            assert False


@pytest.mark.parametrize(
    "explode_infos_fields, explode_infos_fields_expected, remove_fields_not_in_header",
    [
        ("SIFT", ["SIFT"], False),
        ("SIFT,DP,AD", ["SIFT", "DP", "AD"], False),
        ("*", sorted(["NS", "DP", "AA", "CLNSIG", "SIFT"]), False),
        ("DP,*", ["DP", "AA", "CLNSIG", "NS", "SIFT"], False),
        ("*,DP", ["AA", "CLNSIG", "NS", "SIFT", "DP"], False),
        ("DP,DP", ["DP"], False),
        ("DP,SIFT,DP,AA", ["DP", "SIFT", "AA"], False),
        ("DP,*,DP,AA", ["DP", "CLNSIG", "NS", "SIFT", "AA"], False),
        ("NOT_field", ["NOT_field"], False),
        ("NOT_field,DP", ["NOT_field", "DP"], False),
        ("NOT_field", [], True),
        ("NOT_field,DP", ["DP"], True),
        ("NOT_field,*", ["AA", "CLNSIG", "DP", "NS", "SIFT"], True),
        ("AA, DP ,SIFT,NS ", ["AA", "DP", "SIFT", "NS"], False),
        (".*S.*", ["CLNSIG", "NS", "SIFT"], False),
        ("NS,.*S.*", ["NS", "CLNSIG", "SIFT"], False),
        (".*S.*,NS", ["CLNSIG", "SIFT", "NS"], False),
        ("NS,.*S.*,*", ["NS", "CLNSIG", "SIFT", "AA", "DP"], False),
        (".*S.*,NS,*", ["CLNSIG", "SIFT", "NS", "AA", "DP"], False),
        (".*S.*,*,NS", ["CLNSIG", "SIFT", "AA", "DP", "NS"], False),
        ("*,.*S.*,NS", ["AA", "CLNSIG", "DP", "SIFT", "NS"], False),
    ],
)
def test_get_explode_infos_fields(
    explode_infos_fields, explode_infos_fields_expected, remove_fields_not_in_header
):
    """
    The function `test_get_explode_infos_fields()` tests the `get_explode_infos_fields()` method of the
    `Variants` class.
    """

    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # Init files
        input_vcf = tests_data_folder + "/example.vcf.gz"
        output_tsv = f"{tmp_dir}/example.tsv"

        # remove if exists
        remove_if_exists([output_tsv])

        # Create object
        variants = Variants(input=input_vcf, output=output_tsv, load=True)

        # 1 field
        explode_infos_fields_list = variants.get_explode_infos_fields(
            explode_infos_fields,
            remove_fields_not_in_header=remove_fields_not_in_header,
        )
        assert explode_infos_fields_list == explode_infos_fields_expected


def test_export_query():
    """
    This is a test function for exporting data from a VCF file to a TSV file using SQL queries.
    """

    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # Init files
        input_vcf = tests_data_folder + "/example.vcf.gz"
        output_tsv = f"{tmp_dir}/example.tsv"

        # remove if exists
        remove_if_exists([output_tsv])

        # Create object
        variants = Variants(input=input_vcf, output=output_tsv, load=True)

        # Check get_output
        query = 'SELECT "#CHROM", POS, REF, ALT, INFO FROM variants'
        variants.export_output(query=query)
        assert os.path.exists(output_tsv)

        # Check get_output without header
        output_header = output_tsv + ".hdr"
        remove_if_exists([output_tsv, output_tsv + ".hdr"])
        variants.export_output(output_header=output_header, query=query)
        assert os.path.exists(output_tsv) and os.path.exists(output_header)


@pytest.mark.parametrize(
    "database_input_index, database_output_format",
    [
        (database_input_index, database_output_format)
        for database_input_index in [
            "parquet",
            "partition_parquet",
            "vcf",
            "vcf_gz",
            "tsv",
            "csv",
            "example_vcf",
        ]
        for database_output_format in [
            "parquet",
            "partition_parquet",
            "vcf",
            "vcf.gz",
            "tsv",
            "csv",
            "json",
            "bed",
        ]
    ],
)
def test_export_output(database_input_index, database_output_format):
    """
    The function tests the export functionality of a database for various input and output formats.
    """

    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        parquet_partitions = None
        # specific partition_parquet
        if database_output_format in ["partition_parquet"]:
            database_output_format = "parquet"
            parquet_partitions = ["#CHROM"]
        input_database = database_files.get(database_input_index)
        output_database = f"{tmp_dir}/output_database.{database_output_format}"
        output_header = output_database + ".hdr"
        variants = Variants(
            input=input_database,
            output=output_database,
            config=tests_config,
            load=True,
        )
        remove_if_exists([output_database, output_header])
        try:
            assert variants.export_output(
                output_file=output_database,
                output_header=output_header,
                parquet_partitions=parquet_partitions,
            )
            if database_output_format == "vcf":
                try:
                    vcf.Reader(filename=output_database)
                except:
                    assert False
        except:
            assert False


def test_set_get_input():
    """
    This function tests the set_input and get_input methods of the Variants class in Python.
    """

    # Init files
    input_vcf = tests_data_folder + "/example.vcf.gz"
    new_input_vcf = tests_data_folder + "/example.parquet"

    # Create connection
    conn = duckdb.connect(":memory:")

    # Create object
    variants = Variants(conn=conn, input=input_vcf)

    # Check get_input
    assert input_vcf == variants.get_input()

    # Check new input vcf
    variants.set_input(new_input_vcf)
    assert new_input_vcf == variants.get_input()


def test_set_get_output():
    """
    This function tests the set_output and get_output methods of the Variants class in Python.
    """

    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # Init files
        input_vcf = tests_data_folder + "/example.vcf.gz"
        output_vcf = tests_data_folder + f"{tmp_dir}/example.parquet"
        new_output_vcf = tests_data_folder + f"{tmp_dir}/example.parquet"

        # Create connection
        conn = duckdb.connect(":memory:")

        # Create object
        variants = Variants(conn=conn, input=input_vcf, output=output_vcf)

        # Check get_output
        assert output_vcf == variants.get_output()

        # Check new output vcf
        variants.set_output(new_output_vcf)
        assert new_output_vcf == variants.get_output()


def test_set_get_config():
    """
    This function tests the set_config and get_config methods of the Variants class in Python.
    """

    # Init files
    input_vcf = tests_data_folder + "/example.vcf.gz"
    input_config = {"option": "option_value"}
    new_input_config = {"option": "option_value", "new_option": "new_option_value"}

    # Create connection
    conn = duckdb.connect(":memory:")

    # Create object
    variants = Variants(conn=conn, input=input_vcf, config=input_config)

    # Check get_input
    assert input_config == variants.get_config()

    # Check new input vcf
    variants.set_config(new_input_config)
    assert new_input_config == variants.get_config()


def test_set_get_param():
    """
    This function tests the set_param and get_param methods of the Variants class in Python.
    """

    # Init files
    input_vcf = tests_data_folder + "/example.vcf.gz"
    input_param = {"option": "option_value"}
    new_input_param = {"option": "option_value", "new_option": "new_option_value"}

    # Create connection
    conn = duckdb.connect(":memory:")

    # Create object
    variants = Variants(conn=conn, input=input_vcf, param=input_param)

    # Check get_input
    assert input_param == variants.get_param()

    # Check new input vcf
    variants.set_param(new_input_param)
    assert new_input_param == variants.get_param()


def test_get_header():
    """
    This function tests various methods related to getting and setting the header of a VCF file using
    the Variants class in Python.
    """

    # Init files
    input_vcf = tests_data_folder + "/example.vcf.gz"

    # Create connection
    conn = duckdb.connect(":memory:")

    # Create object
    variants = Variants(conn=conn, input=input_vcf)

    # set_header done when vcf object creation

    # Check header VCF
    header_vcf = variants.get_header()
    assert header_vcf.infos != None

    # Check header List and nb
    header_list = variants.get_header(type="list")
    assert header_list != []
    assert len(header_list) == 53

    # check header length
    assert variants.get_header_length() == 52

    # check get_header_columns
    header_columns = variants.get_header_columns().strip()
    header_columns_expected = "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	sample1	sample2	sample3	sample4".strip()
    assert header_columns == header_columns_expected

    # check get_header_columns_as_sql
    header_columns_as_sql = variants.get_header_columns_as_sql().strip()
    header_columns_as_sql_expected = """ "#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","sample1","sample2","sample3","sample4" """.strip()
    assert header_columns_as_sql == header_columns_as_sql_expected

    # check get_header_sample_list
    header_columns_sample_list = variants.get_header_sample_list()
    header_columns_sample_list_expected = ["sample1", "sample2", "sample3", "sample4"]
    assert header_columns_sample_list == header_columns_sample_list_expected


def test_get_header_infos_list():
    """
    The function `test_get_header_infos_list` tests the retrieval of information headers from a VCF file
    using the `Variants` class.
    """

    # Init files
    input_vcf = tests_data_folder + "/example.vcf.gz"

    # Create connection
    conn = duckdb.connect(":memory:")

    # Create object
    variants = Variants(conn=conn, input=input_vcf)

    # set_header done when vcf object creation

    # Check header VCF
    header_vcf = variants.get_header()
    assert header_vcf.infos != None

    # Check header List and nb
    header_list = variants.get_header_infos_list()
    assert header_list == ["NS", "DP", "AA", "CLNSIG", "SIFT"]


def test_get_header_no_samples():
    """
    This function tests various methods related to getting and setting the header of a VCF file when
    there are no samples present.
    """

    # Init files
    input_vcf = tests_data_folder + "/example.no_samples.vcf.gz"

    # Create connection
    conn = duckdb.connect(":memory:")

    # Create object
    variants = Variants(conn=conn, input=input_vcf)

    # set_header done when vcf object creation

    # Check header VCF
    header_vcf = variants.get_header()
    assert header_vcf.infos != None

    # Check header List and nb
    header_list = variants.get_header(type="list")
    assert header_list != []
    assert len(header_list) == 38

    # check header length
    assert variants.get_header_length() == 37

    # check get_header_columns
    header_columns = variants.get_header_columns().strip()
    header_columns_expected = "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO".strip()
    assert header_columns == header_columns_expected

    # check get_header_columns_as_sql
    header_columns_as_sql = variants.get_header_columns_as_sql().strip()
    header_columns_as_sql_expected = (
        """ "#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO" """.strip()
    )
    assert header_columns_as_sql == header_columns_as_sql_expected

    # check get_header_sample_list
    header_columns_sample_list = variants.get_header_sample_list()
    header_columns_sample_list_expected = []
    assert header_columns_sample_list == header_columns_sample_list_expected


def test_get_header_in_config():
    """
    This function tests various methods related to getting the header information from a Variants object
    in Python.
    """

    # Init files
    input_vcf = tests_data_folder + "/example.parquet"
    input_config = {"header_file": tests_data_folder + "/example.parquet.hdr"}

    # Create object
    variants = Variants(input=input_vcf, config=input_config)

    # Check header VCF
    header_vcf = variants.get_header()
    assert header_vcf.infos != None

    # Check header List and nb
    header_list = variants.get_header(type="list")
    assert header_list != []
    assert len(header_list) == 53

    # check header length
    assert variants.get_header_length() == 52

    # check get_header_columns
    header_columns = variants.get_header_columns().strip()
    header_columns_expected = "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	sample1	sample2	sample3	sample4".strip()
    assert header_columns == header_columns_expected

    # check get_header_columns_as_sql
    header_columns_as_sql = variants.get_header_columns_as_sql().strip()
    header_columns_as_sql_expected = """ "#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","sample1","sample2","sample3","sample4" """.strip()
    assert header_columns_as_sql == header_columns_as_sql_expected

    # check get_header_sample_list
    header_columns_sample_list = variants.get_header_sample_list()
    header_columns_sample_list_expected = ["sample1", "sample2", "sample3", "sample4"]
    assert header_columns_sample_list == header_columns_sample_list_expected


def test_load_without_header():
    """
    This function tests if a ValueError is raised when attempting to load a Parquet file without a
    header using the Variants object.
    """

    # Init files
    input_vcf = tests_data_folder + "/example.without_header.parquet"

    # Create object
    try:
        variants = Variants(input=input_vcf)
        assert True
    except:
        assert False


def test_read_vcf_header():
    """
    This function tests the read_vcf_header method of the Variants class by checking if the header list
    is not empty and has a length of 53.
    """

    # Init files
    input_vcf = tests_data_folder + "/example.parquet"
    vcf_header = tests_data_folder + "/example.parquet.hdr"

    # Create connection
    conn = duckdb.connect(":memory:")

    # Create object
    variants = Variants(conn=conn, input=input_vcf)

    # Check read_vcf_header
    with open(vcf_header, "rt") as f:
        header_list = variants.read_vcf_header(f)
    assert header_list != []
    assert len(header_list) == 53


def test_load_when_init():
    """
    This function tests if the Variants object loads data correctly from a VCF file during
    initialization.
    """

    # Init files
    input_vcf = tests_data_folder + "/example.vcf.gz"

    # Create object
    variants = Variants(input=input_vcf, load=True)

    # Check data loaded
    result = variants.get_query_to_df("SELECT count(*) AS count FROM variants")
    nb_variant_in_database = result["count"][0]

    expected_number_of_variants = 7

    assert nb_variant_in_database == expected_number_of_variants


def test_load_format_not_available():
    """
    This is a test function that checks if an input file format is not available.
    """

    # Init files
    input_format = "unknwon"
    input_vcf = tests_folder + f"/data/example.{input_format}"
    input_config = {"header_file": tests_data_folder + "/example.parquet.hdr"}

    # Create object
    with pytest.raises(ValueError) as e:
        variants = Variants(input=input_vcf, config=input_config, load=True)
    assert str(e.value) == f"Input file format '{input_format}' not available"


def test_load_vcf_gz():
    """
    This function tests if a VCF file in gzipped format can be loaded into a Variants object and if the
    expected number of variants is present in the resulting database.
    """

    # Init files
    input_vcf = tests_data_folder + "/example.vcf.gz"

    # Create object
    variants = Variants(input=input_vcf, load=True)

    # Check data loaded
    result = variants.get_query_to_df("SELECT count(*) AS count FROM variants")
    nb_variant_in_database = result["count"][0]

    expected_number_of_variants = 7

    assert nb_variant_in_database == expected_number_of_variants


def test_load_full_unsorted_vcf_gz():
    """
    This function tests if a VCF file can be loaded into a Variants object and the number of variants in
    the file matches the expected number.
    The VCF file is full of various variants type, and unsorted
    """

    # Init files
    input_vcf = tests_data_folder + "/example.full.unsorted.vcf.gz"

    # Create object
    variants = Variants(input=input_vcf, load=True)

    # Check data loaded
    result = variants.get_query_to_df("SELECT count(*) AS count FROM variants")
    nb_variant_in_database = result["count"][0]

    expected_number_of_variants = 36

    assert nb_variant_in_database == expected_number_of_variants


def test_load_vcf_gz_no_samples():
    """
    This function tests if a VCF file with no samples can be loaded into a Variants object and the
    number of variants in the database matches the expected number.
    """

    # Init files
    input_vcf = tests_data_folder + "/example.no_samples.vcf.gz"

    # Create object
    variants = Variants(input=input_vcf, load=True)

    # Check data loaded
    result = variants.get_query_to_df("SELECT count(*) AS count FROM variants")
    nb_variant_in_database = result["count"][0]

    expected_number_of_variants = 10

    assert nb_variant_in_database == expected_number_of_variants


def test_load_parquet():
    """
    This function tests if a Parquet file can be loaded into a Variants object and if the expected
    number of variants is present in the database.
    """

    # Init files
    input_vcf = tests_data_folder + "/example.parquet"

    # Create object
    variants = Variants(input=input_vcf, load=True)

    # Check data loaded
    result = variants.get_query_to_df("SELECT count(*) AS count FROM variants")
    nb_variant_in_database = result["count"][0]

    expected_number_of_variants = 7

    assert nb_variant_in_database == expected_number_of_variants


def test_load_vcf():
    """
    This function tests if a VCF file can be loaded into a Variants object and if the expected number of
    variants is present in the resulting database.
    """

    # Init files
    input_vcf = tests_data_folder + "/example.vcf"

    # Create object
    variants = Variants(input=input_vcf, load=True)

    # Check data loaded
    result = variants.get_query_to_df("SELECT count(*) AS count FROM variants")
    nb_variant_in_database = result["count"][0]

    expected_number_of_variants = 7

    assert nb_variant_in_database == expected_number_of_variants


def test_load_csv():
    """
    This function tests if a CSV file can be loaded into a Variants object and if the expected number of
    variants is present in the database.
    """

    # Init files
    input_vcf = tests_data_folder + "/example.csv"

    # Create object
    variants = Variants(input=input_vcf, load=True)

    # Check data loaded
    result = variants.get_query_to_df("SELECT count(*) AS count FROM variants")
    nb_variant_in_database = result["count"][0]

    expected_number_of_variants = 7

    assert nb_variant_in_database == expected_number_of_variants


def test_load_tsv():
    """
    This function tests if a TSV file can be loaded into a Variants object and if the expected number of
    variants is present in the database.
    """

    # Init files
    input_vcf = tests_data_folder + "/example.tsv"

    # Create object
    variants = Variants(input=input_vcf, load=True)

    # Check data loaded
    result = variants.get_query_to_df("SELECT count(*) AS count FROM variants")
    nb_variant_in_database = result["count"][0]

    expected_number_of_variants = 7

    assert nb_variant_in_database == expected_number_of_variants


def test_load_duckdb():
    """
    This function tests if a DuckDB database containing variant data can be loaded and queried
    correctly.
    """

    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # Init files
        input_vcf = tests_data_folder + "/example.vcf.gz"
        output_duckdb = f"{tmp_dir}/example.duckdb"

        remove_if_exists([output_duckdb])

        duckdb_file = Variants(input=input_vcf, output=output_duckdb, load=True)
        duckdb_file.export_output()

        # Test load duckdb database

        # Init files
        input_vcf = output_duckdb

        # Create object
        variants = Variants(input=input_vcf, load=True)

        # Check data loaded
        result = variants.get_query_to_df("SELECT count(*) AS count FROM variants")
        nb_variant_in_database = result["count"][0]

        expected_number_of_variants = 7

        assert nb_variant_in_database == expected_number_of_variants


def test_get_connexion_db_memory():
    """
    This function tests if the `get_connexion_db()` method of the `Variants` class returns the expected
    value for a memory-based database connection.
    """

    # Init files
    input_vcf = tests_data_folder + "/example.vcf.gz"
    input_config = {"connexion_type": "memory"}

    # Create object
    variants = Variants(input=input_vcf, config=input_config)

    # get connexion_db
    connexion_db = variants.get_connexion_db()

    assert connexion_db == ":memory:"


def test_get_connexion_db_tmpfile():
    """
    This function tests the "get_connexion_db_tmpfile" method of the "Variants" class in Python.
    """

    # Init files
    input_vcf = tests_data_folder + "/example.vcf.gz"
    input_config = {"connexion_type": "tmpfile"}

    # Create object
    variants = Variants(input=input_vcf, config=input_config)

    # get connexion_db
    connexion_db = variants.get_connexion_db()

    assert os.path.exists(connexion_db)

    remove_if_exists(connexion_db)


def test_get_connexion_db_file():
    """
    This function tests the "get_connexion_db" method of a "Variants" object in Python.
    """

    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # Init files
        input_vcf = tests_data_folder + "/example.vcf.gz"
        input_config = {"connexion_type": f"{tmp_dir}/connexion.duckdb"}

        # Remove if exists
        remove_if_exists([f"{tmp_dir}/connexion.duckdb"])

        # Create object
        variants = Variants(input=input_vcf, config=input_config)

        # get connexion_db
        connexion_db = variants.get_connexion_db()

        assert connexion_db == f"{tmp_dir}/connexion.duckdb"


def test_get_table_variants():
    """
    This function tests the get_table_variants method of the Variants class in Python.
    """

    # Init files
    input_vcf = tests_data_folder + "/example.vcf.gz"

    # Create object
    variants = Variants(input=input_vcf)

    # get connexion_db
    table_variants_select = variants.get_table_variants(clause="select")
    assert table_variants_select == "variants"
    table_variants_from = variants.get_table_variants(clause="from")
    assert table_variants_from == "variants as variants"
    table_variants_else = variants.get_table_variants(clause="else")
    assert table_variants_else == "variants"


def test_get_connexion():
    """
    This function tests the ability of a Variants object to establish a database connection and execute
    a query.
    """

    # Init files
    input_vcf = tests_data_folder + "/example.vcf.gz"

    # Create object
    variants = Variants(input=input_vcf)

    # get connexion_db
    connexion = variants.get_connexion()
    result = connexion.query("SELECT 'pass' AS connexion")
    check_connexion = result.df()["connexion"][0] == "pass"

    assert check_connexion


def test_get_connexion_sqlite():
    """
    This function tests the ability to establish a connection to a SQLite database using the Variants
    class.
    """

    # Init files
    input_vcf = tests_data_folder + "/example.vcf.gz"
    input_config = {"connexion_format": "sqlite"}

    # Create object
    variants = Variants(input=input_vcf, config=input_config)

    # get connexion_db
    connexion = variants.get_connexion()
    result = connexion.execute("SELECT 'pass' AS connexion").fetchall()
    assert result[0][0] == "pass"


def test_get_verbose():
    """
    This function tests the "get_verbose" method of the "Variants" class in Python.
    """

    # Init files
    input_vcf = tests_data_folder + "/example.vcf.gz"

    # Create object
    variants = Variants(input=input_vcf)

    # get connexion_db
    verbose = variants.get_verbose()
    assert not verbose

    # config verbose True
    input_config = {"verbose": True}

    # Create object
    variants = Variants(input=input_vcf, config=input_config)

    # get connexion_db
    verbose = variants.get_verbose()
    assert verbose


def test_load_connexion_type_memory():
    """
    This function tests if a Variants object can be created and data can be loaded into memory from a
    VCF file.
    """

    # Init files
    input_vcf = tests_data_folder + "/example.vcf.gz"
    input_config = {"connexion_type": "memory"}

    # Create object
    variants = Variants(input=input_vcf, config=input_config, load=True)

    # Check data loaded
    result = variants.get_query_to_df("SELECT count(*) AS count FROM variants")
    nb_variant_in_database = result["count"][0]

    expected_number_of_variants = 7

    assert nb_variant_in_database == expected_number_of_variants


def test_load_connexion_type_tmpfile():
    """
    This is a unit test function that checks if the number of variants loaded from a VCF file using a
    temporary file connection type matches the expected number.
    """

    # Init files
    input_vcf = tests_data_folder + "/example.vcf.gz"
    input_config = {"connexion_type": "tmpfile"}

    # Create object
    variants = Variants(input=input_vcf, config=input_config, load=True)

    # Check data loaded
    result = variants.get_query_to_df("SELECT count(*) AS count FROM variants")
    nb_variant_in_database = result["count"][0]

    expected_number_of_variants = 7

    assert nb_variant_in_database == expected_number_of_variants


def test_load_connexion_type_file():
    """
    This function tests the loading of a VCF file into a DuckDB database and checks if the expected
    number of variants is present in the database.
    """

    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # Init files
        input_vcf = tests_data_folder + "/example.vcf.gz"
        input_config = {"connexion_type": f"{tmp_dir}/output.duckdb"}

        remove_if_exists([f"{tmp_dir}/output.duckdb"])

        # Create object
        variants = Variants(input=input_vcf, config=input_config, load=True)

        # Check data loaded
        result = variants.get_query_to_df("SELECT count(*) AS count FROM variants")
        nb_variant_in_database = result["count"][0]

        expected_number_of_variants = 7

        assert nb_variant_in_database == expected_number_of_variants


def test_load_connexion_format_sqlite():
    """
    This function tests if a SQLite database is properly loaded with data from a VCF file.
    """

    # Init files
    input_vcf = tests_data_folder + "/example.vcf.gz"
    input_config = {"connexion_format": "sqlite"}

    # Create object
    variants = Variants(input=input_vcf, config=input_config, load=True)

    # Check data loaded
    result = variants.get_query_to_df("SELECT count(*) AS count FROM variants")
    nb_variant_in_database = result["count"][0]

    expected_number_of_variants = 7

    assert nb_variant_in_database == expected_number_of_variants


###
### Export Output
###


def test_export_output_vcf_gz():
    """
    This function tests the export of a VCF file in gzipped format with the pyVCF library.
    """

    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # Init files
        input_vcf = tests_data_folder + "/example.vcf.gz"
        output_vcf = f"{tmp_dir}/example.vcf.gz"

        # remove if exists
        remove_if_exists([output_vcf])

        # Create object
        variants = Variants(input=input_vcf, output=output_vcf, load=True)

        # Check get_output
        variants.export_output()
        assert os.path.exists(output_vcf) and not os.path.exists(output_vcf + ".hdr")

        # Check if VCF is in correct format with pyVCF
        try:
            vcf.Reader(filename=output_vcf)
        except:
            assert False


def test_export_output_vcf_gz_from_full_unsorted():
    """
    Tests the export of a VCF file in compressed format from an unsorted input file.

    This function performs a series of tests on the Variants object to check that it can successfully export a VCF file
    in compressed format from an unsorted input file. It also checks that the exported VCF file is in the correct format
    using the pyVCF library.

    :raises AssertionError: If any of the tests fail.

    :return: None
    """

    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # Init files
        input_vcf = tests_data_folder + "/example.full.unsorted.vcf.gz"
        output_vcf = f"{tmp_dir}/example.vcf.gz"

        # remove if exists
        remove_if_exists([output_vcf])

        # Create object
        variants = Variants(input=input_vcf, output=output_vcf, load=True)

        # Check get_output
        variants.export_output()
        assert os.path.exists(output_vcf) and not os.path.exists(output_vcf + ".hdr")

        # Check if VCF is in correct format with pyVCF
        try:
            vcf.Reader(filename=output_vcf)
        except:
            assert False


def test_export_output_vcf():
    """
    This function tests the export_output method of the Variants class in Python, which exports a VCF
    file and checks if it is in the correct format using pyVCF.
    """

    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # Init files
        input_vcf = tests_data_folder + "/example.vcf.gz"
        output_vcf = f"{tmp_dir}/example.vcf"

        # remove if exists
        remove_if_exists([output_vcf])

        # Create object
        variants = Variants(input=input_vcf, output=output_vcf, load=True)

        # Check get_output
        variants.export_output()
        assert os.path.exists(output_vcf) and not os.path.exists(output_vcf + ".hdr")

        # Check if VCF is in correct format with pyVCF
        vcf.Reader(filename=output_vcf)
        try:
            vcf.Reader(filename=output_vcf)
        except:
            assert False


def test_export_output_parquet():
    """
    This function tests the export of a VCF file to a Parquet file format.
    """

    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # Init files
        input_vcf = tests_data_folder + "/example.vcf.gz"
        output_vcf = f"{tmp_dir}/example.parquet"

        # remove if exists
        remove_if_exists([output_vcf])

        # Create object
        variants = Variants(input=input_vcf, output=output_vcf, load=True)

        # Check get_output
        variants.export_output()
        assert os.path.exists(output_vcf) and os.path.exists(output_vcf + ".hdr")

        # Check get_output without header
        remove_if_exists([output_vcf, output_vcf + ".hdr"])
        variants.export_output(export_header=False)
        assert os.path.exists(output_vcf) and not os.path.exists(output_vcf + ".hdr")


def test_export_output_duckdb():
    """
    This function tests the export_output method of the Variants class in Python's DuckDB library.
    """

    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # Init files
        input_vcf = tests_data_folder + "/example.vcf.gz"
        output_duckdb = f"{tmp_dir}/example.duckdb"

        # remove if exists
        remove_if_exists([output_duckdb])

        # Create object
        variants = Variants(input=input_vcf, output=output_duckdb, load=True)

        # Check get_output
        variants.export_output()
        assert os.path.exists(output_duckdb) and os.path.exists(output_duckdb + ".hdr")

        # remove if exists
        remove_if_exists([output_duckdb, output_duckdb + ".hdr"])

        # Check get_output without header
        variants.export_output(export_header=False)
        assert os.path.exists(output_duckdb) and not os.path.exists(
            output_duckdb + ".hdr"
        )


def test_export_output_tsv():
    """
    This function tests the export_output method of the Variants class in Python, which exports a VCF
    file to a TSV file.
    """

    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # Init files
        input_vcf = tests_data_folder + "/example.vcf.gz"
        output_vcf = f"{tmp_dir}/example.tsv"

        # remove if exists
        remove_if_exists([output_vcf])

        # Create object
        variants = Variants(input=input_vcf, output=output_vcf, load=True)

        # Check get_output
        variants.export_output()
        assert os.path.exists(output_vcf) and os.path.exists(output_vcf + ".hdr")

        # Check get_output without header
        remove_if_exists([output_vcf, output_vcf + ".hdr"])
        variants.export_output(export_header=False)
        assert os.path.exists(output_vcf) and not os.path.exists(output_vcf + ".hdr")


def test_export_output_tsv_gz():
    """
    This function tests the export_output method of the Variants class in Python, which exports a VCF
    file to a TSV file in gzipped format.
    """

    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # Init files
        input_vcf = tests_data_folder + "/example.vcf.gz"
        output_vcf = f"{tmp_dir}/example.tsv.gz"

        # remove if exists
        remove_if_exists([output_vcf])

        # Create object
        variants = Variants(input=input_vcf, output=output_vcf, load=True)

        # Check get_output
        variants.export_output()
        assert os.path.exists(output_vcf) and os.path.exists(output_vcf + ".hdr")

        # Check get_output without header
        remove_if_exists([output_vcf, output_vcf + ".hdr"])
        variants.export_output(export_header=False)
        assert os.path.exists(output_vcf) and not os.path.exists(output_vcf + ".hdr")


def test_export_output_csv():
    """
    This function tests the export_output method of the Variants class in Python by checking if the
    output file exists with and without a header.
    """

    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # Init files
        input_vcf = tests_data_folder + "/example.vcf.gz"
        output_vcf = f"{tmp_dir}/example.csv"

        # remove if exists
        remove_if_exists([output_vcf])

        # Create object
        variants = Variants(input=input_vcf, output=output_vcf, load=True)

        # Check get_output
        variants.export_output()
        assert os.path.exists(output_vcf) and os.path.exists(output_vcf + ".hdr")

        # Check get_output without header
        remove_if_exists([output_vcf, output_vcf + ".hdr"])
        variants.export_output(export_header=False)
        assert os.path.exists(output_vcf) and not os.path.exists(output_vcf + ".hdr")


def test_export_output_tbl():
    """
    This function tests the export_output method of the Variants class in Python.
    """

    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # Init files
        input_vcf = tests_data_folder + "/example.vcf.gz"
        output_vcf = f"{tmp_dir}/example.tbl"

        # remove if exists
        remove_if_exists([output_vcf])

        # Create object
        variants = Variants(input=input_vcf, output=output_vcf, load=True)

        # Check get_output
        variants.export_output()
        assert os.path.exists(output_vcf) and os.path.exists(output_vcf + ".hdr")

        # Check get_output without header
        remove_if_exists([output_vcf, output_vcf + ".hdr"])
        variants.export_output(export_header=False)
        assert os.path.exists(output_vcf) and not os.path.exists(output_vcf + ".hdr")


def test_export_output_tsv_explode_infos():
    """
    This function tests the export of variant information in TSV format with exploded extra information.
    """

    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # Init files
        input_vcf = tests_data_folder + "/example.vcf.gz"
        output_vcf = f"{tmp_dir}/example.tsv"
        param = {"export_extra_infos": True}

        # remove if exists
        remove_if_exists([output_vcf])

        # Create object
        variants = Variants(input=input_vcf, output=output_vcf, load=True, param=param)

        # Explode infos
        variants.explode_infos()

        # Check get_output
        variants.export_output()
        assert os.path.exists(output_vcf) and os.path.exists(output_vcf + ".hdr")

        # Check get_output without header
        remove_if_exists([output_vcf, output_vcf + ".hdr"])
        variants.export_output(export_header=False)
        assert os.path.exists(output_vcf) and not os.path.exists(output_vcf + ".hdr")


def test_export_from_sqlite_output_vcf_gz():
    """
    This function tests the export of a VCF file in gzipped format with the pyVCF library.
    """

    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # Init files
        input_vcf = tests_data_folder + "/example.vcf.gz"
        output_vcf = f"{tmp_dir}/example.vcf.gz"
        input_config = {"connexion_format": "sqlite"}

        # remove if exists
        remove_if_exists([output_vcf])

        # Create object
        variants = Variants(
            input=input_vcf, output=output_vcf, config=input_config, load=True
        )

        # Check get_output
        variants.export_output()
        assert os.path.exists(output_vcf) and not os.path.exists(output_vcf + ".hdr")

        # Check get_output without header
        remove_if_exists([output_vcf, output_vcf + ".hdr"])
        variants.export_output(export_header=False)
        assert os.path.exists(output_vcf) and not os.path.exists(output_vcf + ".hdr")

        # Check if VCF is in correct format with pyVCF
        try:
            vcf.Reader(filename=output_vcf)
        except:
            assert False


def test_export_from_sqlite_output_vcf():
    """
    This function tests the export of a VCF file in gzipped format with the pyVCF library.
    """

    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # Init files
        input_vcf = tests_data_folder + "/example.vcf.gz"
        output_vcf = f"{tmp_dir}/example.vcf"
        input_config = {"connexion_format": "sqlite"}

        # remove if exists
        remove_if_exists([output_vcf])

        # Create object
        variants = Variants(
            input=input_vcf, output=output_vcf, config=input_config, load=True
        )

        # Check get_output
        variants.export_output()
        assert os.path.exists(output_vcf) and not os.path.exists(output_vcf + ".hdr")

        # Check get_output without header
        remove_if_exists([output_vcf, output_vcf + ".hdr"])
        variants.export_output(export_header=False)
        assert os.path.exists(output_vcf) and not os.path.exists(output_vcf + ".hdr")

        # Check if VCF is in correct format with pyVCF
        try:
            vcf.Reader(filename=output_vcf)
        except:
            assert False


def test_export_from_sqlite_output_parquet():
    """
    This function tests the export of a VCF file in gzipped format with the pyVCF library.
    """

    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # Init files
        input_vcf = tests_data_folder + "/example.vcf.gz"
        output_vcf = f"{tmp_dir}/example.parquet"
        input_config = {"connexion_format": "sqlite"}

        # remove if exists
        remove_if_exists([output_vcf])

        # Create object
        variants = Variants(
            input=input_vcf, output=output_vcf, config=input_config, load=True
        )

        # Check get_output
        variants.export_output()
        assert os.path.exists(output_vcf) and os.path.exists(output_vcf + ".hdr")

        # Check get_output without header
        remove_if_exists([output_vcf, output_vcf + ".hdr"])
        variants.export_output(export_header=False)
        assert os.path.exists(output_vcf) and not os.path.exists(output_vcf + ".hdr")


def test_export_from_sqlite_output_tsv():
    """
    This function tests the export of a VCF file in gzipped format with the pyVCF library.
    """

    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # Init files
        input_vcf = tests_data_folder + "/example.vcf.gz"
        output_vcf = f"{tmp_dir}/example.tsv"
        input_config = {"connexion_format": "sqlite"}

        # remove if exists
        remove_if_exists([output_vcf])

        # Create object
        variants = Variants(
            input=input_vcf, output=output_vcf, config=input_config, load=True
        )

        # Check get_output
        variants.export_output()
        assert os.path.exists(output_vcf) and os.path.exists(output_vcf + ".hdr")

        # Check get_output without header
        remove_if_exists([output_vcf, output_vcf + ".hdr"])
        variants.export_output(export_header=False)
        assert os.path.exists(output_vcf) and not os.path.exists(output_vcf + ".hdr")

        # Check if VCF is in correct format with pyVCF
        try:
            vcf.Reader(filename=output_vcf)
        except:
            assert False


def test_export_from_sqlite_output_tsv_gz():
    """
    This function tests the export of a VCF file in gzipped format with the pyVCF library.
    """

    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # Init files
        input_vcf = tests_data_folder + "/example.vcf.gz"
        output_vcf = f"{tmp_dir}/example.tsv.gz"
        input_config = {"connexion_format": "sqlite"}

        # remove if exists
        remove_if_exists([output_vcf])

        # Create object
        variants = Variants(
            input=input_vcf, output=output_vcf, config=input_config, load=True
        )

        # Check get_output
        variants.export_output()
        assert os.path.exists(output_vcf) and os.path.exists(output_vcf + ".hdr")

        # Check get_output without header
        remove_if_exists([output_vcf, output_vcf + ".hdr"])
        variants.export_output(export_header=False)
        assert os.path.exists(output_vcf) and not os.path.exists(output_vcf + ".hdr")

        # Check if VCF is in correct format with pyVCF
        try:
            vcf.Reader(filename=output_vcf)
        except:
            assert False


###
### Explode
###


def test_explode_infos():
    """
    This function tests the explode_infos method of the Variants class in Python.
    """

    # Init files
    input_vcf = tests_data_folder + "/example.vcf.gz"
    param = {"explode": {"explode_infos": True}}

    # Create object
    variants = Variants(input=input_vcf, load=True, param=param)

    # column to check
    column_to_check = "CLNSIG"
    value_to_check = "pathogenic"

    # check column found
    result = variants.execute_query("SELECT * FROM variants LIMIT 0")
    log.debug(f"result={result}")
    assert column_to_check in [col[0] for col in result.description]

    # Check value in column
    result = variants.get_query_to_df(
        f"""SELECT "{column_to_check}" AS column_to_check FROM variants WHERE "#CHROM" = 'chr1' AND POS = 28736 AND REF = 'A' AND ALT = 'C' """
    )
    assert value_to_check == result["column_to_check"][0]

    # Check number of value in column to check
    result = variants.get_query_to_df(
        f"""SELECT "{column_to_check}" AS column_to_check FROM variants WHERE "{column_to_check}" IS NOT NULL """
    )
    assert len(result) == 2


def test_explode_infos_custom():
    """
    This function tests the explode_infos method of the Variants class in Python.
    """

    # Init files
    input_vcf = tests_data_folder + "/example.vcf.gz"
    param = {"explode": {"explode_infos": True, "explode_infos_prefix": "CUSTOM_"}}

    # Create object
    variants = Variants(input=input_vcf, load=True, param=param)

    # column to check
    column_to_check = "CUSTOM_CLNSIG"
    value_to_check = "pathogenic"

    # check column found
    result = variants.execute_query("SELECT * FROM variants LIMIT 0")
    assert column_to_check in [col[0] for col in result.description]

    # Check value in column
    result = variants.get_query_to_df(
        f"""SELECT "{column_to_check}" AS column_to_check FROM variants WHERE "#CHROM" = 'chr1' AND POS = 28736 AND REF = 'A' AND ALT = 'C' """
    )
    assert value_to_check == result["column_to_check"][0]

    # Check number of value in column to check
    result = variants.get_query_to_df(
        f"""SELECT "{column_to_check}" AS column_to_check FROM variants WHERE "{column_to_check}" IS NOT NULL """
    )
    assert len(result) == 2


def test_explode_infos_method():
    """
    This function tests the explode_infos method of the Variants class in Python.
    """

    # Init files
    input_vcf = tests_data_folder + "/example.vcf.gz"

    # Create object
    variants = Variants(input=input_vcf, load=True)

    # Explode infos fields
    variants.explode_infos(prefix="CUSTOM_", fields=["*"])

    # column to check
    column_to_check = "CUSTOM_CLNSIG"
    value_to_check = "pathogenic"

    # check column found
    result = variants.execute_query("SELECT * FROM variants LIMIT 0")
    assert column_to_check in [col[0] for col in result.description]

    # Check value in column
    result = variants.get_query_to_df(
        f"""SELECT "{column_to_check}" AS column_to_check FROM variants WHERE "#CHROM" = 'chr1' AND POS = 28736 AND REF = 'A' AND ALT = 'C' """
    )
    assert value_to_check == result["column_to_check"][0]

    # Check number of value in column to check
    result = variants.get_query_to_df(
        f"""SELECT "{column_to_check}" AS column_to_check FROM variants WHERE "{column_to_check}" IS NOT NULL """
    )
    assert len(result) == 2


def test_explode_infos_no_infos():
    """
    This function tests if the columns in a VCF file remain the same before and after exploding the info
    fields.
    """

    # Init files
    input_vcf = tests_data_folder + "/example.no_samples.vcf.gz"

    # Create object
    variants = Variants(input=input_vcf, load=True)

    # Check column before explode
    result = variants.execute_query("SELECT * FROM variants LIMIT 0")
    columns_before_explode = [col[0] for col in result.description]

    # Explode infos fields
    variants.explode_infos()

    # check column found
    result = variants.execute_query("SELECT * FROM variants LIMIT 0")
    columns_after_explode = [col[0] for col in result.description]

    assert columns_before_explode == columns_after_explode


def test_explode_infos_sqlite():
    """
    This function tests the explode_infos() method and checks if a specific value is present in a column
    of a SQLite database created from a VCF file.
    """

    # Init files
    input_vcf = tests_data_folder + "/example.vcf.gz"
    input_config = {"connexion_format": "sqlite"}

    # Create object
    variants = Variants(input=input_vcf, config=input_config, load=True)

    # Explode infos fields
    variants.explode_infos()

    # Annotation
    variants.annotation()

    # column to check
    column_to_check = "CLNSIG"
    value_to_check = "pathogenic"

    # check column found
    result = variants.execute_query("SELECT * FROM variants LIMIT 0")
    assert column_to_check in [col[0] for col in result.description]

    # Check value in column
    result = variants.get_query_to_df(
        f"""SELECT "{column_to_check}" AS column_to_check FROM variants WHERE "#CHROM" = 'chr1' AND POS = 28736 AND REF = 'A' AND ALT = 'C' """
    )
    assert value_to_check == result["column_to_check"][0]

    # Check number of value in column to check
    result = variants.get_query_to_df(
        f"""SELECT "{column_to_check}" AS column_to_check FROM variants WHERE "{column_to_check}" IS NOT NULL """
    )
    assert len(result) == 2


def test_explode_infos_param_prefix():
    """
    This function tests the functionality of exploding info parameters in a VCF file using Python.
    """

    # Init files
    input_vcf = tests_data_folder + "/example.vcf.gz"
    infos_prefix = "INFO_"
    input_param = {
        "explode": {"explode_infos": True, "explode_infos_prefix": infos_prefix}
    }

    # Create object
    variants = Variants(input=input_vcf, load=True, param=input_param)

    # column to check
    column_to_check = infos_prefix + "CLNSIG"
    value_to_check = "pathogenic"

    # check column found
    result = variants.execute_query("SELECT * FROM variants LIMIT 0")
    assert column_to_check in [col[0] for col in result.description]

    # Check value in column
    result = variants.get_query_to_df(
        f"""SELECT "{column_to_check}" AS column_to_check FROM variants WHERE "#CHROM" = 'chr1' AND POS = 28736 AND REF = 'A' AND ALT = 'C' """
    )
    assert value_to_check == result["column_to_check"][0]

    # Check number of value in column to check
    result = variants.get_query_to_df(
        f"""SELECT "{column_to_check}" AS column_to_check FROM variants WHERE "{column_to_check}" IS NOT NULL """
    )
    assert len(result) == 2


###
### Overview and Stats
###


def test_overview():
    """
    This function tests the get_overview method of the Variants class by creating an object and checking
    if the overview is None.
    """

    # Init files
    input_vcf = tests_data_folder + "/example.vcf.gz"

    # Create object
    variants = Variants(input=input_vcf, load=True)

    # Overview
    overview = variants.get_overview()

    assert overview == None


def test_overview_no_samples():
    """
    This function tests the get_overview method of the Variants class when there are no samples in the
    input VCF file.
    """

    # Init files
    input_vcf = tests_data_folder + "/example.no_samples.vcf.gz"

    # Create object
    variants = Variants(input=input_vcf, load=True)

    # Overview
    overview = variants.get_overview()

    assert overview == None


def test_stats():
    """
    This is a unit test function for the get_stats method of the Variants class in Python.
    """

    # Init files
    input_vcf = tests_data_folder + "/example.vcf.gz"

    # Create object
    variants = Variants(input=input_vcf, load=True)

    # Stats
    stats = variants.get_stats()

    assert isinstance(stats, dict) and len(stats)


def test_stats_no_samples():
    """
    This function tests if the get_stats() method returns None when there are no samples in the input
    VCF file.
    """

    # Init files
    input_vcf = tests_data_folder + "/example.no_samples.vcf.gz"

    # Create object
    variants = Variants(input=input_vcf, load=True)

    # Stats
    stats = variants.get_stats()

    assert isinstance(stats, dict) and len(stats)


###
### No input file
###


def test_no_input_file():
    """
    This function tests the behavior of the Variants class when no input file is provided.
    """

    # Create object
    variant = Variants()

    assert variant.get_input() == None
    assert len(variant.get_header().infos) == 0
    assert len(variant.get_header().formats) == 0
    assert len(variant.get_header().contigs) == 0
    assert variant.get_header_length() == 1
    assert variant.get_header_columns() == "\t".join(vcf_required_columns)
    assert variant.get_header_columns_as_list() == vcf_required_columns
    assert variant.get_header_columns_as_sql() == ",".join(
        f'"{col}"' for col in vcf_required_columns
    )


###
### Query
###


def test_query():
    """
    This function tests the execute_query method of the Variants class in Python by checking if it
    returns the expected results for a given query and for a None query.
    """

    # Create object
    variants = Variants()

    # Query
    query = "SELECT 1 AS query"
    result_query = variants.execute_query(query)
    assert result_query.df()["query"][0] == 1

    # Query none
    result_query = variants.execute_query(None)
    assert result_query == None


def test_get_query_to_df():
    """
    This function tests the get_query_to_df method of the Variants class by running two queries and
    checking their results.
    """

    # Init files
    input_vcf = tests_data_folder + "/example.vcf.gz"

    # Create object
    vcf = Variants(input=input_vcf, load=True)

    # Query
    query = "SELECT 1 AS query"
    result_query = vcf.get_query_to_df(query)
    assert result_query["query"][0] == 1

    # Query
    query = "SELECT * FROM variants"
    result_query = vcf.get_query_to_df(query)
    assert len(result_query) == 7


def test_get_query_to_df_no_samples():
    """
    This function tests the get_query_to_df method of the Variants class when there are no samples in
    the input VCF file.
    """

    # Init files
    input_vcf = tests_data_folder + "/example.no_samples.vcf.gz"

    # Create object
    variants = Variants(input=input_vcf, load=True)

    # Query
    query = "SELECT 1 AS query"
    result_query = variants.get_query_to_df(query)
    assert result_query["query"][0] == 1

    # Query
    query = "SELECT * FROM variants"
    result_query = variants.get_query_to_df(query)
    assert len(result_query) == 10


def test_get_operations_help():
    """
    The function `test_get_operations_help` tests the `get_operations_help` method of the `Variants`
    class.
    """

    expected_operations_help = ["Available calculation operations:"]

    variants = Variants()
    actual_operations_help = variants.get_operations_help()

    print(actual_operations_help[0])

    assert expected_operations_help[0] == actual_operations_help[0]


def test_rename_fields():
    """
    The function `test_rename_fields` renames specified fields in a VCF file and checks if the output
    VCF is in the correct format using pyVCF.
    """

    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # Init files
        input_vcf = tests_data_folder + "/example.annotation_names.vcf"
        output_vcf = f"{tmp_dir}/output.vcf"

        # Create object
        variants = Variants(input=input_vcf, load=True)

        # Fieldst to rename
        fields_to_rename = {
            "CLNSIG": "CLNSIG_renamed",
            "PREFIXCLNSIG": "PREFIXCLNSIG_renamed",
            "DP": "depth",
            "field_not_in_header": "field_not_in_header_renamed",
            "": "",
            "SIFT": None,
            "SPiP_Alt": "SPiP_alternative",
            "SPiP_alternative": None,
        }

        # Rename fields
        fields_renamed = variants.rename_info_fields(fields_to_rename=fields_to_rename)
        assert fields_renamed == {'CLNSIG': 'CLNSIG_renamed', 'PREFIXCLNSIG': 'PREFIXCLNSIG_renamed', 'DP': 'depth', 'SIFT': None, 'SPiP_Alt': 'SPiP_alternative', 'SPiP_alternative': None}
        assert len(variants.get_query_to_df("SELECT INFO FROM variants WHERE INFO LIKE '%SIFT%'")) == 0
        assert len(variants.get_query_to_df("SELECT INFO FROM variants WHERE INFO LIKE '%None=%'")) == 0
        assert len(variants.get_query_to_df("SELECT INFO FROM variants WHERE INFO LIKE '%SPiP%'")) == 0

        # Check if VCF is in correct format with pyVCF
        remove_if_exists([output_vcf])
        variants.export_output(output_file=output_vcf)
        try:
            vcf_obj = vcf.Reader(filename=output_vcf)
        except:
            assert False


def test_rename_fields_to_param_and_export():
    """
    The function `test_rename_fields_to_param_and_export` renames specified fields in a VCF file and
    checks if the output VCF is in the correct format using pyVCF.
    """

    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # Init files
        input_vcf = tests_data_folder + "/example.annotation_names.vcf"
        output_vcf = f"{tmp_dir}/output.test.vcf"

        # Param
        param = {
            "export": {
                "fields_to_rename": {
                    "CLNSIG": "CLNSIG_renamed",
                    "PREFIXCLNSIG": "PREFIXCLNSIG_renamed",
                    "DP": "depth",
                    "field_not_in_header": "field_not_in_header_renamed",
                    "": "",
                    "SIFT": None,
                    "SPiP_Alt": "SPiP_alternative",
                    "SPiP_alternative": None,
                }
            }
        }

        # Create object
        variants = Variants(input=input_vcf, param=param, load=True)

        # Check if VCF is in correct format with pyVCF
        remove_if_exists([output_vcf])
        variants.export_output(output_file=output_vcf)
        try:
            vcf_obj = vcf.Reader(filename=output_vcf)
            assert list(set(vcf_obj.infos.keys())).sort() == ['CLNSIGSUFFIX', 'AA', 'NS', 'PREFIXCLNSIG_renamed', 'CLNSIG_renamed', 'depth'].sort()
        except:
            assert False
