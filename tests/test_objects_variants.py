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


def test_get_explode_infos_fields():
    """
    The function `test_get_explode_infos_fields()` tests the `get_explode_infos_fields()` method of the
    `Variants` class.
    """

    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # Init files
        input_vcf = tests_data_folder + "/example.vcf.gz"
        output_tsv = f"{tmp_dir}/example.tsv"
        fields_all = ["NS", "DP", "AA", "CLNSIG", "SIFT"]
        # ['AA', 'CLNSIG', 'DP', 'NS', 'SIFT'] # Sorted all fields list

        # remove if exists
        remove_if_exists([output_tsv])

        # Create object
        variants = Variants(input=input_vcf, output=output_tsv, load=True)

        # 1 field
        explode_infos_fields = "SIFT"
        explode_infos_fields_list = variants.get_explode_infos_fields(
            explode_infos_fields
        )
        assert explode_infos_fields_list == ["SIFT"]

        # multiple field
        explode_infos_fields = "SIFT,DP,AD"
        explode_infos_fields_list = variants.get_explode_infos_fields(
            explode_infos_fields
        )
        assert explode_infos_fields_list == ["SIFT", "DP", "AD"]

        # * fields
        explode_infos_fields = "*"
        explode_infos_fields_list = variants.get_explode_infos_fields(
            explode_infos_fields
        )
        assert explode_infos_fields_list == sorted(fields_all)

        # 1 field at beginning with all fields
        explode_infos_fields = "DP,*"
        explode_infos_fields_list = variants.get_explode_infos_fields(
            explode_infos_fields
        )
        assert explode_infos_fields_list == ["DP", "AA", "CLNSIG", "NS", "SIFT"]

        # 1 field at the end with all fields
        explode_infos_fields = "*,DP"
        explode_infos_fields_list = variants.get_explode_infos_fields(
            explode_infos_fields
        )
        assert explode_infos_fields_list == ["AA", "CLNSIG", "NS", "SIFT", "DP"]

        # all fields between 2 fields
        explode_infos_fields = "SIFT,*,DP"
        explode_infos_fields_list = variants.get_explode_infos_fields(
            explode_infos_fields
        )
        assert explode_infos_fields_list == ["SIFT", "AA", "CLNSIG", "NS", "DP"]

        # 2 same fields
        explode_infos_fields = "DP,DP"
        explode_infos_fields_list = variants.get_explode_infos_fields(
            explode_infos_fields
        )
        assert explode_infos_fields_list == ["DP"]

        # 2 same fields in a list
        explode_infos_fields = "DP,SIFT,DP,AA"
        explode_infos_fields_list = variants.get_explode_infos_fields(
            explode_infos_fields
        )
        assert explode_infos_fields_list == ["DP", "SIFT", "AA"]

        # 2 same fields in a list conatining *
        explode_infos_fields = "DP,*,DP,AA"
        explode_infos_fields_list = variants.get_explode_infos_fields(
            explode_infos_fields
        )
        assert explode_infos_fields_list == ["DP", "CLNSIG", "NS", "SIFT", "AA"]

        # 1 field not in header and keep it
        explode_infos_fields = "NOT_field"
        explode_infos_fields_list = variants.get_explode_infos_fields(
            explode_infos_fields, remove_fields_not_in_header=False
        )
        assert explode_infos_fields_list == ["NOT_field"]

        # 1 field not in header and keep it, with another field
        explode_infos_fields = "NOT_field,DP"
        explode_infos_fields_list = variants.get_explode_infos_fields(
            explode_infos_fields, remove_fields_not_in_header=False
        )
        assert explode_infos_fields_list == ["NOT_field", "DP"]

        # 1 field not in header and remove it
        explode_infos_fields = "NOT_field"
        explode_infos_fields_list = variants.get_explode_infos_fields(
            explode_infos_fields, remove_fields_not_in_header=True
        )
        assert explode_infos_fields_list == []

        # 1 field not in header and remove it, with another field
        explode_infos_fields = "NOT_field,DP"
        explode_infos_fields_list = variants.get_explode_infos_fields(
            explode_infos_fields, remove_fields_not_in_header=True
        )
        assert explode_infos_fields_list == ["DP"]

        # 1 field not in header and remove it, with ALL fields
        explode_infos_fields = "NOT_field,*"
        explode_infos_fields_list = variants.get_explode_infos_fields(
            explode_infos_fields, remove_fields_not_in_header=True
        )
        assert explode_infos_fields_list == ["AA", "CLNSIG", "DP", "NS", "SIFT"]

        # multiple fields with spaces
        explode_infos_fields = "AA, DP ,SIFT,NS "
        explode_infos_fields_list = variants.get_explode_infos_fields(
            explode_infos_fields, remove_fields_not_in_header=True
        )
        assert explode_infos_fields_list == ["AA", "DP", "SIFT", "NS"]

        # 1 field with pattern
        explode_infos_fields = ".*S.*"
        explode_infos_fields_list = variants.get_explode_infos_fields(
            explode_infos_fields, remove_fields_not_in_header=True
        )
        assert explode_infos_fields_list == ["CLNSIG", "NS", "SIFT"]

        # 1 field with pattern, with one field in pattern at the begenning
        explode_infos_fields = "NS,.*S.*"
        explode_infos_fields_list = variants.get_explode_infos_fields(
            explode_infos_fields, remove_fields_not_in_header=True
        )
        assert explode_infos_fields_list == ["NS", "CLNSIG", "SIFT"]

        # 1 field with pattern, with one field in pattern at the end
        explode_infos_fields = ".*S.*,NS"
        explode_infos_fields_list = variants.get_explode_infos_fields(
            explode_infos_fields, remove_fields_not_in_header=True
        )
        assert explode_infos_fields_list == ["CLNSIG", "SIFT", "NS"]

        # 1 field with pattern, with one field in pattern at the begenning, and all at the end
        explode_infos_fields = "NS,.*S.*,*"
        explode_infos_fields_list = variants.get_explode_infos_fields(
            explode_infos_fields, remove_fields_not_in_header=True
        )
        assert explode_infos_fields_list == ["NS", "CLNSIG", "SIFT", "AA", "DP"]

        # 1 field with pattern, with one field in pattern in the middle, and all at the end
        explode_infos_fields = ".*S.*,NS,*"
        explode_infos_fields_list = variants.get_explode_infos_fields(
            explode_infos_fields, remove_fields_not_in_header=True
        )
        assert explode_infos_fields_list == ["CLNSIG", "SIFT", "NS", "AA", "DP"]

        # 1 field with pattern, all at the middle, with one field in pattern at the end
        explode_infos_fields = ".*S.*,*,NS"
        explode_infos_fields_list = variants.get_explode_infos_fields(
            explode_infos_fields, remove_fields_not_in_header=True
        )
        assert explode_infos_fields_list == ["CLNSIG", "SIFT", "AA", "DP", "NS"]

        # all at the middle, 1 field with pattern (not considered because of all before), with one field in pattern at the end
        explode_infos_fields = "*,.*S.*,NS"
        explode_infos_fields_list = variants.get_explode_infos_fields(
            explode_infos_fields, remove_fields_not_in_header=True
        )
        assert explode_infos_fields_list == ["AA", "CLNSIG", "DP", "SIFT", "NS"]


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


def test_export():
    """
    The function tests the export functionality of a database for various input and output formats.
    """

    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # database input/format
        for database_input_index in [
            "parquet",
            "partition_parquet",
            "vcf",
            "vcf_gz",
            "tsv",
            "csv",
            "example_vcf",
        ]:
            for database_output_format in [
                "parquet",
                "partition_parquet",
                "vcf",
                "vcf.gz",
                "tsv",
                "csv",
                "json",
                "bed",
            ]:
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


def test_set_get_header():
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


def test_set_get_header_no_samples():
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


def test_set_get_header_in_config():
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
    # with pytest.raises(ValueError) as e:
    #     variants = Variants(input=input_vcf)
    # assert str(e.value) == f"No header for file {input_vcf}"


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


# def test_load_psv():
#     """
#     This function tests if a PSV file can be loaded into a Variants object and if the expected number of
#     variants is present in the database.
#     """

#     # Init files
#     input_vcf = tests_data_folder + "/example.psv"

#     # Create object
#     variants = Variants(input=input_vcf, load=True)

#     # Check data loaded
#     result = variants.get_query_to_df("SELECT count(*) AS count FROM variants")
#     nb_variant_in_database = result["count"][0]

#     expected_number_of_variants = 7

#     assert nb_variant_in_database == expected_number_of_variants


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


###
### Prioritization
###


def test_prioritization():
    """
    This is a test function for prioritization of variants in a VCF file using a specified configuration
    and parameter dictionary.
    """

    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # Init files
        input_vcf = tests_data_folder + "/example.vcf.gz"
        output_vcf = f"{tmp_dir}/output.vcf.gz"

        # Construct config dict
        config = {}

        # Construct param dict
        param = {
            "prioritization": {
                "prioritization_config": tests_data_folder
                + "/prioritization_profiles.json",
                "profiles": ["default", "GERMLINE"],
                "pzfields": ["PZFlag", "PZScore", "PZComment", "PZInfos"],
            }
        }

        # Create object
        variants = Variants(
            input=input_vcf, output=output_vcf, load=True, config=config, param=param
        )

        # Prioritization
        variants.prioritization()

        # Check all priorized default profile
        result = variants.get_query_to_df(
            """
            SELECT * FROM variants
            WHERE INFO LIKE '%PZFlag_default=%'
            AND INFO LIKE '%PZScore_default=%'
            AND INFO LIKE '%PZComment_default=%'
            AND INFO LIKE '%PZInfos_default=%'
            """
        )
        assert len(result) == 4

        # Check all priorized GERMLINE profile
        result = variants.get_query_to_df(
            """
            SELECT * FROM variants
            WHERE INFO LIKE '%PZFlag_GERMLINE=%'
            AND INFO LIKE '%PZScore_GERMLINE=%'
            AND INFO LIKE '%PZComment_default=%'
            AND INFO LIKE '%PZInfos_default=%'
            """
        )
        assert len(result) == 4

        # Check all priorized default profile (as default)
        result = variants.get_query_to_df(
            """
            SELECT * FROM variants
            WHERE INFO LIKE '%PZFlag=%'
            AND INFO LIKE '%PZScore=%'
            AND INFO LIKE '%PZComment_default=%'
            AND INFO LIKE '%PZInfos_default=%'
            """
        )
        assert len(result) == 4

        # Check all priorized default profile
        result = variants.get_query_to_df(
            """
            SELECT * FROM variants
            WHERE INFO LIKE '%PZFlag_default=%'
            AND INFO LIKE '%PZScore_default=%'
            """
        )
        assert len(result) == 7

        # Check all priorized GERMLINE profile
        result = variants.get_query_to_df(
            """
            SELECT * FROM variants
            WHERE INFO LIKE '%PZFlag_GERMLINE=%'
            AND INFO LIKE '%PZScore_GERMLINE=%'
            """
        )
        assert len(result) == 7

        # Check all priorized default profile (as default)
        result = variants.get_query_to_df(
            """
            SELECT * FROM variants
            WHERE INFO LIKE '%PZFlag=%'
            AND INFO LIKE '%PZScore=%'
            """
        )
        assert len(result) == 7

        # Check annotation1
        result = variants.get_query_to_df(
            """ SELECT * FROM variants WHERE "#CHROM" = 'chr1' AND POS = 28736 AND REF = 'A' AND ALT = 'C' AND INFO LIKE '%PZScore_default=15%' """
        )
        assert len(result) == 1

        # Check FILTERED
        result = variants.get_query_to_df(
            f""" SELECT INFO FROM variants WHERE INFO LIKE '%FILTERED%' """
        )
        assert len(result) == 1

        # Check if VCF is in correct format with pyVCF
        remove_if_exists([output_vcf])
        variants.export_output()
        try:
            vcf.Reader(filename=output_vcf)
        except:
            assert False


def test_prioritization_full_unsorted():
    """
    This is a test function for prioritization of variants in a VCF file using a specified configuration
    and parameter dictionary.
    Test with a VCF full variants type: SNV, INDEL, MNV, SV
    This VCF is unsorted
    """

    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # Init files
        input_vcf = tests_data_folder + "/example.full.unsorted.vcf.gz"
        output_vcf = f"{tmp_dir}/output.vcf.gz"

        # Construct config dict
        config = {}

        # Construct param dict
        param = {
            "prioritization": {
                "prioritization_config": tests_data_folder
                + "/prioritization_profiles.json",
                "profiles": ["default", "GERMLINE"],
                "pzfields": ["PZFlag", "PZScore", "PZComment", "PZInfos"],
            }
        }

        # Create object
        variants = Variants(
            input=input_vcf, output=output_vcf, load=True, config=config, param=param
        )

        # Prioritization
        variants.prioritization()

        # Check all priorized default profile
        result = variants.get_query_to_df(
            """
            SELECT * FROM variants
            WHERE INFO LIKE '%PZFlag_default=%'
            AND INFO LIKE '%PZScore_default=%'
            AND INFO LIKE '%PZComment_default=%'
            AND INFO LIKE '%PZInfos_default=%'
            """
        )
        assert len(result) == 4

        # Check all priorized GERMLINE profile
        result = variants.get_query_to_df(
            """
            SELECT * FROM variants
            WHERE INFO LIKE '%PZFlag_GERMLINE=%'
            AND INFO LIKE '%PZScore_GERMLINE=%'
            AND INFO LIKE '%PZComment_GERMLINE=%'
            AND INFO LIKE '%PZInfos_GERMLINE=%'
            """
        )
        assert len(result) == 2

        # Check all priorized default profile (as default)
        result = variants.get_query_to_df(
            """
            SELECT * FROM variants
            WHERE INFO LIKE '%PZFlag=%'
            AND INFO LIKE '%PZScore=%'
            AND INFO LIKE '%PZComment=%'
            AND INFO LIKE '%PZInfos=%'
            """
        )
        assert len(result) == 4

        # Check all priorized default profile
        result = variants.get_query_to_df(
            """
            SELECT * FROM variants
            WHERE INFO LIKE '%PZFlag_default=%'
            AND INFO LIKE '%PZScore_default=%'
            """
        )
        assert len(result) == 36

        # Check all priorized GERMLINE profile
        result = variants.get_query_to_df(
            """
            SELECT * FROM variants
            WHERE INFO LIKE '%PZFlag_GERMLINE=%'
            AND INFO LIKE '%PZScore_GERMLINE=%'
            """
        )
        assert len(result) == 36

        # Check all priorized default profile (as default)
        result = variants.get_query_to_df(
            """
            SELECT * FROM variants
            WHERE INFO LIKE '%PZFlag=%'
            AND INFO LIKE '%PZScore=%'
            """
        )
        assert len(result) == 36

        # Check annotation1
        result = variants.get_query_to_df(
            """ SELECT * FROM variants WHERE "#CHROM" = 'chr1' AND POS = 28736 AND REF = 'A' AND ALT = 'C' AND INFO LIKE '%PZScore_default=15%' """
        )
        assert len(result) == 1

        # Check FILTERED
        result = variants.get_query_to_df(
            f""" SELECT INFO FROM variants WHERE INFO LIKE '%FILTERED%' """
        )
        assert len(result) == 1

        # Check if VCF is in correct format with pyVCF
        remove_if_exists([output_vcf])
        variants.export_output()
        try:
            vcf.Reader(filename=output_vcf)
        except:
            assert False


def test_prioritization_varank():
    """
    This is a test function for the prioritization feature of a Python package called "Variants".
    """

    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # Init files
        input_vcf = tests_data_folder + "/example.vcf.gz"
        output_vcf = f"{tmp_dir}/output.vcf.gz"

        # Construct config dict
        config = {}

        # Construct param dict
        param = {
            "prioritization": {
                "prioritization_config": tests_data_folder
                + "/prioritization_profiles.json",
                "profiles": ["default", "GERMLINE"],
                "pzfields": ["PZFlag", "PZScore", "PZComment", "PZInfos"],
                "prioritization_score_mode": "VaRank",
            }
        }

        # Create object
        variants = Variants(
            input=input_vcf, output=output_vcf, load=True, config=config, param=param
        )

        # Prioritization
        variants.prioritization()

        # Check all priorized
        result = variants.get_query_to_df(""" SELECT INFO FROM variants """)
        assert len(result) > 0

        # Check annotation1
        result = variants.get_query_to_df(
            """ SELECT 1 AS count FROM variants WHERE "#CHROM" = 'chr1' AND POS = 28736 AND REF = 'A' AND ALT = 'C' AND INFO LIKE '%PZScore_default=15%' """
        )
        assert len(result) == 1

        # Check if VCF is in correct format with pyVCF
        remove_if_exists([output_vcf])
        variants.export_output()
        try:
            vcf.Reader(filename=output_vcf)
        except:
            assert False


def test_prioritization_all_profiles():
    """
    This function tests if an error is raised when there are no prioritization configuration profiles
    provided.
    """

    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # Init files
        input_vcf = tests_data_folder + "/example.vcf.gz"
        output_vcf = f"{tmp_dir}/output.vcf.gz"

        # Construct config dict
        config = {}

        # Construct param dict
        param = {
            "prioritization": {
                "prioritization_config": tests_data_folder
                + "/prioritization_profiles.json"
            }
        }

        # Create object
        variants = Variants(
            input=input_vcf, output=output_vcf, load=True, config=config, param=param
        )

        # Check if VCF is in correct format with pyVCF
        remove_if_exists([output_vcf])
        variants.export_output()
        try:
            vcf.Reader(filename=output_vcf)
        except:
            assert False

        # # Prioritization fail
        # with pytest.raises(ValueError) as e:
        #     variants.prioritization()
        # assert str(e.value) == f"NO Profiles configuration"


def test_prioritization_profile_not_configured():
    """
    This function tests if an error is raised when there are no prioritization configuration profiles
    provided.
    """

    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # Init files
        input_vcf = tests_data_folder + "/example.vcf.gz"
        output_vcf = f"{tmp_dir}/output.vcf.gz"
        profile_not_defined = "profile_not_defined"

        # Construct config dict
        config = {}

        # Construct param dict
        param = {
            "prioritization": {
                # "prioritization_config": tests_data_folder + "/prioritization_profiles.json"
                "profiles": [profile_not_defined]
            }
        }

        # Create object
        variants = Variants(
            input=input_vcf, output=output_vcf, load=True, config=config, param=param
        )

        # Prioritization fail
        with pytest.raises(ValueError) as e:
            variants.prioritization()
        assert str(e.value) == f"Profile '{profile_not_defined}' NOT configured"


def test_prioritization_no_pzfields():
    """
    This function tests the prioritization method of a Variants object when there are no pzfields
    specified.
    """

    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # Init files
        input_vcf = tests_data_folder + "/example.vcf.gz"
        output_vcf = f"{tmp_dir}/output.vcf.gz"

        # Construct config dict
        config = {}

        # Construct param dict
        param = {
            "prioritization": {
                "prioritization_config": tests_data_folder
                + "/prioritization_profiles.json",
                "profiles": [],
                "pzfields": [],
            }
        }

        # Create object
        variants = Variants(
            input=input_vcf, output=output_vcf, load=True, config=config, param=param
        )

        # Prioritization
        variants.prioritization()

        # Check all priorized
        result = variants.get_query_to_df(
            """ SELECT * FROM variants WHERE INFO LIKE '%PZScore_default=%' """
        )
        assert len(result) == 0

        # Check if VCF is in correct format with pyVCF
        remove_if_exists([output_vcf])
        variants.export_output()
        try:
            vcf.Reader(filename=output_vcf)
        except:
            assert False


def test_prioritization_no_infos():
    """
    This function tests the prioritization method of the Variants class when there is no information
    available.
    """

    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # Init files
        input_vcf = tests_data_folder + "/example.no_samples.vcf.gz"
        output_vcf = f"{tmp_dir}/output.vcf.gz"

        # Construct config dict
        config = {}

        # Construct param dict
        param = {
            "prioritization": {
                "prioritization_config": tests_data_folder
                + "/prioritization_profiles.json",
                "profiles": ["default", "GERMLINE"],
                "pzfields": ["PZFlag", "PZScore", "PZComment", "PZInfos"],
            }
        }

        # Create object
        variants = Variants(
            input=input_vcf, output=output_vcf, load=True, config=config, param=param
        )

        # Prioritization
        variants.prioritization()

        result = variants.get_query_to_df(
            """
            SELECT INFO FROM variants WHERE INFO LIKE '%PZScore=0%' AND INFO LIKE '%PZFlag=PASS%'
            """
        )
        assert len(result) == 10

        # Check if VCF is in correct format with pyVCF
        remove_if_exists([output_vcf])
        variants.export_output()
        try:
            vcf.Reader(filename=output_vcf)
        except:
            assert False


###
### Calculation
###


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


def test_calculation_snpeff_hgvs():
    """
    This is a test function for the calculation of snpeff_hgvs in a VCF file using the Variants class.
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
            """ SELECT * FROM variants WHERE INFO LIKE '%NOMEN=%' """
        )
        assert len(result) == 7

        # Check transcript priority
        result = variants.get_query_to_df(
            """ SELECT * FROM variants WHERE INFO LIKE '%NOMEN=EGFR:NM_001346897%' """
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


def test_calculation_trio():
    """
    This is a test function for the calculation of trio variants in a VCF file using specific
    parameters.
    """

    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # Init files
        input_vcf = tests_data_folder + "/example.vcf.gz"
        output_vcf = f"{tmp_dir}/output.vcf.gz"

        # Construct param dict
        param = {
            "calculation": {
                "calculations": {
                    "trio": {
                        "father": "sample1",
                        "mother": "sample2",
                        "child": "sample3",
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
