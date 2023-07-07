# -*- coding: utf-8 -*-
"""
Tests

Usage:
pytest tests/

Coverage:
coverage run -m pytest . -x -v
coverage report --include=howard/* -m
"""

import logging as log
import os
import sys
import duckdb
import re
import Bio.bgzf as bgzf
import gzip
import pytest

from howard.commons import *
from howard.objects.variants import Variants
from howard.tools.databases import *


# Main tests folder
tests_folder = os.path.dirname(__file__)
tests_data_folder = tests_folder + "/data"
tests_annotations_folder = tests_folder + "/data/annotations"

# Tools folder
tests_tools = "/tools"

# Test config
tests_config = {
  "threads": 1,
  "memory": None,
  "verbosity": "warning",
  "folders": {
    "databases": {
      "root": "",
      "parquet": [f"{tests_folder}/data/annotations"],
      "bcftools": [f"{tests_folder}/data/annotations"],
      "annovar": f"{tests_folder}/data/annotations/annovar",
      "snpeff": f"{tests_folder}/data/annotations/snpeff",
      "varank": f"{tests_folder}/data/annotations/varank"
    }
  },
  "tools": {
    "bcftools": {"bin": "bcftools"},
    "bgzip": {"bin": "bgzip"},
    "snpeff": {"jar": f"{tests_tools}/snpeff/current/bin/snpEff.jar"},
    "java": {"bin": "/usr/bin/java"},
    "annovar": {"bin": f"{tests_tools}/annovar/current/bin/table_annovar.pl"}
  }
}

# Annotation databases
database_files = {
    "parquet" : tests_annotations_folder + "/nci60.parquet",
    "parquet_without_header" : tests_annotations_folder + "/nci60.without_header.parquet",
    "duckdb" : tests_annotations_folder + "/nci60.duckdb",
    "duckdb_no_annotation_table" : tests_annotations_folder + "/nci60.no_annotation_table.duckdb",
    "sqlite" : tests_annotations_folder + "/nci60.sqlite",
    "vcf" : tests_annotations_folder + "/nci60.vcf",
    "vcf_gz" : tests_annotations_folder + "/nci60.vcf.gz",
    "vcf_without_header" : tests_annotations_folder + "/nci60.without_header.vcf",
    "vcf_gz_without_header" : tests_annotations_folder + "/nci60.without_header.vcf.gz",
    "tsv" : tests_annotations_folder + "/nci60.tsv",
    "tsv_alternative_columns" : tests_annotations_folder + "/nci60.alternative_columns.tsv",
    "tsv_failed_columns" : tests_annotations_folder + "/nci60.failed_columns.tsv",
    "tsv_lower_columns" : tests_annotations_folder + "/nci60.lower_columns.tsv",
    "tsv_without_header" : tests_annotations_folder + "/nci60.without_header.tsv",
    "tsv_variants" : tests_annotations_folder + "/nci60.variants.tsv",
    "tsv_gz" : tests_annotations_folder + "/nci60.tsv.gz",
    "csv" : tests_annotations_folder + "/nci60.csv",
    "csv_gz" : tests_annotations_folder + "/nci60.csv.gz",
    "tbl" : tests_annotations_folder + "/nci60.tbl",
    "tbl_gz" : tests_annotations_folder + "/nci60.tbl.gz",
    "json" : tests_annotations_folder + "/nci60.json",
    "json_gz" : tests_annotations_folder + "/nci60.json.gz",
    "bed" : tests_annotations_folder + "/annotation_regions.bed",
    "bed_gz" : tests_annotations_folder + "/annotation_regions.bed.gz",
    "example_vcf" : tests_data_folder + "/example.vcf",

}




def test_export_query():
    """
    This is a test function for exporting data from a VCF file to a TSV file using SQL queries.
    """

    # Init files
    input_vcf = tests_folder + "/data/example.vcf.gz"
    output_tsv = "/tmp/example.tsv"

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

    # database input/format
    for database_input_index in ["parquet", "vcf", "vcf_gz", "tsv", "csv", "tsv_alternative_columns", "example_vcf"]:
        for database_output_format in ["parquet", "vcf", "vcf.gz", "tsv", "csv", "json", "bed"]:
            input_database = database_files.get(database_input_index)
            output_database=f"/tmp/output_database.{database_output_format}"
            output_header=output_database+".hdr"
            variants = Variants(input=input_database, output=output_database, load=True)
            remove_if_exists([output_database,output_header])
            try:
                assert variants.export_output(output_file=output_database, output_header=output_header)
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
    input_vcf = tests_folder + "/data/example.vcf.gz"
    new_input_vcf = tests_folder + "/data/example.parquet"

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

    # Init files
    input_vcf = tests_folder + "/data/example.vcf.gz"
    output_vcf = tests_folder + "/data//tmp/example.parquet"
    new_output_vcf = tests_folder + "/data/tmp/example.parquet"

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
    input_vcf = tests_folder + "/data/example.vcf.gz"
    input_config = { "option": "option_value" }
    new_input_config = { "option": "option_value", "new_option": "new_option_value" }

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
    input_vcf = tests_folder + "/data/example.vcf.gz"
    input_param = { "option": "option_value" }
    new_input_param = { "option": "option_value", "new_option": "new_option_value" }

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
    input_vcf = tests_folder + "/data/example.vcf.gz"

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
    header_columns_sample_list_expected = ['sample1', 'sample2', 'sample3', 'sample4']
    assert header_columns_sample_list == header_columns_sample_list_expected


def test_set_get_header_no_samples():
    """
    This function tests various methods related to getting and setting the header of a VCF file when
    there are no samples present.
    """

    # Init files
    input_vcf = tests_folder + "/data/example.no_samples.vcf.gz"

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
    header_columns_as_sql_expected = """ "#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO" """.strip()
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
    input_vcf = tests_folder + "/data/example.parquet"
    input_config = { "header_file":  tests_folder + "/data/example.parquet.hdr" }

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
    header_columns_sample_list_expected = ['sample1', 'sample2', 'sample3', 'sample4']
    assert header_columns_sample_list == header_columns_sample_list_expected


def test_load_without_header():
    """
    This function tests if a ValueError is raised when attempting to load a Parquet file without a
    header using the Variants object.
    """

    # Init files
    input_vcf = tests_folder + "/data/example.without_header.parquet"

    # Create object
    with pytest.raises(ValueError) as e:
        variants = Variants(input=input_vcf)
    assert str(e.value) == f"No header for file {input_vcf}"


def test_read_vcf_header():
    """
    This function tests the read_vcf_header method of the Variants class by checking if the header list
    is not empty and has a length of 53.
    """

    # Init files
    input_vcf = tests_folder + "/data/example.parquet"
    vcf_header = tests_folder + "/data/example.parquet.hdr"

    # Create connection
    conn = duckdb.connect(":memory:")

    # Create object
    variants = Variants(conn=conn, input=input_vcf)

    # Check read_vcf_header
    with open(vcf_header, 'rt') as f:
        header_list = variants.read_vcf_header(f)
    assert header_list != []
    assert len(header_list) == 53


def test_load_when_init():
    """
    This function tests if the Variants object loads data correctly from a VCF file during
    initialization.
    """

    # Init files
    input_vcf = tests_folder + "/data/example.vcf.gz"

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
    input_config = { "header_file":  tests_folder + "/data/example.parquet.hdr" }

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
    input_vcf = tests_folder + "/data/example.vcf.gz"

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
    input_vcf = tests_folder + "/data/example.full.unsorted.vcf.gz"

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
    input_vcf = tests_folder + "/data/example.no_samples.vcf.gz"

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
    input_vcf = tests_folder + "/data/example.parquet"

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
    input_vcf = tests_folder + "/data/example.vcf"

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
    input_vcf = tests_folder + "/data/example.csv"

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
    input_vcf = tests_folder + "/data/example.tsv"

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
#     input_vcf = tests_folder + "/data/example.psv"

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

    # Create duckdb database

    # Init files
    input_vcf = tests_folder + "/data/example.vcf.gz"
    output_duckdb = "/tmp/example.duckdb"

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
    input_vcf = tests_folder + "/data/example.vcf.gz"
    input_config = { "connexion_type": "memory" }

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
    input_vcf = tests_folder + "/data/example.vcf.gz"
    input_config = { "connexion_type": "tmpfile" }

    # Create object
    variants = Variants(input=input_vcf, config=input_config)

    # get connexion_db
    connexion_db = variants.get_connexion_db()

    assert os.path.exists(connexion_db)


def test_get_connexion_db_file():
    """
    This function tests the "get_connexion_db" method of a "Variants" object in Python.
    """

    # Init files
    input_vcf = tests_folder + "/data/example.vcf.gz"
    input_config = { "connexion_type": "/tmp/connexion.duckdb" }

    # Remove if exists
    remove_if_exists(["/tmp/connexion.duckdb"])

    # Create object
    variants = Variants(input=input_vcf, config=input_config)

    # get connexion_db
    connexion_db = variants.get_connexion_db()

    assert connexion_db == "/tmp/connexion.duckdb"


def test_get_table_variants():
    """
    This function tests the get_table_variants method of the Variants class in Python.
    """

    # Init files
    input_vcf = tests_folder + "/data/example.vcf.gz"

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
    input_vcf = tests_folder + "/data/example.vcf.gz"

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
    input_vcf = tests_folder + "/data/example.vcf.gz"
    input_config = { "connexion_format": "sqlite" }

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
    input_vcf = tests_folder + "/data/example.vcf.gz"

    # Create object
    variants = Variants(input=input_vcf)

    # get connexion_db
    verbose = variants.get_verbose()
    assert not verbose

    # config verbose True
    input_config = { "verbose": True }

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
    input_vcf = tests_folder + "/data/example.vcf.gz"
    input_config = { "connexion_type": "memory" }

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
    input_vcf = tests_folder + "/data/example.vcf.gz"
    input_config = { "connexion_type": "tmpfile" }

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

    # Init files
    input_vcf = tests_folder + "/data/example.vcf.gz"
    input_config = { "connexion_type": "/tmp/output.duckdb" }

    remove_if_exists(["/tmp/output.duckdb"])

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
    input_vcf = tests_folder + "/data/example.vcf.gz"
    input_config = { "connexion_format": "sqlite" }

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

    # Init files
    input_vcf = tests_folder + "/data/example.vcf.gz"
    output_vcf = "/tmp/example.vcf.gz"

    # remove if exists
    remove_if_exists([output_vcf])

    # Create object
    variants = Variants(input=input_vcf, output=output_vcf, load=True)

    # Check get_output
    variants.export_output()
    assert os.path.exists(output_vcf)

    # Check get_output without header
    remove_if_exists([output_vcf])
    variants.export_output(export_header=False)
    assert os.path.exists(output_vcf) and os.path.exists(output_vcf + ".hdr")

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

    # Init files
    input_vcf = tests_folder + "/data/example.full.unsorted.vcf.gz"
    output_vcf = "/tmp/example.vcf.gz"

    # remove if exists
    remove_if_exists([output_vcf])

    # Create object
    variants = Variants(input=input_vcf, output=output_vcf, load=True)

    # Check get_output
    variants.export_output()
    assert os.path.exists(output_vcf)

    # Check get_output without header
    remove_if_exists([output_vcf])
    variants.export_output(export_header=False)
    assert os.path.exists(output_vcf) and os.path.exists(output_vcf + ".hdr")

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

    # Init files
    input_vcf = tests_folder + "/data/example.vcf.gz"
    output_vcf = "/tmp/example.vcf"

    # remove if exists
    remove_if_exists([output_vcf])

    # Create object
    variants = Variants(input=input_vcf, output=output_vcf, load=True)

    # Check get_output
    variants.export_output()
    assert os.path.exists(output_vcf)

    # Check get_output without header
    remove_if_exists([output_vcf])
    variants.export_output(export_header=False)
    assert os.path.exists(output_vcf) and os.path.exists(output_vcf + ".hdr")

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

    # Init files
    input_vcf = tests_folder + "/data/example.vcf.gz"
    output_vcf = "/tmp/example.parquet"

    # remove if exists
    remove_if_exists([output_vcf])

    # Create object
    variants = Variants(input=input_vcf, output=output_vcf, load=True)

    # Check get_output
    variants.export_output()
    assert os.path.exists(output_vcf)

    # Check get_output without header
    remove_if_exists([output_vcf])
    variants.export_output(export_header=False)
    assert os.path.exists(output_vcf) and os.path.exists(output_vcf + ".hdr")


def test_export_output_duckdb():
    """
    This function tests the export_output method of the Variants class in Python's DuckDB library.
    """

    # Init files
    input_vcf = tests_folder + "/data/example.vcf.gz"
    output_duckdb = "/tmp/example.duckdb"

    # remove if exists
    remove_if_exists([output_duckdb])

    # Create object
    variants = Variants(input=input_vcf, output=output_duckdb, load=True)

    # Check get_output
    variants.export_output()
    assert os.path.exists(output_duckdb)

    # remove if exists
    remove_if_exists([output_duckdb])

    # Check get_output without header
    variants.export_output(export_header=False)
    assert os.path.exists(output_duckdb) and os.path.exists(output_duckdb + ".hdr")


def test_export_output_tsv():
    """
    This function tests the export_output method of the Variants class in Python, which exports a VCF
    file to a TSV file.
    """

    # Init files
    input_vcf = tests_folder + "/data/example.vcf.gz"
    output_vcf = "/tmp/example.tsv"

    # remove if exists
    remove_if_exists([output_vcf])

    # Create object
    variants = Variants(input=input_vcf, output=output_vcf, load=True)

    # Check get_output
    variants.export_output()
    assert os.path.exists(output_vcf)

    # Check get_output without header
    remove_if_exists([output_vcf])
    variants.export_output(export_header=False)
    assert os.path.exists(output_vcf) and os.path.exists(output_vcf + ".hdr")


def test_export_output_tsv_gz():
    """
    This function tests the export_output method of the Variants class in Python, which exports a VCF
    file to a TSV file in gzipped format.
    """

    # Init files
    input_vcf = tests_folder + "/data/example.vcf.gz"
    output_vcf = "/tmp/example.tsv.gz"

    # remove if exists
    remove_if_exists([output_vcf])

    # Create object
    variants = Variants(input=input_vcf, output=output_vcf, load=True)

    # Check get_output
    variants.export_output()
    assert os.path.exists(output_vcf)

    # Check get_output without header
    remove_if_exists([output_vcf])
    variants.export_output(export_header=False)
    assert os.path.exists(output_vcf) and os.path.exists(output_vcf + ".hdr")


def test_export_output_csv():
    """
    This function tests the export_output method of the Variants class in Python by checking if the
    output file exists with and without a header.
    """

    # Init files
    input_vcf = tests_folder + "/data/example.vcf.gz"
    output_vcf = "/tmp/example.csv"

    # remove if exists
    remove_if_exists([output_vcf])

    # Create object
    variants = Variants(input=input_vcf, output=output_vcf, load=True)

    # Check get_output
    variants.export_output()
    assert os.path.exists(output_vcf)

    # Check get_output without header
    remove_if_exists([output_vcf])
    variants.export_output(export_header=False)
    assert os.path.exists(output_vcf) and os.path.exists(output_vcf + ".hdr")


def test_export_output_tbl():
    """
    This function tests the export_output method of the Variants class in Python.
    """

    # Init files
    input_vcf = tests_folder + "/data/example.vcf.gz"
    output_vcf = "/tmp/example.tbl"

    # remove if exists
    remove_if_exists([output_vcf])

    # Create object
    variants = Variants(input=input_vcf, output=output_vcf, load=True)

    # Check get_output
    variants.export_output()
    assert os.path.exists(output_vcf)

    # Check get_output without header
    remove_if_exists([output_vcf])
    variants.export_output(export_header=False)
    assert os.path.exists(output_vcf) and os.path.exists(output_vcf + ".hdr")


def test_export_output_tsv_explode_infos():
    """
    This function tests the export of variant information in TSV format with exploded extra information.
    """

    # Init files
    input_vcf = tests_folder + "/data/example.vcf.gz"
    output_vcf = "/tmp/example.tsv"
    param = {
                "export_extra_infos": True
        }

    # remove if exists
    remove_if_exists([output_vcf])

    # Create object
    variants = Variants(input=input_vcf, output=output_vcf, load=True, param=param)

    # Explode infos
    variants.explode_infos()

    # Check get_output
    variants.export_output()
    assert os.path.exists(output_vcf)

    # Check get_output without header
    remove_if_exists([output_vcf])
    variants.export_output(export_header=False)
    assert os.path.exists(output_vcf) and os.path.exists(output_vcf + ".hdr")


def test_export_from_sqlite_output_vcf_gz():
    """
    This function tests the export of a VCF file in gzipped format with the pyVCF library.
    """

    # Init files
    input_vcf = tests_folder + "/data/example.vcf.gz"
    output_vcf = "/tmp/example.vcf.gz"
    input_config = { "connexion_format": "sqlite" }

    # remove if exists
    remove_if_exists([output_vcf])

    # Create object
    variants = Variants(input=input_vcf, output=output_vcf, config=input_config, load=True)

    # Check get_output
    variants.export_output()
    assert os.path.exists(output_vcf)

    # Check get_output without header
    remove_if_exists([output_vcf])
    variants.export_output(export_header=False)
    assert os.path.exists(output_vcf) and os.path.exists(output_vcf + ".hdr")

    # Check if VCF is in correct format with pyVCF
    try:
        vcf.Reader(filename=output_vcf)
    except:
        assert False


def test_export_from_sqlite_output_vcf():
    """
    This function tests the export of a VCF file in gzipped format with the pyVCF library.
    """

    # Init files
    input_vcf = tests_folder + "/data/example.vcf.gz"
    output_vcf = "/tmp/example.vcf"
    input_config = { "connexion_format": "sqlite" }

    # remove if exists
    remove_if_exists([output_vcf])

    # Create object
    variants = Variants(input=input_vcf, output=output_vcf, config=input_config, load=True)

    # Check get_output
    variants.export_output()
    assert os.path.exists(output_vcf)

    # Check get_output without header
    remove_if_exists([output_vcf])
    variants.export_output(export_header=False)
    assert os.path.exists(output_vcf) and os.path.exists(output_vcf + ".hdr")

    # Check if VCF is in correct format with pyVCF
    try:
        vcf.Reader(filename=output_vcf)
    except:
        assert False


def test_export_from_sqlite_output_parquet():
    """
    This function tests the export of a VCF file in gzipped format with the pyVCF library.
    """

    # Init files
    input_vcf = tests_folder + "/data/example.vcf.gz"
    output_vcf = "/tmp/example.parquet"
    input_config = { "connexion_format": "sqlite" }

    # remove if exists
    remove_if_exists([output_vcf])

    # Create object
    variants = Variants(input=input_vcf, output=output_vcf, config=input_config, load=True)

    # Check get_output
    variants.export_output()
    assert os.path.exists(output_vcf)

    # Check get_output without header
    remove_if_exists([output_vcf])
    variants.export_output(export_header=False)
    assert os.path.exists(output_vcf) and os.path.exists(output_vcf + ".hdr")


def test_export_from_sqlite_output_tsv():
    """
    This function tests the export of a VCF file in gzipped format with the pyVCF library.
    """

    # Init files
    input_vcf = tests_folder + "/data/example.vcf.gz"
    output_vcf = "/tmp/example.tsv"
    input_config = { "connexion_format": "sqlite" }

    # remove if exists
    remove_if_exists([output_vcf])

    # Create object
    variants = Variants(input=input_vcf, output=output_vcf, config=input_config, load=True)

    # Check get_output
    variants.export_output()
    assert os.path.exists(output_vcf)

    # Check get_output without header
    remove_if_exists([output_vcf])
    variants.export_output(export_header=False)
    assert os.path.exists(output_vcf) and os.path.exists(output_vcf + ".hdr")

    # Check if VCF is in correct format with pyVCF
    try:
        vcf.Reader(filename=output_vcf)
    except:
        assert False


def test_export_from_sqlite_output_tsv_gz():
    """
    This function tests the export of a VCF file in gzipped format with the pyVCF library.
    """

    # Init files
    input_vcf = tests_folder + "/data/example.vcf.gz"
    output_vcf = "/tmp/example.tsv.gz"
    input_config = { "connexion_format": "sqlite" }

    # remove if exists
    remove_if_exists([output_vcf])

    # Create object
    variants = Variants(input=input_vcf, output=output_vcf, config=input_config, load=True)

    # Check get_output
    variants.export_output()
    assert os.path.exists(output_vcf)

    # Check get_output without header
    remove_if_exists([output_vcf])
    variants.export_output(export_header=False)
    assert os.path.exists(output_vcf) and os.path.exists(output_vcf + ".hdr")

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
    input_vcf = tests_folder + "/data/example.vcf.gz"
    param= {"explode_infos": True}

    # Create object
    variants = Variants(input=input_vcf, load=True, param=param)

    # column to check
    column_to_check = "INFO/CLNSIG"
    value_to_check = "pathogenic"

    # check column found
    result = variants.execute_query("SELECT * FROM variants LIMIT 0")
    assert column_to_check in [col[0] for col in result.description]

    # Check value in column
    result = variants.get_query_to_df(f"""SELECT "{column_to_check}" AS column_to_check FROM variants WHERE "#CHROM" = 'chr1' AND POS = 28736 AND REF = 'A' AND ALT = 'C' """)
    assert value_to_check == result["column_to_check"][0]

    # Check number of value in column to check
    result = variants.get_query_to_df(f"""SELECT "{column_to_check}" AS column_to_check FROM variants WHERE "{column_to_check}" IS NOT NULL """)
    assert len(result) == 2


def test_explode_infos_custom():
    """
    This function tests the explode_infos method of the Variants class in Python.
    """

    # Init files
    input_vcf = tests_folder + "/data/example.vcf.gz"
    param= {"explode_infos": "CUSTOM_"}

    # Create object
    variants = Variants(input=input_vcf, load=True, param=param)

    # column to check
    column_to_check = "CUSTOM_CLNSIG"
    value_to_check = "pathogenic"

    # check column found
    result = variants.execute_query("SELECT * FROM variants LIMIT 0")
    assert column_to_check in [col[0] for col in result.description]

    # Check value in column
    result = variants.get_query_to_df(f"""SELECT "{column_to_check}" AS column_to_check FROM variants WHERE "#CHROM" = 'chr1' AND POS = 28736 AND REF = 'A' AND ALT = 'C' """)
    assert value_to_check == result["column_to_check"][0]

    # Check number of value in column to check
    result = variants.get_query_to_df(f"""SELECT "{column_to_check}" AS column_to_check FROM variants WHERE "{column_to_check}" IS NOT NULL """)
    assert len(result) == 2


def test_explode_infos_method():
    """
    This function tests the explode_infos method of the Variants class in Python.
    """

    # Init files
    input_vcf = tests_folder + "/data/example.vcf.gz"

    # Create object
    variants = Variants(input=input_vcf, load=True)

    # Explode infos fields
    variants.explode_infos(prefix="CUSTOM_")

    # column to check
    column_to_check = "CUSTOM_CLNSIG"
    value_to_check = "pathogenic"

    # check column found
    result = variants.execute_query("SELECT * FROM variants LIMIT 0")
    assert column_to_check in [col[0] for col in result.description]

    # Check value in column
    result = variants.get_query_to_df(f"""SELECT "{column_to_check}" AS column_to_check FROM variants WHERE "#CHROM" = 'chr1' AND POS = 28736 AND REF = 'A' AND ALT = 'C' """)
    assert value_to_check == result["column_to_check"][0]

    # Check number of value in column to check
    result = variants.get_query_to_df(f"""SELECT "{column_to_check}" AS column_to_check FROM variants WHERE "{column_to_check}" IS NOT NULL """)
    assert len(result) == 2


def test_explode_infos_no_infos():
    """
    This function tests if the columns in a VCF file remain the same before and after exploding the info
    fields.
    """

    # Init files
    input_vcf = tests_folder + "/data/example.no_samples.vcf.gz"

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
    input_vcf = tests_folder + "/data/example.vcf.gz"
    input_config = { "connexion_format": "sqlite" }

    # Create object
    variants = Variants(input=input_vcf, config=input_config, load=True)

    # Explode infos fields
    variants.explode_infos()

    # Annotation
    variants.annotation()

    # column to check
    column_to_check = "INFO/CLNSIG"
    value_to_check = "pathogenic"

    # check column found
    result = variants.execute_query("SELECT * FROM variants LIMIT 0")
    assert column_to_check in [col[0] for col in result.description]

    # Check value in column
    result = variants.get_query_to_df(f"""SELECT "{column_to_check}" AS column_to_check FROM variants WHERE "#CHROM" = 'chr1' AND POS = 28736 AND REF = 'A' AND ALT = 'C' """)
    assert value_to_check == result["column_to_check"][0]

    # Check number of value in column to check
    result = variants.get_query_to_df(f"""SELECT "{column_to_check}" AS column_to_check FROM variants WHERE "{column_to_check}" IS NOT NULL """)
    assert len(result) == 2


def test_explode_infos_param_prefix():
    """
    This function tests the functionality of exploding info parameters in a VCF file using Python.
    """

    # Init files
    input_vcf = tests_folder + "/data/example.vcf.gz"
    infos_prefix = "INFO_"
    input_param = {"explode_infos": infos_prefix}

    # Create object
    variants = Variants(input=input_vcf, load=True, param=input_param)

    # column to check
    column_to_check = infos_prefix + "CLNSIG"
    value_to_check = "pathogenic"

    # check column found
    result = variants.execute_query("SELECT * FROM variants LIMIT 0")
    assert column_to_check in [col[0] for col in result.description]

    # Check value in column
    result = variants.get_query_to_df(f"""SELECT "{column_to_check}" AS column_to_check FROM variants WHERE "#CHROM" = 'chr1' AND POS = 28736 AND REF = 'A' AND ALT = 'C' """)
    assert value_to_check == result["column_to_check"][0]

    # Check number of value in column to check
    result = variants.get_query_to_df(f"""SELECT "{column_to_check}" AS column_to_check FROM variants WHERE "{column_to_check}" IS NOT NULL """)
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
    input_vcf = tests_folder + "/data/example.vcf.gz"

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
    input_vcf = tests_folder + "/data/example.no_samples.vcf.gz"

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
    input_vcf = tests_folder + "/data/example.vcf.gz"

    # Create object
    variants = Variants(input=input_vcf, load=True)

    # Stats
    stats = variants.get_stats()

    assert stats == None


def test_stats_no_samples():
    """
    This function tests if the get_stats() method returns None when there are no samples in the input
    VCF file.
    """

    # Init files
    input_vcf = tests_folder + "/data/example.no_samples.vcf.gz"

    # Create object
    variants = Variants(input=input_vcf, load=True)

    # Stats
    stats = variants.get_stats()

    assert stats == None
   

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
    assert variant.get_header_columns_as_sql() == ",".join(f'"{col}"' for col in vcf_required_columns)


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
    input_vcf = tests_folder + "/data/example.vcf.gz"

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
    input_vcf = tests_folder + "/data/example.no_samples.vcf.gz"

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
### Annotation
###

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

    # Init files
    input_vcf = tests_folder + "/data/example.vcf.gz"
    annotation1 = tests_folder + "/data/annotations/nci60.parquet"
    annotation2 = tests_folder + "/data/example.vcf.gz"
    annotation3 = tests_folder + "/data/annotations/refGene.bed.gz"
    output_vcf = "/tmp/output.vcf.gz"

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


def test_annotations_no_samples():
    """
    This function tests the annotation of a VCF file without samples using various annotations
    and checks if the output VCF file is in the correct format.
    """

    # Init files
    input_vcf = tests_folder + "/data/example.no_samples.vcf.gz"
    annotation1 = tests_folder + "/data/annotations/nci60.parquet"
    annotation2 = tests_folder + "/data/example.vcf.gz"
    annotation3 = tests_folder + "/data/annotations/refGene.bed.gz"
    output_vcf = "/tmp/output.vcf.gz"

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
    
    for annotation_format in ["vcf", "vcf.gz", "tsv", "tsv.gz", "csv", "csv.gz", "json", "json.gz", "tbl", "tbl", "parquet", "duckdb"]:

        # Init files
        input_vcf = tests_folder + "/data/example.vcf.gz"
        annotation_parquet = tests_folder + f"/data/annotations/nci60.{annotation_format}"
        output_vcf = "/tmp/output.vcf.gz"

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


def test_annotation_parquet_field_already_in_vcf():
    """
    This function tests if a field already present in a VCF file is not changed during annotation with a
    Parquet file.
    """

    # Init files
    input_vcf = tests_folder + "/data/example.vcf.gz"
    annotation1 = tests_folder + "/data/annotations/nci60.parquet"
    output_vcf = "/tmp/output.vcf.gz"

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

    # Create duckdb database

    # Init files
    annotation_parquet = tests_folder + "/data/annotations/nci60.parquet"
    annotation_duckdb = "/tmp/nci60.duckdb"

    remove_if_exists([annotation_duckdb])

    annotation_database = Variants(input=annotation_parquet, output=annotation_duckdb, load=True)
    annotation_database.export_output()

    # Test annotation with duckdb database

    # Init files
    input_vcf = tests_folder + "/data/example.vcf.gz"
    output_vcf = "/tmp/output.vcf.gz"

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

    # Init files
    input_vcf = tests_folder + "/data/example.vcf.gz"
    annotation_parquet = tests_folder + "/data/annotations/nci60.vcf.gz"
    output_vcf = "/tmp/output.vcf.gz"

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

    # Init files
    input_vcf = tests_folder + "/data/example.vcf.gz"
    annotation_parquet = tests_folder + "/data/annotations/refGene.bed.gz"
    output_vcf = "/tmp/output.vcf.gz"

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

    # Init files
    input_vcf = tests_folder + "/data/example.vcf.gz"
    annotation_annovar = "nci60"
    output_vcf = "/tmp/output.vcf.gz"

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

    # Init files
    input_vcf = tests_folder + "/data/example.full.unsorted.vcf.gz"
    annotation_annovar = "nci60"
    output_vcf = "/tmp/output.vcf.gz"

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

    # Init files
    input_vcf = tests_folder + "/data/example.no_samples.vcf.gz"
    annotation_annovar = "nci60"
    output_vcf = "/tmp/output.vcf.gz"

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

    # Init files
    input_vcf = tests_folder + "/data/example.vcf.gz"
    annotation_annovar = "nci60"
    output_vcf = "/tmp/output.vcf.gz"

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

    # # Check if VCF is in correct format with pyVCF
    # variants.export_output()
    # try:
    #     vcf.Reader(filename=output_vcf)
    # except:
    #     assert False


def test_annotation_quick_annovar():
    """
    This function tests the annotation of a VCF file using Annovar.
    """

    # Init files
    input_vcf = tests_folder + "/data/example.vcf.gz"
    annotation_annovar = "nci60"
    output_vcf = "/tmp/output.vcf.gz"

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

    # Init files
    input_vcf = tests_folder + "/data/example.vcf.gz"
    output_vcf = "/tmp/output.vcf.gz"

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

    # Init files
    input_vcf = tests_folder + "/data/example.full.unsorted.vcf.gz"
    output_vcf = "/tmp/output.vcf.gz"

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
    result = variants.get_query_to_df(""" SELECT INFO, "INFO/ANN" FROM variants """)
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

    # Init files
    input_vcf = tests_folder + "/data/example.no_samples.vcf.gz"
    output_vcf = "/tmp/output.vcf.gz"

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

    # Init files
    input_vcf = tests_folder + "/data/example.vcf.gz"
    annotation_snpeff = "snpeff"
    output_vcf = "/tmp/output.vcf.gz"

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

    # Init files
    input_vcf = tests_folder + "/data/example.vcf.gz"
    output_vcf = "/tmp/output.vcf.gz"

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

    # # Check if VCF is in correct format with pyVCF
    # variants.export_output()
    # try:
    #     vcf.Reader(filename=output_vcf)
    # except:
    #     assert False


def test_annotation_bcftools_sqlite():
    """
    This function tests the annotation of a VCF file using bcftools and SQLite.
    """

    # Init files
    input_vcf = tests_folder + "/data/example.vcf.gz"
    annotation_parquet = tests_folder + "/data/annotations/nci60.vcf.gz"
    output_vcf = "/tmp/output.vcf.gz"

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

    # # Check if VCF is in correct format with pyVCF
    # variants.export_output()
    # try:
    #     vcf.Reader(filename=output_vcf)
    # except:
    #     assert False
    
    
def test_annotation_hgvs():
    """
    The function `test_annotation_hgvs` tests the annotation of a VCF file using bcftools and SQLite.
    """

    # Init files
    input_vcf = tests_folder + "/data/example.vcf.gz"
    output_vcf = "/tmp/output.vcf.gz"

    # Construct config dict
    config = tests_config.copy()

    # Exon Version

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


###
### Prioritization
###

def test_prioritization():
    """
    This is a test function for prioritization of variants in a VCF file using a specified configuration
    and parameter dictionary.
    """

    # Init files
    input_vcf = tests_folder + "/data/example.vcf.gz"
    output_vcf = "/tmp/output.vcf.gz"

    # Construct config dict
    config = {
        "prioritization": {
            "config_profiles": tests_folder + "/data/prioritization_profiles.json"
            }
        }
    
    # Construct param dict
    param = {
                "prioritization": {
                    "profiles": ["default", "GERMLINE"],
                    "pzfields": ["PZFlag", "PZScore", "PZComment", "PZInfos"]
                }
        }

    # Create object
    variants = Variants(input=input_vcf, output=output_vcf, load=True, config=config, param=param)

    # Prioritization
    variants.prioritization()

    # Check all priorized default profile
    result = variants.get_query_to_df("""
        SELECT * FROM variants
        WHERE INFO LIKE '%PZFlag_default=%'
          AND INFO LIKE '%PZScore_default=%'
          AND INFO LIKE '%PZComment_default=%'
          AND INFO LIKE '%PZInfos_default=%'
        """)
    assert len(result) == 4

    # Check all priorized GERMLINE profile
    result = variants.get_query_to_df("""
        SELECT * FROM variants
        WHERE INFO LIKE '%PZFlag_GERMLINE=%'
          AND INFO LIKE '%PZScore_GERMLINE=%'
          AND INFO LIKE '%PZComment_default=%'
          AND INFO LIKE '%PZInfos_default=%'
        """)
    assert len(result) == 4

    # Check all priorized default profile (as default)
    result = variants.get_query_to_df("""
        SELECT * FROM variants
        WHERE INFO LIKE '%PZFlag=%'
          AND INFO LIKE '%PZScore=%'
          AND INFO LIKE '%PZComment_default=%'
          AND INFO LIKE '%PZInfos_default=%'
        """)
    assert len(result) == 4

    # Check all priorized default profile
    result = variants.get_query_to_df("""
        SELECT * FROM variants
        WHERE INFO LIKE '%PZFlag_default=%'
          AND INFO LIKE '%PZScore_default=%'
        """)
    assert len(result) == 7

    # Check all priorized GERMLINE profile
    result = variants.get_query_to_df("""
        SELECT * FROM variants
        WHERE INFO LIKE '%PZFlag_GERMLINE=%'
          AND INFO LIKE '%PZScore_GERMLINE=%'
        """)
    assert len(result) == 7

    # Check all priorized default profile (as default)
    result = variants.get_query_to_df("""
        SELECT * FROM variants
        WHERE INFO LIKE '%PZFlag=%'
          AND INFO LIKE '%PZScore=%'
        """)
    assert len(result) == 7

    # Check annotation1
    result = variants.get_query_to_df(""" SELECT * FROM variants WHERE "#CHROM" = 'chr1' AND POS = 28736 AND REF = 'A' AND ALT = 'C' AND INFO LIKE '%PZScore_default=15%' """)
    assert len(result) == 1

    # Check FILTERED
    result = variants.get_query_to_df(f""" SELECT INFO FROM variants WHERE INFO LIKE '%FILTERED%' """)
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

    # Init files
    input_vcf = tests_folder + "/data/example.full.unsorted.vcf.gz"
    output_vcf = "/tmp/output.vcf.gz"

    # Construct config dict
    config = {
        "prioritization": {
            "config_profiles": tests_folder + "/data/prioritization_profiles.json"
            }
        }
    
    # Construct param dict
    param = {
                "prioritization": {
                    "profiles": ["default", "GERMLINE"],
                    "pzfields": ["PZFlag", "PZScore", "PZComment", "PZInfos"]
                }
        }

    # Create object
    variants = Variants(input=input_vcf, output=output_vcf, load=True, config=config, param=param)

    # Prioritization
    variants.prioritization()

    # Check all priorized default profile
    result = variants.get_query_to_df("""
        SELECT * FROM variants
        WHERE INFO LIKE '%PZFlag_default=%'
          AND INFO LIKE '%PZScore_default=%'
          AND INFO LIKE '%PZComment_default=%'
          AND INFO LIKE '%PZInfos_default=%'
        """)
    assert len(result) == 4

    # Check all priorized GERMLINE profile
    result = variants.get_query_to_df("""
        SELECT * FROM variants
        WHERE INFO LIKE '%PZFlag_GERMLINE=%'
          AND INFO LIKE '%PZScore_GERMLINE=%'
          AND INFO LIKE '%PZComment_GERMLINE=%'
          AND INFO LIKE '%PZInfos_GERMLINE=%'
        """)
    assert len(result) == 2

    # Check all priorized default profile (as default)
    result = variants.get_query_to_df("""
        SELECT * FROM variants
        WHERE INFO LIKE '%PZFlag=%'
          AND INFO LIKE '%PZScore=%'
          AND INFO LIKE '%PZComment=%'
          AND INFO LIKE '%PZInfos=%'
        """)
    assert len(result) == 4

    # Check all priorized default profile
    result = variants.get_query_to_df("""
        SELECT * FROM variants
        WHERE INFO LIKE '%PZFlag_default=%'
          AND INFO LIKE '%PZScore_default=%'
        """)
    assert len(result) == 36

    # Check all priorized GERMLINE profile
    result = variants.get_query_to_df("""
        SELECT * FROM variants
        WHERE INFO LIKE '%PZFlag_GERMLINE=%'
          AND INFO LIKE '%PZScore_GERMLINE=%'
        """)
    assert len(result) == 36

    # Check all priorized default profile (as default)
    result = variants.get_query_to_df("""
        SELECT * FROM variants
        WHERE INFO LIKE '%PZFlag=%'
          AND INFO LIKE '%PZScore=%'
        """)
    assert len(result) == 36

    # Check annotation1
    result = variants.get_query_to_df(""" SELECT * FROM variants WHERE "#CHROM" = 'chr1' AND POS = 28736 AND REF = 'A' AND ALT = 'C' AND INFO LIKE '%PZScore_default=15%' """)
    assert len(result) == 1

    # Check FILTERED
    result = variants.get_query_to_df(f""" SELECT INFO FROM variants WHERE INFO LIKE '%FILTERED%' """)
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

    # Init files
    input_vcf = tests_folder + "/data/example.vcf.gz"
    output_vcf = "/tmp/output.vcf.gz"

    # Construct config dict
    config = {
        "prioritization": {
            "config_profiles": "config/prioritization_profiles.json"
            }
        }

    # Construct param dict
    param = {
                "prioritization": {
                    "profiles": ["default", "GERMLINE"],
                    "pzfields": ["PZFlag", "PZScore", "PZComment", "PZInfos"],
                    "prioritization_score_mode": "VaRank"
                }
        }

    # Create object
    variants = Variants(input=input_vcf, output=output_vcf, load=True, config=config, param=param)

    # Prioritization
    variants.prioritization()

    # Check all priorized
    result = variants.get_query_to_df(""" SELECT INFO FROM variants """)
    assert len(result) > 0

    # Check annotation1
    result = variants.get_query_to_df(""" SELECT 1 AS count FROM variants WHERE "#CHROM" = 'chr1' AND POS = 28736 AND REF = 'A' AND ALT = 'C' AND INFO LIKE '%PZScore_default=15%' """)
    assert len(result) == 1

    # Check if VCF is in correct format with pyVCF
    remove_if_exists([output_vcf])
    variants.export_output()
    try:
        vcf.Reader(filename=output_vcf)
    except:
        assert False


def test_prioritization_no_profiles():
    """
    This function tests if an error is raised when there are no prioritization configuration profiles
    provided.
    """

    # Init files
    input_vcf = tests_folder + "/data/example.vcf.gz"

    # Construct config dict
    config = {
        "prioritization": {
            "config_profiles": None
            }
        }
    # Create object
    variants = Variants(input=input_vcf, load=True, config=config)

    # Prioritization fail
    with pytest.raises(ValueError) as e:
        variants.prioritization()
    assert str(e.value) == f"NO Profiles configuration"


def test_prioritization_no_pzfields():
    """
    This function tests the prioritization method of a Variants object when there are no pzfields
    specified.
    """

    # Init files
    input_vcf = tests_folder + "/data/example.vcf.gz"
    output_vcf = "/tmp/output.vcf.gz"

    # Construct config dict
    config = {
        "prioritization": {
            "config_profiles": "config/prioritization_profiles.json"
            }
        }

    # Construct param dict
    param = {
                "prioritization": {
                    "profiles": [],
                    "pzfields": []
                }
        }

    # Create object
    variants = Variants(input=input_vcf, output=output_vcf, load=True, config=config, param=param)

    # Prioritization
    variants.prioritization()

    # Check all priorized
    result = variants.get_query_to_df(""" SELECT * FROM variants WHERE INFO LIKE '%PZScore_default=%' """)
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

    # Init files
    input_vcf = tests_folder + "/data/example.no_samples.vcf.gz"
    output_vcf = "/tmp/output.vcf.gz"

    # Construct config dict
    config = {
        "prioritization": {
            "config_profiles": tests_folder + "/data/prioritization_profiles.json"
            }
        }
    
    # Construct param dict
    param = {
                "prioritization": {
                    "profiles": ["default", "GERMLINE"],
                    "pzfields": ["PZFlag", "PZScore", "PZComment", "PZInfos"]
                }
        }

    # Create object
    variants = Variants(input=input_vcf, output=output_vcf, load=True, config=config, param=param)

    # Prioritization
    variants.prioritization()

    result = variants.get_query_to_df("""
        SELECT INFO FROM variants WHERE INFO LIKE '%PZScore=0%' AND INFO LIKE '%PZFlag=PASS%'
        """)
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

def test_calculation():
    """
    This function tests the calculation and annotation of genetic variants using input parameters and
    checks if the output VCF file is in the correct format.
    """

    # Init files
    input_vcf = tests_folder + "/data/example.vcf.gz"
    output_vcf = "/tmp/output.vcf.gz"
    input_param = {
            "annotation": {
                "annovar": {
                    "annotations": {
                        "refGene": {
                        "Func_refGene": "location",
                        "Gene_refGene": "gene",
                        "GeneDetail_refGene": "GeneDetail",
                        "ExonicFunc_refGene": "outcome",
                        "AAChange_refGene": "hgvs"
                        }
                    },
                    "options": {
                        "genebase": "-hgvs -splicing_threshold 3 ",
                        "intronhgvs": 10
                    }
                }
            },
            "calculation": {
                "NOMEN": {
                    "options": {
                        "hgvs_field": "hgvs"
                    }
                },
                "middle": None,
                "no_existing_calculation": None
            }
        }

    # Create object
    variants = Variants(input=input_vcf, output=output_vcf, config=tests_config, param=input_param, load=True)

    # Annotation
    variants.annotation()

    # Calculation
    variants.calculation()

    # Check number of NOMEN (2)
    result = variants.get_query_to_df("""SELECT INFO FROM variants WHERE INFO LIKE '%NOMEN=%' """)
    assert len(result) == 2

     # Check number of middle (7)
    result = variants.get_query_to_df("""SELECT INFO FROM variants WHERE INFO LIKE '%middle=%' """)
    assert len(result) == 7

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

    # Init files
    input_vcf = tests_folder + "/data/example.snv.indel.mosaic.vcf"
    output_vcf = "/tmp/output.vcf.gz"

    # Construct param dict
    param = {
        "calculation": {
            "VARTYPE": None
        }
    }

    # Create object
    variants = Variants(conn=None, input=input_vcf, output=output_vcf, param=param, load=True)

    # Calculation
    variants.calculation()

    result = variants.get_query_to_df(""" SELECT * FROM variants WHERE INFO LIKE '%VARTYPE=SNV%' """)
    assert len(result) == 5
    
    result = variants.get_query_to_df(""" SELECT * FROM variants WHERE INFO LIKE '%VARTYPE=INDEL%' """)
    assert len(result) == 1
    
    result = variants.get_query_to_df(""" SELECT * FROM variants WHERE INFO LIKE '%VARTYPE=MOSAIC%' """)
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

    # Init files
    input_vcf = tests_folder + "/data/example.full.unsorted.vcf.gz"
    output_vcf = "/tmp/output.vcf.gz"

    # Construct param dict
    param = {
        "calculation": {
            "VARTYPE": None
        }
    }

    # Create object
    variants = Variants(conn=None, input=input_vcf, output=output_vcf, param=param, load=True)

    # Calculation
    variants.calculation()

    result = variants.get_query_to_df(""" SELECT * FROM variants WHERE INFO LIKE '%VARTYPE=SNV%' """)
    assert len(result) == 4
    
    result = variants.get_query_to_df(""" SELECT * FROM variants WHERE INFO LIKE '%VARTYPE=INDEL%' """)
    assert len(result) == 2
    
    result = variants.get_query_to_df(""" SELECT * FROM variants WHERE INFO LIKE '%VARTYPE=CNV%' """)
    assert len(result) == 3

    result = variants.get_query_to_df(""" SELECT * FROM variants WHERE INFO LIKE '%VARTYPE=INV%' """)
    assert len(result) == 3

    result = variants.get_query_to_df(""" SELECT * FROM variants WHERE INFO LIKE '%VARTYPE=DEL%' """)
    assert len(result) == 3

    result = variants.get_query_to_df(""" SELECT * FROM variants WHERE INFO LIKE '%VARTYPE=INS%' """)
    assert len(result) == 5

    result = variants.get_query_to_df(""" SELECT * FROM variants WHERE INFO LIKE '%VARTYPE=DUP%' """)
    assert len(result) == 6
    
    result = variants.get_query_to_df(""" SELECT * FROM variants WHERE INFO LIKE '%VARTYPE=BND%' """)
    assert len(result) == 7

    result = variants.get_query_to_df(""" SELECT * FROM variants WHERE INFO LIKE '%VARTYPE=MNV%' """)
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

    # Init files
    input_vcf = tests_folder + "/data/example.ann.vcf.gz"
    output_vcf = "/tmp/output.vcf.gz"

    # Construct param dict
    param = {
        "calculation": {
            "snpeff_hgvs": None
        }
    }

    # Create object
    variants = Variants(conn=None, input=input_vcf, output=output_vcf, param=param, load=True)

    # Check if no snpeff_hgvs
    result = variants.get_query_to_df(""" SELECT INFO FROM variants WHERE INFO LIKE '%snpeff_hgvs=%' """)
    assert len(result) == 0

    # Calculation
    variants.calculation()

    # query annotated variant
    result = variants.get_query_to_df(""" SELECT * FROM variants WHERE INFO LIKE '%snpeff_hgvs=%' """)
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

    # Init files
    input_vcf = tests_folder + "/data/example.vcf.gz"
    output_vcf = "/tmp/output.vcf.gz"

    # Construct param dict
    param = {
        "calculation": {
            "snpeff_hgvs": None
        }
    }

    # Create object
    variants = Variants(conn=None, input=input_vcf, output=output_vcf, param=param, load=True)

    # Calculation
    variants.calculation()

    # query annotated variant
    result = variants.get_query_to_df(""" SELECT * FROM variants WHERE INFO LIKE '%snpeff_hgvs=%' """)
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

    # Init files
    input_vcf = tests_folder + "/data/example.snpeff.vcf.gz"
    transcripts_file = tests_folder + "/data/transcripts.tsv"
    output_vcf = "/tmp/output.vcf.gz"

    # Construct param dict
    param = {
        "calculation": {
            "NOMEN": {
                "options": {
                    "hgvs_field": "snpeff_hgvs",
                    "transcripts": transcripts_file
                }
            }
        }
    }

    # Create object
    variants = Variants(conn=None, input=input_vcf, output=output_vcf, param=param, load=True)

    # Calculation
    variants.calculation()

    # query annotated variant
    result = variants.get_query_to_df(""" SELECT * FROM variants WHERE INFO LIKE '%NOMEN=%' """)
    assert len(result) == 7

    # Check transcript priority
    result = variants.get_query_to_df(""" SELECT * FROM variants WHERE INFO LIKE '%NOMEN=EGFR:NM_001346897%' """)
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
    input_vcf = tests_folder + "/data/example.snpeff.vcf.gz"
    transcripts_file = tests_folder + "/data/notranscripts.tsv"

    # Construct param dict
    param = {
        "calculation": {
            "NOMEN": {
                "options": {
                    "hgvs_field": "snpeff_hgvs",
                    "transcripts": transcripts_file
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

    # Init files
    input_vcf = tests_folder + "/data/example.vcf.gz"
    output_vcf = "/tmp/output.vcf.gz"

    # Construct param dict
    param = {
        "calculation": {
            "FINDBYPIPELINE": None
        }
    }

    # Create object
    variants = Variants(conn=None, input=input_vcf, output=output_vcf, param=param, load=True)

    # Calculation
    variants.calculation()

    result = variants.get_query_to_df(""" SELECT INFO FROM variants WHERE INFO LIKE '%findbypipeline%' """)
    assert len(result) == 7

    result = variants.get_query_to_df(""" SELECT * FROM variants WHERE INFO LIKE '%findbypipeline=4/4%' """)
    assert len(result) == 1
    
    result = variants.get_query_to_df(""" SELECT * FROM variants WHERE INFO LIKE '%findbypipeline=3/4%' """)
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

    # Init files
    input_vcf = tests_folder + "/data/example.vcf.gz"
    output_vcf = "/tmp/output.vcf.gz"

    # Construct param dict
    param = {
        "calculation": {
            "FINDBYSAMPLE": None
        }
    }

    # Create object
    variants = Variants(conn=None, input=input_vcf, output=output_vcf, param=param, load=True)

    # Calculation
    variants.calculation()

    result = variants.get_query_to_df(""" SELECT INFO FROM variants WHERE INFO LIKE '%findbysample%' """)
    assert len(result) == 7

    result = variants.get_query_to_df(""" SELECT * FROM variants WHERE INFO LIKE '%findbysample=4/4%' """)
    assert len(result) == 1
    
    result = variants.get_query_to_df(""" SELECT * FROM variants WHERE INFO LIKE '%findbysample=3/4%' """)
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

    # Init files
    input_vcf = tests_folder + "/data/example.vcf.gz"
    output_vcf = "/tmp/output.vcf.gz"

    # Construct param dict
    param = {
        "calculation": {
            "GENOTYPECONCORDANCE": None
        }
    }

    # Create object
    variants = Variants(conn=None, input=input_vcf, output=output_vcf, param=param, load=True)

    # Calculation
    variants.calculation()

    result = variants.get_query_to_df(""" SELECT INFO FROM variants WHERE INFO LIKE '%genotypeconcordance%' """)
    assert len(result) == 7

    result = variants.get_query_to_df(""" SELECT * FROM variants WHERE INFO LIKE '%genotypeconcordance=FALSE%' """)
    assert len(result) == 1
    
    result = variants.get_query_to_df(""" SELECT * FROM variants WHERE INFO LIKE '%genotypeconcordance=TRUE%' """)
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

    # Init files
    input_vcf = tests_folder + "/data/example.vcf.gz"
    output_vcf = "/tmp/output.vcf.gz"

    # Construct param dict
    param = {
        "calculation": {
            "BARCODE": None
        }
    }

    # Create object
    variants = Variants(conn=None, input=input_vcf, output=output_vcf, param=param, load=True)

    # Calculation
    variants.calculation()

    result = variants.get_query_to_df(""" SELECT INFO FROM variants WHERE INFO LIKE '%barcode%' """)
    assert len(result) == 7

    result = variants.get_query_to_df(""" SELECT * FROM variants WHERE INFO LIKE '%barcode=1122%' """)
    assert len(result) == 1

    result = variants.get_query_to_df(""" SELECT * FROM variants WHERE INFO LIKE '%barcode=0111%' """)
    assert len(result) == 1
    
    result = variants.get_query_to_df(""" SELECT * FROM variants WHERE INFO LIKE '%barcode=1011%' """)
    assert len(result) == 4
    
    result = variants.get_query_to_df(""" SELECT * FROM variants WHERE INFO LIKE '%barcode=1101%' """)
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

    # Init files
    input_vcf = tests_folder + "/data/example.vcf.gz"
    output_vcf = "/tmp/output.vcf.gz"

    # Construct param dict
    param = {
        "calculation": {
            "trio": {
                "father": "sample1",
                "mother": "sample2",
                "child": "sample3"
            }
        }
    }

    # Create object
    variants = Variants(conn=None, input=input_vcf, output=output_vcf, param=param, load=True)

    # Calculation
    variants.calculation()

    result = variants.get_query_to_df(""" SELECT INFO FROM variants WHERE INFO LIKE '%trio=recessive%' """)
    assert len(result) == 1

    result = variants.get_query_to_df(""" SELECT * FROM variants WHERE INFO LIKE '%trio=dominant%' """)
    assert len(result) == 5

    result = variants.get_query_to_df(""" SELECT * FROM variants WHERE INFO LIKE '%trio=unknown%' """)
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

    # Init files
    input_vcf = tests_folder + "/data/example.vcf.gz"
    output_vcf = "/tmp/output.vcf.gz"

    # Construct param dict
    param = {
        "calculation": {
            "vaf": None
        }
    }

    # Create object
    variants = Variants(conn=None, input=input_vcf, output=output_vcf, param=param, load=True)

    # Calculation
    variants.calculation()

    result = variants.get_query_to_df(""" SELECT INFO FROM variants WHERE FORMAT LIKE '%:VAF' """)
    assert len(result) == 7

    result = variants.get_query_to_df(""" SELECT * FROM variants WHERE "#CHROM" = 'chr1' AND POS = 28736 AND sample1 LIKE '%:0.279835' """)
    assert len(result) == 1

    result = variants.get_query_to_df(""" SELECT * FROM variants WHERE "#CHROM" = 'chr1' AND POS = 28736 AND sample2 LIKE '%:0.282898' """)
    assert len(result) == 1
    
    result = variants.get_query_to_df(""" SELECT * FROM variants WHERE "#CHROM" = 'chr1' AND POS = 28736 AND sample3 LIKE '%:0.282955' """)
    assert len(result) == 1

    result = variants.get_query_to_df(""" SELECT * FROM variants WHERE "#CHROM" = 'chr1' AND POS = 28736 AND sample4 LIKE '%:0.303819' """)
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

    # Init files
    input_vcf = tests_folder + "/data/example.vcf.gz"
    output_vcf = "/tmp/output.vcf.gz"

    # Construct param dict
    param = {
        "calculation": {
            "vaf": None,
            "vaf_stats": None
        }
    }

    # Create object
    variants = Variants(conn=None, input=input_vcf, output=output_vcf, param=param, load=True)

    # Calculation
    variants.calculation()

    result = variants.get_query_to_df(""" SELECT INFO FROM variants WHERE INFO LIKE '%VAF_stats%' """)
    assert len(result) == 7

    result = variants.get_query_to_df(""" SELECT * FROM variants WHERE "#CHROM" = 'chr1' AND POS = 28736 AND INFO LIKE '%VAF_stats_nb=4%' """)
    assert len(result) == 1

    result = variants.get_query_to_df(""" SELECT * FROM variants WHERE "#CHROM" = 'chr1' AND POS = 28736 AND INFO LIKE '%VAF_stats_min=0.279835%' """)
    assert len(result) == 1
    
    result = variants.get_query_to_df(""" SELECT * FROM variants WHERE "#CHROM" = 'chr1' AND POS = 28736 AND INFO LIKE '%VAF_stats_max=0.303819%' """)
    assert len(result) == 1

    result = variants.get_query_to_df(""" SELECT * FROM variants WHERE "#CHROM" = 'chr1' AND POS = 28736 AND INFO LIKE '%VAF_stats_mean=0.28737675%' """)
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

    # Init files
    input_vcf = tests_folder + "/data/example.vcf.gz"
    output_vcf = "/tmp/output.vcf.gz"

    # Construct param dict
    param = {
        "calculation": {
            "dp_stats": None
        }
    }

    # Create object
    variants = Variants(conn=None, input=input_vcf, output=output_vcf, param=param, load=True)

    # Remove if output file exists
    remove_if_exists([output_vcf])

    # Calculation
    variants.calculation()

    result = variants.get_query_to_df(""" SELECT INFO FROM variants WHERE INFO LIKE '%DP_stats%' """)
    assert len(result) == 7

    result = variants.get_query_to_df(""" SELECT * FROM variants WHERE "#CHROM" = 'chr1' AND POS = 28736 AND INFO LIKE '%DP_stats_nb=4%' """)
    assert len(result) == 1

    result = variants.get_query_to_df(""" SELECT * FROM variants WHERE "#CHROM" = 'chr1' AND POS = 28736 AND INFO LIKE '%DP_stats_min=576.0%' """)
    assert len(result) == 1
    
    result = variants.get_query_to_df(""" SELECT * FROM variants WHERE "#CHROM" = 'chr1' AND POS = 28736 AND INFO LIKE '%DP_stats_max=17664.0%' """)
    assert len(result) == 1

    result = variants.get_query_to_df(""" SELECT * FROM variants WHERE "#CHROM" = 'chr1' AND POS = 28736 AND INFO LIKE '%DP_stats_mean=9158.0%' """)
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

    # Init files
    input_vcf = tests_folder + "/data/example.vcf.gz"
    output_vcf = "/tmp/output.vcf.gz"

    # Construct param dict
    param = {
        "calculation": {
            "variant_id": None
        }
    }

    # Create object
    variants = Variants(conn=None, input=input_vcf, output=output_vcf, param=param, load=True)

    # Remove if output file exists
    remove_if_exists([output_vcf])

    # Calculation
    variants.calculation()

    # Check if all variant have variant_id
    result = variants.get_query_to_df(""" SELECT INFO FROM variants WHERE INFO LIKE '%variant_id%' """)
    assert len(result) == 7

    # Explode info
    variants.explode_infos(prefix="INFO/")

    # Check distinct variant_id
    result = variants.get_query_to_df(""" SELECT distinct "INFO/variant_id" FROM variants """)
    assert len(result) == 7
    
    # Check if VCF is in correct format with pyVCF
    remove_if_exists([output_vcf])
    variants.export_output()
    try:
        vcf.Reader(filename=output_vcf)
    except:
        assert False


def test_get_operations_help():
    operations = { "add": { "name": "addition", "description": "sum two numbers", "available": True, }, "subtract": { "name": "subtraction", "description": "subtract two numbers", "available": False }, }

    expected_operations_help = [
        "Available calculation operations:"
    ]

    variants = Variants()
    actual_operations_help = variants.get_operations_help()

    print(actual_operations_help[0])

    assert expected_operations_help[0] == actual_operations_help[0]
