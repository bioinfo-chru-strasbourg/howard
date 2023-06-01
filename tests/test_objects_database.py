# -*- coding: utf-8 -*-
"""
Tests

Usage:
pytest tests/

Coverage:
coverage run -m pytest . -x -v
coverage report --include=howard/* -m
"""

# import logging as log
import os
# import sys
# import duckdb
# import re
# import Bio.bgzf as bgzf
# import gzip
import pytest

from howard.commons import *
from howard.objects.database import Database


# Main tests folder
tests_folder = os.path.dirname(__file__)

# Annotation subfolder
tests_annotations_folder = tests_folder + "/data/annotations"

# Annotation databases
database_files = {
    "parquet" : tests_annotations_folder + "/nci60.parquet",
    "duckdb" : tests_annotations_folder + "/nci60.duckdb",
    "vcf" : tests_annotations_folder + "/nci60.vcf",
    "vcf_gz" : tests_annotations_folder + "/nci60.vcf.gz",
    "tsv" : tests_annotations_folder + "/nci60.tsv",
    "tsv_alternative_columns" : tests_annotations_folder + "/nci60.alternative_columns.tsv",
    "tsv_failed_columns" : tests_annotations_folder + "/nci60.failed_columns.tsv",
    "tsv_lower_columns" : tests_annotations_folder + "/nci60.lower_columns.tsv",
    "tsv_gz" : tests_annotations_folder + "/nci60.tsv.gz",
    "csv" : tests_annotations_folder + "/nci60.csv",
    "csv_gz" : tests_annotations_folder + "/nci60.csv.gz",
    "psv" : tests_annotations_folder + "/nci60.psv",
    "psv_gz" : tests_annotations_folder + "/nci60.psv.gz",
    "json" : tests_annotations_folder + "/nci60.json",
    "json_gz" : tests_annotations_folder + "/nci60.json.gz",
    "bed" : tests_annotations_folder + "/annotation_regions.bed",
    "bed_gz" : tests_annotations_folder + "/annotation_regions.bed.gz"
}



def test_empty_database():

    # Create Empty database
    try:
        database = Database()
    except:
        assert False
    

def test_get_database():
    """
    This function tests the set_database and get_database methods of the Database class in Python.
    """

    # Init files
    database_file = database_files.get("parquet") #tests_folder + "/data/annotations/nci60.parquet"

    # Create object
    database = Database(database=database_file)

    # Check get_database
    assert database_file == database.get_database()


def test_set_database():
    """
    This function tests the set_database and get_database methods of the Database class in Python.
    """

    # Init files
    database_file = database_files.get("parquet")
    new_database_file = database_files.get("vcf")

    # Create object
    database = Database(database=database_file)

    # Check get_database
    assert database_file == database.get_database()

    # Set new Database
    database.set_database(database=new_database_file)

    # Check get_database
    assert new_database_file == database.get_database()


def test_exists():
    """
    This function test if a database file exists or not
    """

    # Init files
    database_file = database_files.get("parquet")
    database_file_fake = database_files.get("parquet") + "FAKE"

    # Create object
    database = Database()
    
    # Check if database exists
    assert database.exists(database=database_file)

    # Check if database DO NOT exists
    assert not database.exists(database=database_file_fake)

    # Check if database DO NOT exists if None
    assert not database.exists(database=None)


def test_get_format():
    """
    This function test format database
    """

    # Create object
    database = Database()

    # Check parquet
    assert database.get_format(database_files.get("parquet")) == "parquet"

    # Check duckdb
    assert database.get_format(database_files.get("duckdb")) == "duckdb"

    # Check vcf
    assert database.get_format(database_files.get("vcf")) == "vcf"
    assert database.get_format(database_files.get("vcf_gz")) == "vcf"

     # Check tsv
    assert database.get_format(database_files.get("tsv")) == "tsv"
    assert database.get_format(database_files.get("tsv_gz")) == "tsv"

    # Check csv
    assert database.get_format(database_files.get("csv")) == "csv"
    assert database.get_format(database_files.get("csv_gz")) == "csv"

    # Check json
    assert database.get_format(database_files.get("json")) == "json"
    assert database.get_format(database_files.get("json_gz")) == "json"

    # Check bed
    assert database.get_format(database_files.get("bed")) == "bed"
    assert database.get_format(database_files.get("bed_gz")) == "bed"

    # Check None
    assert database.get_format(None) == "unknown"


def test_get_is_compressed():
    """
    This function test if database file are compressed or not
    """

    # Create object
    database = Database()

    # Check parquet
    assert not database.is_compressed(database_files.get("parquet"))

    # Check duckdb
    assert not database.is_compressed(database_files.get("duckdb"))

    # Check vcf
    assert not database.is_compressed(database_files.get("vcf"))
    assert database.is_compressed(database_files.get("vcf_gz"))

     # Check tsv
    assert not database.is_compressed(database_files.get("tsv"))
    assert database.is_compressed(database_files.get("tsv_gz"))

    # Check csv
    assert not database.is_compressed(database_files.get("csv"))
    assert database.is_compressed(database_files.get("csv_gz"))

    # Check json
    assert not database.is_compressed(database_files.get("json"))
    assert database.is_compressed(database_files.get("json_gz"))

    # Check bed
    assert not database.is_compressed(database_files.get("bed"))
    assert database.is_compressed(database_files.get("bed_gz"))

    # Check None
    assert not database.is_compressed(None)


def test_get_type():
    """
    This function test type of database (either variant VCF-like or region BED-like)
    """

    # Create object
    database = Database()

    # Check parquet VCF/BED
    assert database.get_type(database_files.get("parquet")) == "variants"
    assert database.get_type(database_files.get("parquet")) == "variants"

    # Check duckdb variants/regions
    assert database.get_type(database_files.get("duckdb")) == "variants"
    assert database.get_type(database_files.get("duckdb")) == "variants"

    # Check vcf variants/regions
    assert database.get_type(database_files.get("vcf")) == "variants"
    assert database.get_type(database_files.get("vcf_gz")) == "variants"

     # Check tsv
    assert database.get_type(database_files.get("tsv")) == "variants"
    assert database.get_type(database_files.get("tsv_gz")) == "variants"
    assert database.get_type(database_files.get("tsv_alternative_columns")) == "variants"
    assert database.get_type(database_files.get("tsv_failed_columns")) == None
    assert database.get_type(database_files.get("tsv_lower_columns")) == "variants"

    # Check csv
    assert database.get_type(database_files.get("csv")) == "variants"
    assert database.get_type(database_files.get("csv_gz")) == "variants"

    # Check json
    assert database.get_type(database_files.get("json")) == None
    assert database.get_type(database_files.get("json_gz")) == None

    # Check bed
    assert database.get_type(database_files.get("bed")) == "regions"
    assert database.get_type(database_files.get("bed_gz")) == "regions"

    # Check None
    assert database.get_type(None) == None


def test_is_vcf():
    """
    This function test is a database is a vcf (contains all needed columns)
    """

    # Create object
    database = Database(database_files.get("vcf"))

    # Check duckdb
    assert database.is_vcf()

    # Create object
    database = Database()

    # Check parquet VCF/BED
    assert database.is_vcf(database_files.get("parquet"))
    assert database.is_vcf(database_files.get("parquet"))

    # Check duckdb variants/regions
    assert database.is_vcf(database_files.get("duckdb"))
    assert database.is_vcf(database_files.get("duckdb"))

    # Check vcf variants/regions
    assert database.is_vcf(database_files.get("vcf"))
    assert database.is_vcf(database_files.get("vcf_gz"))

     # Check tsv
    assert database.is_vcf(database_files.get("tsv"))
    assert database.is_vcf(database_files.get("tsv_gz"))
    assert database.is_vcf(database_files.get("tsv_alternative_columns"))
    assert not database.is_vcf(database_files.get("tsv_failed_columns"))
    assert database.is_vcf(database_files.get("tsv_lower_columns"))

    # Check csv
    assert database.is_vcf(database_files.get("csv"))
    assert database.is_vcf(database_files.get("csv_gz"))

    # Check json
    assert not database.is_vcf(database_files.get("json"))
    assert not database.is_vcf(database_files.get("json_gz"))

    # Check bed
    assert not database.is_vcf(database_files.get("bed"))
    assert not database.is_vcf(database_files.get("bed_gz"))

    # Check None
    assert not database.is_vcf(None)


def test_get_database_tables():
    """
    This function list tables in a duckdb database
    """

    # Create object
    database = Database(database_files.get("duckdb"))

    # Check duckdb
    assert database.get_database_tables() == ["variants"]

    # Create empty object
    database = Database()

    # Check duckdb
    assert database.get_database_tables(database_files.get("duckdb")) == ["variants"]

    # Check parquet
    assert database.get_database_tables(database_files.get("parquet")) is None

    # Check vcf
    assert database.get_database_tables(database_files.get("vcf")) is None

    # Check None
    assert database.get_database_tables(None) is None


def test_get_database_table():
    """
    This function check which table in a duckdb database contain annotation (variants or regions)
    """

    # Create object
    database = Database(database_files.get("duckdb"))

    # Check duckdb
    assert database.get_database_table() == "variants"

    # Create empty object
    database = Database()

    # Check duckdb
    assert database.get_database_table(database_files.get("duckdb")) == "variants"

    # Check parquet
    assert database.get_database_table(database_files.get("parquet")) is None

    # Check vcf
    assert database.get_database_table(database_files.get("vcf")) is None

    # Check None
    assert database.get_database_table(None) is None


def test_get_sql_from():
    """
    This function get sql from section from a database
    """

    # Create object
    database = Database(database_files.get("duckdb"))

    # Check duckdb
    assert database.get_sql_from() == f"""'{database_files.get("duckdb")}'"""

    # Create empty object
    database = Database()

    # Check duckdb
    assert database.get_sql_from(database_files.get("duckdb")) == f"""'{database_files.get("duckdb")}'"""

    # Check parquet
    assert database.get_sql_from(database_files.get("parquet")) ==  f"""read_parquet('{database_files.get("parquet")}')"""

    # Check vcf
    assert database.get_sql_from(database_files.get("vcf")) == f"""read_csv('{database_files.get("vcf")}', auto_detect=True, delim='\t')"""

    # Check tsv
    assert database.get_sql_from(database_files.get("tsv")) == f"""read_csv('{database_files.get("tsv")}', auto_detect=True, delim='\t')"""

    # Check csv
    assert database.get_sql_from(database_files.get("csv")) == f"""read_csv('{database_files.get("csv")}', auto_detect=True, delim=',')"""

    # Check psv
    assert database.get_sql_from(database_files.get("psv")) == f"""read_csv('{database_files.get("psv")}', auto_detect=True, delim='|')"""

    # Check bed
    assert database.get_sql_from(database_files.get("bed")) == f"""read_csv('{database_files.get("bed")}', auto_detect=True, delim='\t')"""

    # Check None
    assert database.get_sql_from(None) == None


def test_get_colums():
    """
    This function get colums ofa  database
    """

    # Create object
    database = Database(database_files.get("duckdb"))

    # Check duckdb
    assert database.get_columns() == []
    assert database.get_columns(table="variants") == ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'INFO/nci60']

    # Create empty object
    database = Database()

    # Check duckdb
    assert database.get_columns(database_files.get("duckdb")) == []
    assert database.get_columns(database_files.get("duckdb"), table="variants") == ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'INFO/nci60']

    # Check parquet
    assert database.get_columns(database_files.get("parquet")) ==  ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'INFO/nci60']

    # Check vcf
    assert database.get_columns(database_files.get("vcf")) == ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO']

    # Check tsv
    assert database.get_columns(database_files.get("tsv")) == ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO']

    # Check csv
    assert database.get_columns(database_files.get("csv")) == ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO']

    # Check psv
    assert database.get_columns(database_files.get("psv")) == ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO']

    # Check bed
    assert database.get_columns(database_files.get("bed")) == ['#CHROM', 'START', 'END', 'annot1', 'annot2']

    # Check None
    assert database.get_columns(None) == []


def test_get_extra_colums():
    """
    This function get colums ofa  database
    """

    # Create object
    database = Database(database_files.get("duckdb"))

    # Check duckdb
    assert database.get_extra_columns() == ['ID', 'QUAL', 'FILTER', 'INFO', 'INFO/nci60']

    # Create empty object
    database = Database()

    # Check duckdb
    assert database.get_extra_columns(database_files.get("duckdb")) == ['ID', 'QUAL', 'FILTER', 'INFO', 'INFO/nci60']

    # Check parquet
    assert database.get_extra_columns(database_files.get("parquet")) ==  ['ID', 'QUAL', 'FILTER', 'INFO', 'INFO/nci60']

    # Check vcf
    assert database.get_extra_columns(database_files.get("vcf")) == ['ID', 'QUAL', 'FILTER', 'INFO']

    # Check tsv
    assert database.get_extra_columns(database_files.get("tsv")) == ['ID', 'QUAL', 'FILTER', 'INFO']

    # Check csv
    assert database.get_extra_columns(database_files.get("csv")) == ['ID', 'QUAL', 'FILTER', 'INFO']

    # Check psv
    assert database.get_extra_columns(database_files.get("psv")) == ['ID', 'QUAL', 'FILTER', 'INFO']

    # Check bed
    assert database.get_extra_columns(database_files.get("bed")) == ['annot1', 'annot2']

    # Check None
    assert database.get_extra_columns(None) == []

