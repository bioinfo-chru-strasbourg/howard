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
    "parquet_without_header" : tests_annotations_folder + "/nci60.without_header.parquet",
    "duckdb" : tests_annotations_folder + "/nci60.duckdb",
    "duckdb_no_annotation_table" : tests_annotations_folder + "/nci60.no_annotation_table.duckdb",
    "vcf" : tests_annotations_folder + "/nci60.vcf",
    "vcf_gz" : tests_annotations_folder + "/nci60.vcf.gz",
    "vcf_without_header" : tests_annotations_folder + "/nci60.without_header.vcf",
    "vcf_gz_without_header" : tests_annotations_folder + "/nci60.without_header.vcf.gz",
    "tsv" : tests_annotations_folder + "/nci60.tsv",
    "tsv_alternative_columns" : tests_annotations_folder + "/nci60.alternative_columns.tsv",
    "tsv_failed_columns" : tests_annotations_folder + "/nci60.failed_columns.tsv",
    "tsv_lower_columns" : tests_annotations_folder + "/nci60.lower_columns.tsv",
    "tsv_without_header" : tests_annotations_folder + "/nci60.without_header.tsv",
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
    """
    This function test creation of a Database empty
    """

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
    database_file = database_files.get("parquet")

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
    database_file_basename = os.path.basename(database_file)
    databases_folders = [os.path.dirname(database_file)]

    # Create object
    database = Database()

    # Create object
    database = Database(database=database_file)
    assert database_file == database.get_database()

    # Set new Database
    database.set_database(database=new_database_file)
    assert new_database_file == database.get_database()

    # Set None database
    database.set_database(database=None)
    assert None == database.get_database()

    # Set database with base name and without databases folders
    database.set_database(database=database_file_basename)
    assert None == database.get_database()

    # Set database with base name and databases folders
    database.set_database(database=database_file_basename, databases_folders=databases_folders)
    assert database_file == database.get_database()


def test_get_header_from_list():
    """
    This function tests the `get_header_from_list` method of a `Database` class by passing different
    header lists and checking if the expected information is extracted from them.
    """

    # Init files
    default_header_list = [
        '##fileformat=VCFv4.2',
        '#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO'
    ]
    example_header_list = [
        '##fileformat=VCFv4.2',
        '##FILTER=<ID=PASS,Description="All filters passed">',
        '##fileDate=20230212',
        '##source=Example',
        '##reference=hg19',
        '##INFO=<ID=nci60,Number=.,Type=Float,Description="nci60">',
        '##contig=<ID=chr7,length=159138663>',
        '##annotation_list=nci60',
        '#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO'
    ]

    # Create object
    database = Database()

    # Header None
    database_header = database.get_header_from_list(header_list=None)
    assert list(database_header.infos) == []

    # Header default
    database_header = database.get_header_from_list(header_list=default_header_list)
    assert list(database_header.infos) == []

    # Header example
    database_header = database.get_header_from_list(header_list=example_header_list)
    assert list(database_header.infos) == ["nci60"]


def test_get_header_from_file():
    """
    This function tests the `get_header_from_file` method of the `Database` class in Python.
    """

    # Init files
    default_header_list = [
        '##fileformat=VCFv4.2',
        '#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO'
    ]
    default_header_file = "/tmp/header_file.hdr"
    remove_if_exists(default_header_file)
    with open(default_header_file, "w") as f:
        f.write("\n".join(default_header_list))
    example_hearder_file = database_files.get("parquet") + ".hdr"
    example_hearder_inside_file = database_files.get("vcf")
    example_hearder_NO_file = "/tmp/no_file.hdr"
    
    # Create object
    database = Database()

    # Header None
    database_header = database.get_header_from_file(header_file=None)
    assert list(database_header.infos) == []

    # Header default
    database_header = database.get_header_from_file(header_file=default_header_file)
    assert list(database_header.infos) == []

    # Header example
    database_header = database.get_header_from_file(header_file=example_hearder_file)
    assert list(database_header.infos) == ["nci60"]

    # Header example with header within file
    database_header = database.get_header_from_file(header_file=example_hearder_inside_file)
    assert list(database_header.infos) == ["nci60"]

    # Header example with header within file
    database_header = database.get_header_from_file(header_file=example_hearder_NO_file)
    assert list(database_header.infos) == []



def test_get_header():
    """
    This function tests the `get_header` method of the `Database` class in Python.
    """

    # Init files
    default_header_list = [
        '##fileformat=VCFv4.2',
        '#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO'
    ]
    default_header_file = "/tmp/header_file.hdr"
    remove_if_exists(default_header_file)
    with open(default_header_file, "w") as f:
        f.write("\n".join(default_header_list))
    example_hearder_file = database_files.get("parquet") + ".hdr"
    example_hearder_inside_file = database_files.get("vcf")
    example_hearder_NO_file = "/tmp/no_file.hdr"

    # Create object
    database = Database()

    # Header None
    database_header = database.get_header(header_file=None)
    assert not database_header

    # Header default
    database_header = database.get_header(header_file=default_header_file)
    assert database_header
    assert list(database_header.infos) == []

    # Header example
    database_header = database.get_header(header_file=example_hearder_file)
    assert database_header
    assert list(database_header.infos) == ["nci60"]

    # Header example with header within file
    database_header = database.get_header(header_file=example_hearder_inside_file)
    assert database_header
    assert list(database_header.infos) == ["nci60"]

    # Header example with header within file
    database_header = database.get_header(header_file=example_hearder_NO_file)
    assert database_header
    assert list(database_header.infos) == []

    # Header example with header within list default
    database_header = database.get_header(header_list=default_header_list)
    assert database_header
    assert list(database_header.infos) == []

    # Create object with header
    database = Database(database_files.get("vcf"))
    database_header = database.get_header()
    assert list(database_header.infos) == ["nci60"]


def test_read_header_file():
    """
    This function tests the `read_header_file` method of the `Database` class.
    """

    # Init files
    default_header_list = [
        '##fileformat=VCFv4.2',
        '#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO'
    ]
    default_header_file = "/tmp/header_file.hdr"
    remove_if_exists(default_header_file)
    with open(default_header_file, "w") as f:
        # for line in default_header_list:
        #     f.write(line)
        f.write("\n".join(default_header_list))
    example_hearder_file = database_files.get("parquet") + ".hdr"
    example_hearder_inside_file = database_files.get("vcf")
    example_hearder_NO_file = "/tmp/no_file.hdr"
    
    # Create object
    database = Database()

    # Header None
    database_header_list = database.read_header_file(header_file=None)
    assert list(database_header_list) == []

    # Header default
    database_header_list = database.read_header_file(header_file=default_header_file)
    assert [line.strip() for line in list(database_header_list)] == default_header_list

    # Header example
    database_header_list = database.read_header_file(header_file=example_hearder_file)
    assert len(list(database_header_list)) == 37

    # Header example with header within file
    database_header_list = database.read_header_file(header_file=example_hearder_inside_file)
    assert len(list(database_header_list)) == 37

    # Header example with header within file
    database_header_list = database.read_header_file(header_file=example_hearder_NO_file)
    assert list(database_header_list) == []


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


def test_find_database():
    """
    This function test find database file
    """

    # Init files
    database_file = database_files.get("parquet")
    database_file_basename = os.path.basename(database_file)
    databases_folders = [os.path.dirname(database_file)]

    # Create object
    database = Database()

    # Find database with existing file
    assert database.find_database(None) == None

    # Find database with existing file
    assert database.find_database(database=database_file) == database_file

    # Find database with basename and folder list
    assert database.find_database(database=database_file_basename, databases_folders=databases_folders) == database_file
    
    # Find database with basename and folder list multiple
    assert database.find_database(database=database_file_basename, databases_folders=["/not_an_existing_folder"]+databases_folders) == database_file
    assert database.find_database(database=database_file_basename, databases_folders=databases_folders+["/not_an_existing_folder"]) == database_file

    # Find database with basename and folder list without file
    assert database.find_database(database=database_file_basename) == None
    assert database.find_database(database=database_file_basename, databases_folders=["."]) == None
    assert database.find_database(database=database_file_basename, databases_folders=["/not_an_existing_folder"]) == None


def test_set_header():
    """
    This function test set header of a database, depending of its format
    """

    # Set parquet
    database = Database(database=database_files.get("parquet"))
    assert database.header
    assert list(database.header.infos) == ["nci60"]

    # Set parquet with header
    database = Database(database=database_files.get("parquet"), header_file=database_files.get("parquet")+".hdr")
    assert database.header
    assert list(database.header.infos) == ["nci60"]

    # Set parquet without header
    database = Database(database=database_files.get("parquet_without_header"))
    assert not database.header

    # Set vcf
    database = Database(database=database_files.get("vcf"))
    assert database.header
    assert list(database.header.infos) == ["nci60"]

    # Set vcf with header as extra
    database = Database(database=database_files.get("vcf"), header_file=database_files.get("vcf")+".hdr")
    assert database.header
    assert list(database.header.infos) == ["nci60"]

    # Set vcf gz
    database = Database(database=database_files.get("vcf_gz"))
    assert database.header
    assert list(database.header.infos) == ["nci60"]

    # Set vcf without header
    database = Database(database=database_files.get("vcf_without_header"))
    assert database.header
    assert list(database.header.infos) == ["nci60"]

    # Set vcf gz without header
    database = Database(database=database_files.get("vcf_gz_without_header"))
    assert database.header
    assert list(database.header.infos) == ["nci60"]

    # Set tsv without header
    database = Database(database=database_files.get("tsv_without_header"))
    assert database.header
    assert list(database.header.infos) == []

    # Set tsv with header as extra
    database = Database(database=database_files.get("tsv_without_header"), header_file=database_files.get("vcf"))
    assert database.header
    assert list(database.header.infos) == ["nci60"]

    # Set bed 
    database = Database(database=database_files.get("bed"))
    assert database.header
    assert list(database.header.infos) == ['annot1', 'annot2']

    # Set json 
    database = Database(database=database_files.get("json"))
    assert database.header
    assert list(database.header.infos) == ["nci60"]


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

    # Check json
    assert database.get_format(database_files.get("json")) == "json"
    assert database.get_format(database_files.get("json_gz")) == "json"

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

    # Check json
    assert not database.is_compressed(database_files.get("json"))
    assert database.is_compressed(database_files.get("json_gz"))

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

    # Check bed
    assert database.get_type(database_files.get("bed")) == "regions"
    assert database.get_type(database_files.get("bed_gz")) == "regions"

    # Check json
    assert database.get_type(database_files.get("json")) == "variants"
    assert database.get_type(database_files.get("json_gz")) == "variants"

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
    assert database.is_vcf(database_files.get("json"))
    assert database.is_vcf(database_files.get("json_gz"))

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

    # Check duckdb without table anntation but tables
    assert database.get_database_tables(database_files.get("duckdb_no_annotation_table")) == ["variants"]

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

    # Check duckdb without table anntation
    assert database.get_database_table(database_files.get("duckdb_no_annotation_table")) == None

    # Check parquet
    assert database.get_database_table(database_files.get("parquet")) is None

    # Check vcf
    assert database.get_database_table(database_files.get("vcf")) is None

    # Check None
    assert database.get_database_table(None) is None


def test_get_columns():
    """
    The function tests the `get_columns` method of the `Database` class for different database types and
    tables.
    """

    # Create object
    database = Database()

    # Check duckdb
    database = Database()
    assert database.get_columns(database=database_files.get("duckdb"), table="variants") == ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'INFO/nci60']

    # Check duckdb
    database = Database()
    assert database.get_columns(database=database_files.get("duckdb_no_annotation_table"), table="variants") == ['#CHROM', 'NOTPOS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'INFO/nci60']

    # Check duckdb
    database = Database()
    assert database.get_columns(database=database_files.get("duckdb")+".hdr", table="variants") == []

    # Check duckdb
    database = Database()
    assert database.get_columns(database=database_files.get("vcf"), table="variants") == ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO']

    #assert False


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

    # Check vcf gz
    assert database.get_sql_from(database_files.get("vcf_gz")) == f"""read_csv('{database_files.get("vcf_gz")}', auto_detect=True, delim='\t')"""

    # Check tsv
    assert database.get_sql_from(database_files.get("tsv")) == f"""read_csv('{database_files.get("tsv")}', auto_detect=True, delim='\t')"""

    # Check tsv gz
    assert database.get_sql_from(database_files.get("tsv_gz")) == f"""read_csv('{database_files.get("tsv_gz")}', auto_detect=True, delim='\t')"""

    # Check csv
    assert database.get_sql_from(database_files.get("csv")) == f"""read_csv('{database_files.get("csv")}', auto_detect=True, delim=',')"""

    # Check psv
    assert database.get_sql_from(database_files.get("psv")) == f"""read_csv('{database_files.get("psv")}', auto_detect=True, delim='|')"""

    # Check bed
    assert database.get_sql_from(database_files.get("bed")) == f"""read_csv('{database_files.get("bed")}', auto_detect=True, delim='\t')"""

    # Check json
    assert database.get_sql_from(database_files.get("json")) == f"""read_json('{database_files.get("json")}', auto_detect=True)"""

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

