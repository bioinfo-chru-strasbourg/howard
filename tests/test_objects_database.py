# -*- coding: utf-8 -*-
"""
Tests

Usage:
pytest tests/

Coverage:
coverage run -m pytest . -x -v
coverage report --include=howard/* -m
"""

import os
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


def test_get_header_file():
    """
    This function tests the `get_header_file()` method of the `Database` class in Python.
    """

    # Init
    database_file = database_files.get("parquet")
    hearder_file = database_files.get("parquet") + ".hdr"
    database_without_header  = database_files.get("parquet_without_header")
    database_within_header  = database_files.get("vcf_without_header")
    output_header_file = "/tmp/header_file.hdr"

    # With with database and with header
    database = Database(database=database_file, header_file=hearder_file)
    assert database.get_header_file() == hearder_file
    assert database.get_header_file(output_header_file) == output_header_file

    # Check without database and without header file
    database = Database(database=None, header_file=None)
    assert database.get_header_file() == None
    assert database.get_header_file(output_header_file) == output_header_file

    # Check without database and with header file
    database = Database(database=None, header_file=hearder_file)
    assert database.get_header_file() == hearder_file
    assert database.get_header_file(output_header_file) == output_header_file

    # Check with database and with header associated to database
    database = Database(database=database_file, header_file=None)
    assert database.get_header_file() == hearder_file
    assert database.get_header_file(output_header_file) == output_header_file

    # Check with database and without header at all
    database = Database(database=database_without_header, header_file=None)
    assert database.get_header_file() == None
    assert database.get_header_file(output_header_file) == output_header_file

    # Check with database and with header within database
    database = Database(database=database_within_header, header_file=None)
    assert database.get_header_file() == database_within_header
    assert database.get_header_file(output_header_file) == output_header_file


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
    database_format = "parquet"
    database_file_basename_without_extention_format = re.sub(f".{database_format}$", "", database_file_basename)
    databases_folders = [os.path.dirname(database_file)]
    databases_folders_with_subfolder_assembly = [os.path.dirname(databases_folders[0])]
    assembly = os.path.basename(databases_folders[0])

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

    # Find database with basename and folder list and assembly
    assert database.find_database(database=database_file_basename, databases_folders=databases_folders_with_subfolder_assembly, assembly=assembly) == database_file

    # Find database with basename without extension/format and folder list
    assert database.find_database(database=database_file_basename_without_extention_format, databases_folders=databases_folders, format=database_format) == database_file


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
    assert database.header
    assert list(database.header.infos) == ["INFO/nci60"]

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

    # get haeder and set header
    database = Database(database=database_files.get("parquet"))
    header = database.get_header()
    database_header = Database(database=database_files.get("parquet_without_header"), header=header)
    assert list(database_header.header.infos) == ["nci60"]


def test_get_format():
    """
    This function test format database
    """

    # Check format set
    database = Database(database_files.get("parquet"), format="FORMAT")
    assert database.get_format() == "FORMAT"
    assert database.get_format(database_files.get("parquet")) == "FORMAT"

    # Create object
    database = Database()

    # Check parquet
    assert database.get_format(database_files.get("parquet")) == "parquet"

    # Check duckdb
    assert database.get_format(database_files.get("duckdb")) == "duckdb"

    # Check sqlite
    assert database.get_format(database_files.get("sqlite")) == "sqlite"

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

    # Check sqlite
    assert not database.is_compressed(database_files.get("sqlite"))

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
    #assert database.get_type(database_files.get("parquet")) == "variants"

    # Check duckdb variants/regions
    assert database.get_type(database_files.get("duckdb")) == "variants"
    #assert database.get_type(database_files.get("duckdb")) == "variants"

    # Check sqlite variants/regions
    assert database.get_type(database_files.get("sqlite")) == "variants"

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

    # Check duckdb
    assert database.get_database_tables(database_files.get("sqlite")) == ["variants"]

    # Check duckdb without table anntation but tables
    assert database.get_database_tables(database_files.get("duckdb_no_annotation_table")) == ["variants"]

    # Check parquet
    assert database.get_database_tables(database_files.get("parquet")) is None

    # Check vcf
    assert database.get_database_tables(database_files.get("vcf")) is None

    # Check None
    assert database.get_database_tables(None) is None


def test_get_database_basename():
    """
    This is a test function for the `get_database_basename()` method of the `Database` class in Python.
    """

    # Create object
    database = Database(database_files.get("parquet"))
    basename = os.path.basename(database_files.get("parquet"))

    # Check basename
    assert database.get_database_basename() == basename

    # Check basename from database name
    database = Database()
    assert database.get_database_basename(database_files.get("vcf")) == os.path.basename(database_files.get("vcf"))

    # Check basename from NO database
    database = Database()
    assert database.get_database_basename() == None


def test_get_database_dirname():
    """
    This is a test function for the `get_database_dirname()` method of the `Database` class in Python.
    """

    # Create object
    database = Database(database_files.get("parquet"))
    dirname = os.path.dirname(database_files.get("parquet"))

    # Check basename
    assert database.get_database_dirname() == dirname

    # Check basename from database name
    database = Database()
    assert database.get_database_dirname(database_files.get("vcf")) == os.path.dirname(database_files.get("vcf"))

    # Check basename from NO database
    database = Database()
    assert database.get_database_dirname() == None


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

    # Check sqlite
    assert database.get_database_table(database_files.get("sqlite")) == "variants"

    # Check duckdb without table anntation
    assert database.get_database_table(database_files.get("duckdb_no_annotation_table")) == None

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

    # Check sqlite
    database_table = database.get_database_table(database_files.get("sqlite"))
    assert database.get_sql_from(database_files.get("sqlite")) == f"""(SELECT * FROM sqlite_scan('{database_files.get("sqlite")}', '{database_table}'))"""

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



def test_get_sql_database_link():
    """
    This function get sql_database_link section from a database
    """

    # Check None
    database = Database()
    sql_database_link = database.get_sql_database_link(database=None)
    sql_database_attach = database.get_sql_database_attach(database=None, output="query")
    assert not sql_database_link
    assert not sql_database_attach    

    # Check duckdb as init
    database = Database(database=database_files.get("duckdb"))
    sql_database_link = database.get_sql_database_link()
    sql_database_attach = database.get_sql_database_attach(output="query")
    assert sql_database_link
    assert sql_database_attach
    if sql_database_attach:
        database.conn.query(sql_database_attach)
    assert len(database.conn.query(f""" SELECT * FROM {sql_database_link} """).df())

    # Check parquet as init
    database = Database(database=database_files.get("parquet"))
    sql_database_link = database.get_sql_database_link()
    sql_database_attach = database.get_sql_database_attach(output="query")
    assert sql_database_link
    assert not sql_database_attach
    if sql_database_attach:
        database.conn.query(sql_database_attach)
    assert len(database.conn.query(f""" SELECT * FROM {sql_database_link} """).df())

    # Check duckdb
    database = Database()
    sql_database_link = database.get_sql_database_link(database=database_files.get("duckdb"))
    sql_database_attach = database.get_sql_database_attach(database=database_files.get("duckdb"), output="query")
    assert sql_database_link
    assert sql_database_attach
    if sql_database_attach:
        database.conn.query(sql_database_attach)
    assert len(database.conn.query(f""" SELECT * FROM {sql_database_link} """).df())

    # Check sqlite
    database = Database()
    sql_database_link = database.get_sql_database_link(database=database_files.get("sqlite"))
    sql_database_attach = database.get_sql_database_attach(database=database_files.get("sqlite"), output="query")
    assert sql_database_link
    assert sql_database_attach
    if sql_database_attach:
        database.conn.query(sql_database_attach)
    assert len(database.conn.query(f""" SELECT * FROM {sql_database_link} """).df())

    # Check parquet
    database = Database()
    sql_database_link = database.get_sql_database_link(database=database_files.get("parquet"))
    sql_database_attach = database.get_sql_database_attach(database=database_files.get("parquet"), output="query")
    assert sql_database_link
    assert not sql_database_attach
    if sql_database_attach:
        database.conn.query(sql_database_attach)
    assert len(database.conn.query(f""" SELECT * FROM {sql_database_link} """).df())

    # Check vcf
    database = Database()
    sql_database_link = database.get_sql_database_link(database=database_files.get("vcf"))
    sql_database_attach = database.get_sql_database_attach(database=database_files.get("vcf"), output="query")
    assert sql_database_link
    assert not sql_database_attach
    if sql_database_attach:
        database.conn.query(sql_database_attach)
    assert len(database.conn.query(f""" SELECT * FROM {sql_database_link} """).df())

    # Check vcf gz
    database = Database()
    sql_database_link = database.get_sql_database_link(database=database_files.get("vcf_gz"))
    sql_database_attach = database.get_sql_database_attach(database=database_files.get("vcf_gz"), output="query")
    assert sql_database_link
    assert not sql_database_attach
    if sql_database_attach:
        database.conn.query(sql_database_attach)
    assert len(database.conn.query(f""" SELECT * FROM {sql_database_link} """).df())

    # Check tsv
    database = Database()
    sql_database_link = database.get_sql_database_link(database=database_files.get("tsv"))
    sql_database_attach = database.get_sql_database_attach(database=database_files.get("tsv"), output="query")
    assert sql_database_link
    assert not sql_database_attach
    if sql_database_attach:
        database.conn.query(sql_database_attach)
    assert len(database.conn.query(f""" SELECT * FROM {sql_database_link} """).df())

    # Check csv
    database = Database()
    sql_database_link = database.get_sql_database_link(database=database_files.get("csv"))
    sql_database_attach = database.get_sql_database_attach(database=database_files.get("csv"), output="query")
    assert sql_database_link
    assert not sql_database_attach
    if sql_database_attach:
        database.conn.query(sql_database_attach)
    assert len(database.conn.query(f""" SELECT * FROM {sql_database_link} """).df())

    # Check psv
    database = Database()
    sql_database_link = database.get_sql_database_link(database=database_files.get("psv"))
    sql_database_attach = database.get_sql_database_attach(database=database_files.get("psv"), output="query")
    assert sql_database_link
    assert not sql_database_attach
    if sql_database_attach:
        database.conn.query(sql_database_attach)
    assert len(database.conn.query(f""" SELECT * FROM {sql_database_link} """).df())

    # Check json
    database = Database()
    sql_database_link = database.get_sql_database_link(database=database_files.get("json"))
    sql_database_attach = database.get_sql_database_attach(database=database_files.get("json"), output="query")
    assert sql_database_link
    assert not sql_database_attach
    if sql_database_attach:
        database.conn.query(sql_database_attach)
    assert len(database.conn.query(f""" SELECT * FROM {sql_database_link} """).df())


def test_get_columns():
    """
    The function tests the ability of a database object to retrieve columns from various file formats.
    """

    # Create object
    database = Database(database_files.get("duckdb"))

    # Check duckdb
    assert database.get_columns() == ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'INFO/nci60']
    assert database.get_columns(table="variants") == ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'INFO/nci60']

    # Create empty object
    database = Database()

    # Check duckdb
    assert database.get_columns(database=database_files.get("duckdb")) == ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'INFO/nci60']
    assert database.get_columns(database=database_files.get("duckdb"), table="variants") == ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'INFO/nci60']
    assert database.get_columns(database=database_files.get("duckdb_no_annotation_table")) == []
    assert database.get_columns(database=database_files.get("duckdb_no_annotation_table"), table="variants") == ['#CHROM', 'NOTPOS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'INFO/nci60']

    # Check sqlite
    assert database.get_columns(database=database_files.get("sqlite")) == ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'INFO/nci60']
    assert database.get_columns(database=database_files.get("sqlite"), table="variants") == ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'INFO/nci60']

    # Check parquet
    assert database.get_columns(database=database_files.get("parquet")) ==  ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'INFO/nci60']

    # Check vcf
    assert database.get_columns(database=database_files.get("vcf")) == ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO']

    # Check vcf without header
    assert database.get_columns(database=database_files.get("vcf_without_header")) == ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO']

    # Check tsv
    assert database.get_columns(database=database_files.get("tsv")) == ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO']

    # Check csv
    assert database.get_columns(database=database_files.get("csv")) == ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO']

    # Check psv
    assert database.get_columns(database=database_files.get("psv")) == ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO']

    # Check bed
    assert database.get_columns(database=database_files.get("bed")) == ['#CHROM', 'START', 'END', 'annot1', 'annot2']

    # Check None
    assert database.get_columns(database=None) == []


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

    # Check sqlite
    assert database.get_extra_columns(database_files.get("sqlite")) == ['ID', 'QUAL', 'FILTER', 'INFO', 'INFO/nci60']

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


def test_get_header_from_columns():
    """
    The function tests the `get_header_from_columns` method of the `Database` class for different database types and
    tables.
    """

    # Check duckdb
    database = Database()

    # duckdb with associated header file
    assert list(database.get_header_from_columns(database=database_files.get("duckdb")).infos) == ["INFO/nci60"]
    assert list(database.get_header_from_columns(database=database_files.get("duckdb"), header_columns=None).infos) == ["INFO/nci60"]
    assert list(database.get_header_from_columns(database=database_files.get("duckdb"), header_columns=[]).infos) == ["INFO/nci60"]

    # Parquet without header and a annotation column
    assert list(database.get_header_from_columns(database=database_files.get("parquet_without_header")).infos) == ["INFO/nci60"]

    # TSV without header and no extra annotation columns
    assert list(database.get_header_from_columns(database=database_files.get("tsv_without_header")).infos) == []

    # VCF
    assert list(database.get_header_from_columns(database=database_files.get("vcf")).infos) == []

    # VCF without external header
    assert list(database.get_header_from_columns(database=database_files.get("vcf_without_header")).infos) == []

    #assert False


def test_get_annotations():
    """
    The function tests the `get_annotations` method of the `Database` class for different database types and
    tables.
    """

    # Check duckdb
    database = Database()

    # no database
    assert database.get_annotations(database=None) == None

    # duckdb with associated header file
    assert list(database.get_annotations(database=database_files.get("duckdb"))) == ["nci60"]

    # Parquet without header and a annotation column
    assert list(database.get_annotations(database=database_files.get("parquet_without_header"))) == ["INFO/nci60"]

    # TSV without header and no extra annotation columns
    assert list(database.get_annotations(database=database_files.get("tsv_without_header"))) == []

    # VCF
    assert list(database.get_annotations(database=database_files.get("vcf"))) == ["nci60"]

    # VCF without external header
    assert list(database.get_annotations(database=database_files.get("vcf_without_header"))) == ["nci60"]


def test_query():
    """
    The function tests the query method of a Database object.
    """

    # Create object
    database = Database()
    
    # Query empty
    assert database.query() == None

    # Query database
    database = Database(database=database_files.get("parquet"))
    assert len(database.query(query=f"""SELECT * FROM {database.get_sql_from()}"""))


def test_get_needed_columns():
    """
    
    """

    # Init
    variants_needed_columns_empty = {
        'chromosome': None,
        'position': None,
        'reference': None,
        'alternative': None
    }
    variants_needed_columns_only_CHROM = {
        'chromosome': "CHROM",
        'position': None,
        'reference': None,
        'alternative': None
    }
    variants_needed_columns_only_ALL = {
        'chromosome': "#CHROM",
        'position': "POS",
        'reference': "REF",
        'alternative': "ALT"
    }
    vcf_needed_columns_empty = {
        'chromosome': None,
        'position': None,
        'reference': None,
        'alternative': None,
        'info': None,
    }
    vcf_needed_columns_only_CHROM = {
        'chromosome': "CHROM",
        'position': None,
        'reference': None,
        'alternative': None,
        'info': None,
    }
    vcf_needed_columns_only_ALL = {
        'chromosome': "#CHROM",
        'position': "POS",
        'reference': "REF",
        'alternative': "ALT",
        'info': "INFO",
    }
    regions_needed_columns_empty = {
        'chromosome': None,
        'start': None,
        'end': None
    }
    regions_needed_columns_only_CHROM = {
        'chromosome': "CHROM",
        'start': None,
        'end': None
    }
    regions_needed_columns_only_ALL = {
        'chromosome': "#CHROM",
        'start': "START",
        'end': "END"
    }

    # Create object
    database = Database()

    # empty
    assert database.get_needed_columns(database_columns = [], database_type = None) == {}

    # Variants No columns
    assert database.get_needed_columns(database_columns = [], database_type = "variants") == variants_needed_columns_empty

    # Variants Only CHROM
    assert database.get_needed_columns(database_columns = ["CHROM"], database_type = "variants") == variants_needed_columns_only_CHROM

    # Variants ALL
    assert database.get_needed_columns(database_columns = ["#CHROM", "POS", "REF", "ALT"], database_type = "variants") == variants_needed_columns_only_ALL

    # Variants INFO MORE
    assert database.get_needed_columns(database_columns = ["#CHROM", "POS", "REF", "ALT", "INFO", "MORE"], database_type = "variants") == variants_needed_columns_only_ALL

    # VCF No columns
    assert database.get_needed_columns(database_columns = [], database_type = "vcf") == vcf_needed_columns_empty

    # VCF Only CHROM
    assert database.get_needed_columns(database_columns = ["CHROM"], database_type = "vcf") == vcf_needed_columns_only_CHROM

    # VCF ALL
    assert database.get_needed_columns(database_columns = ["#CHROM", "POS", "REF", "ALT", "INFO"], database_type = "vcf") == vcf_needed_columns_only_ALL

    # VCF MORE
    assert database.get_needed_columns(database_columns = ["#CHROM", "POS", "REF", "ALT", "INFO", "MORE"], database_type = "vcf") == vcf_needed_columns_only_ALL

    # Regions No columns
    assert database.get_needed_columns(database_columns = [], database_type = "regions") == regions_needed_columns_empty

    # Regions Only CHROM
    assert database.get_needed_columns(database_columns = ["CHROM"], database_type = "regions") == regions_needed_columns_only_CHROM

    # Regions ALL
    assert database.get_needed_columns(database_columns = ["#CHROM", "START", "END"], database_type = "regions") == regions_needed_columns_only_ALL

    # Regions MORE
    assert database.get_needed_columns(database_columns = ["#CHROM", "START", "END", "MORE"], database_type = "regions") == regions_needed_columns_only_ALL
