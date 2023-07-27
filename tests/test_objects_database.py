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
from howard.objects.variants import Variants
from howard.objects.database import Database
from test_needed import *



def test_database_as_conn():
    """
    This function test creation of a Database empty
    """

    with TemporaryDirectory(dir=".") as tmp_dir:

        # Create conn with variants in a table by loading a VCF with Variants object
        variants = Variants(input=database_files.get("example_vcf"), load=True)

        # Create database object with conn
        try:
            database = Database(conn=variants.get_connexion())
        except:
            assert False

        # Check get_conn
        assert database.get_conn() == variants.get_connexion()

        # Create database object using conn, variants object header, and point to table "variants"
        try:
            database = Database(database=variants.get_connexion(), header=variants.get_header(), table="variants")
        except:
            assert False

        # Check table variants is queryable
        assert database.query(query="SELECT * FROM variants")

        # Check if some variants
        assert len(database.query(query="SELECT * FROM variants"))

        # get_database_tables
        assert database.get_database_tables() == ["variants"]

        # Check get_database_basename
        assert database.get_database_basename() == None
        
        # Check is compressed
        assert not database.is_compressed()

        # Check export
        output_database = f"{tmp_dir}/output_database.vcf"
        remove_if_exists([output_database])
        try:
            assert database.export(output_database=output_database)
            assert database.query(database=output_database, query=f"""{database.get_sql_database_link(database=output_database)}""")
        except:
            assert False

        # Check export duckdb
        output_database = f"{tmp_dir}/output_database.duckdb"
        remove_if_exists([output_database])
        assert database.export(output_database=output_database)


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

    with TemporaryDirectory(dir=".") as tmp_dir:

        # Init files
        default_header_list = [
            '##fileformat=VCFv4.2',
            '#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO'
        ]
        default_header_file = f"{tmp_dir}/header_file.hdr"
        remove_if_exists(default_header_file)
        with open(default_header_file, "w") as f:
            f.write("\n".join(default_header_list))
        example_hearder_file = database_files.get("parquet") + ".hdr"
        example_hearder_inside_file = database_files.get("vcf")
        example_hearder_NO_file = f"{tmp_dir}/no_file.hdr"
        
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


def test_get_header_infos_list():
    """
    This function tests the `get_header` method of the `Database` class in Python.
    """

    # Init files

    # Create object
    database = Database()

    # Header None
    database_header_infos_list = database.get_header_infos_list()
    assert database_header_infos_list == []

    # Header parquet
    database_header_infos_list = database.get_header_infos_list(database=database_files.get("parquet"))
    assert database_header_infos_list == ["nci60"]

    # Create conn with variants in a table by loading a Parquet with Variants object
    # Header list is columns of the table
    variants = Variants(input=database_files.get("parquet"), load=True)
    database = Database(conn=variants.get_connexion())
    database_header_infos_list = database.get_header_infos_list(database=database.get_conn())
    assert database_header_infos_list == ["INFO/nci60"]


def test_get_header():
    """
    This function tests the `get_header` method of the `Database` class in Python.
    """

    with TemporaryDirectory(dir=".") as tmp_dir:

        # Init files
        default_header_list = [
            '##fileformat=VCFv4.2',
            '#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO'
        ]
        default_header_file = f"{tmp_dir}/header_file.hdr"
        remove_if_exists(default_header_file)
        with open(default_header_file, "w") as f:
            f.write("\n".join(default_header_list))
        example_hearder_file = database_files.get("parquet") + ".hdr"
        example_hearder_inside_file = database_files.get("vcf")
        example_hearder_NO_file = f"{tmp_dir}/no_file.hdr"

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

    with TemporaryDirectory(dir=".") as tmp_dir:

        # Init files
        default_header_list = [
            '##fileformat=VCFv4.2',
            '#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO'
        ]
        default_header_file = f"{tmp_dir}/header_file.hdr"
        remove_if_exists(default_header_file)
        with open(default_header_file, "w") as f:
            # for line in default_header_list:
            #     f.write(line)
            f.write("\n".join(default_header_list))
        example_hearder_file = database_files.get("parquet") + ".hdr"
        example_hearder_inside_file = database_files.get("vcf")
        example_hearder_NO_file = f"{tmp_dir}/no_file.hdr"
        
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

    with TemporaryDirectory(dir=".") as tmp_dir:

        # Init
        database_file = database_files.get("parquet")
        hearder_file = database_files.get("parquet") + ".hdr"
        database_without_header  = database_files.get("parquet_without_header")
        database_within_header  = database_files.get("vcf_without_header")
        output_header_file = f"{tmp_dir}/header_file.hdr"

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

    # Check duckdb
    assert database.get_type(database_files.get("duckdb")) == "variants"

    # Check sqlite
    assert database.get_type(database_files.get("sqlite")) == "variants"

    # Check vcf
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

    # Check parquet
    assert not database.is_vcf(database_files.get("parquet"))

    # Check duckdb
    assert not database.is_vcf(database_files.get("duckdb"))

    # Check vcf
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


def test_is_compressed():
    """
    This function test is a database is a vcf (contains all needed columns)
    """

    # Create object
    database = Database(database_files.get("vcf"))

    # Check duckdb
    assert not database.is_compressed()

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
    assert not database.is_compressed(database_files.get("tsv_alternative_columns"))
    assert not database.is_compressed(database_files.get("tsv_failed_columns"))
    assert not database.is_compressed(database_files.get("tsv_lower_columns"))

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
    assert database.get_sql_from(database_files.get("vcf")) == f"""read_csv('{database_files.get("vcf")}', names=['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO'], auto_detect=True, delim='\t')"""

    # Check vcf gz
    assert database.get_sql_from(database_files.get("vcf_gz")) == f"""read_csv('{database_files.get("vcf_gz")}', names=['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO'], auto_detect=True, delim='\t')"""

    # Check tsv
    assert database.get_sql_from(database_files.get("tsv")) == f"""read_csv('{database_files.get("tsv")}', names=['#CHROM', 'POS', 'REF', 'ALT', 'ID', 'QUAL', 'FILTER', 'INFO'], auto_detect=True, delim='\t')"""

    # Check tsv gz
    assert database.get_sql_from(database_files.get("tsv_gz")) == f"""read_csv('{database_files.get("tsv_gz")}', names=['#CHROM', 'POS', 'REF', 'ALT', 'ID', 'QUAL', 'FILTER', 'INFO'], auto_detect=True, delim='\t')"""

    # Check csv
    assert database.get_sql_from(database_files.get("csv")) == f"""read_csv('{database_files.get("csv")}', names=['#CHROM', 'POS', 'REF', 'ALT', 'ID', 'QUAL', 'FILTER', 'INFO'], auto_detect=True, delim=',')"""

    # Check tbl
    assert database.get_sql_from(database_files.get("tbl")) == f"""read_csv('{database_files.get("tbl")}', names=['#CHROM', 'POS', 'REF', 'ALT', 'ID', 'QUAL', 'FILTER', 'INFO'], auto_detect=True, delim='|')"""

    # Check bed
    assert database.get_sql_from(database_files.get("bed")) == f"""read_csv('{database_files.get("bed")}', names=['#CHROM', 'START', 'END', 'annot1', 'annot2'], auto_detect=True, delim='\t')"""

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

    # Check tbl
    database = Database()
    sql_database_link = database.get_sql_database_link(database=database_files.get("tbl"))
    sql_database_attach = database.get_sql_database_attach(database=database_files.get("tbl"), output="query")
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
    assert database.get_columns(database=database_files.get("parquet")) ==  ['#CHROM', 'POS', 'REF', 'ALT', 'ID', 'QUAL', 'FILTER', 'INFO', 'INFO/nci60']

    # Check vcf
    assert database.get_columns(database=database_files.get("vcf")) == ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO']

    # Check vcf without header
    assert database.get_columns(database=database_files.get("vcf_without_header")) == ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO']

    # Check tsv
    assert database.get_columns(database=database_files.get("tsv")) == ['#CHROM', 'POS', 'REF', 'ALT', 'ID', 'QUAL', 'FILTER', 'INFO']

    # Check csv
    assert database.get_columns(database=database_files.get("csv")) == ['#CHROM', 'POS', 'REF', 'ALT', 'ID', 'QUAL', 'FILTER', 'INFO']

    # Check tbl
    assert database.get_columns(database=database_files.get("tbl")) == ['#CHROM', 'POS', 'REF', 'ALT', 'ID', 'QUAL', 'FILTER', 'INFO']

    # Check bed
    assert database.get_columns(database=database_files.get("bed")) == ['#CHROM', 'START', 'END', 'annot1', 'annot2']

    # Check refgene
    assert database.get_columns(database=database_files.get("refgene")) == ['#CHROM', 'START', 'END', 'symbol', 'transcript', 'strand']

    # Check refgene with an extra header file
    assert database.get_columns(database=database_files.get("refgene"),header_file=database_files.get("refgene")+".hdr") == ['#CHROM', 'START', 'END', 'symbol', 'transcript', 'strand']

    # Check refgene without header
    assert database.get_columns(database=database_files.get("refgene_without_header")) == ['#CHROM', 'START', 'END', 'column3', 'column4', 'column5']

    # Check refgene gz
    assert database.get_columns(database=database_files.get("refgene_gz"),header_file=database_files.get("refgene_gz")+".hdr") == ['#CHROM', 'START', 'END', 'symbol', 'transcript', 'strand']

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

    # Check tbl
    assert database.get_extra_columns(database_files.get("tbl")) == ['ID', 'QUAL', 'FILTER', 'INFO']

    # Check bed
    assert database.get_extra_columns(database_files.get("bed")) == ['annot1', 'annot2']

    # Check None
    assert database.get_extra_columns(None) == []


def test_find_column():
    """
    This is a test function for a Python program that checks the functionality of a method called
    "find_column" in a database.
    """

    # Check duckdb
    database = Database(database=database_files.get("parquet"))

    # find column
    assert database.find_column() == "INFO"

    assert database.find_column(column="nci60") == "INFO/nci60"
    assert database.find_column(column="nci60", prefixes=["INFO/"]) == "INFO/nci60"
    assert database.find_column(column="nci60", prefixes=["PREFIX/"]) == "INFO"

    assert database.find_column(column="INFO/nci60") == "INFO/nci60"
    assert database.find_column(column="INFO/nci60", prefixes=["INFO/"]) == "INFO/nci60"
    assert database.find_column(column="INFO/nci60", prefixes=["PREFIX/"]) == "INFO/nci60"

    assert database.find_column(column="OTHER/nci60", prefixes=["PREFIX/"]) == None


    assert database.find_column(column="column") == None
    assert database.find_column(column="column", prefixes=["INFO/"]) == None
    assert database.find_column(column="column", prefixes=["PREFIX/"]) == None

    # Empty Database
    database = Database()
    assert database.find_column(database=database_files.get("parquet"), column="nci60") == "INFO/nci60"
    assert database.find_column(database=database_files.get("parquet"), column="INFO/nci60") == "INFO/nci60"
    assert database.find_column(database=database_files.get("parquet"), column="nci60", prefixes=["INFO/"]) == "INFO/nci60"
    assert database.find_column(database=database_files.get("parquet"), column="nci60", prefixes=["OTHER/"]) == "INFO"
    assert database.find_column(database=database_files.get("parquet"), column="column") == None


def test_map_columns():
    """
    The function tests the `map_columns` method of a `Database` class in Python.
    """

    # Check duckdb
    database = Database(database=database_files.get("parquet"))

    # find column
    assert database.map_columns() == {}

    assert database.map_columns(columns=["nci60"]) == {"nci60": "INFO/nci60"}
    assert database.map_columns(columns=["nci60"], prefixes=["INFO/"]) == {"nci60": "INFO/nci60"}
    assert database.map_columns(columns=["nci60"], prefixes=["PREFIX/"]) == {"nci60": "INFO"}

    assert database.map_columns(columns=["nci60", "QUAL"]) == {"nci60": "INFO/nci60", "QUAL": "QUAL"}
    assert database.map_columns(columns=["nci60", "QUAL"], prefixes=["INFO/"]) == {"nci60": "INFO/nci60", "QUAL": "QUAL"}
    assert database.map_columns(columns=["nci60", "QUAL"], prefixes=["PREFIX/"]) == {"nci60": "INFO", "QUAL": "QUAL"}


    assert database.map_columns(columns=["INFO/nci60"]) == {"INFO/nci60": "INFO/nci60"}
    assert database.map_columns(columns=["INFO/nci60"], prefixes=["INFO/"]) == {"INFO/nci60": "INFO/nci60"}
    assert database.map_columns(columns=["INFO/nci60"], prefixes=["PREFIX/"]) == {"INFO/nci60": "INFO/nci60"}

    assert database.map_columns(columns=["column"]) == {"column": None}
    assert database.map_columns(columns=["column"], prefixes=["INFO/"]) == {"column": None}
    assert database.map_columns(columns=["column"], prefixes=["PREFIX/"]) == {"column": None}

    assert database.map_columns(columns=["column", "QUAL"]) == {"column": None, "QUAL": "QUAL"}
    assert database.map_columns(columns=["column", "QUAL"], prefixes=["INFO/"]) == {"column": None, "QUAL": "QUAL"}
    assert database.map_columns(columns=["column", "QUAL"], prefixes=["PREFIX/"]) == {"column": None, "QUAL": "QUAL"}

    # Empty Database
    database = Database()
    assert database.map_columns(database=database_files.get("parquet"), columns=["nci60"]) == {"nci60": "INFO/nci60"}
    assert database.map_columns(database=database_files.get("parquet"), columns=["INFO/nci60"]) == {"INFO/nci60": "INFO/nci60"}
    assert database.map_columns(database=database_files.get("parquet"), columns=["nci60"], prefixes=["INFO/"]) == {"nci60": "INFO/nci60"}
    assert database.map_columns(database=database_files.get("parquet"), columns=["nci60"], prefixes=["OTHER/"]) == {"nci60": "INFO"}
    assert database.map_columns(database=database_files.get("parquet"), columns=["column"]) == {"column": None}

    assert database.map_columns(database=database_files.get("parquet"), columns=["nci60", "QUAL"]) == {"nci60": "INFO/nci60", "QUAL": "QUAL"}


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
        '#CHROM': None,
        'POS': None,
        'REF': None,
        'ALT': None
    }
    variants_needed_columns_only_CHROM = {
        '#CHROM': "CHROM",
        'POS': None,
        'REF': None,
        'ALT': None
    }
    variants_needed_columns_only_ALL = {
        '#CHROM': "#CHROM",
        'POS': "POS",
        'REF': "REF",
        'ALT': "ALT"
    }
    vcf_needed_columns_empty = {
        '#CHROM': None,
        'POS': None,
        'ID': None,
        'REF': None,
        'ALT': None,
        'QUAL': None,
        'FILTER': None,
        'INFO': None,
    }
    vcf_needed_columns_only_CHROM = {
        '#CHROM': "CHROM",
        'POS': None,
        'ID': None,
        'REF': None,
        'ALT': None,
        'QUAL': None,
        'FILTER': None,
        'INFO': None,
    }
    vcf_needed_columns_only_ALL = {
        '#CHROM': "#CHROM",
        'POS': "POS",
        'ID': None,
        'REF': "REF",
        'ALT': "ALT",
        'QUAL': None,
        'FILTER': None,
        'INFO': "INFO",
    }
    regions_needed_columns_empty = {
        '#CHROM': None,
        'START': None,
        'END': None
    }
    regions_needed_columns_only_CHROM = {
        '#CHROM': "CHROM",
        'START': None,
        'END': None
    }
    regions_needed_columns_only_ALL = {
        '#CHROM': "#CHROM",
        'START': "START",
        'END': "END"
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


def test_export():
    """
    The function tests the export functionality of a database for various input and output formats.
    """

    with TemporaryDirectory(dir=".") as tmp_dir:

        # No database input
        database = Database()
        output_database=f"{tmp_dir}/output_database.no_input.parquet"
        assert not database.export(output_database)

        # database input/format
        for database_input_index in ["parquet", "partition_parquet", "vcf", "vcf_gz", "tsv", "csv", "tbl", "tsv_alternative_columns", "tsv_variants", "json", "example_vcf", "bed"]:
            for database_output_format in ["duckdb", "parquet", "partition_parquet", "vcf", "vcf.gz", "tsv", "csv", "tbl", "json", "bed"]:
                parquet_partitions = None
                # specific partition_parquet
                if database_output_format in ["partition_parquet"]:
                    database_output_format = "parquet"
                    parquet_partitions = ["#CHROM"]
                input_database = database_files.get(database_input_index)
                database = Database(database_files.get(database_input_index))
                output_database=f"{tmp_dir}/output_database.{database_output_format}"
                output_header=output_database+".hdr"
                remove_if_exists([output_database,output_header])

                try:
                    assert database.export(output_database=output_database, output_header=output_header, parquet_partitions=parquet_partitions)
                    if database.get_sql_database_attach(database=output_database):
                        database.query(database=output_database, query=f"""{database.get_sql_database_attach(database=output_database)}""")
                    assert database.query(database=output_database, query=f"""{database.get_sql_database_link(database=output_database)}""")
                except:
                    assert False

