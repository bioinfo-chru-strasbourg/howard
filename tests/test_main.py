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

from howard.objects.variants import Variants


tests_folder = os.path.dirname(__file__)


def remove_if_exists(filepath):
    if os.path.exists(filepath):
        os.remove(filepath)


def test_load():
    """
    This function tests that the input VCF file is correctly loaded into the Variants object
    """
    input_vcf = tests_folder + "/data/example.vcf.gz"
    vcf = Variants(input=input_vcf)
    input_vcf_test = vcf.get_input()
    assert input_vcf_test == input_vcf


def test_export_vcf():
    """
    It loads a VCF file into a DuckDB database, and then exports it back to a VCF file
    """
    input_vcf = tests_folder + "/data/example.vcf.gz"
    output_vcf = "/tmp/output.vcf"
    conn = duckdb.connect(":memory:")
    vcf = Variants(conn=conn, input=input_vcf, output=output_vcf)
    vcf.load_data()
    remove_if_exists(output_vcf)
    vcf.export_output()
    assert os.path.exists(output_vcf)


def test_export_vcf_gz():
    """
    It loads a VCF file into a DuckDB database, and then exports it back to a VCF file
    """
    input_vcf = tests_folder + "/data/example.vcf.gz"
    output_vcf = "/tmp/output.vcf.gz"
    conn = duckdb.connect(":memory:")
    vcf = Variants(conn=conn, input=input_vcf, output=output_vcf)
    vcf.load_data()
    remove_if_exists(output_vcf)
    vcf.export_output()
    assert os.path.exists(output_vcf)


def test_export_parquet():
    """
    It loads a VCF file into a DuckDB database, and then exports it back to a VCF file
    """
    input_vcf = tests_folder + "/data/example.vcf.gz"
    output_vcf = "/tmp/output.parquet"
    conn = duckdb.connect(":memory:")
    vcf = Variants(conn=conn, input=input_vcf, output=output_vcf)
    vcf.load_data()
    remove_if_exists(output_vcf)
    vcf.export_output()
    assert os.path.exists(output_vcf)


def test_export_header():
    """
    It loads a VCF file into a DuckDB database, and then exports it back to a VCF file and check header file
    """

    input_vcf = tests_folder + "/data/example.vcf.gz"
    output_vcf = "/tmp/output.vcf.gz"
    output_vcf_header = "/tmp/output.vcf.gz.hdr"
    conn = duckdb.connect(":memory:")
    vcf = Variants(conn=conn, input=input_vcf, output=output_vcf)
    vcf.load_data()
    remove_if_exists(output_vcf)
    vcf.export_output()
    assert os.path.exists(output_vcf_header)


def test_query():
    """
    This function connects to a duckdb database, executes a query, and returns the result
    """
    conn = duckdb.connect(":memory:")
    result = conn.execute("SELECT 1 AS count").df()["count"][0]
    print(result)
    assert result == 1


def test_get_table_variants_select_vcf_gz():
    input_vcf = tests_folder + "/data/example.vcf.gz"
    conn = duckdb.connect(":memory:")
    vcf = Variants(conn=conn, input=input_vcf)
    vcf.load_data()
    variants_table = vcf.get_table_variants(clause="select")
    assert variants_table == "variants"


def test_get_table_variants_from_vcf_gz():
    input_vcf = tests_folder + "/data/example.vcf.gz"
    conn = duckdb.connect(":memory:")
    vcf = Variants(conn=conn, input=input_vcf)
    vcf.load_data()
    variants_table = vcf.get_table_variants(clause="from")
    assert variants_table == "variants as variants"


def test_get_table_variants_select_parquet():
    input_vcf = tests_folder + "/data/example.parquet"
    conn = duckdb.connect(":memory:")
    vcf = Variants(conn=conn, input=input_vcf)
    vcf.load_data()
    variants_table = vcf.get_table_variants(clause="select")
    assert variants_table == "variants"


def test_get_table_variants_from_parquet():
    input_vcf = tests_folder + "/data/example.parquet"
    conn = duckdb.connect(":memory:")
    vcf = Variants(conn=conn, input=input_vcf)
    vcf.load_data()
    variants_table = vcf.get_table_variants(clause="from")
    assert variants_table == "variants as variants"

def test_get_table_variants_from_parquet_ro():
    input_vcf = tests_folder + "/data/example.parquet"
    conn = duckdb.connect(":memory:")
    config = {"access": "RO"}
    vcf = Variants(conn=conn, input=input_vcf, config=config)
    vcf.load_data()
    variants_table = vcf.get_table_variants(clause="from")
    assert variants_table == f"'{input_vcf}' as variants"

