# -*- coding: utf-8 -*-
"""
Tests

Usage:
pytest tests/
"""
from __future__ import division
from __future__ import print_function

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


def test_export():
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


def test_query():
    """
    This function connects to a duckdb database, executes a query, and returns the result
    """
    conn = duckdb.connect(":memory:")
    result = conn.execute("SELECT 1 AS count").df()["count"][0]
    print(result)
    assert result == 1

if __name__ == "__main__":
    test_load()
