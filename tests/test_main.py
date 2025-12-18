# -*- coding: utf-8 -*-
"""
Tests

Usage:
pytest tests/

Coverage:
coverage run -m pytest . -x -v
#coverage run -m pytest --mypy . -x -v --ignore-missing-imports
coverage report --include=howard/* -m
"""

import pytest  # type: ignore
import os
import duckdb  # type: ignore
from tempfile import TemporaryDirectory

# from howard.functions.commons import *
from howard.functions.commons import remove_if_exists
from howard.objects.variants import Variants


tests_folder = os.path.dirname(__file__)


def test_load():
    """
    This function tests that the input VCF file is correctly loaded into the Variants object
    """
    input_vcf = tests_folder + "/data/example.vcf.gz"
    vcf = Variants(input=input_vcf)
    input_vcf_test = vcf.get_input()
    assert input_vcf_test == input_vcf


@pytest.mark.parametrize(
    "input, output, header",
    [
        (tests_folder + "/data/example.vcf.gz", "output.tsv.gz", "output.tsv.gz.hdr"),
        (tests_folder + "/data/example.vcf.gz", "output.vcf", None),
        (tests_folder + "/data/example.vcf.gz", "output.vcf.gz", None),
        (tests_folder + "/data/example.vcf.gz", "output.parquet", "output.parquet.hdr"),
    ],
)
def test_export(input, output, header):
    """
    It loads a VCF file into a DuckDB database, and then exports it back to a VCF file
    """

    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        output = os.path.join(tmp_dir, output)
        vcf = Variants(input=input, output=output)
        vcf.load_data()
        remove_if_exists([output])
        vcf.export_output()
        assert os.path.exists(output)
        if header is not None:
            header = os.path.join(tmp_dir, header)
            assert os.path.exists(header)


@pytest.mark.parametrize(
    "input_vcf, clause, variants_table, config",
    [
        (tests_folder + "/data/example.vcf.gz", None, "variants", {}),
        (tests_folder + "/data/example.vcf.gz", "select", "variants", {}),
        (tests_folder + "/data/example.vcf.gz", "from", "variants as variants", {}),
        (tests_folder + "/data/example.parquet", None, "variants", {}),
        (tests_folder + "/data/example.parquet", "select", "variants", {}),
        (tests_folder + "/data/example.parquet", "from", "variants as variants", {}),
        (
            tests_folder + "/data/example.parquet",
            "from",
            f"'{tests_folder + '/data/example.parquet'}' as variants",
            {"access": "RO"},
        ),
    ],
)
def test_get_table_variants(input_vcf, clause, variants_table, config):
    vcf = Variants(input=input_vcf, config=config)
    vcf.load_data()
    if clause:
        result = vcf.get_table_variants(clause=clause)
    else:
        result = vcf.get_table_variants()
    assert variants_table == result


def test_query():
    """
    This function connects to a duckdb database, executes a query, and returns the result
    """
    conn = duckdb.connect(":memory:")
    result = conn.execute("SELECT 1 AS count").df()["count"][0]
    assert result == 1
    conn.close()
