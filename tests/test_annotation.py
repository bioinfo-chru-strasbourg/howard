# -*- coding: utf-8 -*-
"""
Tests

Usage:
pytest tests/

Coverage:
coverage run -m pytest
"""
from __future__ import division
from __future__ import print_function


import logging as log
import os
import sys
import duckdb
import re
import Bio.bgzf as bgzf
import gzip

from howard.objects.variants import Variants


tests_folder = os.path.dirname(__file__)



def remove_if_exists(filepath):
    if os.path.exists(filepath):
        os.remove(filepath)


def test_remove_if_exists():

    filename = "/tmp/output.test"

    fhandle = open(filename, 'a')
    try:
        os.utime(filename, None)
    finally:
        fhandle.close()

    created = os.path.exists(filename)

    remove_if_exists(filename)

    deleted = not os.path.exists(filename)

    assert created and deleted



def test_annotation_parquet():

    # Init files
    input_vcf = tests_folder + "/data/example.vcf.gz"
    annotation_parquet = tests_folder + "/data/annotations/nci60.parquet"
    output_vcf = "/tmp/output.vcf.gz"
    # Create connection
    conn = duckdb.connect(":memory:")
    # Construct param dict
    param = {"annotation": {"parquet": {annotation_parquet: {"INFO": None}}}}
    # Create object
    vcf = Variants(conn=conn, input=input_vcf, output=output_vcf, param=param)
    # Load data
    vcf.load_data()
    # Remove if output file exists
    remove_if_exists(output_vcf)
    # Annotation
    vcf.annotation()
    # query annotated variant
    result = vcf.execute_query("SELECT 1 AS count FROM variants WHERE \"#CHROM\" = 'chr7' AND POS = 55249063 AND REF = 'G' AND ALT = 'A'")
    length = len(result.df())
    
    assert length

