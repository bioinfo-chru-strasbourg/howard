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


tests_folder = os.path.dirname(__file__)




    

def test_set_get_input():

    # Init files
    input_vcf = tests_folder + "/data/example.vcf.gz"
    new_input_vcf = tests_folder + "/data/example.parquet"

    # Create connection
    conn = duckdb.connect(":memory:")

    # Create object
    vcf = Variants(conn=conn, input=input_vcf)

    # Check get_input
    check_get_input = input_vcf == vcf.get_input()

    # Check new input vcf
    vcf.set_input(new_input_vcf)
    check_set_input = new_input_vcf == vcf.get_input()

    assert check_get_input and check_set_input


def test_set_get_output():

    # Init files
    input_vcf = tests_folder + "/data/example.vcf.gz"
    output_vcf = tests_folder + "/data//tmp/example.parquet"
    new_output_vcf = tests_folder + "/data/tmp/example.parquet"

    # Create connection
    conn = duckdb.connect(":memory:")

    # Create object
    vcf = Variants(conn=conn, input=input_vcf, output=output_vcf)

    # Check get_output
    check_get_output = output_vcf == vcf.get_output()

    # Check new output vcf
    vcf.set_output(new_output_vcf)
    check_set_output = new_output_vcf == vcf.get_output()

    assert check_get_output and check_set_output


def test_set_get_config():

    # Init files
    input_vcf = tests_folder + "/data/example.vcf.gz"
    input_config = { "option": "option_value" }
    new_input_config = { "option": "option_value", "new_option": "new_option_value" }

    # Create connection
    conn = duckdb.connect(":memory:")

    # Create object
    vcf = Variants(conn=conn, input=input_vcf, config=input_config)

    # Check get_input
    check_get_config = input_config == vcf.get_config()

    # Check new input vcf
    vcf.set_config(new_input_config)
    check_set_config = new_input_config == vcf.get_config()

    assert check_get_config and check_set_config


def test_set_get_param():

    # Init files
    input_vcf = tests_folder + "/data/example.vcf.gz"
    input_param = { "option": "option_value" }
    new_input_param = { "option": "option_value", "new_option": "new_option_value" }

    # Create connection
    conn = duckdb.connect(":memory:")

    # Create object
    vcf = Variants(conn=conn, input=input_vcf, param=input_param)

    # Check get_input
    check_get_param = input_param == vcf.get_param()

    # Check new input vcf
    vcf.set_param(new_input_param)
    check_set_param = new_input_param == vcf.get_param()

    assert check_get_param and check_set_param


def test_set_get_header():

    # Init files
    input_vcf = tests_folder + "/data/example.vcf.gz"

    # Create connection
    conn = duckdb.connect(":memory:")

    # Create object
    vcf = Variants(conn=conn, input=input_vcf)

    # set_header done when vcf object creation

    # Check header VCF
    header_vcf = vcf.get_header()
    check_header_vcf = header_vcf.infos != None

    # Check header List and nb
    header_list = vcf.get_header(type="list")
    check_header_list = header_list != []
    check_header_nb = len(header_list) == 53

    # check header length
    check_header_length = vcf.get_header_length() == 52

    # check get_header_columns
    header_columns = vcf.get_header_columns().strip()
    header_columns_expected = "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	sample1	sample2	sample3	sample4".strip()
    check_header_columns = header_columns == header_columns_expected
    

    # check get_header_columns_as_sql
    header_columns_as_sql = vcf.get_header_columns_as_sql().strip()
    header_columns_as_sql_expected = """ "#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","sample1","sample2","sample3","sample4" """.strip()
    check_header_columns_as_sql = header_columns_as_sql == header_columns_as_sql_expected

    # check get_header_sample_list
    header_columns_sample_list = vcf.get_header_sample_list()
    header_columns_sample_list_expected = ['sample1', 'sample2', 'sample3', 'sample4']
    check_header_columns_sample_list = header_columns_sample_list == header_columns_sample_list_expected

    
    # check full
    check_full = check_header_vcf and check_header_list and check_header_nb and check_header_length and check_header_columns and check_header_columns_as_sql and check_header_columns_sample_list
    #check_full = check_header_vcf and check_header_list and check_header_columns and check_header_columns_as_sql and check_header_columns_sample_list

    assert check_full


def test_set_get_header_in_config():

    # Init files
    input_vcf = tests_folder + "/data/example.parquet"
    input_config = { "header_file":  tests_folder + "/data/example.parquet.hdr" }

    # Create object
    vcf = Variants(input=input_vcf, config=input_config)

    # Check header VCF
    header_vcf = vcf.get_header()
    check_header_vcf = header_vcf.infos != None

    # Check header List and nb
    header_list = vcf.get_header(type="list")
    check_header_list = header_list != []
    check_header_nb = len(header_list) == 53

    # check header length
    check_header_length = vcf.get_header_length() == 52

    # check get_header_columns
    header_columns = vcf.get_header_columns().strip()
    header_columns_expected = "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	sample1	sample2	sample3	sample4".strip()
    check_header_columns = header_columns == header_columns_expected
    

    # check get_header_columns_as_sql
    header_columns_as_sql = vcf.get_header_columns_as_sql().strip()
    header_columns_as_sql_expected = """ "#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","sample1","sample2","sample3","sample4" """.strip()
    check_header_columns_as_sql = header_columns_as_sql == header_columns_as_sql_expected

    # check get_header_sample_list
    header_columns_sample_list = vcf.get_header_sample_list()
    header_columns_sample_list_expected = ['sample1', 'sample2', 'sample3', 'sample4']
    check_header_columns_sample_list = header_columns_sample_list == header_columns_sample_list_expected

    # check full
    check_full = check_header_vcf and check_header_list and check_header_nb and check_header_length and check_header_columns and check_header_columns_as_sql and check_header_columns_sample_list
    #check_full = check_header_vcf and check_header_list and check_header_columns and check_header_columns_as_sql and check_header_columns_sample_list

    assert check_full


def test_load_without_header():

    # Init files
    input_vcf = tests_folder + "/data/example.without_header.parquet"

    # Create object
    with pytest.raises(ValueError) as e:
        vcf = Variants(input=input_vcf)
    assert str(e.value) == f"No header for file {input_vcf}"


def test_read_vcf_header():

    # Init files
    input_vcf = tests_folder + "/data/example.parquet"
    vcf_header = tests_folder + "/data/example.parquet.hdr"

    # Create connection
    conn = duckdb.connect(":memory:")

    # Create object
    vcf = Variants(conn=conn, input=input_vcf)

    # Check read_vcf_header
    with open(vcf_header, 'rt') as f:
        header_list = vcf.read_vcf_header(f)
    check_header_list = header_list != []
    check_header_list_length = len(header_list) == 53

    # check full
    check_full = check_header_list and check_header_list_length
    #check_full = check_header_list

    assert check_full


def test_load_when_init():

    # Init files
    input_vcf = tests_folder + "/data/example.vcf.gz"

    # Create object
    vcf = Variants(input=input_vcf, load=True)

    # Check data loaded
    result = vcf.execute_query("SELECT count(*) AS count FROM variants")
    nb_variant_in_database = result.df()["count"][0]
    #length = result.df()["count"][0]

    expected_number_of_variants = 7

    assert nb_variant_in_database == expected_number_of_variants


def test_load_format_not_available():

    # Init files
    input_format = "unknwon"
    input_vcf = tests_folder + f"/data/example.{input_format}"
    input_config = { "header_file":  tests_folder + "/data/example.parquet.hdr" }

    # Create object
    with pytest.raises(ValueError) as e:
        vcf = Variants(input=input_vcf, config=input_config, load=True)
    assert str(e.value) == f"Input file format '{input_format}' not available"


def test_load_vcf_gz():

    # Init files
    input_vcf = tests_folder + "/data/example.vcf.gz"

    # Create object
    vcf = Variants(input=input_vcf)

    # Load data
    vcf.load_data()

    # Check data loaded
    result = vcf.execute_query("SELECT count(*) AS count FROM variants")
    nb_variant_in_database = result.df()["count"][0]
    #length = result.df()["count"][0]

    expected_number_of_variants = 7

    assert nb_variant_in_database == expected_number_of_variants


def test_load_parquet():

    # Init files
    input_vcf = tests_folder + "/data/example.parquet"

    # Create connection
    #conn = duckdb.connect(":memory:")

    # Create object
    #vcf = Variants(conn=conn, input=input_vcf)
    vcf = Variants(input=input_vcf)

    # Load data
    vcf.load_data()

    # Check data loaded
    result = vcf.execute_query("SELECT count(*) AS count FROM variants")
    nb_variant_in_database = result.df()["count"][0]
    #length = result.df()["count"][0]

    expected_number_of_variants = 7

    assert nb_variant_in_database == expected_number_of_variants


def test_load_vcf():

    # Init files
    input_vcf = tests_folder + "/data/example.vcf"

    # Create connection
    #conn = duckdb.connect(":memory:")

    # Create object
    #vcf = Variants(conn=conn, input=input_vcf)
    vcf = Variants(input=input_vcf)

    # Load data
    vcf.load_data()

    # Check data loaded
    result = vcf.execute_query("SELECT count(*) AS count FROM variants")
    nb_variant_in_database = result.df()["count"][0]
    #length = result.df()["count"][0]

    expected_number_of_variants = 7

    assert nb_variant_in_database == expected_number_of_variants


def test_load_csv():

    # Init files
    input_vcf = tests_folder + "/data/example.csv"

    # Create connection
    #conn = duckdb.connect(":memory:")

    # Create object
    #vcf = Variants(conn=conn, input=input_vcf)
    vcf = Variants(input=input_vcf)

    # Load data
    vcf.load_data()

    # Check data loaded
    result = vcf.execute_query("SELECT count(*) AS count FROM variants")
    nb_variant_in_database = result.df()["count"][0]
    #length = result.df()["count"][0]

    expected_number_of_variants = 7

    assert nb_variant_in_database == expected_number_of_variants


def test_load_tsv():

    # Init files
    input_vcf = tests_folder + "/data/example.tsv"

    # Create connection
    #conn = duckdb.connect(":memory:")

    # Create object
    #vcf = Variants(conn=conn, input=input_vcf)
    vcf = Variants(input=input_vcf)

    # Load data
    vcf.load_data()

    # Check data loaded
    result = vcf.execute_query("SELECT count(*) AS count FROM variants")
    nb_variant_in_database = result.df()["count"][0]
    #length = result.df()["count"][0]

    expected_number_of_variants = 7

    assert nb_variant_in_database == expected_number_of_variants


def test_load_psv():

    # Init files
    input_vcf = tests_folder + "/data/example.psv"

    # Create connection
    #conn = duckdb.connect(":memory:")

    # Create object
    #vcf = Variants(conn=conn, input=input_vcf)
    vcf = Variants(input=input_vcf)

    # Load data
    vcf.load_data()

    # Check data loaded
    result = vcf.execute_query("SELECT count(*) AS count FROM variants")
    nb_variant_in_database = result.df()["count"][0]
    #length = result.df()["count"][0]

    expected_number_of_variants = 7

    assert nb_variant_in_database == expected_number_of_variants


def test_load_duckdb():

    # Init files
    input_vcf = tests_folder + "/data/example.duckdb"

    # Create object
    vcf = Variants(input=input_vcf)

    # Load data
    vcf.load_data()

    # Check data loaded
    result = vcf.execute_query("SELECT count(*) AS count FROM variants")
    nb_variant_in_database = result.df()["count"][0]

    expected_number_of_variants = 7

    assert nb_variant_in_database == expected_number_of_variants


def test_get_connexion_db_memory():

    # Init files
    input_vcf = tests_folder + "/data/example.vcf.gz"
    input_config = { "connexion_type": "memory" }

    # Create object
    vcf = Variants(input=input_vcf, config=input_config)

    # get connexion_db
    connexion_db = vcf.get_connexion_db()

    assert connexion_db == ":memory:"


def test_get_connexion_db_tmpfile():

    # Init files
    input_vcf = tests_folder + "/data/example.vcf.gz"
    input_config = { "connexion_type": "tmpfile" }

    # Create object
    vcf = Variants(input=input_vcf, config=input_config)

    # get connexion_db
    connexion_db = vcf.get_connexion_db()

    assert os.path.exists(connexion_db)


def test_get_connexion_db_file():

    # Init files
    input_vcf = tests_folder + "/data/example.vcf.gz"
    input_config = { "connexion_type": "/tmp/connexion.duckdb" }

    # Create object
    vcf = Variants(input=input_vcf, config=input_config)

    # get connexion_db
    connexion_db = vcf.get_connexion_db()

    assert connexion_db == "/tmp/connexion.duckdb"


def test_get_table_variants():

    # Init files
    input_vcf = tests_folder + "/data/example.vcf.gz"

    # Create object
    vcf = Variants(input=input_vcf)

    # get connexion_db
    table_variants_select = vcf.get_table_variants(clause="select")
    table_variants_select_check = table_variants_select == "variants"
    table_variants_from = vcf.get_table_variants(clause="from")
    table_variants_from_check = table_variants_select == "variants"
    table_variants_else = vcf.get_table_variants(clause="else")
    table_variants_else_check = table_variants_select == "variants"

    assert table_variants_select_check and table_variants_from_check and table_variants_else_check


def test_get_connexion():

    # Init files
    input_vcf = tests_folder + "/data/example.vcf.gz"

    # Create object
    vcf = Variants(input=input_vcf)

    # get connexion_db
    connexion = vcf.get_connexion()
    result = connexion.query("SELECT 'pass' AS connexion")
    check_connexion = result.df()["connexion"][0] == "pass"
    
    assert check_connexion


def test_get_verbose():

    # Init files
    input_vcf = tests_folder + "/data/example.vcf.gz"
    #input_config = { "verbose": True }

    # Create object
    vcf = Variants(input=input_vcf)

    # get connexion_db
    verbose = vcf.get_verbose()
    check_verbose_false = not verbose

    # config verbose True
    input_config = { "verbose": True }

    # Create object
    vcf = Variants(input=input_vcf, config=input_config)

    # get connexion_db
    verbose = vcf.get_verbose()
    check_verbose_true = verbose
    
    assert check_verbose_false and check_verbose_true


def test_load_connexion_type_memory():

    # Init files
    input_vcf = tests_folder + "/data/example.vcf.gz"
    input_config = { "connexion_type": "memory" }

    # Create object
    vcf = Variants(input=input_vcf, config=input_config)

    # Load data
    vcf.load_data()

    # Check data loaded
    result = vcf.execute_query("SELECT count(*) AS count FROM variants")
    nb_variant_in_database = result.df()["count"][0]

    expected_number_of_variants = 7

    assert nb_variant_in_database == expected_number_of_variants


def test_load_connexion_type_tmpfile():

    # Init files
    input_vcf = tests_folder + "/data/example.vcf.gz"
    input_config = { "connexion_type": "tmpfile" }

    # Create object
    vcf = Variants(input=input_vcf, config=input_config)

    # Load data
    vcf.load_data()

    # Check data loaded
    result = vcf.execute_query("SELECT count(*) AS count FROM variants")
    nb_variant_in_database = result.df()["count"][0]

    expected_number_of_variants = 7

    assert nb_variant_in_database == expected_number_of_variants


def test_load_connexion_type_file():

    # Init files
    input_vcf = tests_folder + "/data/example.vcf.gz"
    input_config = { "connexion_type": "/tmp/output.duckdb" }

    remove_if_exists(["/tmp/output.duckdb"])

    # Create object
    vcf = Variants(input=input_vcf, config=input_config)

    # Load data
    vcf.load_data()

    # Check data loaded
    result = vcf.execute_query("SELECT count(*) AS count FROM variants")
    nb_variant_in_database = result.df()["count"][0]

    expected_number_of_variants = 7

    assert nb_variant_in_database == expected_number_of_variants


def test_load_connexion_format_sqlite():

    # Init files
    input_vcf = tests_folder + "/data/example.vcf.gz"
    input_config = { "connexion_format": "sqlite" }

    # Create object
    vcf = Variants(input=input_vcf, config=input_config)

    # Load data
    vcf.load_data()

    # Check data loaded
    result = vcf.get_query_to_df("SELECT count(*) AS count FROM variants")
    nb_variant_in_database = result["count"][0]

    expected_number_of_variants = 7

    assert nb_variant_in_database == expected_number_of_variants


def test_export_output_bcf():

    # Init files
    input_vcf = tests_folder + "/data/example.vcf.gz"
    output_vcf = "/tmp/example.bcf"

    # remove if exists
    remove_if_exists([output_vcf])

    # Create object
    vcf = Variants(input=input_vcf, output=output_vcf, load=True)

    # Check get_output
    vcf.export_output()
    assert os.path.exists(output_vcf)

    # Check get_output without header
    remove_if_exists([output_vcf])
    vcf.export_output(export_header=False)
    assert os.path.exists(output_vcf) and os.path.exists(output_vcf + ".hdr")


def test_export_output_vcf_gz():

    # Init files
    input_vcf = tests_folder + "/data/example.vcf.gz"
    output_vcf = "/tmp/example.vcf.gz"

    # remove if exists
    remove_if_exists([output_vcf])

    # Create object
    vcf = Variants(input=input_vcf, output=output_vcf, load=True)

    # Check get_output
    vcf.export_output()
    assert os.path.exists(output_vcf)

    # Check get_output without header
    remove_if_exists([output_vcf])
    vcf.export_output(export_header=False)
    assert os.path.exists(output_vcf) and os.path.exists(output_vcf + ".hdr")


def test_export_output_vcf():

    # Init files
    input_vcf = tests_folder + "/data/example.vcf.gz"
    output_vcf = "/tmp/example.vcf"

    # remove if exists
    remove_if_exists([output_vcf])

    # Create object
    vcf = Variants(input=input_vcf, output=output_vcf, load=True)

    # Check get_output
    vcf.export_output()
    assert os.path.exists(output_vcf)

    # Check get_output without header
    remove_if_exists([output_vcf])
    vcf.export_output(export_header=False)
    assert os.path.exists(output_vcf) and os.path.exists(output_vcf + ".hdr")


def test_export_output_parquet():

    # Init files
    input_vcf = tests_folder + "/data/example.vcf.gz"
    output_vcf = "/tmp/example.parquet"

    # remove if exists
    remove_if_exists([output_vcf])

    # Create object
    vcf = Variants(input=input_vcf, output=output_vcf)

    # Load data
    vcf.load_data()

    # Check get_output
    vcf.export_output()
    assert os.path.exists(output_vcf)

    # Check get_output without header
    remove_if_exists([output_vcf])
    vcf.export_output(export_header=False)
    assert os.path.exists(output_vcf) and os.path.exists(output_vcf + ".hdr")


def test_export_output_duckdb():

    # Init files
    input_vcf = tests_folder + "/data/example.vcf.gz"
    output_vcf = "/tmp/example.duckdb"

    # remove if exists
    remove_if_exists([output_vcf])

    # Create object
    vcf = Variants(input=input_vcf, output=output_vcf)

    # Load data
    vcf.load_data()

    # Check get_output
    vcf.export_output()
    assert os.path.exists(output_vcf)

    # Check get_output without header
    vcf.export_output(export_header=False)
    assert os.path.exists(output_vcf) and os.path.exists(output_vcf + ".hdr")


def test_export_output_tsv():

    # Init files
    input_vcf = tests_folder + "/data/example.vcf.gz"
    output_vcf = "/tmp/example.tsv"

    # remove if exists
    remove_if_exists([output_vcf])

    # Create object
    vcf = Variants(input=input_vcf, output=output_vcf)

    # Load data
    vcf.load_data()

    # Check get_output
    vcf.export_output()
    assert os.path.exists(output_vcf)

    # Check get_output without header
    remove_if_exists([output_vcf])
    vcf.export_output(export_header=False)
    assert os.path.exists(output_vcf) and os.path.exists(output_vcf + ".hdr")


def test_export_output_tsv_gz():

    # Init files
    input_vcf = tests_folder + "/data/example.vcf.gz"
    output_vcf = "/tmp/example.tsv.gz"

    # remove if exists
    remove_if_exists([output_vcf])

    # Create object
    vcf = Variants(input=input_vcf, output=output_vcf)

    # Load data
    vcf.load_data()

    # Check get_output
    vcf.export_output()
    assert os.path.exists(output_vcf)

    # Check get_output without header
    remove_if_exists([output_vcf])
    vcf.export_output(export_header=False)
    assert os.path.exists(output_vcf) and os.path.exists(output_vcf + ".hdr")


def test_export_output_csv():

    # Init files
    input_vcf = tests_folder + "/data/example.vcf.gz"
    output_vcf = "/tmp/example.csv"

    # remove if exists
    remove_if_exists([output_vcf])

    # Create object
    vcf = Variants(input=input_vcf, output=output_vcf)

    # Load data
    vcf.load_data()

    # Check get_output
    vcf.export_output()
    assert os.path.exists(output_vcf)

    # Check get_output without header
    remove_if_exists([output_vcf])
    vcf.export_output(export_header=False)
    assert os.path.exists(output_vcf) and os.path.exists(output_vcf + ".hdr")


def test_export_output_psv():

    # Init files
    input_vcf = tests_folder + "/data/example.vcf.gz"
    output_vcf = "/tmp/example.psv"

    # remove if exists
    remove_if_exists([output_vcf])

    # Create object
    vcf = Variants(input=input_vcf, output=output_vcf)

    # Load data
    vcf.load_data()

    # Check get_output
    vcf.export_output()
    assert os.path.exists(output_vcf)

    # Check get_output without header
    remove_if_exists([output_vcf])
    vcf.export_output(export_header=False)
    assert os.path.exists(output_vcf) and os.path.exists(output_vcf + ".hdr")


def test_export_output_vcf_explode_infos():

    # Init files
    input_vcf = tests_folder + "/data/example.vcf.gz"
    output_vcf = "/tmp/example.tsv"
    param = {
                "export_extra_infos": True
        }

    # remove if exists
    remove_if_exists([output_vcf])

    # Create object
    vcf = Variants(input=input_vcf, output=output_vcf, load=True, param=param)

    # Explode infos
    vcf.explode_infos()

    # Check get_output
    vcf.export_output()
    assert os.path.exists(output_vcf)

    # Check get_output without header
    remove_if_exists([output_vcf])
    vcf.export_output(export_header=False)
    assert os.path.exists(output_vcf) and os.path.exists(output_vcf + ".hdr")


def test_prioritization_varank():

    # Init files
    input_vcf = tests_folder + "/data/example.vcf.gz"

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
    vcf = Variants(input=input_vcf, load=True, config=config, param=param)

    # Prioritization
    vcf.prioritization()

    # Check all priorized
    result = vcf.execute_query(""" SELECT INFO FROM variants """).df()
    assert len(result) > 0


    # Check annotation1
    result = vcf.execute_query(""" SELECT 1 AS count FROM variants WHERE "#CHROM" = 'chr1' AND POS = 28736 AND REF = 'A' AND ALT = 'C' AND INFO LIKE '%PZScore_default=15%' """).df()
    assert len(result) == 1


def test_prioritization_no_profiles():

    # Init files
    input_vcf = tests_folder + "/data/example.vcf.gz"

    # Construct config dict
    config = {
        "prioritization": {
            "config_profiles": None
            }
        }
    # Create object
    vcf = Variants(input=input_vcf, load=True, config=config)

    # Prioritization fail
    with pytest.raises(ValueError) as e:
        vcf.prioritization()
    assert str(e.value) == f"NO Profiles configuration"


def test_prioritization_no_pzfields():

    # Init files
    input_vcf = tests_folder + "/data/example.vcf.gz"

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
    vcf = Variants(input=input_vcf, load=True, config=config, param=param)

    # Prioritization
    vcf.prioritization()

    # Check all priorized
    result = vcf.execute_query(""" SELECT count(*) AS count FROM variants WHERE INFO LIKE '%PZScore_default=%' """).df()
    check_priorization = False
    if len(result["count"]):
        if result["count"][0] == 0:
            check_priorization = True
    assert check_priorization


def test_annotations():

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
    vcf = Variants(conn=None, input=input_vcf, output=output_vcf, param=param, load=True)

    # Remove if output file exists
    remove_if_exists([output_vcf])

    # Annotation
    vcf.annotation()

    # check param
    param_input = vcf.get_param()
    expected_param = {'annotations': param_annotations,
                      'annotation': {
                        'parquet': {'annotations': {annotation1: {'INFO': None}}},
                        'bcftools': {'annotations': {annotation2: {'CLNSIG': 'CLNSIG_new'}, annotation3: {'symbol': 'gene'}}}}
                    }

    check_param = param_input and expected_param

    # Check annotation1
    result = vcf.execute_query("SELECT 1 AS count FROM variants WHERE \"#CHROM\" = 'chr1' AND POS = 28736 AND REF = 'A' AND ALT = 'C' AND INFO LIKE '%CLNSIG_new=%'").df()
    check_annotation1 = False
    if len(result["count"]):
        if result["count"][0] == 1:
            check_annotation1 = True

    # Check annotation2
    result = vcf.execute_query("SELECT 1 AS count FROM variants WHERE \"#CHROM\" = 'chr7' AND POS = 55249063 AND REF = 'G' AND ALT = 'A' AND INFO = 'DP=125;nci60=0.66;gene=EGFR,EGFR-AS1'").df()
    check_annotation2 = False
    if len(result["count"]):
        if result["count"][0] == 1:
            check_annotation2 = True
    
    # check full
    check_full = check_param and check_annotation1 and check_annotation2

    assert check_full


def test_annotation_parquet():

    # Init files
    input_vcf = tests_folder + "/data/example.vcf.gz"
    annotation_parquet = tests_folder + "/data/annotations/nci60.parquet"
    output_vcf = "/tmp/output.vcf.gz"

    # Construct param dict
    param = {"annotation": {"parquet": {"annotations": {annotation_parquet: {"INFO": None}}}}}

    # Create object
    vcf = Variants(conn=None, input=input_vcf, output=output_vcf, param=param, load=True)

    # Remove if output file exists
    remove_if_exists([output_vcf])

    # Annotation
    vcf.annotation()

    # query annotated variant
    result = vcf.execute_query("SELECT 1 AS count FROM variants WHERE \"#CHROM\" = 'chr7' AND POS = 55249063 AND REF = 'G' AND ALT = 'A' AND INFO = 'DP=125;nci60=0.66'")
    length = len(result.df())
    
    assert length == 1


def test_annotation_parquet_field_already_in_vcf():

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
    vcf = Variants(conn=None, input=input_vcf, output=output_vcf, param=param, load=True)

    # Remove if output file exists
    remove_if_exists([output_vcf])

    # Annotation
    vcf.annotation()

    # check param
    param_input = vcf.get_param()
    expected_param = {'annotations': param_annotations,
                      'annotation': {
                        'parquet': {'annotations': {annotation1: {"nci60": "DP"}}}
                      }
                    }

    assert param_input and expected_param

    # Check annotation not changed
    result = vcf.execute_query("SELECT 1 AS count FROM variants WHERE \"#CHROM\" = 'chr7' AND POS = 55249063 AND REF = 'G' AND ALT = 'A' AND INFO = 'DP=125'").df()
    assert len(result) == 1


def test_annotation_duckdb():

    # Init files
    input_vcf = tests_folder + "/data/example.vcf.gz"
    annotation_parquet = tests_folder + "/data/annotations/nci60.duckdb"
    output_vcf = "/tmp/output.vcf.gz"

    # Construct param dict
    param = {"annotation": {"parquet": {"annotations": {annotation_parquet: {"INFO": None}}}}}

    # Create object
    vcf = Variants(conn=None, input=input_vcf, output=output_vcf, param=param, load=True)

    # Remove if output file exists
    remove_if_exists([output_vcf])

    # Annotation
    vcf.annotation()

    # query annotated variant
    result = vcf.execute_query("SELECT 1 AS count FROM variants WHERE \"#CHROM\" = 'chr7' AND POS = 55249063 AND REF = 'G' AND ALT = 'A' AND INFO = 'DP=125;nci60=0.66'")
    length = len(result.df())
    
    assert length == 1


def test_annotation_bcftools():

    # Init files
    input_vcf = tests_folder + "/data/example.vcf.gz"
    annotation_parquet = tests_folder + "/data/annotations/nci60.vcf.gz"
    output_vcf = "/tmp/output.vcf.gz"

    # Construct param dict
    param = {"annotation": {"bcftools": {"annotations":  {annotation_parquet: {"INFO": None}}}}}

    # Create object
    vcf = Variants(conn=None, input=input_vcf, output=output_vcf, param=param, load=True)

    # Remove if output file exists
    remove_if_exists([output_vcf])

    # Annotation
    vcf.annotation()

    # query annotated variant
    result = vcf.execute_query("""SELECT 1 AS count FROM variants WHERE "#CHROM" = 'chr7' AND POS = 55249063 AND REF = 'G' AND ALT = 'A' AND INFO = 'DP=125;nci60=0.66'""")
    length = len(result.df())
    
    assert length == 1


def test_annotation_bcftools_bed():

    # Init files
    input_vcf = tests_folder + "/data/example.vcf.gz"
    annotation_parquet = tests_folder + "/data/annotations/refGene.bed.gz"
    output_vcf = "/tmp/output.vcf.gz"

    # Construct param dict
    param = {"annotation": {"bcftools": {"annotations":  {annotation_parquet: {"symbol": None}}}}}

    # Create object
    vcf = Variants(conn=None, input=input_vcf, output=output_vcf, param=param, load=True)

    # Remove if output file exists
    remove_if_exists([output_vcf])

    # Annotation
    vcf.annotation()

    # TEST
    print(vcf.execute_query("""SELECT "#CHROM", POS, REF, ALT, INFO FROM variants """).df())

    # query annotated variant
    result = vcf.execute_query("""SELECT 1 AS count FROM variants WHERE "#CHROM" = 'chr7' AND POS = 55249063 AND REF = 'G' AND ALT = 'A' AND INFO = 'DP=125;symbol=EGFR,EGFR-AS1'""")
    length = len(result.df())
    
    assert length == 1


def test_annotation_annovar():

    # Init files
    input_vcf = tests_folder + "/data/example.vcf.gz"
    annotation_annovar = "nci60"
    output_vcf = "/tmp/output.vcf.gz"

    # Construct param dict
    param = {"annotation": {"annovar": {"annotations":  {annotation_annovar: {"INFO": None}}}}}

    # Create object
    vcf = Variants(conn=None, input=input_vcf, output=output_vcf, param=param, load=True)

    # Remove if output file exists
    remove_if_exists([output_vcf])

    # Annotation
    vcf.annotation()

    # query annotated variant
    result = vcf.execute_query("""SELECT 1 AS count FROM variants WHERE "#CHROM" = 'chr7' AND POS = 55249063 AND REF = 'G' AND ALT = 'A' AND INFO = 'DP=125;nci60=0.66'""")
    length = len(result.df())
    
    assert length == 1


def test_annotation_quick_annovar():

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
    vcf = Variants(conn=None, input=input_vcf, output=output_vcf, param=param, load=True)

    # Remove if output file exists
    remove_if_exists([output_vcf])

    # Annotation
    vcf.annotation()

    # query annotated variant
    result = vcf.execute_query("""SELECT 1 AS count FROM variants WHERE "#CHROM" = 'chr7' AND POS = 55249063 AND REF = 'G' AND ALT = 'A' AND INFO = 'DP=125;nci60=0.66'""")
    length = len(result.df())
    
    assert length == 1


def test_annotation_snpeff():

    # Init files
    input_vcf = tests_folder + "/data/example.vcf.gz"
    annotation_snpeff = "snpeff"
    output_vcf = "/tmp/output.vcf.gz"

    # Construct param dict
    param = {"annotation": {"snpeff": {"options": "-lof -hgvs -oicr -noShiftHgvs -spliceSiteSize 3 "}}}

    # Create object
    vcf = Variants(conn=None, input=input_vcf, output=output_vcf, param=param, load=True)

    # Remove if output file exists
    remove_if_exists([output_vcf])

    # Annotation
    vcf.annotation()

    # query annotated variant
    result = vcf.execute_query(""" SELECT 1 AS count FROM variants """)
    length = len(result.df())
    
    assert length == 7


def test_annotation_quick_snpeff():

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
    vcf = Variants(conn=None, input=input_vcf, output=output_vcf, param=param, load=True)

    # Remove if output file exists
    remove_if_exists([output_vcf])

    # Annotation
    vcf.annotation()

    # query annotated variant
    result = vcf.execute_query(""" SELECT 1 AS count FROM variants """)
    length = len(result.df())
    
    assert length == 7


def test_annotation_bcftools_sqlite():

    # Init files
    input_vcf = tests_folder + "/data/example.vcf.gz"
    annotation_parquet = tests_folder + "/data/annotations/nci60.vcf.gz"
    output_vcf = "/tmp/output.vcf.gz"

    # Construct config dict
    config = {"connexion_format": "sqlite"}

    # Construct param dict
    param = {"annotation": {"bcftools": {"annotations":  {annotation_parquet: {"INFO": None}}}}}

    # Create object
    vcf = Variants(conn=None, input=input_vcf, output=output_vcf, param=param, config=config, load=True)

    # Remove if output file exists
    remove_if_exists([output_vcf])

    # Annotation
    vcf.annotation()

    # query annotated variant
    # result = vcf.execute_query("""SELECT 1 AS count FROM variants WHERE "#CHROM" = 'chr7' AND POS = 55249063 AND REF = 'G' AND ALT = 'A' AND INFO = 'DP=125;nci60=0.66'""")
    # length = len(result.df())

    result = vcf.get_query_to_df("""SELECT 1 AS count FROM variants WHERE "#CHROM" = 'chr7' AND POS = 55249063 AND REF = 'G' AND ALT = 'A' AND INFO = 'DP=125;nci60=0.66'""")
    length = len(result)
    
    assert length == 1


def test_explode_infos():

    # Init files
    input_vcf = tests_folder + "/data/example.vcf.gz"
    output_vcf = "/tmp/output.vcf.gz"

    # Create object
    vcf = Variants(input=input_vcf, output=output_vcf)

    # Load data
    vcf.load_data()

    # Explode infos fields
    vcf.explode_infos()

    # Remove if output file exists
    remove_if_exists([output_vcf])

    # Annotation
    vcf.annotation()

    # column to check
    column_to_check = "INFO/CLNSIG"
    value_to_check = "pathogenic"

    # check column found
    result = vcf.execute_query("SELECT * FROM variants LIMIT 0")
    assert column_to_check in [col[0] for col in result.description]

    # Check value in column
    result = vcf.execute_query(f"""SELECT "{column_to_check}" AS column_to_check FROM variants WHERE "#CHROM" = 'chr1' AND POS = 28736 AND REF = 'A' AND ALT = 'C' """)
    assert value_to_check == result.df()["column_to_check"][0]


def test_explode_infos_param_prefix():

    # Init files
    input_vcf = tests_folder + "/data/example.vcf.gz"
    output_vcf = "/tmp/output.vcf.gz"
    infos_prefix = "INFO_"
    input_param = {"explode_infos": infos_prefix}

    # Create object
    vcf = Variants(input=input_vcf, output=output_vcf, load=True, param=input_param)

    # Explode infos fields
    #vcf.explode_infos(prefix=infos_prefix)

    # Remove if output file exists
    remove_if_exists([output_vcf])

    # Annotation
    vcf.annotation()

    # column to check
    column_to_check = infos_prefix + "CLNSIG"
    value_to_check = "pathogenic"

    # check column found
    result = vcf.execute_query("SELECT * FROM variants LIMIT 0")
    assert column_to_check in [col[0] for col in result.description]

    # Check value in column
    result = vcf.execute_query(f"""SELECT "{column_to_check}" AS column_to_check FROM variants WHERE "#CHROM" = 'chr1' AND POS = 28736 AND REF = 'A' AND ALT = 'C' """)
    assert value_to_check == result.df()["column_to_check"][0]


def test_overview():

    # Init files
    input_vcf = tests_folder + "/data/example.vcf.gz"

    # Create object
    vcf = Variants(input=input_vcf)

    # Load data
    vcf.load_data()

    # Overview
    overview = vcf.get_overview()

    assert overview == None


def test_stats():

    # Init files
    input_vcf = tests_folder + "/data/example.vcf.gz"

    # Create object
    vcf = Variants(input=input_vcf, load=True)

    # Stats
    stats = vcf.get_stats()

    assert stats == None


def test_query():

    # Init files
    input_vcf = tests_folder + "/data/example.vcf.gz"
    output_query = "/tmp/output.tsv"

    # Create object
    vcf = Variants(input=input_vcf, load=True)

    # Query
    query = "SELECT 1 AS query"
    result_query = vcf.execute_query("SELECT 1 AS query")
    assert result_query.df()["query"][0] == 1

    # Query none
    result_query = vcf.execute_query(None)
    assert result_query == None

    # Remove if output file exists
    remove_if_exists([output_query])

    # Output query
    vcf.export_output(output_file=output_query, query=query, export_header=False)
    assert os.path.exists(output_query)


def test_calculation():

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
    vcf = Variants(input=input_vcf, output=output_vcf, param=input_param, load=True)

    # Annotation
    vcf.annotation()

    # Calculation
    vcf.calculation()

    # Check number of NOMEN (2)
    result = vcf.get_query_to_df(f"""SELECT INFO FROM variants WHERE INFO LIKE '%NOMEN=%' """)
    print(result)
    
    assert len(result) == 2

     # Check number of middle (7)
    result = vcf.get_query_to_df(f"""SELECT INFO FROM variants WHERE INFO LIKE '%middle=%' """)
    assert len(result) == 7


def test_calculation_snpeff_hgvs():

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
    vcf = Variants(conn=None, input=input_vcf, output=output_vcf, param=param, load=True)

    # Check if no snpeff_hgvs
    result = vcf.get_query_to_df(f""" SELECT INFO FROM variants WHERE INFO LIKE '%snpeff_hgvs=%' """)
    assert len(result) == 0

    # Remove if output file exists
    remove_if_exists([output_vcf])

    # Calculation
    vcf.calculation()

    # query annotated variant
    result = vcf.get_query_to_df(f""" SELECT * FROM variants WHERE INFO LIKE '%snpeff_hgvs=%' """)
    assert len(result) == 7
    

def test_calculation_snpeff_hgvs_no_ann():

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
    vcf = Variants(conn=None, input=input_vcf, output=output_vcf, param=param, load=True)

    # Remove if output file exists
    remove_if_exists([output_vcf])

    # Calculation
    vcf.calculation()

    # query annotated variant
    result = vcf.get_query_to_df(f""" SELECT * FROM variants WHERE INFO LIKE '%snpeff_hgvs=%' """)
    assert len(result) == 0
    

def test_calculation_snpeff_hgvs_transcripts():

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
    vcf = Variants(conn=None, input=input_vcf, output=output_vcf, param=param, load=True)

    # Remove if output file exists
    remove_if_exists([output_vcf])

    # Calculation
    vcf.calculation()

    # query annotated variant
    result = vcf.get_query_to_df(f""" SELECT * FROM variants WHERE INFO LIKE '%NOMEN=%' """)
    assert len(result) == 7

    # Check transcript priority
    result = vcf.get_query_to_df(f""" SELECT * FROM variants WHERE INFO LIKE '%NOMEN=EGFR:NM_001346897%' """)
    assert len(result) == 1


def test_calculation_snpeff_hgvs_notranscripts():

    # Init files
    input_vcf = tests_folder + "/data/example.snpeff.vcf.gz"
    transcripts_file = tests_folder + "/data/notranscripts.tsv"
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
    vcf = Variants(conn=None, input=input_vcf, output=output_vcf, param=param, load=True)

    # Remove if output file exists
    remove_if_exists([output_vcf])

    # Calculation
    with pytest.raises(ValueError) as e:
        vcf.calculation()
    assert str(e.value) == f"Transcript file '{transcripts_file}' does NOT exist"

    
def test_calculation_vartype():

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
    vcf = Variants(conn=None, input=input_vcf, output=output_vcf, param=param, load=True)

    # Remove if output file exists
    remove_if_exists([output_vcf])

    # Calculation
    vcf.calculation()

    result = vcf.get_query_to_df(f""" SELECT * FROM variants WHERE INFO LIKE '%VARTYPE=SNV%' """)
    assert len(result) == 5
    
    result = vcf.get_query_to_df(f""" SELECT * FROM variants WHERE INFO LIKE '%VARTYPE=INDEL%' """)
    assert len(result) == 1
    
    result = vcf.get_query_to_df(f""" SELECT * FROM variants WHERE INFO LIKE '%VARTYPE=MOSAIC%' """)
    assert len(result) == 1
    

def test_prioritization():

    # Init files
    input_vcf = tests_folder + "/data/example.vcf.gz"

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
    vcf = Variants(input=input_vcf, load=True, config=config, param=param)

    # Prioritization
    vcf.prioritization()

    # Check all priorized default profile
    result = vcf.execute_query("""
        SELECT * FROM variants
        WHERE INFO LIKE '%PZFlag_default=%'
          AND INFO LIKE '%PZScore_default=%'
          AND INFO LIKE '%PZComment_default=%'
          AND INFO LIKE '%PZInfos_default=%'
        """).df()
    assert len(result) == 7

    # Check all priorized GERMLINE profile
    result = vcf.execute_query("""
        SELECT * FROM variants
        WHERE INFO LIKE '%PZFlag_GERMLINE=%'
          AND INFO LIKE '%PZScore_GERMLINE=%'
          AND INFO LIKE '%PZComment_GERMLINE=%'
          AND INFO LIKE '%PZInfos_GERMLINE=%'
        """).df()
    assert len(result) == 7

    # Check all priorized default profile (as default)
    result = vcf.execute_query("""
        SELECT * FROM variants
        WHERE INFO LIKE '%PZFlag=%'
          AND INFO LIKE '%PZScore=%'
          AND INFO LIKE '%PZComment=%'
          AND INFO LIKE '%PZInfos=%'
        """).df()
    assert len(result) == 7

    # Check annotation1
    result = vcf.execute_query(""" SELECT * FROM variants WHERE "#CHROM" = 'chr1' AND POS = 28736 AND REF = 'A' AND ALT = 'C' AND INFO LIKE '%PZScore_default=15%' """).df()
    assert len(result) == 1

    # Check FILTERED
    result = vcf.execute_query(f""" SELECT INFO FROM variants WHERE INFO LIKE '%FILTERED%' """).df()
    assert len(result) == 1

    
def test_calculation_findbypipeline():

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
    vcf = Variants(conn=None, input=input_vcf, output=output_vcf, param=param, load=True)

    # Remove if output file exists
    remove_if_exists([output_vcf])

    # Calculation
    vcf.calculation()

    result = vcf.get_query_to_df(f""" SELECT INFO FROM variants WHERE INFO LIKE '%findbypipeline%' """)
    assert len(result) == 7

    result = vcf.get_query_to_df(f""" SELECT * FROM variants WHERE INFO LIKE '%findbypipeline=4/4%' """)
    assert len(result) == 1
    
    result = vcf.get_query_to_df(f""" SELECT * FROM variants WHERE INFO LIKE '%findbypipeline=3/4%' """)
    assert len(result) == 6
    
    
def test_calculation_genotype_concordance():

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
    vcf = Variants(conn=None, input=input_vcf, output=output_vcf, param=param, load=True)

    # Remove if output file exists
    remove_if_exists([output_vcf])

    # Calculation
    vcf.calculation()

    result = vcf.get_query_to_df(f""" SELECT INFO FROM variants WHERE INFO LIKE '%genotypeconcordance%' """)
    assert len(result) == 7

    result = vcf.get_query_to_df(f""" SELECT * FROM variants WHERE INFO LIKE '%genotypeconcordance=FALSE%' """)
    assert len(result) == 1
    
    result = vcf.get_query_to_df(f""" SELECT * FROM variants WHERE INFO LIKE '%genotypeconcordance=TRUE%' """)
    assert len(result) == 6


def test_calculation_barcode():

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
    vcf = Variants(conn=None, input=input_vcf, output=output_vcf, param=param, load=True)

    # Remove if output file exists
    remove_if_exists([output_vcf])

    # Calculation
    vcf.calculation()

    result = vcf.get_query_to_df(f""" SELECT INFO FROM variants """)
    print(result)

    result = vcf.get_query_to_df(f""" SELECT INFO FROM variants WHERE INFO LIKE '%barcode%' """)
    assert len(result) == 7

    result = vcf.get_query_to_df(f""" SELECT * FROM variants WHERE INFO LIKE '%barcode=1122%' """)
    assert len(result) == 1

    result = vcf.get_query_to_df(f""" SELECT * FROM variants WHERE INFO LIKE '%barcode=0111%' """)
    assert len(result) == 1
    
    result = vcf.get_query_to_df(f""" SELECT * FROM variants WHERE INFO LIKE '%barcode=1011%' """)
    assert len(result) == 4
    
    result = vcf.get_query_to_df(f""" SELECT * FROM variants WHERE INFO LIKE '%barcode=1101%' """)
    assert len(result) == 1
    
