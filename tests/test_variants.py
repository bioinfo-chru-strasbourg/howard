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
    variants = Variants(conn=conn, input=input_vcf)

    # Check get_input
    assert input_vcf == variants.get_input()

    # Check new input vcf
    variants.set_input(new_input_vcf)
    assert new_input_vcf == variants.get_input()


def test_set_get_output():

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

    # Init files
    input_vcf = tests_folder + "/data/example.without_header.parquet"

    # Create object
    with pytest.raises(ValueError) as e:
        variants = Variants(input=input_vcf)
    assert str(e.value) == f"No header for file {input_vcf}"


def test_read_vcf_header():

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

    # Init files
    input_format = "unknwon"
    input_vcf = tests_folder + f"/data/example.{input_format}"
    input_config = { "header_file":  tests_folder + "/data/example.parquet.hdr" }

    # Create object
    with pytest.raises(ValueError) as e:
        variants = Variants(input=input_vcf, config=input_config, load=True)
    assert str(e.value) == f"Input file format '{input_format}' not available"


def test_load_vcf_gz():

    # Init files
    input_vcf = tests_folder + "/data/example.vcf.gz"

    # Create object
    variants = Variants(input=input_vcf, load=True)

    # Check data loaded
    result = variants.get_query_to_df("SELECT count(*) AS count FROM variants")
    nb_variant_in_database = result["count"][0]

    expected_number_of_variants = 7

    assert nb_variant_in_database == expected_number_of_variants


def test_load_vcf_gz_no_samples():

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

    # Init files
    input_vcf = tests_folder + "/data/example.tsv"

    # Create object
    variants = Variants(input=input_vcf, load=True)

    # Check data loaded
    result = variants.get_query_to_df("SELECT count(*) AS count FROM variants")
    nb_variant_in_database = result["count"][0]

    expected_number_of_variants = 7

    assert nb_variant_in_database == expected_number_of_variants


def test_load_psv():

    # Init files
    input_vcf = tests_folder + "/data/example.psv"

    # Create object
    variants = Variants(input=input_vcf, load=True)

    # Check data loaded
    result = variants.get_query_to_df("SELECT count(*) AS count FROM variants")
    nb_variant_in_database = result["count"][0]

    expected_number_of_variants = 7

    assert nb_variant_in_database == expected_number_of_variants


def test_load_duckdb():

    # Init files
    input_vcf = tests_folder + "/data/example.duckdb"

    # Create object
    variants = Variants(input=input_vcf, load=True)

    # Check data loaded
    result = variants.get_query_to_df("SELECT count(*) AS count FROM variants")
    nb_variant_in_database = result["count"][0]

    expected_number_of_variants = 7

    assert nb_variant_in_database == expected_number_of_variants


def test_get_connexion_db_memory():

    # Init files
    input_vcf = tests_folder + "/data/example.vcf.gz"
    input_config = { "connexion_type": "memory" }

    # Create object
    variants = Variants(input=input_vcf, config=input_config)

    # get connexion_db
    connexion_db = variants.get_connexion_db()

    assert connexion_db == ":memory:"


def test_get_connexion_db_tmpfile():

    # Init files
    input_vcf = tests_folder + "/data/example.vcf.gz"
    input_config = { "connexion_type": "tmpfile" }

    # Create object
    variants = Variants(input=input_vcf, config=input_config)

    # get connexion_db
    connexion_db = variants.get_connexion_db()

    assert os.path.exists(connexion_db)


def test_get_connexion_db_file():

    # Init files
    input_vcf = tests_folder + "/data/example.vcf.gz"
    input_config = { "connexion_type": "/tmp/connexion.duckdb" }

    # Create object
    variants = Variants(input=input_vcf, config=input_config)

    # get connexion_db
    connexion_db = variants.get_connexion_db()

    assert connexion_db == "/tmp/connexion.duckdb"


def test_get_table_variants():

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

def test_export_output_bcf():

    # Init files
    input_vcf = tests_folder + "/data/example.vcf.gz"
    output_vcf = "/tmp/example.bcf"

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

    # # Check if VCF is in correct format with pyVCF
    # try:
    #     vcf.Reader(filename=output_vcf)
    # except:
    #     assert False


def test_export_output_vcf_gz():

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


def test_export_output_vcf():

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

    # Init files
    input_vcf = tests_folder + "/data/example.vcf.gz"
    output_vcf = "/tmp/example.duckdb"

    # remove if exists
    remove_if_exists([output_vcf])

    # Create object
    variants = Variants(input=input_vcf, output=output_vcf, load=True)

    # Check get_output
    variants.export_output()
    assert os.path.exists(output_vcf)

    # Check get_output without header
    variants.export_output(export_header=False)
    assert os.path.exists(output_vcf) and os.path.exists(output_vcf + ".hdr")


def test_export_output_tsv():

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


def test_export_output_psv():

    # Init files
    input_vcf = tests_folder + "/data/example.vcf.gz"
    output_vcf = "/tmp/example.psv"

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


def test_export_output_tcf_explode_infos():

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


###
### Explode
###

def test_explode_infos():

    # Init files
    input_vcf = tests_folder + "/data/example.vcf.gz"

    # Create object
    variants = Variants(input=input_vcf, load=True)

    # Explode infos fields
    variants.explode_infos()

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


def test_explode_infos_no_infos():

    # Init files
    input_vcf = tests_folder + "/data/example.no_samples.vcf.gz"
    #output_vcf = "/tmp/output.vcf.gz"

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

    # Init files
    input_vcf = tests_folder + "/data/example.vcf.gz"

    # Create object
    variants = Variants(input=input_vcf, load=True)

    # Overview
    overview = variants.get_overview()

    assert overview == None


def test_overview_no_samples():

    # Init files
    input_vcf = tests_folder + "/data/example.no_samples.vcf.gz"

    # Create object
    variants = Variants(input=input_vcf, load=True)

    # Overview
    overview = variants.get_overview()

    assert overview == None


def test_stats():

    # Init files
    input_vcf = tests_folder + "/data/example.vcf.gz"

    # Create object
    variants = Variants(input=input_vcf, load=True)

    # Stats
    stats = variants.get_stats()

    assert stats == None


def test_stats_no_samples():

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

    # Init files
    #input_vcf = tests_folder + "/data/example.vcf.gz"

    # Create object
    #vcf = Variants(input=input_vcf, load=True)
    variant = Variants()

    assert variant.get_input() == None
    assert variant.get_header() == None
    assert variant.get_header_length() == 0
    assert variant.get_header_columns() == ""
    assert variant.get_header_columns_as_list() == []
    assert variant.get_header_columns_as_sql() == ""


###
### Query
###

def test_query():

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


def test_annotation_parquet():

    # Init files
    input_vcf = tests_folder + "/data/example.vcf.gz"
    annotation_parquet = tests_folder + "/data/annotations/nci60.parquet"
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

    # Init files
    input_vcf = tests_folder + "/data/example.vcf.gz"
    annotation_parquet = tests_folder + "/data/annotations/nci60.duckdb"
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
    assert len(result) == 1

    # Check if VCF is in correct format with pyVCF
    variants.export_output()
    try:
        vcf.Reader(filename=output_vcf)
    except:
        assert False


def test_annotation_bcftools():

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

    # Init files
    input_vcf = tests_folder + "/data/example.vcf.gz"
    annotation_annovar = "nci60"
    output_vcf = "/tmp/output.vcf.gz"

    # Construct param dict
    param = {"annotation": {"annovar": {"annotations":  {annotation_annovar: {"INFO": None}}}}}

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


def test_annotation_annovar_no_samples():

    # Init files
    input_vcf = tests_folder + "/data/example.no_samples.vcf.gz"
    annotation_annovar = "nci60"
    output_vcf = "/tmp/output.vcf.gz"

    # Construct param dict
    param = {"annotation": {"annovar": {"annotations":  {annotation_annovar: {"INFO": None}}}}}

    # Create object
    variants = Variants(conn=None, input=input_vcf, output=output_vcf, param=param, load=True)

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

    # Init files
    input_vcf = tests_folder + "/data/example.vcf.gz"
    annotation_annovar = "nci60"
    output_vcf = "/tmp/output.vcf.gz"

    # Construct config dict
    config = {"connexion_format": "sqlite"}

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


def test_annotation_snpeff():

    # Init files
    input_vcf = tests_folder + "/data/example.vcf.gz"
    output_vcf = "/tmp/output.vcf.gz"

    # Construct param dict
    param = {"annotation": {"snpeff": {"options": "-lof -hgvs -oicr -noShiftHgvs -spliceSiteSize 3 "}}}

    # Create object
    variants = Variants(conn=None, input=input_vcf, output=output_vcf, param=param, load=True)

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


def test_annotation_snpeff_no_samples():

    # Init files
    input_vcf = tests_folder + "/data/example.no_samples.vcf.gz"
    output_vcf = "/tmp/output.vcf.gz"

    # Construct param dict
    param = {"annotation": {"snpeff": {"options": "-lof -hgvs -oicr -noShiftHgvs -spliceSiteSize 3 "}}}

    # Create object
    variants = Variants(conn=None, input=input_vcf, output=output_vcf, param=param, load=True)

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
    variants = Variants(conn=None, input=input_vcf, output=output_vcf, param=param, load=True)

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

    # Init files
    input_vcf = tests_folder + "/data/example.vcf.gz"
    output_vcf = "/tmp/output.vcf.gz"

    # Construct config dict
    config = {"connexion_format": "sqlite"}

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

    # Init files
    input_vcf = tests_folder + "/data/example.vcf.gz"
    annotation_parquet = tests_folder + "/data/annotations/nci60.vcf.gz"
    output_vcf = "/tmp/output.vcf.gz"

    # Construct config dict
    config = {"connexion_format": "sqlite"}

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
    

###
### Prioritization
###

def test_prioritization():

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
    assert len(result) == 7

    # Check all priorized GERMLINE profile
    result = variants.get_query_to_df("""
        SELECT * FROM variants
        WHERE INFO LIKE '%PZFlag_GERMLINE=%'
          AND INFO LIKE '%PZScore_GERMLINE=%'
          AND INFO LIKE '%PZComment_GERMLINE=%'
          AND INFO LIKE '%PZInfos_GERMLINE=%'
        """)
    assert len(result) == 7

    # Check all priorized default profile (as default)
    result = variants.get_query_to_df("""
        SELECT * FROM variants
        WHERE INFO LIKE '%PZFlag=%'
          AND INFO LIKE '%PZScore=%'
          AND INFO LIKE '%PZComment=%'
          AND INFO LIKE '%PZInfos=%'
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


def test_prioritization_varank():

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
        SELECT INFO FROM variants WHERE INFO NOT IN ('','.') AND INFO NOT NULL
        """)
    assert len(result) == 0

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
    variants = Variants(input=input_vcf, output=output_vcf, param=input_param, load=True)

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

    # Init files
    input_vcf = tests_folder + "/data/example.full.vcf"
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
    assert len(result) == 3
    
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

    # Check if VCF is in correct format with pyVCF
    remove_if_exists([output_vcf])
    variants.export_output()
    try:
        vcf.Reader(filename=output_vcf)
    except:
        assert False


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
 