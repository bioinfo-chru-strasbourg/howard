# -*- coding: utf-8 -*-
"""
Tests

Usage:
pytest tests/

Coverage:
coverage run -m pytest tests/test_commons.py -x -vv --log-cli-level=INFO --capture=tee-sys
coverage report --include=howard/* -m
"""

import logging as log
import os
import sys
from tempfile import TemporaryDirectory
import duckdb
import re
import Bio.bgzf as bgzf
import gzip
import pytest
import pandas as pd
from pandas.testing import assert_frame_equal
from unittest.mock import patch

from howard.objects.variants import Variants
from howard.functions.commons import *
from test_needed import *


def test_detect_column_type():
    """
    The function `test_detect_column_type` contains test cases for detecting the type of data in a
    column using the `detect_column_type` function.
    """

    # Test case 1: Numeric values
    column = pd.Series([1, 2, 3, 4])
    expected_output = "DOUBLE"
    assert detect_column_type(column) == expected_output

    # Test case 2: Datetime values
    column = pd.to_datetime(["2022-01-01", "2022-01-02"])
    expected_output = "DATETIME"
    assert detect_column_type(column) == expected_output

    # Test case 3: Boolean values
    column = pd.Series([True, False, True])
    expected_output = "BOOLEAN"
    assert detect_column_type(column) == expected_output

    # Test case 4: Mixed values
    column = pd.Series([1, "2022-01-01", True, "Hello"])
    expected_output = "VARCHAR"
    assert detect_column_type(column) == expected_output


def test_explode_annotation_format():

    # Test case 1: Basic input
    annotation = "A|B|C,D|E|F"
    uniquify = False
    output_format = "fields"
    prefix = "ANN_"
    header = ["Allele", "Annotation", "Annotation_Impact"]
    expected_output = "ANN_Allele=A,D;ANN_Annotation=B,E;ANN_AnnotationImpact=C,F"
    assert (
        explode_annotation_format(annotation, uniquify, output_format, prefix, header)
        == expected_output
    )

    # Test case 2: Uniquify and JSON format
    annotation = "A|B|C,D|E|F"
    uniquify = True
    output_format = "JSON"
    prefix = "ANN_"
    header = ["Allele", "Annotation", "Annotation_Impact"]
    expected_output = '{"0":{"Allele":"A","Annotation":"B","Annotation_Impact":"C"},"1":{"Allele":"D","Annotation":"E","Annotation_Impact":"F"}}'
    assert (
        explode_annotation_format(annotation, uniquify, output_format, prefix, header)
        == expected_output
    )


def test_basic_functionality():
    # Test the basic functionality with default separators and clearing
    result = params_string_to_dict("app:param1=value1:param2=value2+value3")
    expected = {"param1": "value1", "param2": "value2,value3"}
    assert result == expected


def test_basic_without_header():
    # Test the basic functionality with default separators and clearing, without header
    result = params_string_to_dict("param1=value1:param2=value2+value3", header=False)
    expected = {"param1": "value1", "param2": "value2,value3"}
    assert result == expected


def test_custom_separators():
    # Test with custom separators
    result = params_string_to_dict(
        "app;param1=value1;param2=value2 value3",
        param_sep=";",
        val_clear={"+": ",", " ": "_"},
    )
    expected = {"param1": "value1", "param2": "value2_value3"}
    assert result == expected


def test_empty_values():
    # Test with empty values
    result = params_string_to_dict("app:param1=:param2=")
    expected = {"param1": None, "param2": None}
    log.debug(f"result={result}")
    log.debug(f"expected={expected}")
    assert result == expected


def test_no_clearing():
    # Test with no clearing of values
    result = params_string_to_dict(
        "app:param1=value1+value2:param2=value3 value4", val_clear={}
    )
    expected = {"param1": "value1+value2", "param2": "value3 value4"}
    assert result == expected


def test_invalid_input():
    # Test with invalid input
    result = params_string_to_dict("=value1:param2=value2")
    expected = {"param2": "value2"}
    assert result == expected


def test_help():

    # Init
    help_json_file = os.path.join(folder_main, "docs", "json", "help.param.json")

    # Help JSON to MD
    help_content_md = help_generation_from_json(
        help_json_file=help_json_file,
        output_type="markdown",
        title="test",
        code_type="json",
    )
    assert help_content_md != ""

    # Help JSON to HTML
    help_content_html = help_generation_from_json(
        help_json_file=help_json_file,
        output_type="html",
        title="test",
        code_type="json",
    )
    assert help_content_html != ""

    # Help
    from howard.tools.tools import arguments, commands_arguments, shared_arguments

    setup_cfg = f"{main_folder}/../../setup.cfg"
    arguments_dict = {
        "arguments": arguments,
        "commands_arguments": commands_arguments,
        "shared_arguments": shared_arguments,
    }
    help_content = help_generation(
        arguments_dict=arguments_dict, setup=setup_cfg, output_type="markdown"
    )
    assert help_content != ""

    # Gooey argument
    argument_gooey = get_argument_gooey(
        arguments=arguments, arg=list(arguments.keys())[0]
    )
    assert argument_gooey != ""


def test_identical_with_identical_files():
    with tempfile.NamedTemporaryFile(
        mode="w", delete=False
    ) as f1, tempfile.NamedTemporaryFile(mode="w", delete=False) as f2:
        f1.write("## Header\n")
        f1.write("Data Line 1\n")
        f1.write("Data Line 2\n")
        f2.write("## Header\n")
        f2.write("Data Line 1\n")
        f2.write("Data Line 2\n")
        f1.close()
        f2.close()

        assert identical([f1.name, f2.name])

        os.unlink(f1.name)
        os.unlink(f2.name)


def test_identical_with_different_header():
    with tempfile.NamedTemporaryFile(
        mode="w", delete=False
    ) as f1, tempfile.NamedTemporaryFile(mode="w", delete=False) as f2:
        f1.write("## Header\n")
        f1.write("Data Line 1\n")
        f1.write("Data Line 2\n")
        f2.write("## Another Header\n")  # Different
        f2.write("Data Line 1\n")
        f2.write("Data Line 2\n")
        f1.close()
        f2.close()

        assert identical([f1.name, f2.name])

        os.unlink(f1.name)
        os.unlink(f2.name)


def test_identical_with_different_content():
    with tempfile.NamedTemporaryFile(
        mode="w", delete=False
    ) as f1, tempfile.NamedTemporaryFile(mode="w", delete=False) as f2:
        f1.write("## Header\n")
        f1.write("Data Line 1\n")
        f1.write("Data Line 2\n")
        f2.write("## Header\n")
        f2.write("Data Line 1\n")
        f2.write("Data Line 3\n")  # Different
        f1.close()
        f2.close()

        assert not identical([f1.name, f2.name])

        os.unlink(f1.name)
        os.unlink(f2.name)


def test_transcripts_file_to_df():
    """
    The `test_transcripts_file_to_df` function tests the functionality of the `transcripts_file_to_df`
    function with different input files.
    """

    # Input transcript file
    transcripts_file = tests_data_folder + "/transcripts.tsv"
    transcripts_file_empty = tests_data_folder + "/transcripts.empty.tsv"
    transcripts_file_with_header = tests_data_folder + "/transcripts.with_header.tsv"
    transcripts_file_with_comments = (
        tests_data_folder + "/transcripts.with_comments.tsv"
    )

    df = transcripts_file_to_df(transcripts_file=transcripts_file)
    assert list(df.columns) == ["transcript", "gene"]
    assert len(df) == 4

    df = transcripts_file_to_df(transcripts_file=transcripts_file_empty)
    assert list(df.columns) == ["transcript", "gene"]
    assert len(df) == 0

    df = transcripts_file_to_df(transcripts_file="file_not_exist")
    log.debug(f"df={df}")
    assert list(df.columns) == ["transcript", "gene"]
    assert len(df) == 0

    df = transcripts_file_to_df(transcripts_file=transcripts_file_with_header)
    assert list(df.columns) == ["transcript", "gene"]
    assert len(df) == 4

    df = transcripts_file_to_df(transcripts_file=transcripts_file_with_comments)
    assert list(df.columns) == ["transcript", "gene"]
    assert len(df) == 4


def test_get_bin():
    """
    The `test_get_bin` function tests the `get_bin` function with different configurations and
    inputs.
    """

    # Test java
    config = {"tools": {"java": {"bin": "java"}}}
    java_bin = get_bin(
        tool="java",
        bin="java",
        bin_type="bin",
        config=config,
        default_folder="/usr/bin",
    )
    assert java_bin == "/usr/bin/java"

    # Config
    tool_name = "snpeff"
    tool_type = "bin"
    tool_bin = "snpEff.jar"
    tool_bin_path = tests_config.get("tools").get("snpeff")
    config = {"tools": {tool_name: {tool_type: tool_bin_path}}}

    # Test with config and tool name
    assert get_bin(tool=tool_name, config=config) == tool_bin_path

    # Test with config and tool bin
    assert os.path.basename(get_bin(bin=tool_bin, config=config)) == tool_bin

    # Test with config and tool name but config bin does NOT exists
    config_bad = {"tools": {tool_name: "~/howard/tools/snpeff/5.1d/bin/NOT_snpEff.jar"}}
    assert get_bin(tool=tool_name, config=config_bad) == None

    # Test with config with tool type, and tool name but config bin does NOT exists
    config_bad = {"tools": {tool_name: {"jar": tool_bin_path}}}
    assert get_bin(tool=tool_name, config=config_bad) == tool_bin_path

    # Test with config and tool name and tool bin, but config bin does NOT exists, and search with tool bin in default tool folder
    assert (
        os.path.basename(get_bin(tool=tool_name, bin=tool_bin, config=config))
        == tool_bin
    )

    # Test with no config and tool name
    assert get_bin(tool=tool_name) == tool_bin_path

    # Test with no config and tool bin
    assert os.path.basename(get_bin(bin=tool_bin)) == tool_bin

    # Test with no config
    assert get_bin() == None

    # Test with no config, tool bin and default tool folder
    assert (
        os.path.basename(
            get_bin(bin=tool_bin, default_folder=full_path("~/howard/tools"))
        )
        == tool_bin
    )


def test_get_bin_command():
    """ """

    # Test docker
    snpeff_bin_path = tests_config.get("tools").get("snpeff")
    config = {
        "tools": {
            "docker": {"bin": "docker"},
            "java": {"bin": "java"},
            "snpeff": {"jar": snpeff_bin_path},
            "bcftools": {
                "bin": "bcftools",
                "docker": {
                    "image": "howard:0.11.0",
                    "entrypoint": "bcftools",
                    "options": None,
                    "command": None,
                },
            },
        },
        "threads": 12,
        "memory": "40g",
        "tmp": "/tmp",
    }
    param = {
        "threads": 2,
        "memory": "16g",
        "tmp": "/tmp/howard",
    }

    # Test command bcftools found with bin
    tool_command = get_bin_command(tool="bcftools", config=config, param=param)
    assert os.path.basename(tool_command) == "bcftools"

    # Test command bcftools found with docker (specified)
    tool_command = get_bin_command(
        tool="bcftools", bin_type="docker", config=config, param=param
    )
    assert (
        "run  --rm  -v /tmp/howard:/tmp/howard  --cpus=2  --memory=16g  --entrypoint='bcftools'  howard:0.11.0 "
        in tool_command
    )

    # Test command bcftools found with docker with added options
    tool_command = get_bin_command(
        tool="bcftools",
        bin_type="docker",
        config=config,
        param=param,
        add_options="-v /host/path/to/mount:/inner/path_to/mount",
    )
    assert (
        "run  --rm  -v /tmp/howard:/tmp/howard  --cpus=2  --memory=16g  --entrypoint='bcftools'  -v /host/path/to/mount:/inner/path_to/mount  howard:0.11.0 "
        in tool_command
    )

    # Test command snpeff found with java/jar (added options for java optional)
    tool_command = get_bin_command(tool="snpeff", config=config, param=param)
    assert tool_command.endswith("snpEff.jar")
    assert "java  -XX:ParallelGCThreads=2  -XX:MaxHeapSize=16G  -jar " in tool_command


def test_download_file():
    """
    The `test_download_file` function downloads a file from a specified URL using the Aria download
    manager and checks if the downloaded file matches the expected file.
    """

    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        url = "https://raw.githubusercontent.com/bioinfo-chru-strasbourg/howard/master/README.md"
        dest = os.path.join(tmp_dir, "file")
        threads = 2

        log.debug("Download by Aria")
        remove_if_exists([dest])
        download_file(url=url, dest_file_path=dest)
        assert os.path.exists(dest)

        log.debug("Download by Aria with threads")
        remove_if_exists([dest])
        download_file(url=url, dest_file_path=dest, threads=threads)
        assert os.path.exists(dest)

        log.debug("Download by Request (no aria)")
        remove_if_exists([dest])
        download_file(url=url, dest_file_path=dest, try_aria=False)
        assert os.path.exists(dest)

        log.debug("Download by Request (no aria, with chunk)")
        remove_if_exists([dest])
        download_file(url=url, dest_file_path=dest, try_aria=False, chunk_size=1024)
        assert os.path.exists(dest)


def test_get_compression_type():
    """
    The function `test_get_compression_type` tests the `get_compression_type` function with different
    file types and asserts that the returned compression type matches the expected value.
    """

    assert get_compression_type(database_files.get("vcf")) == "none"

    assert get_compression_type(database_files.get("vcf_gz")) == "bgzip"

    assert get_compression_type(database_files.get("example_vcf_gzip")) == "gzip"

    assert get_compression_type(database_files.get("bcf")) == "unknown"

    assert get_compression_type(database_files.get("parquet")) == "none"


def test_concat_and_compress_files_vcf():
    """
    The function `test_concat_and_compress_files_vcf` tests the `concat_and_compress_files` function
    with different input and output file types.
    """

    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # Input file as plain text, output as plain text
        input_files = [database_files.get("vcf")]
        output_file = os.path.join(tmp_dir, "test1")
        compression_type = "none"
        assert concat_and_compress_files(
            input_files=input_files,
            output_file=output_file,
            compression_type=compression_type,
        )

        # Input file as plain text, output as bgzip
        input_files = [database_files.get("vcf")]
        output_file = os.path.join(tmp_dir, "test2")
        compression_type = "bgzip"
        assert concat_and_compress_files(
            input_files=input_files,
            output_file=output_file,
            compression_type=compression_type,
        )

        # Input file as plain text, output as gzip
        input_files = [database_files.get("vcf")]
        output_file = os.path.join(tmp_dir, "test3")
        compression_type = "gzip"
        assert concat_and_compress_files(
            input_files=input_files,
            output_file=output_file,
            compression_type=compression_type,
        )


def test_concat_and_compress_files_vcf_gz():
    """
    The function `test_concat_and_compress_files_vcf_gz` tests the `concat_and_compress_files` function
    with different compression types.
    """

    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # Input file as bgzip, output as plain text
        input_files = [database_files.get("vcf_gz")]
        output_file = os.path.join(tmp_dir, "test4")
        compression_type = "none"
        assert concat_and_compress_files(
            input_files=input_files,
            output_file=output_file,
            compression_type=compression_type,
        )

        # Input file as bgzip, output as bgzip
        input_files = [database_files.get("vcf_gz")]
        output_file = os.path.join(tmp_dir, "test5")
        compression_type = "bgzip"
        assert concat_and_compress_files(
            input_files=input_files,
            output_file=output_file,
            compression_type=compression_type,
        )

        # Input file as bgzip, output as gzip
        input_files = [database_files.get("vcf_gz")]
        output_file = os.path.join(tmp_dir, "test6")
        compression_type = "gzip"
        assert concat_and_compress_files(
            input_files=input_files,
            output_file=output_file,
            compression_type=compression_type,
        )


def test_concat_and_compress_files_vcf_gzip():
    """
    The function `test_concat_and_compress_files_vcf_gzip` tests the `concat_and_compress_files`
    function with different compression types.
    """

    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # Input file as gzip, output as plain text
        input_files = [database_files.get("example_vcf_gzip")]
        output_file = os.path.join(tmp_dir, "test7")
        compression_type = "none"
        assert concat_and_compress_files(
            input_files=input_files,
            output_file=output_file,
            compression_type=compression_type,
        )

        # Input file as gzip, output as bgzip
        input_files = [database_files.get("example_vcf_gzip")]
        output_file = os.path.join(tmp_dir, "test8")
        compression_type = "bgzip"
        assert concat_and_compress_files(
            input_files=input_files,
            output_file=output_file,
            compression_type=compression_type,
        )

        # Input file as gzip, output as gzip
        input_files = [database_files.get("example_vcf_gzip")]
        output_file = os.path.join(tmp_dir, "test9")
        compression_type = "gzip"
        assert concat_and_compress_files(
            input_files=input_files,
            output_file=output_file,
            compression_type=compression_type,
        )


def test_concat_and_compress_files_sort_and_index():
    """
    The function `test_concat_and_compress_files_sort_and_index` tests the `concat_and_compress_files`
    function by compressing, sorting, and indexing input files.
    """

    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # Compress and sort and index
        input_files = [database_files.get("example_vcf_gzip")]
        output_file = os.path.join(tmp_dir, "test10")
        compression_type = "bgzip"
        assert concat_and_compress_files(
            input_files=input_files,
            output_file=output_file,
            compression_type=compression_type,
            sort=True,
            index=True,
        )

        # Compress and sort and index gzip
        input_files = [database_files.get("example_vcf_gzip")]
        output_file = os.path.join(tmp_dir, "test11")
        compression_type = "gzip"
        assert concat_and_compress_files(
            input_files=input_files,
            output_file=output_file,
            compression_type=compression_type,
            sort=True,
            index=True,
        )

        # Compress and sort and index none
        input_files = [database_files.get("example_vcf_gzip")]
        output_file = os.path.join(tmp_dir, "test12")
        compression_type = "none"
        assert concat_and_compress_files(
            input_files=input_files,
            output_file=output_file,
            compression_type=compression_type,
            sort=True,
            index=True,
        )


def test_concat_and_compress_files_multiple_inputs():
    """
    The function `test_concat_and_compress_files_multiple_inputs()` tests the
    `concat_and_compress_files()` function with multiple input files and different compression types.
    """

    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # Multiple Input files, output as none
        input_files = [
            database_files.get("vcf"),
            database_files.get("example_vcf_gzip"),
            database_files.get("vcf_gz"),
        ]
        output_file = os.path.join(tmp_dir, "test13")
        compression_type = "none"
        assert concat_and_compress_files(
            input_files=input_files,
            output_file=output_file,
            compression_type=compression_type,
        )

        # Multiple Input files, output as gzip
        input_files = [
            database_files.get("vcf"),
            database_files.get("example_vcf_gzip"),
            database_files.get("vcf_gz"),
        ]
        output_file = os.path.join(tmp_dir, "test14")
        compression_type = "gzip"
        assert concat_and_compress_files(
            input_files=input_files,
            output_file=output_file,
            compression_type=compression_type,
        )

        # Multiple Input files, output as bgzip
        input_files = [
            database_files.get("vcf"),
            database_files.get("example_vcf_gzip"),
            database_files.get("vcf_gz"),
        ]
        output_file = os.path.join(tmp_dir, "test15")
        compression_type = "bgzip"
        assert concat_and_compress_files(
            input_files=input_files,
            output_file=output_file,
            compression_type=compression_type,
        )


def test_get_file_compressed():
    """
    This function tests whether a given file is compressed or not, based on its file extension.
    """

    # Test pour un fichier compressé .gz
    assert get_file_compressed("testfile.gz") == True

    # Test pour un fichier compressé .bcf
    assert get_file_compressed("testfile.bcf") == False

    # Test pour un fichier non compressé .vcf
    assert get_file_compressed("testfile.vcf") == False

    # Test pour un fichier non compressé .tsv
    assert get_file_compressed("testfile.tsv") == False

    # Test pour un fichier compressé .csv.gz
    assert get_file_compressed("testfile.csv.gz") == True

    # Test pour un fichier compressé .parquet
    assert get_file_compressed("testfile.parquet") == False

    # Test pour un fichier compressé .parquet.gz
    assert get_file_compressed("testfile.parquet.gz") == True


def test_get_file_format():
    """
    The function tests the `get_file_format` function for various file formats and asserts that the
    output is correct.
    """

    # Test pour un fichier .vcf
    assert get_file_format("testfile.vcf") == "vcf"

    # Test pour un fichier .tsv.gz
    assert get_file_format("testfile.tsv.gz") == "tsv"

    # Test pour un fichier .csv
    assert get_file_format("testfile.csv") == "csv"

    # Test pour un fichier .bcf.gz
    assert get_file_format("testfile.bcf") == "bcf"

    # Test pour un fichier .bcf.gz
    assert get_file_format("testfile.bcf.gz") == "bcf"

    # Test pour un fichier .parquet
    assert get_file_format("testfile.parquet") == "parquet"

    # Test pour un fichier .txt.gz
    assert get_file_format("testfile.txt.gz") == "txt"


def test_remove_if_exists():
    """
    This function tests the functionality of the "remove_if_exists" function by creating a file,
    attempting to delete it using "remove_if_exists", and checking if the file was successfully deleted.
    """

    # filename
    filename = "/tmp/output.test"

    # create filename
    fhandle = open(filename, "a")
    try:
        os.utime(filename, None)
    finally:
        fhandle.close()
    created = os.path.exists(filename)

    # remove filename
    remove_if_exists([filename])

    # check delete
    deleted = not os.path.exists(filename)

    assert created and deleted


def test_remove_if_exists_complete():
    """
    This is a test function for the "remove_if_exists" function in Python that tests its ability to
    remove existing files and handle non-existent files.
    """

    # Crée un fichier pour tester la fonction
    test_path = "test_dir"
    os.makedirs(test_path, exist_ok=True)
    test_file_path = os.path.join(test_path, "test_file.txt")
    with open(test_file_path, "w") as f:
        f.write("Test file content")

    # Teste la fonction avec un fichier qui existe
    remove_if_exists(test_file_path)
    assert not os.path.exists(test_file_path)

    # Teste la fonction avec un fichier qui n'existe pas
    remove_if_exists("nonexistent_file.txt")

    # Teste la fonction avec une liste de fichiers dont certains existent
    filepaths = [test_file_path, "nonexistent_file.txt"]
    remove_if_exists(filepaths)
    assert not os.path.exists(test_file_path)


def test_set_log_level():
    """
    This is a unit test for a function called "set_log_level" that sets the verbosity level and checks
    if it was set correctly.
    """

    # define verbosity
    verbosity = "info"

    # set verbosity
    result_verbosity = set_log_level(verbosity)

    assert verbosity == result_verbosity


def test_set_log_level_error():
    """
    This function tests if an error is raised when an unknown verbosity level is passed to the
    set_log_level function in Python.
    """

    # define verbosity
    verbosity = "not_a_level"

    # check verbosity
    with pytest.raises(ValueError) as e:
        set_log_level(verbosity)
    assert str(e.value) == "Unknown verbosity level:" + verbosity


def test_split_interval_either():
    """
    This is a test function that checks if an error is raised when neither step nor ncuts is provided in
    the split_interval function.
    """

    start = 0
    end = 1000
    step = None
    ncuts = None
    with pytest.raises(ValueError) as e:
        split_interval(start, end, step=step, ncuts=ncuts)
    assert str(e.value) == "Either step or ncuts must be provided"


def test_split_interval_only():
    """
    This is a test function that checks if an error is raised when both step and ncuts are provided as
    arguments in the split_interval function.
    """

    start = 0
    end = 1000
    step = 100
    ncuts = 4
    with pytest.raises(ValueError) as e:
        split_interval(start, end, step=step, ncuts=ncuts)
    assert str(e.value) == "Only one of step or ncuts must be provided"


def test_split_interval_step():
    """
    This is a test function that checks if the output of the split_interval function matches the
    expected output for a given set of input parameters.
    """

    start = 0
    end = 1000
    step = 100
    ncuts = None
    split = split_interval(start, end, step=step, ncuts=ncuts)
    split_expectd = [0, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000]

    assert split == split_expectd


def test_split_interval_ncuts():
    """
    This is a test function that checks if the split_interval function correctly splits an interval into
    a specified number of cuts.
    """

    start = 0
    end = 1000
    step = None
    ncuts = 4
    split = split_interval(start, end, step=step, ncuts=ncuts)
    split_expectd = [0, 250, 500, 750, 1000]

    assert split == split_expectd


def test_merged_regions():
    """
    The function tests whether the merge_regions function correctly merges a list of genomic regions.
    """

    # Define a list of genomic regions
    regions = [
        ("chr1", 100, 200),
        ("chr1", 150, 300),
        ("chr2", 500, 600),
        ("chr2", 550, 650),
        ("chr3", 800, 900),
    ]
    # Call the merge_regions function
    merged_regions = merge_regions(regions)
    # Define merged regions expected
    merged_regions_expected = [
        ("chr1", 100, 300),
        ("chr2", 500, 650),
        ("chr3", 800, 900),
    ]

    assert merged_regions == merged_regions_expected


def test_create_where_clause():
    """
    This function tests the create_where_clause function by defining a list of merged regions and a
    table, calling the create_where_clause function, and comparing the output to an expected where
    clause.
    """

    # Define a list of merged regions
    merged_regions = [
        ("chr1", 100, 300),
        ("chr2", 500, 650),
        ("chr3", 800, 900),
        ("chr3", 1000, 1200),
    ]

    # define table
    table = "variants"

    # Call the create_where_clause function
    where_clause = create_where_clause(merged_regions, table=table)

    # Create expected where clause
    where_clause_expected = """ ( variants."#CHROM" = 'chr1' AND (   (variants.POS >= 100 AND variants.POS <= 300)  ) )   OR  ( variants."#CHROM" = 'chr2' AND (   (variants.POS >= 500 AND variants.POS <= 650)  ) )   OR  ( variants."#CHROM" = 'chr3' AND (   (variants.POS >= 800 AND variants.POS <= 900)   OR  (variants.POS >= 1000 AND variants.POS <= 1200)  ) ) """

    assert where_clause.strip() == where_clause_expected.strip()


def test_command():
    """
    This function tests the execution of a command and compares its output to an expected result.
    """

    # Creation of a command
    cmd = "echo 'Command'"

    # Execute command
    command_output = command(cmd)

    # Expected result
    results_expected = "Command"

    assert command_output == results_expected


def test_run_parallel_commands():
    """
    This function tests the execution of multiple commands in parallel and compares the results to the
    expected output.
    """

    # Creation of a list of command
    commands = ["echo 'Command 1'", "echo 'Command 2'", "echo 'Command 3'"]

    # Execute commands in parallele
    results = run_parallel_commands(commands, 2)

    # Expected results
    results_expected = ["Command 1", "Command 2", "Command 3"]

    assert results == results_expected


def test_run_parallel_functions():
    """
    This function tests the functionality of running multiple functions in parallel using a specified
    number of threads.
    """

    # List of functions (examples)
    functions = [example_function(1, "hello"), example_function(2, "world")]
    threads = 2

    # Launch functions
    results = run_parallel_functions(functions, threads)

    # Number of result expected
    expected_output_length = 2

    assert len(results) == expected_output_length


def test_example_function():
    """
    This is a test function that checks if the output of the example_function matches the expected
    result.
    """

    # Launch functions
    result = example_function(1, "hello")

    # result expected
    expected_result = [1, "hello"]

    assert result == expected_result


def test_get_index():
    """
    The function `test_get_index()` tests the `get_index()` function with various inputs.
    """

    # Test with a list of values
    values = ["a", "b", "c", "d"]
    assert get_index("a", values) == 0
    assert get_index("d", values) == 3
    assert get_index("e", values) == -1

    # Test with an empty list
    assert get_index("a", []) == -1
    assert get_index("", []) == -1

    # Test with a None element
    assert get_index(None, values) == -1


def test_find_nomen_full():
    """
    The function `test_find_nomen_full()` tests the `find_nomen()` function with various input cases and
    expected outputs.
    """

    # Test case 1
    hgvs = "NM_001637.3:c.1582G>T"
    transcripts = ["NM_001637.3"]
    expected_output = {
        "NOMEN": "NM_001637:c.1582G>T",
        "CNOMEN": "c.1582G>T",
        "RNOMEN": None,
        "NNOMEN": None,
        "PNOMEN": None,
        "TVNOMEN": "NM_001637.3",
        "TNOMEN": "NM_001637",
        "TPNOMEN": None,
        "TPVNOMEN": None,
        "VNOMEN": "3",
        "ENOMEN": None,
        "GNOMEN": None,
    }
    assert find_nomen(hgvs, transcripts=transcripts) == expected_output

    # Test case 2
    hgvs = "NM_001637.3:c.1582G>T,NM_001637.3:c.1583G>T"
    transcripts = ["NM_001637.3"]
    expected_output = {
        "NOMEN": "NM_001637:c.1582G>T",
        "CNOMEN": "c.1582G>T",
        "RNOMEN": None,
        "NNOMEN": None,
        "PNOMEN": None,
        "TVNOMEN": "NM_001637.3",
        "TNOMEN": "NM_001637",
        "TPNOMEN": None,
        "TPVNOMEN": None,
        "VNOMEN": "3",
        "ENOMEN": None,
        "GNOMEN": None,
    }
    assert find_nomen(hgvs, transcripts=transcripts) == expected_output

    # Test case 3
    hgvs = "NM_001637.3:c.1582G>T,NM_001637.3:c.1583G>T,NM_001637.2:c.1582G>T:p.G12D"
    transcripts = ["NM_001637.2", "NM_001637.3"]
    expected_output = {
        "NOMEN": "NM_001637:c.1582G>T:p.G12D",
        "CNOMEN": "c.1582G>T",
        "RNOMEN": None,
        "NNOMEN": None,
        "PNOMEN": "p.G12D",
        "TVNOMEN": "NM_001637.2",
        "TNOMEN": "NM_001637",
        "TPNOMEN": None,
        "TPVNOMEN": None,
        "VNOMEN": "2",
        "ENOMEN": None,
        "GNOMEN": None,
    }
    assert find_nomen(hgvs, transcripts=transcripts) == expected_output

    # Test case 4
    hgvs = "Gene1:exon12:n.1582G>T:NR_001637.3"
    transcripts = []
    expected_output = {
        "NOMEN": "Gene1:NR_001637:exon12:n.1582G>T",
        "CNOMEN": None,
        "RNOMEN": None,
        "NNOMEN": "n.1582G>T",
        "PNOMEN": None,
        "TVNOMEN": "NR_001637.3",
        "TNOMEN": "NR_001637",
        "TPNOMEN": None,
        "TPVNOMEN": None,
        "VNOMEN": "3",
        "ENOMEN": "exon12",
        "GNOMEN": "Gene1",
    }
    assert find_nomen(hgvs, transcripts=transcripts) == expected_output

    # Test case 5
    hgvs = "Gene1:exon12:r.1582G>T"
    transcripts = []
    expected_output = {
        "NOMEN": "Gene1:exon12:r.1582G>T",
        "CNOMEN": None,
        "RNOMEN": "r.1582G>T",
        "NNOMEN": None,
        "PNOMEN": None,
        "TVNOMEN": None,
        "TNOMEN": None,
        "TPNOMEN": None,
        "TPVNOMEN": None,
        "VNOMEN": None,
        "ENOMEN": "exon12",
        "GNOMEN": "Gene1",
    }
    assert find_nomen(hgvs, transcripts=transcripts) == expected_output


def test_find():
    """
    This is a test function for the "find" function, which tests if the function can correctly locate
    files in a given directory and its subdirectories.
    """
    # Crée un fichier pour tester la fonction
    test_path = "test_dir"
    os.makedirs(test_path, exist_ok=True)
    test_file_path = os.path.join(test_path, "test_file.txt")
    with open(test_file_path, "w") as f:
        f.write("Test file content")

    # Teste la fonction avec un nom de fichier qui n'existe pas
    assert find("nonexistent_file.txt", test_path) == ""

    # Teste la fonction avec un nom de fichier qui existe dans le dossier de base
    assert find("test_file.txt", test_path) == test_file_path

    # Teste la fonction avec un nom de fichier qui existe dans un sous-dossier
    sub_dir_path = os.path.join(test_path, "sub_dir")
    os.makedirs(sub_dir_path, exist_ok=True)
    sub_dir_file_path = os.path.join(sub_dir_path, "sub_dir_file.txt")
    with open(sub_dir_file_path, "w") as f:
        f.write("Sub-directory file content")
    assert find("sub_dir_file.txt", test_path) == sub_dir_file_path

    # Nettoie le dossier de test
    os.remove(test_file_path)
    os.remove(sub_dir_file_path)
    os.rmdir(sub_dir_path)
    os.rmdir(test_path)


def test_find_all():
    """
    This function tests the find_all function by creating temporary directories and files and checking
    if the function returns the correct paths.
    """
    # Create a temporary directory structure
    with tempfile.TemporaryDirectory() as tmpdir:
        # Create some files with the name 'test_file' in different directories
        open(os.path.join(tmpdir, "test_file"), "a").close()
        os.makedirs(os.path.join(tmpdir, "subdir"))
        open(os.path.join(tmpdir, "subdir", "test_file"), "a").close()

        # Test that find_all returns the correct paths
        assert find_all("test_file", tmpdir) == [
            os.path.join(tmpdir, "test_file"),
            os.path.join(tmpdir, "subdir", "test_file"),
        ]


def test_find_genome():
    """
    This is a test function for the find_genome function, which checks if the function can correctly
    find an existing genome file and handle errors when a genome file does not exist.
    """

    # Default genome filenaùe
    genome_filename = "hg19.fa"

    # Create existing genome
    tmp_genome = NamedTemporaryFile(prefix="genome.", suffix=".fa", dir="/tmp")
    tmp_genome_name = tmp_genome.name
    # call the function to find the genome file
    genome_path_found = find_genome(genome_path=tmp_genome_name, file=genome_filename)
    # check if the genome file was found
    assert os.path.exists(genome_path_found) and genome_path_found == tmp_genome_name

    # Either genome in the system or not

    # create a temporary directory
    with tempfile.TemporaryDirectory() as tmpdir:
        # specify a non-existent path for the genome file
        genome_path_nonexistent = os.path.join(tmpdir, "nonexistent_genome.fa")
        # call the function to find the genome file
        error = None
        try:
            genome_path_found = find_genome(
                genome_path=genome_path_nonexistent, file=genome_filename
            )
        except:
            with pytest.raises(ValueError) as e:
                genome_path_found = find_genome(
                    genome_path=genome_path_nonexistent, file=genome_filename
                )
                error = str(e.value)

        # check if the genome file was found
        assert (
            os.path.exists(genome_path_found)
            and genome_path_found != genome_path_nonexistent
        ) or error == None


def test_findbypipeline():
    """
    The function `test_findbypipeline()` tests the `findbypipeline()` function with three different test
    cases.
    """

    # Test case 1: Sample with GT
    data = {
        "FORMAT": "GT:DP:AD",
        "S1": "0/1:20:10,10",
        "S2": "0/0:20:20,0",
        "S3": "./.:20:0,0",
    }
    row = pd.Series(data)
    expected_result = "1/3"
    assert findbypipeline(row, ["S1", "S2", "S3"]) == expected_result

    # Test case 2: No sample/pipeline
    data = {
        "FORMAT": "GT:DP:AD",
        "S1": "0/1:20:10,10",
        "S2": "0/0:20:20,0",
        "S3": "./.:20:0,0",
    }
    row = pd.Series(data)
    expected_result = "0/0"
    assert findbypipeline(row, []) == expected_result

    # Test case 3: All samples have missing genotype
    data = {
        "FORMAT": "GT:DP:AD",
        "S1": "./.:20:0,0",
        "S2": "./.:20:0,0",
        "S3": "./.:20:0,0",
    }
    row = pd.Series(data)
    expected_result = "0/3"
    assert findbypipeline(row, ["S1", "S2", "S3"]) == expected_result


def test_genotypeconcordance():
    """
    The function tests the concordance of genotypes between a given variant and a list of samples.
    """

    # Define test data
    test_dataframe = pd.DataFrame(
        {
            "CHROM": ["chr1"],
            "POS": [100],
            "ID": ["rs123"],
            "REF": ["A"],
            "ALT": ["T"],
            "QUAL": [30],
            "FILTER": ["PASS"],
            "INFO": ["AF=0.5"],
            "FORMAT": ["GT:DP"],
            "Sample1": ["0/1:10"],
            "Sample2": ["0/0:20"],
            "Sample3": ["1/1:30"],
        }
    )

    # Test case 1: No samples
    assert genotypeconcordance(test_dataframe.iloc[0], []) == "0/0"

    # Test case 2: All samples have the same genotype
    assert genotypeconcordance(test_dataframe.iloc[0], ["Sample2", "Sample2"]) == "TRUE"

    # Test case 3: At least one sample has a different genotype
    assert (
        genotypeconcordance(test_dataframe.iloc[0], ["Sample1", "Sample3"]) == "FALSE"
    )

    # Test case 4: Some samples have null or unknown genotypes
    assert (
        genotypeconcordance(
            test_dataframe.iloc[0], ["Sample1", "Sample2", "Sample3", "Sample4"]
        )
        == "FALSE"
    )

    # Test case 5: All samples have null or unknown genotypes
    assert (
        genotypeconcordance(test_dataframe.iloc[0], ["Sample4", "Sample5"]) == "FALSE"
    )

    # Test case 6: All samples have the same null or unknown genotype
    assert (
        genotypeconcordance(test_dataframe.iloc[0], ["Sample2", "Sample4", "Sample5"])
        == "TRUE"
    )


def test_genotype_compression():
    """
    This is a test function that checks if the genotype compression function returns the expected output
    for various input test cases.
    """
    # Test cases with expected output
    test_cases = {
        "0/0": "0",
        "0/1": "01",
        "1/1": "1",
        "1|2": "12",
        "./.": "0",
        ".|1": "01",
        "1|1|1": "1",
        ".|.|.": "0",
        "1/2|2": "12",
        ".": "0",
        "": "",
    }

    for genotype, expected_output in test_cases.items():
        assert genotype_compression(genotype) == expected_output


def test_genotype_barcode():
    """
    This is a unit test for the function `genotype_barcode()` which tests various input genotype cases
    and their expected output.
    """
    test_cases = {
        "1|1": "2",
        "0|0": "0",
        "1|0": "1",
        "1/1": "2",
        "0/0": "0",
        "1/0": "1",
        "1/2": "1",
        "0/2": "1",
        "./.": "0",
        ".|.": "0",
        "": "?",
    }
    for input_genotype, expected_output in test_cases.items():
        assert genotype_barcode(input_genotype) == expected_output


def test_barcode():
    """
    The function tests the barcode function with different test cases using a test dataset.
    """

    test_data = pd.DataFrame(
        {
            "CHROM": ["chr1", "chr2"],
            "POS": [100, 200],
            "REF": ["A", "A"],
            "ALT": ["T", "T"],
            "FORMAT": ["GT", "GT"],
            "sample1": ["0/0", "0/0"],
            "sample2": ["0/1", "1|0"],
            "sample3": ["./.", "./0"],
            "sample4": ["1/1", "2/2"],
        }
    )

    # Case 1: empty samples list
    assert barcode(test_data.iloc[0], []) == ""

    # Case 2: all samples have same barcode
    assert barcode(test_data.iloc[0], ["sample1", "sample2", "sample4"]) == "012"

    # Case 3: some samples have same barcode, some have different barcode
    assert (
        barcode(test_data.iloc[1], ["sample1", "sample2", "sample3", "sample4"])
        == "0102"
    )


def test_vaf_normalization():
    """
    This is a unit test function for testing the vaf_normalization function on a sample dataset.
    """

    sample_data = pd.DataFrame(
        {
            "FORMAT": ["GT:AD:DP:GQ:PL", "GT:FREQ", "GT:DP4"],
            "Sample1": ["0/1:10,5:15:99:255,0,255", "0/1:50.0%", "0/1:6,4,3,2"],
            "Sample2": ["1/1:.:.", "0/0", "0/1:4,2,2,1"],
            "Sample3": ["./.:.:.", "./.", ""],
        }
    )

    expected_output = pd.DataFrame(
        {
            "FORMAT": ["GT:AD:DP:GQ:PL", "GT:FREQ", "GT:DP4"],
            "Sample1": [
                "0/1:10,5:15:99:255,0,255:0.333333",
                "0/1:50.0%:0.5",
                "0/1:6,4,3,2:0.333333",
            ],
            "Sample2": ["1/1:.:.:.:.:.", "0/0:.:.", "0/1:4,2,2,1:0.333333"],
            "Sample3": ["./.:.:.:.:.:.", "./.:.:.", "./.:.:."],
        }
    )

    actual_output = sample_data.copy()
    for sample in sample_data.columns[1:]:
        actual_output[sample] = sample_data.apply(
            lambda x: vaf_normalization(x, sample), axis=1
        )

    assert_frame_equal(actual_output, expected_output)


def test_genotype_stats():
    """
    The function tests the genotype statistics of a given dataset.
    """

    test_data = pd.DataFrame(
        {
            "FORMAT": "GT:AD:DP:GQ:PL:VAF",
            "Sample1": "0/1:0,10:10:20:255,0,255:0.5",
            "Sample2": "0/0:5,0:5:15:255,0,255:0",
            "Sample3": "1/1:0,5:5:10:0,255,255:1",
        },
        index=[0],
    )

    expected_output = {
        "VAF_stats_nb": 2,
        "VAF_stats_list": "0.5:1.0",
        "VAF_stats_min": 0.5,
        "VAF_stats_max": 1.0,
        "VAF_stats_mean": 0.75,
        "VAF_stats_mediane": 0.75,
        "VAF_stats_stdev": 0.3535533905932738,
    }

    # test for all samples
    output = genotype_stats(test_data.iloc[0], ["Sample1", "Sample2", "Sample3"], "VAF")
    assert output == expected_output

    # test for one sample only
    output = genotype_stats(test_data.iloc[0], ["Sample1"], "VAF")
    expected_output = {
        "VAF_stats_nb": 1,
        "VAF_stats_list": "0.5",
        "VAF_stats_min": 0.5,
        "VAF_stats_max": 0.5,
        "VAF_stats_mean": 0.5,
        "VAF_stats_mediane": 0.5,
        "VAF_stats_stdev": None,
    }
    assert output == expected_output

    # test for empty samples
    output = genotype_stats(test_data.iloc[0], [], "VAF")
    expected_output = {
        "VAF_stats_nb": 0,
        "VAF_stats_list": None,
        "VAF_stats_min": None,
        "VAF_stats_max": None,
        "VAF_stats_mean": None,
        "VAF_stats_mediane": None,
        "VAF_stats_stdev": None,
    }
    assert output == expected_output

    output = genotype_stats(test_data.iloc[0], ["Sample1"], "invalid")
    expected_output = {
        "invalid_stats_nb": 0,
        "invalid_stats_list": None,
        "invalid_stats_min": None,
        "invalid_stats_max": None,
        "invalid_stats_mean": None,
        "invalid_stats_mediane": None,
        "invalid_stats_stdev": None,
    }
    assert output == expected_output


def test_trio():
    """
    The function "test_trio" is not defined and therefore cannot be summarized.
    """
    # Define a sample DataFrame
    df = pd.DataFrame(
        {
            "CHROM": ["1", "1"],
            "POS": [100, 200],
            "ID": [".", "."],
            "REF": ["A", "C"],
            "ALT": ["T", "G"],
            "QUAL": [10.0, 20.0],
            "FILTER": ["PASS", "PASS"],
            "INFO": ["AF=0.5", "AF=0.2"],
            "FORMAT": "GT:DP:AD:GQ:PL",
            "sample1": ["0/1:10:5,5:99:10,20,30", "0/0:15:15,0:99:0,30,40"],
            "sample2": ["1/1:20:0,20:99:50,0,60", "0/0:25:10,15:99:20,0,80"],
            "sample3": ["1/1:30:0,30:99:70,0,80", "0/1:10:10,0:99:0,10,20"],
        }
    )

    # Test case 1: Empty samples list
    assert trio(df.iloc[0], []) == ""

    # Test case 2: Non-empty samples list with dominant variant
    assert trio(df.iloc[0], ["sample1", "sample2", "sample3"]) == "recessive"

    # Test case 3: Non-empty samples list with recessive variant
    assert trio(df.iloc[0], ["sample1", "sample1", "sample1"]) == "dominant"

    # Test case 4: Non-empty samples list with recessive variant
    assert trio(df.iloc[1], ["sample1", "sample2", "sample3"]) == "denovo"

    # Test case 5: Non-existent sample name
    assert trio(df.iloc[0], ["sample1", "non_existent_sample", "sample2"]) == "unknown"


def test_load_duckdb_extension():
    """
    The function tests whether the `load_duckdb_extension` function can successfully load a specified
    extension in a DuckDB connection.
    """

    # Create connexion
    conn = duckdb.connect()

    # tests
    assert load_duckdb_extension(conn, ["sqlite_scanner"])
    assert load_duckdb_extension(conn, ["json"])
    assert not load_duckdb_extension(conn, ["not_an_extension"])


def test_get_duckdb_extension_file():
    """
    The function "test_get_duckdb_extension_file" tests the "get_duckdb_extension_file" function with
    the argument "sqlite_scanner".
    """

    # Create connexion
    conn = duckdb.connect()

    assert get_duckdb_extension_file("sqlite_scanner", conn=conn)


def test_clean_annotation_field_basic_alphanumeric():
    assert clean_annotation_field("HelloWorld") == "HelloWorld"


def test_clean_annotation_field_with_special_characters():
    assert clean_annotation_field("Hello, World!") == "HelloWorld"


def test_clean_annotation_field_with_allowed_characters():
    assert clean_annotation_field("Hello-World", char_allowed=["-"]) == "Hello-World"


def test_clean_annotation_field_empty_string():
    assert clean_annotation_field("") == ""


def test_clean_annotation_field_no_allowed_characters():
    assert clean_annotation_field("Hello@World#2023") == "HelloWorld2023"


def test_clean_annotation_field_all_characters_removed():
    assert clean_annotation_field("!!!") == ""


def test_clean_annotation_field_non_alphanumeric_with_allowed_chars():
    assert clean_annotation_field("Test123!@#", char_allowed=["!"]) == "Test123!"
