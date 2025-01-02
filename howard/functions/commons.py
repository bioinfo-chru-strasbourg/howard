import io
import multiprocessing
import os
from pathlib import Path
import platform
import re
import statistics
import string
import subprocess
import sys
from tempfile import NamedTemporaryFile
import tempfile
import typing
import duckdb
import json
import argparse
import pandas as pd
import vcf
import logging as log
import shutil
import urllib.request
import zipfile
import gzip
import requests
import fnmatch
import ast

import random

import pgzip
import mgzip
import bgzip

import pysam
import yaml
import pysam.bcftools

import signal
from contextlib import contextmanager

from configparser import ConfigParser, NoSectionError, NoOptionError

from importlib.metadata import metadata

from shutil import which

file_folder = os.path.dirname(__file__)

# plugin subfolder
subfolder_plugins = "plugins"

# Main folder
folder_main = os.path.abspath(os.path.join(file_folder, "../.."))
folder_config = os.path.abspath(os.path.join(folder_main, "config"))
folder_user_home = os.path.abspath(os.path.expanduser("~"))
folder_howard_home = os.path.join(folder_user_home, "howard")
folder_plugins = os.path.join(folder_main, "plugins")

comparison_map = {
    "gt": ">",
    "gte": ">=",
    "lt": "<",
    "lte": "<=",
    "equals": "=",
    "contains": "SIMILAR TO",
}


code_type_map = {"Integer": 0, "String": 1, "Float": 2, "Flag": 3}


code_type_map_to_sql = {
    "Integer": "INTEGER",
    "String": "VARCHAR",
    "Float": "FLOAT",
    "Flag": "VARCHAR",
}


file_format_delimiters = {"vcf": "\t", "tsv": "\t", "csv": ",", "psv": "|", "bed": "\t"}

file_format_allowed = list(file_format_delimiters.keys()) + [
    "json",
    "parquet",
    "duckdb",
]

file_compressed_format = ["gz", "bgz"]


vcf_required_release = "##fileformat=VCFv4.2"
vcf_required_columns = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"]

vcf_required = [vcf_required_release, "\t".join(vcf_required_columns)]

# Tools
DEFAULT_TOOLS_FOLDER = os.path.join(folder_howard_home, "tools")

DEFAULT_TOOLS_BIN = {
    "bcftools": {"bin": "bcftools"},
    "bgzip": {"bin": "bgzip"},
    "java": {"bin": "java"},
    "snpeff": {"jar": "~/howard/tools/snpeff/current/bin/snpEff.jar"},
    "annovar": {"perl": "~/howard/tools/annovar/current/bin/table_annovar.pl"},
    "exomiser": {"jar": "~/howard/tools/exomiser/current/bin/exomiser.jar"},
    "docker": {"bin": "docker"},
    "splice": {
        "docker": {
            "image": "bioinfochrustrasbourg/splice:0.2.1",
            "entrypoint": "/bin/bash",
            "options": None,
            "command": None,
        }
    },
}

# URL
DEFAULT_ANNOVAR_URL = "http://www.openbioinformatics.org/annovar/download"
DEFAULT_REFSEQ_URL = "http://hgdownload.soe.ucsc.edu/goldenPath"
DEFAULT_DBNSFP_URL = "https://dbnsfp.s3.amazonaws.com"
DEFAULT_EXOMISER_URL = "http://data.monarchinitiative.org/exomiser"
DEFAULT_EXOMISER_REMM_URL = "https://kircherlab.bihealth.org/download/ReMM"
DEFAULT_EXOMISER_CADD_URL = "https://kircherlab.bihealth.org/download/CADD"
DEFAULT_ALPHAMISSENSE_URL = "https://storage.googleapis.com/dm_alphamissense"
DEFAULT_DBSNP_URL = "https://ftp.ncbi.nih.gov/snp/archive"


# Databases default folder
DEFAULT_DATABASE_FOLDER = os.path.join(folder_howard_home, "databases")
DEFAULT_ANNOTATIONS_FOLDER = f"{DEFAULT_DATABASE_FOLDER}/annotations/current"
DEFAULT_GENOME_FOLDER = f"{DEFAULT_DATABASE_FOLDER}/genomes/current"
DEFAULT_SNPEFF_FOLDER = f"{DEFAULT_DATABASE_FOLDER}/snpeff/current"
DEFAULT_ANNOVAR_FOLDER = f"{DEFAULT_DATABASE_FOLDER}/annovar/current"
DEFAULT_REFSEQ_FOLDER = f"{DEFAULT_DATABASE_FOLDER}/refseq/current"
DEFAULT_DBNSFP_FOLDER = f"{DEFAULT_DATABASE_FOLDER}/dbnsfp/current"
DEFAULT_EXOMISER_FOLDER = f"{DEFAULT_DATABASE_FOLDER}/exomiser/current"
DEFAULT_DBSNP_FOLDER = f"{DEFAULT_DATABASE_FOLDER}/exomiser/dbsnp"
DEFAULT_SPLICE_FOLDER = f"{DEFAULT_DATABASE_FOLDER}/splice"
DEFAULT_SPLICEAI_FOLDER = f"{DEFAULT_DATABASE_FOLDER}/spliceai"
DEFAULT_SPIP_FOLDER = f"{DEFAULT_DATABASE_FOLDER}/spip"

# Data default folder
DEFAULT_DATA_FOLDER = os.path.join(folder_howard_home, "data")

# Deefault Assembly
DEFAULT_ASSEMBLY = "hg19"

# DuckDB extension
DUCKDB_EXTENSION = f"{file_folder}/duckdb_extension"

# Variables
MACHIN_LIST = {"amd64": "amd64", "arm64": "arm64"}

# bcftools format allowed
BCFTOOLS_FORMAT = ["vcf", "bed"]

LOG_FORMAT = "#[%(asctime)s] [%(levelname)s] %(message)s"

CODE_TYPE_MAP = {"Integer": 0, "String": 1, "Float": 2, "Flag": 3}

GENOTYPE_MAP = {None: ".", -1: "A", -2: "G", -3: "R"}

DTYPE_LIMIT_AUTO = 10000

DEFAULT_CHUNK_SIZE = 1024 * 1024


def remove_if_exists(filepaths: list) -> None:
    """
    The function removes a file if it exists at the specified filepath(s).

    :param filepaths: A list of file paths that you want to check for existence and remove if they exist
    :type filepaths: list
    """

    # If input is string
    if isinstance(filepaths, str):
        filepaths = [filepaths]

    for filepath in filepaths:

        if filepath:

            # full path
            filepath = full_path(filepath)

            # If path exists
            if os.path.exists(filepath):

                # Folder
                if os.path.isdir(filepath):
                    shutil.rmtree(filepath)

                # File
                else:
                    os.remove(filepath)


def set_log_level(verbosity: str, log_file: str = None) -> str:
    """
    It sets the log level of the Python logging module

    :param verbosity: The level of verbosity
    """

    verbosity = verbosity.lower()
    configs = {
        "debug": log.DEBUG,
        "info": log.INFO,
        "warning": log.WARNING,
        "error": log.ERROR,
        "critical": log.CRITICAL,
        "notset": log.NOTSET,
    }
    if verbosity not in configs.keys():
        raise ValueError("Unknown verbosity level:" + verbosity)

    log.basicConfig(
        filename=log_file,
        encoding="utf-8",
        format=LOG_FORMAT,
        datefmt="%Y-%m-%d %H:%M:%S",
        level=configs[verbosity],
    )

    return verbosity


def split_interval(start: int, end: int, step: int = None, ncuts: int = None):
    """
    It takes a start and end value, and either a step size or a number of cuts, and returns a list of
    values that split the interval into equal-sized pieces

    :param start: the start of the interval
    :param end: the end of the interval
    :param step: the step size between each cut
    :param ncuts: number of cuts to make
    :return: A list of numbers.
    """
    if step is None and ncuts is None:
        raise ValueError("Either step or ncuts must be provided")
    if step is not None and ncuts is not None:
        raise ValueError("Only one of step or ncuts must be provided")
    if step is not None:
        return list(range(start, end, step)) + [end]
    if ncuts is not None:
        step = (end - start) / ncuts
        return [start + i * step for i in range(ncuts + 1)]


def merge_regions(regions: list) -> list:
    """
    It takes a list of genomic regions and returns a list of genomic regions where overlapping regions
    have been merged

    :param regions: A list of tuples representing genomic regions with the values of the chrom, start
    and end columns
    :return: A list of tuples representing the merged regions with the values of the columns chrom,
    start and end.
    """

    merged_regions = []

    if regions:

        # Sort regions by chromosomes and first position
        sorted_regions = sorted(regions, key=lambda x: (x[0], x[1]))

        # Init current region
        current_region = sorted_regions[0]

        # Fetch sorted regions
        for region in sorted_regions[1:]:
            # If current region overlap next region, merge both
            if current_region[0] == region[0] and current_region[2] >= region[1]:
                current_region = (
                    current_region[0],
                    current_region[1],
                    max(current_region[2], region[2]),
                )
            # Else, add current region to merged regions list, and next region
            else:
                merged_regions.append(current_region)
                current_region = region

        # Add last region to merged regions list
        merged_regions.append(current_region)

    return merged_regions


def create_where_clause(merged_regions: list, table: str = "variants") -> str:
    """
    It takes a list of merged regions and returns a SQL WHERE clause that can be used to filter variants
    in a SQL table

    :param merged_regions: a list of tuples representing the merged regions with the values of the
    chrom, start and end columns
    :param table: The name of the table to query, defaults to variants (optional)
    :return: A dictionary with the chromosome as key and the where clause as value.
    """

    where_clause = " "
    where_clause_chrom = {}
    for i, region in enumerate(merged_regions):
        chrom = region[0]
        start = region[1]
        stop = region[2]
        if chrom not in where_clause_chrom:
            where_clause_chrom[chrom] = ""
            where_clause_chrom_sep = ""
        else:
            where_clause_chrom_sep = " OR "
        where_clause_chrom[
            chrom
        ] += f" {where_clause_chrom_sep} ({table}.POS >= {start} AND {table}.POS <= {stop}) "

    nb_chrom = 0
    where_clause_sep = ""
    for chrom in where_clause_chrom:
        nb_chrom += 1
        if nb_chrom > 1:
            where_clause_sep = " OR "
        where_clause += f" {where_clause_sep} ( {table}.\"#CHROM\" = '{chrom}' AND ( {where_clause_chrom[chrom]} ) ) "

    return where_clause


def command(command: str) -> str:
    """
    It runs a command in the shell and waits for it to finish

    :param command: The command to run
    :return: The return value is the exit status of the process.
    """

    output = subprocess.check_output(command, shell=True)
    return output.decode("utf-8").strip()


def run_parallel_commands(commands: list, threads: int = 1) -> list:
    """
    It takes a list of commands and a number of threads, and runs the commands in parallel

    :param commands: a list of commands to run
    :param threads: The number of threads to use
    :return: A list of results from the commands.
    """
    pool = multiprocessing.Pool(threads)
    results = []
    for cmd in commands:
        results.append(
            pool.apply_async(command, args=(cmd,), error_callback=lambda e: print(e))
        )
    pool.close()
    pool.join()
    return [result.get().strip() for result in results]


def run_parallel_functions(functions: list, threads: int = 1) -> list:
    """
    It takes a list of functions and a number of threads, and runs the functions in parallel using the
    number of threads specified

    :param functions: a list of functions to run in parallel
    :param threads: The number of threads to use
    :return: A list of multiprocessing.pool.ApplyResult objects.
    """
    pool = multiprocessing.Pool(threads)
    results = []
    for func in functions:
        results.append(
            pool.apply_async(func, args=(1, "hello"), error_callback=lambda e: print(e))
        )
    pool.close()
    pool.join()
    return results


def example_function(num, word):
    """
    `example_function` takes in a number and a word and returns a list of the number and the word

    :param num: a number
    :param word: a string
    :return: [num, word]
    """
    return [num, word]


def find(name: str, path: str) -> str:
    """
    It recursively walks the directory tree starting at the given path, and returns the first file it
    finds with the given name

    :param name: The name of the file you're looking for
    :param path: The path to search for the file
    :return: The path to the file.
    """
    for root, dirs, files in os.walk(path):
        if name in files:
            return os.path.join(root, name)
    return ""


def find_all(name: str, path: str) -> list:
    """
    "Walk the directory tree starting at path, and for each regular file with the name name, append its
    full path to the result list."

    The os.walk function is a generator that yields a 3-tuple containing the name of a directory, a list
    of its subdirectories, and a list of the files in that directory. The name of the directory is a
    string, and the lists of subdirectories and files are lists of strings

    :param name: The name of the file you're looking for
    :param path: The path to search in
    :return: A list of all the files in the directory that have the name "name"
    """

    # result
    result = []

    # Full path
    path = full_path(path)

    for root, dirs, files in os.walk(path):
        for filename in fnmatch.filter(files, name):
            if os.path.exists(os.path.join(root, filename)):
                result.append(os.path.join(root, filename))
    return result


def find_genome(genome_path: str, assembly: str = None, file: str = None) -> str:
    """
    The `find_genome` function checks if a genome file exists at the specified path, and if not, it
    tries to find it using the provided assembly name or file name.

    :param genome_path: The path to the genome file
    :type genome_path: str
    :param assembly: The `assembly` parameter is a string that represents the name of the genome
    assembly. It is used to search for the genome file with the specified assembly name in the
    `genome_dir` directory. If a genome file with the assembly name is found, its path is returned
    :type assembly: str
    :param file: The `file` parameter is the name of the genome file that you want to find
    :type file: str
    :return: the path to the genome file.
    """

    # Full path
    genome_path = full_path(genome_path)

    # check genome
    if os.path.exists(genome_path) and not os.path.isdir(genome_path):
        return genome_path
    else:
        log.debug(f"Genome warning: Try to find genome in '{genome_path}'...")
        genome_dir = genome_path
        genome_path = ""
        # Try to find genome
        if file and find_all(file, genome_dir):
            genome_path = find_all(file, genome_dir)[0]
        elif assembly and find_all(assembly + ".fa", genome_dir):
            genome_path = find_all(assembly + ".fa", genome_dir)[0]
    return genome_path


def find_file_prefix(
    input_file: str = None, prefix: str = None, folder: str = None, assembly: str = None
) -> str:
    """
    The function `find_file_prefix` is used to find a specific file based on input parameters such as
    input file, folder, and assembly.

    :param input_file: The input file is the file that you want to find the prefix for. It can be a file
    path or just the file name if it is in the current directory
    :type input_file: str
    :param folder: The `folder` parameter is a string that represents the directory where the file is
    located
    :type folder: str
    :param assembly: The "assembly" parameter is a string that represents the assembly version of the
    file you are looking for. It is used to search for files with the specific assembly version in their
    filename
    :type assembly: str
    :return: the path of the output file.
    """

    # Output
    output_file = None

    # Full path
    input_file = full_path(input_file)
    folder = full_path(folder)

    if input_file and os.path.exists(input_file):
        output_file = input_file
    else:
        # Find in specific assembly folder
        if find_all(f"{prefix}.txt", f"{folder}/{assembly}"):
            output_file = find_all(f"{prefix}.txt", f"{folder}/{assembly}")[0]
        # Find with assembly in filename
        elif find_all(f"{prefix}.{assembly}.txt", folder):
            output_file = find_all(f"{prefix}.{assembly}.txt", folder)[0]
        # Find within the entire folder
        elif find_all(f"{prefix}.txt", folder):
            output_file = find_all(f"{prefix}.txt", folder)[0]

    return output_file


def find_nomen(
    hgvs: pd.DataFrame,
    transcript: str = None,
    transcripts: list = [],
    transcripts_source_order: list = None,
    pattern=None,
    transcripts_len: int = None,
) -> dict:
    """
    The function `find_nomen` takes a HGVS string and a list of transcripts, parses the HGVS string, and
    returns a dictionary with the best NOMEN based on specified patterns.

    :param hgvs: The `hgvs` parameter is a DataFrame containing the HGVS strings to parse. It seems like
    the function is designed to process multiple HGVS strings at once. You can pass this DataFrame to
    the function for processing. If you have a specific DataFrame that you would like to use, please
    provide it
    :type hgvs: pd.DataFrame
    :param transcript: The `transcript` parameter in the `find_nomen` function is used to specify a
    single transcript to use for ranking. It is a string that represents the transcript. If provided,
    this transcript will be used along with the transcripts from the `transcripts` list to determine the
    best NOMEN
    :type transcript: str
    :param transcripts: Transcripts are a list of transcripts to use for ranking in the `find_nomen`
    function. You can provide a list of transcripts that you want to consider when constructing the
    NOMEN for a given HGVS string
    :type transcripts: list
    :param transcripts_source_order: The `transcripts_source_order` parameter is a list that specifies
    the order in which different sources of transcripts should be considered. In the provided function,
    the default order is `["column", "file"]`, which means that transcripts from a column in the input
    data will be considered first, followed by
    :type transcripts_source_order: list
    :param pattern: The `pattern` parameter in the `find_nomen` function is used to specify the format
    in which the NOMEN should be constructed. By default, the pattern is set to
    "GNOMEN:TNOMEN:ENOMEN:CNOMEN:RNOMEN:NNOMEN
    :return: The function `find_nomen` returns a dictionary containing the following keys:
    - NOMEN
    - CNOMEN
    - RNOMEN
    - NNOMEN
    - PNOMEN
    - TVNOMEN
    - TNOMEN
    - TPVNOMEN
    - TPNOMEN
    - VNOMEN
    - ENOMEN
    - GNOMEN
    """

    # Result structure
    empty_nomen_dict = {
        "NOMEN": None,
        "CNOMEN": None,
        "RNOMEN": None,
        "NNOMEN": None,
        "PNOMEN": None,
        "TVNOMEN": None,
        "TNOMEN": None,
        "TPVNOMEN": None,
        "TPNOMEN": None,
        "VNOMEN": None,
        "ENOMEN": None,
        "GNOMEN": None,
    }
    nomen_dict = empty_nomen_dict.copy()

    # Transcripts list
    transcripts_list = transcripts.copy()

    # Transcript as a list (create dict with rank, slow)
    if isinstance(transcripts_list, list):
        transcripts_list_tmp = {}
        for rank, source in enumerate(transcripts_list, start=1):
            transcripts_list_tmp[source] = rank
        transcripts_list = transcripts_list_tmp

    # Pattern
    if pattern is None:
        pattern = "GNOMEN:TNOMEN:ENOMEN:CNOMEN:RNOMEN:NNOMEN:PNOMEN"

    # Transcripts source order and rank
    if transcripts_source_order is None:
        transcripts_source_order = ["column", "file"]
    transcripts_source_order_rank = {}
    for rank, source in enumerate(transcripts_source_order, start=1):
        transcripts_source_order_rank[source] = rank

    # Transcripts file length
    if transcripts_len is None:
        transcripts_len = len(transcripts_list)

    # Check if transcript in column transcript, and add to transcripts list/rank
    transcript_shift = 0
    if transcript is not None and str(transcript) != "nan":
        transcript_list = str(transcript).split(",")
        if transcripts_source_order_rank.get(
            "column", 0
        ) < transcripts_source_order_rank.get("file", 0):
            transcript_shift = len(transcript_list)
        for rank, transcript in enumerate(transcript_list, start=1):
            if transcripts_source_order_rank.get(
                "column", 0
            ) < transcripts_source_order_rank.get("file", 0):
                transcripts_list[transcript] = rank - len(transcript_list)
            elif transcript not in transcripts:
                transcripts_list[transcript] = rank + transcripts_len

    # If hgvs not empty
    if str(hgvs) != "nan":

        hgvs_split = str(hgvs).split(",")

        nomen_score_max = 0

        for one_hgvs in hgvs_split:
            one_hgvs_split = one_hgvs.split(":")

            one_nomen_score = 0
            one_nomen_dict = empty_nomen_dict.copy()

            for one_hgvs_infos in one_hgvs_split:

                if re.match(r"^[NX][MR]_(.*)$", one_hgvs_infos):
                    # Transcript with version
                    one_nomen_dict["TVNOMEN"] = one_hgvs_infos
                    one_nomen_score += 1
                    # Split transcript
                    one_hgvs_infos_split = one_hgvs_infos.split(".")
                    # Transcript
                    one_nomen_dict["TNOMEN"] = one_hgvs_infos_split[0]
                    # Transcript version
                    if len(one_hgvs_infos_split) > 1:
                        one_nomen_dict["VNOMEN"] = one_hgvs_infos_split[1]
                    # NOMEN Score
                    if re.match(r"^NM_(.*)$", one_hgvs_infos) or re.match(
                        r"^NM_(.*)$", one_hgvs_infos
                    ):
                        one_nomen_score += 2
                    elif re.match(r"^NR_(.*)$", one_hgvs_infos):
                        one_nomen_score += 1
                    # NOMEN with default transcript
                    if (
                        one_nomen_dict["TVNOMEN"] in transcripts_list
                        or one_nomen_dict["TNOMEN"] in transcripts_list
                    ):
                        rank = (
                            max(
                                transcripts_list.get(
                                    one_nomen_dict["TVNOMEN"], -transcript_shift - 1
                                ),
                                transcripts_list.get(
                                    one_nomen_dict["TNOMEN"], -transcript_shift - 1
                                ),
                            )
                            + transcript_shift
                        )
                        if rank >= 0:
                            one_nomen_score += 100 * (len(transcripts_list) - rank + 1)

                elif re.match(r"^[NX]P_(.*)$", one_hgvs_infos):
                    # Transcript Protein with version
                    one_nomen_dict["TPVNOMEN"] = one_hgvs_infos
                    one_nomen_score += 1
                    # Split transcript
                    one_hgvs_infos_split = one_hgvs_infos.split(".")
                    # Transcript Protein
                    one_nomen_dict["TPNOMEN"] = one_hgvs_infos_split[0]

                elif (
                    re.match(r"^c\.(.*)$", one_hgvs_infos)
                    or re.match(r"^g\.(.*)$", one_hgvs_infos)
                    or re.match(r"^m\.(.*)$", one_hgvs_infos)
                ):
                    one_nomen_dict["CNOMEN"] = one_hgvs_infos
                    one_nomen_score += 1
                elif re.match(r"^n\.(.*)$", one_hgvs_infos):
                    one_nomen_dict["NNOMEN"] = one_hgvs_infos
                    one_nomen_score += 1
                elif re.match(r"^r\.(.*)$", one_hgvs_infos):
                    one_nomen_dict["RNOMEN"] = one_hgvs_infos
                    one_nomen_score += 1
                elif re.match(r"^p\.(.*)$", one_hgvs_infos):
                    one_nomen_dict["PNOMEN"] = one_hgvs_infos
                    one_nomen_score += 1
                elif re.match(r"^exon(.*)$", one_hgvs_infos):
                    one_nomen_dict["ENOMEN"] = one_hgvs_infos
                    one_nomen_score += 1
                else:
                    one_nomen_dict["GNOMEN"] = one_hgvs_infos

            if one_nomen_score > nomen_score_max:
                nomen_dict = one_nomen_dict.copy()
                nomen_score_max = one_nomen_score

        # Contruct NOMEN from pattern
        nomen = []
        for n in pattern.split(":"):
            if nomen_dict.get(n, None):
                nomen.append(nomen_dict.get(n, None))
        nomen_dict["NOMEN"] = ":".join(nomen)

    return nomen_dict


def explode_annotation_format(
    annotation: str = "",
    uniquify: bool = False,
    output_format: str = "fields",
    prefix: str = "ANN_",
    header: list = [
        "Allele",
        "Annotation",
        "Annotation_Impact",
        "Gene_Name",
        "Gene_ID",
        "Feature_Type",
        "Feature_ID",
        "Transcript_BioType",
        "Rank",
        "HGVS.c",
        "HGVS.p",
        "cDNA.pos / cDNA.length",
        "CDS.pos / CDS.length",
        "AA.pos / AA.length",
        "Distance",
        "ERRORS / WARNINGS / INFO",
    ],
) -> str:
    """
    The `explode_annotation_format` function takes an annotation string and formats it into a specified
    output format with optional customization parameters.

    :param annotation: The `annotation` parameter is a string containing multiple annotations separated
    by commas and pipe symbols. Each annotation consists of different fields separated by pipe symbols.
    For example, an annotation string could look like this: "A|B|C,D|E|F"
    :type annotation: str
    :param uniquify: The `uniquify` parameter in the `explode_annotation_format` function is a boolean
    flag that determines whether to keep only unique values for each annotation field. If set to `True`,
    only unique values will be retained for each field before joining them together. If set to `False`,
    all values, defaults to False
    :type uniquify: bool (optional)
    :param output_format: The `output_format` parameter specifies the format in which you want the
    output to be generated. The function supports two output formats: "fields" and "JSON". If you choose
    "fields", the output will be a string with annotations separated by semicolons. If you choose
    "JSON", the, defaults to fields
    :type output_format: str (optional)
    :param prefix: The `prefix` parameter in the `explode_annotation_format` function is used to specify
    the prefix that will be added to each annotation field when generating the exploded annotation
    string. In the provided function, the default prefix value is set to "ANN_". You can customize this
    prefix value to suit your specific, defaults to ANN_
    :type prefix: str (optional)
    :param header: The `header` parameter in the `explode_annotation_format` function is a list of
    column names that will be used to create a DataFrame from the input annotation string. Each element
    in the `header` list corresponds to a specific field in the annotation data
    :type header: list
    :return: The function `explode_annotation_format` returns a string that contains the exploded and
    formatted annotation information based on the input parameters provided. The format of the returned
    string depends on the `output_format` parameter. If `output_format` is set to "JSON", the function
    returns a JSON-formatted string. Otherwise, it returns a string with annotations formatted based on
    the other parameters such as `uniquify
    """

    # Split annotation ann values
    annotation_infos = [x.split("|") for x in annotation.split(",")]

    # Create dictionary
    annotation_dict = {header[i]: [] for i in range(len(header))}
    for info in annotation_infos:
        for i in range(len(header)):
            annotation_dict[header[i]].append(info[i])

    # Fetch each annotations
    if output_format.upper() == "JSON":

        # Transpose dict
        annotation_explode = {
            i: {f"{prefix}{key}": value[i] for key, value in annotation_dict.items()}
            for i in range(len(next(iter(annotation_dict.values()))))
        }

    else:
        ann_list = []
        for key, values in annotation_dict.items():
            if uniquify:
                values = set(values)
            ann_list_infos = ",".join(values)
            if ann_list_infos:
                header_clean = "".join(char for char in key if char.isalnum())
                ann_list.append(f"{prefix}{header_clean}={ann_list_infos}")

        # join list
        annotation_explode = ";".join(ann_list)

    return annotation_explode


def extract_snpeff_hgvs(
    snpeff: str = "",
    header: list = [
        "Allele",
        "Annotation",
        "Annotation_Impact",
        "Gene_Name",
        "Gene_ID",
        "Feature_Type",
        "Feature_ID",
        "Transcript_BioType",
        "Rank",
        "HGVS.c",
        "HGVS.p",
        "cDNA.pos / cDNA.length",
        "CDS.pos / CDS.length",
        "AA.pos / AA.length",
        "Distance",
        "ERRORS / WARNINGS / INFO",
    ],
) -> str:
    """
    This function extracts HGVS annotations from a given snpEff annotation string and returns them as a
    comma-separated string.

    :param snpeff: The `snpeff` parameter is a string that contains annotations for genetic variants in
    a specific format. It is used as input to extract HGVS notation for the variants
    :type snpeff: str
    :param header: The header parameter is a list of column names that will be used to create a pandas
    DataFrame from the snpeff string input. It is used to extract specific information from the snpeff
    annotations
    :type header: str
    :return: a string that contains the HGVS annotations extracted from the input SNPEff annotation
    string.
    """

    # Split snpeff ann values
    snpeff_infos = [x.split("|") for x in snpeff.split(",")]

    # Create Dataframe
    snpeff_dict = {}
    for i in range(len(header)):
        snpeff_dict[header[i]] = [x[i] for x in snpeff_infos]
    df = pd.DataFrame.from_dict(snpeff_dict, orient="index").transpose()

    # Fetch each annotations
    hgvs_list = []
    for i, row in df.iterrows():

        # Catch values
        gene_id = row["Gene_ID"]
        feature_id = row["Feature_ID"]
        rank = row["Rank"]
        hgvs_c = row["HGVS.c"]
        hgvs_p = row["HGVS.p"]

        # Concatenate with ":" if not empty
        values = []
        if gene_id != "":
            values.append(gene_id)
        if feature_id != "":
            values.append(feature_id)
        if rank != "":
            values.append("exon" + rank.split("/")[0])
        if hgvs_c != "":
            values.append(hgvs_c)
        if hgvs_p != "":
            values.append(hgvs_p)
        hgvs = ":".join(values)

        # Add to list
        hgvs_list.append(hgvs)

    # join list
    snpeff_hgvs = ",".join(hgvs_list)

    return snpeff_hgvs


def explode_snpeff_ann(
    snpeff: str = "",
    uniquify: bool = False,
    output_format: str = "fields",
    prefix: str = "ANN_",
    header: list = [
        "Allele",
        "Annotation",
        "Annotation_Impact",
        "Gene_Name",
        "Gene_ID",
        "Feature_Type",
        "Feature_ID",
        "Transcript_BioType",
        "Rank",
        "HGVS.c",
        "HGVS.p",
        "cDNA.pos / cDNA.length",
        "CDS.pos / CDS.length",
        "AA.pos / AA.length",
        "Distance",
        "ERRORS / WARNINGS / INFO",
    ],
) -> str:
    """
    The `explode_snpeff_ann` function takes a string of SNPEff annotations, splits and processes them
    based on specified parameters, and returns the processed annotations in a specified output format.

    :param snpeff: The `snpeff` parameter is a string containing annotations separated by commas. Each
    annotation is further divided into different fields separated by pipes (|)
    :type snpeff: str
    :param uniquify: The `uniquify` parameter in the `explode_snpeff_ann` function is a boolean flag
    that determines whether to keep only unique values for each annotation field or not. If `uniquify`
    is set to `True`, only unique values will be kept for each annotation field. If, defaults to False
    :type uniquify: bool (optional)
    :param output_format: The `output_format` parameter in the `explode_snpeff_ann` function specifies
    the format in which the output will be generated. The function supports two output formats: "fields"
    and "JSON", defaults to fields
    :type output_format: str (optional)
    :param prefix: The `prefix` parameter in the `explode_snpeff_ann` function is used to specify the
    prefix that will be added to each annotation field in the output. For example, if the prefix is set
    to "ANN_", then the output annotations will be formatted as "ANN_Annotation=example_annotation,
    defaults to ANN_
    :type prefix: str (optional)
    :param header: The `header` parameter in the `explode_snpeff_ann` function is a list of strings that
    represent the column names or fields for the output data. These strings include information such as
    allele, annotation, gene name, gene ID, feature type, transcript biotype, and various other details
    related
    :type header: list
    :return: The function `explode_snpeff_ann` returns a string that contains the exploded and formatted
    SNPEff annotations based on the input parameters provided. The specific format of the returned
    string depends on the `output_format`, `uniquify`, and other parameters specified in the function.
    """

    # Split snpeff ann values
    snpeff_infos = [x.split("|") for x in snpeff.split(",")]

    # Create Dataframe
    snpeff_dict = {}
    for i in range(len(header)):
        if output_format.upper() in ["JSON"]:
            header_clean = header[i]
        else:
            header_clean = "".join(char for char in header[i] if char.isalnum())
        snpeff_dict[header_clean] = [x[i] for x in snpeff_infos]
    df = pd.DataFrame.from_dict(snpeff_dict, orient="index").transpose()

    # Fetch each annotations
    if output_format.upper() in ["JSON"]:
        snpeff_ann_explode = df.transpose().to_json()
    else:
        ann_list = []
        for annotation in df:
            if uniquify:
                ann_list_infos = ",".join(df[annotation].unique())
            else:
                ann_list_infos = ",".join(df[annotation])
            if ann_list_infos:
                ann_list.append(f"{prefix}{annotation}={ann_list_infos}")

        # join list
        snpeff_ann_explode = ";".join(ann_list)

    return snpeff_ann_explode


def get_index(value, values: list = []) -> int:
    """
    The function returns the index of a given value in a list, or -1 if the value is not in the list.

    :param value: The value to search for in the list
    :param values: The parameter "values" is a list of values in which we want to find the index of a
    specific value. It is an optional parameter with a default value of an empty list
    :type values: list
    :return: The function `get_index` returns the index of the first occurrence of the `value` parameter
    in the `values` list. If the `value` parameter is not found in the `values` list, the function
    returns -1.
    """
    try:
        return values.index(value)
    except ValueError:
        # If element does not exist, return -1
        return -1


def get_file_format(filename: str = None) -> str:
    """
    It takes a filename and returns the file format

    :param filename: the name of the file you want to get the format of
    :type filename: str
    :return: The file format of the file.
    """
    if filename:
        filename_name, filename_extension = os.path.splitext(filename)
        filename_format = filename_extension.replace(".", "")
        if filename_format in file_compressed_format:
            _, filename_name_extension = os.path.splitext(filename_name)
            filename_format = filename_name_extension.replace(".", "")
    else:
        filename_format = "unknown"
    return filename_format


def findbypipeline(df, samples: list = []):
    """
    This function takes a dataframe and a list of samples, and returns the number of pipelines found in
    the samples that have a non-null GT value.

    :param df: The input dataframe containing genetic variant information
    :param samples: The `samples` parameter is a list of strings representing the names of the
    samples/pipelines to be searched for in the input dataframe `df`
    :type samples: list
    :return: a string in the format of "nb_pipeline_find/nb_pipeline", where nb_pipeline_find is the
    number of pipelines in the input list samples that have a non-null GT value in the input dataframe
    df, and nb_pipeline is the total number of pipelines in the input list samples. If the input list
    samples is empty, the function returns "0/0".
    """

    # format
    format_fields = df["FORMAT"].split(":")

    # no sample/pipeline
    if not samples:
        return "0/0"

    # init
    nb_pipeline = len(samples)
    nb_pipeline_find = 0

    # For each sample/pipeline
    for sample in samples:

        # Split snpeff ann values
        sample_infos = df[sample].split(":")

        # Create Dataframe
        sample_dict = {}
        for i in range(len(format_fields)):
            if len(sample_infos) > i:
                sample_dict[format_fields[i]] = sample_infos[i]

        # Check if GT not null
        if sample_dict["GT"].replace("0", ".") not in ["", ".", "./.", ".|."]:
            nb_pipeline_find += 1

    return f"{nb_pipeline_find}/{nb_pipeline}"


def genotypeconcordance(df, samples: list = []):
    """
    The function checks the genotype concordance of a given list of samples in a dataframe.

    :param df: The input dataframe containing genetic variant information, including genotype
    information for each sample/pipeline
    :param samples: The parameter "samples" is a list of sample/pipeline names that are present in the
    input dataframe "df". These samples/pipelines have genotype information that will be used to
    calculate genotype concordance
    :type samples: list
    :return: a string that indicates whether the genotypes of the specified samples in the input
    dataframe are concordant or not. The string is either "TRUE" or "FALSE", depending on whether all
    the specified samples have the same genotype or not.
    """

    # format
    format_fields = df["FORMAT"].split(":")

    # no sample/pipeline
    if not samples:
        return "0/0"

    # init
    genotype_list = {}

    # For each sample/pipeline
    for sample in samples:

        if sample in df:

            # Split snpeff ann values
            sample_infos = df[sample].split(":")

            # Create Dataframe
            sample_dict = {}
            for i in range(len(format_fields)):
                if len(sample_infos) > i:
                    sample_dict[format_fields[i]] = sample_infos[i]

            # Check if GT not null
            # genotype_list[sample_dict["GT"]] = 1
            if sample_dict["GT"] not in ["", ".", "./.", ".|."]:
                genotype_list[sample_dict["GT"]] = 1

    return str(len(genotype_list) == 1).upper()


def genotype_compression(genotype: str = "") -> str:
    """
    The function takes a genotype string, replaces dots with zeros, removes non-digit characters, sorts
    and removes duplicates, and returns the compressed genotype string.

    :param genotype: The input genotype as a string. It is a DNA sequence that contains genetic
    information
    :type genotype: str
    :return: The function `genotype_compression` returns a compressed version of the input genotype
    string. The compressed string has all dots replaced with 0s, all non-digit characters removed, and
    duplicates removed and sorted. The compressed string is returned as a string.
    """

    genotype_compressed = "".join(
        sorted(set(re.sub(r"\D", "", genotype.replace(".", "0"))))
    )

    return genotype_compressed


def genotype_barcode(genotype: str = "") -> str:
    """
    This function takes a genotype string and compresses it, then returns a barcode string based on the
    length and content of the compressed genotype.

    :param genotype: The genotype parameter is a string that represents a genetic sequence or code
    :type genotype: str
    :return: The function `genotype_barcode` returns a string representing the barcode for a given
    genotype. The barcode can be "0", "1", "2", or "?" depending on the length and content of the
    compressed genotype string.
    """

    genotype_compressed = genotype_compression(genotype)
    if len(genotype_compressed) == 1:
        if genotype_compressed == "0":
            barcode = "0"
        else:
            barcode = "2"
    elif len(genotype_compressed) > 1:
        barcode = "1"
    else:
        barcode = "?"

    return barcode


def barcode(df, samples: list = []):
    """
    Generates a barcode based on the genotype of the specified samples.

    :param df: A pandas DataFrame containing the genetic data.
    :type df: pandas.DataFrame

    :param samples: A list of sample names to use for generating the barcode.
    :type samples: list(str)

    :return: A barcode string based on the genotype of the specified samples.
    :rtype: str
    """
    # format
    format_fields = df["FORMAT"].split(":")

    # no sample/pipeline
    if not samples:
        return ""

    # init
    barcode = []

    # For each sample/pipeline
    for sample in samples:

        if sample in df:

            # Split snpeff ann values
            sample_infos = df[sample].split(":")

            # Create Dataframe
            sample_dict = {}
            for i in range(len(format_fields)):
                if len(sample_infos) > i:
                    sample_dict[format_fields[i]] = sample_infos[i]

            # generate barcode
            barcode.append(genotype_barcode(sample_dict["GT"]))

    return "".join(barcode)


def trio(df, samples: list = []):
    """
    The function trio(df, samples:list = []) determines the type of variant (denovo, dominant, or
    recessive) in a trio based on the barcode generated from the samples.

    :param df: The input dataframe containing genetic variant information
    :param samples: A list of sample IDs to be used in the analysis
    :type samples: list
    :return: The function `trio` returns a string that represents the type of variant in a trio
    analysis, which can be "denovo", "dominant", "recessive", or "unknown".
    """

    # no sample/pipeline
    if not samples:
        return ""

    # init
    trio_barcode = barcode(df, samples)

    # switcher
    switcher = {
        "001": "denovo",
        ("011", "101", "111", "021", "201", "121", "211"): "dominant",
        ("112", "212", "122", "222"): "recessive",
    }

    trio_variant_type = "unknown"
    for case in switcher:
        if trio_barcode in case:
            trio_variant_type = switcher[case]
            break

    return trio_variant_type


def vaf_normalization(row, sample: str) -> str:
    """
    This function takes in a row of data and a sample name, extracts the genotype information for that
    sample, calculates the variant allele frequency (VAF) from the genotype information, and adds the
    VAF to the genotype information before returning it.

    :param row: The input row of a pandas DataFrame containing information about a genetic variant
    :param sample: The parameter "sample" is a string representing the name of the sample for which we
    want to calculate the VAF (Variant Allele Frequency). It is used to extract the genotype information
    for that particular sample from the input row
    :type sample: str
    :return: a string that represents the genotype information for a given sample with an added "VAF"
    field that represents the variant allele frequency.
    """

    # format
    format_fields = row["FORMAT"].split(":")

    # Sample genotype
    sample_genotype = row[sample]

    # No genotype
    if not sample_genotype:
        sample_genotype = "./."

    # Split samples values values
    sample_genotype_infos = sample_genotype.split(":")

    # Create Dataframe
    sample_genotype_dict = {}
    for i in range(len(format_fields)):
        if len(sample_genotype_infos) > i:
            sample_genotype_dict[format_fields[i]] = sample_genotype_infos[i]
        else:
            sample_genotype_dict[format_fields[i]] = "."

    # Find VAF
    # Default VAF
    vaf = "."
    # VAF from FREQ
    if "FREQ" in sample_genotype_dict:
        if sample_genotype_dict["FREQ"] != ".":
            vaf_freq = sum(
                map(float, sample_genotype_dict["FREQ"].replace("%", "").split(","))
            )
            if vaf_freq:
                vaf = round(vaf_freq / 100, 4)
    # VAF from DP4
    elif "DP4" in sample_genotype_dict:
        if sample_genotype_dict["DP4"] != ".":
            dp4_split = sample_genotype_dict["DP4"].split(",")
            if dp4_split != ["."]:
                dp4_dp = sum(map(int, dp4_split))
                pd4_alt = sum(map(int, dp4_split[2:]))
                vaf = round(pd4_alt / dp4_dp, 6) if dp4_dp else "."
    # VAF from AD
    elif "AD" in sample_genotype_dict:
        if sample_genotype_dict["AD"] != ".":
            ad_split = sample_genotype_dict["AD"].split(",")
            if ad_split != ["."]:
                ad_dp = sum(map(int, ad_split))
                ad_alt = sum(map(int, ad_split[1:]))
                vaf = round(ad_alt / ad_dp, 6) if ad_dp else "."

    # add vaf into genotype
    sample_genotype_dict["VAF"] = vaf

    # Create new genotype info
    genotype_with_vaf = ":".join([str(v) for v in sample_genotype_dict.values()])

    return genotype_with_vaf


def genotype_stats(df, samples: list = [], info: str = "VAF"):
    """
    This function computes statistics on a specified information field (e.g. VAF) for a given set of
    samples in a pandas dataframe.

    :param df: The input dataframe containing variant information
    :param samples: The list of sample/pipeline names for which to compute the genotype statistics. If
    empty, the function will return an empty dictionary
    :type samples: list
    :param info: The parameter "info" is a string that represents the type of information being analyzed
    in the function. In this case, it is used to compute statistics on the Variant Allele Frequency
    (VAF) of genetic variants, defaults to VAF
    :type info: str (optional)
    :return: a dictionary containing statistics related to a specified information field (default is
    "VAF") for a given set of samples in a pandas DataFrame. The statistics include the number of
    values, a list of values, minimum value, maximum value, mean, median, and standard deviation. If no
    samples are specified, an empty dictionary is returned.
    """

    # format
    format_fields = df["FORMAT"].split(":")

    # init
    vaf_stats = {
        info + "_stats_nb": 0,
        info + "_stats_list": None,
        info + "_stats_min": None,
        info + "_stats_max": None,
        info + "_stats_mean": None,
        info + "_stats_mediane": None,
        info + "_stats_stdev": None,
    }

    # no sample/pipeline
    if not samples:
        return vaf_stats

    # init
    vaf_list = []

    # For each sample/pipeline
    for sample in samples:

        # Split snpeff ann values
        sample_infos = df[sample].split(":")

        # Create Dataframe
        sample_dict = {}
        for i in range(len(format_fields)):
            if len(sample_infos) > i:
                sample_dict[format_fields[i]] = sample_infos[i]

        # Check if GT not null
        if info in sample_dict:
            try:
                vaf_float = float(sample_dict[info])
            except:
                vaf_float = None
            if vaf_float:
                vaf_list.append(vaf_float)

    vaf_stats[info + "_stats_nb"] = len(vaf_list)
    if vaf_list:
        vaf_stats[info + "_stats_list"] = ":".join([str(x) for x in vaf_list])

    # Compute min, max, mean and median only if the list is not empty
    if vaf_list:
        vaf_stats[info + "_stats_min"] = min(vaf_list)
        vaf_stats[info + "_stats_max"] = max(vaf_list)
        vaf_stats[info + "_stats_mean"] = statistics.mean(vaf_list)
        vaf_stats[info + "_stats_mediane"] = statistics.median(vaf_list)

    # Check if there are at least 2 values in the list before computing variance or stdev
    if len(vaf_list) >= 2:
        vaf_stats[info + "_stats_stdev"] = statistics.stdev(vaf_list)

    return vaf_stats


def extract_file(file_path: str, path: str = None, threads: int = 1):
    """
    The function extracts a compressed file in .zip or .gz format based on the file path provided.

    :param file_path: The file path parameter is a string that represents the path to a file that needs
    to be extracted. The function checks if the file has a ".zip" or ".gz" extension and extracts it
    accordingly
    :type file_path: str
    :param path: The `path` parameter is an optional string that represents the directory where the
    extracted files will be saved. If no `path` is provided, the function will use the directory of the
    `file_path` as the extraction destination
    :type path: str
    :param threads: The `threads` parameter is an optional parameter that specifies the number of
    threads to use for extraction. By default, it is set to 1, meaning the extraction will be done using
    a single thread, defaults to 1
    :type threads: int (optional)
    """

    if file_path.endswith(".zip"):
        if not path:
            path = os.path.dirname(file_path)
        with zipfile.ZipFile(file_path, "r") as zip_ref:
            zip_ref.extractall(path)
    elif file_path.endswith(".gz"):
        concat_and_compress_files(
            input_files=[file_path],
            output_file=file_path[:-3],
            compression_type="none",
            threads=threads,
        )


def download_file(
    url: str,
    dest_file_path: str,
    chunk_size: int = 1024 * 1024,
    try_aria: bool = True,
    aria_async_dns: bool = False,
    threads: int = 1,
    quiet: bool = True,
):
    """
    The `download_file` function is a Python function that downloads a file from a given URL and saves
    it to a specified destination file path in chunks.

    :param url: The `url` parameter is the URL of the file you want to download. It should be a string
    that represents the complete URL, including the protocol (e.g., "http://example.com/file.txt")
    :type url: str
    :param dest_file_path: The `dest_file_path` parameter is the path where the downloaded file will be
    saved. It should be a string representing the file path, including the file name and extension. For
    example, if you want to save the file as "myfile.txt" in the current directory, you can set `dest
    :type dest_file_path: str
    :param chunk_size: The `chunk_size` parameter determines the size of each chunk of data that is
    downloaded at a time. In this case, the default value is set to 1 MB, which means that the file will
    be downloaded in chunks of 1 MB at a time. This parameter can be adjusted according to
    :type chunk_size: int
    :param try_aria: The `try_aria` parameter is a boolean value that determines whether to use the
    Aria2c command-line tool for downloading the file. If set to `True`, the function will attempt to
    download the file using Aria2c. If set to `False`, the function will use the, defaults to True
    :type try_aria: bool (optional)
    :param aria_async_dns: The `aria_async_dns` parameter is a boolean value that determines whether to
    use asynchronous DNS resolution with Aria2c. If set to `True`, Aria2c will use asynchronous DNS
    resolution, which can improve download performance. If set to `False`, Aria2c will use synchronous,
    defaults to False
    :type aria_async_dns: bool (optional)
    :param threads: The `threads` parameter specifies the number of threads to be used for downloading
    the file. It determines the number of simultaneous connections that will be made to download the
    file. By default, it is set to 1, which means that only one connection will be made at a time.
    Increasing the value, defaults to 1
    :type threads: int (optional)
    :param quiet: The `quiet` parameter is a boolean value that determines whether to suppress the
    output of the download process. If set to `True`, the output will be suppressed. If set to `False`,
    the output will be displayed. By default, it is set to `True`, defaults to True
    :type quiet: bool (optional)
    :return: a boolean value indicating whether the file was successfully downloaded and saved to the
    specified destination file path.
    """

    # Full path
    dest_file_path = full_path(dest_file_path)

    # Create folder if not exists
    if not os.path.exists(os.path.dirname(dest_file_path)):
        Path(os.path.dirname(dest_file_path)).mkdir(parents=True, exist_ok=True)

    try:

        if try_aria:

            # Aria options
            aria_async_dns_option = str(aria_async_dns).lower()
            if quiet and log.root.level >= 20:
                aria_quiet_option = " --quiet "
                aria_redirect_option = " 2>/dev/null "
            else:
                aria_quiet_option = ""
                aria_redirect_option = ""

            # Aria command
            aria_command = f"aria2c -c -s {threads} -x {threads} -j 1 {url} -d {os.path.dirname(dest_file_path)} -o {os.path.basename(dest_file_path)} {aria_quiet_option}"

            # Launch command
            # Try with --async-dns option
            try:
                output = os.system(
                    f"{aria_command} --async-dns={aria_async_dns_option} {aria_redirect_option}"
                )
                if output:
                    assert False
            except:
                output = os.system(aria_command)

            # Test output file
            if output or not (
                os.path.exists(dest_file_path) and os.stat(dest_file_path).st_size > 0
            ):
                assert False

        else:

            assert False

    except:

        # Request
        with requests.get(url, stream=True) as r:
            r.raise_for_status()

            # Create a temporary file
            tmp_file_path = dest_file_path + ".tmp"

            # Open the temporary file for writing in binary mode
            with open(tmp_file_path, "wb") as f:

                # Download the file by chunks
                for chunk in r.iter_content(chunk_size=chunk_size):

                    # Write the chunk to the temporary file
                    f.write(chunk)

        # Move the temporary file to the final destination
        shutil.move(tmp_file_path, dest_file_path)

    return os.path.exists(dest_file_path)


def whereis_bin(bin_file: str) -> str:
    """ """

    if isinstance(bin_file, str):
        whereis_bin = bin_file
    else:
        return None

    # Switch to which
    if whereis_bin and which(whereis_bin):
        whereis_bin = which(whereis_bin)

    # Full path
    whereis_bin = full_path(whereis_bin)

    if whereis_bin and (
        os.path.exists(whereis_bin) or (whereis_bin and which(whereis_bin))
    ):
        return whereis_bin
    else:
        return None


def get_bin(
    bin: str = None,
    tool: str = None,
    bin_type: str = "bin",
    config: dict = {},
    default_folder: str = DEFAULT_TOOLS_FOLDER,
    output_type: str = "bin",
) -> str:
    """
    The `get_bin` function retrieves the path to a specified binary file from a configuration dictionary
    or searches for it in the file system if it is not specified in the configuration.

    :param bin: The `bin` parameter is a string or a pattern that represents the name of the binary file (e.g.,
    `snpEff.jar`, `exomiser-cli*.jar`) that you want to retrieve the path for
    :type bin: str
    :param tool: The `tool` parameter is a string that represents the name of the tool. It is used to
    retrieve the path to the tool's binary file
    :type tool: str
    :param bin_type: The `bin_type` parameter is a string that specifies the type of binary file to
    search for in the config dict (e.g., `jar`, `bin`). In this case, the default value is "bin". A value "jar" indicates that the function is searching
    for a JAR file. Defaults to bin
    :type bin_type: str (optional)
    :param config: A dictionary containing configuration information for the snpEff tool, including the
    path to the snpEff jar file. If no configuration is provided, an empty dictionary is used
    :type config: dict
    :param default_folder: The `default_folder` parameter is a string that represents the default folder
    where the tool binaries are located. If the `bin_file` is not found in the configuration dictionary
    or in the file system, the function will search for it in this default folder
    :type default_folder: str
    :return: the path to the snpEff.jar file. If the file is not found, it returns None.
    """

    # Full path
    default_folder = full_path(default_folder)
    # log.debug(f"default_folder={default_folder}")

    # Config
    config_tool = config.get("tools", {}).get(tool)

    # Allowed dict conf
    tool_dict_conf_type_allowed = ["jar", "java", "docker"]

    # Find bin_file
    bin_file = None
    bin_file_type = bin_type

    # For input config for tool and default config (if not provided)
    for conf in [config_tool, DEFAULT_TOOLS_BIN.get(tool, {})]:

        # If config not empty and bin_file not already found
        if conf and not bin_file:

            # If tool config is a dict (simple dict with type of more complex such as docker)
            if isinstance(conf, dict) and bin_type:

                # If type found in dict (simple dict such as
                # {'bin': '/path/to/tool'})
                # {'jar': '/path/to/tool.jar'})
                if bin_type in conf:
                    bin_file_new = conf.get(bin_type)
                    # If bin is a env binary or a existing file (see whereis function)
                    if whereis_bin(bin_file_new):
                        bin_file = whereis_bin(bin_file_new)
                    # If specific bin type is asked (e.g. docker, which is a dict not a string)
                    elif bin_type in tool_dict_conf_type_allowed and isinstance(
                        bin_file_new, dict
                    ):
                        if output_type in ["bin", "dict"]:
                            return bin_file_new
                        elif output_type in ["type"]:
                            return bin_type

                # If bin not found as a specific input type
                if not bin_file:

                    # For each other type in conf (if any)
                    for bin_type_conf in conf:

                        # Bin file/config and bin type
                        bin_file = conf.get(bin_type_conf)
                        bin_file_type = bin_type_conf

                        # If type is allowed to be return (e.g. 'docker', 'java')
                        if bin_file_type in tool_dict_conf_type_allowed:

                            # return if output type is for a dict (for futher configuration, e.g. for 'docker', 'java')
                            if output_type in ["dict"] and isinstance(bin_file, dict):
                                return bin_file
                            # else:
                            #     return None

            # If config is a path to binary/file as a string
            elif isinstance(conf, str):
                bin_file = conf

    # Return not "bin"
    if output_type in ["type"]:
        return bin_file_type
    elif output_type in ["dict"] and isinstance(bin_file, dict):
        return bin_file

    # Check bin_file validation
    bin_file = whereis_bin(bin_file)

    # Config - check tools
    if not bin_file:

        # Try to find bin file
        try:
            bin_file = find_all(bin, default_folder)[0]
        except:
            return None

        # Full path
        bin_file = full_path(bin_file)

        # Check if found
        if not os.path.exists(bin_file):
            return None

    return bin_file


def get_bin_command(
    bin: str = None,
    tool: str = None,
    bin_type: str = "bin",
    config: dict = {},
    param: dict = {},
    default_folder: str = DEFAULT_TOOLS_FOLDER,
    add_options: str = None,
) -> str:
    """
    The function `get_bin_command` generates a command based on the tool type (jar, java, docker) and
    specified parameters.

    :param bin: The `bin` parameter in the `get_bin_command` function is used to specify the binary
    executable file that you want to run. It is a string that represents the path or name of the binary
    file. If you provide this parameter, the function will attempt to locate the binary file based on
    the
    :type bin: str
    :param tool: The `tool` parameter in the `get_bin_command` function represents the name of the tool
    for which you want to retrieve the command. It is used to identify the specific tool for which the
    command is being generated
    :type tool: str
    :param bin_type: The `bin_type` parameter in the `get_bin_command` function specifies the type of
    binary executable that the tool uses. It can have values like "bin", "jar", "java", "docker", etc.,
    depending on the type of tool being executed. The function uses this parameter to determine,
    defaults to bin
    :type bin_type: str (optional)
    :param config: The `config` parameter in the `get_bin_command` function is a dictionary that holds
    configuration settings for the tool being used. It can include various settings such as paths,
    environment variables, or any other configuration options needed for the tool to run properly
    :type config: dict
    :param param: The `param` parameter in the `get_bin_command` function is a dictionary that contains
    additional parameters or configurations for the tool being executed. These parameters can be used to
    customize the behavior or settings of the tool when generating the command for execution. The
    function uses the `param` dictionary along with the
    :type param: dict
    :param default_folder: The `default_folder` parameter in the `get_bin_command` function is used to
    specify the default folder where the tools are located. If a specific folder is not provided when
    calling the function, it will default to the value of `DEFAULT_TOOLS_FOLDER`
    :type default_folder: str
    :param add_options: The `add_options` parameter in the `get_bin_command` function allows you to pass
    additional options or arguments to the command being constructed based on the tool type. These
    additional options can be specific configurations, flags, or any other parameters that you want to
    include in the final command. When provided,
    :type add_options: str
    :return: The `get_bin_command` function returns a string representing the command to execute a
    specific tool based on the provided parameters. The returned command can be either a Java command
    for running a JAR file or a Docker command for running a Docker image/container. If the tool type is
    not Java or Docker, it returns the default tool bin.
    """

    # Tool bin and type
    tool_bin = get_bin(
        bin=bin,
        tool=tool,
        bin_type=bin_type,
        config=config,
        default_folder=default_folder,
    )
    tool_bin_type = get_bin(
        bin=bin,
        tool=tool,
        bin_type=bin_type,
        config=config,
        output_type="type",
        default_folder=default_folder,
    )

    # Threads and memory
    threads = get_threads(config=config, param=param)
    memory = extract_memory_in_go(get_memory(config=config, param=param))
    tmp = get_tmp(config=config, param=param)
    # Java and jar command
    if tool_bin_type in ["jar", "java"]:

        # Find java bin
        java_bin = get_bin(tool="java", config=config)

        # Java options
        java_options = ""

        # Threads
        if threads:
            java_options += f" -XX:ParallelGCThreads={threads} "

        # Memory
        if memory:
            java_options += f" -XX:MaxHeapSize={memory}G "

        # Additional options
        if add_options:
            java_options += f" {add_options} "

        # Java Command for tool
        tool_command = f"{java_bin} {java_options} -jar {tool_bin}"
        # Return tool java command
        return tool_command

    # Docker command
    elif tool_bin_type in ["docker"]:

        # Find docker bin
        docker_bin = get_bin(tool="docker", config=config)

        # Docker param for tool
        tool_bin_dict = get_bin(
            tool=tool, config=config, bin_type="docker", output_type="dict"
        )
        tool_bin_docker_image = tool_bin_dict.get("image", None)
        tool_bin_docker_entrypoint = tool_bin_dict.get("entrypoint", None)
        tool_bin_docker_command = tool_bin_dict.get("command", "")
        tool_bin_docker_options = tool_bin_dict.get("options", "")
        if not tool_bin_docker_image:
            msg_error = f"Error in Docker command for tool '{tool}'"
            log.error(msg_error)
            raise ValueError(msg_error)

        # Create default param
        tool_bin_docker_params = ""
        if (
            config.get("tools", {})
            .get(tool, {})
            .get("docker", {})
            .get("config", {})
            .get("notremove", "")
            is True
        ):
            tool_bin_docker_params += ""
        elif config.get("docker", {}).get("notremove", "") is True:
            tool_bin_docker_params += ""
        else:
            tool_bin_docker_params += " --rm "

        # Temp folder by default
        if (
            config.get("tools", {})
            .get(tool, {})
            .get("docker", {})
            .get("config", {})
            .get("tmp", "")
            is True
        ):
            tool_bin_docker_params += f" -v {tmp}:{tmp} "
        elif config.get("docker", {}).get("tmp", "") is True:
            tool_bin_docker_params += f" -v {tmp}:{tmp} "
        else:
            tool_bin_docker_params += ""

        # Threads
        if threads:
            tool_bin_docker_params += f" --cpus={threads} "

        # Memory
        if memory:
            tool_bin_docker_params += f" --memory={memory}g "

        def mount_folder(config):
            add_options = ""
            for databases, path in (
                config.get("folders", {}).get("databases", {}).items()
            ):
                if path and isinstance(path, str):
                    add_options += f" -v {path}:{path}:ro"
                elif path and isinstance(path, list):
                    for pathlist in path:
                        add_options += f" -v {pathlist}:{pathlist}:ro"
            return add_options

        # Init docker param
        if tool_bin_docker_options:
            tool_bin_docker_params += f" {tool_bin_docker_options} "

        # Mount folder only
        if (
            config.get("tools", {})
            .get(tool, {})
            .get("docker", {})
            .get("config", {})
            .get("automount", "")
            is False
        ):
            tool_bin_docker_params += mount_folder(config)
        elif (
            config.get("tools", {})
            .get(tool, {})
            .get("docker", {})
            .get("config", {})
            .get("automount", "")
            is True
        ):
            tool_bin_docker_params += docker_automount()
        elif config.get("docker", {}).get("automount", "") is False:
            tool_bin_docker_params += mount_folder(config)
        elif config.get("docker", {}).get("automount", "") is True:
            tool_bin_docker_params += docker_automount()
        else:
            tool_bin_docker_params += mount_folder(config)

        # Entrypoint
        if tool_bin_docker_entrypoint:
            tool_bin_docker_params += f" --entrypoint='{tool_bin_docker_entrypoint}' "

        # Additional options
        if add_options:
            tool_bin_docker_params += f" {add_options} "

        # Additioanl command
        if not tool_bin_docker_command:
            tool_bin_docker_command = ""

        # Docker command
        tool_command = f"{docker_bin} run {tool_bin_docker_params} {tool_bin_docker_image} {tool_bin_docker_command}"

        # Return tool docker command
        return tool_command

    else:
        return tool_bin


def get_tmp(config: dict = {}, param: dict = None, default_tmp: str = "/tmp") -> str:
    """
    The `get_tmp` function returns the value of the "tmp" parameter from either the `param` dictionary,
    `config` dictionary, or a default value "/tmp".

    :param config: Config is a dictionary that contains configuration settings for the function. It is
    an optional parameter with a default value of an empty dictionary. It can be used to provide
    additional configuration settings to the function `get_tmp`
    :type config: dict
    :param param: The `param` parameter is a dictionary containing parameters that can be passed to the
    function `get_tmp`. It can include various key-value pairs, but in this context, the function
    specifically looks for the key "tmp" within the `param` dictionary to determine the temporary path
    value. If the "
    :type param: dict
    :param default_tmp: The `default_tmp` parameter in the `get_tmp` function is a string that
    represents the default path for temporary files. If the "tmp" key is not found in the `param`
    dictionary or the `config` dictionary, the function will return this `default_tmp` value, which is,
    defaults to /tmp
    :type default_tmp: str (optional)
    :return: The function `get_tmp` returns the value of the "tmp" key from the `param` dictionary if it
    exists. If the "tmp" key is not found in the `param` dictionary, it returns the value of the "tmp"
    key from the `config` dictionary. If neither key is found in `param` or `config`, it returns the
    default value "/tmp".
    """

    # Tmp in param or config
    tmp_param = param.get("tmp", config.get("tmp", default_tmp))
    if not tmp_param:
        tmp_param = default_tmp

    # Return tmp
    return tmp_param


def get_threads(config: dict = {}, param: dict = {}) -> int:
    """
    This Python function retrieves the number of threads to use based on input parameters and system
    configuration.

    :param config: The `config` parameter is a dictionary that contains configuration settings for the
    function `get_threads`. It can be used to provide default values for the number of threads to use in
    the function
    :type config: dict
    :param param: The `param` parameter is a dictionary that may contain the key "threads" which
    specifies the number of threads to use. If the "threads" key is not present in the `param`
    dictionary, the function will look for the "threads" key in the `config` dictionary. If neither
    :type param: dict
    :return: The function `get_threads` returns the number of threads to be used based on the input
    parameters.
    """

    # Config
    if not config:
        config = {}

    # Param
    if not param:
        param = {}

    # Threads system
    nb_threads = os.cpu_count()

    # Threads in param or config
    threads_param = int(param.get("threads", config.get("threads", 0)))

    # Check threads in params
    if threads_param:
        input_thread = threads_param
    else:
        input_thread = None

    # Check threads value
    if not threads_param or int(threads_param) <= 0:
        threads = nb_threads
    else:
        threads = int(input_thread)

    # Return threads
    return threads


def get_memory(config: dict = {}, param: dict = None) -> str:
    """
    The `get_memory` function retrieves memory information using psutil and calculates a default memory
    value based on total memory, with the option to specify a custom memory value.

    :param config: The `config` parameter is a dictionary that may contain configuration settings for
    the function `get_memory`. It is used to provide default values or settings for the function
    :type config: dict
    :param param: The `param` parameter is a dictionary that may contain a key "memory" which represents
    the amount of memory to be used. If the "memory" key is not present in the `param` dictionary, the
    function will try to retrieve the value from the `config` dictionary using the key "
    :type param: dict
    :return: The function `get_memory` returns a string representing the amount of memory to be used.
    This memory value is calculated based on the total memory available on the system, with a default
    value set to 80% of the total memory. The function first checks if a specific memory value is
    provided in the `param` dictionary, and if not, it looks for a default value in the `config`
    """

    import psutil

    # Memory system
    mem = psutil.virtual_memory()
    mem_total = mem.total / 1024 / 1024 / 1024
    mem_default = int(mem_total * 0.8)

    # Threads in param or config
    if config:
        memory_config = config.get("memory", None)
    else:
        memory_config = None
    if param:
        memory_param = param.get("memory", memory_config)
    else:
        memory_param = memory_config

    # Check memory value
    if mem_default < 1:
        mem_default = 1
    if memory_param:
        memory = memory_param
    else:
        memory = f"{mem_default}G"

    return memory


def extract_float_from_str(text: str = "") -> float:
    """
    The function `extract_float_from_str` extracts a float value from a given string input.

    :param text: The `extract_float_from_str` function is designed to extract a floating-point number
    from a given string input. The function uses a regular expression to find the first occurrence of a
    floating-point number in the input string and returns it as a float
    :type text: str
    :return: The function `extract_float_from_str` returns a float value extracted from the input text
    string. If a float value is found in the text, it is returned as a float. If no float value is
    found, it returns `None`.
    """

    value_in_float = re.findall(r"[-+]?\d*\.*\d+", text)

    if len(value_in_float):
        return float(value_in_float[0])
    else:
        return None


def extract_memory_in_go(
    memory_str, default_value: int = 1, default_unit: str = "G"
) -> int:
    """
    The `extract_memory_in_go` function converts a memory size string in the format FLOAT[kMG] to an
    integer value in Go memory units.

    :param memory_str: The `memory_str` parameter should be a string representing a memory value with a
    unit suffix in the format FLOAT[kMG]. For example, it could be "1G", "512M", or "2k"
    :param default: The `default` parameter in the `extract_memory_in_go` function is used to specify a
    default integer value if the conversion of the memory size string fails or if the value cannot be
    extracted from the input string. If no valid value can be extracted from the input string, the
    function will return the, defaults to 1
    :type default: int (optional)
    :return: The `extract_memory_in_go` function is returning an integer value representing the memory
    size in Go units based on the input memory string provided.
    """

    # Convert in upper ccase
    memory_str = str(memory_str).upper()

    # Coefficient
    coefficient = {
        "K": 1 / (1024 * 1024),
        "M": 1 / 1024,
        "G": 1,
    }

    try:
        # Extract in float
        value = extract_float_from_str(text=memory_str)
        if not value:
            value = default_value

        # Extract unit
        unit = re.findall(r"[?([A-Za-z]+", memory_str)
        if not unit:
            unit = default_unit
        else:
            unit = unit[0]

        # Convert in int Go
        if unit in coefficient:
            return int(value * coefficient[unit])
        else:
            raise ValueError(
                f"Invalid unit '{unit}'. Use 'k', 'm' ou 'g' i upper or lower case."
            )
    except ValueError:
        raise ValueError(
            f"Invalid memory format '{memory_str}'. Use format FLOAT[kMG]."
        )


def concat_file(input_files: list, output_file: str) -> bool:
    """
    This function concatenates multiple input files into a single output file.

    :param input_files: A list of file paths to the input files that need to be concatenated
    :type input_files: list
    :param output_file: The parameter "output_file" is a string that represents the name of the file
    that will be created by the function and will contain the concatenated content of all the input
    files
    :type output_file: str
    :return: a boolean value indicating whether the output file was successfully created or not. It
    checks if the output file exists using the `os.path.exists()` function and returns `True` if it
    exists and `False` otherwise.
    """

    # Full path
    output_file = full_path(output_file)

    with open(output_file, "w") as outfile:
        for input_file in input_files:
            input_file = full_path(input_file)
            with open(input_file) as infile:
                for line in infile:
                    outfile.write(line)

    return os.path.exists(output_file)


def compress_file(input_file: str, output_file: str) -> bool:
    """
    This function compresses a file using the BGZF compression algorithm.

    :param input_file: The path and name of the input file that needs to be compressed
    :type input_file: str
    :param output_file: The output_file parameter is a string that represents the name and path of
    the file where the compressed data will be written
    :type output_file: str
    """

    # Full path
    output_file = full_path(output_file)

    # Concat and compress
    concat_and_compress_files(input_files=[input_file], output_file=output_file)

    # Return
    return os.path.exists(output_file)


def get_compression_type(filepath: str) -> str:
    """
    The function `get_compression_type` determines the compression type of a file based on its first few
    bytes.

    :param filepath: The `filepath` parameter is a string that represents the path to the file for which
    we want to determine the compression type
    :type filepath: str
    :return: The function `get_compression_type` returns a string indicating the compression type of the
    file specified by the `filepath` parameter. The possible return values are "gzip" if the file is
    compressed using gzip, "bgzip" if the file is compressed using bgzip, "unknown" if the compression
    type is unknown, and "none" if the file is not compressed.
    """

    try:
        with open(filepath, "rb") as test_f:
            # Firsts bits
            bit_1 = test_f.read(1)
            bit_2 = test_f.read(1)
            bit_3 = test_f.read(1)
            bit_4 = test_f.read(1)
            # If bit is compress
            if bit_1 == b"\x1f" and bit_2 == b"\x8b" and bit_3 == b"\x08":
                # If bit 4 == x04 is compress 'bgzip' type
                if bit_4 == b"\x04":
                    return "bgzip"
                # else is compress 'gzip' type
                else:
                    return "gzip"
            # If no compress type
            else:
                return "none"
    except:
        return "unknown"


def get_file_compressed(filename: str = None) -> bool:
    """
    This function takes a filename as input and returns True if the file is compressed (in bgzip) and False if it
    is not

    :param filename: the name of the file to be checked
    :type filename: str
    :return: A boolean value.
    """

    # Full path
    filename = full_path(filename)

    compression_type = None
    if filename:
        # If filename exists, try to detect compression within content (first bits)
        if os.path.exists(filename):
            compression_type = get_compression_type(filename)
        # If filename does not exists, use extension to determine if it is compressed
        else:
            filename_name, filename_extension = os.path.splitext(filename)
            compress_format = filename_extension.replace(".", "")
            if compress_format in file_compressed_format:
                compression_type = "bgzip"
            else:
                compression_type = "none"
    return compression_type in ["gzip", "bgzip"]


def concat_into_infile(
    input_files: list,
    compressed_file: object,
    compression_type: str = "none",
    threads: int = 1,
    block: int = 10**6,
) -> bool:
    """
    The function `concat_into_infile` concatenates multiple input files into a compressed output file,
    with support for different compression types and multi-threading.

    :param input_files: A list of input file paths that need to be concatenated into the compressed file
    :type input_files: list
    :param compressed_file: The `compressed_file` parameter is an object that represents the file where
    the concatenated contents of the input files will be written. It is expected to be a file object
    that has write capabilities
    :type compressed_file: object
    :param compression_type: The `compression_type` parameter specifies the type of compression to be
    used for the output file. The default value is "none", which means no compression will be applied.
    Other possible values include "bgzip" and "gzip", which indicate that the output file should be
    compressed using the bgzip and, defaults to none
    :type compression_type: str (optional)
    :param threads: The "threads" parameter specifies the number of threads to use for compression or
    decompression. It determines how many parallel processes can be executed simultaneously, which can
    help improve performance when dealing with large files or multiple files, defaults to 1
    :type threads: int (optional)
    :param block: The `block` parameter is used to specify the size of the block when reading the input
    files. It is set to `10 ** 6`, which means 1 million bytes. This parameter determines how much data
    is read from the input files at a time
    :type block: int
    :return: a boolean value, specifically `True`.
    """

    # Output file compressions type
    if compression_type in ["none"]:
        open_type = "t"
    else:
        open_type = "b"

    # Open input files
    for input_file in input_files:

        # Full path
        input_file = full_path(input_file)

        # Input file compression type
        input_compression_type = get_compression_type(input_file)
        if input_compression_type in ["bgzip"]:
            with open(input_file, "rb") as raw:
                with bgzip.BGZipReader(
                    raw, num_threads=threads, raw_read_chunk_size=block
                ) as infile:
                    if open_type in ["t"]:
                        compressed_file.write(str(infile.read(), "utf-8"))
                    else:
                        shutil.copyfileobj(infile, compressed_file)
        elif input_compression_type in ["gzip"]:
            # See https://pypi.org/project/mgzip/
            with mgzip.open(input_file, "r" + open_type, thread=threads) as infile:
                shutil.copyfileobj(infile, compressed_file)
        elif input_compression_type in ["none"]:
            with open(input_file, "r" + open_type) as infile:
                shutil.copyfileobj(infile, compressed_file)
        else:
            raise ValueError(f"Input file compression type unknown: {input_file}")

    return True


def concat_and_compress_files(
    input_files: list,
    output_file: str,
    compression_type: str = "bgzip",
    threads: int = 1,
    memory: int = 1,
    block: int = 10**6,
    compression_level: int = 6,
    sort: bool = False,
    index: bool = False,
) -> bool:
    """
    The function `concat_and_compress_files` takes a list of input files, an output file name, and
    optional parameters for compression type, number of threads, block size, compression level, sorting,
    and indexing, and concatenates and compresses the input files into the output file.

    :param input_files: A list of input file paths that need to be concatenated and compressed
    :type input_files: list
    :param output_file: The `output_file` parameter is a string that specifies the path and name of the
    output file that will be created after concatenating and compressing the input files
    :type output_file: str
    :param compression_type: The `compression_type` parameter specifies the type of compression to be
    applied to the output file. It can take one of three values: "bgzip", "gzip", or "none", defaults to
    bgzip
    :type compression_type: str (optional)
    :param threads: The `threads` parameter specifies the number of threads to use for compression and
    decompression. It determines the level of parallelism in the compression process, allowing for
    faster execution when multiple threads are used, defaults to 1
    :type threads: int (optional)
    :param memory: The `memory` parameter specifies the amount of max memory (in Gb) to use for sorting.
    defaults to 1
    :type memory: int (optional)
    :param block: The `block` parameter specifies the size of the block used for reading and writing
    data during compression. It is set to a default value of 10^6 (1 million) bytes
    :type block: int
    :param compression_level: The `compression_level` parameter determines the level of compression to
    be used when compressing the output file. It is an integer value ranging from 0 to 9, where 0
    indicates no compression and 9 indicates maximum compression. The higher the compression level, the
    smaller the resulting compressed file size, defaults to 6
    :type compression_level: int (optional)
    :param sort: The `sort` parameter is a boolean flag that determines whether the output file should
    be sorted or not. If `sort` is set to `True`, the output file will be sorted using
    `pysam.bcftools.sort` before renaming it. If `sort` is set to `False, defaults to False
    :type sort: bool (optional)
    :param index: The `index` parameter is a boolean flag that determines whether or not to index the
    output file after concatenation and compression. If `index` is set to `True`, the output file will
    be indexed using the `pysam.tabix_index` function with the preset "vcf". Make sure VCF is sorted.
    Defaults to False
    :type index: bool (optional)
    :return: a boolean value indicating whether the output file exists or not.
    """

    # Prevent compression type not available
    if compression_type not in ["bgzip", "gzip", "none"]:
        compression_type = "bgzip"

    # Output file compressions type
    if compression_type in ["none"]:
        open_type = "t"
    else:
        open_type = "b"

    # Full path
    output_file = full_path(output_file)

    # Tmp file
    output_file_tmp = output_file + "." + str(random.randrange(1000)) + ".tmp"

    if compression_type in ["gzip"]:
        # See https://pypi.org/project/mgzip/
        with pgzip.open(
            output_file_tmp,
            "w" + open_type,
            thread=threads,
            blocksize=threads * block,
            compresslevel=compression_level,
        ) as compressed_file:
            concat_into_infile(
                input_files,
                compressed_file,
                compression_type=compression_type,
                threads=threads,
                block=block,
            )
    elif compression_type in ["bgzip"]:
        with open(output_file_tmp, "w" + open_type) as compressed_file_raw:
            with bgzip.BGZipWriter(
                compressed_file_raw, num_threads=threads
            ) as compressed_file:
                concat_into_infile(
                    input_files,
                    compressed_file,
                    compression_type=compression_type,
                    threads=threads,
                    block=block,
                )
    elif compression_type in ["none"]:
        with open(output_file_tmp, "w" + open_type) as compressed_file:
            concat_into_infile(
                input_files,
                compressed_file,
                compression_type=compression_type,
                threads=threads,
                block=block,
            )

    # Output file
    if sort:
        # Sort with pysam
        try:
            pysam.bcftools.sort(
                f"-Oz{compression_level}",
                "-o",
                output_file,
                "-T",
                output_file_tmp,
                output_file_tmp,
                "-m",
                f"{memory}G",
                threads=threads,
                catch_stdout=False,
            )
            # Remove tmp file
            os.remove(output_file_tmp)
        except:
            raise ValueError(f"Output file sorting failed: {output_file_tmp}")

    else:
        # Rename tmp file
        os.rename(output_file_tmp, output_file)

    # Index file
    if index:
        # Index file with tabix
        try:
            pysam.tabix_index(output_file, preset="vcf", force=True)
        except:
            raise ValueError(f"Output file indexing failed: {output_file}")

    # Return output file
    return os.path.exists(output_file)


def get_plateform_name_from_duckdb(conn: duckdb.DuckDBPyConnection) -> str:
    """
    The function `get_plateform_name_from_duckdb` returns the platform information from a DuckDB connection.

    :param conn: The `conn` parameter is an instance of the `DuckDBPyConnection` class from the `duckdb`
    module. It represents a connection to a DuckDB database
    :type conn: duckdb.DuckDBPyConnection
    :return: the platform information from the DuckDB connection.
    """

    return conn.query(f"PRAGMA platform").df()["platform"][0]


def get_duckdb_extension_file(
    extension_name: str, conn: duckdb.DuckDBPyConnection, download: bool = True
) -> str:
    """
    This function returns the file path of a DuckDB extension based on the extension name and platform.

    :param extension_name: The name of the DuckDB extension file that is being requested
    :type extension_name: str
    :return: a string that represents the file path of a DuckDB extension file. The file path is
    constructed using the constant `DUCKDB_EXTENSION`, the platform name obtained from the
    `get_plateform_name_from_duckdb()` function, and the extension name passed as an argument to the function.
    """

    # Init
    release_version_number = duckdb.__version__
    platform_name = get_plateform_name_from_duckdb(conn)

    # File
    extension_file_path = (
        f"v{release_version_number}/{platform_name}/{extension_name}.duckdb_extension"
    )
    extension_file = f"{DUCKDB_EXTENSION}/{extension_file_path}"
    extension_file = full_path(extension_file)
    extension_file_gz = f"{extension_file}.gz"
    url = f"http://extensions.duckdb.org/{extension_file_path}.gz"

    # Download and extract if not exists
    if not os.path.exists(extension_file) and download:
        # Create folder if not exists
        if not os.path.exists(os.path.dirname(extension_file)):
            Path(os.path.dirname(extension_file)).mkdir(parents=True, exist_ok=True)
        # Download extension if not exists
        if not os.path.exists(extension_file_gz):
            try:
                download_file(url, extension_file_gz)
            except:
                log.error(f"Fail download '{url}'")
        # Uncompress extention file
        if os.path.exists(extension_file_gz):
            with mgzip.open(extension_file_gz, "rb", thread=1) as infile:
                with open(extension_file, "wb") as compressed_file:
                    shutil.copyfileobj(infile, compressed_file)

    if os.path.exists(extension_file):
        return extension_file
    else:
        return None


def load_duckdb_extension(
    conn: duckdb.DuckDBPyConnection, duckdb_extensions: list
) -> bool:
    """
    This function loads DuckDB extensions into a connection object and returns a boolean indicating
    whether all extensions were successfully loaded.

    :param conn: duckdb.DuckDBPyConnection object representing a connection to a DuckDB database
    :type conn: duckdb.DuckDBPyConnection
    :param duckdb_extensions: A list of strings representing the names of the DuckDB extensions to be
    loaded
    :type duckdb_extensions: list
    :return: a boolean value indicating whether all the specified DuckDB extensions were successfully
    loaded or not.
    """

    loaded = True

    for extension_name in duckdb_extensions:
        duckdb_extension_file = get_duckdb_extension_file(
            extension_name, conn=conn, download=True
        )
        try:
            if duckdb_extension_file and os.path.exists(duckdb_extension_file):
                # Try loading extension by file
                conn.query(
                    f"INSTALL '{duckdb_extension_file}'; LOAD {extension_name}; "
                )
            else:
                try:
                    with time_limit(60):
                        # Try loading extension by name
                        conn.query(
                            f"INSTALL '{extension_name}'; LOAD {extension_name}; "
                        )
                except:
                    loaded = False
        except:
            try:
                with time_limit(60):
                    # Try loading extension by name
                    conn.query(f"INSTALL '{extension_name}'; LOAD {extension_name}; ")
            except:
                loaded = False

    return loaded


class TimeoutException(Exception):
    pass


def time_limit(seconds):
    def signal_handler(signum, frame):
        raise TimeoutException("Timed out!")

    signal.signal(signal.SIGALRM, signal_handler)
    signal.alarm(seconds)
    try:
        yield
    finally:
        signal.alarm(0)


def duckdb_execute(query: str, threads: int = 1) -> bool:
    """
    The `duckdb_execute` function executes a query using the DuckDB database engine and returns a
    boolean indicating whether the query was successful or not.

    :param query: The `query` parameter is a string that represents the SQL query you want to execute in
    DuckDB. It can be any valid SQL statement, such as SELECT, INSERT, UPDATE, DELETE, etc
    :type query: str
    :param threads: The "threads" parameter specifies the number of threads to use for executing the
    query. By default, it is set to 1, meaning that the query will be executed using a single thread,
    defaults to 1
    :type threads: int (optional)
    :return: The function `duckdb_execute` returns a boolean value. It returns `True` if the query
    execution is successful, and `False` if it is not successful.
    """

    conn = duckdb.connect(config={"threads": threads})
    conn.execute("SET max_expression_depth TO 10000")
    if conn.execute(query):
        conn.close()
        return True
    else:
        conn.close()
        return False


def genome_build_switch(assembly: str) -> str:
    """
    The `genome_build_switch` function takes an assembly name as input and returns a new
    assembly name if a different version of the same genome is available, otherwise it returns
    None.

    :param assembly: The `assembly` parameter is a string that represents the name or identifier
    of a genome assembly
    :type assembly: str
    :return: The function `genome_build_switch` returns a string.
    """

    import genomepy

    genome_list = genomepy.search(assembly, exact=False)

    if genome_list:
        for row in genome_list:
            tax_id = row[3]
            break
    else:
        log.error(f"Assembly '{assembly}' NOT found")
        return None

    genome_list = genomepy.search(assembly, exact=False)
    for row in genome_list:
        new_assembly = row[0].split(".")[0]
        new_tax_id = row[3]
        if new_assembly != assembly and tax_id == new_tax_id:
            return new_assembly

    return None


# get argument
def get_argument(
    arguments: dict = {},
    arg: str = "",
    required: bool = False,
    remove_infos: list = ["gooey", "extra"],
    add_metavar: bool = False,
) -> dict:
    """
    The `get_argument` function retrieves information about a specific argument from a dictionary, and
    can also set its "required" status.

    :param arguments: A dictionary containing information about the arguments passed to a function or
    method
    :type arguments: dict
    :param arg: The `arg` parameter is a string that represents the name of the argument that you want
    to retrieve information for
    :type arg: str
    :param required: The `required` parameter is a boolean value that determines whether the argument is
    required or not. If set to True, the function will return an empty dictionary if the argument is not
    found in the `arguments` dictionary. If set to False (default), the function will still return an
    empty dictionary if, defaults to False
    :type required: bool (optional)
    :param remove_infos: The `remove_infos` parameter is a list that contains the names of specific
    information that you want to remove from the argument dictionary. In the code, it is used to remove
    specific argument information such as "gooey" from the `arg_infos` dictionary
    :type remove_infos: list
    :return: a dictionary containing information about a specific argument, specified by the `arg`
    parameter. If the argument is found in the `arguments` dictionary, the function returns a dictionary
    containing the information about that argument. If the argument is not found, an empty dictionary is
    returned. The `required` parameter is used to specify whether the argument is required or not, and
    this information is added to
    """

    if arg in arguments:
        arg_infos = arguments.get(arg, {}).copy()
        for arg_info in remove_infos:
            arg_infos.pop(arg_info, None)
        if required != None:
            arg_infos["required"] = required
        if add_metavar and "metavar" not in arg_infos:
            arg_infos["metavar"] = arg.replace("_", " ")
        return arg_infos
    else:
        return {}


# get_argument_gooey
def get_argument_gooey(arguments: dict = {}, arg: str = ""):
    """
    The function `get_argument_gooey` takes an argument and returns the corresponding widget and options
    for the Gooey library in Python.

    :param arg: The `arg` parameter is a string that represents the name of the argument you want to
    retrieve information for
    :type arg: str
    :return: The function `get_argument_gooey` returns two values: `widget` and `options`.
    """

    # Init
    argument = get_argument(arguments=arguments, arg=arg, remove_infos=[])
    argument_type = argument.get("type", None)
    gooey_argument = argument.get("gooey", {})

    # Widget
    widget = None
    if gooey_argument.get("widget", None):
        widget = gooey_argument.get("widget")
    else:
        if str(argument_type) == "FileType('r')":
            widget = "FileChooser"
        elif str(argument_type) == "FileType('w')":
            widget = "FileSaver"

    # options
    options = gooey_argument.get("options", {})

    # Return
    return widget, options


# get argument
def get_argument_to_mk(arg: str, argument: dict = {}, mode: str = "mk") -> str:
    """
    The function `get_argument_to_mk` generates a formatted string containing information about a
    command line argument, which can be output in either Markdown or HTML format.

    :param arg: The `arg` parameter is a string that represents the name of the argument. It is used to
    generate the header and text for the argument
    :type arg: str
    :param argument: The `argument` parameter is a dictionary that contains information about the
    argument. It has the following keys:
    :type argument: dict
    :param mode: The `mode` parameter is used to specify the format of the output. It can have two
    possible values: "mk" or "html". If "mk" is specified, the output will be formatted using Markdown
    syntax. If "html" is specified, the output will be formatted using HTML syntax, defaults to mk
    :type mode: str (optional)
    :return: a formatted string that provides information about a command line argument. The format of
    the string depends on the value of the `mode` parameter. If `mode` is set to "html", the string is
    formatted as an HTML `<pre>` block. Otherwise, the string is formatted as a Markdown code block. The
    string includes the argument name, metavariable, help text, required
    """

    from html import escape

    text = ""

    # Option info
    metavar = argument.get("metavar", arg)
    help = argument.get("help", None)
    required = argument.get("required", None)
    choices = argument.get("choices", None)
    default = argument.get("default", None)
    action = argument.get("action", None)

    # header
    text_header = f"--{arg}"
    if not action:
        text_header += f"=<{metavar}>"
    if choices:
        text_header += f" {choices}"
    if default:
        text_header += f" (default: {default})"
    if required:
        text_header += " | required"

    # text
    if mode == "html":
        text += f"<pre>"
        text += escape(text_header)
        text += "\n"
        text += escape(str(help))
        text += "\n\n"
        text += f"</pre>"
    else:
        text += "\n<small>\n"
        text += "\n> ```\n"
        text += f">     {text_header}"
        text += "\n> \n> "
        text += str(help).strip().replace("\n", "\n> ")
        text += "\n> ```\n"
        text += "\n</small>\n"
        text += "\n"
        text += ""

    return text


def help_generation_from_dict(
    element: str,
    help_dict: dict,
    previous: str = "",
    output_type: str = "markdown",
    level: int = 1,
    table: str = "",
    generate_table: bool = False,
    code_type: str = "",
    auto_default: bool = True,
    previous_sections: bool = False,
):
    """
    The `help_generation_from_dict` function generates help documentation from a dictionary input,
    supporting markdown and HTML output formats with specific sections like "__help", "__format",
    "__default", "__examples", "__code", and "__examples_code".

    :param element: The `element` parameter in the `help_generation_from_dict` function is a string that
    represents the current element or key in the dictionary for which help documentation is being
    generated. It is the specific key or element within the dictionary that you want to generate help
    documentation for
    :type element: str
    :param help_dict: The `help_dict` parameter in the `help_generation_from_dict` function is a
    dictionary that contains the help documentation for various elements or keys. This dictionary
    structure allows for organizing and storing information related to each element, such as help text,
    formatting details, default values, and examples. The function processes
    :type help_dict: dict
    :param previous: The `previous` parameter in the `help_generation_from_dict` function is used to
    keep track of the previous elements in the hierarchy. It is a string that represents the path to the
    current element being processed. This parameter helps in maintaining the correct hierarchy level
    when generating help documentation for nested elements in a
    :type previous: str
    :param output_type: The `output_type` parameter in the `help_generation_from_dict` function
    specifies the type of output format that you want the generated help documentation to be in. It can
    take two possible values: "markdown" or "html". By default, the output type is set to markdown,
    defaults to markdown
    :type output_type: str (optional)
    :param level: The `level` parameter in the `help_generation_from_dict` function is used to keep
    track of the depth or level of recursion in the generation process. It starts at 1 for the initial
    call and increments by 1 for each level of recursion into sub-elements. This parameter helps in
    formatting the, defaults to 1
    :type level: int (optional)
    :param table: The `table` parameter in the `help_generation_from_dict` function is used to store the
    table of contents for the generated help documentation. It is a string that contains the formatted
    table of contents with links to different sections or elements within the documentation. This table
    helps users navigate through the documentation easily
    :type table: str
    :param generate_table: The `generate_table` parameter in the `help_generation_from_dict` function is
    a boolean flag that determines whether the function should generate a table of contents for the help
    documentation. When set to `True`, the function will include a table of contents in the output based
    on the hierarchy of elements in the, defaults to False
    :type generate_table: bool (optional)
    :param code_type: The `code_type` parameter in the `help_generation_from_dict` function specifies
    the type of code examples that will be included in the generated help documentation. It defaults to
    "json", meaning that the code examples provided in the "__examples_code" section of the dictionary
    will be in JSON format
    :type code_type: str
    :param auto_default: The `auto_default` parameter in the `help_generation_from_dict` function is a
    boolean flag that determines whether the function should automatically populate certain sections of
    the help documentation based on the information available in the dictionary and the element's
    arguments. When set to `True`, the function will automatically fill in sections, defaults to True
    :type auto_default: bool (optional)
    :param previous_sections: The `previous_sections` parameter in the `help_generation_from_dict`
    function is a boolean flag that determines whether the function should include previous sections in
    the hierarchy when generating help documentation for nested elements. When set to `True`, the
    function will maintain the previous sections in the hierarchy path, helping to provide, defaults to
    False
    :type previous_sections: bool (optional)
    :return: The `help_generation_from_dict` function returns the generated help documentation based on
    the input `help_dict` dictionary. The output is formatted based on the specified `output_type`
    (either "markdown" or "html") and includes sections such as "__help", "__format", "__default", and
    "__examples" if they are present in the `help_dict`.
    """

    from howard.tools.tools import arguments

    element_argument_infos = arguments.get(
        element, arguments.get(element.replace("_", "-"), {})
    )

    # Level marker
    level_marker = "::"

    # Output
    output = ""
    output_table = ""

    # Specific help sections
    section_list = [
        "__help",
        "__type",
        "__choices",
        "__format",
        "__default",
        "__default_eval",
        "__examples",
        "__code",
        "__examples_code",
    ]

    # Previous level
    if previous and previous_sections:
        prefix = f"{previous}"
        previous = f"{prefix}{level_marker}"
    else:
        prefix = f"{element}"
        previous = ""

    if "__code_type" in help_dict:
        code_type = help_dict.get("__code_type", code_type)

    if "__auto" in help_dict:
        auto_default = help_dict.get("__auto", True)

    if not help_dict:
        auto_default = True

    if auto_default and level > 1:

        # default
        if "__help" not in help_dict:
            if "help" in element_argument_infos:
                help_dict["__help"] = element_argument_infos.get("help", None)  #
                # Clean help
                help_dict["__help"] = format_arg_help(help_dict["__help"])

        # default
        if "__default" not in help_dict:
            if "default" in element_argument_infos:
                help_dict["__default"] = element_argument_infos.get("default", None)

        # type
        if "__type" not in help_dict:
            if "type" in element_argument_infos:
                help_dict["__type"] = element_argument_infos.get("type", None).__name__

        # format
        if "__format" not in help_dict:
            if "format" in element_argument_infos:
                help_dict["__format"] = element_argument_infos.get("format", None)

        # choices
        if "__choices" not in help_dict:
            if "choices" in element_argument_infos:
                help_dict["__choices"] = element_argument_infos.get("choices", None)

        # format
        if "__format" not in help_dict:
            if (
                "extra" in element_argument_infos
                and "format" in element_argument_infos.get("extra", {})
            ):
                help_dict["__format"] = element_argument_infos.get("extra", {}).get(
                    "format", None
                )

        # example (code)
        if (
            "__examples" not in help_dict
            and "__examples_code" not in help_dict
            and "__code" not in help_dict
        ):
            if (
                "extra" in element_argument_infos
                and "examples" in element_argument_infos.get("extra", {})
            ):
                help_dict["__examples"] = element_argument_infos.get("extra", {}).get(
                    "examples", None
                )

    # If no section "__help" (mandatory)
    if "__help" not in help_dict:
        return output

    # For each help section
    for section in section_list:

        if section in help_dict:

            # If content not empty
            if section in help_dict:

                # element content
                help_dict_content = help_dict.get(section, None)

                # Break section
                section_break = ""

                # Output variables
                output_header = ""
                output_header_table = ""
                output_content = ""

                # Markdown
                if output_type == "markdown":

                    # Level
                    level_md = "#" * level

                    # line break
                    if section in ["__code", "__examples_code", "__examples"]:
                        line_break = "\n"
                    else:
                        line_break = "\n\n"

                    # Help content from type
                    if isinstance(help_dict_content, str):
                        help_md = help_dict_content.replace("\n", line_break)
                    # elif isinstance(help_dict_content, list) and section not in ["__choices"]:
                    elif (
                        isinstance(help_dict_content, list)
                        or isinstance(help_dict_content, set)
                    ) and section not in ["__choices"]:
                        help_md = line_break.join(help_dict_content)
                        section_break = "\n"
                    else:
                        help_md = help_dict_content

                    if section in ["__code", "__examples_code", "__examples"]:

                        # Format help dict
                        # if isinstance(help_md, set):
                        #     help_md = list(help_md)
                        if isinstance(help_md, str):
                            help_md_dict = {}
                            example_header = ""
                            for example_line in help_md.split("\n"):
                                if example_line.startswith("#"):
                                    example_header = example_line
                                else:

                                    if example_header not in help_md_dict:
                                        help_md_dict[example_header] = ""
                                        example_sep = ""
                                    help_md_dict[example_header] += (
                                        example_sep + example_line
                                    )
                                    example_sep = "\n"
                        elif isinstance(help_md, dict):
                            help_md_dict = help_md
                        else:
                            help_md_dict = help_md

                        # Format examples
                        help_md = ""
                        for example in help_md_dict:
                            example_code = help_md_dict.get(example, "")
                            if isinstance(example_code, list):
                                example_code = "\n".join(example_code)
                            example = re.sub(r"^#*\s*", "", example)
                            help_md += f"""\n\n> {example}\n"""
                            # Clean example code
                            if (
                                len(example_code) > 0
                                and code_type.upper() == "JSON"
                                and example_code.strip()[0] != "{"
                            ):
                                example_code = example_code.replace("\n", "\n   ")
                                example_code = f"""{{\n   {example_code}\n}}"""
                            example_code = example_code.replace("\n", "\n> ")
                            # Add example code
                            help_md += (
                                f"""\n> ```{code_type}\n> {example_code}\n> ```"""
                            )

                    if section in ["__default"]:
                        if help_md in [""]:
                            help_md = " "
                        help_md = f"```{help_md}```"

                    if section in ["__default_eval"]:
                        help_md = "```" + eval(help_md) + "```"

                    if section in ["__type"]:
                        help_md = f"```{help_md}```"

                    if section in ["__format"]:
                        help_md = f"```{help_md}```"

                    if section in ["__choices"]:
                        help_md = f"```{help_md}```"

                    # Output variables
                    output_header += f"{level_md} {prefix}\n\n"
                    output_header_table += (
                        f"#{level_md} Table of contents\n\n{table}\n\n"
                    )
                    output_content += f"{help_md}\n\n"

                # HTML
                elif output_type == "html":

                    # Level
                    level_html = f"H{level}"

                    # Help content from type
                    if isinstance(help_dict_content, str):
                        help_html = help_dict_content.replace("\n", "<br>")
                    elif isinstance(help_dict_content, list):
                        help_dict_content = map(str, help_dict_content)
                        if section in ["__code", "__examples_code", "__examples"]:
                            help_html = "\n".join(help_dict_content)
                            section_break = ""
                        else:
                            help_html = "<br>".join(help_dict_content)
                            section_break = "<br>"
                    elif isinstance(help_dict_content, dict):
                        if section in ["__code", "__examples_code", "__examples"]:
                            for example_header in help_dict_content:
                                help_html = f"{example_header}"
                                if isinstance(
                                    help_dict_content.get(example_header, []), list
                                ):
                                    help_html += "\n".join(
                                        help_dict_content.get(example_header, [])
                                    )
                                else:
                                    help_html += help_dict_content.get(
                                        example_header, ""
                                    )
                    else:
                        help_html = help_dict_content

                    if section in ["__code", "__examples_code", "__examples"]:
                        help_html = f"<xmp>{help_html}</xmp>"

                    # Output variables
                    output_header += (
                        f"<{level_html} id='{prefix}'>{prefix}</{level_html}>"
                    )
                    output_header_table += f"""
                    <H{level+1}>Table of contents</{level_html}>
                    {table}
                    <br>
                    """
                    output_content += f"{help_html}<br>"

                # Output section
                output_section = re.sub(r"^__", "", section).split("_")[0].capitalize()

                # Output construction
                if section == "__help":
                    output += output_header
                    output += output_content
                    if table and level == 1:
                        output += output_header_table
                else:
                    output_section = (
                        re.sub(r"^__", "", section).split("_")[0].capitalize()
                    )
                    output += f"{output_section}: {section_break}{output_content}"

    # Output table construction
    if output_type == "markdown":
        output_table += (
            "   " * (level - 1)
            + f"- [{element}](#{prefix.lower().replace(' ','-').replace(':','')})\n"
        )
    elif output_type == "html":
        output_table += (
            "&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;" * (level - 1)
            + f"- <a href='#{prefix}'>{element}</a><BR>"
        )

    # Recusivity to other sub-element
    for element in help_dict:

        # If element is a sub-element
        if element not in section_list and isinstance(
            help_dict.get(element, None), dict
        ):

            # output_tabl
            output_table += help_generation_from_dict(
                element=element,
                help_dict=help_dict.get(element, None),
                previous=f"{previous}{element}",
                output_type=output_type,
                level=level + 1,
                table=output_table,
                generate_table=True,
                code_type=code_type,
            )

            # Add sub-element output
            output += help_generation_from_dict(
                element=element,
                help_dict=help_dict.get(element, None),
                previous=f"{previous}{element}",
                output_type=output_type,
                level=level + 1,
                table=output_table,
                code_type=code_type,
            )

    # Return output
    if generate_table:
        return output_table
    else:
        return output


def help_generation_from_json(
    help_json_file: str,
    output_type: str = "markdown",
    title="Help",
    code_type: str = "",
    include_toc: bool = False,
):
    """
    The `help_generation_from_json` function reads a JSON file containing help information, converts it
    into a specified output format, and returns the generated help content.

    :param help_json_file: The `help_json_file` parameter is a string that should contain the file path
    to the JSON file from which help information will be extracted. This JSON file likely contains
    structured data that will be used to generate the help content
    :type help_json_file: str
    :param output_type: The `output_type` parameter in the `help_generation_from_json` function
    specifies the format in which the generated help content will be output. By default, it is set to
    "markdown", which means the help content will be formatted using Markdown syntax. However, you can
    also specify other output formats such, defaults to markdown
    :type output_type: str (optional)
    :param title: The `title` parameter in the `help_generation_from_json` function is a string that
    represents the title of the help documentation that will be generated. It is used to provide a title
    for the help content to make it more organized and informative. By default, the title is set to
    "Help", defaults to Help (optional)
    :param code_type: The `code_type` parameter in the `help_generation_from_json` function is used to
    specify the type of code examples that will be included in the generated help content. This
    parameter allows you to define the format or language of the code examples to be displayed alongside
    the help information extracted from the JSON
    :type code_type: str
    :param include_toc: The `include_toc` parameter in the `help_generation_from_json` function is a
    boolean flag that determines whether a table of contents (TOC) should be included in the generated
    help content. If `include_toc` is set to `True`, a table of contents will be generated based,
    defaults to False
    :type include_toc: bool (optional)
    :return: The function `help_generation_from_json` returns the generated help content based on the
    information stored in the JSON file provided as input.
    """

    # Read JSON file
    help_dict = {}
    with open(help_json_file) as f:
        help_dict = yaml.safe_load(f)

    # Generate help table of content
    if include_toc:
        output_table = help_generation_from_dict(
            element=title,
            help_dict=help_dict,
            previous="",
            output_type=output_type,
            level=1,
            generate_table=True,
            code_type=code_type,
        )
    else:
        output_table = None

    # Generate help
    output = help_generation_from_dict(
        element=title,
        help_dict=help_dict,
        previous="",
        output_type=output_type,
        level=1,
        table=output_table,
        code_type=code_type,
    )

    # Return help
    return output


class RawTextArgumentDefaultsHelpFormatter(
    argparse.ArgumentDefaultsHelpFormatter, argparse.RawTextHelpFormatter
):
    pass


def help_header(setup: str = None) -> str:
    """
    The `help_header` function generates a header for the help documentation based on the metadata
    information provided in the setup file.

    :param setup: The `setup` parameter is a string that represents the path to a configuration file.
    This file contains metadata about the program, such as its name, version, description, and long
    description content type
    :type setup: str
    :return: The function `help_header` returns a string that represents the header for the help
    documentation. The header includes the program name, version, authors, and description.

    """

    # Config Parser
    # setup = "/tmp/config"
    if os.path.isfile(setup):
        cf = ConfigParser()
        cf.read(setup)
        prog_name = cf.get("metadata", "name", fallback="Unknown Program")
        prog_version = cf.get("metadata", "version", fallback="0.0.0")
        prog_authors = cf.get("metadata", "author", fallback="0.0.0")
        prog_description = cf.get(
            "metadata", "description", fallback="No description available."
        )
    else:
        # meta
        try:
            meta = metadata("howard-ann")
            prog_name = meta.get("Name")
            prog_version = meta.get("Version")
            prog_authors = meta.get("author")
            prog_description = meta.get("Summary")
        # Default
        except Exception as e:
            print(f"Erreur : {e}")
            prog_name = "HOWARD"
            prog_version = "0.0.0"
            prog_authors = "ALB"
            prog_description = "HOWARD - Highly Open Workflow for Annotation & Ranking toward genomic variant Discovery"

    # Logo
    import pyfiglet

    ascii_logo = pyfiglet.figlet_format(
        prog_name.split("-")[0].upper()
    )  # , font="slant")

    # Header
    header = f"""{ascii_logo}{prog_name.split('-')[0].upper()}::{prog_version} [{prog_authors}]\n{prog_description}\n"""

    # Return
    return header


def help_generation(
    arguments_dict: dict = {},
    parser=None,
    output_type: str = "parser",
):
    """
    The `help_generation` function generates a parser object for command-line arguments, as well as
    markdown or HTML help documentation for those arguments.

    :param arguments_dict: A dictionary containing the arguments for the function. It has three keys:
    :type arguments_dict: dict
    :param parser: The `parser` parameter is an instance of the `argparse.ArgumentParser` class. It is
    used to define the command-line interface and parse the command-line arguments. If no `parser` is
    provided, a new instance of `argparse.ArgumentParser` will be created
    :param setup: The `setup` parameter is a string that represents the path to a configuration file.
    This file contains metadata about the program, such as its name, version, description, and long
    description content type
    :type setup: str
    :param output_type: The `output_type` parameter determines the format of the output. It can be one
    of the following values:, defaults to parser
    :type output_type: str (optional)
    :return: The function `help_generation` returns different outputs based on the value of the
    `output_type` parameter.
    """

    # Arguments
    arguments = arguments_dict.get("arguments", {})
    commands_arguments = arguments_dict.get("commands_arguments", {})
    shared_arguments = arguments_dict.get("shared_arguments", {})

    # Parser default
    if not parser:
        parser = argparse.ArgumentParser()

    # Parser information
    parser.prog = metadata("howard-ann").get("Name").split("-")[0]
    parser.description = ""
    parser.epilog = (
        "\nUsage examples:\n"
        + """   howard process --input=tests/data/example.vcf.gz --output=/tmp/example.annotated.vcf.gz --param=config/param.json \n"""
        + """   howard annotation --input=tests/data/example.vcf.gz --output=/tmp/example.howard.vcf.gz --annotations='tests/databases/annotations/current/hg19/dbnsfp42a.parquet,tests/databases/annotations/current/hg19/gnomad211_genome.parquet' \n"""
        + """   howard calculation --input=tests/data/example.full.vcf --output=/tmp/example.calculation.tsv --calculations='vartype' \n"""
        + """   howard prioritization --input=tests/data/example.vcf.gz --output=/tmp/example.prioritized.vcf.gz --prioritization_config=config/prioritization_profiles.json --prioritizations='default,GERMLINE' \n"""
        + """   howard query --input=tests/data/example.vcf.gz --explode_infos --query='SELECT "#CHROM", POS, REF, ALT, "DP", "CLNSIG", sample2, sample3 FROM variants WHERE "DP" >= 50 OR "CLNSIG" NOT NULL ORDER BY "CLNSIG" DESC, "DP" DESC' \n"""
        + """   howard stats --input=tests/data/example.vcf.gz \n"""
        + """   howard convert --input=tests/data/example.vcf.gz --output=/tmp/example.tsv --explode_infos && cat /tmp/example.tsv \n"""
    )
    parser.formatter_class = argparse.RawTextHelpFormatter

    # Optionals
    parser._optionals.title = "Shared arguments"

    # Sub parser
    subparsers = parser.add_subparsers(title="Tools", dest="command")

    # Help options
    options_md = ""
    options_html = ""

    # Options for markdown
    options_md += f"# HOWARD Help"
    options_md += "\n\n"
    options_md += f"## Introduction"
    options_md += "\n\n"
    options_md += "<!--TOC-->"
    options_md += "\n\n"
    options_md += parser.description.replace("\n", "\n\n")
    options_md += re.sub(
        r"> $", "", parser.epilog.replace("\n", "\n\n").replace("   ", "> ")
    )
    options_md += "\n\n"

    # Options for HTML
    options_html += f"<H1>HOWARD Help</h1>\n"
    options_html += "<p>" + parser.description.replace("\n", "<br>") + "</p>"
    options_html += parser.epilog.replace("\n", "<br>").replace(
        "   ", "&nbsp;&nbsp;&nbsp;"
    )

    # Create commands arguments
    for command in commands_arguments:

        # Command description
        command_description = commands_arguments[command].get("description", "")
        command_help = commands_arguments[command].get("help", "")
        command_epilog = commands_arguments[command].get("epilog", "")

        # Command parser
        command_parser = subparsers.add_parser(
            command,
            description=command_description,
            help=command_help,
            epilog=command_epilog,
            # formatter_class=argparse.RawTextHelpFormatter
            # formatter_class=argparse.ArgumentDefaultsHelpFormatter
            formatter_class=RawTextArgumentDefaultsHelpFormatter,
        )

        # Markdown
        options_md += f"## {command.upper()} tool\n"
        options_md += command_description.replace("\n", "\n\n")
        options_md += "\n\n"
        options_md += re.sub(
            r"> $", "", command_epilog.replace("\n", "\n\n").replace("   ", "> ")
        )
        options_md += "\n\n"

        # HTML
        options_html += f"<H2>{command.upper()}</H2>\n"
        options_html += "<p>" + command_description.replace("\n", "<br>") + "</p>"
        options_html += command_epilog.replace("\n", "<br>").replace(
            "   ", "&nbsp;&nbsp;&nbsp;"
        )

        # Gooey - add metavar
        if output_type == "gooey":
            add_metavar = True
            command_group_suffix = ""
        else:
            add_metavar = False
            command_group_suffix = " options"

        # Main args
        command_parser._optionals.title = "Options"
        if "main" in commands_arguments[command]["groups"]:
            group = "main"
            options_md += f"### Main options\n"
            options_html += f"<H3>Main options</H3>\n"
            for arg in commands_arguments[command]["groups"][group]:
                required = commands_arguments[command]["groups"][group][arg]
                argument = get_argument(
                    arguments=arguments.copy(),
                    arg=arg,
                    required=required,
                    add_metavar=add_metavar,
                )
                if output_type == "gooey":
                    widget, options = get_argument_gooey(arguments=arguments, arg=arg)
                    argument["widget"] = widget
                    argument["gooey_options"] = options
                    if argument.get("help", "") in ["==SUPPRESS=="]:
                        argument["help"] = arg
                    argument["help"] = format_arg_help(
                        argument.get("help", ""), str(argument.get("default", None))
                    )
                command_parser.add_argument(f"--{arg}", **argument)
                options_md += get_argument_to_mk(arg, argument)
                options_html += get_argument_to_mk(arg, argument, mode="html")

        for group in commands_arguments[command]["groups"]:
            if group != "main":
                options_md += f"### {group}\n"
                options_html += f"<H3>{group}</H3>\n"
                command_group = command_parser.add_argument_group(
                    f"{group}{command_group_suffix}"
                )
                for arg in commands_arguments[command]["groups"][group]:
                    required = commands_arguments[command]["groups"][group][arg]
                    argument = get_argument(
                        arguments=arguments.copy(),
                        arg=arg,
                        required=required,
                        add_metavar=add_metavar,
                    )
                    if output_type == "gooey":
                        widget, options = get_argument_gooey(
                            arguments=arguments, arg=arg
                        )
                        argument["widget"] = widget
                        argument["gooey_options"] = options
                        if argument.get("help", "") in ["==SUPPRESS=="]:
                            argument["help"] = arg
                        argument["help"] = format_arg_help(
                            argument.get("help", ""), str(argument.get("default", None))
                        )
                    command_group.add_argument(f"--{arg}", **argument)
                    options_md += get_argument_to_mk(arg, argument)
                    options_html += get_argument_to_mk(arg, argument, mode="html")

        # Shared arguments
        shared_group = command_parser.add_argument_group("Shared options")
        for arg in shared_arguments:
            argument = get_argument(
                arguments=arguments, arg=arg, required=False, add_metavar=add_metavar
            )
            if output_type == "gooey":
                widget, options = get_argument_gooey(arguments=arguments, arg=arg)
                argument["widget"] = widget
                argument["gooey_options"] = options
                if argument.get("help", "") in ["==SUPPRESS=="]:
                    argument["help"] = arg
                argument["help"] = format_arg_help(
                    argument.get("help", ""), str(argument.get("default", None))
                )
            shared_group.add_argument(f"--{arg}", **argument)

        options_md += "\n\n"

    # Shared arguments for help files
    options_md += f"## Shared arguments\n"
    options_html += f"<H2>Shared arguments</H2>\n".upper()
    for arg in shared_arguments:
        required = False
        argument = get_argument(arguments=arguments.copy(), arg=arg, required=required)
        options_md += get_argument_to_mk(arg, argument)
        options_html += get_argument_to_mk(arg, argument, mode="html")

    # Output
    if output_type in ["parser", "gooey"]:
        return parser
    elif output_type == "markdown":
        return options_md
    elif output_type == "html":
        return options_html
    else:
        return parser


def format_arg_help(help_message: str, default_value: object = None) -> str:
    """
    The function `format_arg_help` formats a help message for a function argument, including a default
    value if provided.

    :param help_message: The `help_message` parameter is a string that contains the description or help
    message for a function or method argument. It provides information about the purpose or usage of the
    argument
    :type help_message: str
    :param default_value: The `default_value` parameter in the `format_arg_help` function is an optional
    parameter that specifies a default value for the argument being described in the help message. If a
    default value is provided, it will be included in the formatted help message to indicate the default
    value for that argument
    :type default_value: object
    :return: The function `format_arg_help` returns a formatted help message with a default value
    appended at the end if provided.
    """

    help_return = help_message

    random_string = "".join(
        random.choices(string.ascii_uppercase + string.digits, k=10)
    )
    help_return = (
        re.sub(r"\n\s*-", random_string, help_return)
        .replace("\n", " ")
        .replace(random_string, "\n-")
    )
    if default_value:
        help_return += "\n(default: " + str(default_value) + ")"

    return help_return


def bed_sort(input: str, output: str) -> str:
    """
    The `bed_sort` function reads a tab-separated input file, sorts the data based on columns 0, 1, and
    2 in ascending order, and writes the sorted data to a tab-separated output file.

    :param input: The `input` parameter is the path to the input file that contains the data to be
    sorted. This file should be in a tab-separated format
    :type input: str
    :param output: The `output` parameter is a string that specifies the path and filename of the output
    file where the sorted data will be saved
    :type output: str
    """

    # Read input BED (uncompressed)
    # data = pd.read_csv(input, sep="\t", dtype={"chr": "string", "start": "int", "end": "int"}, header=None)
    data = pd.read_csv(input, sep="\t", header=None, low_memory=False)

    # Sort BED regions
    data_sorted = data.sort_values(
        [0, 1, 2], axis=0, ascending=True, na_position="first"
    )

    # Write output BED (uncompressed)
    data_sorted.to_csv(output, sep="\t", header=0, index=0)


def full_path(path: str) -> str:
    """
    The function `full_path` takes a path string as input and returns the full expanded path.

    :param path: The `full_path` function takes a string `path` as input and returns the full path by
    expanding the user's home directory in the path if it is not None
    :type path: str
    :return: The function `full_path` is returning the expanded version of the input `path` using
    `os.path.expanduser(path)`. This function expands the `~` character in the path to the user's home
    directory.
    """

    if isinstance(path, str):
        return os.path.abspath(os.path.expandvars(os.path.expanduser(path)))
    else:
        return path


def get_default_argument(arguments_dict: dict, argument: str):
    """
    The function `get_default_argument` retrieves the default value of a specified argument from a
    dictionary of arguments.

    :param arguments_dict: The `arguments_dict` parameter is a dictionary that contains information
    about arguments
    :type arguments_dict: dict
    :param argument: The `get_default_argument` function takes in two parameters:
    :type argument: str
    :return: The function is attempting to return the default value of a specific argument from a
    dictionary of arguments. However, there is a mistake in the code. The correct key to access the
    argument's default value should be "argument" instead of "arguments". Therefore, the function will
    return the default value of the specified argument if it exists, otherwise it will return None.
    """

    return arguments_dict.get("arguments", {}).get(argument, {}).get("default", None)


def set_param(
    param: dict,
    args: argparse,
    arguments_dict: dict,
    argument: str,
    section: list = None,
) -> dict:
    """
    The function `set_param` takes input arguments and adds them to a dictionary based on certain
    conditions.

    :param param: The `param` parameter is a dictionary that stores configuration parameters or
    settings. It is used to collect and store various arguments and their values based on the conditions
    specified in the `set_param` function
    :type param: dict
    :param args: The `args` parameter in the `set_param` function is likely an instance of the
    `argparse.Namespace` class, which is typically used to store the command-line arguments parsed by
    the `argparse` module in Python. It contains the values of the arguments provided by the user when
    the script
    :type args: argparse
    :param arguments_dict: The `arguments_dict` parameter seems to be a dictionary that likely contains
    information about arguments and their default values. This dictionary is used in the function
    `set_param` to determine whether a specific argument should be included in the `param` dictionary
    based on certain conditions
    :type arguments_dict: dict
    :param argument: The `argument` parameter in the `set_param` function represents the specific
    argument that you want to set in the `param` dictionary. It is the key that will be used to store
    the value in the dictionary
    :type argument: str
    :param section: The `section` parameter in the `set_param` function is used to specify a section
    within the `param` dictionary where the argument value should be stored. If a `section` is provided,
    the argument value will be stored under that section in the `param` dictionary. If no `section
    :type section: str
    :return: the updated `param` dictionary after setting the specified argument value based on the
    conditions provided in the function.
    """

    # Argument value
    value = vars(args).get(argument, None)

    # Sections
    if section:
        if isinstance(section, str):
            sections = section.split(":")
        else:
            sections = section
    else:
        sections = []

    # Check if to include in param
    if argument in args and (
        not arguments_dict
        or value
        not in [get_default_argument(arguments_dict=arguments_dict, argument=argument)]
    ):
        sections.append(argument)
        param = add_value_into_dict(dict_tree=param, sections=sections, value=value)

    return param


def add_value_into_dict(dict_tree: dict, sections: list = [], value=None):
    """
    The function `add_value_into_dict` adds a value into a dictionary tree based on the provided
    sections.

    :param dict_tree: The `dict_tree` parameter is a dictionary representing a tree structure. It serves
    as the starting point for adding a value based on the provided sections
    :type dict_tree: dict
    :param sections: The `sections` parameter in the `add_value_into_dict` function represents a list of
    sections corresponding to successive keys in the dictionary. These sections are used to traverse the
    dictionary tree and determine the location where the value should be added. Each element in the
    `sections` list corresponds to a key in
    :type sections: list
    :param value: The `value` parameter in the `add_value_into_dict` function represents the value that
    you want to add into the dictionary tree at the specified location determined by the `sections`
    list. This value can be of any data type (e.g., int, str, list, dict, etc.)
    :return: The function `add_value_into_dict` returns the updated dictionary tree after adding the
    value based on the given sections.
    """

    # Pointer to traverse the tree
    current_node = dict_tree

    # If sections is empty, add the value to the root of the dictionary tree
    if not sections:
        return current_node

    # Traverse the tree based on sections to find the appropriate location for the value
    for section in sections[:-1]:
        # If the section does not exist yet in the current dictionary, create it
        if section not in current_node:
            current_node[section] = {}
        # Move down the tree following the current section
        current_node = current_node[section]

    # Add the value to the last section
    last_section = sections[-1]
    current_node[last_section] = value

    return dict_tree


def load_param(args: argparse) -> dict:
    """
    The function `load_param` takes command line arguments and returns a dictionary containing
    parameters loaded from a file or as JSON.

    :param args: It seems like the code snippet you provided is a function named `load_param` that takes
    an argument `args` of type `argparse` and returns a dictionary. The function is intended to load
    parameters from a file or a string
    :type args: argparse
    :return: A dictionary containing the loaded parameters is being returned.
    """

    param = {}
    if "param" in args:
        if isinstance(args.param, str) and os.path.exists(full_path(args.param)):
            with open(full_path(args.param)) as param_file:
                param = yaml.safe_load(param_file)

        else:
            param = json.loads(args.param)

    return param


def load_config_args(args):
    """
    The function `load_config_args` takes in arguments, extracts specific keys from them, and loads
    parameters in JSON format.

    :param args: The `load_config_args` function takes in an `args` object as input. This `args` object
    seems to contain various configuration parameters that the function will use to load and return
    specific values
    :return: The function `load_config_args` returns the variables `arguments_dict`, `setup_cfg`,
    `config`, and `param`.
    """

    # Arguments dict
    if "arguments_dict" in args:
        arguments_dict = args.arguments_dict
    else:
        arguments_dict = None

    # Setup config
    if "setup_cfg" in args:
        setup_cfg = args.setup_cfg
    else:
        setup_cfg = None

    # Config
    if "config" in args:
        config = args.config
    else:
        config = {}

    # Load parameters in JSON format
    param = load_param(args)

    return arguments_dict, setup_cfg, config, param


def load_args(
    param: dict,
    args: argparse,
    arguments_dict: dict,
    command: str = None,
    arguments_list: dict = {},
    strict: bool = False,
    section_prefix: list = [],
) -> dict:
    """
    The `load_args` function processes arguments based on specified parameters and conditions, raising
    an error if a specified argument is not found.

    :param param: The `param` parameter in the `load_args` function is a dictionary that stores the
    arguments and their values. It is used to keep track of the arguments that have been loaded or
    processed during the argument parsing process
    :type param: dict
    :param args: The `args` parameter in the `load_args` function is an instance of the
    `argparse.ArgumentParser` class from the `argparse` module in Python. This object is used to parse
    command-line arguments and options. It contains information about the arguments passed to the script
    when it was executed
    :type args: argparse
    :param arguments_dict: The `arguments_dict` parameter in the `load_args` function is a dictionary
    that likely contains information about the arguments expected by the script. It may include details
    such as the argument names, their corresponding sections, and any additional parameters related to
    each argument. This dictionary is used within the `load_args
    :type arguments_dict: dict
    :param command: The `command` parameter in the `load_args` function is a string that represents a
    specific command or action for which arguments need to be loaded. This parameter is used to identify
    the command-specific arguments that should be processed during argument parsing
    :type command: str
    :param arguments_list: The `arguments_list` parameter in the `load_args` function is a dictionary
    that contains the names of arguments that are expected to be present in the `args` object. This list
    is used to specify which arguments should be processed by the function `load_args` during the
    argument parsing process
    :type arguments_list: dict
    :param strict: The `strict` parameter in the `load_args` function is a boolean flag that determines
    whether an error should be raised if an argument specified in the `arguments_list` list is not found
    in the `args` object. If `strict` is set to `True`, an error will be raised, defaults to False
    :type strict: bool (optional)
    :return: The function `load_args` is returning a dictionary named `param` after processing the
    arguments based on the input parameters and conditions specified in the function.
    """

    # Variables
    arguments_list_to_load = {}
    param_section_not_found = get_random(N=16)

    # List from command
    if command:
        command_infos = arguments_dict.get("commands_arguments", {}).get(command, {})
        for command_group in command_infos.get("groups", {}):
            for command_argument in command_infos.get("groups").get(command_group):
                command_group_clean = command_group.replace(" ", "_").lower()
                if command_group_clean in ["main"]:
                    command_group_clean = None
                command_argument = command_argument.replace("-", "_")
                arguments_list_to_load[command_argument] = command_group_clean

    # Add arguments from arguments
    if arguments_list:
        for argument in arguments_list:
            arguments_list_to_load[argument] = arguments_list.get(argument)

    # Load arguments
    if arguments_list_to_load:
        for argument in arguments_list_to_load:
            if argument in args:
                section = (
                    arguments_dict.get("arguments", {})
                    .get(argument, {})
                    .get("extra", {})
                    .get("param_section", param_section_not_found)
                )
                if section in [param_section_not_found]:
                    section = arguments_list_to_load.get(argument, None)
                # section str to list
                if isinstance(section, str):
                    section = section.split(":")
                else:
                    section = []
                if section_prefix:
                    section = section_prefix + section
                param = set_param(
                    param=param,
                    args=args,
                    arguments_dict=arguments_dict,
                    section=section,
                    argument=argument,
                )

            elif strict:
                msg_error = f"Argument '{argument}' not found in list of aguments"
                log.error(msg_error)
                raise ValueError(msg_error)

    return param


def get_random(N: int = 10) -> str:
    """
    The function `get_random` generates a random string of uppercase letters and digits with a default
    length of 10.

    :param N: The parameter `N` in the `get_random` function represents the length of the random string
    that will be generated. By default, if no value is provided for `N`, it will generate a random
    string of length 10 consisting of uppercase letters and digits, defaults to 10
    :type N: int (optional)
    :return: A random string of length N consisting of uppercase letters and digits.
    """

    return "".join(random.choices(string.ascii_uppercase + string.digits, k=N))


def transcripts_file_to_df(
    transcripts_file: str, column_names: list = ["transcript", "gene"]
) -> pd.DataFrame:
    """
    This Python function reads a transcripts file into a pandas DataFrame, filtering out comment lines.

    :param transcripts_file: The `transcripts_file` parameter is a string that represents the file path
    to a file containing transcript information. This function is designed to read the contents of this
    file and convert it into a pandas DataFrame. The file is expected to be tab-separated with two
    columns: "transcript" and "gene
    :type transcripts_file: str
    :param column_names: The `column_names` parameter is a list that specifies the column names expected
    in the transcripts file. By default, it is set to `["transcript", "gene"]`, indicating that the file
    should have two columns named "transcript" and "gene". If the actual column names in the
    :type column_names: list
    :return: A pandas DataFrame containing transcript and gene information read from the specified file
    after filtering out comment lines is being returned.
    """

    # Full path
    transcripts_file = full_path(transcripts_file)

    # Transcript dataframe
    transcripts_dataframe = pd.DataFrame(columns=column_names)

    # If file exists
    if transcripts_file and os.path.exists(transcripts_file):

        # Data
        data = []

        # Check column names length by reading all lines
        with open(transcripts_file, "r") as f:
            for line in f:
                # Split line by tab
                values = line.strip().split("\t")
                data.append(values)

        # Convert into DataFrame
        transcripts_dataframe_nb_columns = pd.DataFrame(data)

        # If dataframe not with the same len than input name columns
        if len(transcripts_dataframe_nb_columns.columns) != len(column_names):
            column_names_new = []
            for i in range(len(transcripts_dataframe_nb_columns.columns)):
                if len(column_names) > i:
                    column_names_new.append(column_names[i])
                else:
                    column_names_new.append(f"column_{i+1}")
            column_names = column_names_new

        # Create dataframe
        transcripts_dataframe = pd.read_csv(
            transcripts_file,
            sep="\t",
            header=None,
            names=column_names,
            index_col=False,
        )

        # First column
        first_column = column_names[0]

        # Filter comment on lines
        transcripts_dataframe = transcripts_dataframe[
            transcripts_dataframe[first_column].str.contains("^[^#]+")
        ]

    # Return
    return transcripts_dataframe


def identical(
    vcf_list: typing.List[str], begin: str = "##", line_strip: bool = True
) -> bool:
    """
    The `identical` function compares the contents of multiple VCF files to determine if they are
    identical.

    :param vcf_list: The `vcf_list` parameter is a list of file paths to VCF (Variant Call Format) files
    that you want to compare for identity. The function reads the contents of these files and checks if
    they are identical based on the specified conditions
    :type vcf_list: typing.List[str]
    :param begin: The `begin` parameter in the `identical` function is used to specify a string that
    indicates the beginning of a line in the input files. If a line in the input file starts with the
    specified `begin` string, it will be skipped and not included in the comparison process. By default,
    defaults to ##
    :type begin: str (optional)
    :param line_strip: The `line_strip` parameter in the `identical` function is a boolean flag that
    determines whether each line read from the input files should be stripped of leading and trailing
    whitespaces before being compared. If `line_strip` is set to `True`, each line will be stripped
    using the `strip, defaults to True
    :type line_strip: bool (optional)
    :return: The function `identical` is returning a boolean value. It returns `True` if all the lines
    in the VCF files provided in the `vcf_list` are identical, and `False` otherwise.
    """

    vcfs_lines = []
    k = 0
    for vcf in vcf_list:
        vcfs_lines.append([])
        with open(vcf, "r") as f:
            for l in f:
                if not l.startswith(begin) or begin == "":
                    if line_strip:
                        l = l.strip()
                    vcfs_lines[k].append(l)
        k += 1

    for i in range(len(vcf_list) - 1):
        if vcfs_lines[i] != vcfs_lines[i + 1]:
            return False
    return True


def check_docker_image_exists(image_with_tag: str) -> bool:
    """
    Checks if a Docker image with a specific tag exists in the local repository.

    :param image_with_tag: Image name with tag (e.g., "image:version")
    :return: True if the image exists, False otherwise
    """
    image_name, image_tag = image_with_tag.split(":")
    try:
        # Run the `docker images` command and capture the output
        result = subprocess.run(
            ["docker", "images", "--format", "{{json . }}"],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
        # Check for errors
        if result.returncode != 0:
            log.warning(f"Error running Docker command: {result.stderr}")
            return False

        # Search for the image with the specified tag in the command output
        for image in result.stdout.split("\n"):
            if image:
                image_dict = ast.literal_eval(image)
                if (
                    image_dict.get("Repository") == image_name
                    and image_dict.get("Tag") == image_tag
                ):
                    log.debug(f"Find locally {image_with_tag}")
                    return True
        return False
    except Exception as e:
        log.warning(f"Docker image check: {e}")
        return False


def params_string_to_dict(
    params: str,
    param_sep: str = ":",
    var_val_sep: str = "=",
    val_clear: dict = {"+": ",", " ": ""},
    header: bool = True,
) -> dict:
    """
    The `params_string_to_dict` function in Python converts a string of parameters into a dictionary
    using specified separators and clears certain characters from the parameter values.

    :param params: The `params` parameter in the `params_string_to_dict` function is a string of
    parameters that you want to convert into a dictionary. It contains the information you want to parse
    and organize into key-value pairs
    :type params: str
    :param param_sep: The `param_sep` parameter in the `params_string_to_dict` function is used to
    specify the separator that separates different parameters in the input string `params`. By default,
    the `param_sep` is set to ":" in the function definition. This means that the function expects the
    parameters in the input, defaults to :
    :type param_sep: str (optional)
    :param var_val_sep: The `var_val_sep` parameter in the `params_string_to_dict` function is used to
    specify the separator between the variable and value in the input string `params`. By default, it is
    set to `"="`, which means that the function expects the format of each parameter in the `params`,
    defaults to =
    :type var_val_sep: str (optional)
    :param val_clear: The `val_clear` parameter in the `params_string_to_dict` function is a dictionary
    that contains key-value pairs used to clear specific characters from the parameter values before
    storing them in the resulting dictionary
    :type val_clear: dict
    :param header: The `header` parameter in the `params_string_to_dict` function is a boolean flag that
    determines whether the input string `params` has a header that should be skipped when processing the
    parameters. If `header` is set to `True`, the function will start processing parameters from the
    second line onwards, defaults to True
    :type header: bool (optional)
    :return: The function `params_string_to_dict` returns a dictionary containing the parameters
    extracted from the input string `params`.
    """

    # Params dict
    params_dict = {}

    # Header
    if header:
        start = 1
    else:
        start = 0

    # Params
    params_split = params.split(param_sep)
    for params_option in params_split[start:]:
        if params_option != "":
            params_option_var_val = params_option.split(var_val_sep)
            params_option_var = params_option_var_val[0].strip()
            params_option_val = params_option_var_val[1].strip()
            if params_option_val is not None:
                if params_option_val in [""]:
                    params_option_val = None
                else:
                    # clear value
                    for c in val_clear.items():
                        params_option_val = params_option_val.replace(c[0], c[1])
                params_dict[params_option_var] = params_option_val
    return params_dict


def determine_value_type(
    value: str, sep: str = ";", skip_null: list = ["", "."]
) -> str:
    """
    The function `determine_value_type` determines the type of a given value in a string format,
    handling lists of values separated by a specified separator and skipping specified null-like
    values.

    :param value: The `value` parameter in the `determine_value_type` function is the input value
    that you want to determine the type of. It can be a string containing one or more values
    separated by a specified separator (default is ';')
    :type value: str
    :param sep: The `sep` parameter in the `determine_value_type` function is used to specify the
    separator character that is used to split the input `value` string into individual values. By
    default, the separator is set to ";", but you can change it to a different character if needed,
    defaults to ;
    :type sep: str (optional)
    :param skip_null: The `skip_null` parameter in the `determine_value_type` function is a list
    that contains values that should be skipped during the type determination process. These values
    are considered as null-like or empty values and are not taken into account when determining the
    type of the given value
    :type skip_null: list
    :return: The function `determine_value_type` returns a string indicating the type of the given
    value. The possible return values are:
    - "VARCHAR" if the value contains at least one non-numeric character
    - "DOUBLE" if the value contains at least one floating-point number
    - "BIGINT" if the value contains only integers
    - None if the value is empty or does not match any
    """

    # Split the value by sep (dafult ';') to handle lists of values
    values = str(value).split(sep)

    has_float = False
    has_int = False
    has_string = False

    for val in values:
        val = val.strip()

        # Skip empty or null-like values
        if skip_null and val in skip_null:
            continue

        # Check if value is integer
        if re.match(r"^-?\d+$", val):
            has_int = True
        # Check if value is float
        elif re.match(r"^-?\d*\.\d+$", val):
            has_float = True
        else:
            has_string = True
            return "VARCHAR"  # force return VARCHAR to speed up

    if has_string:
        return "VARCHAR"
    elif has_float:
        return "DOUBLE"
    elif has_int:
        return "BIGINT"
    else:
        return None  # Default to None if no identifiable type is found


def determine_column_types(values_list: list) -> str:
    """
    The function `determine_column_types` analyzes a list of values to determine the predominant
    data type among VARCHAR, DOUBLE, and BIGINT.

    :param values_list: It seems like you have provided the code snippet for a function that
    determines the type of values in a list, but you have not provided the actual values_list that
    the function will operate on. If you provide me with the values_list, I can help you test the
    function and see how it determines the
    :type values_list: list
    :return: the type of the column based on the types of values present in the input list. It will
    return "VARCHAR" if the list contains any string values, "DOUBLE" if it contains any float
    values, "BIGINT" if it contains any integer values, and "VARCHAR" if none of the specific types
    are found.
    """

    has_float = False
    has_int = False
    has_string = False

    for value in values_list:
        value_type = determine_value_type(value)

        if value_type == "VARCHAR":
            has_string = True
            return "VARCHAR"  # force return VARCHAR to speed up
        elif value_type == "DOUBLE":
            has_float = True
        elif value_type == "BIGINT":
            has_int = True

    if has_string:
        return "VARCHAR"
    elif has_float:
        return "DOUBLE"
    elif has_int:
        return "BIGINT"
    else:
        return "VARCHAR"  # Default to VARCHAR if no identifiable type is found


def detect_column_type(column) -> str:
    """
    The function `detect_column_type` determines the type of a given column in a DataFrame as either
    DATETIME, BOOLEAN, DOUBLE, or VARCHAR.

    :param column: The function `detect_column_type` takes a column as input and determines its data
    type based on certain conditions. The conditions are as follows:
    :return: The function `detect_column_type` returns a string indicating the type of data in the input
    column. The possible return values are "DATETIME", "BOOLEAN", "DOUBLE", or "VARCHAR" based on the
    conditions checked in the function.
    """

    from pandas.api.types import is_datetime64_any_dtype as is_datetime

    if len(column) == 0:
        return "VARCHAR"
    elif is_datetime(column):
        return "DATETIME"
    elif column.dropna().apply(lambda x: str(x).lower() in ["true", "false"]).all():
        return "BOOLEAN"
    elif pd.to_numeric(column, errors="coerce").notnull().all():
        return "DOUBLE"
    else:
        try:
            pd.to_numeric(column)
            return "DOUBLE"
        except:
            return "VARCHAR"


def determine_column_number(values_list: list) -> str:
    """ """

    for value in values_list:
        if ";" in str(value):
            return "."

    return "1"


def clean_annotation_field(name: str = "", char_allowed: list = None) -> str:
    """
    The `clean_annotation_field` function removes characters from a string that are not alphanumeric or
    in a specified list.

    :param name: The `name` parameter is a string that represents the input text that you want to clean.
    It typically contains annotations or other text that you want to process
    :type name: str
    :param char_allowed: The `char_allowed` parameter is a list that contains characters that are
    allowed to remain in the `name` string after cleaning. Any character in the `name` string that is
    not alphanumeric and not in the `char_allowed` list will be removed during the cleaning process
    :type char_allowed: list
    :return: The function `clean_annotation_field` returns a cleaned version of the `name` string, where
    only alphanumeric characters and characters from the `char_allowed` list are kept.
    """

    # Init
    if char_allowed is None:
        char_allowed = []

    # Convert char_allowed to a set for faster membership testing
    char_allowed_set = set(char_allowed)

    return "".join(
        char for char in name if (char.isalnum() or char in char_allowed_set)
    )


def docker_automount() -> str:
    """
     Add needed volume to the tool container, check first if we are already inside one otherwise return empty string
    :param containerid: for other linux distribution catch container mount from container ID
    :return: string containing volume to add
    """
    if not os.path.exists("/.dockerenv"):
        log.warning("Not inside docker container, block automount option")
        return ""

    container_id = command(
        r"cat /proc/self/cgroup | grep 'docker' | sed 's/^.*\///' | tail -n1"
    )
    if not container_id:
        container_id = command("cat /proc/1/cpuset | cut -d/ -f3")
    if not container_id:
        raise ValueError("Can't get container ID from within container EXIT")

    mounts_new = ""
    mounts = json.loads(command(f"docker inspect -f json {container_id}"))[0]["Mounts"]
    for volume in mounts:
        if "sock" not in volume.get("Source") and "tmp" not in volume.get("Source"):
            mounts_new += f" -v {volume.get('Source')}:{volume.get ('Destination')}:{volume.get('Mode')}"
    return mounts_new
