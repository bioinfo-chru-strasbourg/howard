import io
import multiprocessing
import os
from pathlib import Path
import platform
import re
import statistics
import subprocess
import sys
from tempfile import NamedTemporaryFile
import tempfile
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

import random

import pgzip
import mgzip
import bgzip

import pysam
import pysam.bcftools

import signal
from contextlib import contextmanager


file_folder = os.path.dirname(__file__)


# Main folder
folder_main = os.path.abspath(os.path.join(file_folder, ".."))
folder_config = os.path.abspath(os.path.join(folder_main, "config"))


comparison_map = {
    "gt": ">",
    "gte": ">=",
    "lt": "<",
    "lte": "<=",
    "equals": "=",
    "contains": "SIMILAR TO"
}


code_type_map = {
    "Integer": 0,
    "String": 1,
    "Float": 2,
    "Flag": 3
}


code_type_map_to_sql = {
    "Integer": "INTEGER",
    "String": "VARCHAR",
    "Float": "FLOAT",
    "Flag": "VARCHAR"
}


file_format_delimiters = {
    "vcf": "\t",
    "tsv": "\t",
    "csv": ",",
    "psv": "|",
    "bed": "\t"
}

file_format_allowed = list(file_format_delimiters.keys()) + ["json", "parquet", "duckdb"]

file_compressed_format = ["gz", "bgz"]


vcf_required_release = '##fileformat=VCFv4.2'
vcf_required_columns = ['#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO']

vcf_required = [
            vcf_required_release,
            "\t".join(vcf_required_columns)
            ]

# Tools
DEFAULT_TOOLS_FOLDER = "/tools"
DEFAULT_SNPEFF_BIN = f"{DEFAULT_TOOLS_FOLDER}/snpeff/5.1d/bin/snpEff.jar"

# URL
DEFAULT_ANNOVAR_URL = "http://www.openbioinformatics.org/annovar/download"
DEFAULT_REFSEQ_URL = "http://hgdownload.soe.ucsc.edu/goldenPath"
DEFAULT_DBNSFP_URL = "https://dbnsfp.s3.amazonaws.com"
DEFAULT_EXOMISER_URL = "http://data.monarchinitiative.org/exomiser"
DEFAULT_ALPHAMISSENSE_URL = "https://storage.googleapis.com/dm_alphamissense"
DEFAULT_DBSNP_URL = "https://ftp.ncbi.nih.gov/snp/archive"


# Databases default folder
DEFAULT_DATABASE_FOLDER = "/databases"
DEFAULT_ANNOTATIONS_FOLDER = f"{DEFAULT_DATABASE_FOLDER}/annotations/current"
DEFAULT_GENOME_FOLDER = f"{DEFAULT_DATABASE_FOLDER}/genomes/current"
DEFAULT_SNPEFF_FOLDER = f"{DEFAULT_DATABASE_FOLDER}/snpeff/current"
DEFAULT_ANNOVAR_FOLDER = f"{DEFAULT_DATABASE_FOLDER}/annovar/current"
DEFAULT_REFSEQ_FOLDER = f"{DEFAULT_DATABASE_FOLDER}/refseq/current"
DEFAULT_DBNSFP_FOLDER = f"{DEFAULT_DATABASE_FOLDER}/dbnsfp/current"
DEFAULT_EXOMISER_FOLDER = f"{DEFAULT_DATABASE_FOLDER}/exomiser/current"
DEFAULT_DBSNP_FOLDER = f"{DEFAULT_DATABASE_FOLDER}/exomiser/dbsnp"


# Deefault Assembly
DEFAULT_ASSEMBLY = "hg19"

# DuckDB extension
DUCKDB_EXTENSION = f"{file_folder}/duckdb_extension"

# Variables
MACHIN_LIST = {
    "amd64": "amd64",
    "arm64": "arm64"
}

# bcftools format allowed 
BCFTOOLS_FORMAT = ["vcf", "bed"]

LOG_FORMAT = "#[%(asctime)s] [%(levelname)s] %(message)s"

CODE_TYPE_MAP = {
            "Integer": 0,
            "String": 1,
            "Float": 2,
            "Flag": 3
        }

GENOTYPE_MAP = {
            None: ".",
            -1: "A",
            -2: "G",
            -3: "R"
        }

DTYPE_LIMIT_AUTO = 10000

DEFAULT_CHUNK_SIZE = 1024*1024

def remove_if_exists(filepaths: list) -> None:
    """
    The function removes a file if it exists at the specified filepath(s).
    
    :param filepaths: A list of file paths that you want to check for existence and remove if they exist
    :type filepaths: list
    """
    if type(filepaths) is str:
        filepaths = [filepaths]
    for filepath in filepaths:
        if os.path.exists(filepath):
            if os.path.isdir(filepath):
                #os.rmdir(filepath)
                shutil.rmtree(filepath)
            else:
                os.remove(filepath)


def set_log_level(verbosity: str, log_file:str = None) -> str:
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
        encoding='utf-8', 
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
        return [start + i*step for i in range(ncuts+1)]


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
                current_region = (current_region[0], current_region[1], max(
                    current_region[2], region[2]))
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
        where_clause_chrom[chrom] += f" {where_clause_chrom_sep} ({table}.POS >= {start} AND {table}.POS <= {stop}) "

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
    return output.decode('utf-8').strip()


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
        results.append(pool.apply_async(command, args=(
            cmd,), error_callback=lambda e: print(e)))
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
        results.append(pool.apply_async(func, args=(1, "hello"),
                       error_callback=lambda e: print(e)))
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
    
    result = []
    for root, dirs, files in os.walk(path):
        for filename in fnmatch.filter(files, name):
            if os.path.exists(os.path.join(root, filename)):
                result.append(os.path.join(root, filename))
        # if name in files:
        #     result.append(os.path.join(root, name))
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
        elif assembly and find_all(assembly+".fa", genome_dir):
            genome_path = find_all(assembly+".fa", genome_dir)[0]
    return genome_path


def find_file_prefix(input_file:str = None, prefix:str = None, folder:str = None, assembly:str = None) -> str:
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

    output_file = None
    if os.path.exists(input_file):
        output_file = input_file
    else:
        #refgene_prefix = "refGene"
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


def find_nomen(hgvs: str = "", pattern="GNOMEN:TNOMEN:ENOMEN:CNOMEN:RNOMEN:NNOMEN:PNOMEN", transcripts: list = []) -> dict:
    """
    > This function takes a HGVS string and a list of transcripts and returns a dictionary with the best
    NOMEN for each HGVS string

    :param hgvs: The HGVS string to parse
    :type hgvs: str
    :param pattern: This is the pattern that you want to use to construct the NOMEN. The default is
    "GNOMEN:TNOMEN:ENOMEN:CNOMEN:RNOMEN:NNOMEN:PNOMEN". This means that the NOMEN will be constructed by
    joining, defaults to GNOMEN:TNOMEN:ENOMEN:CNOMEN:RNOMEN:NNOMEN:PNOMEN (optional)
    :param transcripts: list of transcripts to use for ranking
    :type transcripts: list
    :return: A dictionary with the following keys:
        NOMEN
        CNOMEN
        RNOMEN
        NNOMEN
        PNOMEN
        TVNOMEN
        TNOMEN
        VNOMEN
        ENOMEN
        GNOMEN
    """

    empty_nomen_dict = {
        "NOMEN": None,
        "CNOMEN": None,
        "RNOMEN": None,
        "NNOMEN": None,
        "PNOMEN": None,
        "TVNOMEN": None,
        "TNOMEN": None,
        "VNOMEN": None,
        "ENOMEN": None,
        "GNOMEN": None,
    }

    nomen_dict = empty_nomen_dict.copy()

    if hgvs != "nan":

        hgvs_split = str(hgvs).split(',')

        nomen_score_max = 0

        for one_hgvs in hgvs_split:
            # print(one_hgvs)
            one_hgvs_split = one_hgvs.split(':')
            # print(one_hgvs_split)

            one_nomen_score = 0
            one_nomen_dict = empty_nomen_dict.copy()

            for one_hgvs_infos in one_hgvs_split:

                if re.match(r"^[NX][MRP]_(.*)$", one_hgvs_infos):
                    # Transcript with version
                    one_nomen_dict["TVNOMEN"] = one_hgvs_infos
                    one_nomen_score += 1
                    # Split transcript
                    one_hgvs_infos_split = one_hgvs_infos.split('.')
                    # Transcript
                    one_nomen_dict["TNOMEN"] = one_hgvs_infos_split[0]
                    # Transcript version
                    if len(one_hgvs_infos_split) > 1:
                        one_nomen_dict["VNOMEN"] = one_hgvs_infos_split[1]
                    # NOMEN Score
                    if re.match(r"^NM_(.*)$", one_hgvs_infos) or re.match(r"^NM_(.*)$", one_hgvs_infos):
                        one_nomen_score += 2
                    elif re.match(r"^NR_(.*)$", one_hgvs_infos):
                        one_nomen_score += 1
                    # NOMEN with default transcript
                    if one_nomen_dict["TVNOMEN"] in transcripts or one_nomen_dict["TNOMEN"] in transcripts:
                        rank = max( get_index(one_nomen_dict["TVNOMEN"], transcripts), get_index(one_nomen_dict["TNOMEN"], transcripts))
                        if rank >= 0:
                            one_nomen_score += (100 *
                                                (len(transcripts) - rank))

                elif re.match(r"^c\.(.*)$", one_hgvs_infos) or re.match(r"^g\.(.*)$", one_hgvs_infos) or re.match(r"^m\.(.*)$", one_hgvs_infos):
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
        for n in pattern.split(':'):
            if nomen_dict.get(n, None):
                nomen.append(nomen_dict.get(n, None))
        nomen_dict["NOMEN"] = ":".join(nomen)

    return nomen_dict


def extract_snpeff_hgvs(snpeff:str = "", header:str = ['Allele', 'Annotation', 'Annotation_Impact', 'Gene_Name', 'Gene_ID', 'Feature_Type', 'Feature_ID', 'Transcript_BioType', 'Rank', 'HGVS.c', 'HGVS.p', 'cDNA.pos / cDNA.length', 'CDS.pos / CDS.length', 'AA.pos / AA.length', 'Distance', 'ERRORS / WARNINGS / INFO']) -> str:
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

    snpeff_hgvs = ""

    if snpeff != "nan":
        snpeff_hgvs = "snpeff_hgvs_list"

    # Split snpeff ann values
    snpeff_infos = [x.split('|') for x in snpeff.split(",")]

    # Create Dataframe
    snpeff_dict = {}
    for i in range(len(header)):
        snpeff_dict[header[i]] = [x[i] for x in snpeff_infos]
    df = pd.DataFrame.from_dict(snpeff_dict, orient='index').transpose()

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
        # Si l'élément n'existe pas dans la liste, renvoyer un index négatif pour qu'il soit ignoré dans le calcul
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
            filename_name_name, filename_name_extension = os.path.splitext(
                filename_name)
            filename_format = filename_name_extension.replace(".", "")
    else:
        filename_format = "unknown"
    return filename_format


def findbypipeline(df, samples:list = []):
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
            if len(sample_infos)>i:
                sample_dict[format_fields[i]] = sample_infos[i]
        
        # Check if GT not null
        if sample_dict["GT"].replace("0",".") not in ['','.','./.','.|.']:
            nb_pipeline_find += 1

    return f"{nb_pipeline_find}/{nb_pipeline}"


def genotypeconcordance(df, samples:list = []):
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
                if len(sample_infos)>i:
                    sample_dict[format_fields[i]] = sample_infos[i]
            
            # Check if GT not null
            #genotype_list[sample_dict["GT"]] = 1
            if sample_dict["GT"] not in ['','.','./.','.|.']:
                genotype_list[sample_dict["GT"]] = 1

    return str(len(genotype_list) == 1).upper()


def genotype_compression(genotype:str = "") -> str:
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

    genotype_compressed = genotype
    genotype_compressed = re.sub(r'\.', '0', genotype_compressed)
    genotype_compressed = re.sub(r'\D', '', genotype_compressed)
    genotype_compressed = ''.join(sorted(set(genotype_compressed)))

    return genotype_compressed


def genotype_barcode(genotype:str = "") -> str:
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


def barcode(df, samples:list = []):
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
                if len(sample_infos)>i:
                    sample_dict[format_fields[i]] = sample_infos[i]
            
            # generate barcode
            barcode.append(genotype_barcode(sample_dict["GT"]))

    return "".join(barcode)


def trio(df, samples:list = []):
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
        ("011","101","111","021","201","121","211"): "dominant",
        ("112","212","122","222"): "recessive"
    }

    trio_variant_type = "unknown"
    for case in switcher:
        if trio_barcode in case:
            trio_variant_type = switcher[case]
            break

    return trio_variant_type


def vaf_normalization(row, sample:str) -> str:
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
        if len(sample_genotype_infos)>i:
            sample_genotype_dict[format_fields[i]] = sample_genotype_infos[i]
        else:
            sample_genotype_dict[format_fields[i]] = "."

    # Find VAF
    # Default VAF
    vaf = "."
    # VAF from FREQ
    if "FREQ" in sample_genotype_dict:
        if sample_genotype_dict["FREQ"] != ".":
            vaf_freq = sum(map(float, sample_genotype_dict["FREQ"].replace('%', '').split(',')))
            if vaf_freq:
                vaf = round(vaf_freq / 100, 4)
    # VAF from DP4
    elif "DP4" in sample_genotype_dict:
        if sample_genotype_dict["DP4"] != ".":
            dp4_split = sample_genotype_dict["DP4"].split(",")
            if dp4_split != ['.']:
                dp4_dp = sum(map(int, dp4_split))
                pd4_alt = sum(map(int, dp4_split[2:]))
                vaf = round(pd4_alt / dp4_dp, 6) if dp4_dp else "."
    # VAF from AD
    elif "AD" in sample_genotype_dict:
        if sample_genotype_dict["AD"] != ".":
            ad_split = sample_genotype_dict["AD"].split(",")
            if ad_split != ['.']:
                ad_dp = sum(map(int, ad_split))
                ad_alt = sum(map(int, ad_split[1:]))
                vaf = round(ad_alt / ad_dp, 6) if ad_dp else "."

    # add vaf into genotype
    sample_genotype_dict["VAF"] = vaf

    # Create new genotype info
    genotype_with_vaf = ":".join([str(v) for v in sample_genotype_dict.values()])

    return genotype_with_vaf


def genotype_stats(df, samples:list = [], info:str = "VAF"):
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
        info+"_stats_nb": 0,
        info+"_stats_list": None,
        info+"_stats_min": None,
        info+"_stats_max": None,
        info+"_stats_mean": None,
        info+"_stats_mediane": None,
        info+"_stats_stdev": None
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
            if len(sample_infos)>i:
                sample_dict[format_fields[i]] = sample_infos[i]
        
        # Check if GT not null
        if info in sample_dict:
            try:
                vaf_float = float(sample_dict[info])
            except:
                vaf_float = None
            if vaf_float:
                vaf_list.append(vaf_float)

    vaf_stats[info+"_stats_nb"] = len(vaf_list)
    if vaf_list:
        vaf_stats[info+"_stats_list"] = ":".join([str(x) for x in vaf_list])

    # Compute min, max, mean and median only if the list is not empty
    if vaf_list:
        vaf_stats[info+"_stats_min"] = min(vaf_list)
        vaf_stats[info+"_stats_max"] = max(vaf_list)
        vaf_stats[info+"_stats_mean"] = statistics.mean(vaf_list)
        vaf_stats[info+"_stats_mediane"] = statistics.median(vaf_list)

    # Check if there are at least 2 values in the list before computing variance or stdev
    if len(vaf_list) >= 2:
        vaf_stats[info+"_stats_stdev"] = statistics.stdev(vaf_list)

    return vaf_stats


def extract_file(file_path:str, path:str = None, threads:int = 1):
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

    if file_path.endswith('.zip'):
        if not path:
            path = os.path.dirname(file_path)
        with zipfile.ZipFile(file_path, 'r') as zip_ref:
            zip_ref.extractall(path)
    elif file_path.endswith('.gz'):
        concat_and_compress_files(input_files=[file_path], output_file=file_path[:-3], compression_type="none", threads=threads)


def download_file(url:str, dest_file_path:str, chunk_size:int = 1024*1024, try_aria:bool = True, threads:int = 1):
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
    :param threads: The `threads` parameter specifies the number of threads to be used for downloading
    the file. It determines the number of simultaneous connections that will be made to download the
    file. By default, it is set to 1, which means that only one connection will be made at a time.
    Increasing the value, defaults to 1
    :type threads: int (optional)
    :return: a boolean value indicating whether the file was successfully downloaded and saved to the
    specified destination file path.
    """

    # Create folder if not exists
    if not os.path.exists(os.path.dirname(dest_file_path)):
        Path(os.path.dirname(dest_file_path)).mkdir(parents=True, exist_ok=True)

    try:

        if try_aria:

            # Aria command
            aria_command = f"aria2c -c -s {threads} -x {threads} -j 1 {url} -d {os.path.dirname(dest_file_path)} -o {os.path.basename(dest_file_path)}"

            # Launch command
            output = os.system(aria_command)

        else:

            assert False

    except:

        # Request
        with requests.get(url, stream=True) as r:
            r.raise_for_status()

            # Create a temporary file
            tmp_file_path = dest_file_path + '.tmp'

            # Open the temporary file for writing in binary mode
            with open(tmp_file_path, 'wb') as f:

                # Download the file by chunks
                for chunk in r.iter_content(chunk_size=chunk_size):

                    # Write the chunk to the temporary file
                    f.write(chunk)

        # Move the temporary file to the final destination
        shutil.move(tmp_file_path, dest_file_path)

    return os.path.exists(dest_file_path)


def get_bin(bin:str = None, tool:str = None, bin_type:str = "jar", config:dict = {}, default_folder:str = DEFAULT_TOOLS_FOLDER):
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
    search for in the config dict (e.g., `jar`, `sh`). In this case, the default value is "jar", which indicates that the function is searching
    for a JAR file, defaults to jar
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

    # Config - snpEff
    bin_file = config.get("tools", {}).get(
        tool, {}).get(bin_type, None)

    # Config - check tools
    if not bin_file or not os.path.exists(bin_file):
        
        # Try to find bin file
        try:
            bin_file = find_all(bin, default_folder)[0]
        except:
            return None

        # Check if found
        if not os.path.exists(bin_file):
            return None

    return bin_file


def get_snpeff_bin(config:dict = {}):
    """
    This function retrieves the path to the snpEff.jar file from a configuration dictionary or searches
    for it in the file system if it is not specified in the configuration.
    
    :param config: A dictionary containing configuration information for the snpEff tool, including the
    path to the snpEff jar file. If no configuration is provided, an empty dictionary is used
    :type config: dict
    :return: the path to the snpEff.jar file, either from the configuration dictionary or by searching
    for it in the file system. If the file is not found, it returns None.
    """

    # Config - snpEff
    snpeff_jar = config.get("tools", {}).get(
        "snpeff", {}).get("jar", None)

    # Config - check tools
    if not snpeff_jar or not os.path.exists(snpeff_jar):
        # Try to find snpEff.jar
        try:
            snpeff_jar = find_all('snpEff.jar', '/')[0]
        except:
            return None
        # Check if found
        if not os.path.exists(snpeff_jar):
            return None
            
    return snpeff_jar



def concat_file(input_files:list, output_file:str) -> bool:
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

    with open(output_file, 'w') as outfile:
        for input_file in input_files:
            with open(input_file) as infile:
                for line in infile:
                    outfile.write(line)

    return os.path.exists(output_file)


def compress_file(input_file:str, output_file:str) -> bool:
    """
    This function compresses a file using the BGZF compression algorithm.
    
    :param input_file: The path and name of the input file that needs to be compressed
    :type input_file: str
    :param output_file: The output_file parameter is a string that represents the name and path of
    the file where the compressed data will be written
    :type output_file: str
    """

    concat_and_compress_files(input_files=[input_file], output_file=output_file)

    return os.path.exists(output_file)


def get_compression_type(filepath:str) -> str:
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
        with open(filepath, 'rb') as test_f:
            # Firsts bits
            bit_1 = test_f.read(1)
            bit_2 = test_f.read(1)
            bit_3 = test_f.read(1)
            bit_4 = test_f.read(1)
            # If bit is compress
            if bit_1 == b'\x1f' and bit_2 == b'\x8b' and bit_3 == b'\x08':
                # If bit 4 == x04 is compress 'bgzip' type
                if bit_4 == b'\x04':
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



def concat_into_infile(input_files:list, compressed_file:object, compression_type:str = "none", threads:int = 1, block:int = 10 ** 6) -> bool:
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
    if compression_type in ['none']:
        open_type = "t"
    else:
        open_type = "b"

    # Open input files
    for input_file in input_files:

        # Input file compression type
        input_compression_type = get_compression_type(input_file)
        if input_compression_type in ['bgzip']:
            with open(input_file, 'rb') as raw:
                with bgzip.BGZipReader(raw, num_threads=threads, raw_read_chunk_size=block) as infile:
                    if open_type in ['t']:
                        compressed_file.write(str(infile.read(), 'utf-8'))
                    else:
                        shutil.copyfileobj(infile, compressed_file)
        elif input_compression_type in ['gzip']:
            # See https://pypi.org/project/mgzip/
            with mgzip.open(input_file, 'r'+open_type, thread=threads) as infile:
                shutil.copyfileobj(infile, compressed_file)
        elif input_compression_type in ['none']:
            with open(input_file, 'r'+open_type) as infile:
                shutil.copyfileobj(infile, compressed_file)
        else:
            raise ValueError(f"Input file compression type unknown: {input_file}")

    return True


def concat_and_compress_files(input_files: list, output_file: str, compression_type:str = "bgzip", threads:int = 1, block:int = 10 ** 6, compression_level:int = 6, sort:bool = False, index:bool = False) -> bool:
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
    if compression_type not in ['bgzip', 'gzip', 'none']:
        compression_type = 'bgzip'

    # Output file compressions type
    if compression_type in ['none']:
        open_type = "t"
    else:
        open_type = "b"

    output_file_tmp = output_file+"."+str(random.randrange(1000))+".tmp"

    if compression_type in ['gzip']:
        # See https://pypi.org/project/mgzip/
        with pgzip.open(output_file_tmp, 'w'+open_type, thread=threads, blocksize=threads * block, compresslevel=compression_level) as compressed_file:
            concat_into_infile(input_files, compressed_file, compression_type=compression_type, threads=threads, block=block)
    elif compression_type in ['bgzip']:
        with open(output_file_tmp, 'w'+open_type) as compressed_file_raw:
            with bgzip.BGZipWriter(compressed_file_raw, num_threads=threads) as compressed_file:
                concat_into_infile(input_files, compressed_file, compression_type=compression_type, threads=threads, block=block)
    elif compression_type in ['none']:
        with open(output_file_tmp, 'w'+open_type) as compressed_file:
            concat_into_infile(input_files, compressed_file, compression_type=compression_type, threads=threads, block=block)     

    # Output file
    if sort:
        # Sort with pysam
        try:
            pysam.bcftools.sort(f"-Oz{compression_level}", "-o", output_file, "-T", output_file_tmp, output_file_tmp, threads=threads, catch_stdout=False)
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


def get_plateform_name_from_duckdb(conn:duckdb.DuckDBPyConnection) -> str:
    """
    The function `get_plateform_name_from_duckdb` returns the platform information from a DuckDB connection.
    
    :param conn: The `conn` parameter is an instance of the `DuckDBPyConnection` class from the `duckdb`
    module. It represents a connection to a DuckDB database
    :type conn: duckdb.DuckDBPyConnection
    :return: the platform information from the DuckDB connection.
    """

    return conn.query(f"PRAGMA platform").df()["platform"][0]    


def get_duckdb_extension_file(extension_name:str, conn:duckdb.DuckDBPyConnection, download:bool = True) -> str:
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
    extension_file_path = f'v{release_version_number}/{platform_name}/{extension_name}.duckdb_extension'
    extension_file = f'{DUCKDB_EXTENSION}/{extension_file_path}'
    extension_file_gz = f'{DUCKDB_EXTENSION}/{extension_file_path}.gz'
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
            with mgzip.open(extension_file_gz, 'rb', thread=1) as infile:
                with open(extension_file, 'wb') as compressed_file:
                    shutil.copyfileobj(infile, compressed_file)

    if os.path.exists(extension_file):
        return extension_file
    else:
        return None


def load_duckdb_extension(conn:duckdb.DuckDBPyConnection, duckdb_extensions:list) -> bool:
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
        duckdb_extension_file = get_duckdb_extension_file(extension_name, conn=conn, download=True)
        try:
            if duckdb_extension_file and os.path.exists(duckdb_extension_file):
                # Try loading extension by file
                conn.query(f"INSTALL '{duckdb_extension_file}'; LOAD {extension_name}; ")
            else:
                loaded = False
        except:
            try:
                with time_limit(60):
                    # Try loading extension by name
                    conn.query(f"INSTALL '{extension_name}'; LOAD {extension_name}; ")
            except:
                loaded = False

    return loaded


class TimeoutException(Exception): pass

def time_limit(seconds):
    def signal_handler(signum, frame):
        raise TimeoutException("Timed out!")
    signal.signal(signal.SIGALRM, signal_handler)
    signal.alarm(seconds)
    try:
        yield
    finally:
        signal.alarm(0)


def duckdb_execute(query:str, threads:int = 1) -> bool:
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

    conn = duckdb.connect(config={"threads":threads})
    conn.execute("SET max_expression_depth TO 10000") 
    if conn.execute(query):
        conn.close()
        return True
    else:
        conn.close()
        return False


def genome_build_switch(assembly:str) -> str:
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
