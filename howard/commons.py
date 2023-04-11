import io
import multiprocessing
import os
import re
import statistics
import subprocess
from tempfile import NamedTemporaryFile
import tempfile
import duckdb
import json
import argparse
import Bio.bgzf as bgzf
import pandas as pd
# import pyarrow as pa
# import pyarrow.parquet as pq
import vcf
import logging as log
import shutil
import urllib.request
import zipfile
import gzip
import requests


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
    "bcf": "\t",
    "tsv": "\t",
    "csv": ",",
    "psv": "|"
}


vcf_required_release = '##fileformat=VCFv4.2'
vcf_required_columns = ['#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO']

vcf_required = [
            vcf_required_release,
            "\t".join(vcf_required_columns)
            ]

default_snpeff_bin = "/tools/snpeff/5.1d/bin/snpEff.jar"

default_annovar_url = "http://www.openbioinformatics.org/annovar/download/"


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
            os.remove(filepath)


def set_log_level(verbosity: str) -> str:
    """
    It sets the log level of the Python logging module

    :param verbosity: The level of verbosity
    """
    configs = {
        "debug": log.DEBUG,
        "info": log.INFO,
        "warning": log.WARNING,
        "error": log.ERROR,
        "critical": log.CRITICAL,
    }
    if verbosity not in configs.keys():
        raise ValueError("Unknown verbosity level:" + verbosity)
    log.basicConfig(
        format="#[%(asctime)s] [%(levelname)s] %(message)s",
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


def run_parallel_functions_new(functions: list, threads:int = 1) -> None:
    """
    Runs a list of functions in parallel using multiprocessing.Pool.

    :param functions: A list of functions to run.
    :param threads: Number of threads to use.
    """
    with multiprocessing.Pool(threads) as p:
        p.map(lambda f: f(), functions)


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
        if name in files:
            result.append(os.path.join(root, name))
    return result


def get_bgzip(threads: int = 1, level: int = 1):
    """
    It checks if the bgzip command supports the --threads option, and if it does, it adds it to the
    command

    :param threads: number of threads to use for compression, defaults to 1 (optional)
    :param level: Compression level, defaults to 1 (optional)
    :return: The command to use for bgzip or gzip.
    """
    command_gzip = f" bgzip -c "
    # Check threads in bgzip command (error in macos)
    result_command_bgzip = subprocess.run(
        "bgzip --help 2>&1 | grep 'threads'", shell=True, stdout=subprocess.PIPE)
    if not result_command_bgzip.returncode:
        command_gzip += f" --threads={threads} --compress-level={level} "
    else:
        command_gzip = f" gzip -c -{level} "
    return command_gzip


def find_genome(genome_path: str, genome: str = "hg19.fa"):
    """
    It checks if the genome file exists, and if not, it tries to find it

    :param genome_path: the path to the genome file
    :param genome: the name of the genome file, defaults to hg19.fa (optional)
    :return: The path to the genome file.
    """
    # check genome
    if not os.path.exists(genome_path):
        log.warning(f"Genome warning: no genome '{genome}'. Try to find...")
        # Try to find genome
        try:
            genome_path = find_all(genome, '/databases')[0]
        except:
            log.error(f"Genome failed: no genome '{genome}'")
            raise ValueError(f"Genome failed: no genome '{genome}'")
    return genome_path


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
        if filename_format in ["gz"]:
            filename_name_name, filename_name_extension = os.path.splitext(
                filename_name)
            filename_format = filename_name_extension.replace(".", "")
    else:
        filename_format = "unknown"
    return filename_format


def get_file_compressed(filename: str = None) -> bool:
    """
    This function takes a filename as input and returns True if the file is compressed (in bgzip) and False if it
    is not

    :param filename: the name of the file to be checked
    :type filename: str
    :return: A boolean value.
    """
    compress_format = None
    if filename:
        filename_name, filename_extension = os.path.splitext(filename)
        compress_format = filename_extension.replace(".", "")
    return compress_format in ["gz", "bcf"]


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


def extract_file(file_path:str):
    """
    The function extracts a compressed file in .zip or .gz format based on the file path provided.
    
    :param file_path: The file path parameter is a string that represents the path to a file that needs
    to be extracted. The function checks if the file has a ".zip" or ".gz" extension and extracts it
    accordingly
    :type file_path: str
    """
    """
    The function extracts a compressed file if it is in .zip or .gz format.
    
    :param file_path: The file path parameter is a string that represents the path to a file that needs
    to be extracted. The function checks if the file has a ".zip" or ".gz" extension and extracts it
    accordingly
    :type file_path: str
    """
    if file_path.endswith('.zip'):
        with zipfile.ZipFile(file_path, 'r') as zip_ref:
            zip_ref.extractall(os.path.dirname(file_path))
    elif file_path.endswith('.gz'):
        with gzip.open(file_path, 'rb') as f_in:
            with open(file_path[:-3], 'wb') as f_out:
                f_out.write(f_in.read())


def download_file(url:str, dest_file_path:str, chunk_size:int = 1024*1024):
    """
    This function downloads a file from a given URL and saves it to a specified destination file path in
    chunks.
    
    :param url: The URL of the file to be downloaded
    :param dest_file_path: The path where the downloaded file will be saved
    :param chunk_size: The size of each chunk of data to be downloaded at a time. The default value is 1
    MB
    """
    # Create a temporary file
    tmp_file_path = dest_file_path + '.tmp'
    with requests.get(url, stream=True) as r:
        r.raise_for_status()
        # Open the temporary file for writing in binary mode
        with open(tmp_file_path, 'wb') as f:
            # Download the file by chunks
            for chunk in r.iter_content(chunk_size=chunk_size):
                # Write the chunk to the temporary file
                f.write(chunk)
    # Move the temporary file to the final destination
    shutil.move(tmp_file_path, dest_file_path)


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
