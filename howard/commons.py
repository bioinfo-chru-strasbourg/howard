import io
import multiprocessing
import os
import re
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


def remove_if_exists(filepath: str) -> None:
    """
    If the filepath exists, remove it

    :param filepath: The path to the file you want to remove
    """
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
            genome_path = find_all(genome, '/')[0]
        except:
            log.error(f"Genome failed: no genome '{genome}'")
            raise ValueError(f"Genome failed: no genome '{genome}'")
    return genome_path


def find_nomen(hgvs:str = "", pattern = "GNOMEN:TNOMEN:ENOMEN:CNOMEN:RNOMEN:NNOMEN:PNOMEN", transcripts:list = []) -> dict:
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
            #print(one_hgvs)
            one_hgvs_split = one_hgvs.split(':')
            #print(one_hgvs_split)

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
                    if len(one_hgvs_infos_split)>1:
                        one_nomen_dict["VNOMEN"] = one_hgvs_infos_split[1]
                    # NOMEN Score 
                    if re.match(r"^NM_(.*)$", one_hgvs_infos) or re.match(r"^NM_(.*)$", one_hgvs_infos):
                        one_nomen_score += 2
                    elif re.match(r"^NR_(.*)$", one_hgvs_infos):
                        one_nomen_score += 1
                    # NOMEN with default transcript
                    if one_nomen_dict["TVNOMEN"] in transcripts or one_nomen_dict["TNOMEN"] in transcripts:
                        rank = max(get_index(one_nomen_dict["TVNOMEN"]), get_index(one_nomen_dict["TNOMEN"])) + 1
                        if rank >= 0:
                            one_nomen_score += (100 * (len(transcripts) - rank) )

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
            if nomen_dict.get(n,None):
                nomen.append(nomen_dict.get(n,None))
        nomen_dict["NOMEN"] = ":".join(nomen)

    return nomen_dict



def get_index(value, values:list = []) -> int:
    try:
        return values.index(value)
    except ValueError:
        # Si l'élément n'existe pas dans la liste, renvoyer un index négatif pour qu'il soit ignoré dans le calcul
        return -1
