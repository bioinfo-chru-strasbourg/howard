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


def add_numbers(a, b):
    return a + b

def multiply_numbers(a, b):
    return a * b


def set_log_level(verbosity):
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
        raise ValueError(
            "Unknown verbosity level:" + verbosity +
            "\nPlease use any in:" + configs.keys()
        )
    log.basicConfig(
        format="#[%(asctime)s] [%(levelname)s] %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        level=configs[verbosity],
    )


def split_interval(start, end, step=None, ncuts=None):
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
        return [start + i*step for i in range(ncuts+1)] + [end]


def merge_regions(regions):
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


def create_where_clause(merged_regions, table="variants"):
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


def run_parallel_commands(commands, threads):
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
    return results


def command(command):
    """
    It runs a command in the shell and waits for it to finish
    
    :param command: The command to run
    :return: The return value is the exit status of the process.
    """
    return subprocess.Popen(command, shell=True).wait()


def run_parallel_functions(functions, threads):
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


def function_query(obj, query):
    """
    It takes a query and executes it on the object
    
    :param obj: The object that you want to query
    :param query: The query to be executed
    :return: the result of the query.
    """
    return obj.execute_query(query)
