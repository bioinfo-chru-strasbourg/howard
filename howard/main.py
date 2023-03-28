#!/usr/bin/env python

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
import vcf
import logging as log
import sys

from howard.objects.variants import Variants
from howard.objects.annotation import Annotation
from howard.commons import *


# Usage
# python -m pip install -e .
# howard --input=my.vcf.gz --output=my.output.vcf --annotations=my.annotations.vcf.gz --stats --overview_footer
# python -m howard.main --input=my.vcf.gz --output=my.output.vcf --annotations=my.annotations.vcf.gz --stats --overview_footer


# Main function
def main() -> None:
    """
    It loads a VCF file in multiple format (VCF, parquet, DB), and process, query, export data
    """

    
    # group.add_argument('--input', help='input file')
    # group.add_argument('--query', help='query string')

    # Args
    parser = argparse.ArgumentParser(
        description="Load a VCF file in multiple format (VCF, parquet, DB), and process, query, export data")
    parser.add_argument(
        "--input", help="Input file path (format: vcf, vcf.gz, parquet or db) Required", required=False)
    parser.add_argument(
        "--output", help="Output file path (format: vcf, vcf.gz, parquet or db)", required=False)
    parser.add_argument(
        "--config", help="Configuration file (format: JSON) (default: {})", default="{}")
    parser.add_argument(
        "--param", help="Parameters file (format: JSON) (default: {})", default="{}")
    parser.add_argument(
        "--query", help="Query (format: SQL) (default: None) (example: 'SELECT * FROM variants LIMIT 5')", default=None)
    parser.add_argument(
        "--output_query", help="Output Query file (format: VCF, TSV, Parquet...)", default=None)
    parser.add_argument(
        "--annotations", help="Quick annotation with databases file (format: list of files) (default: null)", default=None)
    parser.add_argument(
        "--calculations", help="Quick calculations (format: list of calculation operations) (example: NOMEN) (default: null)", default=None)
    parser.add_argument(
        "--prioritizations", help="Quick prioritization (format: file with profiles) (default: null)", default=None)
    parser.add_argument(
        "--threads", help="Number of threads. Will be added/replace to config file. (format: Integer) (default: null)", default=None)
    parser.add_argument("--overview", "--overview_header",
                        help="Overview after loading data", action="store_true")
    parser.add_argument(
        "--overview_footer", help="Overview before data processing", action="store_true")
    parser.add_argument("--stats", "--stats_header",
                        help="Statistics after loading data", action="store_true")
    parser.add_argument(
        "--stats_footer", help="Statistics before data processing", action="store_true")
    parser.add_argument("--verbose", help="Verbose", action="store_true")
    parser.add_argument("--debug", help="Debug", action="store_true")
    parser.add_argument(
        "--verbosity", help="Verbosity level: CRITICAL, ERROR, WARNING, INFO, DEBUG or NOTSET", required=False, default="warning")
    args = parser.parse_args()

    # Check args
    if not (args.input or args.query):
        parser.error('At least one of --input or --query is required.')

    # Verbosity
    # Verbose
    if args.verbose:
        args.verbosity = "info"
    # Debug
    if args.debug:
        args.verbosity = "debug"
    # Overview and Stats verbosity
    if args.overview or args.overview_footer or args.stats or args.stats_footer:
        args.verbosity = "info"

    # Logging
    set_log_level(args.verbosity)

    log.info("Start")

    # Load configuration in JSON format
    if os.path.exists(args.config):
        with open(args.config) as config_file:
            config = json.load(config_file)
    else:
        config = json.loads(args.config)

    # add to config
    config["verbosity"] = args.verbosity
    if args.threads:
        config["threads"] = args.threads

    # Load parameters in JSON format
    if os.path.exists(args.param):
        with open(args.param) as param_file:
            param = json.load(param_file)
    else:
        param = json.loads(args.param)


    # Connexion config
    connexion_config = {}
    if config.get("threads", None):
        connexion_config["threads"] = config.get("threads")
    if config.get("memory_limit", None):
        connexion_config["memory_limit"] = config.get("memory_limit")
    # if config.get("duckdb_compression",None):
    #     connexion_config["compression"] = config.get("duckdb_compression") # 'lz4'


    # Create VCF object
    if args.input:
        vcfdata_obj = Variants(None, args.input, args.output, config, param)

        params = vcfdata_obj.get_param()

        # Quick Annotation
        if args.annotations:
            annotation_file_list = args.annotations.split(",")
            log.info(f"Quick Annotation Files: {annotation_file_list}")
            param_quick_annotations = param.get("annotations",{})
            for annotation_file in annotation_file_list:
                param_quick_annotations[annotation_file] = {"INFO": None}
            params["annotations"] = param_quick_annotations

        # Quick calculations
        if args.calculations:
            calculations_list= args.calculations.split(",")
            log.info(f"Quick Calculations list: {calculations_list}")
            param_quick_calculations = param.get("calculation",{})
            for calculation_operation in calculations_list:
                param_quick_calculations[calculation_operation] = {}
            params["calculation"] = param_quick_calculations

        # Quick prioritization
        if args.prioritizations:
            config_profiles= args.prioritizations
            log.info(f"Quick Prioritization Config file: {config_profiles}")
            param_quick_prioritizations = param.get("prioritization",{})
            param_quick_prioritizations["config_profiles"] = config_profiles
            params["prioritization"] = param_quick_prioritizations
        
        vcfdata_obj.set_param(params)
            

        # Load data from input file
        vcfdata_obj.load_data()

        # Overview
        if args.overview:
            vcfdata_obj.get_overview()

        # Stats
        if args.stats:
            vcfdata_obj.get_stats()

        # Annotation
        if vcfdata_obj.get_param().get("annotations", None) or vcfdata_obj.get_param().get("annotation", None):
            vcfdata_obj.annotation()

        # Calculation
        if vcfdata_obj.get_param().get("calculations", None) or vcfdata_obj.get_param().get("calculation", None):
            vcfdata_obj.calculation()

        # Prioritization
        if vcfdata_obj.get_param().get("prioritizations", None) or vcfdata_obj.get_param().get("prioritization", None):
            vcfdata_obj.prioritization()

        # Output
        if vcfdata_obj.get_output():
            log.info("Exporting...")
            vcfdata_obj.export_output(export_header=True)

        # Overview footer
        if args.overview_footer:
            vcfdata_obj.get_overview()

        # Stats footer
        if args.stats_footer:
            vcfdata_obj.get_overview()

        # Query
        if args.query or param.get("query", None):
            log.info("Querying...")
            if args.query:
                query = args.query
                #result = vcfdata_obj.execute_query(args.query)
            elif param.get("query", None):
                query = param.get("query", None)
            result = vcfdata_obj.execute_query(query)
            print(result.df())

            # Output Query
            if args.output_query or param.get("output_query", None):
                log.info("Output Querying...")
                if args.output_query:
                    output_query = args.output_query
                elif param.get("output_query", None):
                    output_query = param.get("output_query", None)
                vcfdata_obj.export_output(output_file=output_query, query=query, export_header=True)

        # Close connexion
        vcfdata_obj.close_connexion()

    else:

        conn = duckdb.connect(":memory:", config=connexion_config)

        # Query
        if args.query or param.get("query", None):
            log.info("Querying...")
            if args.query:
                result_dataframe = conn.execute(args.query).df()
            elif param.get("query", None):
                result_dataframe = conn.execute(param.get("query", None)).df()
            print(result_dataframe)

            # Output Query
            if args.output_query or param.get("output_query", None):
                log.info("Output Querying (TSV only)...")
                if args.output_query:
                    output_query = args.output_query
                elif param.get("output_query", None):
                    output_query = param.get("output_query", None)
                result_dataframe.to_csv(output_query, sep='\t', index=False)

        # Close connexion
        conn.close()


    log.info("End")



if __name__ == '__main__':
    main()


    
