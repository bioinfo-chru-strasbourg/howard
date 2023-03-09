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



# Main function
def main():
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
        "--annotation", help="Quick annotation with a database file (format: file) (default: null)", default=None)
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

        # Connexion
        if vcfdata_obj.get_input_format() in ["db", "duckdb"]:
            connexion_db = vcfdata_obj.get_input()
            vcfdata_obj.set_output(args.input)
        elif vcfdata_obj.get_output_format() in ["db", "duckdb"]:
            connexion_db = vcfdata_obj.get_output()
            # vcfdata_obj.set_output(None)
        elif vcfdata_obj.get_connexion_type() in ["memory", None]:
            connexion_db = ":memory:"
        elif vcfdata_obj.get_connexion_type() in ["tmpfile"]:
            tmp_name = tempfile.mkdtemp(prefix=vcfdata_obj.get_prefix(
            ), dir=vcfdata_obj.get_tmp_dir(), suffix=".db")
            connexion_db = f"{tmp_name}/tmp.db"
        elif vcfdata_obj.get_connexion_type() != "":
            connexion_db = vcfdata_obj.get_connexion_type()
        else:
            connexion_db = ":memory:"
        
        conn = duckdb.connect(connexion_db, config=connexion_config)

        # Quick Annotation
        if args.annotation:
            if os.path.exists(args.annotation):
                log.info(f"Quick Annotation File {args.annotation}")
                quick_annotation_file = args.annotation
                quick_annotation_name, quick_annotation_extension = os.path.splitext(
                    args.annotation)
                quick_annotation_format = quick_annotation_extension.replace(
                    ".", "")
                if quick_annotation_format in ["parquet", "duckdb"]:
                    param_quick_annotation = {
                        "annotation": {
                            "parquet": {
                                "annotations": {
                                    f"{quick_annotation_file}": {
                                        "INFO": None
                                    }
                                }
                            }
                        }
                    }
                elif quick_annotation_format in ["gz"]:
                    param_quick_annotation = {
                        "annotation": {
                            "bcftools": {
                                "annotations": {
                                    f"{quick_annotation_file}": {
                                        "INFO": None
                                    }
                                }
                            }
                        }
                    }
                else:
                    log.error(
                        f"Quick Annotation File {args.annotation} - format {quick_annotation_format} not supported yet")
                    raise ValueError(
                        f"Quick Annotation File {args.annotation} - format {quick_annotation_format} not supported yet"
                    )
                vcfdata_obj.set_param(param_quick_annotation)
            else:
                log.error(
                    f"Quick Annotation File {args.annotation} does NOT exist")
            # return
            # vcfdata_obj.get_overview()

        # Load data from input file
        vcfdata_obj.load_data()

        # Overview
        if args.overview:
            vcfdata_obj.get_overview()

        # Stats
        if args.stats:
            vcfdata_obj.get_stats()

        # Annotation
        # if param.get("annotation",None):
        if vcfdata_obj.get_param().get("annotation", None):
            vcfdata_obj.annotation()

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
                result = vcfdata_obj.execute_query(args.query)
            elif param.get("query", None):
                result = vcfdata_obj.execute_query(param.get("query", None))
            print(result.df())

    else:

        conn = duckdb.connect(":memory:", config=connexion_config)

        # Query
        if args.query or param.get("query", None):
            log.info("Querying...")
            if args.query:
                result = conn.execute(args.query)
            elif param.get("query", None):
                result = conn.execute(param.get("query", None))
            print(result.df())

    
    log.info("End")



if __name__ == '__main__':
    main()
    # my_variants = Variants('Alice')
    # my_variants.say_hello()

    # my_annotation = Annotation('Alice')
    # my_annotation.say_hello()

    # result = add_numbers(3, 4)
    # print(result)

    
