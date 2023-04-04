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
from howard.tools.download import *


# Usage
# python -m pip install -e .
# howard analysis --input=my.vcf.gz --output=my.output.vcf --annotations=my.annotations.vcf.gz --stats --overview_footer
# python -m howard.main --input=my.vcf.gz --output=my.output.vcf --annotations=my.annotations.vcf.gz --stats --overview_footer


def analysis(args) -> None:

    log.info("Start")

    config = args.config

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

    # Create VCF object
    if args.input:
        vcfdata_obj = Variants(None, args.input, args.output, config, param)

        params = vcfdata_obj.get_param()

        # Quick Annotation
        if args.annotations:
            #annotation_file_list = args.annotations.split(",")
            annotation_file_list = [value for val in args.annotations for value in val.split(',')]
            log.info(f"Quick Annotation Files: {annotation_file_list}")
            param_quick_annotations = param.get("annotations",{})
            for annotation_file in annotation_file_list:
                param_quick_annotations[annotation_file] = {"INFO": None}
            params["annotations"] = param_quick_annotations

        # Quick calculations
        if args.calculations:
            #calculations_list= args.calculations.split(",")
            calculations_list= [value for val in args.calculations for value in val.split(',')]
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


def download(args) -> None:

    log.info("Start")

    databases_download(args)

    log.info("End")

# Main function
def main() -> None:
    """
    It loads a VCF file in multiple format (VCF, parquet, DB), and process, query, export data
    """

    # Main parser
    parser = argparse.ArgumentParser(
        prog="howard",
        description="""howard annotates and prioritizes genetic variations, calculates and normalizes annotations, translates vcf format and generates variants statistics""",
        #usage="howard [<shared-args>]",
        epilog="Examples:\n"
            """howard analysis --input=input.vcf.gz --output=output.tsv \n"""
            """howard analysis --input=input.vcf.gz --query="SELECT * FROM variants WHERE REF = 'A' AND POS < 100000" \n""",
        formatter_class=argparse.RawTextHelpFormatter
        )
    #"""howard download --download-url http://..."""

    parser._optionals.title = "Shared arguments"

    subparsers = parser.add_subparsers(title="Commands", dest='command')

    # Sub-command analysis
    analysis_parser = subparsers.add_parser(
        'analysis',
        description="""howard analysis command manage genetic variations to:\n- query genetic variants and annotations\n- annotates genetic variants with annotation files and tools\n- prioritizes variants with profiles (list of citeria) to calculate scores and flags\n- calculates and normalizes annotations\n- translates into various formats\n- generates variants statistics""",
        #parents=[parser],
        help='Load a VCF file in multiple format (VCF, parquet, DB), and process, query, export data',
        add_help=True,
        formatter_class=argparse.RawTextHelpFormatter
    )

    analysis_parser.add_argument(
        "--input",
        metavar="FILE",
        help="""Input file path\nFormat: BCF, VCF, TSVv, CSV, PSV, Parquet or duckDB\nFiles can be compressesd (e.g. vcf.gz, tsv.gz)""",
        required=False
    )
    analysis_parser.add_argument(
        "--output",
        metavar="FILE",
        help="""Output file path\nFormat: BCF, VCF, TSVv, CSV, PSV, Parquet or duckDB\nFiles can be compressesd (e.g. vcf.gz, tsv.gz)""",
        required=False
    )
    analysis_parser.add_argument(
        "--param",
        metavar="JSON",
        help="""Parameters file\nFormat: JSON\nDefault: {}""",
        default="{}"
    )
    analysis_parser.add_argument(
        "--query",
        metavar="QUERY",
        help="""Query in SQL format\nFormat: SQL\nExample: 'SELECT * FROM variants LIMIT 5'""",
        default=None
    )
    analysis_parser.add_argument(
        "--output_query",
        metavar="FILE",
        help="""Output Query file\nformat: VCF, TSV, Parquet...""",
        default=None
    )
    analysis_parser.add_argument(
        "--annotations",
        metavar="ANNOTATION",
        help="""Quick annotation with databases file\nFormat: list of files in Parquet, VCF, BED\nFor snpeff annotation, use keyword 'snpeff'\nFor Annovar annotation, use keyword 'annovar' with annovar code (e.g. 'annovar:refGene', 'annovar:cosmic70')""",
        nargs='+',
        default=None
    )
    analysis_parser.add_argument(
        "--calculations",
        metavar="OPERATION",
        help="Quick calculations (format: list of calculation operations) (example: NOMEN) (default: null)",
        nargs='+',
        default=None
    )
    analysis_parser.add_argument(
        "--prioritizations",
        metavar="PROFILES",
        help="Quick prioritization (format: file with profiles) (default: null)",
        default=None
    )
    analysis_parser.add_argument(
        "--overview", "--overview_header",
        help="Overview after loading data",
        action="store_true"
    )
    analysis_parser.add_argument(
        "--overview_footer",
        help="Overview before data processing",
        action="store_true"
    )
    analysis_parser.add_argument(
        "--stats", "--stats_header",
        help="Statistics after loading data",
        action="store_true"
    )
    analysis_parser.add_argument(
        "--stats_footer",
        help="Statistics before data processing",
        action="store_true"
    )

    # Sub-command download
    download_parser = subparsers.add_parser(
        'download',
        description="""Dowload databases and needed files for howar and associated tools""",
        #parents=[parser],
        help='Dowload databases and needed files for howard and associated tools',
        add_help=True,
        formatter_class=argparse.RawTextHelpFormatter
    )

    download_parser.add_argument(
        "--assembly",
        metavar="ASSEMBLY",
        help="""Assembly to download\nDefault: 'hg19'""",
        nargs='+',
        required=False,
        default=["hg19"]
    )
    download_parser.add_argument(
        "--download-annovar",
        metavar="FOLDER",
        help="Download Annovar databases within Annovar folder",
        required=False
    )
    download_parser.add_argument(
        "--download-annovar-files",
        metavar="FILE",
        help=f"""Download Annovar databases for a list of Annovar file code (see Annovar Doc)\nDefault: All available files\nExample: refGene,gnomad211_exome,cosmic70,clinvar_202*,nci60\nNote: refGene will be at leaset downloaded\nNote2: Only file that not exists or with a different size will be downloaded""",
        nargs='+',
        required=False
    )
    download_parser.add_argument(
        "--download-annovar-url",
        metavar="URL",
        help=f"""Download Annovar databases URL (see Annovar Doc)\nDefault: 'http://www.openbioinformatics.org/annovar/download/'""",
        required=False,
        default="http://www.openbioinformatics.org/annovar/download/"
    )
    download_parser.add_argument(
        "--download-snpeff",
        metavar="FOLDER",
        help="Download snpEff databases within snpEff folder",
        required=False
    )

    

    # Shared Arguments
    
    parser.add_argument(
        "--config",
        metavar="JSON",
        help="""Configuration file\nFormat: JSON\nDefault: {}""",
        default="{}"
    )
    parser.add_argument(
        "--threads",
        metavar="INTEGER",
        help="""Number of threads (replace config)""",
        default=None
    )
    parser.add_argument(
        "--verbosity",
        metavar="LEVEL",
        help="""Verbosity level (CRITICAL, ERROR, WARNING, INFO, DEBUG or NOTSET)\nDefault: INFO""",
        required=False,
        default="warning"
    )
    parser.add_argument(
        "--quiet",
        help=argparse.SUPPRESS,
        action="store_true"
    )
    parser.add_argument(
        "--verbose",
        help=argparse.SUPPRESS,
        action="store_true"
    )
    parser.add_argument(
        "--debug",
        help=argparse.SUPPRESS,
        action="store_true"
    )

    # Parse args
    args, remaining = parser.parse_known_args()

    # Verbosity
    # Default
    args.verbosity = "info"
    # Quiet
    if args.quiet:
        args.verbosity = "warning"
    # Verbose
    if args.verbose:
        args.verbosity = "info"
    # Debug
    if args.debug:
        args.verbosity = "debug"

    # Logging
    set_log_level(args.verbosity)

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

    # Change config
    args.config = config

    # Check args
    if not args.command:
        parser.print_help()
        return
    elif args.command == "analysis":
        if not (args.input or args.query):
            analysis_parser.error('At least one of --input or --query is required.')
        else:  
            analysis(args)
    elif args.command == "download":
        if not args.download_annovar and not args.download_snpeff:
            download_parser.error('At least one database to download is required.')
        else:
            download(args)



if __name__ == '__main__':
    main()


    
