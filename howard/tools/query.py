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
from tabulate import tabulate

from howard.objects.variants import Variants
from howard.objects.database import Database
from howard.commons import *
from howard.tools.databases import *



def query(args:argparse) -> None:
    """
    This Python function loads and queries data from a VCF file based on user input and exports the
    results.
    
    :param args: args is an object that contains the arguments passed to the function. It is likely a
    Namespace object created by parsing command line arguments using argparse
    :type args: argparse
    """

    log.info("Start")

    # Config infos
    if "arguments_dict" in args:
        arguments_dict = args.arguments_dict
    else:
        arguments_dict = None
    if "setup_cfg" in args:
        setup_cfg = args.setup_cfg
    else:
        setup_cfg = None
    config = args.config

    # Load parameters in JSON format
    param = {}
    if "param" in args:
        if isinstance(args.param, str) and os.path.exists(full_path(args.param)):
            with open(full_path(args.param)) as param_file:
                param = json.load(param_file)
        else:
            param = json.loads(args.param)

    vcfdata_obj = Variants(None, args.input, args.output, config, param)

    param = vcfdata_obj.get_param()

    # Explode Infos
    if args.explode_infos and "explode_infos" not in param:
        param["explode_infos"] = args.explode_infos
    if args.explode_infos_prefix and "explode_infos_prefix" not in param:
        param["explode_infos_prefix"] = args.explode_infos_prefix
    if args.explode_infos_fields and "explode_infos_fields" not in param:
        param["explode_infos_fields"] = args.explode_infos_fields

    if not param.get("explode_infos",False):
        config["access"] = "RO"

    vcfdata_obj.set_param(param)
    vcfdata_obj.set_config(config)

    # Load
    if args.input:
        vcfdata_obj.load_data()

    # include_header
    if "include_header" in args and "include_header" not in param:
        param["header_in_output"] = args.include_header

    # query_limit
    query_limit=10
    if "query_limit" in args and "query_limit" not in param:
        param["query_limit"] = int(args.query_limit)
        query_limit = int(args.query_limit)

    # query_print_mode
    query_print_mode=None
    if "query_print_mode" in args and "query_print_mode" not in param:
        query_print_mode = args.query_print_mode

    # Explode infos
    if param.get("explode_infos",False):
        vcfdata_obj.explode_infos()

    # Query
    if args.query or param.get("query", None):
        query = ""
        if args.query:
            query = args.query
        elif param.get("query", None):
            query = param.get("query")
        if args.output:
            log.info("Exporting Querying...")
            vcfdata_obj.export_output(query=query, export_header=True)
        else:
            log.info("Querying...")
            if query_print_mode in ["markdown"]:
                print(vcfdata_obj.get_query_to_df(query, limit=query_limit).to_markdown())
            elif query_print_mode in ["tabulate"]:
                print(tabulate(vcfdata_obj.get_query_to_df(query, limit=query_limit), headers='keys', tablefmt='psql'))
            else:
                print(vcfdata_obj.get_query_to_df(query, limit=query_limit))
    else:
        # Parser
        parser = help_generation(arguments_dict=arguments_dict, setup=setup_cfg, output_type="parser")
        parser.print_help()
        print("")
        log.error(f"No query provided")
        raise ValueError(f"No query provided")
            

    log.info("End")
