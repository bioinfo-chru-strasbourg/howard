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

    config = args.config

    param = {}

    vcfdata_obj = Variants(None, args.input, args.output, config, param)

    params = vcfdata_obj.get_param()

    # Explode Infos
    if args.explode_infos:
        params["explode_infos"] = args.explode_infos
        params["explode_infos_prefix"] = args.explode_infos_prefix
        params["explode_infos_fields"] = args.explode_infos_fields
    else:
        config["access"] = "RO"

    vcfdata_obj.set_param(params)
    vcfdata_obj.set_config(config)

    # Load
    if args.input:
        vcfdata_obj.load_data()

    # include_header
    if "include_header" in args and args.include_header:
        params["header_in_output"] = args.include_header

    # query_limit
    query_limit=10
    if "query_limit" in args and args.query_limit:
        params["query_limit"] = int(args.query_limit)
        query_limit = int(args.query_limit)

    # query_print_mode
    query_print_mode=None
    if "query_print_mode" in args and args.query_print_mode:
        query_print_mode = args.query_print_mode

    # Explode infos
    if args.explode_infos:
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
            

    log.info("End")
