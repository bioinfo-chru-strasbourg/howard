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
from howard.functions.commons import *
from howard.functions.databases import *



def process(args:argparse) -> None:
    """
    The "process" function processes input arguments, loads parameters in JSON format, creates a VCF
    object, performs quick annotations, calculations, prioritizations, and queries, exports output, and
    closes the connection.
    
    :param args: args is a variable that contains the arguments passed to the function "process". It is
    assumed to be an object with several attributes, including "config", "param", "input", "output",
    "annotations", "calculations", "prioritizations", and "query". These attributes are used to
    :type args: argparse
    """

    log.info("Start")

    # Load config args
    arguments_dict, setup_cfg, config, param = load_config_args(args)

    # Create variants object
    vcfdata_obj = Variants(input=args.input, output=args.output, config=config, param=param)

    # Get Config and Params
    config = vcfdata_obj.get_config()
    param = vcfdata_obj.get_param()

    # Load args into param
    param = load_args(param=param, args=args, arguments_dict=arguments_dict, command="process", strict=False)

    # Re-Load Config and Params
    vcfdata_obj.set_param(param)
    vcfdata_obj.set_config(config)

    # Load data
    vcfdata_obj.load_data()

    # Annotation
    vcfdata_obj.annotation_hgvs()
    vcfdata_obj.annotation()
    vcfdata_obj.calculation()
    vcfdata_obj.prioritization()
    
    # Query
    if param.get("query", {}).get("query", None):

        log.info("Querying...")

        # Parameters
        query = param.get("query", {}).get("query", None)
        query_limit = param.get("query", {}).get("query_limit", None)
        query_print_mode = param.get("query", {}).get("query_print_mode", None)

        # Print query
        if query_print_mode in ["markdown"]:
            print(vcfdata_obj.get_query_to_df(query, limit=query_limit).to_markdown())
        elif query_print_mode in ["tabulate"]:
            print(tabulate(vcfdata_obj.get_query_to_df(query, limit=query_limit), headers='keys', tablefmt='psql'))
        else:
            print(vcfdata_obj.get_query_to_df(query, limit=query_limit))

    # Export
    vcfdata_obj.export_output(query=param.get("query", {}).get("query", None))
    #vcfdata_obj.export_output()

    # Close connexion
    vcfdata_obj.close_connexion()

    log.info("End")