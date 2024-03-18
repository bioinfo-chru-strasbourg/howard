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
from howard.objects.database import Database
from howard.functions.commons import *
from howard.functions.databases import *



def stats(args:argparse) -> None:
    """
    The stats() function takes in arguments, loads data from an input file, gets statistics on the data,
    and closes the connection.
    
    :param args: args is a parameter that is passed to the function stats(). It is likely an object or a
    dictionary that contains various arguments or parameters that are needed by the function to perform
    its tasks. Some of the arguments that may be included in args are input file path, configuration
    settings, and other parameters that are
    :type args: argparse
    """

    log.info("Start")

    # Load config args
    arguments_dict, setup_cfg, config, param = load_config_args(args)

    # Create variants object
    vcfdata_obj = Variants(input=args.input, config=config, param=param)

    # Get Config and Params
    config = vcfdata_obj.get_config()
    param = vcfdata_obj.get_param()

    # Load args into param
    param = load_args(param=param, args=args, arguments_dict=arguments_dict, command="stats", strict=False)
    
    # Access
    config["access"] = "RO"

    # Re-Load Config and Params
    vcfdata_obj.set_param(param)
    vcfdata_obj.set_config(config)

    # Load data
    vcfdata_obj.load_data()

    # Parameters
    stats_md = param.get("stats",{}).get("stats_md",None)
    stats_json = param.get("stats",{}).get("stats_json",None)

    # Stats
    vcfdata_obj.print_stats(output_file=stats_md, json_file=stats_json)

    # Close connexion
    vcfdata_obj.close_connexion()

    log.info("End")



