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
from howard.commons import *
from howard.tools.databases import *



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
    config["access"] = "RO"

    # Load parameters in JSON format
    param = {}
    if "param" in args:
        if isinstance(args.param, str) and os.path.exists(full_path(args.param)):
            with open(full_path(args.param)) as param_file:
                param = json.load(param_file)
        else:
            param = json.loads(args.param)

    # MarkDown file
    stats_md=None
    if "stats_md" in args and args.stats_md:
        if isinstance(args.stats_md, str):
            stats_md = args.stats_md
        else:
            stats_md = args.stats_md.name
        param["stats_md"] = stats_md

    # JSON file
    stats_json=None
    if "stats_json" in args and args.stats_json:
        if isinstance(args.stats_json, str):
            stats_json = args.stats_json
        else:
            stats_json = args.stats_json.name
        param["stats_json"] = stats_json

    # Create VCF object
    if args.input:
        vcfdata_obj = Variants(None, args.input, config=config, param=param)

        # Load data from input file
        vcfdata_obj.load_data()

        # Stats
        vcfdata_obj.print_stats(output_file=stats_md, json_file=stats_json)

        # Close connexion
        vcfdata_obj.close_connexion()

    log.info("End")



