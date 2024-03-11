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



def convert(args:argparse) -> None:
    """
    The `convert` function converts a VCF file to a different format and can optionally explode info
    fields.
    
    :param args: `args` is a parameter passed to the `convert` function, likely an object or dictionary
    containing various arguments needed for the function to perform its task. These arguments could
    include things like input and output file paths, configuration settings, and other parameters
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
    if args.explode_infos:
        param["explode_infos"] = args.explode_infos
        param["explode_infos_prefix"] = args.explode_infos_prefix
        param["explode_infos_fields"] = args.explode_infos_fields
        
    else:
        config["access"] = "RO"

    # order_by
    if "order_by" in args and args.order_by:
        param["order_by"] = args.order_by

    # include_header
    if "include_header" in args and args.include_header:
        param["header_in_output"] = args.include_header

    # parquet_partitions
    if "parquet_partitions" in args and args.parquet_partitions:
        param["parquet_partitions"] = args.parquet_partitions.split(",")

    vcfdata_obj.set_param(param)
    vcfdata_obj.set_config(config)

    # Create VCF object
    if args.input:

        # Create
        vcfdata_obj = Variants(None, args.input, args.output, config=config, param=param, load=True)

        # Output
        vcfdata_obj.export_output()

        # Close connexion
        vcfdata_obj.close_connexion()

    log.info("End")



