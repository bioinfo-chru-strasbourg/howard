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



def prioritization(args:argparse) -> None:
    """
    The function performs prioritization on a VCF file based on user-specified configurations and
    exports the results.
    
    :param args: args is an object that contains the command line arguments passed to the script. It is
    used to configure the behavior of the script and to provide input and output file paths, as well as
    other parameters needed for the execution of the script
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
    param = load_args(param=param, args=args, arguments_dict=arguments_dict, command="prioritization", strict=False)
    
    # Re-Load Config and Params
    vcfdata_obj.set_param(param)
    vcfdata_obj.set_config(config)

    # Load data
    if vcfdata_obj.get_input():
        vcfdata_obj.load_data()

    log.debug(f"param={param}")

    # if param.get("calculation",{}).get("show_calculations",False):
    #     operations_config_file = param.get("calculation").get("calculation_config")
    #     log.debug(f"operations_config_file={operations_config_file}")
    #     for help_line in vcfdata_obj.get_operations_help(operations_config_file=operations_config_file):
    #         log.info(help_line)
    #     exit()

    # Annotation
    vcfdata_obj.prioritization()

    # Export
    vcfdata_obj.export_output()

    # Close connexion
    vcfdata_obj.close_connexion()

    log.info("End")



