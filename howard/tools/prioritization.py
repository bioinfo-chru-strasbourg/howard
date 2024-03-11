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

    # Create VCF object
    #if args.input:
    if args.prioritizations or param.get("prioritizations", None) or param.get("prioritization", None):
        vcfdata_obj = Variants(None, args.input, args.output, config, param)

        param = vcfdata_obj.get_param()

        # Quick prioritization
        if args.prioritizations:
            if isinstance(args.prioritizations, str):
                config_profiles = args.prioritizations
            else:
                config_profiles = args.prioritizations.name
            #config_profiles= args.prioritizations
            log.info(f"Quick Prioritization Config file: {config_profiles}")
            param_quick_prioritizations = param.get("prioritization",{})
            param_quick_prioritizations["config_profiles"] = config_profiles
            param["prioritization"] = param_quick_prioritizations

            # Profiles
            if args.profiles:
                param["prioritization"]["profiles"] = [value for value in args.profiles.split(',')]

            # PZFields
            if args.pzfields:
                param["prioritization"]["pzfields"] = [value for value in args.pzfields.split(',')]

            # Default
            if args.default_profile:
                param["prioritization"]["default_profile"] = args.default_profile

            # Score Mode
            if args.prioritization_score_mode:
                param["prioritization"]["prioritization_score_mode"] = args.prioritization_score_mode

            # Config profiles
            if config_profiles in args and args.config_profiles:
                if isinstance(args.config_profiles, str):
                    config_profiles_file = args.config_profiles
                else:
                    config_profiles_file = args.config_profiles.name
                param["prioritization"]["config_profiles"] = config_profiles_file

        vcfdata_obj.set_param(param)
            
        # Load data from input file
        vcfdata_obj.load_data()

        # Prioritization
        if vcfdata_obj.get_param().get("prioritizations", None) or vcfdata_obj.get_param().get("prioritization", None):
            vcfdata_obj.prioritization()

        # Export
        if vcfdata_obj.get_output():
            log.info("Exporting...")
            vcfdata_obj.export_output(export_header=True)

        # Close connexion
        vcfdata_obj.close_connexion()

    else:
        # Parser
        parser = help_generation(arguments_dict=arguments_dict, setup=setup_cfg, output_type="parser")
        parser.print_help()
        print("")
        log.error(f"No prioritizations provided")
        raise ValueError(f"No prioritizations provided")

    log.info("End")



