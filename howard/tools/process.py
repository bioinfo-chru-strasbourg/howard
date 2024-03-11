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
    #if os.path.exists(args.param):
    if isinstance(args.param, str) and os.path.exists(full_path(args.param)):
        with open(full_path(args.param)) as param_file:
            param = json.load(param_file)
    else:
        param = json.loads(args.param)

    # Create VCF object
    if args.input:
        vcfdata_obj = Variants(None, args.input, args.output, config, param)

        param = vcfdata_obj.get_param()

        # Quick HGVS
        if "hgvs" in args and args.hgvs:
            log.info(f"Quick HGVS: {args.hgvs}")
            if not param.get("hgvs", None):
                param["hgvs"] = {}
            for option in args.hgvs.split(","):
                option_var_val = option.split(":")
                option_var = option_var_val[0]
                if len(option_var_val) > 1:
                    option_val = option_var_val[1]
                else:
                    option_val = "True"
                if option_val.upper() in ["TRUE"]:
                    option_val = True
                elif option_val.upper() in ["FALSE"]:
                    option_val = False
                param["hgvs"][option_var] = option_val

        # Quick Annotation
        if args.annotations:
            annotation_file_list = [value for value in args.annotations.split(',')]
            log.info(f"Quick Annotation: {annotation_file_list}")
            param_quick_annotations = param.get("annotations",{})
            for annotation_file in annotation_file_list:
                param_quick_annotations[annotation_file] = {"INFO": None}
            param["annotations"] = param_quick_annotations

        # Quick calculations
        if args.calculations:
            calculations_list= [value for value in args.calculations.split(',')]
            log.info(f"Quick Calculations: {calculations_list}")
            param_quick_calculations = param.get("calculation",{})
            for calculation_operation in calculations_list:
                param_quick_calculations[calculation_operation] = {}
            param["calculation"] = param_quick_calculations

        # Quick prioritization
        if args.prioritizations:
            if isinstance(args.prioritizations, str):
                config_profiles = args.prioritizations
            else:
                config_profiles = args.prioritizations.name
            #config_profiles= args.prioritizations
            log.info(f"Quick Prioritization: {config_profiles}")
            param_quick_prioritizations = param.get("prioritization",{})
            param_quick_prioritizations["config_profiles"] = config_profiles
            param["prioritization"] = param_quick_prioritizations

        # Quick query
        if args.query:
            param["query"] = args.query

        # Explode infos
        param["explode_infos"] = args.explode_infos
        param["explode_infos_prefix"] = args.explode_infos_prefix
        param["explode_infos_fields"] = args.explode_infos_fields

        # include_header
        if "include_header" in args and args.include_header:
            param["header_in_output"] = args.include_header

        # Set param
        vcfdata_obj.set_param(param)

        # Load data from input file
        vcfdata_obj.load_data()

        # Annotation
        if vcfdata_obj.get_param().get("hgvs", None):
            vcfdata_obj.annotation_hgvs()

        # Annotation
        if vcfdata_obj.get_param().get("annotations", None) or vcfdata_obj.get_param().get("annotation", None):
            vcfdata_obj.annotation()

        # Calculation
        if vcfdata_obj.get_param().get("calculations", None) or vcfdata_obj.get_param().get("calculation", None):
            vcfdata_obj.calculation()

        # Prioritization
        if vcfdata_obj.get_param().get("prioritizations", None) or vcfdata_obj.get_param().get("prioritization", None):
            vcfdata_obj.prioritization()

        # param["explode_infos"]
        if param.get("explode_infos", False):
            vcfdata_obj.explode_infos()

        # Query
        if param.get("query", None):
            log.info("Querying...")
            query = param.get("query", None)
            print(vcfdata_obj.get_query_to_df(query))

            # Output Query
            log.info("Exporting Querying...")
            vcfdata_obj.export_output(query=query, export_header=True)

        # Export Ouptut
        elif vcfdata_obj.get_output():
            log.info("Exporting...")
            vcfdata_obj.export_output(export_header=True)

        # Close connexion
        vcfdata_obj.close_connexion()

    log.info("End")



