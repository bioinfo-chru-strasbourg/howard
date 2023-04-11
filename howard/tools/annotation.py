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
from howard.tools.databases import *



def annotation(args) -> None:

    log.info("Start")

    config = args.config

    # # Load parameters in JSON format
    # if os.path.exists(args.param):
    #     with open(args.param) as param_file:
    #         param = json.load(param_file)
    # else:
    #     param = json.loads(args.param)

    param = {}

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

        # # Quick calculations
        # if args.calculations:
        #     #calculations_list= args.calculations.split(",")
        #     calculations_list= [value for val in args.calculations for value in val.split(',')]
        #     log.info(f"Quick Calculations list: {calculations_list}")
        #     param_quick_calculations = param.get("calculation",{})
        #     for calculation_operation in calculations_list:
        #         param_quick_calculations[calculation_operation] = {}
        #     params["calculation"] = param_quick_calculations

        # # Quick prioritization
        # if args.prioritizations:
        #     config_profiles= args.prioritizations
        #     log.info(f"Quick Prioritization Config file: {config_profiles}")
        #     param_quick_prioritizations = param.get("prioritization",{})
        #     param_quick_prioritizations["config_profiles"] = config_profiles
        #     params["prioritization"] = param_quick_prioritizations
        
        vcfdata_obj.set_param(params)
            

        # Load data from input file
        vcfdata_obj.load_data()

        # Overview
        # if args.overview:
        #     vcfdata_obj.get_overview()

        # # Stats
        # if args.stats:
        #     vcfdata_obj.get_stats()

        # Annotation
        if vcfdata_obj.get_param().get("annotations", None) or vcfdata_obj.get_param().get("annotation", None):
            vcfdata_obj.annotation()

        # # Calculation
        # if vcfdata_obj.get_param().get("calculations", None) or vcfdata_obj.get_param().get("calculation", None):
        #     vcfdata_obj.calculation()

        # # Prioritization
        # if vcfdata_obj.get_param().get("prioritizations", None) or vcfdata_obj.get_param().get("prioritization", None):
        #     vcfdata_obj.prioritization()

        # Export
        if vcfdata_obj.get_output():
            log.info("Exporting...")
            vcfdata_obj.export_output(export_header=True)

        # Close connexion
        vcfdata_obj.close_connexion()

    log.info("End")



