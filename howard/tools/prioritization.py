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



def prioritization(args) -> None:

    log.info("Start")

    config = args.config

    param = {}

    # Create VCF object
    if args.input:
        vcfdata_obj = Variants(None, args.input, args.output, config, param)

        params = vcfdata_obj.get_param()

        # Quick prioritization
        if args.prioritizations:
            config_profiles= args.prioritizations
            log.info(f"Quick Prioritization Config file: {config_profiles}")
            param_quick_prioritizations = param.get("prioritization",{})
            param_quick_prioritizations["config_profiles"] = config_profiles
            params["prioritization"] = param_quick_prioritizations

            # Profiles
            if args.profiles:
                params["prioritization"]["profiles"] = [value for val in args.profiles for value in val.split(',')]

            # PZFields
            if args.pzfields:
                params["prioritization"]["pzfields"] = [value for val in args.pzfields for value in val.split(',')]

            # Default
            if args.default_profile:
                params["prioritization"]["default_profile"] = args.default_profile

            # Score Mode
            if args.prioritization_score_mode:
                params["prioritization"]["prioritization_score_mode"] = args.prioritization_score_mode

        
        vcfdata_obj.set_param(params)
            
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

    log.info("End")



