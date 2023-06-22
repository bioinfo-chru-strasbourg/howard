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



def calculation(args:argparse) -> None:
    """
    This function performs calculations on VCF data based on user input and exports the results.
    
    :param args: The `args` parameter is a command line argument parser object that contains the
    arguments passed to the script when it was executed
    :type args: argparse
    """

    log.info("Start")

    config = args.config

    param = {}

    # Show available calculations
    if args.show_calculations:
        vcfdata_obj = Variants()
        for help_line in vcfdata_obj.get_operations_help():
            log.info(help_line)

    # Create VCF object
    elif args.input and args.output and args.calculations:

        vcfdata_obj = Variants(None, args.input, args.output, config, param)

        params = vcfdata_obj.get_param()

        # Quick calculations
        if args.calculations:
            calculations_list= [value for value in args.calculations.split(',')]
            log.info(f"Quick Calculations list: {calculations_list}")
            param_quick_calculations = param.get("calculation",{})
            for calculation_operation in calculations_list:
                param_quick_calculations[calculation_operation.upper()] = {}
            params["calculation"] = param_quick_calculations

        # HGVS Field
        if args.hgvs_field and "NOMEN" in params["calculation"]:
            if "options" not in params["calculation"]["NOMEN"]:
                params["calculation"]["NOMEN"]["options"] = {}
            params["calculation"]["NOMEN"]["options"]["hgvs_field"] = args.hgvs_field

        # HGVS Transcripts
        if args.transcripts and "NOMEN" in params["calculation"]:
            if "options" not in params["calculation"]["NOMEN"]:
                params["calculation"]["NOMEN"]["options"] = {}
            params["calculation"]["NOMEN"]["options"]["transcripts"] = args.transcripts

        # TRIO pedigree
        if args.trio_pedigree and "TRIO" in params["calculation"]:
            trio_pedigree = {}
            # Load trio_pedigree in JSON format
            if os.path.exists(args.trio_pedigree):
                with open(args.trio_pedigree) as trio_pedigree_file:
                    trio_pedigree = json.load(trio_pedigree_file)
            else:
                trio_pedigree = json.loads(args.trio_pedigree)
            params["calculation"]["TRIO"] = trio_pedigree

        vcfdata_obj.set_param(params)

        # Load data from input file
        vcfdata_obj.load_data()

        # Calculation
        if vcfdata_obj.get_param().get("calculations", None) or vcfdata_obj.get_param().get("calculation", None):
            vcfdata_obj.calculation()

        # Export
        if vcfdata_obj.get_output():
            log.info("Exporting...")
            vcfdata_obj.export_output(export_header=True)

        # Close connexion
        vcfdata_obj.close_connexion()

    # If no arguments
    else:

        log.info("""The following arguments are required:""")
        log.info("""   --input, --output, --calculations""")
        log.info("""   --show_calculations""")

    log.info("End")



