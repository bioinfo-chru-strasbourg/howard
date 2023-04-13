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



def calculation(args) -> None:

    log.info("Start")

    config = args.config

    param = {}

    # Create VCF object
    if args.input:
        vcfdata_obj = Variants(None, args.input, args.output, config, param)

        if args.show_calculations:
            for help_line in vcfdata_obj.get_operations_help():
                log.info(help_line)
            return 

        params = vcfdata_obj.get_param()

        # Quick calculations
        if args.calculations:
            calculations_list= [value for val in args.calculations for value in val.split(',')]
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

    log.info("End")



