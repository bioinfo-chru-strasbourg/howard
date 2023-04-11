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

    param = {}


    # Create VCF object
    if args.input:
        vcfdata_obj = Variants(None, args.input, args.output, config, param)

        params = vcfdata_obj.get_param()

        # Quick Annotation
        if args.annotations:
            annotation_file_list = [value for val in args.annotations for value in val.split(',')]
            log.info(f"Quick Annotation Files: {annotation_file_list}")
            param_quick_annotations = param.get("annotations",{})
            for annotation_file in annotation_file_list:
                param_quick_annotations[annotation_file] = {"INFO": None}
            params["annotations"] = param_quick_annotations
        
        vcfdata_obj.set_param(params)
            

        # Load data from input file
        vcfdata_obj.load_data()

        # Annotation
        if vcfdata_obj.get_param().get("annotations", None) or vcfdata_obj.get_param().get("annotation", None):
            vcfdata_obj.annotation()

        # Export
        if vcfdata_obj.get_output():
            log.info("Exporting...")
            vcfdata_obj.export_output(export_header=True)

        # Close connexion
        vcfdata_obj.close_connexion()

    log.info("End")



