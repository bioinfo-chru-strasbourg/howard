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



def annotation(args:argparse) -> None:
    """
    The `annotation` function performs annotation on a VCF file based on specified parameters and
    exports the annotated data.
    
    :param args: The `args` parameter is likely an object or dictionary containing various arguments
    passed to the `annotation` function. It is not clear from the code snippet what specific arguments
    are expected or required
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
    if args.annotations or param.get("annotations", None) or param.get("annotation", None):
        vcfdata_obj = Variants(None, args.input, args.output, config, param)

        param = vcfdata_obj.get_param()

        # Prapare annotation dict
        if not param.get("annotation", None):
            param["annotation"] = {}
        if not param.get("annotation", {}).get("options", None):
            param["annotation"]["options"] = {}

        # Quick Annotation
        if args.annotations:
            annotation_file_list = [value for value in args.annotations.split(',')]
            log.info(f"Quick Annotation Files: {annotation_file_list}")
            param_quick_annotations = param.get("annotations",{})
            for annotation_file in annotation_file_list:
                param_quick_annotations[annotation_file] = {"INFO": None}
            param["annotations"] = param_quick_annotations
        
        if args.annotations_update:
            param["annotation"]["options"]["annotations_update"] = True

        if args.annotations_append:
            param["annotation"]["options"]["annotations_append"] = True

        vcfdata_obj.set_param(param)
        
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

    else:
        # Parser
        parser = help_generation(arguments_dict=arguments_dict, setup=setup_cfg, output_type="parser")
        parser.print_help()
        print("")
        log.error(f"No annotations provided")
        raise ValueError(f"No annotations provided")

    log.info("End")



