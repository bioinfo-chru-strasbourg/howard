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

    # Load config args
    arguments_dict, setup_cfg, config, param = load_config_args(args)

    # Create variants object
    vcfdata_obj = Variants(input=args.input, output=args.output, config=config, param=param)

    # Get Config and Params
    config = vcfdata_obj.get_config()
    param = vcfdata_obj.get_param()

    # Load args into param
    param = load_args(param=param, args=args, arguments_dict=arguments_dict, command="annotation", strict=False)

    # Re-Load Config and Params
    vcfdata_obj.set_param(param)
    vcfdata_obj.set_config(config)

    # Load data
    vcfdata_obj.load_data()

    # Annotation
    vcfdata_obj.annotation()

    # Export
    vcfdata_obj.export_output()

    # Close connexion
    vcfdata_obj.close_connexion()

    log.info("End")