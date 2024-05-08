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


def hgvs(args: argparse) -> None:
    """
    The `hgvs` function takes command line arguments, creates a VCF object, sets parameters and
    configurations, loads data from an input file, performs annotation using HGVS notation, exports the
    output, and closes the connection.

    :param args: The `args` parameter is of type `argparse.Namespace` and is used to parse command line
    arguments. It contains the following attributes:
    :type args: argparse
    """

    log.info("Start")

    # Load config args
    arguments_dict, setup_cfg, config, param = load_config_args(args)

    # Create variants object
    vcfdata_obj = Variants(
        input=args.input, output=args.output, config=config, param=param
    )

    # Get Config and Params
    config = vcfdata_obj.get_config()
    param = vcfdata_obj.get_param()

    # Load args into param
    param = load_args(
        param=param,
        args=args,
        arguments_dict=arguments_dict,
        command="hgvs",
        strict=False,
    )

    # Re-Load Config and Params
    vcfdata_obj.set_param(param)
    vcfdata_obj.set_config(config)

    # Load data
    vcfdata_obj.load_data()

    # Prioritization
    vcfdata_obj.annotation_hgvs()

    # Export
    vcfdata_obj.export_output(export_header=True)

    # Close connexion
    vcfdata_obj.close_connexion()

    log.info("End")
