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



def stats(args) -> None:

    log.info("Start")

    config = args.config

    param = {}

    # Create VCF object
    if args.input:
        vcfdata_obj = Variants(None, args.input, config=config, param=param)

        # Load data from input file
        vcfdata_obj.load_data()

        # Stats
        vcfdata_obj.get_stats()

        # Close connexion
        vcfdata_obj.close_connexion()

    log.info("End")



