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



def hgvs(args:argparse) -> None:
    """
    The `hgvs` function takes command line arguments, creates a VCF object, sets parameters and
    configurations, loads data from an input file, performs annotation using HGVS notation, exports the
    output, and closes the connection.
    
    :param args: The `args` parameter is of type `argparse.Namespace` and is used to parse command line
    arguments. It contains the following attributes:
    :type args: argparse
    """

    log.info("Start")

    config = args.config

    params = {}
    
    # Create VCF object
    if args.input:
        vcfdata_obj = Variants(None, args.input, args.output, config, params)
        
        # Params
        params = vcfdata_obj.get_param()
        config = vcfdata_obj.get_config()

        if not params.get("hgvs", None):
            params["hgvs"] = {}

        if "use_exon" in args:
            params["hgvs"]["use_exon"] = args.use_exon

        if "use_gene" in args:
            params["hgvs"]["use_gene"] = args.use_gene

        if "use_protein" in args:
            params["hgvs"]["use_protein"] = args.use_protein

        if "add_protein" in args:
            params["hgvs"]["add_protein"] = args.add_protein

        if "full_format" in args:
            params["hgvs"]["full_format"] = args.full_format
        
        if "use_version" in args:
            params["hgvs"]["use_version"] = args.use_version

        if "codon_type" in args:
            params["hgvs"]["codon_type"] = args.codon_type

        if "refgene" in args:
            params["hgvs"]["refgene"] = args.refgene

        if "refseqlink" in args:
            params["hgvs"]["refseqlink"] = args.refseqlink

        if "assembly" in args:
            params["assembly"] = args.assembly

        if "genomes_folder" in args and args.genomes_folder:
            if "folders" not in config:
                config["folders"] = {}
            if "databases" not in config["folders"]:
                config["folders"]["databases"] = {}
            config["folders"]["databases"]["genomes"] = args.genomes_folder

        if "refseq_folder" in args and args.refseq_folder:
            if "folders" not in config:
                config["folders"] = {}
            if "databases" not in config["folders"]:
                config["folders"]["databases"] = {}
            config["folders"]["databases"]["refseq"] = args.refseq_folder

        vcfdata_obj.set_param(params)
        vcfdata_obj.set_config(config)
        
        # Load data from input file
        vcfdata_obj.load_data()

        # Prioritization
        #if vcfdata_obj.get_param().get("hgvs", None):
        vcfdata_obj.annotation_hgvs()
        
        # Export
        if vcfdata_obj.get_output():
            #log.info("Exporting...")
            vcfdata_obj.export_output(export_header=True)

        # Close connexion
        vcfdata_obj.close_connexion()

    log.info("End")



