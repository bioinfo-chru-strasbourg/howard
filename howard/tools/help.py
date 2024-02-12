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
from configparser import ConfigParser

from howard.objects.variants import Variants
from howard.objects.database import Database
from howard.commons import *
from howard.tools.databases import *
from howard.tools.tools import *


main_folder = os.path.dirname(__file__)


def help(args:argparse) -> None:
    """
    The `help` function generates help documentation in various formats (parser, Markdown, HTML) based
    on the provided arguments and setup configuration.
    
    :param args: The `args` parameter is of type `argparse.Namespace`. It is used to pass command-line
    arguments to the `help` function. The `argparse` module provides a way to parse command-line
    arguments and generate help messages. The `Namespace` object holds the values of the command-line
    arguments
    :type args: argparse
    """

    # Config infos
    arguments_dict = args.arguments_dict
    setup_cfg = args.setup_cfg
    
    parser = help_generation(arguments_dict=arguments_dict, setup=setup_cfg, output_type="parser")
    parser.print_help()
    print("")

    log.info("Start")

    # MarkDown file
    if "help_md" in args and args.help_md:
        help_file = args.help_md.name
        log.info(f"Help - generate Markdown help file ['{help_file}']")
        help_content = help_generation(arguments_dict=arguments_dict, setup=setup_cfg, output_type="markdown")
        f = open(help_file, "w")
        f.write(help_content)
        f.close()

    # HTML file
    if "help_html" in args and args.help_html:
        help_file = args.help_html.name
        log.info(f"Help - generate HTML help file ['{help_file}']")
        help_content = help_generation(arguments_dict=arguments_dict, setup=setup_cfg, output_type="html")
        f = open(help_file, "w")
        f.write(help_content)
        f.close()

    log.info("End")
