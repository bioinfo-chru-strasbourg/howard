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


# Import Gooey only if exists
import imp
try:
    imp.find_module('gooey')
    from gooey import Gooey, GooeyParser
except ImportError:
    found = False

# Images folder
main_folder = os.path.dirname(__file__)
image_dir = main_folder + "/../../images/"

# Gooey config
@Gooey(advanced=True,          # toggle whether to show advanced config or not 
    language="english",  # Translations configurable via json
    auto_start=False,           # skip config screens all together
    #target="",     # Explicitly set the subprocess executable arguments
    program_name='HOWARD',       # Defaults to script name
    program_description='HOWARD Graphical User Interface',       # Defaults to ArgParse Description
    default_size=(1200, 800),   # starting size of the GUI
    required_cols=2,           # number of columns in the "Required" section
    optional_cols=2,           # number of columns in the "Optional" section
    dump_build_config=False,   # Dump the JSON Gooey uses to configure itself
    load_build_config=None,    # Loads a JSON Gooey-generated configuration
    monospace_display=True,   # Uses a mono-spaced font in the output screen
    navigation="SIDEBAR",    # SIDEBAR TABBED
    show_sidebar=True,
    sidebar_title="Commands",
    tabbed_groups=True,
    image_dir=image_dir
)

# GUI function
def gui(args:argparse) -> None:
    """
    The `gui` function generates a graphical user interface (GUI) for a Python script using the
    `argparse` module and the `Gooey` library.
    
    :param args: The `args` parameter is of type `argparse`, which is a module in Python used for
    parsing command-line arguments. It is used to define the arguments that the program accepts and to
    generate help messages. In this code, it seems that `args` is an object that contains information
    about the
    :type args: argparse
    """

    # Config infos
    arguments_dict = args.arguments_dict
    setup_cfg = args.setup_cfg

    # Parser Gooey
    parser = GooeyParser()
    parser_gooey = help_generation(arguments_dict=arguments_dict, parser=parser, setup=setup_cfg, output_type="gooey")
    parser_gooey.print_help()
    print("")

    log.info("Start")

    args = parser_gooey.parse_args()

    log.info("End")