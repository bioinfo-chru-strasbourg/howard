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
import psutil
import markdown
from configparser import ConfigParser

from howard.objects.variants import Variants
from howard.objects.database import Database
from howard.commons import *

from howard.tools.tools import *


# Usage
# python -m pip install -e .
# howard --help
# howard query --help
# howard analysis --input=my.vcf.gz --output=my.output.vcf --annotations=my.annotations.vcf.gz --stats --overview_footer
# howard gui
# python -m howard.main --input=my.vcf.gz --output=my.output.vcf --annotations=my.annotations.vcf.gz --stats --overview_footer


main_folder = os.path.dirname(__file__)

# Main function
def main() -> None:
    """
    It loads a VCF file in multiple format (VCF, parquet, DB), and process, query, export data
    """

    # Config infos
    setup_cfg = f'{main_folder}/../setup.cfg'
    arguments_dict = {
        "arguments": arguments,
        "commands_arguments": commands_arguments,
        "shared_arguments": shared_arguments 
    }
    
    # Generate parser
    parser = argparse.ArgumentParser()
    parser = help_generation(arguments_dict=arguments_dict, parser=parser, setup=setup_cfg, output_type="parser")

    # Parse args
    args = parser.parse_args()

    # Quiet
    if "quiet" in args and args.quiet:
        args.verbosity = "warning"
    # Verbose
    if "verbose" in args and args.verbose:
        args.verbosity = "info"
    # Debug
    if "debug" in args and args.debug:
        args.verbosity = "debug"
    # verbosity
    if "verbosity" not in args:
        args.verbosity = "info"
    # log
    if "log" not in args:
        args.log = None

    # Config infos
    args.arguments_dict = arguments_dict
    args.setup_cfg = setup_cfg

    # Logging
    set_log_level(args.verbosity, args.log)

    # Threads
    nb_threads = os.cpu_count()
    if "threads" in args and args.threads:
        threads = args.threads
    else:
        threads = nb_threads
    if threads == -1:
        threads = nb_threads

    # Memory
    mem = psutil.virtual_memory()
    mem_total = mem.total / 1024 / 1024 / 1024
    mem_default = int(mem_total * 0.8)
    if mem_default < 1:
        mem_default = 1
    if "memory" in args and args.memory:
        memory = args.memory
    else:
        memory = f"{mem_default}G"

    # chunk_size
    chunk_size = DEFAULT_CHUNK_SIZE
    if "chunk_size" in args and args.chunk_size:
        chunk_size = int(args.chunk_size)

    # Temporary folder
    tmp = None
    if "tmp" in args and args.tmp:
        tmp = args.tmp
    
    # duckdb settings
    duckdb_settings = None
    if "duckdb_settings" in args and args.duckdb_settings:
        duckdb_settings = args.duckdb_settings
        
    # Assembly
    if "assembly" in args and args.assembly:
        assembly = args.assembly
    else:
        assembly = "hg19"

    # Load configuration in JSON format
    if "config" in args:
        if os.path.exists(args.config):
            with open(args.config) as config_file:
                config = json.load(config_file)
        else:
            config = json.loads(args.config)
    else:
        config = {}
    
    # add to config
    config["verbosity"] = args.verbosity
    config["threads"] = threads
    config["memory"] = memory
    config["chunk_size"] = chunk_size
    config["tmp"] = tmp
    config["duckdb_settings"] = duckdb_settings
    config["assembly"] = assembly

    # Change config
    args.config = config
    log.debug(f"config: {config}")

    # Command eval
    if not args.command:
        parser.print_help()
    else:
        if args.command == "gui" and not tool_gui_enable:
            #from gooey import Gooey, GooeyParser
            log.error("""HOWARD GUI disabled""")
            log.error("""ModuleNotFoundError: No module named 'gooey'""")
            log.error("""Install module 'gooey': pip install gooey""")
            log.error("""Or install requirements: pip install -r requirements-gui.txt""")
            raise ValueError("""HOWARD GUI disabled""")
        command_function = commands_arguments[args.command]["function"]
        log.debug(f"Command/Tool: {command_function}")
        eval(f"{command_function}(args)")

if __name__ == '__main__':
    main()

