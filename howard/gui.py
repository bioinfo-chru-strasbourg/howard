from gooey import Gooey, GooeyParser
#from message import display_message

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

from howard.objects.variants import Variants
from howard.objects.database import Database
from howard.commons import *

from howard.tools.tools import *


@Gooey(advanced=True,          # toggle whether to show advanced config or not 
       language="english",  # Translations configurable via json
       #auto_start=True,           # skip config screens all together
       #target="",     # Explicitly set the subprocess executable arguments
       program_name='HOWARD',       # Defaults to script name
       program_description='Graphical User Interface of HOWARD',       # Defaults to ArgParse Description
       default_size=(800, 600),   # starting size of the GUI
       # required_cols=2,           # number of columns in the "Required" section
       # optional_cols=2,           # number of columns in the "Optional" section
       dump_build_config=False,   # Dump the JSON Gooey uses to configure itself
       load_build_config=None,    # Loads a JSON Gooey-generated configuration
       #monospace_display=False,   # Uses a mono-spaced font in the output screen
       navigation="TABBED",    # SIDEBAR TABBED
       show_sidebar=True,
       sidebar_title="Commands",
       #tabbed_groups=True,
)

def main() -> None:
    """
    It loads a VCF file in multiple format (VCF, parquet, DB), and process, query, export data
    """

    # Main parser
    desc = "Graphical User Interface of HOWARD ?"
    parser = GooeyParser(description=desc)

    parser._optionals.title = "Shared arguments"

    subparsers = parser.add_subparsers(title="Tools", dest='command')

    # Create commands arguments
    for command in commands_arguments:
        command_parser = subparsers.add_parser(
            command,
            description = commands_arguments[command].get("description",""),
            help = commands_arguments[command].get("help",""),
            epilog = commands_arguments[command].get("epilog",""),
            formatter_class=argparse.RawTextHelpFormatter
        )
        # Main args
        command_parser._optionals.title = "Main Options"
        if "main" in commands_arguments[command]["groups"]:
            for arg in commands_arguments[command]["groups"]["main"]:
                required = commands_arguments[command]["groups"]["main"][arg]
                arg_infos = get_argument(arguments=arguments, arg=arg, required=required)
                widget, options = get_argument_gooey(arg=arg)
                command_parser.add_argument(f"--{arg}", **arg_infos, widget=widget, gooey_options=options)

        for group in commands_arguments[command]["groups"]:
            if group != "main":
                command_group = command_parser.add_argument_group(f"{group} options")
                for arg in commands_arguments[command]["groups"][group]:
                    required = commands_arguments[command]["groups"][group][arg]
                    arg_infos = get_argument(arguments=arguments, arg=arg, required=required)
                    widget, options = get_argument_gooey(arg=arg)
                    command_group.add_argument(f"--{arg}", **arg_infos, widget=widget, gooey_options=options)

        # Shared arguments
        shared_group = command_parser.add_argument_group('Shared options')
        for arg in shared_arguments:
            arg_infos = get_argument(arguments=arguments, arg=arg, required=required)
            widget, options = get_argument_gooey(arg=arg)
            #shared_group.add_argument(f"--{arg}", **get_argument(arguments=arguments, arg=arg, required=False))
            shared_group.add_argument(f"--{arg}", **arg_infos, widget=widget, gooey_options=options)


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
        return
    else:
        command_function = commands_arguments[args.command]["function"]
        log.debug(f"Command/Tool: {command_function}")
        eval(f"{command_function}(args)")


if __name__ == '__main__':
    main()