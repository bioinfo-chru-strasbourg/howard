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

from howard.tools.tools import *


# Usage
# python -m pip install -e .
# howard analysis --input=my.vcf.gz --output=my.output.vcf --annotations=my.annotations.vcf.gz --stats --overview_footer
# python -m howard.main --input=my.vcf.gz --output=my.output.vcf --annotations=my.annotations.vcf.gz --stats --overview_footer


# Main function
def main() -> None:
    """
    It loads a VCF file in multiple format (VCF, parquet, DB), and process, query, export data
    """

    # Main parser
    parser = argparse.ArgumentParser(
        prog="howard",
        description="""howard annotates and prioritizes genetic variations, calculates and normalizes annotations, translates vcf format and generates variants statistics""",
        #usage="howard [<shared-args>]",
        epilog="Examples:\n"
            """   howard process --input=input.vcf.gz --output=output.tsv \n"""
            """   howard query --input=input.vcf.gz --query="SELECT * FROM variants WHERE REF = 'A' AND POS < 100000" \n""",
        formatter_class=argparse.RawTextHelpFormatter
        )

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
        command_parser._optionals.title = "Options"
        if "main" in commands_arguments[command]["groups"]:
            for arg in commands_arguments[command]["groups"]["main"]:
                required = commands_arguments[command]["groups"]["main"][arg]
                command_parser.add_argument(f"--{arg}", **get_argument(arguments=arguments, arg=arg, required=required))

        for group in commands_arguments[command]["groups"]:
            if group != "main":
                command_group = command_parser.add_argument_group(f"{group} options")
                for arg in commands_arguments[command]["groups"][group]:
                    required = commands_arguments[command]["groups"][group][arg]
                    command_group.add_argument(f"--{arg}", **get_argument(arguments=arguments, arg=arg, required=required))

        # Shared arguments
        shared_group = command_parser.add_argument_group('Shared options')
        for arg in shared_arguments:
            shared_group.add_argument(f"--{arg}", **get_argument(arguments=arguments, arg=arg, required=False))


    # Parse args
    args, remaining = parser.parse_known_args()

    # Verbosity
    # Default
    args.verbosity = "info"
    # Quiet
    if "quiet" in args and args.quiet:
        args.verbosity = "warning"
    # Verbose
    if "verbose" in args and args.verbose:
        args.verbosity = "info"
    # Debug
    if "debug" in args and args.debug:
        args.verbosity = "debug"

    # Logging
    set_log_level(args.verbosity)

    # Threads
    if "threads" in args and args.threads:
        threads = args.threads
    else:
        threads = 1

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


    # Change config
    args.config = config


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


    
