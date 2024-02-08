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
# howard analysis --input=my.vcf.gz --output=my.output.vcf --annotations=my.annotations.vcf.gz --stats --overview_footer
# python -m howard.main --input=my.vcf.gz --output=my.output.vcf --annotations=my.annotations.vcf.gz --stats --overview_footer


main_folder = os.path.dirname(__file__)

# Main function
def main() -> None:
    """
    It loads a VCF file in multiple format (VCF, parquet, DB), and process, query, export data
    """
    
    # Config Parser
    cf = ConfigParser()
    cf.read(f'{main_folder}/../setup.cfg')
    prog_name = cf['metadata']['name']
    prog_version = cf['metadata']['version']
    prog_description = cf['metadata']['description']
    prog_long_description_content_type = cf['metadata']['long_description_content_type']
    # prog_author = cf['metadata']['author']
    # prog_author_email = cf['metadata']['author_email']

    # Main parser
    parser = argparse.ArgumentParser(
        prog="howard",
        description=f"""{prog_name.upper()}:{prog_version}\n"""
                    f"""{prog_description}\n"""
                    f"""{prog_long_description_content_type}""",
        epilog="Usage examples:\n"
            """   howard process --input=tests/data/example.vcf.gz --output=/tmp/example.annotated.vcf.gz --param=config/param.json \n"""
            """   howard annotation --input=tests/data/example.vcf.gz --output=/tmp/example.howard.vcf.gz --annotations='tests/databases/annotations/hg19/dbnsfp42a.parquet,tests/databases/annotations/hg19/gnomad211_genome.parquet' \n"""
            """   howard calculation --input=tests/data/example.full.vcf --output=/tmp/example.calculation.tsv --calculations='vartype' \n"""
            """   howard prioritization --input=tests/data/example.vcf.gz --output=/tmp/example.prioritized.vcf.gz --prioritizations=config/prioritization_profiles.json --profiles='default,GERMLINE' \n"""
            """   howard query --input=tests/data/example.vcf.gz --explode_infos --query='SELECT "#CHROM", POS, REF, ALT, "DP", "CLNSIG", sample2, sample3 FROM variants WHERE "DP" >= 50 OR "CLNSIG" NOT NULL ORDER BY "CLNSIG" DESC, "DP" DESC' \n"""
            """   howard stats --input=tests/data/example.vcf.gz \n"""
            """   howard convert --input=tests/data/example.vcf.gz --output=/tmp/example.tsv --explode_infos && cat /tmp/example.tsv \n"""
            ,
        formatter_class=argparse.RawTextHelpFormatter
        )
    
    parser._optionals.title = "Shared arguments"

    subparsers = parser.add_subparsers(title="Tools", dest='command')

    options_md = ""
    options_html = ""
    
    options_md += f"# HOWARD Help"
    options_md += "\n\n"
    options_md += parser.description.replace("\n","\n\n")
    options_md += re.sub(r'> $',"",parser.epilog.replace("\n","\n\n").replace("   ","> "))
    options_md += "\n\n"

    options_html += f"<H1>HOWARD Help</h1>\n"
    options_html += "<p>" + parser.description.replace("\n","<br>") + "</p>"
    options_html += parser.epilog.replace("\n","<br>").replace("   ","&nbsp;&nbsp;&nbsp;")
        

    # Create commands arguments
    for command in commands_arguments:
        
        command_description = commands_arguments[command].get("description","")
        command_help = commands_arguments[command].get("help","")
        command_epilog = commands_arguments[command].get("epilog","")
        
        command_parser = subparsers.add_parser(
            command,
            description = command_description,
            help = command_help,
            epilog = command_epilog,
            formatter_class=argparse.RawTextHelpFormatter
        )

        # Markdown
        options_md += f"## {command.upper()} tool\n"
        options_md += command_description.replace("\n","\n\n")
        options_md += "\n\n"
        options_md += re.sub(r'> $',"",command_epilog.replace("\n","\n\n").replace("   ","> "))
        options_md += "\n\n"

        # HTML
        options_html += f"<H2>{command.upper()}</H2>\n"
        options_html += "<p>" + command_description.replace("\n","<br>") + "</p>"
        options_html += command_epilog.replace("\n","<br>").replace("   ","&nbsp;&nbsp;&nbsp;")
        
        
        # Main args
        command_parser._optionals.title = "Options"
        if "main" in commands_arguments[command]["groups"]:
            group = "main"
            options_md += f"### Main options\n"
            options_html += f"<H3>Main options</H3>\n"
            for arg in commands_arguments[command]["groups"][group]:
                required = commands_arguments[command]["groups"][group][arg]
                argument = get_argument(arguments=arguments.copy(), arg=arg, required=required)
                command_parser.add_argument(f"--{arg}", **argument)
                options_md += get_argument_to_mk(arg, argument)
                options_html += get_argument_to_mk(arg, argument, mode="html")

        for group in commands_arguments[command]["groups"]:
            if group != "main":
                options_md += f"### {group}\n"
                options_html += f"<H3>{group}</H3>\n"
                command_group = command_parser.add_argument_group(f"{group} options")
                for arg in commands_arguments[command]["groups"][group]:
                    required = commands_arguments[command]["groups"][group][arg]
                    argument = get_argument(arguments=arguments.copy(), arg=arg, required=required)
                    command_group.add_argument(f"--{arg}", **argument)
                    options_md += get_argument_to_mk(arg, argument)
                    options_html += get_argument_to_mk(arg, argument, mode="html")

        # Shared arguments
        shared_group = command_parser.add_argument_group('Shared options')
        for arg in shared_arguments:
            shared_group.add_argument(f"--{arg}", **get_argument(arguments=arguments, arg=arg, required=False))

        #options_mk += "> " + re.sub(r'> $',"",command_epilog.replace("\n","\n\n> "))
        options_md += "\n\n"


    options_md += f"## Shared arguments\n"
    options_html += f"<H2>Shared arguments</H2>\n".upper()
    for arg in shared_arguments:
        required = False
        argument = get_argument(arguments=arguments.copy(), arg=arg, required=required)
        options_md += get_argument_to_mk(arg, argument)
        options_html += get_argument_to_mk(arg, argument, mode="html")
 

    # Parse args
    args, remaining = parser.parse_known_args()

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
    elif args.command == "help":
        if "help_md" in args and args.help_md is not None:
            help_file = args.help_md.name
            f = open(help_file, "w")
            f.write(options_md)
            f.close()
        if "help_html" in args and args.help_html is not None:
            help_file = args.help_html.name
            f = open(help_file, "w")
            f.write(options_html)
            f.close()
        if "help_md" not in args and "help_html" not in args:
            parser.print_help()
        return
    else:
        command_function = commands_arguments[args.command]["function"]
        log.debug(f"Command/Tool: {command_function}")
        eval(f"{command_function}(args)")


if __name__ == '__main__':
    main()

