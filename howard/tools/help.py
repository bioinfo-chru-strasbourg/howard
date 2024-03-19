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
from howard.functions.commons import *
from howard.functions.databases import *
from howard.tools.tools import *


main_folder = os.path.dirname(__file__)


def help(args: argparse) -> None:
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
    if "arguments_dict" in args:
        arguments_dict = args.arguments_dict
    else:
        arguments_dict = None
    if "setup_cfg" in args:
        setup_cfg = args.setup_cfg
    else:
        setup_cfg = None

    # Parser
    parser = help_generation(
        arguments_dict=arguments_dict, setup=setup_cfg, output_type="parser"
    )
    parser.print_help()
    print("")

    log.info("Start")

    # Help JSON input
    if "help_json_input" in args and args.help_json_input:
        help_json_file = args.help_json_input
    else:
        help_json_file = None

    # Help JSON input title
    if "help_json_input_title" in args and args.help_json_input_title:
        help_json_input_title = args.help_json_input_title
    else:
        help_json_input_title = ""

    # Help example code type
    if "code_type" in args and args.code_type:
        code_type = args.code_type
    else:
        code_type = ""
    log.debug(f"code_type={code_type}")

    # MarkDown file
    if "help_md" in args and args.help_md:

        # Help file
        help_file = args.help_md
        log.info(f"Help - generate Markdown help file ['{help_file}']")

        # If Help input JSON file
        if help_json_file:
            log.info(f"Help -     from JSON help file ['{help_json_file}']")
            help_content = help_generation_from_json(
                help_json_file=help_json_file,
                output_type="markdown",
                title=help_json_input_title,
                code_type=code_type,
            )

        # Help from options
        else:
            help_content = help_generation(
                arguments_dict=arguments_dict, setup=setup_cfg, output_type="markdown"
            )

        # Write file
        f = open(help_file, "w")
        f.write(help_content)
        f.close()

        # Generate Table Of Content (if marker <!--TOC-->)
        import md_toc

        toc = md_toc.build_toc(help_file)
        md_toc.write_string_on_file_between_markers(help_file, toc, "<!--TOC-->")

    # HTML file
    if "help_html" in args and args.help_html:

        # Help file
        help_file = args.help_html
        log.info(f"Help - generate HTML help file ['{help_file}']")

        # If Help input JSON file
        if help_json_file:
            log.info(f"Help -     from JSON help file ['{help_json_file}']")
            help_content = help_generation_from_json(
                help_json_file=help_json_file,
                output_type="html",
                title=help_json_input_title,
            )

        # Help from options
        else:
            help_content = help_generation(
                arguments_dict=arguments_dict, setup=setup_cfg, output_type="html"
            )

        # Write file
        f = open(help_file, "w")
        f.write(help_content)
        f.close()

    log.info("End")
