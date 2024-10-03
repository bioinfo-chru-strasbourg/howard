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
import pypandoc

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

    # Help MD input
    if "help_md_input" in args and args.help_md_input:
        help_md_file = args.help_md_input
    else:
        help_md_file = None

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
    if ("help_md" in args and args.help_md) or help_json_file or help_md_file:

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

        # If Help input MD file
        elif help_md_file:
            log.info(f"Help -     from MD help file ['{help_md_file}']")
            help_file = help_md_file
            with open(help_md_file, "r") as file:
                help_content = file.read()

        # Help from options
        else:
            help_content = help_generation(
                arguments_dict=arguments_dict, setup=setup_cfg, output_type="markdown"
            )

        # Clean MD content
        help_content = help_content.replace(os.path.expanduser("~"), "~")

        with TemporaryDirectory() as tmp_dir:
            # help md tmp
            help_file_tmp = os.path.join(tmp_dir, "tmp.md")

            # Write file
            if help_md_file:
                toc = False
                shift_heading_level_by = 0
            else:
                toc = True
                shift_heading_level_by = -1

            # Write file
            f = open(help_file_tmp, "w")
            f.write(help_content)
            f.close()

            # # Generate Table Of Content (if marker <!--TOC-->)
            # if not help_md_file:
            #     import md_toc

            #     toc = md_toc.api.build_toc(help_file_tmp)
            #     md_toc.api.write_string_on_file_between_markers(
            #         help_file_tmp, toc, "<!--TOC-->"
            #     )

            # shutil.copy(help_file_tmp, help_file + ".tmp.md")

            # Generate MarkDown
            log.debug(f"help_file_tmp={help_file_tmp}")
            if help_file:
                pdoc_args = [
                    "-s",
                    f"--table-of-contents={toc}",
                    "-N",
                    f"--shift-heading-level-by={shift_heading_level_by}",
                    # "--fenced_code_attributes=true",
                ]
                pypandoc.convert_file(
                    help_file_tmp,
                    to="gfm",
                    format="markdown-smart+fenced_code_attributes",
                    outputfile=help_file,
                    extra_args=pdoc_args,
                )

            if "help_html" in args and args.help_html:
                help_file_html = full_path(args.help_html)
                log.info(f"Help - generate Markdown help file ['{help_file_html}']")
                if not os.path.exists(os.path.dirname(help_file_html)):
                    Path(os.path.dirname(help_file_html)).mkdir(
                        parents=True, exist_ok=True
                    )

                # Tmp md file
                help_file_tmp_html = help_file_tmp + ".html.md"
                with open(help_file_tmp, "r", encoding="utf-8") as f:
                    content = f.read()
                content_changed = re.sub(r"\.md([\)\"#])", r".html\1", content)
                with open(help_file_tmp_html, "w", encoding="utf-8") as fichier:
                    fichier.write(content_changed)

                # pdoc args
                pdoc_args = [
                    "-s",
                    f"--table-of-contents={toc}",
                    "-V",
                    "toc-title:Contents",
                    "-N",
                    f"--shift-heading-level-by={shift_heading_level_by}",
                    "--ascii=true",
                    "--metadata",
                    f"title={help_json_input_title}",
                    "-V",
                    "maxwidth:1000px",
                ]

                # Generate doc
                pypandoc.convert_file(
                    help_file_tmp_html,
                    to="html5",
                    format="markdown",
                    outputfile=help_file_html,
                    extra_args=pdoc_args,
                )
                # command = f"pandoc -s --from=markdown-smart --table-of-contents=true --to=html5 --metadata title='{help_json_input_title}' -o {help_file_html} {help_file}"
                # try:
                #     run_parallel_commands([command])
                # except Exception as inst:
                #     log.error(
                #         "Python pandoc package need to be installed ('pip install pandoc')"
                #     )
                #     log.error(inst)

            if "help_pdf" in args and args.help_pdf:
                help_file_pdf = full_path(args.help_pdf)
                log.info(f"Help - generate Markdown help file ['{help_file_pdf}']")
                if not os.path.exists(os.path.dirname(help_file_pdf)):
                    Path(os.path.dirname(help_file_pdf)).mkdir(
                        parents=True, exist_ok=True
                    )

                # Tmp md file
                help_file_tmp_pdf = help_file_tmp + ".pdf.md"
                with open(help_file_tmp, "r", encoding="utf-8") as f:
                    content = f.read()
                content_changed = re.sub(r"\.md([\)\"#])", r".pdf\1", content)
                with open(help_file_tmp_pdf, "w", encoding="utf-8") as fichier:
                    fichier.write(content_changed)

                pdoc_args = [
                    "-s",
                    f"--table-of-contents={toc}",
                    "-V",
                    "toc-title:Contents",
                    "-N",
                    "-V",
                    "geometry:margin=2cm",
                    "-V",
                    "fontsize:10pt",
                    f"--shift-heading-level-by={shift_heading_level_by}",
                    "--metadata",
                    f"title={help_json_input_title}",
                ]

                pypandoc.convert_file(
                    help_file_tmp_pdf,
                    to="pdf",
                    format="markdown-smart",
                    outputfile=help_file_pdf,
                    extra_args=pdoc_args,
                )

                # command = f"pandoc -s --from=markdown-smart --table-of-contents=true -V geometry:margin=3cm --to=pdf --metadata title='{help_json_input_title}' -o {help_file_pdf} {help_file_tmp}"
                # try:
                #     run_parallel_commands([command])
                # except Exception as inst:
                #     log.error(
                #         "Python pandoc package need to be installed ('pip install pandoc')"
                #     )
                #     log.error(inst)

            # Remove tmp
            # remove_if_exists(filepaths=[help_file_tmp])

            # if "help_html" in args and args.help_html:
            #     help_file_html = full_path(args.help_html)
            #     log.info(f"Help - generate Markdown help file ['{help_file_html}']")
            #     if not os.path.exists(os.path.dirname(help_file_html)):
            #         Path(os.path.dirname(help_file_html)).mkdir(parents=True, exist_ok=True)
            #     command = f"pandoc -s --from=markdown-smart --table-of-contents=true --to=html5 --metadata title='{help_json_input_title}' -o {help_file_html} {help_file}"
            #     try:
            #         run_parallel_commands([command])
            #     except Exception as inst:
            #         log.error(
            #             "Python pandoc package need to be installed ('pip install pandoc')"
            #         )
            #         log.error(inst)

            # if "help_pdf" in args and args.help_pdf:
            #     help_file_pdf = full_path(args.help_pdf)
            #     log.info(f"Help - generate Markdown help file ['{help_file_pdf}']")
            #     if not os.path.exists(os.path.dirname(help_file_pdf)):
            #         Path(os.path.dirname(help_file_pdf)).mkdir(parents=True, exist_ok=True)
            #     command = f"pandoc -s --from=markdown-smart --table-of-contents=true -V geometry:margin=3cm --to=pdf --metadata title='{help_json_input_title}' -o {help_file_pdf} {help_file}"
            #     try:
            #         run_parallel_commands([command])
            #     except Exception as inst:
            #         log.error(
            #             "Python pandoc package need to be installed ('pip install pandoc')"
            #         )
            #         log.error(inst)

            # # Generate Table Of Content (if marker <!--TOC-->)
            # if not help_md_file:
            #     import md_toc

            #     toc = md_toc.api.build_toc(help_file)
            #     md_toc.api.write_string_on_file_between_markers(
            #         help_file, toc, "<!--TOC-->"
            #     )

            # # PYPANDOC
            # import pypandoc

            # # pdoc_args = ["--mathjax", "+smart"]
            # pdoc_args = []

            # # With an input file: it will infer the input format from the filename
            # output = pypandoc.convert_file(
            #     help_file, "pdf", outputfile="docs/pdf/pypandoc.pdf", extra_args=pdoc_args
            # )
            # output = pypandoc.convert_file(
            #     help_file,
            #     "html",
            #     outputfile="docs/html/pypandoc.html",
            #     extra_args=pdoc_args,
            # )

    log.debug(f"main_folder={main_folder}")

    log.info("End")
