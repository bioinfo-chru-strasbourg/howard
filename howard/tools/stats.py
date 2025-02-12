import argparse
import logging as log

from howard.functions.commons import load_args, load_config_args
from howard.objects.variants import Variants


def stats(args: argparse) -> None:
    """
    The stats() function takes in arguments, loads data from an input file, gets statistics on the data,
    and closes the connection.

    :param args: args is a parameter that is passed to the function stats(). It is likely an object or a
    dictionary that contains various arguments or parameters that are needed by the function to perform
    its tasks. Some of the arguments that may be included in args are input file path, configuration
    settings, and other parameters that are
    :type args: argparse
    """

    log.info("Start")

    # Load config args
    arguments_dict, _, config, param = load_config_args(args)

    # Create variants object
    vcfdata_obj = Variants(input=args.input, config=config, param=param)

    # Get Config and Params
    config = vcfdata_obj.get_config()
    param = vcfdata_obj.get_param()

    # Load args into param
    param = load_args(
        param=param,
        args=args,
        arguments_dict=arguments_dict,
        command="stats",
        strict=False,
    )

    # Access
    config["access"] = "RO"

    # Re-Load Config and Params
    vcfdata_obj.set_param(param)
    vcfdata_obj.set_config(config)

    # Load data
    vcfdata_obj.load_data()

    # Parameters
    stats_stdout = param.get("stats", {}).get("stats_stdout", False)
    stats_md = param.get("stats", {}).get("stats_md", None)
    stats_json = param.get("stats", {}).get("stats_json", None)
    stats_html = param.get("stats", {}).get("stats_html", None)
    stats_pdf = param.get("stats", {}).get("stats_pdf", None)
    annotations_stats = param.get("stats", {}).get("annotations_stats", False)
    queries = param.get("stats", {}).get("queries", None)
    queries_view = param.get("stats", {}).get("queries_view", None)

    # Force if no stats file
    if (
        stats_md is None
        and stats_json is None
        and stats_html is None
        and stats_pdf is None
    ):
        stats_stdout = True

    # Stats
    vcfdata_obj.print_stats(
        stdout=stats_stdout,
        output_file=stats_md,
        json_file=stats_json,
        html_file=stats_html,
        pdf_file=stats_pdf,
        annotations_stats=annotations_stats,
        queries=queries,
        queries_view=queries_view,
    )

    # Log
    log.info("End")

    # Return Variants object
    return vcfdata_obj
