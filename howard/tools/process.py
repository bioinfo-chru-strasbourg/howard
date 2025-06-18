import argparse
import logging as log
from tabulate import tabulate  # type: ignore

from howard.functions.commons import load_args, load_config_args
from howard.objects.variants import Variants

# from howard.functions.commons import *
# from howard.functions.databases import *


def process(args: argparse) -> None:
    """
    The "process" function processes input arguments, loads parameters in JSON format, creates a VCF
    object, performs quick annotations, calculations, prioritizations, and queries, exports output, and
    closes the connection.

    :param args: args is a variable that contains the arguments passed to the function "process". It is
    assumed to be an object with several attributes, including "config", "param", "input", "output",
    "annotations", "calculations", "prioritizations", and "query". These attributes are used to
    :type args: argparse
    """

    log.info("Start")

    # Load config args
    arguments_dict, _, config, param = load_config_args(args)

    # Create variants object
    vcfdata_obj = Variants(
        input=args.input, output=args.output, config=config, param=param
    )

    # Get Config and Params
    config = vcfdata_obj.get_config()
    param = vcfdata_obj.get_param()

    # Load args into param
    param = load_args(
        param=param,
        args=args,
        arguments_dict=arguments_dict,
        command="process",
        strict=False,
    )

    # Re-Load Config and Params
    vcfdata_obj.set_param(param)
    vcfdata_obj.set_config(config)

    # Load data
    vcfdata_obj.load_data()

    # Annotation HGVS
    vcfdata_obj.annotation_hgvs()

    # Annotations
    vcfdata_obj.annotation()

    # Calculations
    vcfdata_obj.calculation()

    # Prioritization
    vcfdata_obj.prioritization()

    # Explode infos
    if param.get("explode_infos", {}).get("explode_infos", False):
        vcfdata_obj.explode_infos()

    # Query
    if param.get("query", {}).get("query", None):

        log.info("Querying...")

        # Parameters
        query = param.get("query", {}).get("query", None)
        query_limit = param.get("query", {}).get("query_limit", None)
        query_print_mode = param.get("query", {}).get("query_print_mode", None)

        # Print query
        if query_print_mode in ["markdown"]:
            print(vcfdata_obj.get_query_to_df(query, limit=query_limit).to_markdown())
        elif query_print_mode in ["tabulate"]:
            print(
                tabulate(
                    vcfdata_obj.get_query_to_df(query, limit=query_limit),
                    headers="keys",
                    tablefmt="psql",
                )
            )
        else:
            print(vcfdata_obj.get_query_to_df(query, limit=query_limit))

    # Export
    vcfdata_obj.export_output(query=param.get("query", {}).get("query", None))
    # vcfdata_obj.export_output()

    # Log
    log.info("End")

    # Return Variants object
    return vcfdata_obj
