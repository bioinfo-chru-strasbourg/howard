import argparse
import logging as log
from tabulate import tabulate  # type: ignore

from howard.functions.commons import load_args, load_config_args
from howard.objects.variants import Variants


def query(args: argparse) -> None:
    """
    This Python function loads and queries data from a VCF file based on user input and exports the
    results.

    :param args: args is an object that contains the arguments passed to the function. It is likely a
    Namespace object created by parsing command line arguments using argparse
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
        command="query",
        strict=False,
    )

    # Access
    if config.get("access", None) is None:
        input_format = vcfdata_obj.get_input_format()
        if param.get("explode", {}).get("explode_infos", False) or input_format not in [
            "duckdb",
            "parquet",
        ]:
            config["access"] = "RW"
        else:
            config["access"] = "RO"

    # Re-Load Config and Params
    vcfdata_obj.set_param(param)
    vcfdata_obj.set_config(config)

    # Load data
    if vcfdata_obj.get_input():
        vcfdata_obj.load_data()
        vcfdata_obj.load_header()
        vcfdata_obj.create_annotations_view(
            view="variants_view",
            view_type="view",
            view_mode="explore",
            info_prefix_column="",
            fields_needed_all=True,
            info_struct_column="INFOS",
            sample_struct_column="SAMPLES",
            detect_type_list=True,
        )

    # Query
    if param.get("query", {}).get("query", None):

        log.info("Querying...")

        # Explode
        if param.get("explode", {}).get("explode_infos", False):
            vcfdata_obj.explode_infos()

        # Parameters
        query = param.get("query", {}).get("query", None)
        query_limit = param.get("query", {}).get("query_limit", None)
        query_print_mode = param.get("query", {}).get("query_print_mode", None)

        # Print query
        if query_print_mode is not None:
            if query_print_mode.lower() in ["dataframe"]:
                print(vcfdata_obj.get_query_to_df(query, limit=query_limit))
            elif query_print_mode.lower() in ["markdown"]:
                print(
                    vcfdata_obj.get_query_to_df(query, limit=query_limit).to_markdown()
                )
            elif query_print_mode.lower() in ["tabulate"]:
                print(
                    tabulate(
                        vcfdata_obj.get_query_to_df(query, limit=query_limit),
                        headers="keys",
                        tablefmt="psql",
                    )
                )
            elif query_print_mode.lower() in ["no", "none", "null", "disabled"]:
                log.info("Query print mode disabled")
            else:
                print(vcfdata_obj.get_query_to_df(query, limit=query_limit))
        # if not output/export
        elif not vcfdata_obj.get_output():
            print(vcfdata_obj.get_query_to_df(query, limit=query_limit))

    # Export
    if vcfdata_obj.get_output():
        vcfdata_obj.export_output(
            query=param.get("query", {}).get("query", None), export_header=True
        )

    # Log
    log.info("End")

    # Return variants object
    return vcfdata_obj
