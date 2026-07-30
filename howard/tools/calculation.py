import argparse
import logging as log

from howard.functions.commons import load_param_and_config
from howard.tools.process import process


def calculation(args: argparse) -> None:
    """
    This function performs calculations on VCF data based on user input and exports the results.

    :param args: The `args` parameter is a command line argument parser object that contains the
    arguments passed to the script when it was executed
    :type args: argparse
    """

    # Load args, param, config and vcfdata_obj
    _, _, param, vcfdata_obj = load_param_and_config(args=args, command="calculation", strict=False, load_data=False)


    # Show calculation
    if param.get("calculation", {}).get("show_calculations", False) or param.get(
        "calculation", {}
    ).get("show_calculations_md", None):

        # Operations config file
        operations_config_file = param.get("calculation", {}).get("calculation_config")

        # Show calculation in log
        if param.get("calculation", {}).get("show_calculations", False):
            for help_line in vcfdata_obj.get_operations_help(
                operations_config_file=operations_config_file
            ):
                log.info(help_line)

        # Show calculation in markdown
        if param.get("calculation", {}).get("show_calculations_md", None):
            for help_line in vcfdata_obj.get_operations_help(
                operations_config_file=operations_config_file, format="markdown"
            ):
                log.info(help_line)

        exit()


    # Calculation
    vcfdata_obj = process(args=args, tools=["calculation"])

    # Return Variants object
    return vcfdata_obj
