import argparse
import logging as log

from howard.functions.commons import load_param_and_config
from howard.tools.process import process


def prioritization(args: argparse) -> None:
    """
    The function performs prioritization on a VCF file based on user-specified configurations and
    exports the results.

    :param args: args is an object that contains the command line arguments passed to the script. It is
    used to configure the behavior of the script and to provide input and output file paths, as well as
    other parameters needed for the execution of the script
    :type args: argparse
    """

    # Load args, param, config and vcfdata_obj
    #arguments_dict, config, param, vcfdata_obj = load_param_and_config(args=args, command="annotation", strict=False, load_data=False)

    vcfdata_obj = process(args=args, tools=["prioritization"])

    return vcfdata_obj
