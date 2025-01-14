import argparse
import logging as log

from howard.functions.commons import load_args, load_config_args
from howard.objects.variants import Variants


def convert(args: argparse) -> None:
    """
    The `convert` function converts a VCF file to a different format and can optionally explode info
    fields.

    :param args: `args` is a parameter passed to the `convert` function, likely an object or dictionary
    containing various arguments needed for the function to perform its task. These arguments could
    include things like input and output file paths, configuration settings, and other parameters
    :type args: argparse
    """

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
        command="convert",
        strict=False,
    )

    # Access
    if not param.get("explode", {}).get("explode_infos", False):
        config["access"] = "RO"

    # Re-Load Config and Params
    vcfdata_obj.set_param(param)
    vcfdata_obj.set_config(config)

    # Load data
    vcfdata_obj.load_data()

    # Export
    vcfdata_obj.export_output()

    # Log
    log.info("End")

    # Return Variants object
    return vcfdata_obj
