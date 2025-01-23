import argparse
import logging as log

from howard.functions.commons import load_args, load_config_args
from howard.objects.variants import Variants


def annotation(args: argparse) -> None:
    """
    The `annotation` function performs annotation on a VCF file based on specified parameters and
    exports the annotated data.

    :param args: The `args` parameter is likely an object or dictionary containing various arguments
    passed to the `annotation` function. It is not clear from the code snippet what specific arguments
    are expected or required
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
        command="annotation",
        strict=False,
    )

    # Re-Load Config and Params
    vcfdata_obj.set_param(param)
    vcfdata_obj.set_config(config)

    # Load data
    vcfdata_obj.load_data()

    # Annotation
    vcfdata_obj.annotation()

    # Export
    vcfdata_obj.export_output()

    # Log
    log.info("End")

    # Return Variants object
    return vcfdata_obj
