import argparse
import logging as log

from howard.functions.commons import load_param_and_config
from howard.objects.variants import Variants
from howard.tools.process import process


def annotation(args: argparse) -> None:
    """
    The `annotation` function performs annotation on a VCF file based on specified parameters and
    exports the annotated data.

    :param args: The `args` parameter is likely an object or dictionary containing various arguments
    passed to the `annotation` function. It is not clear from the code snippet what specific arguments
    are expected or required
    :type args: argparse
    """

    # Load args, param, config and vcfdata_obj
    #arguments_dict, config, param, vcfdata_obj = load_param_and_config(args=args, command="annotation", strict=False, load_data=False)

    vcfdata_obj = process(args=args, tools=["annotation"])

    return vcfdata_obj
