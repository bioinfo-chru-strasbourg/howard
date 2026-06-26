import argparse
import logging as log

from howard.functions.commons import load_param_and_config


def hgvs(args: argparse) -> None:
    """
    The `hgvs` function takes command line arguments, creates a VCF object, sets parameters and
    configurations, loads data from an input file, performs annotation using HGVS notation, exports the
    output, and closes the connection.

    :param args: The `args` parameter is of type `argparse.Namespace` and is used to parse command line
    arguments. It contains the following attributes:
    :type args: argparse
    """

    log.info("Start")

    # Load args, param, config and vcfdata_obj
    _, _, _, vcfdata_obj = load_param_and_config(args=args, command="hgvs", strict=False, load_data=True)

    # HGVS annotation
    vcfdata_obj.annotation_hgvs()

    # Export
    vcfdata_obj.export_output(export_header=True)

    # Log
    log.info("End")

    # Return Variants object
    return vcfdata_obj
