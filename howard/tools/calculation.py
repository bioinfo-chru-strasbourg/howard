import argparse
import logging as log

from howard.functions.commons import load_param_and_config
#from howard.objects.variants import Variants
from howard.tools.process import process


# def load_param_and_config(args: argparse, load_data:bool = False) -> tuple:
#     """
#     The `load_param_and_config` function loads configuration and parameter settings for a given set of arguments.

#     :param args: The `args` parameter is an instance of the `argparse.Namespace` class, which contains the command-line
#     arguments passed to the script. It is used to access the values of the arguments and options specified by the user
#     when running the script
#     :type args: argparse.Namespace
#     :return: a tuple containing four elements: `arguments_dict`, `config`, `param`, and `vcfdata_obj`.
#     """

#     from howard.objects.variants import Variants

#     # Load config args
#     arguments_dict, _, config, param = load_config_args(args)

#     log.debug(f"param={param}")

#     # Create variants object
#     vcfdata_obj = Variants(
#         input=args.input, output=args.output, config=config, param=param
#     )

#     # Get Config and Params
#     config = vcfdata_obj.get_config()
#     param = vcfdata_obj.get_param()

#     log.debug(f"param={param}")

#     # Load args into param
#     param = load_args(
#         param=param,
#         args=args,
#         arguments_dict=arguments_dict,
#         command="calculation",
#         strict=False,
#     )

#     log.debug(f"param2={param}")

#     # Re-Load Config and Params
#     vcfdata_obj.set_param(param)
#     vcfdata_obj.set_config(config)

#     # Load data
#     if load_data and vcfdata_obj.get_input():
#         vcfdata_obj.load_data()

#     return arguments_dict, config, param, vcfdata_obj


def calculation(args: argparse) -> None:
    """
    This function performs calculations on VCF data based on user input and exports the results.

    :param args: The `args` parameter is a command line argument parser object that contains the
    arguments passed to the script when it was executed
    :type args: argparse
    """

    # vcfdata_obj = process(args=args, tools=["calculation"])

    # return vcfdata_obj

    # log.info("Start")

    # # Load config args
    # arguments_dict, _, config, param = load_config_args(args)

    # log.debug(f"param={param}")

    # # Create variants object
    # vcfdata_obj = Variants(
    #     input=args.input, output=args.output, config=config, param=param
    # )

    # # Get Config and Params
    # config = vcfdata_obj.get_config()
    # param = vcfdata_obj.get_param()

    # log.debug(f"param={param}")

    # # Load args into param
    # param = load_args(
    #     param=param,
    #     args=args,
    #     arguments_dict=arguments_dict,
    #     command="calculation",
    #     strict=False,
    # )

    # log.debug(f"param2={param}")

    # # Re-Load Config and Params
    # vcfdata_obj.set_param(param)
    # vcfdata_obj.set_config(config)

    # # Load data
    # if vcfdata_obj.get_input():
    #     vcfdata_obj.load_data()

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
    #vcfdata_obj.calculation(operations_config_file=operations_config_file)

    # # Export
    # if vcfdata_obj.get_input() or vcfdata_obj.get_output():
    #     vcfdata_obj.export_output()

    # # Log
    # log.info("End")

    # Return Variants object
    return vcfdata_obj
