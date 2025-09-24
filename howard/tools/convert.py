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
    config["access"] = config.get("access", "RO")

    # Init
    param["explode"] = param.get("explode", {})

    # Re-Load Config and Params
    vcfdata_obj.set_param(param)
    vcfdata_obj.set_config(config)

    # Determine view type and mode (either "table" or "view", either "explore" or "full")
    view_type = "view"
    view_mode = "explore"

    # Output format
    output_format = vcfdata_obj.get_output_format()

    # Load data
    vcfdata_obj.load_data()

    # Explode Infos

    # Init
    query = None

    # If input format is vcf, explode_infos is set to False
    if output_format in ["vcf"]:
        param["explode"]["explode_infos"] = False

    # If explode infos is set to True, create annotation view and set query
    elif param.get("explode", {}).get("explode_infos", False):

        # Fields to explode
        fields = param.get("explode", {}).get("explode_infos_fields", None)
        if fields is not None:
            fields = fields.split(",")

        # Prefix
        info_prefix_column = param.get("explode", {}).get("explode_infos_prefix", "")

        # Create annotation view with infos from explode_infos_fields
        annotation_view_name = "variants_view_export"
        annotation_view_name = vcfdata_obj.create_annotations_view(
            view=annotation_view_name,
            view_type=view_type,
            view_mode=view_mode,
            info_prefix_column=info_prefix_column,
            fields=fields,
            fields_not_exists=False,
            fields_needed_all=True,
            fields_forced_as_varchar=True,
            # info_struct_column=None,
            # sample_struct_column=None,
            detect_type_list=False,
        )
        query = f"SELECT * FROM {annotation_view_name}"
        param["explode"]["explode_infos"] = False

    # Export
    vcfdata_obj.export_output(query=query)

    # Log
    log.info("End")

    # Return Variants object
    return vcfdata_obj
