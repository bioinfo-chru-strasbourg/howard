import argparse
import logging as log

from howard.functions.commons import load_param_and_config


def convert(args: argparse) -> None:
    """
    The `convert` function converts a VCF file to a different format and can optionally explode info
    fields.

    :param args: `args` is a parameter passed to the `convert` function, likely an object or dictionary
    containing various arguments needed for the function to perform its task. These arguments could
    include things like input and output file paths, configuration settings, and other parameters
    :type args: argparse
    """

    # Log
    log.info("Start")

    # Load args, param, config and vcfdata_obj
    _, config, param, vcfdata_obj = load_param_and_config(args=args, command="convert", strict=False, load_data=False)

    # Access
    config["access"] = config.get("access", "RO") or "RO"

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
        # log.debug(f"Query for export: {query}")
        param["explode"]["explode_infos"] = False

    # Export
    vcfdata_obj.export_output(query=query)

    # Log
    log.info("End")

    # Return Variants object
    return vcfdata_obj
