import argparse
import logging as log
from tabulate import tabulate  # type: ignore

from howard.functions.commons import load_args, load_config_args
from howard.objects.variants import Variants


def sort(args: argparse) -> None:
    """
    This Python function loads and sort variants from a VCF file based on user input and exports the
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

    # Access
    input_format = vcfdata_obj.get_input_format()
    if param.get("explode", {}).get("explode_infos", False) or not input_format in [
        "duckdb",
        "parquet",
    ]:
        access = "RW"
    else:
        access = "RO"
    config["access"] = access

    # Load args into param
    param = load_args(
        param=param,
        args=args,
        arguments_dict=arguments_dict,
        command="filter",
        strict=False,
    )

    # Load data
    if vcfdata_obj.get_input():
        vcfdata_obj.load_data()
        # vcfdata_obj.load_header()
        # view_name = "variants_view"
        # vcfdata_obj.create_annotations_view(
        #     view=view_name,
        #     view_type="view",
        #     view_mode="explore",
        #     info_prefix_column="",
        #     fields_needed_all=True,
        #     info_struct_column="INFOS",
        #     sample_struct_column="SAMPLES",
        #     detect_type_list=True,
        # )

    # Filtering
    log.info("Sorting...")

    # # Filter
    # filter = param.get("filters", {}).get("filter", None)

    # # Columns
    # columns = vcfdata_obj.get_header_columns_as_list()

    # # Samples
    # samples_param = param.get("filters", {}).get("samples", None)
    # samples = []
    # if not (samples_param is None or samples_param.strip() == ""):

    #     # Check samples in file
    #     samples_in_file = vcfdata_obj.get_header_sample_list(check=True)

    #     for s in samples_param.split(","):
    #         # Check if sample in file
    #         if s.strip() in samples_in_file:
    #             samples.append(s.strip())
    #         else:
    #             log.warning(f"Sample '{s.strip()}' not in file")

    #     if len(samples):
    #         # Remove samples from columns if not selected
    #         for s in samples_in_file:
    #             if s not in samples:
    #                 columns.remove(s)

    # Sort contigs
    vcfdata_obj.sort_contigs()

    # variants table
    table_variants = vcfdata_obj.get_table_variants()

    # Create case clause
    case_clause = ""
    for i, chrom in enumerate(vcfdata_obj.get_header().contigs):
        case_clause += f"""    WHEN "#CHROM" = '{chrom}' THEN {i + 1}\n"""

    # Create case clause order by
    if case_clause != "":
        case_clause_order_by = f"""
            ORDER BY 
            CASE
                {case_clause}
            END
        """
    else:
        case_clause_order_by = ""

    # Create sort query
    query_sort = f"""
        SELECT *
        FROM {table_variants}
        {case_clause_order_by}
    """

    # Export
    vcfdata_obj.export_output(query=query_sort, export_header=True)

    # Log
    log.info("End")

    # Return variants object
    return vcfdata_obj
