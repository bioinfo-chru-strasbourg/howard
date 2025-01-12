import argparse
import logging as log
from pathlib import Path
import os
import pandas as pd

from howard.functions.commons import load_args, load_config_args, full_path
from howard.objects.variants import Variants

# Arguments
arguments = {
    "add_variants_view": {
        "help": """Create a sheet with all INFO fields exploded.\n""",
        "action": "store_true",
        "default": False,
    },
    "add_header": {
        "help": """Create a sheet with all INFO fields header descritions.\n""",
        "action": "store_true",
        "default": False,
    },
}

# Command
commands_arguments = {
    "to_excel": {
        "function": "to_excel",
        "description": """Convert VCF file to Excel '.xlsx' format.\n""",
        "help": """Convert VCF file to Excel '.xlsx' format""",
        "epilog": """Usage examples:\n"""
        """   howard to_excel --input=tests/data/example.vcf.gz --output=/tmp/example.xlsx --add_variants_view\n"""
        """    \n""",
        "groups": {
            "main": {"input": True, "output": True},
            "Add": {"add_variants_view": False, "add_header": False},
        },
    }
}


# Function to write dataframe by batch into a Excel file as writer
def write_to_excel_in_batches(
    vcfdata_obj, query, table, writer, batch_size: int = None
):
    """
    Function to write dataframe by batch into a Excel file as writer
    """

    if batch_size is not None:

        offset = 0
        while True:
            # Read batch
            query_with_limit = f"{query} LIMIT {batch_size} OFFSET {offset}"
            dataframe = vcfdata_obj.get_query_to_df(query_with_limit)

            if dataframe.empty:
                break

            # Write batch
            if offset == 0:
                dataframe.to_excel(writer, sheet_name=table, index=False)
            else:
                dataframe.to_excel(
                    writer,
                    sheet_name=table,
                    index=False,
                    header=False,
                    startrow=offset + 1,
                )

            offset += batch_size

    else:

        # write dataframe
        dataframe = vcfdata_obj.get_query_to_df(query)
        dataframe.to_excel(writer, sheet_name=table, index=False)


# To Excel
def to_excel(
    vcfdata_obj: Variants, output: str, tables: list = None, batch_size: int = None
) -> list:
    """
    Generate Excel file
    """

    # Output fiel and folder
    output = full_path(output)
    output_folder = os.path.dirname(output)

    # Create folder if not  exists
    if not os.path.exists(output_folder):
        Path(output_folder).mkdir(parents=True, exist_ok=True)

    # Find all tables
    if tables is None:
        tables = list(vcfdata_obj.get_query_to_df(f"SHOW tables")["name"])

    # Write all tables into a Excel file
    with pd.ExcelWriter(output) as writer:

        # for each table
        for table in tables:

            # Query
            query = f"SELECT * FROM {table}"

            # Log
            log.info(f"Add sheet '{table}'")

            # Batch
            write_to_excel_in_batches(
                vcfdata_obj, query, table, writer, batch_size=batch_size
            )

    return tables


# Main function
def main(args: argparse) -> None:

    log.info("START")

    log.debug(f"Input file: {args.input}")

    # Load config args
    arguments_dict, _, config, param = load_config_args(args)

    # Create variants object
    vcfdata_obj = Variants(
        input=args.input, output=args.output, config=config, param=param
    )

    # Get Config and Params
    # config = vcfdata_obj.get_config()
    param = vcfdata_obj.get_param()
    log.debug(f"param={param}")

    # Load args into param
    param = load_args(
        param=param,
        args=args,
        arguments_dict=arguments_dict,
        command="to_excel",
        strict=False,
    )

    log.debug(f"param={param}")

    # Re-Load Config and Params
    # vcfdata_obj.set_config(config)
    vcfdata_obj.set_param(param)

    # Load data
    vcfdata_obj.load_data()

    if param.get("add", {}).get("add_header", False):

        # Log
        log.info("Add VCF header")

        # Load header
        try:
            vcfdata_obj.load_header()
        except:
            log.debug("error in header table creation")

    if param.get("add", {}).get("add_variants_view", False):

        # Log
        log.info("Add Variants view")

        # Table variants

        table_variants = vcfdata_obj.get_table_variants()

        try:
            # Create variants view
            vcfdata_obj.create_annotations_view(
                table=table_variants,
                fields_needed_all=True,
                detect_type_list=False,
                fields_forced_as_varchar=True,
                info_prefix_column="",
            )
        except:
            log.debug("error in annotations view creation")

    # Output
    output = vcfdata_obj.get_output()

    # Batch size
    # chunk_size = config.get("chunk_size", None)
    chunk_size = None

    # Convert
    sheets = to_excel(
        vcfdata_obj=vcfdata_obj, output=vcfdata_obj.get_output(), batch_size=chunk_size
    )

    # Lof
    log.debug(f"Excel file '{os.path.basename(output)}' created, with sheets {sheets}")

    log.info("END")

    return vcfdata_obj
