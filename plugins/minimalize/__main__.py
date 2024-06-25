import argparse
import logging as log

from howard.functions.commons import load_args, load_config_args, vcf_required_columns
from howard.objects.variants import Variants

# Arguments
arguments = {
    "minimalize_info": {
        "help": """Minimalize INFO field (e.g. '.' value).\n""",
        "action": "store_true",
        "default": False,
    },
    "minimalize_id": {
        "help": """Minimalize ID field (e.g. '.' value).\n""",
        "action": "store_true",
        "default": False,
    },
    "minimalize_qual": {
        "help": """Minimalize QUAL field (e.g. '.' value).\n""",
        "action": "store_true",
        "default": False,
    },
    "minimalize_filter": {
        "help": """Minimalize FILTER field (e.g. '.' value).\n""",
        "action": "store_true",
        "default": False,
    },
    "minimalize_samples": {
        "help": """Minimalize samples to keep only genotypes (i.e. 'GT').\n""",
        "action": "store_true",
        "default": False,
    },
    "remove_samples": {
        "help": """Remove all samples to keep only variants.\n""",
        "action": "store_true",
        "default": False,
    },
}

# Command
commands_arguments = {
    "minimalize": {
        "function": "minimalize",
        "description": """Minimalize a VCF file consists in put missing value ('.') on INFO/Tags, ID, QUAL or FILTER fields. Options can also minimalize samples (keep only GT) or remove all samples. INFO/tags can by exploded before minimalize to keep tags into separated columns (useful for Parquet or TSV format to constitute a database).\n""",
        "help": """Minimalize a VCF file, such as removing INFO/Tags or samples""",
        "epilog": """Usage examples:\n"""
        """   howard minimalize --input=tests/data/example.vcf.gz --output=/tmp/example.minimal.vcf.gz --minimalize_info --minimalize_filter --minimalize_qual --minimalize_id --minimalize_samples\n"""
        """   howard minimalize --input=tests/data/example.vcf.gz --output=/tmp/example.minimal.tsv --remove_samples --explode_infos --minimalize_info\n"""
        """    \n""",
        "groups": {
            "main": {"input": True, "output": True, "param": False},
            "Minimalize": {
                "minimalize_info": False,
                "minimalize_id": False,
                "minimalize_qual": False,
                "minimalize_filter": False,
                "minimalize_samples": False,
                "remove_samples": False,
            },
            "Explode": {
                "explode_infos": False,
                "explode_infos_prefix": False,
                "explode_infos_fields": False,
            },
            "Export": {
                "include_header": False,
                "order_by": False,
                "parquet_partitions": False,
            },
        },
    }
}


# Missing value
# QUAL as NULL because integer format (?!?)
missing_value = {
    "INFO": "'.'",
    "ID": "'.'",
    "FILTER": "'.'",
    "QUAL": "NULL",
}


# Minimalize
def minimalize(vcfdata_obj: Variants, field: str):

    if field in ["samples"]:
        header_columns = vcfdata_obj.get_header_columns_as_list()
        update_query_set = []
        if "FORMAT" in header_columns:
            update_query_set.append("FORMAT='GT' ")
        for sample in header_columns:
            if sample not in vcf_required_columns + ["FORMAT"]:
                update_query_set.append(
                    f""" "{sample}"=string_split(CAST("{sample}" AS VARCHAR), ':')[1] """
                )
        update_query = "UPDATE variants SET " + ", ".join(update_query_set)
    else:
        update_query = f"UPDATE variants SET {field.upper()}={missing_value.get(field.upper(), None)}"
        log.info(f"Minimalization {field.upper()}...")
    vcfdata_obj.execute_query(update_query)


# Remove
def remove(vcfdata_obj: Variants, field: str):

    if field in ["samples"]:
        header_columns = vcfdata_obj.get_header_columns_as_list()
        for sample in header_columns:
            if sample not in vcf_required_columns:
                if sample != "FORMAT":
                    log.info(f"Removing '{sample}' sample...")
                vcfdata_obj.execute_query(f"ALTER TABLE variants DROP {sample};")


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

    # Load args into param
    param = load_args(
        param=param,
        args=args,
        arguments_dict=arguments_dict,
        command="minimalize",
        strict=False,
    )

    log.debug(f"param={param}")

    # Re-Load Config and Params
    # vcfdata_obj.set_config(config)
    vcfdata_obj.set_param(param)

    # Load data
    vcfdata_obj.load_data()

    param_minimalize = param.get("minimalize", {})
    for minimalize_action in param_minimalize:
        if param_minimalize.get(minimalize_action):
            minimalize_action_infos = minimalize_action.split("_")
            minimalize_action_type = minimalize_action_infos[0]
            minimalize_action_target = minimalize_action_infos[1]
            eval(
                f"{minimalize_action_type}(vcfdata_obj=vcfdata_obj,field='{minimalize_action_target}')"
            )

    # Export
    vcfdata_obj.export_output()

    log.info("END")
