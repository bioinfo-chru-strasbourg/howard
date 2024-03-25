import argparse
import logging as log
import vcf
import pandas as pd

from howard.functions.commons import load_args, load_config_args, code_type_map
from howard.objects.variants import Variants

# # Arguments
arguments = {}


# Command
commands_arguments = {
    "genebe": {
        "function": "genebe",
        "description": """GeneBe annotation using REST API (see https://genebe.net/).\n""",
        "help": """GeneBe annotation using REST API""",
        "epilog": """Usage examples:\n"""
        """   howard genebe --input=tests/data/example.vcf.gz --output=/tmp/example.geenbe.vcf.gz\n"""
        """    \n""",
        "groups": {
            "main": {"input": True, "output": True, "param": False, "assembly": False},
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

# dtype and vcf header map
TYPES = {
    "int": "Integer",
    "int64": "Integer",
    "float": "Float",
    "float64": "Float",
    "object": "String",
}


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

    # Re-Load Config and Params
    # vcfdata_obj.set_config(config)
    vcfdata_obj.set_param(param)

    # Load data
    vcfdata_obj.load_data()

    # GeneBe
    import genebe as gnb

    # Create dataframe for GeneBe
    df = vcfdata_obj.get_query_to_df(
        'SELECT "#CHROM" as chr, POS as pos, REF as ref, ALT as alt FROM variants'
    )

    # Annotate with GeneBe
    log.info("GeneBe annotation...")
    annotated = gnb.annotate_dataframe_variants(
        df,
        use_ensembl=False,
        use_refseq=True,
        genome=args.assembly,
        flatten_consequences=True,
    )

    # Header pointer
    header = vcfdata_obj.get_header()

    query_concat = []

    # For each column/annotation
    for column in annotated.columns:

        # Only new columns
        if column not in ["chr", "pos", "ref", "alt"]:

            # Debug
            log.debug(f"GeneBe Annotation '{column}'")

            # Header type
            header_type = TYPES.get(str(annotated[column].dtype))

            # Add column in Header
            header.infos[column] = vcf.parser._Info(
                column,
                ".",
                header_type,
                f"GeneBe '{column}' annotation",
                "GeneBe",
                "Unknown",
                code_type_map[header_type],
            )

            # Create query concat to INFO
            query_concat.append(
                f"""CASE
                        WHEN annotated.{column} NOT NULL
                        THEN concat(';{column}=', annotated.{column})
                        ELSE ''
                    END"""
            )

    # Query Update
    query_update = f"""
        UPDATE variants
        SET INFO = regexp_replace(concat(INFO, {','.join(query_concat)}), '^;', '')
        FROM annotated
        WHERE variants."#CHROM"=annotated.chr
          AND variants.POS=annotated.pos
          AND variants.REF=annotated.ref
          AND variants.ALT=annotated.alt
        """
    log.debug(query_update)
    vcfdata_obj.execute_query(query_update)

    ff = vcfdata_obj.get_query_to_df("SELECT * FROM variants")
    log.debug(ff)

    # Export
    vcfdata_obj.export_output()

    log.info("END")
