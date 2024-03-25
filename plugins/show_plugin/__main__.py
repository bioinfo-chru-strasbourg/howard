import argparse
import logging as log
from howard.objects.variants import Variants

# Arguments
arguments = {
    "show": {
        "help": """show variants in input file.\n""",
        "action": "store_true",
        "default": False,
    },
    "limit": {
        "help": """Limit of output variants.\n""" """Either '2' or '5' lines.\n""",
        "default": 2,
        "type": int,
        "choices": [2, 5],
    },
}

# Command
commands_arguments = {
    "show_plugin": {
        "function": "show_plugin",
        "description": """Show variants in an input file.\n""",
        "help": """Short description of the plugin1""",
        "epilog": """Usage examples:\n"""
        """   howard show_plugin --input=tests/data/example.vcf.gz --output=/tmp/example.minimal.vcf.gz  --show --limit=5 \n"""
        """   howard show_plugin --input=tests/data/example.vcf.gz --output=/tmp/example.minimal.tsv  --show \n"""
        """    \n""",
        "groups": {
            "main": {"input": True, "output": False, "param": False},
            "Options": {"show": False, "limit": False},
        },
    }
}


# Main function
def main(args: argparse) -> None:
    """
    Query input VCF file and show result
    """

    # Log
    log.info("START")

    # Debug
    log.debug(f"Input file: {args.input}")
    log.debug(f"Output file: {args.output}")

    # Load variants file
    variants_obj = Variants(input=args.input, output=args.output, load=True)

    # Create query and show results
    query = f"SELECT * FROM variants LIMIT {args.limit}"
    if args.show:
        df = variants_obj.get_query_to_df(query)
        log.info(df)

    # Export
    variants_obj.export_output(query=query, export_header=True)

    log.info("END")
