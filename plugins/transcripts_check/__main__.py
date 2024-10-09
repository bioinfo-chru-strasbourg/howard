import argparse
import json
import logging as log
import os
from pathlib import Path

from howard.functions.commons import load_args, load_config_args, transcripts_file_to_df
from howard.objects.variants import Variants
from howard.tools.tools import PathType

# Arguments
arguments = {
    "transcripts_expected": {
        "metavar": "List of transcripts (file)",
        "help": """File with a list of transcripts in first column.\n""",
        "required": False,
        "default": None,
        "type": PathType(exists=True, type="file"),
        "gooey": {
            "widget": "FileChooser",
            "options": {
                "wildcard": "TSV file (*.tsv)|*.tsv|" "All files (*)|*",
            },
        },
    },
    "transcripts_missing": {
        "metavar": "List of missing transcripts (file)",
        "help": """File with a list of missing transcripts in first column.\n""",
        "required": False,
        "default": None,
        "type": PathType(exists=None, type=None),
        "gooey": {
            "widget": "FileSaver",
            "options": {
                "wildcard": "TSV file (*.tsv)|*.tsv|" "All files (*)|*",
            },
        },
    },
}

# Command
commands_arguments = {
    "transcripts_check": {
        "function": "transcripts_check",
        "description": """Check if a transcript list is present in a generated transcript table from a input VCF file.\n""",
        "help": """Check transcript list in transcript table""",
        "epilog": """Usage examples:\n"""
        """   howard transcripts_check --input=plugins/transcripts_check/tests/data/example.ann.transcripts.vcf.gz --param=plugins/transcripts_check/tests/data/param.transcripts.json --transcripts_expected=plugins/transcripts_check/tests/data/transcripts.tsv --stats=/tmp/transcripts.stats.json --transcripts_missing=/tmp/transcripts.missing.tsv\n"""
        """    \n""",
        "groups": {
            "main": {
                "input": True,
                "param": True,
                "transcripts_expected": True,
                "transcripts_missing": False,
                "stats_json": False,
            }
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
def check_transcript(
    vcfdata_obj: Variants, transcripts_file: str, transcripts_table: str
) -> dict:

    # Init
    stats = {}

    # Transcript list dataframe
    transcripts_list = transcripts_file_to_df(
        transcripts_file=transcripts_file, column_names=["transcript"]
    )

    # Query and dataframe of existing transcripts
    check_query_transcripts = f"""
        SELECT string_split(transcripts.transcript, '.')[1] AS 'transcript_id', string_split(transcripts.transcript, '.')[2] AS 'transcript_ver', 
        FROM {transcripts_table} as transcripts
        GROUP BY transcript_id, transcript_ver
    """
    df_transcripts = vcfdata_obj.get_query_to_df(query=check_query_transcripts)

    # Query and dataframe of expecting transcripts
    check_query_transcripts_list = """
        SELECT string_split(transcripts_list.transcript, '.')[1] AS 'transcript_id', string_split(transcripts_list.transcript, '.')[2] AS 'transcript_ver', 
        FROM transcripts_list
        GROUP BY transcript_id, transcript_ver
    """
    df_transcripts_list = vcfdata_obj.get_query_to_df(
        query=check_query_transcripts_list
    )

    # Intersection
    check_query_intersection = f"""
        (
            {check_query_transcripts}
        )
        INTERSECT
        (
            {check_query_transcripts_list}
        )
    """
    df_intersection = vcfdata_obj.get_query_to_df(query=check_query_intersection)

    # Union
    check_query_union = f"""
        (
            {check_query_transcripts}
        )
        UNION
        (
            {check_query_transcripts_list}
        )
    """
    df_union = vcfdata_obj.get_query_to_df(query=check_query_union)

    # Missing
    check_query_missing = f"""
        (
            {check_query_transcripts_list}
        )
        EXCEPT
        (
            {check_query_transcripts}
        )
    """
    df_missing = vcfdata_obj.get_query_to_df(query=check_query_missing)

    # Missing list
    check_query_missing_list = f"""
        SELECT concat(
            transcript_id,
            CASE
                WHEN transcript_ver IS NOT NULL
                THEN concat ('.', transcript_ver)
                ELSE ''
            END
            ) AS 'transcript'
        FROM
        (
            (
                {check_query_transcripts_list}
            )
            EXCEPT
            (
                {check_query_transcripts}
            )
        )
        GROUP BY transcript
    """
    df_missing_list = vcfdata_obj.get_query_to_df(query=check_query_missing_list)

    # Stats
    stats = {
        "available": len(df_transcripts),
        "list": len(df_transcripts_list),
        "intersection": len(df_intersection),
        "union": len(df_union),
        "percent": len(df_intersection) / len(df_transcripts_list),
        "missing": len(df_missing),
        "missing_list": list(df_missing_list["transcript"]),
    }

    # Log
    log.debug(stats)

    return stats


# Main function
def main(args: argparse) -> None:

    log.info("START")

    # Load config args
    arguments_dict, _, config, param = load_config_args(args)

    # Create variants object
    vcfdata_obj = Variants(input=args.input, config=config, param=param)

    # Get Config and Params
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
    vcfdata_obj.set_param(param)

    # Load data
    vcfdata_obj.load_data()

    # Generate transcript table
    log.info("Create transcript table...")
    transcripts_table = vcfdata_obj.create_transcript_view()
    log.info(f"Create transcript table: '{transcripts_table}'")

    # Check transcripts table
    log.info("Generate stats...")
    stats = check_transcript(
        vcfdata_obj=vcfdata_obj,
        transcripts_file=args.transcripts_expected,
        transcripts_table=transcripts_table,
    )

    # Transcript missing file
    transcripts_missing = args.transcripts_missing

    # Create transcript missing object
    transcripts_missing_list = "\n".join(stats.get("missing_list", [])) + "\n"

    # Stats JSON file
    stats_json = args.stats_json

    # Create json object
    json_object = json.dumps(stats, indent=4)

    # Transcripts missing
    if args.transcripts_missing:

        # Log
        log.info(f"Export transcripts missing : '{transcripts_missing}'")

        # Create folder if not exists
        if not os.path.exists(os.path.dirname(transcripts_missing)):
            Path(os.path.dirname(transcripts_missing)).mkdir(
                parents=True, exist_ok=True
            )

        # Write output file
        with open(transcripts_missing, "w") as outfile:
            outfile.write(transcripts_missing_list)

    # Stats in JSON (or stdout)
    if args.stats_json:

        # Log
        log.info(f"Export stats in JSON format: '{stats_json}'")

        # Create folder if not exists
        if not os.path.exists(os.path.dirname(stats_json)):
            Path(os.path.dirname(stats_json)).mkdir(parents=True, exist_ok=True)

        # Write output file
        with open(stats_json, "w") as outfile:
            outfile.write(json_object)

    else:

        # Log
        log.info(f"Stats:")

        # Print stats
        print(json_object)

    log.info("END")
