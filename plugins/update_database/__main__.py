import argparse
import logging as log
import sys
import os
from howard.functions.commons import DEFAULT_DATABASE_FOLDER
import multiprocess as mp
import concurrent.futures

sys.path.append(os.path.join(os.path.dirname(__file__)))
from plugins.update_database import clinvar, gnomad, cadd, omim
from plugins.update_database import utils

arguments = {
    "databases_folder": {
        "help": """Path of HOWARD database folder.\n""",
        "type": str,
        "default": DEFAULT_DATABASE_FOLDER,
    },
    "database": {
        "help": """Which database to update.\n""",
        "type": str,
        "choices": ["clinvar", "gnomad", "CADD", "omim"],
    },
    "data_folder": {
        "help": """Path of data needed to update database.\n""",
        "type": str,
    },
    "update_config": {
        "help": """Path of json configuration file.\n""",
        "default": os.path.join(os.path.dirname(__file__), "config", "update_databases.json"),
        "type": str,
    },
    "current_folder": {
        "help": """Path of json configuration file.\n""",
        "type": str,
        "default": "current",
    },
    "refseq": {
        "help": """Path of refseq file.\n""",
        "type": str,
    },
}

# Command
commands_arguments = {
    "update_database": {
        "function": "update_database",
        "description": """Update HOWARD database\n""",
        "help": """Update HOWARD database""",
        "epilog": """Usage examples:\n"""
        """   howard update_database --database clinvar --databases_folder /home1/DB/HOWARD --update_config update_databases.json  \n"""
        """    \n""",
        "groups": {
            "main": {"param": False},
            "Update_database": {
                "databases_folder": False,
                "database": False,
                "data_folder": False,
                "update_config": False,
                "current_folder": False,
                "refseq": False,
            },
            "Options": {"show": False, "limit": False},
        },
    }
}

# from plugins.update_database.ucsc import Ucsc


# Main function
def main(args: argparse) -> None:
    """hist | grep
    Query input VCF file and show result
    """

    # Log
    log.info("START")
    if args.database == "clinvar":
        log.info("Update Clinvar")
        clinvar.Clinvar(
            database=args.database,
            databases_folder=args.databases_folder,
            config_json=args.update_config,
            current_folder=args.current_folder,
        ).update_clinvar()

    elif args.database == "gnomad":
        log.info("Update Gnomad")
        gnom = gnomad.Gnomad(
            database=args.database,
            databases_folder=args.databases_folder,
            config_json=args.update_config,
            current_folder=args.current_folder,
            data_folder=args.data_folder,
        )

        with concurrent.futures.ProcessPoolExecutor(max_workers=6) as executor:
            futures = {
                executor.submit(gnom, file): os.path.join(gnom.data_folder, file)
                for file in os.listdir(gnom.data_folder)
                if file.endswith(".bgz")
                and not os.path.exists(
                    os.path.join(gnom.data_folder, file).replace(".vcf.bgz", ".parsed.vcf.gz")
                )
            }
            for future in concurrent.futures.as_completed(futures):
                try:
                    result = future.result()
                    log.info(result)
                except Exception as e:
                    log.error(f"Error processing {futures[future]}: {e}")
        gnom.update_gnomad()

    elif args.database == "CADD":
        cadd_input = [
            os.path.join(args.data_folder, file)
            for file in os.listdir(args.data_folder)
            if file.endswith(".tsv.gz")
        ][0]
        if cadd_input:
            log.info(f"CADD input {cadd_input}")
            input_args = cadd.update_cadd(
                cadd_input,
                os.path.join(args.data_folder, "processing"),
                os.path.join(args.data_folder, f"CADD.generated.{utils.now()}.partition.parquet"),
            )
        # Start processing
        with mp.Pool(10) as p:
            result = p.starmap(cadd.create_vcf_chunks, input_args)
            for r in result:
                log.debug(f"Generated {r} files")

    elif args.database == "omim":
        o = omim.Omim(
            database=args.database,
            databases_folder=args.databases_folder,
            config_json=args.update_config,
            current_folder=args.current_folder,
            refseq=args.refseq,
            data_folder=args.data_folder,
        )
        o.update_omim()

    # Debug
    log.info("END")
