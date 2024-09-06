#!/usr/bin/env python

import argparse
import datetime
from functools import partial
import itertools
import multiprocessing
import os
import subprocess
import pyarrow.parquet as pq
import pyarrow as pa
from pyarrow import csv
import duckdb
import pandas as pd
import Bio.bgzf as bgzf
import numpy as np
import concurrent.futures
from multiprocessing import Pool, cpu_count
import dask.dataframe as dd
import logging as log
import fnmatch
import glob

import os
import requests
import shutil
import zipfile
import gzip
import pandas as pd
from typing import List
from tempfile import TemporaryDirectory

from jproperties import Properties  # jproperties 2.1.1


from howard.functions.commons import *
from howard.objects.variants import *
from howard.functions.databases import *
from howard.functions.from_annovar import *
from howard.functions.from_extann import from_extann


def databases(args: argparse) -> None:
    """
    The function downloads databases and logs the start and end of the process.

    :param args: The "args" parameter is likely an object or dictionary containing various arguments or
    options related to the "databases" function. Without more context, it's difficult to say exactly
    what these arguments might be, but they could include things like the names or locations of
    databases to download, authentication credentials, or
    :type args: argparse
    """

    log.info("Start")

    # Load config args
    arguments_dict, setup_cfg, config, param = load_config_args(args)

    # Load args into param
    param = load_args(
        param=param,
        args=args,
        arguments_dict=arguments_dict,
        command="databases",
        strict=False,
        section_prefix=["databases"],
    )

    # Assembly
    assemblies = param.get("databases", {}).get(
        "assemblies", param.get("databases", {}).get("assembly", args.assembly)
    )

    if assemblies:
        if isinstance(assemblies, str):
            assemblies = [value for value in assemblies.split(",")]
        param = add_value_into_dict(
            dict_tree=param, sections=["databases", "assemblies"], value=assemblies
        )

    # Genome folder
    genomes_folder = param.get("databases", {}).get("genomes_folder", None)
    if param.get("databases", {}).get("genomes_folder", None) and not isinstance(
        param.get("databases", {}).get("genomes_folder", None), str
    ):
        genomes_folder = param.get("databases", {}).get("genomes_folder", None).name
    if (
        "config" in args
        and args.config
        and args.config.get("folders", {}).get("databases", {}).get("genomes")
        and param.get("databases", {}).get("genomes_folder", None)
        == DEFAULT_REFSEQ_FOLDER
    ):
        genomes_folder = (
            args.config.get("folders", {}).get("databases", {}).get("genomes")
        )
    param = add_value_into_dict(
        dict_tree=param, sections=["databases", "genomes_folder"], value=genomes_folder
    )

    # Annovar files
    annovar_files = (
        param.get("databases", {})
        .get("annovar", {})
        .get("download_annovar_files", None)
    )
    if annovar_files and isinstance(annovar_files, str):
        annovar_files = [value for value in annovar_files.split(",")]
        param = add_value_into_dict(
            dict_tree=param,
            sections=["databases", "annovar", "download_annovar_files"],
            value=annovar_files,
        )

    # refSeq files
    refseq_files = (
        param.get("databases", {}).get("refseq", {}).get("download_refseq_files", None)
    )
    if refseq_files and isinstance(refseq_files, str):
        refseq_files = [value for value in refseq_files.split(",")]
        param = add_value_into_dict(
            dict_tree=param,
            sections=["databases", "refseq", "download_refseq_files"],
            value=refseq_files,
        )

    # Exomise Application Properties
    exomiser_application_properties = (
        param.get("databases", {})
        .get("exomiser", {})
        .get("download_exomiser_application_properties", None)
    )
    if exomiser_application_properties and not isinstance(
        exomiser_application_properties, str
    ):
        exomiser_application_properties = exomiser_application_properties.name
        param = add_value_into_dict(
            dict_tree=param,
            sections=[
                "databases",
                "exomiser",
                "download_exomiser_application_properties",
            ],
            value=exomiser_application_properties,
        )

    # dbSNP releases
    download_dbsnp_releases = (
        param.get("databases", {}).get("dbsnp", {}).get("download_dbsnp_releases", None)
    )

    if download_dbsnp_releases and isinstance(download_dbsnp_releases, str):
        download_dbsnp_releases = [
            value for value in download_dbsnp_releases.split(",")
        ]
        param = add_value_into_dict(
            dict_tree=param,
            sections=["databases", "dbsnp", "download_dbsnp_releases"],
            value=download_dbsnp_releases,
        )

    # Convert HGMD file
    convert_hgmd_file = (
        param.get("databases", {}).get("hgmd", {}).get("convert_hgmd_file", None)
    )
    if convert_hgmd_file and not isinstance(convert_hgmd_file, str):
        convert_hgmd_file = convert_hgmd_file.name
        param = add_value_into_dict(
            dict_tree=param,
            sections=[
                "databases",
                "hgmd",
                "convert_hgmd_file",
            ],
            value=convert_hgmd_file,
        )

    # Threads
    threads = get_threads(config=config, param=param)

    # Memory
    memory = extract_memory_in_go(get_memory(config=config, param=param))

    # Param
    if "generate_param" in args and args.generate_param:
        generate_databases_param(args=args, assemblies=assemblies)
        return None

    # Param database
    param_databases = param.get("databases", {})

    if not param_databases.get("assemblies", None):
        msg_error = "Not assemblies defined"
        log.error(msg_error)
        raise ValueError(msg_error)
    else:
        assemblies = param_databases.get("assemblies", [])

    # Log
    log.debug(f"param_databases={param_databases}")

    # Genomes
    if param_databases.get("genomes", {}).get("download_genomes", None):
        log.debug(f"Download Genomes")
        param_databases_genomes = param_databases.get("genomes", {})
        databases_download_genomes(
            assemblies=assemblies,
            genomes_folder=param_databases_genomes.get("download_genomes"),
            provider=param_databases_genomes.get("download_genomes_provider"),
            contig_regex=param_databases_genomes.get("download_genomes_contig_regex"),
            threads=threads,
        )

    # Annovar
    if param_databases.get("annovar", {}).get("download_annovar", None):
        log.debug(f"Download Annovar databases")
        param_databases_annovar = param_databases.get("annovar", {})
        databases_download_annovar(
            folder=param_databases_annovar.get("download_annovar"),
            files=param_databases_annovar.get("download_annovar_files"),
            assemblies=assemblies,
            annovar_url=param_databases_annovar.get("download_annovar_url"),
            threads=threads,
        )

    # snpEff
    if param_databases.get("snpeff", {}).get("download_snpeff", None):
        log.debug(f"Download snpEff databases")
        param_databases_snpeff = param_databases.get("snpeff", {})
        databases_download_snpeff(
            folder=param_databases_snpeff.get("download_snpeff"),
            assemblies=assemblies,
            config=config,
            threads=threads,
        )

    # refSeq
    if param_databases.get("refseq", {}).get("download_refseq", None):
        log.debug(f"Download refSeq databases")
        param_databases_refseq = param_databases.get("refseq", {})
        databases_download_refseq(
            assemblies=assemblies,
            refseq_folder=param_databases_refseq.get("download_refseq", None),
            refseq_url=param_databases_refseq.get("download_refseq_url", None),
            refseq_prefix=param_databases_refseq.get("download_refseq_prefix", None),
            refseq_files=param_databases_refseq.get("download_refseq_files", None),
            refseq_format_file=param_databases_refseq.get(
                "download_refseq_format_file", None
            ),
            include_utr_5=param_databases_refseq.get(
                "download_refseq_include_utr5", None
            ),
            include_utr_3=param_databases_refseq.get(
                "download_refseq_include_utr3", None
            ),
            include_chrM=param_databases_refseq.get(
                "download_refseq_include_chrM", None
            ),
            include_non_canonical_chr=param_databases_refseq.get(
                "download_refseq_include_non_canonical_chr", None
            ),
            include_non_coding_transcripts=param_databases_refseq.get(
                "download_refseq_include_non_coding_transcripts", None
            ),
            include_transcript_ver=param_databases_refseq.get(
                "download_refseq_include_transcript_version", None
            ),
            threads=threads,
            memory=memory,
        )

    # dbNSFP
    if param_databases.get("dbnsfp", {}).get("download_dbnsfp", None):
        log.debug(f"Download dbNSFP")
        param_databases_dbnsfp = param_databases.get("dbnsfp", {})
        databases_download_dbnsfp(
            assemblies=assemblies,
            dbnsfp_folder=param_databases_dbnsfp.get("download_dbnsfp", None),
            dbnsfp_url=param_databases_dbnsfp.get("download_dbnsfp_url", None),
            dbnsfp_release=param_databases_dbnsfp.get("download_dbnsfp_release", None),
            parquet_size=param_databases_dbnsfp.get(
                "download_dbnsfp_parquet_size", None
            ),
            generate_sub_databases=param_databases_dbnsfp.get(
                "download_dbnsfp_subdatabases", None
            ),
            generate_parquet_file=param_databases_dbnsfp.get(
                "download_dbnsfp_parquet", None
            ),
            generate_vcf_file=param_databases_dbnsfp.get("download_dbnsfp_vcf", None),
            not_generate_files_all=param_databases_dbnsfp.get(
                "download_dbnsfp_no_files_all", None
            ),
            add_info=param_databases_dbnsfp.get("download_dbnsfp_add_info", None),
            uniquify=param_databases_dbnsfp.get("download_dbnsfp_uniquify", False),
            row_group_size=param_databases_dbnsfp.get(
                "download_dbnsfp_row_group_size", None
            ),
            genomes_folder=param.get("databases", {}).get("genomes_folder", None),
            threads=threads,
            memory=memory,
        )

    # AlphaMissense
    if param_databases.get("alphamissense", {}).get("download_alphamissense", None):
        log.debug(f"Download AlphaMissense")
        param_databases_alphamissense = param_databases.get("alphamissense", {})
        databases_download_alphamissense(
            assemblies=assemblies,
            alphamissense_folder=param_databases_alphamissense.get(
                "download_alphamissense", None
            ),
            alphamissense_url=param_databases_alphamissense.get(
                "download_alphamissense_url", None
            ),
            threads=threads,
        )

    # Exomiser
    if param_databases.get("exomiser", {}).get("download_exomiser", None):
        log.debug(f"Download Exomiser")
        param_databases_exomiser = param_databases.get("exomiser", {})
        databases_download_exomiser(
            assemblies=assemblies,
            exomiser_folder=param_databases_exomiser.get("download_exomiser", None),
            exomiser_application_properties=param_databases_exomiser.get(
                "download_exomiser_application_properties", None
            ),
            exomiser_url=param_databases_exomiser.get("download_exomiser_url", None),
            exomiser_release=param_databases_exomiser.get(
                "download_exomiser_release", None
            ),
            exomiser_phenotype_release=param_databases_exomiser.get(
                "download_exomiser_phenotype_release", None
            ),
            exomiser_remm_release=param_databases_exomiser.get(
                "download_exomiser_remm_release", None
            ),
            exomiser_remm_url=param_databases_exomiser.get(
                "download_exomiser_remm_url", None
            ),
            exomiser_cadd_release=param_databases_exomiser.get(
                "download_exomiser_cadd_release", None
            ),
            exomiser_cadd_url=param_databases_exomiser.get(
                "download_exomiser_cadd_url", None
            ),
            exomiser_cadd_url_snv_file=param_databases_exomiser.get(
                "download_exomiser_cadd_url_snv_file", None
            ),
            exomiser_cadd_url_indel_file=param_databases_exomiser.get(
                "download_exomiser_cadd_url_indel_file", None
            ),
            threads=threads,
        )

    # dbSNP
    if param_databases.get("dbsnp", {}).get("download_dbsnp", None):
        log.debug(f"Download dbSNP")
        param_databases_dbsnp = param_databases.get("dbsnp", {})
        databases_download_dbsnp(
            assemblies=assemblies,
            dbsnp_folder=param_databases_dbsnp.get("download_dbsnp", None),
            dbsnp_releases=param_databases_dbsnp.get("download_dbsnp_releases", None),
            dbsnp_release_default=param_databases_dbsnp.get(
                "download_dbsnp_release_default", None
            ),
            dbsnp_url=param_databases_dbsnp.get("download_dbsnp_url", None),
            dbsnp_url_files=param_databases_dbsnp.get("download_dbsnp_url_files", None),
            dbsnp_url_files_prefix=param_databases_dbsnp.get(
                "download_dbsnp_url_files_prefix", None
            ),
            dbsnp_assemblies_map=param_databases_dbsnp.get(
                "download_dbsnp_assemblies_map", None
            ),
            dbsnp_vcf=param_databases_dbsnp.get("download_dbsnp_vcf", None),
            dbsnp_parquet=param_databases_dbsnp.get("download_dbsnp_parquet", None),
            genomes_folder=param.get("databases", {}).get("genomes_folder", None),
            threads=threads,
            memory=memory,
        )

    # HGMD
    if param_databases.get("hgmd", {}).get("convert_hgmd", None):
        log.debug(f"Convert HGMD")
        param_databases_hgmd = param_databases.get("hgmd", {})
        databases_download_hgmd(
            assemblies=assemblies,
            hgmd_folder=param_databases_hgmd.get("convert_hgmd", None),
            hgmd_file=param_databases_hgmd.get("convert_hgmd_file", None),
            output_basename=param_databases_hgmd.get("convert_hgmd_basename", None),
            genomes_folder=param.get("databases", {}).get("genomes_folder", None),
            threads=threads,
            memory=memory,
        )

    # from Annovar
    if args.input_annovar and args.output_annovar and args.genome:
        log.debug(f"Convert Annovar")
        from_annovar(args=args)

    # from ExtAnn
    if args.input_extann and args.output_extann:
        log.debug("Convert ExtAnn")
        from_extann(args=args)

    log.info("End")
