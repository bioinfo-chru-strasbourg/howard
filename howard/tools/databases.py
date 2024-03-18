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

from jproperties import Properties # jproperties 2.1.1


from howard.functions.commons import *
from howard.objects.variants import *
from howard.functions.databases import *




def databases(args:argparse) -> None:
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
    param = load_args(param={}, args=args, arguments_dict=arguments_dict, command="databases", strict=False, section_prefix=["databases"])
    #log.debug(f"param={param}")

    # Assembly
    assemblies = [value for value in args.assembly.split(',')]
    param = add_value_into_dict(dict_tree=param, sections=["databases", "assemblies"], value=assemblies)

    # Genome folder
    genomes_folder = param.get("databases", {}).get("genomes_folder", None)
    if param.get("databases", {}).get("genomes_folder", None) and not isinstance(param.get("databases", {}).get("genomes_folder", None), str):
        genomes_folder = param.get("databases", {}).get("genomes_folder", None).name
    if "config" in args and args.config and args.config.get("folders",{}).get("databases",{}).get("genomes") and param.get("databases", {}).get("genomes_folder", None) == DEFAULT_REFSEQ_FOLDER:
        genomes_folder = args.config.get("folders",{}).get("databases",{}).get("genomes")
    param = add_value_into_dict(dict_tree=param, sections=["databases", "genomes_folder"], value=genomes_folder)


    # Threads
    nb_threads = os.cpu_count()
    if "threads" in args:
        input_thread = args.threads
    else:
        input_thread = None
    if not input_thread or int(input_thread) <= 0:
        threads = nb_threads
    else:
        threads = int(input_thread)

    # Param
    if "generate_param" in args and args.generate_param:
        generate_databases_param(args=args, assemblies=assemblies)
        return None


    # Param database
    param_databases = param.get("databases", {})
    

    log.debug(f"param={param}")
    # Genomes
    if param_databases.get("genomes",{}).get("download_genomes"):
        log.debug(f"Download Genomes")
        if assemblies:
            databases_download_genomes(
                assemblies=param_databases.get("assemblies",[]), 
                genomes_folder=param_databases.get("genomes",{}).get("download_genomes"),
                provider=param_databases.get("genomes",{}).get("download_genomes_provider"),
                contig_regex=param_databases.get("genomes",{}).get("download_genomes_contig_regex"),
                threads=threads
                )

    # Annovar
    if args.download_annovar:
        log.debug(f"Download Annovar databases")
        if args.download_annovar_files:
            files = [value for value in args.download_annovar_files.split(',')]
        else:
            files = []
        databases_download_annovar(
            folder=args.download_annovar,
            files=files,
            assemblies = assemblies,
            annovar_url=args.download_annovar_url,
            threads=threads
            )

    # snpEff
    if args.download_snpeff:
        log.debug(f"Download snpEff databases")
        databases_download_snpeff(
            folder=args.download_snpeff,
            assemblies = assemblies,
            config=args.config,
            threads=threads
            )
        
    # refSeq
    if args.download_refseq:
        log.debug(f"Download refSeq databases")
        if args.download_refseq_files:
            files = [value for value in args.download_refseq_files.split(',')]
        else:
            files = []
        databases_download_refseq(
            assemblies = assemblies,
            refseq_folder=args.download_refseq,
            refseq_url=args.download_refseq_url,
            refseq_prefix=args.download_refseq_prefix,
            refseq_files=files,
            refseq_format_file=args.download_refseq_format_file,
            include_utr_5=args.download_refseq_include_utr5,
            include_utr_3=args.download_refseq_include_utr3,
            include_chrM=args.download_refseq_include_chrM,
            include_non_canonical_chr=args.download_refseq_include_non_canonical_chr,
            include_non_coding_transcripts=args.download_refseq_include_non_coding_transcripts,
            include_transcript_ver=args.download_refseq_include_transcript_version,
            threads=threads
            )

    # dbNSFP
    if args.download_dbnsfp:
        log.debug(f"Download dbNSFP")
        databases_download_dbnsfp(
            assemblies = assemblies,
            dbnsfp_folder=args.download_dbnsfp,
            dbnsfp_url=args.download_dbnsfp_url,
            dbnsfp_release=args.download_dbnsfp_release,
            threads=threads,
            parquet_size=args.download_dbnsfp_parquet_size,
            generate_sub_databases=args.download_dbnsfp_subdatabases,
            generate_parquet_file=args.download_dbnsfp_parquet,
            generate_vcf_file=args.download_dbnsfp_vcf,
            not_generate_files_all=args.download_dbnsfp_no_files_all,
            add_info=args.download_dbnsfp_add_info,
            row_group_size=args.download_dbnsfp_row_group_size,
            genomes_folder=param.get("databases", {}).get("genomes_folder", None)
            )

    # AlphaMissense
    if args.download_alphamissense:
        log.debug(f"Download AlphaMissense")
        databases_download_alphamissense(
            assemblies = assemblies,
            alphamissense_folder=args.download_alphamissense,
            alphamissense_url=args.download_alphamissense_url,
            threads=threads
            )
    
    # Exomiser
    if args.download_exomiser:
        log.debug(f"Download Exomiser")
        if "download_exomiser_application_properties" in args and args.download_exomiser_application_properties:
            if isinstance(args.download_exomiser_application_properties, str):
                download_exomiser_application_properties = args.download_exomiser_application_properties
            else:
                download_exomiser_application_properties = args.download_exomiser_application_properties.name
        else:
            download_exomiser_application_properties=None
        databases_download_exomiser(
            assemblies = assemblies,
            exomiser_folder=args.download_exomiser,
            exomiser_application_properties=download_exomiser_application_properties,
            exomiser_url=args.download_exomiser_url,
            exomiser_release=args.download_exomiser_release,
            exomiser_phenotype_release=args.download_exomiser_phenotype_release,
            exomiser_remm_release=args.download_exomiser_remm_release,
            exomiser_remm_url=args.download_exomiser_remm_url,
            exomiser_cadd_release=args.download_exomiser_cadd_release,
            exomiser_cadd_url=args.download_exomiser_cadd_url,
            exomiser_cadd_url_snv_file=args.download_exomiser_cadd_url_snv_file,
            exomiser_cadd_url_indel_file=args.download_exomiser_cadd_url_indel_file,
            threads=threads
            )

    # dbSNP
    if args.download_dbsnp:
        log.debug(f"Download dbSNP")
        if args.download_dbsnp_releases:
            dbsnp_releases = [value for value in args.download_dbsnp_releases.split(',')]
        else:
            dbsnp_releases = []
        databases_download_dbsnp(
            assemblies = assemblies,
            dbsnp_folder=args.download_dbsnp,
            dbsnp_releases=dbsnp_releases,
            dbsnp_release_default=args.download_dbsnp_release_default,
            dbsnp_url=args.download_dbsnp_url,
            dbsnp_url_files=args.download_dbsnp_url_files,
            dbsnp_url_files_prefix=args.download_dbsnp_url_files_prefix,
            dbsnp_assemblies_map=args.download_dbsnp_assemblies_map,
            genomes_folder=param.get("databases", {}).get("genomes_folder", None),
            threads=threads,
            dbsnp_vcf=args.download_dbsnp_vcf,
            dbsnp_parquet=args.download_dbsnp_parquet
            )

    # HGMD
    if args.convert_hgmd:
        log.debug(f"Convert HGMD")
        if "convert_hgmd_file" in args and args.convert_hgmd_file:
            if isinstance(args.convert_hgmd_file, str):
                convert_hgmd_file = args.convert_hgmd_file
            else:
                convert_hgmd_file = args.convert_hgmd_file.name
        databases_download_hgmd(
            assemblies = assemblies,
            hgmd_folder=args.convert_hgmd,
            hgmd_file=convert_hgmd_file,
            output_basename=args.convert_hgmd_basename,
            genomes_folder=param.get("databases", {}).get("genomes_folder", None),
            threads=threads
            )

    log.info("End")