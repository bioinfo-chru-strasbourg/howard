#!/usr/bin/env python

import argparse
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


from howard.commons import *


def query_and_concatenate_columns(parquet_file: str, output_file: str, columns: list):
    """
    This function performs an SQL query on a large Parquet file and concatenates multiple columns (if not empty),
    including the column name in the concatenation.

    :param parquet_file: The path to the Parquet file
    :type parquet_file: str
    :param output_file: The path to the output file where the concatenated data will be written
    :type output_file: str
    :param columns: The list of columns to concatenate
    :type columns: list
    """
    parquet_reader = pq.ParquetFile(parquet_file)
    num_rows = parquet_reader.metadata.num_rows

    chunk_size = 100000  # Adjust the chunk size as per your memory constraints

    with pd.HDFStore(output_file, mode='w') as store:
        for i in range(0, num_rows, chunk_size):
            table = parquet_reader.read_row_group(0, columns=columns, skip_rows=i, num_rows=chunk_size).to_pandas()
            table = table.apply(lambda row: ';'.join([f"{col}={row[col]}" for col in columns if pd.notnull(row[col])]), axis=1)
            store.append('data', table)

    # Export the concatenated data to a CSV file
    concatenated_data = pd.read_hdf(output_file, 'data')
    concatenated_data.to_csv(output_file + '.csv', index=False)


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

    databases_download(args)

    log.info("End")


def databases_download(args:argparse) -> None:
    """
    The `databases_download` function downloads genome, Annovar, and snpEff databases based on
    user-specified arguments.
    
    :param args: The `args` parameter is an object of the `argparse` module that contains the input
    arguments for the `databases_download` function. These arguments are used to determine which genome,
    Annovar, and snpEff databases to download
    :type args: argparse
    """

    log.debug(f"Args {args}")

    # Assembly
    assemblies = [value for value in args.assembly.split(',')]

    # Genomes
    if args.download_genomes:
        log.debug(f"Download Genomes")
        if assemblies:
            databases_download_genomes(
                assemblies=assemblies, 
                genome_folder=args.download_genomes,
                provider=args.download_genomes_provider,
                contig_regex=args.download_genomes_contig_regex,
                threads=args.threads
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
            annovar_url=args.download_annovar_url
            )

    # snpEff
    if args.download_snpeff:
        log.debug(f"Download snpEff databases")
        databases_download_snpeff(
            folder=args.download_snpeff,
            assemblies = assemblies,
            config=args.config
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
            include_transcript_ver=args.download_refseq_include_transcript_version
            )

    # dbNSFP
    if args.download_dbnsfp:
        log.debug(f"Download dbNSFP")
        databases_download_dbnsfp(
            assemblies = assemblies,
            dbnsfp_folder=args.download_dbnsfp,
            dbnsfp_url=args.download_dbnsfp_url,
            dbnsfp_release=args.download_dbnsfp_release,
            threads=args.threads,
            #nb_data_files=args.download_dbnsfp_nb_data_files,
            parquet_size=args.download_dbnsfp_parquet_size,
            generate_sub_databases=args.download_dbnsfp_subdatabases,
            generate_parquet_file=args.download_dbnsfp_parquet,
            generate_vcf_file=args.download_dbnsfp_vcf,
            #generate_vcf_file_all=args.download_dbnsfp_vcf_all,
            genomes_folder=args.genomes_folder,
            )



def databases_download_annovar(folder:str = None, files:list = None, assemblies:list = ["hg19"], annovar_url:str = "http://www.openbioinformatics.org/annovar/download") -> None:
    """
    This function downloads and extracts Annovar databases for specified assemblies and files.
    
    :param folder: The folder where the Annovar databases will be downloaded to
    :type folder: str
    :param files: The `files` parameter is a list of specific Annovar database files to download. If not
    provided, only the mandatory files will be downloaded. If set to "ALL", all available files will be
    downloaded
    :type files: list
    :param assemblies: A list of genome assemblies for which Annovar databases will be downloaded.
    Default is ["hg19"]
    :type assemblies: list
    :param annovar_url: The URL where Annovar databases can be downloaded from, defaults to
    http://www.openbioinformatics.org/annovar/download
    :type annovar_url: str (optional)
    """

    log.info(f"Download Annovar databases {assemblies}")

    # Minimum files to download
    files_minimum = ["refGene*"]

    for assembly in assemblies:

        log.debug(f"Download Annovar for assembly '{assembly}'")

        folder_assembly = f"{folder}/{assembly}"

        if not os.path.exists(folder_assembly):
            os.makedirs(folder_assembly)
        
        log.debug(f"Download Annovar databases in folder '{folder_assembly}'")

        if not files:
            log.debug("Only mandatory files will be downloaded")
            patterns_file_to_check = files_minimum
        elif files == ["ALL"]:
            log.debug("All files will be downloaded")
            patterns_file_to_check = None
        else:
            patterns_file_to_check = list(set(files+files_minimum))
            log.debug(f"Following files will be downloaded: {patterns_file_to_check}")

        log.debug(f"Download list of Annovar files from Annovar URL")

        avdblist_file = f"{assembly}_avdblist.txt"
        avdblist_url_file = f"{annovar_url}/{avdblist_file}"
        avdblist_folder_file = f"{folder_assembly}/{avdblist_file}"
        log.debug(f"Download list of Annovar files {avdblist_file} from Annovar URL {avdblist_url_file} to Annovar folder {avdblist_folder_file}")
        download_file(avdblist_url_file, avdblist_folder_file)

        if not os.path.exists(avdblist_folder_file):

            log.error(f"Download list of Annovar files from Annovar URL: {avdblist_url_file}")
            raise ValueError(f"Download list of Annovar files from Annovar URL: {avdblist_url_file}")
        
        else:

            # Open file
            with open(avdblist_folder_file, "r") as f:
                lines = f.readlines()
            # Dataframe
            data = []
            for line in lines:
                line = line.strip().split("\t")
                d = {
                    "file": line[0],
                    "version": line[1],
                    "size": line[2]
                }
                data.append(d)
            df = pd.DataFrame(data)

            # Create list of files associated to Annovar code
            files_to_download = []
            if not patterns_file_to_check:
                files_to_download = list(df["file"])
            else:
                for file in patterns_file_to_check:
                    files_to_download.append(f"{assembly}_{file}.txt.gz")
                    files_to_download.append(f"{assembly}_{file}.txt.idx.gz")
                    files_to_download.append(f"{assembly}_{file}.fa.gz")
                    files_to_download.append(f"{assembly}_{file}.zip")

            # Found if files exists in folder
            list_files_in_folder = []
            list_files_in_folder_size = {}
            # Loop through each file in the list
            for file in files_to_download:
                # Get a list of files in the folder that match the pattern
                files_in_folder = glob.glob(os.path.join(folder_assembly, file))
                if files_in_folder:
                    # If the file exists, add it
                    list_files_in_folder += [os.path.basename(f) for f in files_in_folder]
                    for f in files_in_folder:
                        list_files_in_folder_size[os.path.basename(f)] = os.path.getsize(f)
            # List of files in folder
            log.debug(f"List of file in Annovar folder: {list_files_in_folder}")
            log.debug(f"List of file in Annovar folder with size: {list_files_in_folder_size}")
            
            # Found if files exists in url
            matched_files = []
            #for pattern in list_files_in_folder:
            for pattern in files_to_download:
                matched_files += fnmatch.filter(df['file'], pattern)
            # deduplicate
            matched_files = list(set(matched_files))
            # List of files available
            log.debug(f"List of file available in Annovar URL: {matched_files}")

            # Find files in url but not in folder
            files_to_download = []
            for file in matched_files:
                file_url_size = next((d for d in data if fnmatch.fnmatch(d["file"], file)), None).get("size",0)
                file_folder_size = list_files_in_folder_size.get(file,0)
                if file not in list_files_in_folder or int(file_url_size) != int(file_folder_size):
                    files_to_download.append(file)

            # List of files to download
            log.debug(f"List of file to download: {files_to_download}")

            # Dowload and extract files
            if files_to_download:
                log.info(f"Download Annovar databases {[assembly]}")
                for file in files_to_download:
                    log.info(f"Download Annovar databases {[assembly]} - file '{file}' Downloading...")
                    file_url = os.path.join(annovar_url, file)
                    file_path = os.path.join(folder_assembly, file)

                    log.debug(f"Download file {file} from {file_url} to {file_path}...")
                    download_file(file_url, file_path)
                    log.debug(f"Extract file {file} to {folder_assembly}...")
                    extract_file(file_path)
            else:
                log.info(f"Download Annovar databases {[assembly]} - already exists")


def databases_download_snpeff(folder:str = None, assemblies:list = ["hg19"], config:dict = {}) -> None:
    """
    This function downloads and extracts snpEff databases for specified genome assemblies.
    
    :param folder: The folder where the snpEff databases will be downloaded and stored
    :type folder: str
    :param assemblies: The assemblies parameter is a list of genome assemblies for which the snpEff
    databases need to be downloaded
    :type assemblies: list
    :param config: The `config` parameter is a dictionary that contains information about the tools and
    their configurations. It is used to retrieve the path to the Java binary and the path to the snpEff
    binary
    :type config: dict
    """

    log.info(f"Download snpEff databases {assemblies}")

    # Java bin and snpEff jar
    snpeff_bin = get_snpeff_bin(config=config)
    java_bin = config.get("tools", {}).get("java", {}).get("bin", "java")

    # database list
    snpeff_databases_list = "snpeff_databases.list"
    snpeff_databases_list_path = os.path.join(folder,snpeff_databases_list)
    
    # create folder if not exists
    if folder:
        if not os.path.exists(folder):
            os.makedirs(folder)

    # For each assembly
    for assembly in assemblies:

        # Destination folder
        folder_assembly = os.path.join(folder, assembly)

        # Check if destination folder exists (already downloaded)
        if not os.path.exists(folder_assembly):

            # Download list of databases if file does not exists
            if not os.path.exists(snpeff_databases_list_path):
                snpeff_command_list_databases = f""" {java_bin} -jar {snpeff_bin} databases > {snpeff_databases_list_path} """
                log.info(f"snpEff databases downloading - list of databases...")
                log.debug(f"snpEff databases downloading - list of databases: {snpeff_command_list_databases}")
                command(snpeff_command_list_databases)
            
            # Open list of databases and create list of url
            snpeff_list_databases = {}
            with open(snpeff_databases_list_path, "r") as f:
                for line in f:
                    cols = line.strip().split()
                    snpeff_list_databases[cols[0]] = []
                    # Last 2 infos on file. Fix from snpEff between release 5.0 and 5.1
                    snpeff_list_databases[cols[0]].append(cols[-2].replace(",","").replace("]","").replace("[",""))
                    snpeff_list_databases[cols[0]].append(cols[-1].replace(",","").replace("]","").replace("[",""))

            # Strat download
            log.info(f"Download snpEff databases for assembly '{assembly}'...")
            #print(snpeff_list_databases.keys())
            # Try to download files
            file_path = None
            for file_url in snpeff_list_databases[assembly]:
                # File to be downloaded
                file_path = os.path.join(folder, os.path.basename(file_url))
                # Check if file already downloaded
                if os.path.exists(file_path):
                    break
                else:
                    # try to download
                    try:
                        log.debug(f"Download snpEff '{file_url}'...")
                        download_file(file_url, file_path)
                    # If fail, just pass to next url
                    except:
                        log.debug(f"Download snpEff '{file_url}' failed")

            # If download file exists
            if file_path is not None and os.path.exists(file_path):
                log.debug(f"Extract file {file_path} to {folder}...")
                # Extract file
                extract_file(file_path)
                # Move to destination folder
                extracted_folder = os.path.join(folder, "data", assembly)
                shutil.move(extracted_folder, folder_assembly)

        else:

            log.info(f"Download snpEff databases {[assembly]} - already exists")


def databases_download_genomes(assemblies: list, genome_folder: str = None, provider:str = "UCSC", contig_regex:str = None, threads:int = 1) -> None:
    """
    This function downloads genome assemblies using genomepy package with options to specify genome
    folder, provider, contig regex, and number of threads.
    
    :param assemblies: a list of genome assembly names to download
    :type assemblies: list
    :param genome_folder: The folder where the downloaded genome files will be saved. If no folder is
    specified, the default folder will be used
    :type genome_folder: str
    :param provider: The provider parameter specifies the source of the genome data. In this case, the
    default provider is set to "UCSC", which refers to the University of California, Santa Cruz Genome
    Browser. Other possible providers could include NCBI or Ensembl, defaults to UCSC
    :type provider: str (optional)
    :param contig_regex: The contig_regex parameter is a regular expression used to filter the contigs
    (chromosomes or scaffolds) to be downloaded for a given genome assembly. It allows users to download
    only a subset of the available contigs, based on their names or other characteristics. If
    contig_regex is not specified
    :type contig_regex: str
    :param threads: The "threads" parameter specifies the number of threads (parallel processes) to use
    for downloading the genomes. This can speed up the process if the computer has multiple cores or
    processors. The default value is 1, meaning that the download will be done using a single thread,
    defaults to 1
    :type threads: int (optional)
    :return: None is being returned.
    """

    log.info(f"Download Genomes {assemblies}")

    import genomepy

    if not genome_folder:
        genome_folder = DEFAULT_GENOME_FOLDER

    if os.path.exists(genome_folder):
        installed_genomes = genomepy.list_installed_genomes(genomes_dir=genome_folder)
    else:
        installed_genomes = []

    for assembly in assemblies:
        if assembly in installed_genomes:
            log.info(f"Download Genomes '{[assembly]}' - already exists")
        else:
            log.info(f"Download Genomes '{[assembly]}' downloading...")
            genomepy.install_genome(assembly, annotation=False, provider=provider, genomes_dir=genome_folder, threads=threads, regex=contig_regex)

    return None


def databases_download_refseq(assemblies:list, refseq_folder:str = None, refseq_url:str = None, refseq_prefix:str = "ncbiRefSeq", refseq_files:List = ["ncbiRefSeq.txt", "ncbiRefSeqLink.txt"], refseq_format_file:str = "ncbiRefSeq.txt", refseq_format_file_output:str = None, include_utr_5:bool = True, include_utr_3:bool = True, include_chrM:bool = True, include_non_canonical_chr:bool = True, include_non_coding_transcripts:bool = True, include_transcript_ver:bool = True) -> dict:
    """
    The `databases_download_refseq` function downloads RefSeq files for a list of assemblies and returns
    a dictionary of installed RefSeq files for each assembly.
    
    :param assemblies: A list of assemblies for which the RefSeq files need to be downloaded. Each
    assembly is represented as a string
    :type assemblies: list
    :param refseq_folder: The `refseq_folder` parameter is a string that specifies the folder where the
    RefSeq files will be downloaded and stored. If this parameter is not provided, a default folder will
    be used
    :type refseq_folder: str
    :param refseq_url: The `refseq_url` parameter is a string that represents the URL where the RefSeq
    files can be downloaded from
    :type refseq_url: str
    :param refseq_prefix: The `refseq_prefix` parameter is a string that specifies the prefix for the
    downloaded RefSeq files. By default, it is set to "ncbiRefSeq". This prefix is used to identify the
    RefSeq files for each assembly. For example, if the prefix is set to "ncbi, defaults to ncbiRefSeq
    :type refseq_prefix: str (optional)
    :param refseq_files: The `refseq_files` parameter is a list of filenames that need to be downloaded
    for each assembly. The default value is `["ncbiRefSeq.txt", "ncbiRefSeqLink.txt"]`, but you can
    provide your own list of filenames if needed
    :type refseq_files: List
    :param refseq_format_file: The `refseq_format_file` parameter is a string that specifies the
    filename of the RefSeq file that needs to be formatted. This file will be used as input for the
    `databases_format_refseq` function. By default, the value is set to "ncbiRefSeq.txt", defaults to
    ncbiRefSeq.txt
    :type refseq_format_file: str (optional)
    :param refseq_format_file_output: The `refseq_format_file_output` parameter is a string that
    specifies the output file path for the formatted RefSeq file. This file will be generated by the
    `databases_format_refseq` function and will contain the formatted RefSeq data. If this parameter is
    not provided, the formatted RefSeq
    :type refseq_format_file_output: str
    :param include_utr_5: A boolean parameter that specifies whether to include the 5' untranslated
    region (UTR) in the downloaded RefSeq files. If set to True, the 5' UTR will be included. If set to
    False, the 5' UTR will be excluded, defaults to True
    :type include_utr_5: bool (optional)
    :param include_utr_3: The `include_utr_3` parameter is a boolean that specifies whether to include
    the 3' untranslated region (UTR) in the downloaded RefSeq files. If set to `True`, the 3' UTR will
    be included. If set to `False`, the 3', defaults to True
    :type include_utr_3: bool (optional)
    :param include_chrM: A boolean parameter that determines whether to include the mitochondrial
    chromosome (chrM) in the downloaded RefSeq files. If set to True, the chrM will be included; if set
    to False, it will be excluded, defaults to True
    :type include_chrM: bool (optional)
    :param include_non_canonical_chr: The `include_non_canonical_chr` parameter is a boolean value that
    determines whether or not to include non-canonical chromosomes in the downloaded RefSeq files. If
    set to `True`, non-canonical chromosomes will be included. If set to `False`, non-canonical
    chromosomes will be excluded, defaults to True
    :type include_non_canonical_chr: bool (optional)
    :param include_non_coding_transcripts: The parameter `include_non_coding_transcripts` is a boolean
    flag that determines whether non-coding transcripts should be included in the downloaded RefSeq
    files. If set to `True`, non-coding transcripts will be included. If set to `False`, non-coding
    transcripts will be excluded, defaults to True
    :type include_non_coding_transcripts: bool (optional)
    :param include_transcript_ver: The `include_transcript_ver` parameter is a boolean value that
    determines whether to include the transcript version in the downloaded RefSeq files. If set to
    `True`, the transcript version will be included. If set to `False`, the transcript version will be
    excluded, defaults to True
    :type include_transcript_ver: bool (optional)
    :return: The function `databases_download_refseq` returns a dictionary `installed_refseq` which
    contains information about the downloaded RefSeq files for each assembly. The keys of the dictionary
    are the assembly names, and the values are lists of the installed RefSeq files for each assembly.
    """

    # Log
    log.info(f"Download refSeq databases {assemblies}")

    # Default refSeq Folder
    if not refseq_folder:
        refseq_folder = DEFAULT_REFSEQ_FOLDER

    # Default refSeq URL
    if not refseq_url:
        refseq_url = DEFAULT_REFSEQ_URL

    # Create folder if not exists
    if not os.path.exists(refseq_folder):
        os.makedirs(refseq_folder)

    # Installed refSeq files
    installed_refseq = {}

    for assembly in assemblies:

        # Create folder if not exists
        if not os.path.exists(f"{refseq_folder}/{assembly}"):
            os.makedirs(f"{refseq_folder}/{assembly}")

        # Strat download needed files
        installed_refseq[assembly] = []
        if os.path.exists(refseq_folder):
            log.debug(f"Download refSeq databases {assemblies} - '{assembly}'")
            existing_files_path = f"{refseq_folder}/{assembly}"
            existing_files = glob.glob(rf'{existing_files_path}/{refseq_prefix}*', recursive=True)
            # For refSeq files to download
            for refseq_file in refseq_files:
                new_refseq_file = False
                # If refSeq file exists
                if f"{existing_files_path}/{refseq_file}" in existing_files:
                    log.info(f"Download refSeq databases ['{assembly}'] - '{refseq_file}' already exists")
                    installed_refseq[assembly].append(refseq_file)
                # If refSeq file need to be downloaded
                else:
                    # files to download
                    file_url = f"{refseq_url}/{assembly}/database/{refseq_file}.gz"
                    file_path = f"{refseq_folder}/{assembly}/{refseq_file}.gz"
                    # try to download
                    try:
                        log.info(f"Download refSeq databases ['{assembly}'] - '{refseq_file}' downloading...")
                        # Download file
                        download_file(file_url, file_path)
                        # Extract file
                        extract_file(file_path)
                        # add to installed files
                        installed_refseq[assembly].append(refseq_file)
                        new_refseq_file = True
                    # If fail, just pass to next url
                    except:
                        log.debug(f"Download refSeq databases {assemblies} - '{assembly}' - '{refseq_file}' downloading failed")
                        raise ValueError(f"Download refSeq databases {assemblies} - '{assembly}' - '{refseq_file}' downloading failed")
                # format refSeq
                file_path = f"{refseq_folder}/{assembly}/{refseq_file}"
                if refseq_format_file_output:
                    file_path_bed = refseq_format_file_output
                else:
                    file_path_bed = re.sub(r'txt$', 'bed', file_path)
                if refseq_file == refseq_format_file and (new_refseq_file or not os.path.exists(file_path_bed)):
                    log.info(f"Download refSeq databases ['{assembly}'] - '{refseq_file}' formatting...")
                    databases_format_refseq(refseq_file=file_path, output_file=file_path_bed, include_utr_5=include_utr_5, include_utr_3=include_utr_3, include_chrM=include_chrM, include_non_canonical_chr=include_non_canonical_chr, include_non_coding_transcripts=include_non_coding_transcripts, include_transcript_ver=include_transcript_ver)

    # Log
    log.debug(f"installed_refseq: {installed_refseq}")

    return installed_refseq


def databases_format_refseq(refseq_file:str, output_file:str, include_utr_5:bool = True, include_utr_3:bool = True, include_chrM:bool = True, include_non_canonical_chr:bool = True, include_non_coding_transcripts:bool = True, include_transcript_ver:bool = True) -> str:
    """
    The function `databases_format_refseq` takes a RefSeq file as input and formats it according to
    specified criteria, such as including or excluding certain features, and outputs the formatted file.
    
    :param refseq_file: The `refseq_file` parameter is a string that represents the path to the input
    RefSeq file. This file contains information about gene annotations, including chromosome, start and
    end positions, strand, and other details
    :type refseq_file: str
    :param output_file: The output file is the name of the file where the formatted data will be written
    :type output_file: str
    :param include_utr_5: The parameter `include_utr_5` determines whether to include the 5' UTR
    (untranslated region) in the output. If set to `True`, the 5' UTR will be included. If set to
    `False`, the 5' UTR will be excluded, defaults to True
    :type include_utr_5: bool (optional)
    :param include_utr_3: The parameter `include_utr_3` determines whether to include the 3' UTR
    (untranslated region) in the output. If set to `True`, the 3' UTR will be included. If set to
    `False`, the 3' UTR will be excluded, defaults to True
    :type include_utr_3: bool (optional)
    :param include_chrM: A boolean parameter that determines whether to include transcripts from the
    mitochondrial chromosome (chrM or chrMT) in the output file. If set to True, transcripts from the
    mitochondrial chromosome will be included. If set to False, transcripts from the mitochondrial
    chromosome will be excluded, defaults to True
    :type include_chrM: bool (optional)
    :param include_non_canonical_chr: The parameter `include_non_canonical_chr` determines whether or
    not to include non-canonical chromosomes in the output. If set to `True`, non-canonical chromosomes
    will be included. If set to `False`, non-canonical chromosomes will be excluded, defaults to True
    :type include_non_canonical_chr: bool (optional)
    :param include_non_coding_transcripts: The parameter `include_non_coding_transcripts` determines
    whether non-coding transcripts should be included in the output file. If set to `True`, non-coding
    transcripts will be included. If set to `False`, non-coding transcripts will be excluded, defaults
    to True
    :type include_non_coding_transcripts: bool (optional)
    :param include_transcript_ver: The parameter `include_transcript_ver` determines whether to include
    the transcript version in the output file. If set to `True`, the transcript version will be included
    in the output file. If set to `False`, the transcript version will be removed from the output file,
    defaults to True
    :type include_transcript_ver: bool (optional)
    :return: the path of the output file.
    """

    # Open refSeq file
    with open(refseq_file, "r") as fi:

        # Open refSeq output file
        with open(output_file, "w") as fo:

            # For each line
            for l in fi:

                # Header
                if not l.startswith("#"):

                    # line
                    l = l.strip().split()
                    
                    # check non-coding transcripts or unwanted chr
                    condition = (
                        (include_non_coding_transcripts or l[1].startswith("NM_"))
                        and (
                            (include_non_canonical_chr or re.match(r'^chr([0-9]*|X|Y)$', l[2]))
                            or (include_chrM and l[2] in {"chrM", "chrMT"})
                        )
                    )

                    # include transcript
                    if condition:

                        # Transcript version
                        if not include_transcript_ver:
                            l[1] = l[1].split('.')[0]
                        
                        # Exon to write (depending on utr)
                        # packing exons together to be able to number them (need to start from the end if negative strand, if using IGV as standard)
                        exon_starts = l[9][:-1].split(",")
                        exon_ends = l[10][:-1].split(",")
                        exons_to_write = []
                        for i in range(len(exon_starts)):
                            
                            if include_utr_5 and include_utr_3:
                                exons_to_write.append([l[2], exon_starts[i], exon_ends[i], l[12], l[1], l[3]])
                    
                            elif not include_utr_5 and include_utr_3:
                                if int(exon_starts[i]) < int(l[6]):
                                    exons_to_write.append([l[2], l[6], exon_ends[i], l[12], l[1], l[3]])
                                else:
                                    exons_to_write.append([l[2], exon_starts[i], exon_ends[i], l[12], l[1], l[3]])
                    
                            elif include_utr_5 and not include_utr_3:
                                if int(exon_ends[i]) > int(l[7]):
                                    exons_to_write.append([l[2], exon_starts[i], l[7], l[12], l[1], l[3]])
                                else:
                                    exons_to_write.append([l[2], exon_starts[i], exon_ends[i], l[12], l[1], l[3]])
                    
                            elif not include_utr_5 and not include_utr_3:
                                if int(exon_starts[i]) >= int(l[6]) and int(exon_ends[i]) <= int(l[7]):
                                    exons_to_write.append([l[2], exon_starts[i], exon_ends[i], l[12], l[1], l[3]])
                                elif int(exon_starts[i]) < int(l[6]) and int(l[6]) <= int(exon_ends[i]) <= int(l[7]):
                                    exons_to_write.append([l[2], l[6], exon_ends[i], l[12], l[1], l[3]])
                                elif int(l[6])<= int(exon_starts[i]) <= int(l[7]) and int(exon_ends[i]) > int(l[7]):
                                    exons_to_write.append([l[2], exon_starts[i], l[7], l[12], l[1], l[3]])
                                elif int(exon_starts[i]) <= int(l[6]) <= int(exon_ends[i]) and int(exon_starts[i]) <= int(l[7]) <= int(exon_ends[i]):
                                    exons_to_write.append([l[2], l[6], l[7], l[12], l[1], l[3]])

                        # Write exon (depending on strand)
                        if l[3] == "-":
                            exon_num = len(exons_to_write)
                            exon_num_incrementation = -1
                        else:
                            exon_num = 1
                            exon_num_incrementation = 1

                        for i in range(len(exons_to_write)):
                            exons_to_write[i].append("exon"+str(exon_num))
                            exon_num += exon_num_incrementation
                            fo.write("\t".join(exons_to_write[i])+"\n")
                
    return output_file



def databases_download_dbnsfp(assemblies:list, dbnsfp_folder:str = None, dbnsfp_url:str = None, dbnsfp_release:str = "4.4a", threads:int = None, nb_data_files:int = None, parquet_size:int = 100, generate_parquet_file:bool = False, generate_sub_databases:bool = False, generate_vcf_file:bool = False, generate_vcf_file_all:bool = True, genomes_folder:str = None) -> bool:
    """
    The function `databases_download_dbnsfp` is used to download and process dbNSFP databases for specified
    genome assemblies.
    
    :param assemblies: A list of genome assemblies for which to download and process dbNSFP data. Each
    assembly should be specified as a string
    :type assemblies: list
    :param dbnsfp_folder: The `dbnsfp_folder` parameter is a string that specifies the folder where the
    dbNSFP database files are located. If this parameter is not provided, the function will attempt to
    download the dbNSFP database files from the `dbnsfp_url` parameter
    :type dbnsfp_folder: str
    :param dbnsfp_url: The URL where the dbNSFP database files can be downloaded from
    :type dbnsfp_url: str
    :param dbnsfp_release: The version of the dbNSFP database to be used. The default value is "4.4a",
    defaults to 4.4a
    :type dbnsfp_release: str (optional)
    :param threads: The `threads` parameter specifies the number of threads to use for parallel
    processing. It determines how many tasks can be executed simultaneously. Increasing the number of
    threads can potentially speed up the execution time of the function, especially if there are
    multiple cores available on the machine
    :type threads: int (optional)
    :param nb_data_files: The parameter "nb_data_files" is used to specify the number of data files to
    be processed. It is an optional parameter and its value should be an integer
    :type nb_data_files: int (optional)
    :param parquet_size: The parameter "parquet_size" is used to specify the maximum size (Mb) of data files in parquet folder. It is an optional parameter and its value should be an integer
    :type parquet_size: int (optional)
    :param generate_parquet_file: A boolean flag indicating whether to generate a parquet file or not,
    defaults to False
    :type generate_parquet_file: bool (optional)
    :param generate_sub_databases: A boolean parameter that determines whether to generate sub-databases
    or not. If set to True, the function will generate sub-databases based on the assemblies provided.
    If set to False, the function will not generate sub-databases, defaults to False
    :type generate_sub_databases: bool (optional)
    :param generate_vcf_file: A boolean flag indicating whether to generate a VCF file or not,
    defaults to False
    :type generate_vcf_file: bool (optional)
    :param generate_vcf_file_all: A boolean flag indicating whether to generate a ALL database VCF file or not,
    defaults to True
    :type generate_vcf_file_all: bool (optional)
    :param genomes_folder: A string that specifies where are genomes
    :type genomes_folder: str (optional)
    :return: bool as success or not
    """

    def get_database_files(pattern:str) -> List:
        """
        The function `get_database_files` returns a list of database files that match a given pattern,
        sorted by file size in descending order.
        
        :param pattern: The pattern parameter is a string that represents the file pattern to search
        for. It can include wildcards such as "*", which matches any number of characters, and "?",
        which matches any single character. The pattern is used to search for files in a specific
        directory or directories
        :type pattern: str
        :return: a list of database files that match the given pattern.
        """
        
        database_files = sorted(glob.glob(pattern), key=os.path.getsize, reverse=True)
        return database_files
    

    def get_columns_structure(database_file:str, sample_size:int = 1000000, threads:int = 1) -> dict:
        """
        The function `get_columns_structure` creates a view of a database file and retrieves the structure
        of its columns.
        
        :param database_file: The `database_file` parameter is a string that represents the path to the
        database file that you want to retrieve the column structure from
        :type database_file: str
        :param sample_size: The `sample_size` parameter is an optional parameter that specifies the
        number of rows to sample from the database file in order to determine the structure of the
        columns. By default, it is set to 1,000,000. This means that the function will read a sample of
        1,000, defaults to 1000000
        :type sample_size: int (optional)
        :return: a dictionary called `columns_structure` which contains the column names and their
        corresponding data types from the database file.
        """
        
        # Check columns structure
        log.info(f"Download dbNSFP {assemblies} - Check database structure...")

        # Create view and retrive structure
        query = f"""
                    CREATE VIEW view_for_structure AS
                    (SELECT *
                    FROM read_csv_auto('{database_file}', compression=gzip, ALL_VARCHAR=0, delim='\t', nullstr='.', sample_size={sample_size})
                )
            """
        #log.debug(query)
        db_structure = duckdb.connect(config={"threads":threads})
        db_structure.execute(query)
        res_structure = db_structure.query("PRAGMA table_info('view_for_structure');")

        # Constructure structure dict
        columns_structure = {}
        for column in res_structure.df().iterrows():
            column_name = column[1].iloc[1]
            data_type = column[1].iloc[2]
            columns_structure[column_name] = data_type

        # Log
        log.info(f"Download dbNSFP {assemblies} - Check database structure - {len(columns_structure)} columns found")

        db_structure.close()

        if not columns_structure:
            log.error(f"Download dbNSFP {assemblies} - Check database structure - No database structure found")
            raise ValueError(f"Download dbNSFP {assemblies} - Check database structure - No database structure found")
        else:
            log.info(f"Download dbNSFP {assemblies} - Check database structure - Database structure generated")

        return columns_structure
    

    def clean_name(name:str) -> str:
        """
        The clean_name function replaces hyphens and plus signs in a string with underscores.
        
        :param name: The `name` parameter is a string that represents a name
        :type name: str
        :return: The function `clean_name` returns a string.
        """

        mapping_table = str.maketrans({'-': '_', '+': '_'})
        return name.translate(mapping_table)


    def get_columns_select_clause(columns_structure:dict, assembly:str = "hg19", for_parquet:bool = False, sub_database:str = None, null_if_no_annotation:bool = True, print_log:bool = True) -> dict:
        """
        The `get_columns_select_clause` function generates the SELECT and WHERE clauses for a SQL query
        based on the input columns structure, assembly, sub-database, and null_if_no_annotation flag.
        
        :param columns_structure: A dictionary that represents the structure of the columns. Each key is
        a column name and the corresponding value is the column type
        :type columns_structure: dict
        :param assembly: The `assembly` parameter specifies the genome assembly version. It can take
        values like "hg19", "hg38", or "hg18". By default, it is set to "hg19", defaults to hg19
        :type assembly: str (optional)
        :param for_parquet: The `for_parquet` parameter is a boolean flag that determines whether the
        SQL query is being generated for a Parquet file. If `for_parquet` is set to `True`, the function
        will include the position and REF/ALT columns in the SELECT clause, regardless of whether they
        are specified, defaults to False
        :type for_parquet: bool (optional)
        :param sub_database: The `sub_database` parameter is an optional parameter that specifies a
        sub-database within the main database. It is used to filter the columns based on the specified
        sub-database. If a sub-database is provided, only columns that belong to that sub-database or
        start with the sub-database
        :type sub_database: str
        :param null_if_no_annotation: The `null_if_no_annotation` parameter is a boolean parameter that
        determines whether to return null if there are no annotations found. If set to True, the
        function will return null if there are no annotations. If set to False, the function will return
        an empty string if there are no annotations, defaults to True
        :type null_if_no_annotation: bool (optional)
        :param log: The `log` parameter is a boolean flag that determines whether to enable logging or
        not. If set to `True`, the function will log debug messages during execution. If set to `False`,
        logging will be disabled, defaults to True
        :type log: bool (optional)
        :return: The `get_columns_select_clause` function returns a dictionary containing the "select"
        and "where" clauses for a SQL query. The "select" clause includes the columns to be selected in
        the query, while the "where" clause includes the conditions for filtering the data.
        Additionally, the dictionary also includes a "annotations" key that contains a dictionary of the
        selected annotations and their corresponding column types.
        """
        

        columns_select_position = {}
        columns_select_ref_alt = {}
        columns_select_annotations = {}
        columns_where_annotations = {}
        columns_list_annotations = {}
        variant_position_columns =  {
            "hg38": {
                    "#chr": {
                        "alias": "#CHROM",
                        "type": "VARCHAR",
                        "prefix": "chr"
                    },
                    "pos(1-based)": {
                        "alias": "POS",
                        "type": "BIGINT"
                    },
                },
            "hg19": {
                    "hg19_chr": {
                        "alias": "#CHROM",
                        "type": "VARCHAR",
                        "prefix": "chr"
                    },
                    "hg19_pos(1-based)": {
                        "alias": "POS",
                        "type": "BIGINT"
                    },
                },
            "hg18": {
                    "hg18_chr": {
                        "alias": "#CHROM",
                        "type": "VARCHAR",
                        "prefix": "chr"
                    },
                    "hg18_pos(1-based)": {
                        "alias": "POS",
                        "type": "BIGINT"
                    },
                }
        }

        for column in columns_structure:
            
            # REF/ALT columns
            if column in ["ref", "alt"]:
                column_alias = column.upper()
                column_key = f""" "{column}" AS "{column_alias}" """
                columns_select_ref_alt[column_key] = "VARCHAR"

            # Keep assembly position columns
            elif column in variant_position_columns.get(assembly, {}).keys():

                column_prefix = variant_position_columns.get(assembly, {}).get(column,{}).get("prefix","")
                column_alias = variant_position_columns.get(assembly, {}).get(column,{}).get("alias",None)
                column_type = variant_position_columns.get(assembly, {}).get(column,{}).get("type",None)
                if column_prefix:
                    column_prefix_concat = f""" '{column_prefix}' ||  """
                else:
                    column_prefix_concat = ""
                column_key = f""" {column_prefix_concat} "{column}" AS "{column_alias}" """
                columns_select_position[column_key] = column_type
                # Force type on structure
                columns_structure[column] = variant_position_columns.get(assembly, {}).get(column,{}).get("type",None)
            
            else:
                column_to_remove = False
                # Remove other assemby position columns
                for other_assembly in variant_position_columns.keys():
                    if column in variant_position_columns.get(other_assembly, {}).keys():
                        # Force type on structure
                        column_type = variant_position_columns.get(other_assembly, {}).get(column,{}).get("type",None)
                        columns_structure[column] = column_type
                        column_to_remove = True
                
                # other columns - annotations
                if not column_to_remove:
                    if sub_database is None or sub_database in ['ALL'] or column == sub_database or column.startswith(f"{sub_database}_"):
                        column_alias = clean_name(column)
                        if for_parquet:
                            column_key = f""" "{column_alias}" """
                            column_where = f""" "{column_alias}" IS NOT NULL """
                        else:
                            if column == column_alias:
                                column_key = f""" "{column}" """
                            else:
                                column_key = f""" "{column}" AS "{column_alias}" """
                            column_where = f""" "{column}" IS NOT NULL """
                        columns_select_annotations[column_key] = columns_structure[column]
                        columns_where_annotations[column_where] = columns_structure[column]
                        columns_list_annotations[column_alias] = columns_structure[column]

        # Force position and ref/alt if parquet
        if for_parquet:
            columns_select_position = {'"#CHROM"': '#CHROM', '"POS"': 'POS'}
            columns_select_ref_alt = {'"REF"': 'REF', '"ALT"': 'ALT'}

        # Create select clause
        columns_select = list(columns_select_position.keys()) + list(columns_select_ref_alt.keys()) + list(columns_select_annotations.keys())
        columns_select_clause = ", ".join(columns_select)

        # Create where clause
        columns_where = list(columns_where_annotations.keys())
        columns_where_clause = " OR ".join(columns_where)

        # log
        if print_log:
            log.debug(f"""Download dbNSFP ['{assembly}'] - Database '{sub_database}' - Check database structure - {len(list(columns_select_annotations.keys()))} annotations found""")

        # If no annotation found (usually for sub_database ref, alt, pos...)
        if null_if_no_annotation and not len(list(columns_select_annotations.keys())):
            columns_select_clause = None

        # Create output
        clauses = {
            "select": columns_select_clause,
            "where": columns_where_clause,
            "annotations": columns_list_annotations
        }

        return clauses


    def get_annotation_description(readme:str) -> dict:
        """
        The function `get_annotation_description` reads a README file and extracts annotation
        descriptions from it, storing them in a dictionary.
        
        :param readme: The `readme` parameter is a string that represents the file path to a README
        file. This function reads the contents of the README file and extracts information about
        annotation descriptions
        :type readme: str
        """

        with open(readme, "r") as f:
            lines = f.readlines()
        line_in_variant = False
        annotation_dict = {}
        for line in lines:
            if line.startswith("Columns of dbNSFP_variant"):
                line_in_variant = True
            elif line.startswith("\n"):
                line_in_variant = False
            elif line_in_variant:
                line_split = line.split("\t")
                if line_split[0] != "":
                    line_num = int(line_split[0])
                    line_name = line_split[1].split(":")[0]
                    line_desc = ":".join(line_split[1].split(":")[1:]).strip()
                else:
                    if line_split[2][0].isupper() or not line_split[2][0].isalnum():
                        sep = "."
                    else:
                        sep = ""
                    line_desc += f"""{sep} {line_split[2].strip()}"""
                annotation_dict[line_name] = line_desc.replace('"', "'")

        return annotation_dict


    # Log
    log.info(f"Download dbNSFP {assemblies}")

    # genomes folder
    if not genomes_folder:
        genomes_folder = DEFAULT_GENOME_FOLDER

    # Default refSeq Folder
    if not dbnsfp_folder:
        dbnsfp_folder = DEFAULT_DBNSFP_FOLDER

    # Default refSeq URL
    if not dbnsfp_url:
        dbnsfp_url = DEFAULT_DBNSFP_URL

    # Create folder if not exists
    if not os.path.exists(dbnsfp_folder):
        Path(dbnsfp_folder).mkdir(parents=True, exist_ok=True)

    # Files for dbNSFP
    # https://dbnsfp.s3.amazonaws.com/dbNSFP4.4a.zip
    dbnsfp_zip = f"dbNSFP{dbnsfp_release}.zip"
    dbnsfp_readme = f"dbNSFP{dbnsfp_release}.readme.txt"
    dbnsfp_zip_url = os.path.join(dbnsfp_url, dbnsfp_zip)
    dbnsfp_zip_dest = os.path.join(dbnsfp_folder, dbnsfp_zip)
    dbnsfp_readme_dest = os.path.join(dbnsfp_folder, dbnsfp_readme)
    dbnsfp_zip_dest_folder = os.path.dirname(dbnsfp_zip_dest)

    # Download dbNSFP
    if not os.path.exists(dbnsfp_zip_dest):
        log.info(f"Download dbNSFP {assemblies} - Download '{dbnsfp_zip}'...")
        download_file(dbnsfp_zip_url, dbnsfp_zip_dest)
    else:
        log.info(f"Download dbNSFP {assemblies} - Database '{dbnsfp_zip}' already exists")

    # Extract dbNSFP
    if not os.path.exists(dbnsfp_readme_dest):
        log.info(f"Download dbNSFP {assemblies} - Extract '{dbnsfp_zip}'...")
        extract_file(dbnsfp_zip_dest)
    else:
        log.info(f"Download dbNSFP {assemblies} - Database '{dbnsfp_zip}' already extracted")
    


    def get_header(header_file:str, dbnsfp_readme_dest:str = None, columns_structure:dict = {}, annotations:list = [], readme_annotations_description:dict = None) -> bool:
        """
        The function `get_header` generates a VCF header based on provided annotations and writes it to
        a file.
        
        :param header_file: The path to the output file where the VCF header will be written
        :type header_file: str
        :param dbnsfp_readme_dest: The `dbnsfp_readme_dest` parameter is a string that represents the
        destination file path for the dbNSFP readme file. This file contains the description of the
        annotations present in the dbNSFP database
        :type dbnsfp_readme_dest: str
        :param columns_structure: The `columns_structure` parameter is a dictionary that maps annotation
        names to their corresponding column types in a database. The keys of the dictionary are the
        annotation names, and the values are the column types
        :type columns_structure: dict
        :param annotations: The `annotations` parameter is a list of strings that represents the
        annotations to be included in the VCF header
        :type annotations: list
        :return: a boolean value indicating whether the header file was successfully created.
        """

        # Default VCF header
        default_header_list = [
                '##fileformat=VCFv4.2',
                '#CHROM\tPOS\tREF\tALT\t' + "\t".join(annotations)
                ]
        header_vcf = vcf.Reader(io.StringIO("\n".join(default_header_list)))

        code_type_map = {
                "Integer": 0,
                "String": 1,
                "Float": 2,
                "Flag": 3
            }
        code_type_map_from_sql = {
            "BIGINT": "Integer",
            "VARCHAR": "String",
            "DOUBLE": "Float",
            "BOOLEAN": "Integer"
        }

        # Extract annotation description
        if not readme_annotations_description:
            readme_annotations_description = get_annotation_description(dbnsfp_readme_dest)

        # Number of annotations
        #log.debug(f"""{len(readme_annotations_description)} annotations found """)

        # Add annotations to header
        for annotation in readme_annotations_description:
            annotation_name = clean_name(annotation)

            if annotation_name in annotations:
                #log.debug(f"{annotation}: {readme_annotations_description[annotation]}")
                vcf_header_infos_number = "."
                vcf_header_infos_type = code_type_map_from_sql.get(columns_structure.get(annotation,"VARCHAR"),"String")
                vcf_header_infos_description = readme_annotations_description[annotation] + f" [dbNSFP{dbnsfp_release}]"
                vcf_header_infos_source = dbnsfp_release
                vcf_header_infos_version = dbnsfp_release

                header_vcf.infos[annotation_name] = vcf.parser._Info(
                    annotation_name,
                    vcf_header_infos_number,
                    vcf_header_infos_type,
                    vcf_header_infos_description,
                    vcf_header_infos_source,
                    vcf_header_infos_version,
                    code_type_map[vcf_header_infos_type]
                )

        f = open(header_file, 'w')
        vcf.Writer(f, header_vcf)
        f.close()

        return os.path.exists(header_file)


    def size_to_mb(size:int) -> int:
        """
        The function `size_to_mb` converts a given size in megabytes to bytes.
        
        :param size: The size parameter represents the size of a partition in megabytes
        :type size: int
        :return: the partition size in megabytes as an integer.
        """

        if size:
            byte = 1
            kilobytes = byte * 1024
            megabytes = kilobytes * 1024
            return int(size * megabytes)
        else:
            return None


    # Init
    partition_size = size_to_mb(parquet_size)
    sample_size = 1000000
    database_files = []
    columns_structure = {}
    readme_annotations_description = None
    if not threads:
        threads = os.cpu_count()
    if not nb_data_files:
        nb_data_files = threads

    # Create duckdb connexion
    db = duckdb.connect(config={"threads":threads})

    parquet_partition_generated_list = []
    parquet_partition_already_generated_list = []
    parquet_generated_list = []
    parquet_already_generated_list = []
    vcf_generated_list = []
    vcf_already_generated_list = []

    for assembly in assemblies:

        output_assembly = f"{dbnsfp_folder}/{assembly}"
        if os.path.exists(output_assembly):
            log.info(f"Download dbNSFP ['{assembly}'] - Parquet folder already exists")
            continue
        else:
            # Create output chromosome
            Path(output_assembly).mkdir(parents=True, exist_ok=True)
  
        # Check database files
        log.info(f"Download dbNSFP {assemblies} - Check database files...")
        database_files = get_database_files(f"{dbnsfp_zip_dest_folder}/dbNSFP{dbnsfp_release}_variant.chr*.gz")
        log.debug(f"Download dbNSFP {assemblies} - Check database files - Found {len(database_files)} files")
        if not database_files:
            log.error(f"Download dbNSFP {assemblies} - Check database files - No files found")
            raise ValueError(f"Download dbNSFP {assemblies} - Check database files - No files found")

        # Database structure
        if not columns_structure:
            columns_structure = get_columns_structure(database_file=database_files[0], sample_size=sample_size, threads=threads)

        #log.info(f"Download dbNSFP ['{assembly}']")

        # Output prefix
        output_prefix = f"{dbnsfp_folder}/{assembly}/dbNSFP{dbnsfp_release}"

        # Sub database
        sub_database = "ALL"

        # Create clauses (select and where)
        columns_clauses = get_columns_select_clause(columns_structure=columns_structure, assembly=assembly, sub_database=sub_database)
        columns_select_clause = columns_clauses.get("select")
        columns_annotations = columns_clauses.get("annotations")
        
        # If clause select not empty
        if columns_select_clause:

            log.info(f"""Download dbNSFP ['{assembly}']""")
            log.info(f"""Download dbNSFP ['{assembly}'] - Parquet folder generation""")

            # Process database files
            for database_file in sorted(database_files, key=str.casefold, reverse=False):

                log.info(f"""Download dbNSFP ['{assembly}'] - Database '{sub_database}' - Process file '{os.path.basename(database_file)}'...""")

                # Chromosome
                log.debug(f"""Download dbNSFP ['{assembly}'] - Database '{sub_database}' - Process file '{os.path.basename(database_file)}' - Check chomosome...""")
                query = f"""
                            SELECT "#chr"
                            FROM read_csv_auto('{database_file}', compression=gzip, ALL_VARCHAR=1, delim='\t', nullstr='.', SAMPLE_SIZE=10)
                            WHERE \"#chr\" IS NOT NULL
                            LIMIT 1
                    """
                res = db.query(query).df()
                chromosome = "chr" + res["#chr"][0]
                log.debug(f"""Download dbNSFP ['{assembly}'] - Database '{sub_database}' - Process file '{os.path.basename(database_file)}' - Found chromosome '{chromosome}'""")

                # Output
                output = f"{output_prefix}.{sub_database}.partition.parquet"
                if not os.path.exists(output):
                    Path(output).mkdir(parents=True, exist_ok=True)
                output_parquet = os.path.join(output, f"#CHROM={chromosome}")

                # Query
                query = f"""
                        COPY (
                            SELECT {columns_select_clause}
                            FROM read_csv_auto(
                                '{database_file}',
                                header=1,
                                delim='\t',
                                compression=gzip,
                                nullstr='.',
                                columns={columns_structure},
                                types={columns_structure}
                                )
                            WHERE \"#CHROM\" IS NOT NULL
                            )
                        TO '{output_parquet}' WITH (FORMAT PARQUET, PER_THREAD_OUTPUT TRUE)
                    """
                
                # Log
                log.debug(f"""Download dbNSFP ['{assembly}'] - Database '{sub_database}' - Process file '{os.path.basename(database_file)}' - Write parquet files...""")

                if os.path.exists(f"{output_parquet}"):
                    log.debug(f"""Download dbNSFP ['{assembly}'] - Database '{sub_database}' - Process file '{os.path.basename(database_file)}' - Parquet file already exists""")

                else:

                    # Query
                    db_copy = duckdb.connect(config={"threads":nb_data_files})
                    db_copy.query(query)
                    db_copy.close()

                    # # Split parquet to max size
                    if partition_size:
                        log.debug(f"""Download dbNSFP ['{assembly}'] - Database '{sub_database}' - Process file '{os.path.basename(database_file)}' - Split parquet files...""")
                        for file in os.listdir(output_parquet):
                            file_path = os.path.join(output_parquet,file)
                            file_size = os.stat(file_path).st_size
                            file_base = file.replace(".parquet","")
                            log.debug(f"""Download dbNSFP ['{assembly}'] - Database '{sub_database}' - Process file '{os.path.basename(database_file)}' - Split parquet files ['{file}']...""")
                            if file_size > partition_size:
                                schema = pq.ParquetFile(file_path).schema_arrow
                                dd.read_parquet(file_path).repartition(partition_size=partition_size*8).to_parquet(output_parquet, schema=schema, name_function=lambda x: f"{file_base}-{x}.parquet")
                                os.remove(file_path)
                            else:
                                new_file_name = os.path.join(output_parquet,file.replace(".parquet","-0.parquet"))
                                os.rename(file_path, new_file_name)


                    # Log
                    log.debug(f"""Download dbNSFP ['{assembly}'] - Database '{sub_database}' - Process file '{os.path.basename(database_file)}' - Write parquet file - done.""")

            # Header
            header_file = f"{output}.hdr"
            if not os.path.exists(header_file):
                if not readme_annotations_description:
                    readme_annotations_description = get_annotation_description(dbnsfp_readme_dest)
                if get_header(header_file=header_file, dbnsfp_readme_dest=dbnsfp_readme_dest, columns_structure=columns_structure, annotations=columns_annotations.keys(), readme_annotations_description=readme_annotations_description):
                    log.debug(f"Download dbNSFP ['{assembly}'] - Database '{sub_database}' - Write header")
                else:
                    log.error(f"Download dbNSFP ['{assembly}'] - Database '{sub_database}' - Write header failed")
                    raise ValueError(f"Download dbNSFP ['{assembly}'] - Database '{sub_database}' - Write header failed")
                

    if generate_sub_databases:

        log.info(f"""Download dbNSFP ['{assembly}']""")
        log.info(f"""Download dbNSFP ['{assembly}'] - Parquet folders generation""")

        for assembly in assemblies:

            # Find database structure
            if not database_files:
                database_files = get_database_files(f"{dbnsfp_zip_dest_folder}/dbNSFP{dbnsfp_release}_variant.chr*.gz")
            if not columns_structure:
                columns_structure = get_columns_structure(database_file=database_files[0], sample_size=10, threads=threads)

            output_prefix = f"{dbnsfp_folder}/{assembly}/dbNSFP{dbnsfp_release}"
            parquet_all_annotation = f"{output_prefix}.ALL.partition.parquet/*/*.parquet"

            sub_databases_structure = {}

            # Find sub databases
            for column in columns_structure:
                sub_database = column.split("_")[0]
                if sub_database not in sub_databases_structure:
                    sub_databases_structure[sub_database] = {}
                sub_databases_structure[sub_database][column] = columns_structure[column]

            # Chromosomes
            chromosomes = []

            # for each sub database
            for sub_database in sorted(set(sub_databases_structure), key=str.casefold, reverse=False):
                
                # Clauses
                columns_clauses = get_columns_select_clause(columns_structure=columns_structure, assembly=assembly, for_parquet=True, sub_database=sub_database, null_if_no_annotation=True, print_log=False)
                columns_select_clause = columns_clauses.get("select")
                columns_where_clause = columns_clauses.get("where")
                columns_annotations = columns_clauses.get("annotations")

                # Clean sub database name
                sub_database_name = clean_name(sub_database)

                # Parquet folder
                parquet_sub_database_annotation = f"{output_prefix}.{sub_database_name}.partition.parquet"

                if columns_select_clause:

                    # Parquet Folder
                    if os.path.exists(f"{parquet_sub_database_annotation}"):
                        log.debug(f"""Download dbNSFP ['{assembly}'] - Database '{sub_database}' - Parquet folder ['{sub_database_name}'] already exists""")
                        parquet_partition_already_generated_list.append(sub_database)
                    else:
                        log.info(f"""Download dbNSFP ['{assembly}'] - Database '{sub_database}' - Parquet folder ['{sub_database_name}'] generation...""")
                        
                        db_copy = duckdb.connect(config={"threads":nb_data_files})

                        if not chromosomes:
                            query_chromosomes = f"""
                                SELECT distinct "#CHROM"
                                    FROM read_parquet('{parquet_all_annotation}')
                            """
                            res = db_copy.query(query_chromosomes)
                            chromosomes = sorted(list(res.df()["#CHROM"]))

                        # Generate folder by chromosome
                        for chromosome in chromosomes:
                            parquet_sub_database_annotation_chromosome = os.path.join(parquet_sub_database_annotation, f"#CHROM={chromosome}")
                            Path(parquet_sub_database_annotation_chromosome).mkdir(parents=True, exist_ok=True)
                            query_copy = f"""
                                COPY (
                                    SELECT {columns_select_clause}
                                        FROM read_parquet('{parquet_all_annotation}')
                                        WHERE "#CHROM" in ('{chromosome}')
                                            AND ({columns_where_clause})
                                    )
                                TO '{parquet_sub_database_annotation_chromosome}' WITH (FORMAT PARQUET, PER_THREAD_OUTPUT TRUE)
                                """

                            db_copy.query(query_copy)

                            # Split parquet to max size
                            if partition_size:
                                #log.debug(f"""Download dbNSFP ['{assembly}'] - Database '{sub_database}' - Parquet folder ['{sub_database_name}'] generation - Split parquet files...""")
                                for file in os.listdir(parquet_sub_database_annotation_chromosome):
                                    log.debug(f"""Download dbNSFP ['{assembly}'] - Database '{sub_database}' - Parquet folder ['{sub_database_name}'] generation - Split parquet files ['{file}']...""")
                                    file_path = os.path.join(parquet_sub_database_annotation_chromosome,file)
                                    file_size = os.stat(file_path).st_size
                                    file_base = file.replace(".parquet","")
                                    if file_size > partition_size:
                                        schema = pq.ParquetFile(file_path).schema_arrow
                                        dd.read_parquet(file_path).repartition(partition_size=partition_size*8).to_parquet(parquet_sub_database_annotation_chromosome, schema=schema, name_function=lambda x: f"{file_base}-{x}.parquet")
                                        os.remove(file_path)
                                    else:
                                        new_file_name = os.path.join(parquet_sub_database_annotation_chromosome,file.replace(".parquet","-0.parquet"))
                                        os.rename(file_path, new_file_name)

                        
                        db_copy.close()
                        parquet_partition_generated_list.append(sub_database)

                        log.debug(f"""Download dbNSFP ['{assembly}'] - Database '{sub_database}' - Parquet folder ['{sub_database_name}'] generation - done.""")

                    # Header
                    header_file = f"{parquet_sub_database_annotation}.hdr"
                    if not os.path.exists(header_file):
                        if not readme_annotations_description:
                            readme_annotations_description = get_annotation_description(dbnsfp_readme_dest)
                        if get_header(header_file=header_file, dbnsfp_readme_dest=dbnsfp_readme_dest, columns_structure=columns_structure, annotations=columns_annotations.keys(), readme_annotations_description=readme_annotations_description):
                            log.debug(f"Download dbNSFP ['{assembly}'] - Database '{sub_database}' - Parquet folder ['{sub_database_name}'] generation - Write header")
                        else:
                            log.error(f"Download dbNSFP ['{assembly}'] - Database '{sub_database}' - Parquet folder ['{sub_database_name}'] generation - Write header")
                            raise ValueError(f"Download dbNSFP ['{assembly}'] - Database '{sub_database}' - Parquet folder ['{sub_database_name}'] generation - Write header")

    
    # Full DB
    if generate_parquet_file or generate_vcf_file:

        # Init
        genomes_sizes = {}

        for assembly in assemblies:

            log.info(f"Download dbNSFP ['{assembly}']")
            log.info(f"Download dbNSFP ['{assembly}'] - Parquet/VCF files generation")

            output_prefix = f"{dbnsfp_folder}/{assembly}/dbNSFP{dbnsfp_release}"
            partition_parquet_folder_all = f"{output_prefix}.ALL.partition.parquet"
            partition_parquet_folders = sorted(glob.glob(f"{output_prefix}.*.partition.parquet"), key=os.path.basename, reverse=True)
            if partition_parquet_folders:
                partition_parquet_folders.remove(partition_parquet_folder_all)

            for partition_parquet_folder in [partition_parquet_folder_all] + sorted(partition_parquet_folders, key=str.casefold, reverse=False):

                # sub database
                sub_database = os.path.basename(partition_parquet_folder).split(".")[-3]

                # Input Parquet folder
                input_partition_parquet = f"{output_prefix}.{sub_database}.partition.parquet"
                input_partition_parquet_header = f"{input_partition_parquet}.hdr"

                # Input Parquet files
                input_parquet_files = sorted(glob.glob(f"{input_partition_parquet}/*/*parquet"), key=os.path.normcase, reverse=False)

                # Parquet file
                if generate_parquet_file:

                    # Output
                    output_parquet = f"{output_prefix}.{sub_database}.parquet"
                    output_header_file = f"{output_parquet}.hdr"

                    # Parquet file
                    if not os.path.exists(output_parquet):

                        if input_parquet_files:

                            # Log
                            log.info(f"""Download dbNSFP ['{assembly}'] - Database '{sub_database}' - Parquet file...""")

                            # Generate Parquet files
                            #log.debug(pq.ParquetFile(input_parquet_files[0]).schema_arrow)
                            with pq.ParquetWriter(output_parquet, schema=pq.ParquetFile(input_parquet_files[0]).schema_arrow) as writer:
                                chromosome_folder_previous = ""
                                for parquet_file in input_parquet_files:
                                    chromosome_folder = f'{os.path.basename(os.path.dirname(parquet_file))}'
                                    if chromosome_folder != chromosome_folder_previous:
                                        chromosome_folder_previous = chromosome_folder
                                        log.debug(f"""Download dbNSFP ['{assembly}'] - Database '{sub_database}' - Parquet file - Process '{os.path.basename(os.path.dirname(parquet_file))}/*' files...""")
                                    log.debug(f"""Download dbNSFP ['{assembly}'] - Database '{sub_database}' - Parquet file - Process '{os.path.basename(os.path.dirname(parquet_file))}/{os.path.basename(parquet_file)}' file...""")
                                    #log.debug(pq.ParquetFile(parquet_file).schema_arrow)
                                    writer.write_table(pq.read_table(parquet_file))
                            parquet_generated_list.append(sub_database)

                        else:
                            log.error(f"""Download dbNSFP ['{assembly}'] - Database '{sub_database}' - Parquet file - No parquet files found""")
                            raise ValueError(f"""Download dbNSFP ['{assembly}'] - Database '{sub_database}' - Parquet file - No parquet files found""")

                    else:

                        # Log
                        log.debug(f"""Download dbNSFP ['{assembly}'] - Database '{sub_database}' - Parquet file already generated""")

                        parquet_already_generated_list.append(sub_database)


                    # Header
                    if not os.path.exists(output_header_file):
                        if os.path.exists(input_partition_parquet_header):
                            shutil.copy(input_partition_parquet_header, output_header_file)
                            log.debug(f"Download dbNSFP ['{assembly}'] - Database '{sub_database}' - Parquet file - Write header")
                        else:
                            log.error(f"Download dbNSFP ['{assembly}'] - Database '{sub_database}' - Parquet file - Write header failed")
                            raise ValueError(f"Download dbNSFP ['{assembly}'] - Database '{sub_database}' - Parquet file - Write header failed")

                if generate_vcf_file:
                    
                    if sub_database in ['ALL'] and not generate_vcf_file_all:
                        continue

                    # Output VCF
                    output_vcf = f"{output_prefix}.{sub_database}.vcf.gz"
                    header_file = f"{output_vcf}.hdr"

                    # Chromosome sizes
                    if not genomes_sizes:
                        genomes_sizes_file = os.path.join(genomes_folder, assembly, f"{assembly}.fa.sizes")
                        if os.path.exists(genomes_sizes_file):
                            with open(genomes_sizes_file, "r") as f:
                                for line in f:
                                    genomes_sizes[line.split("\t")[0]] = line.split("\t")[1].strip()
                        else:
                            genomes_sizes = None

                    # VCF file
                    if not os.path.exists(output_vcf):

                        if input_parquet_files:

                            # Init
                            list_for_vcf = []

                            # Log
                            log.info(f"""Download dbNSFP ['{assembly}'] - Database '{sub_database}' - VCF file...""")

                            # Find database structure
                            if not database_files:
                                database_files = get_database_files(f"{dbnsfp_zip_dest_folder}/dbNSFP{dbnsfp_release}_variant.chr*.gz")
                            if not columns_structure:
                                columns_structure = get_columns_structure(database_file=database_files[0], sample_size=10, threads=threads)
                            columns_clauses = get_columns_select_clause(columns_structure=columns_structure, assembly=assembly, for_parquet=True, sub_database=sub_database, null_if_no_annotation=True, print_log=False)
                            columns_annotations = columns_clauses.get("annotations")

                            # Generate columns concat for INFO column
                            column_INFO = ""
                            i = 0
                            for column_annotation in columns_annotations:
                                if i < 1000000:
                                    column_INFO += f"""
                                    CASE
                                        WHEN "{column_annotation}" IS NOT NULL
                                        THEN '{column_annotation}=' || replace("{column_annotation}", ';', ',') || ';'
                                        ELSE ''
                                    END ||
                                    """
                                i += 1

                            if column_INFO:
                                
                                # Generate VCF files within tmp folder
                                with TemporaryDirectory(dir='/tmp') as tmp_dir:
                                    
                                    # Init
                                    file_num = -1
                                
                                    # Create duckDB connexion
                                    db_copy = duckdb.connect(config={"threads":threads})

                                    # Set max expression depth for big sub databases (e.g. gnomAD, ALL)
                                    db_copy.execute("SET max_expression_depth TO 10000")

                                    # Log
                                    log.debug(f"""Download dbNSFP ['{assembly}'] - Database '{sub_database}' - VCF file - Check chomosomes and number of variants...""")

                                    # Check chomosomes and total number of variants
                                    query_chromosomes = f"""
                                                SELECT distinct "#CHROM" AS chromosome, count(*) AS nb_variants
                                                FROM read_parquet('{input_partition_parquet}/*/*parquet')
                                                GROUP BY "#CHROM"
                                            """
                                    chromosomes = sorted(list(db_copy.query(query_chromosomes).df()["chromosome"]))
                                    nb_variants = sum(db_copy.query(query_chromosomes).df()["nb_variants"])
                                    
                                    # Log
                                    log.debug(f"""Download dbNSFP ['{assembly}'] - Database '{sub_database}' - VCF file - Found {len(chromosomes)} chromosomes""")
                                    log.debug(f"""Download dbNSFP ['{assembly}'] - Database '{sub_database}' - VCF file - Found {nb_variants} variants""")

                                    log.debug(f"""Download dbNSFP ['{assembly}'] - Database '{sub_database}' - VCF file - Generate files...""")
                                    
                                    for parquet_file in input_parquet_files:
                                        
                                        # Check number of variant for parquet file
                                        query_nb_variants = f"""
                                                    SELECT count(*) AS nb_variants
                                                    FROM read_parquet('{parquet_file}')
                                                """
                                        parquet_file_nb_variants = sum(db_copy.query(query_nb_variants).df()["nb_variants"])

                                        # Log
                                        log.debug(f"""Download dbNSFP ['{assembly}'] - Database '{sub_database}' - VCF file - Process '{os.path.basename(os.path.dirname(parquet_file))}/{os.path.basename(parquet_file)}' file [{parquet_file_nb_variants} variants]...""")

                                        # Query for parquet file
                                        query_copy = f"""
                                                    SELECT "#CHROM", POS, '.' AS ID, REF, ALT, 0 AS QUAL, 'PASS' AS FILTER, regexp_replace(replace({column_INFO} '','"', ''), ';$', '') INFO
                                                    FROM read_parquet('{parquet_file}')
                                                """
                                        # CSV generation with PyArrow
                                        if parquet_file_nb_variants:

                                            # Options
                                            batch_size = 102400
                                            write_options = csv.WriteOptions(include_header=False, delimiter="\t", quoting_style="none", batch_size=batch_size)

                                            # Query
                                            res = db_copy.execute(query_copy)

                                            # Fetch rows
                                            record_batch_reader = res.fetch_record_batch(rows_per_batch=nb_variants)
                                            chunk = record_batch_reader.read_next_batch()
                                            if len(chunk) > 0:
                                                i += 1
                                                # VCF file name
                                                vcf_file = os.path.join(tmp_dir,f'variants.{file_num}.tsv')
                                                # Write VCF file
                                                csv.write_csv(chunk, vcf_file, write_options=write_options)
                                                # Append to list of VCF file
                                                list_for_vcf.append(vcf_file)

                                    # Log
                                    log.debug(f"""Download dbNSFP ['{assembly}'] - Database '{sub_database}' - VCF file - {len(list_for_vcf)} files generated""")

                                    # Header
                                    if not os.path.exists(header_file):

                                        # Init
                                        header_list = []

                                        # Temporary header file for INFO fields
                                        header_file_tmp = os.path.join(tmp_dir,f"header.hdr")
                                        if os.path.exists(input_partition_parquet_header):
                                            shutil.copy(input_partition_parquet_header, header_file_tmp)
                                            log.debug(f"Download dbNSFP ['{assembly}'] - Database '{sub_database}' - VCF file - Write header")
                                        else:
                                            log.error(f"Download dbNSFP ['{assembly}'] - Database '{sub_database}' - VCF file - Write header failed")
                                            raise ValueError(f"Download dbNSFP ['{assembly}'] - Database '{sub_database}' - VCF file - Write header failed")
                                        
                                        # Read header file to dict
                                        with open(header_file_tmp, "r") as f:
                                            for line in f:
                                                header_list.append(line.strip())

                                        # Remove last line with #CHROM
                                        header_list = header_list[:-1]

                                        # Append FILTER PASS
                                        header_list.append('##FILTER=<ID=PASS,Description="All filters passed">')

                                        # Add chromosomes
                                        for chromosome in chromosomes:
                                            header_list.append(f"##contig=<ID={chromosome},length={genomes_sizes.get(chromosome,0)},assembly={assembly}>")

                                        # Add last line with #CHROM
                                        header_list.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO")

                                        # Write header
                                        with open(header_file, 'w') as fp:
                                            for item in header_list:
                                                # write each item on a new line
                                                fp.write("%s\n" % item)

                                        # Insert header to lists of files to concat
                                        list_for_vcf.insert(0, header_file)

                                    # Log
                                    log.debug(f"Download dbNSFP ['{assembly}'] - Database '{sub_database}' - VCF file - Concate and Compress files...")
                                    concat_and_compress_files(input_files=list_for_vcf, output_file=output_vcf)

                                    db_copy.close()

                                    vcf_generated_list.append(sub_database)

                            else:
                                log.error(f"""Download dbNSFP ['{assembly}'] - Database '{sub_database}' - VCF file - No INFO columns found""")
                                raise ValueError(f"""Download dbNSFP ['{assembly}'] - Database '{sub_database}' - VCF file - No INFO columns found""")

                        else:
                            log.error(f"""Download dbNSFP ['{assembly}'] - Database '{sub_database}' - VCF file - No parquet files found""")
                            raise ValueError(f"""Download dbNSFP ['{assembly}'] - Database '{sub_database}' - VCF file - No parquet files found""")
                                
                        # if subdatabase in ['ALL']:
                        #     exit()
                    
                    else:
                        vcf_already_generated_list.append(sub_database)

    # Log
    if parquet_partition_generated_list:
        log.info(f"""Download dbNSFP ['{assembly}'] - {len(set(parquet_partition_generated_list))} Parquet folders generated {list(set(parquet_partition_generated_list))}""")
    if parquet_partition_already_generated_list:
        log.info(f"""Download dbNSFP ['{assembly}'] - {len(set(parquet_partition_already_generated_list))} Parquet folders already generated {list(set(parquet_partition_already_generated_list))} """)
    if parquet_generated_list:
        log.info(f"""Download dbNSFP ['{assembly}'] - {len(set(parquet_generated_list))} Parquet files generated {list(set(parquet_generated_list))}""")
    if parquet_already_generated_list:
        log.info(f"""Download dbNSFP ['{assembly}'] - {len(set(parquet_already_generated_list))} Parquet files already generated {list(set(parquet_already_generated_list))} """)
    if vcf_generated_list:
        log.info(f"""Download dbNSFP ['{assembly}'] - {len(set(vcf_generated_list))} VCF files generated {list(set(vcf_generated_list))}""")
    if vcf_already_generated_list:
        log.info(f"""Download dbNSFP ['{assembly}'] - {len(set(vcf_already_generated_list))} VCF files already generated {list(set(vcf_already_generated_list))} """)


    # Close duckdb connexion
    db.close()

    return True

