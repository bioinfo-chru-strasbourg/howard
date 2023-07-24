#!/usr/bin/env python

import argparse
from functools import partial
import itertools
import multiprocessing
import os
import subprocess
import pyarrow.parquet as pq
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


from howard.commons import *


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
                contig_regex=args.download_genomes_contig_regex
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

