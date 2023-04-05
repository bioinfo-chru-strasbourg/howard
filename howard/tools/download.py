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



def databases_download(args) -> None:

    log.debug(f"Args {args}")

    # Assembly
    assemblies = [value for val in args.assembly for value in val.split(',')]

    # Annovar
    if args.download_annovar:
        log.debug(f"Download Annovar databases")
        if args.download_annovar_files:
            files = [value for val in args.download_annovar_files for value in val.split(',')]
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


def databases_download_annovar(folder:str = None, files:list = None, assemblies:list = ["hg19"], annovar_url:str = "http://www.openbioinformatics.org/annovar/download/") -> None:

    log.info(f"Download Annovar databases {assemblies}")

    # Minimum files to download
    #files_minimum = ["refGene.txt.gz","refGene.txt.idx.gz","refGeneMrna.fa.gz"]
    files_minimum = ["refGene*"]

    for assembly in assemblies:

        log.debug(f"Download Annovar for assembly '{assembly}'")

        if folder:
            if not os.path.exists(folder):
                os.makedirs(folder)
        
        log.debug(f"Download Annovar databases in folder '{folder}'")

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
        avdblist_url_file = f"{annovar_url}{avdblist_file}"
        avdblist_folder_file = f"{folder}/{avdblist_file}"
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
                files_in_folder = glob.glob(os.path.join(folder, file))
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
                    file_path = os.path.join(folder, file)

                    log.debug(f"Download file {file} from {file_url} to {file_path}...")
                    download_file(file_url, file_path)
                    log.debug(f"Extract file {file} to {folder}...")
                    extract_file(file_path)
            else:
                log.info(f"Download Annovar databases {[assembly]} already exists")


def databases_download_snpeff(folder:str = None, assemblies:list = ["hg19"], config:dict = {}) -> None:

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

            log.info(f"Database snpEff databases {[assembly]} already exists")

