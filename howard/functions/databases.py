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



def generate_databases_param(args:argparse, assemblies:list = []):
    """
    This function generates database parameters based on specified arguments and assemblies.
    
    :param args: The `args` parameter in the `generate_databases_param` function is expected to be an
    instance of the `argparse` module, which is commonly used for parsing command-line arguments. This
    parameter is used to retrieve various arguments and options provided by the user when running the
    script or program
    :type args: argparse
    :param assemblies: The `assemblies` parameter is a list containing the assemblies for which
    databases will be generated. The function `generate_databases_param` takes various arguments using
    the `argparse` module and generates database parameters based on these inputs. If the
    `generate_param` argument is provided and set to True,
    :type assemblies: list
    :return: None
    """

    # Param
    if "generate_param" in args and args.generate_param:
        log.debug(f"Generate param config")

        if isinstance(args.generate_param, str):
            generate_param_file = args.generate_param
        else:
            generate_param_file = args.generate_param.name

        if "generate_param_description" in args and args.generate_param_description:
            if isinstance(args.generate_param_description, str):
                generate_param_description_file = args.generate_param_description
            else:
                generate_param_description_file = args.generate_param_description.name
        else:
            generate_param_description_file = None

        # Check assembly
        if len(assemblies) == 1:
            assembly = assemblies[0]
        else:
            log.debug(f"Error: choose one uniq assembly '{assemblies}'")
            raise ValueError(f"Error: choose one uniq assembly '{assemblies}'")
        
        # Databases releases
        if args.generate_param_releases:
            generate_param_releases = [value for value in args.generate_param_releases.split(',')]
        else:
            generate_param_releases = []

        # Databases formats
        if args.generate_param_formats:
            generate_param_formats = [value for value in args.generate_param_formats.split(',')]
        else:
            generate_param_formats = []

        # databases infos
        databases_infos_dict = databases_infos(
            database_folder_releases=generate_param_releases,
            assembly = assembly,
            database_formats=generate_param_formats,
            config=args.config
            )
        databases_param_stats = databases_param(
            databases_infos_dict=databases_infos_dict,
            bcftools_preference=args.generate_param_bcftools,
            output=generate_param_file,
            output_description=generate_param_description_file
            )

        # Show param (3 levels of param)
        for databases_param_stat in databases_param_stats:
            if not isinstance(databases_param_stats[databases_param_stat], dict):
                log.info(f"{databases_param_stat}: {databases_param_stats[databases_param_stat]}")
            else:
                log.info(f"{databases_param_stat}")
                for databases_param_substat in databases_param_stats[databases_param_stat]:
                    if not isinstance(databases_param_stats[databases_param_stat][databases_param_substat], dict):
                        log.info(f"   {databases_param_substat}: {databases_param_stats[databases_param_stat][databases_param_substat]}")
                    else:
                        log.info(f"   {databases_param_substat}")
                        for databases_param_subsubstat in databases_param_stats[databases_param_stat][databases_param_substat]:
                            log.info(f"      {databases_param_subsubstat}: {databases_param_stats[databases_param_stat][databases_param_substat][databases_param_subsubstat]}")

        return None


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


def databases_infos(database_folders:list = [], database_folder_releases:list = ["current"], assembly:str = "hg19", database_formats:list = None, config:dict = {}) -> dict:
    """
    The `databases_infos` function scans database folders and retrieves information about the databases
    found, including their folder, release, assembly, subdatabase, format, header, and parameters.
    
    :param database_folders: A list of folders where the databases are located
    :type database_folders: list
    :param database_folder_releases: A list of specific releases of the database folders to include in
    the search. If None, all releases will be included
    :type database_folder_releases: list
    :param assembly: The `assembly` parameter is a string that specifies the assembly version of the
    databases to be searched. It is used to filter the databases based on their assembly version. The
    default value is "hg19", defaults to hg19
    :type assembly: str (optional)
    :param database_formats: The `database_formats` parameter is a list that specifies the formats of
    the databases to include in the results. If this parameter is not provided or is set to `None`, all
    database formats will be included
    :type database_formats: list
    :param config: The `config` parameter is a dictionary that contains configuration settings for the
    function. It has the following structure:
    :type config: dict
    :return: The `databases_infos` function returns a dictionary containing information about the
    databases found in the specified database folders. The keys of the dictionary are the paths to the
    database files, and the values are dictionaries containing the following information: folder,
    release, assembly, subdatabase, format, header, and parameters.
    """

    # Init Database object
    from howard.objects.database import Database
    vcfdata_obj = Database()

    # Init results
    infos = {}

    # Check databases folders from config
    if not database_folders:
        config_database_folders = config.get("folders",{}).get("databases",{}).get("annotations",[])
        for config_database_folder in config_database_folders:
            config_database_folder_dirname = os.path.dirname(re.sub(r"/*$", "", config_database_folder))
            log.debug(f"config_database_folder: {config_database_folder}")
            log.debug(f"   config_database_folder_dirname: {config_database_folder_dirname}")
            database_folders.append(config_database_folder_dirname)

    # Find all databases folders
    if not database_folders:
        database_folders = []
        default_database_folder = DEFAULT_DATABASE_FOLDER
        database_folders_list = os.listdir(default_database_folder)
        for database_folder in database_folders_list:
            if os.path.isdir(os.path.join(default_database_folder, database_folder)):
                database_folders.append(os.path.join(default_database_folder, database_folder))

    log.debug(f"database_folders: {database_folders}")

    # Scan database folders
    for database_folder in database_folders:
        
        # Database folder name
        database_folder_name = os.path.basename(database_folder)

        # Scan all folder with allowed structure in database folder (<release>/<assembly>)
        for root, dirs, files in os.walk(f"{database_folder}", followlinks=True):
            
            # database structure and level
            database_subfolder = root.replace(database_folder, '')
            database_subfolder_structure = database_subfolder.split("/")[1:]
            level = database_subfolder.count(os.sep)

            # Level as correct structure
            if level == 2:

                # If "(<release>/<assembly>" is a folder
                if os.path.isdir(root):

                    # Database infos (release and assembly)
                    database_release = database_subfolder_structure[0]
                    database_assembly = database_subfolder_structure[1]

                    # Filter with input release and assembly
                    if (not database_folder_releases or database_release in database_folder_releases) and (database_assembly == assembly):

                        # List files within the folder
                        list_of_files = glob.glob(os.path.join(root,"**"), recursive=True)

                        # List file as a database (to check)
                        for database in list_of_files:

                            # Check if file exists (ERROR woth glob.glob???)
                            database = full_path(database)
                            if os.path.exists(database):

                                # Format of the database
                                database_format = vcfdata_obj.get_format(database=database)

                                # Check if it is a database (header exists as file or within vcf)
                                # Filter with input database format
                                if (os.path.isfile(f"{database}.hdr") or (os.path.isfile(database) and database.endswith(("vcf","vcf.gz","vcf.bgz")))) and (not database_formats or database_format in database_formats):

                                    # Database subfolder (allow files in subfolder, such as database versions)
                                    database_subfolder = os.path.dirname(database.replace(f"{root}/", ''))

                                    # Database header
                                    database_header = vcfdata_obj.get_header(database=database).infos
                                    
                                    # Init
                                    header_infos_dict = {}
                                    header_infos_dict_param = {}

                                    # List header infos 
                                    for info in list(database_header):

                                        # Init
                                        header_infos_dict[info] = {}
                                        header_infos_dict_param[info] = info

                                        # ID
                                        header_infos_dict[info]["id"] = info

                                        # num
                                        if database_header[info].num in GENOTYPE_MAP.keys():
                                            header_infos_dict[info]["Number"] = GENOTYPE_MAP.get(database_header[info].num)
                                        else:
                                            header_infos_dict[info]["Number"] = database_header[info].num

                                        # type
                                        if database_header[info].type:
                                            header_infos_dict[info]["Type"] = database_header[info].type
                                        else:
                                            header_infos_dict[info]["Type"] = "."

                                        # desc
                                        if database_header[info].desc != None:
                                            header_infos_dict[info]["Description"] = database_header[info].desc
                                        else:
                                            header_infos_dict[info]["Description"] = ""

                                    # Add database to dict
                                    if header_infos_dict:
                                        infos[database] = {
                                            "folder": database_folder_name,
                                            "release": database_release,
                                            "assembly": database_assembly,
                                            "subdatabase": database_subfolder,
                                            "format": database_format,
                                            "header": header_infos_dict,
                                            "param": header_infos_dict_param
                                        }

                                        # Log
                                        log.debug(f"---")
                                        log.debug(f"database: {database}")
                                        log.debug(json.dumps({
                                            "folder": database_folder_name,
                                            "release": database_release,
                                            "assembly": database_assembly,
                                            "subdatabase": database_subfolder,
                                            "format": database_format,
                                            "header": len(header_infos_dict),
                                            "param": len(header_infos_dict_param)
                                        }, indent=4))
                        
    return infos


def databases_param(databases_infos_dict:dict, output:str = None, output_description:str = None, bcftools_preference:bool = False) -> dict:
    """
    The `databases_param` function takes in a dictionary of database information, an optional output
    file path, and a boolean flag for bcftools preference, and returns a dictionary containing the
    parameters for parquet and bcftools annotations.
    
    :param databases_infos_dict: A dictionary containing information about databases. Each key in the
    dictionary represents the name of a database, and the corresponding value is another dictionary
    containing information about the database, such as its format and parameters
    :type databases_infos_dict: dict
    :param output: The `output` parameter is a string that specifies the path and filename of the output
    file where the generated JSON object will be written. If this parameter is not provided or is set to
    `None`, the JSON object will not be written to a file
    :type output: str
    :param output_description: The `output_description` parameter is a string that specifies the path
    and filename of the output file where the description of the databases will be written. If this
    parameter is not provided or is set to `None`, the description will not be written to a file
    :type output_description: str
    :param bcftools_preference: The `bcftools_preference` parameter is a boolean flag that determines
    whether to prioritize databases in the BCFTOOLS format. If `bcftools_preference` is set to `True`,
    databases in the BCFTOOLS format will be given priority over other formats. If `bcftools, defaults
    to False
    :type bcftools_preference: bool (optional)
    :return: The function `databases_param` returns a dictionary object named "param_stats_show".
    """
    
    # Full path
    output = full_path(output)
    output_description = full_path(output_description)

    # Init
    param = {
                "annotation": {}
    }
    param_parquet = {}
    param_bcftools = {}
    param_stats = {
        "assembly": [],
        "databases": [],
        "releases": [],
        "formats": [],
        "fields": []
    }
    param_stats_show = {}

    # If None (not empty) dict in input, check by default
    if databases_infos_dict == None:
        databases_infos_dict = databases_infos()

    # List databases
    for database_infos in databases_infos_dict:

        database_format = databases_infos_dict[database_infos].get("format", None)
        
        # If preference for bcftools
        if bcftools_preference and database_format in BCFTOOLS_FORMAT:
            param_bcftools[database_infos] = databases_infos_dict[database_infos].get("param")
        else:
            param_parquet[database_infos] = databases_infos_dict[database_infos].get("param")
        
        # Stats
        param_stats["assembly"].append(databases_infos_dict[database_infos].get("assembly"))
        param_stats["databases"].append(database_infos)
        param_stats["releases"].append(databases_infos_dict[database_infos].get("release"))
        param_stats["formats"].append(databases_infos_dict[database_infos].get("format"))
        param_stats["fields"] += list(databases_infos_dict[database_infos].get("header").keys())

    # Stats
    param_stats_show["Assembly"] = param_stats["assembly"][0]
    param_stats_show["Number of databases"] = len(set(param_stats["databases"])) 
    param_stats_show["Releases of databases"] = list(set(param_stats["releases"]))
    param_stats_show["Formats of databases"] = list(set(param_stats["formats"]))
    param_stats_show["Number of Annotation fields"] = len(set(param_stats["fields"])) 

    # Add parquet annotations to param
    if param_parquet:
        param["annotation"]["parquet"] = {
            "annotations": param_parquet
        }

    # Add bcftools annotations to param
    if bcftools_preference and param_bcftools:
        param["annotation"]["bcftools"] = {
            "annotations": param_bcftools
        }

    # Output
    if output:

        # Create json object
        json_object = json.dumps(param, indent=4)
        
        # Create folder if not exists
        if not os.path.exists(os.path.dirname(output)):
            Path(os.path.dirname(output)).mkdir(parents=True, exist_ok=True)

        # Write output file
        with open(output, "w") as outfile:
            outfile.write(json_object)

    # Output description
    if output_description:

        # Create json object
        json_object = json.dumps(databases_infos_dict, indent=4)
        
        # Create folder if not exists
        if not os.path.exists(os.path.dirname(output_description)):
            Path(os.path.dirname(output_description)).mkdir(parents=True, exist_ok=True)

        # Write output file
        with open(output_description, "w") as outfile:
            outfile.write(json_object)

    return param_stats_show


def databases_download_annovar(folder:str = None, files:list = None, assemblies:list = ["hg19"], annovar_url:str = "http://www.openbioinformatics.org/annovar/download", threads:int = 1) -> None:
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
    :param threads: The "threads" parameter specifies the number of threads (parallel processes) to use
    for download and extract/uncompress files. Default: 1
    :type threads: int (optional)
    """

    log.info(f"Download Annovar databases {assemblies}")

    # Full path
    folder = full_path(folder)

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
        download_file(avdblist_url_file, avdblist_folder_file, threads=threads)

        if not os.path.exists(avdblist_folder_file):
            log.error(f"Download list of Annovar files from Annovar URL: {avdblist_url_file}")
            log.error(f"Annovar Database list: {avdblist_folder_file}")
            raise FileNotFoundError(f"Annovar Database list: {avdblist_folder_file}")
        
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
                    download_file(file_url, file_path, threads=threads)
                    log.debug(f"Extract file {file} to {folder_assembly}...")
                    extract_file(file_path=file_path, threads=threads)
            else:
                log.info(f"Download Annovar databases {[assembly]} - already exists")


def databases_download_snpeff(folder:str = None, assemblies:list = ["hg19"], config:dict = {}, threads:int = 1) -> None:
    """
    The `databases_download_snpeff` function downloads and extracts snpEff databases for specified
    genome assemblies.
    
    :param folder: The `folder` parameter is a string that specifies the folder where the snpEff
    databases will be downloaded and stored. If the folder does not exist, it will be created
    :type folder: str
    :param assemblies: The `assemblies` parameter is a list of genome assemblies for which the snpEff
    databases need to be downloaded. It specifies the genome assemblies for which you want to download
    the snpEff databases. For example, if you want to download the snpEff databases for the human genome
    assembly hg
    :type assemblies: list
    :param config: The `config` parameter is a dictionary that contains information about the tools and
    their configurations. It is used to retrieve the path to the Java binary and the path to the snpEff
    binary
    :type config: dict
    :param threads: The `threads` parameter specifies the number of threads to be used for downloading
    the snpEff databases. It determines the parallelism of the download process, allowing multiple files
    to be downloaded simultaneously, defaults to 1
    :type threads: int (optional)
    """

    log.info(f"Download snpEff databases {assemblies}")

    # Full path
    folder = full_path(folder)

    # Java bin and snpEff jar
    snpeff_bin = get_bin(tool="snpeff", bin="snpEff.jar", bin_type="jar", config=config, default_folder=f"{DEFAULT_TOOLS_FOLDER}/snpeff")
    java_bin = get_bin(tool="java", bin="java", bin_type="bin", config=config, default_folder="/usr/bin")

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
                log.info(f"Download snpEff databases {assemblies} - list of databases...")
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

            if not snpeff_list_databases:
                log.error(f"Download snpEff databases {[assembly]} - list of databases empty - check file '{snpeff_databases_list_path}'")
                raise ValueError(f"Download snpEff databases {[assembly]} - list of databases empty - check file '{snpeff_databases_list_path}'")

            # Start download
            log.info(f"Download snpEff databases {[assembly]} - downloading...")
            # Try to download files
            file_path = None
            for file_url in snpeff_list_databases.get(assembly,[]):
                # File to be downloaded
                file_path = os.path.join(folder, os.path.basename(file_url))
                # Check if file already downloaded
                if os.path.exists(file_path):
                    break
                else:
                    # try to download
                    try:
                        log.debug(f"Download snpEff '{file_url}'...")
                        # Download file
                        if download_file(file_url, file_path, threads=threads):
                            break
                        else:
                            log.error(f"Download snpEff '{file_url}' failed")
                    # If fail, just pass to next url
                    except:
                        log.error(f"Download snpEff '{file_url}' failed")

            # If download file exists
            if file_path is not None and os.path.exists(file_path):
                log.info(f"Download snpEff databases {[assembly]} - extracting...")
                log.debug(f"Extract file {file_path} to {folder}...")
                # Extract file
                extract_file(file_path)
                # Move to destination folder
                extracted_folder = os.path.join(folder, "data", assembly)
                shutil.move(extracted_folder, folder_assembly)
            else:
                log.error(f"Download snpEff '{file_url}' failed")
                raise ValueError(f"Download snpEff '{file_url}' failed")

        else:

            log.info(f"Download snpEff databases {[assembly]} - already exists")


def databases_download_genomes(assemblies: list, genomes_folder: str = DEFAULT_GENOME_FOLDER, provider:str = "UCSC", contig_regex:str = None, threads:int = 1) -> None:
    """
    This function downloads genome assemblies using genomepy package with options to specify genome
    folder, provider, contig regex, and number of threads.
    
    :param assemblies: a list of genome assembly names to download
    :type assemblies: list
    :param genomes_folder: The folder where the downloaded genome files will be saved. If no folder is
    specified, the default folder will be used
    :type genomes_folder: str
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

    if not genomes_folder:
        genomes_folder = DEFAULT_GENOME_FOLDER

    # Full path
    genomes_folder = full_path(genomes_folder)

    if os.path.exists(genomes_folder):
        installed_genomes = genomepy.list_installed_genomes(genomes_dir=genomes_folder)
    else:
        installed_genomes = []

    for assembly in assemblies:
        if assembly in installed_genomes:
            log.info(f"Download Genomes {[assembly]} - already exists")
        else:
            log.info(f"Download Genomes {[assembly]} downloading...")
            genomepy.install_genome(assembly, annotation=False, provider=provider, genomes_dir=genomes_folder, threads=threads, regex=contig_regex)

    return None


def databases_download_refseq(assemblies:list, refseq_folder:str = None, refseq_url:str = None, refseq_prefix:str = "ncbiRefSeq", refseq_files:List = ["ncbiRefSeq.txt", "ncbiRefSeqLink.txt"], refseq_format_file:str = "ncbiRefSeq.txt", refseq_format_file_output:str = None, include_utr_5:bool = True, include_utr_3:bool = True, include_chrM:bool = True, include_non_canonical_chr:bool = True, include_non_coding_transcripts:bool = True, include_transcript_ver:bool = True, threads:int = 1) -> dict:
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
    :param include_utr_3: The `include_utr_3` parameter is a boolean value that specifies whether to
    include the 3' untranslated region (UTR) in the downloaded RefSeq files. If set to `True`, the 3'
    UTR will be included. If set to `False`, the 3, defaults to True
    :type include_utr_3: bool (optional)
    :param include_chrM: The `include_chrM` parameter is a boolean value that determines whether to
    include the mitochondrial chromosome (chrM) in the downloaded RefSeq files. If set to True, the chrM
    will be included; if set to False, it will be excluded, defaults to True
    :type include_chrM: bool (optional)
    :param include_non_canonical_chr: The `include_non_canonical_chr` parameter is a boolean value that
    determines whether or not to include non-canonical chromosomes in the downloaded RefSeq files. If
    set to `True`, non-canonical chromosomes will be included. If set to `False`, non-canonical
    chromosomes will be excluded, defaults to True
    :type include_non_canonical_chr: bool (optional)
    :param include_non_coding_transcripts: The `include_non_coding_transcripts` parameter is a boolean
    flag that determines whether non-coding transcripts should be included in the downloaded RefSeq
    files. If set to `True`, non-coding transcripts will be included. If set to `False`, non-coding
    transcripts will be excluded, defaults to True
    :type include_non_coding_transcripts: bool (optional)
    :param include_transcript_ver: The `include_transcript_ver` parameter is a boolean value that
    determines whether to include the transcript version in the downloaded RefSeq files. If set to
    `True`, the transcript version will be included. If set to `False`, the transcript version will be
    excluded, defaults to True
    :type include_transcript_ver: bool (optional)
    :param threads: The `threads` parameter specifies the number of threads to use for downloading and
    extracting the RefSeq files. It determines the level of parallelism in the download and extraction
    process. By default, it is set to 1, which means that the download and extraction will be performed
    sequentially. If you want, defaults to 1
    :type threads: int (optional)
    :return: The function `databases_download_refseq` returns a dictionary `installed_refseq` which
    contains information about the downloaded RefSeq files for each assembly. The keys of the dictionary
    are the assembly names, and the values are lists of the installed RefSeq files for each assembly.
    """

    # Log
    log.info(f"Download refSeq databases {assemblies}")

    # Default refSeq Folder
    if not refseq_folder:
        refseq_folder = DEFAULT_REFSEQ_FOLDER

    # Full path
    refseq_folder = full_path(refseq_folder)

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
                        download_file(file_url, file_path, threads=threads)
                        # Extract file
                        extract_file(file_path, threads=threads)
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
                    file_path_bed = f"{refseq_format_file_output}.uncompress.txt"
                    file_path_bed_compressed = refseq_format_file_output
                else:
                    file_path_bed = re.sub(r'txt$', 'bed', file_path)
                    file_path_bed_compressed = f"{file_path_bed}.gz"
                file_path_bed_hdr = f"{file_path_bed}.hdr"
                file_path_bed_compressed_hdr = f"{file_path_bed_compressed}.hdr"
                file_path_bed_basename = os.path.basename(file_path_bed)
                file_path_bed_compressed_basename = os.path.basename(file_path_bed_compressed)
                if refseq_file == refseq_format_file and (new_refseq_file or not os.path.exists(file_path_bed)):
                    log.info(f"Download refSeq databases ['{assembly}'] - '{refseq_file}' formatting to '{file_path_bed_basename}'...")
                    # Format refSeq in BED format
                    databases_format_refseq(refseq_file=file_path, output_file=file_path_bed, include_utr_5=include_utr_5, include_utr_3=include_utr_3, include_chrM=include_chrM, include_non_canonical_chr=include_non_canonical_chr, include_non_coding_transcripts=include_non_coding_transcripts, include_transcript_ver=include_transcript_ver, sort=True, header=True, header_first_line=True)
                    # compress file
                    log.info(f"Download refSeq databases ['{assembly}'] - '{file_path_bed_basename}' compressing to '{file_path_bed_compressed_basename}'...")
                    concat_and_compress_files(input_files=[file_path_bed], output_file=file_path_bed_compressed, threads=threads, compression_type="bgzip", sort=False, index=False)
                    shutil.copy(file_path_bed_hdr, file_path_bed_compressed_hdr)
                    

    # Log
    log.debug(f"installed_refseq: {installed_refseq}")

    return installed_refseq


def databases_format_refseq(refseq_file:str, output_file:str, include_utr_5:bool = True, include_utr_3:bool = True, include_chrM:bool = True, include_non_canonical_chr:bool = True, include_non_coding_transcripts:bool = True, include_transcript_ver:bool = True, sort:bool = False, header:bool = False, header_first_line:bool = True) -> str:
    """
    The `databases_format_refseq` function takes a RefSeq file as input, formats it according to
    specified criteria, and outputs the formatted file.
    
    :param refseq_file: The `refseq_file` parameter is a string that represents the path to the input
    RefSeq file. This file contains information about gene annotations, including chromosome, start and
    end positions, strand, and other details
    :type refseq_file: str
    :param output_file: The `output_file` parameter is a string that represents the name of the file
    where the formatted data will be written
    :type output_file: str
    :param include_utr_5: The `include_utr_5` parameter is a boolean that determines whether to include
    the 5' UTR (untranslated region) in the output file. If set to `True`, the 5' UTR will be included.
    If set to `False`, the 5' U, defaults to True
    :type include_utr_5: bool (optional)
    :param include_utr_3: A boolean parameter that determines whether to include the 3' UTR
    (untranslated region) in the output. If set to True, the 3' UTR will be included. If set to False,
    the 3' UTR will be excluded, defaults to True
    :type include_utr_3: bool (optional)
    :param include_chrM: The `include_chrM` parameter is a boolean that determines whether to include
    transcripts from the mitochondrial chromosome (chrM or chrMT) in the output file. If set to True,
    transcripts from the mitochondrial chromosome will be included. If set to False, transcripts from
    the mitochondrial chromosome will be excluded, defaults to True
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
    :param include_transcript_ver: The `include_transcript_ver` parameter determines whether to include
    the transcript version in the output file. If set to `True`, the transcript version will be included
    in the output file. If set to `False`, the transcript version will be removed from the output file.
    The default value is `True, defaults to True
    :type include_transcript_ver: bool (optional)
    :param sort: The `sort` parameter determines whether to sort the output file in ascending order
    based on the chromosome and start position. If set to `True`, the file will be sorted. If set to
    `False`, the file will not be sorted. The default value is `False`, defaults to False
    :type sort: bool (optional)
    :param header: The `header` parameter is a boolean that determines whether to include a header line
    in the output file. If set to `True`, a header line will be included. If set to `False`, no header
    line will be included. The default value is `False`, defaults to False
    :type header: bool (optional)
    :param header_first_line: The `header_first_line` parameter is a boolean that determines whether to
    include the header line as the first line in the output file. If set to `True`, the header line will
    be included as the first line. If set to `False`, the header line will not be included as the first,
    defaults to True
    :type header_first_line: bool (optional)
    :return: The function `databases_format_refseq` returns the path of the output file.
    """

    # Header
    header_line = '#CHROM\tSTART\tEND\tname\ttranscript\tstrand\texon'

    # Open refSeq file
    with open(refseq_file, "r") as fi:

        # Open refSeq output file
        with open(output_file, "w") as fo:

            # Write haeder line
            if header_first_line:
                fo.write(header_line+"\n")

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
    
    # Sort BED file
    if sort:
        bed_sort(input=output_file, output=output_file)

    # Create header
    if header:

        # Header file
        header_file = f"{output_file}.hdr"

        # VCF header dict
        header_dict = {
            "name": {
                "number": ".",
                "type": "String",
                "description": "Name of the refSeq Gene",
                "source": "refSeq",
                "version": "refSeq"
            },
            "transcript": {
                "number": ".",
                "type": "String",
                "description": "Name of the refSeq Transcript",
                "source": "refSeq",
                "version": "refSeq"
            },
            "strand": {
                "number": ".",
                "type": "String",
                "description": "DNA strand orientation",
                "source": "refSeq",
                "version": "refSeq"
            },
            "exon": {
                "number": ".",
                "type": "String",
                "description": "Name of the refSeq Exon",
                "source": "refSeq",
                "version": "refSeq"
            },
        }

        # Create header
        header_list = [
                '##fileformat=VCFv4.2',
                header_line
                ]
        header_vcf = vcf.Reader(io.StringIO("\n".join(header_list)))

        # List of header annotations
        for annotation_name in header_dict:
            header_vcf.infos[annotation_name] = vcf.parser._Info(
                        annotation_name,
                        header_dict[annotation_name]["number"],
                        header_dict[annotation_name]["type"],
                        header_dict[annotation_name]["description"],
                        header_dict[annotation_name]["source"],
                        header_dict[annotation_name]["version"],
                        code_type_map[header_dict[annotation_name]["type"]]
                    )

        # Write file
        f = open(header_file, 'w')
        vcf.Writer(f, header_vcf)
        f.close()

    return output_file


def databases_download_dbnsfp(assemblies:list, dbnsfp_folder:str = None, dbnsfp_url:str = None, dbnsfp_release:str = "4.4a", threads:int = None, parquet_size:int = 100, generate_parquet_file:bool = False, generate_sub_databases:bool = False, generate_vcf_file:bool = False, not_generate_files_all:bool = False, genomes_folder:str = None, add_info:bool = False, row_group_size:int = 100000) -> bool:
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
    :param not_generate_files_all: A boolean flag indicating to not generate database Parquet/VCF file for the entire database,
    defaults to False
    :type not_generate_files_all: bool (optional)
    :param genomes_folder: A string that specifies where are genomes
    :type genomes_folder: str (optional)
    :param add_info: Add INFO column in Parquet folder/file
    :type add_info: bool (optional)
    :param row_group_size: Row group size to generate parquet folder and file (see duckDB doc)
    :type row_group_size: int (optional)
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
        :param threads: The `threads` parameter specifies the number of threads to use for parallel
        processing. It determines how many tasks can be executed simultaneously. Increasing the number of
        threads can potentially speed up the execution time of the function, especially if there are
        multiple cores available on the machine
        :type threads: int (optional)
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


    def get_columns_select_clause(columns_structure:dict, assembly:str = "hg19", for_parquet:bool = False, sub_database:str = None, null_if_no_annotation:bool = True, print_log:bool = True, add_info:bool = False) -> dict:
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
        :param print_log: The `log` parameter is a boolean flag that determines whether to enable logging or
        not. If set to `True`, the function will log debug messages during execution. If set to `False`,
        logging will be disabled, defaults to True
        :type print_log: bool (optional)
        :param add_info: Add INFO column in Parquet folder/file
        :type add_info: bool (optional)
        :return: The `get_columns_select_clause` function returns a dictionary containing the "select"
        and "where" clauses for a SQL query. The "select" clause includes the columns to be selected in
        the query, while the "where" clause includes the conditions for filtering the data.
        Additionally, the dictionary also includes a "annotations" key that contains a dictionary of the
        selected annotations and their corresponding column types.
        """
        

        columns_select_position = {}
        columns_select_ref_alt = {}
        columns_select_info = {}
        columns_select_annotations = {}
        columns_where_annotations = {}
        columns_list_annotations = {}
        columns_info_annotations = {}
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
                    column_alias = clean_name(column)
                    if sub_database is None or sub_database in ['ALL'] or column.upper() == sub_database or column_alias.upper().startswith(f"{sub_database}_"):
                        if for_parquet:
                            column_key = f""" "{column_alias}" """
                            column_where = f""" "{column_alias}" IS NOT NULL """
                            column_info_key = f"""
                                CASE
                                    WHEN "{column_alias}" IS NOT NULL
                                    THEN 
                                        concat(
                                            '{column_alias}=',
                                            "{column_alias}",
                                            ';'
                                        )
                                    ELSE ''
                                END
                            """
                        else:
                            # columns for each annotation
                            column_key = f"""
                                list_aggregate(list_distinct(array_filter(string_split("{column}", ';'), x -> x != '.')), 'string_agg', ',') AS "{column_alias}"
                            """
                            # columns for INFO clumn
                            column_info_key = f"""
                                CASE
                                    WHEN len(list_distinct(array_filter(string_split("{column}", ';'), x -> x != '.'))) > 0
                                    THEN 
                                        concat(
                                            '{column_alias}=',
                                            list_aggregate(list_distinct(array_filter(string_split("{column}", ';'), x -> x != '.')), 'string_agg', ','),
                                            ';'
                                        )
                                    ELSE ''
                                END
                            """
                            column_where = f""" "{column}" IS NOT NULL """
                        columns_info_annotations[column_info_key] = None
                        columns_select_annotations[column_key] = columns_structure[column]
                        columns_where_annotations[column_where] = columns_structure[column]
                        columns_list_annotations[column_alias] = columns_structure[column]

        # Force position and ref/alt if parquet
        if for_parquet:
            columns_select_position = {'"#CHROM"': '#CHROM', '"POS"': 'POS'}
            columns_select_ref_alt = {'"REF"': 'REF', '"ALT"': 'ALT'}

        # DEVEL
        if add_info:
            info_concat = f"""
                regexp_replace(
                    replace(
                        concat({", ".join(list(columns_info_annotations))}),
                        '"', ''),
                    ';$', '')
                AS "INFO"
            """
            columns_select_info[info_concat] = "INFO"
            columns_list_annotations["INFO"] = "VARCHAR"

        # Create select clause
        columns_select = list(columns_select_position.keys()) + list(columns_select_ref_alt.keys()) + list(columns_select_info.keys()) + list(columns_select_annotations.keys())
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

    # Full path
    genomes_folder = full_path(genomes_folder)

    # Default refSeq Folder
    if not dbnsfp_folder:
        dbnsfp_folder = DEFAULT_DBNSFP_FOLDER

    # Full path
    dbnsfp_folder = full_path(dbnsfp_folder)

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
        download_file(dbnsfp_zip_url, dbnsfp_zip_dest, threads=threads)
    else:
        log.info(f"Download dbNSFP {assemblies} - Database '{dbnsfp_zip}' already exists")

    # Extract dbNSFP
    if not os.path.exists(dbnsfp_readme_dest):
        log.info(f"Download dbNSFP {assemblies} - Extract '{dbnsfp_zip}'...")
        extract_file(dbnsfp_zip_dest, threads=threads)
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

        # Full path
        header_file = full_path(header_file)

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
            return 0


    # Init
    partition_size = size_to_mb(parquet_size)
    sample_size = 1000000
    database_files = []
    columns_structure = {}
    readme_annotations_description = None
    
    if not threads or int(threads) <= 0:
        threads = os.cpu_count()
    threads = int(threads)

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
        columns_clauses = get_columns_select_clause(columns_structure=columns_structure, assembly=assembly, sub_database=sub_database, add_info=add_info)
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
                        TO '{output_parquet}' WITH (FORMAT PARQUET, PER_THREAD_OUTPUT TRUE, ROW_GROUP_SIZE {row_group_size})
                    """
                
                # Log
                log.debug(f"""Download dbNSFP ['{assembly}'] - Database '{sub_database}' - Process file '{os.path.basename(database_file)}' - Write parquet files...""")

                if os.path.exists(f"{output_parquet}"):
                    log.debug(f"""Download dbNSFP ['{assembly}'] - Database '{sub_database}' - Process file '{os.path.basename(database_file)}' - Parquet file already exists""")

                else:

                    # Query
                    db_copy = duckdb.connect(config={"threads":threads})
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
            if not columns_structure and len(database_files):
                columns_structure = get_columns_structure(database_file=database_files[0], sample_size=10, threads=threads)

            output_prefix = f"{dbnsfp_folder}/{assembly}/dbNSFP{dbnsfp_release}"
            parquet_all_annotation = f"{output_prefix}.ALL.partition.parquet/*/*.parquet"

            sub_databases_structure = {}

            # Find sub databases
            for column in columns_structure:
                sub_database = clean_name(column).split("_")[0].upper()
                if sub_database not in sub_databases_structure:
                    sub_databases_structure[sub_database] = {}
                sub_databases_structure[sub_database][column] = columns_structure[column]

            # Chromosomes
            chromosomes = []

            # for each sub database
            for sub_database in sorted(set(sub_databases_structure), key=str.casefold, reverse=False):
                
                # Clauses
                columns_clauses = get_columns_select_clause(columns_structure=columns_structure, assembly=assembly, for_parquet=True, sub_database=sub_database, null_if_no_annotation=True, print_log=False, add_info=add_info)
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
                        log.debug(f"""Download dbNSFP ['{assembly}'] - Database '{sub_database}' - Parquet folder already exists""")
                        parquet_partition_already_generated_list.append(sub_database)
                    else:
                        log.info(f"""Download dbNSFP ['{assembly}'] - Database '{sub_database}' - Parquet folder...""")
                        
                        db_copy = duckdb.connect(config={"threads":threads})

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
                                TO '{parquet_sub_database_annotation_chromosome}' WITH (FORMAT PARQUET, PER_THREAD_OUTPUT TRUE, ROW_GROUP_SIZE {row_group_size})
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

                # Skip ALL sub database
                if sub_database in ['ALL'] and not_generate_files_all:
                    continue

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
                            list_of_queries = []
                            list_for_vcf = []

                            # Log
                            log.info(f"""Download dbNSFP ['{assembly}'] - Database '{sub_database}' - VCF file...""")

                            # Find database structure
                            if not database_files:
                                database_files = get_database_files(f"{dbnsfp_zip_dest_folder}/dbNSFP{dbnsfp_release}_variant.chr*.gz")
                            if not columns_structure:
                                columns_structure = get_columns_structure(database_file=database_files[0], sample_size=10, threads=threads)
                            columns_clauses = get_columns_select_clause(columns_structure=columns_structure, assembly=assembly, for_parquet=True, sub_database=sub_database, null_if_no_annotation=True, print_log=False, add_info=add_info)
                            columns_annotations = columns_clauses.get("annotations")

                            # Generate columns concat for INFO column
                            column_info = []

                            # INFO column check
                            query_info_check = f"""
                                    SELECT *
                                    FROM read_parquet('{input_parquet_files[0]}')
                                """
                            columns_info_check = list(db.query(query_info_check).df())
                            
                            # If INFO column exists
                            if "INFO" in columns_info_check:
                                column_info.append(""" "INFO" """)

                            else:
                                for column_annotation in columns_annotations:
                                    column_info.append(f"""
                                    CASE
                                        WHEN "{column_annotation}" IS NOT NULL
                                        THEN concat('{column_annotation}=', replace("{column_annotation}", ';', ','), ';')
                                        ELSE ''
                                    END
                                    """)

                            if column_info:
                                
                                # Generate VCF files within tmp folder
                                with TemporaryDirectory(dir=output_assembly) as tmp_dir:
                                    
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
                                        file_num += 1
                                        vcf_file_gz = os.path.join(tmp_dir,f'variants.{file_num}.tsv.gz')
                                        query_copy = f"""
                                                COPY (
                                                    SELECT "#CHROM", POS, '.' AS ID, REF, ALT, 0 AS QUAL, 'PASS' AS FILTER, regexp_replace(replace(concat({", ".join(column_info)}),'"', ''), ';$', '') AS INFO
                                                    FROM read_parquet('{parquet_file}')
                                                    )
                                                TO '{vcf_file_gz}' WITH (FORMAT CSV, DELIM '\t', HEADER 0, QUOTE '', COMPRESSION 'gzip')
                                                """
                                        # Query
                                        list_of_queries.append(query_copy)
                                        #db_copy.execute(query_copy)
                                        list_for_vcf.append(vcf_file_gz)

                                    # Log
                                    log.debug(f"""Download dbNSFP ['{assembly}'] - Database '{sub_database}' - VCF file - {len(list_for_vcf)} files process in parallel [{threads} threads]...""")
                                    with Pool(processes=threads) as pool:
                                        pool.map(duckdb_execute, list_of_queries)

                                    # Log
                                    log.debug(f"""Download dbNSFP ['{assembly}'] - Database '{sub_database}' - VCF file - {len(list_for_vcf)} files generated""")

                                    # Header

                                    # Init
                                    header_list = []

                                    # Temporary header file for INFO fields
                                    header_file_tmp = os.path.join(tmp_dir,"header.hdr")
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
                                    concat_and_compress_files(input_files=list_for_vcf, output_file=output_vcf, sort=True, index=True)

                                    db_copy.close()

                                    vcf_generated_list.append(sub_database)

                            else:
                                log.error(f"""Download dbNSFP ['{assembly}'] - Database '{sub_database}' - VCF file - No INFO columns found""")
                                raise ValueError(f"""Download dbNSFP ['{assembly}'] - Database '{sub_database}' - VCF file - No INFO columns found""")

                        else:
                            log.error(f"""Download dbNSFP ['{assembly}'] - Database '{sub_database}' - VCF file - No parquet files found""")
                            raise ValueError(f"""Download dbNSFP ['{assembly}'] - Database '{sub_database}' - VCF file - No parquet files found""")
                                
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


def databases_download_alphamissense(assemblies:list, alphamissense_folder:str = DEFAULT_ANNOTATIONS_FOLDER, alphamissense_url:str = DEFAULT_ALPHAMISSENSE_URL, threads:int = None) -> bool:
    """
    The `databases_download_alphamissense` function downloads and converts AlphaMissense databases for a
    list of assemblies.
    
    :param assemblies: `assemblies` is a list of assemblies for which the AlphaMissense database needs
    to be downloaded. Each assembly represents a specific genome or genetic sequence
    :type assemblies: list
    :param alphamissense_folder: The `alphamissense_folder` parameter is a string that specifies the
    folder where the AlphaMissense files will be downloaded and stored. It is set to
    `DEFAULT_ANNOTATIONS_FOLDER` by default, which is likely a predefined constant or variable in your
    code
    :type alphamissense_folder: str
    :param alphamissense_url: The `alphamissense_url` parameter is a string that specifies the URL where
    the AlphaMissense files are located. It is used to construct the download URL for each assembly's
    AlphaMissense file
    :type alphamissense_url: str
    :param threads: The `threads` parameter is an optional parameter that specifies the number of
    threads to use for the conversion process. It determines the level of parallelism when converting
    the AlphaMissense TSV file to the Parquet format. If not specified, the default value will be used
    :type threads: int
    :return: The function `databases_download_alphamissense` returns a boolean value `True`.
    """

    from howard.objects.database import Database

    # Log
    log.info(f"Download AlphaMissense {assemblies}")

    # Full path
    alphamissense_folder = full_path(alphamissense_folder)

    for assembly in assemblies:

        # Files for AlphaMissense
        alphamissense_tsv = f"AlphaMissense_{assembly}.tsv.gz"
        alphamissense_tsv_url = os.path.join(alphamissense_url, alphamissense_tsv)
        alphamissense_tsv_dest = os.path.join(alphamissense_folder, assembly, alphamissense_tsv)
        alphamissense_tsv_dest_folder = os.path.dirname(alphamissense_tsv_dest)
        alphamissense_parquet = f"AlphaMissense.parquet"
        alphamissense_parquet_dest = os.path.join(alphamissense_folder, assembly, alphamissense_parquet)

        # Create folder if not exists
        if not os.path.exists(alphamissense_tsv_dest_folder):
            Path(alphamissense_tsv_dest_folder).mkdir(parents=True, exist_ok=True)

        log.debug(f"{alphamissense_tsv}, {alphamissense_tsv_url}, {alphamissense_tsv_dest}, {alphamissense_tsv_dest_folder}")

        # Download AlphaMissense
        if not os.path.exists(alphamissense_tsv_dest):
            log.info(f"Download AlphaMissense {assemblies} - Download '{alphamissense_tsv}'...")
            download_file(alphamissense_tsv_url, alphamissense_tsv_dest, threads=threads)
        else:
            log.info(f"Download AlphaMissense {assemblies} - Database '{alphamissense_tsv}' already exists")

        # Convert to Parquet
        if os.path.exists(alphamissense_tsv_dest):
            if not os.path.exists(alphamissense_parquet_dest):
                log.info(f"Download AlphaMissense {assemblies} - Convert to '{alphamissense_parquet}'...")
                alphamissense_database = Database(database=alphamissense_tsv_dest, conn_config={"threads":threads})
                alphamissense_database.export(output_database=alphamissense_parquet_dest, output_header=alphamissense_parquet_dest+".hdr", threads=threads)
            else:
                log.info(f"Download AlphaMissense {assemblies} - Database '{alphamissense_parquet}' already exists")
        else:
            log.error(f"Download AlphaMissense {assemblies} - Database '{alphamissense_tsv}' DOES NOT exists")
    
    return True



def databases_download_exomiser(assemblies:list, exomiser_folder:str = DEFAULT_EXOMISER_FOLDER, exomiser_application_properties:str = None, exomiser_url:str = DEFAULT_EXOMISER_URL, exomiser_release:str = None, exomiser_phenotype_release:str = None, exomiser_remm_release:str = None, exomiser_remm_url:str = "https://kircherlab.bihealth.org/download/ReMM", exomiser_cadd_release:str = None, exomiser_cadd_url:str = "https://kircherlab.bihealth.org/download/CADD", exomiser_cadd_url_snv_file:str = "whole_genome_SNVs.tsv.gz", exomiser_cadd_url_indel_file:str = "InDels.tsv.gz", threads:int = 1) -> bool:
    """
    The `databases_download_exomiser` function downloads and sets up the Exomiser database for the
    specified assemblies.
    
    :param assemblies: A list of assemblies for which to download Exomiser databases. Each assembly is a
    string representing a genome build, such as "GRCh37" or "GRCh38"
    :type assemblies: list
    :param exomiser_folder: The `exomiser_folder` parameter is a string that specifies the folder where
    the Exomiser databases will be downloaded and stored. If the folder does not exist, it will be
    created
    :type exomiser_folder: str
    :param exomiser_application_properties: The `exomiser_application_properties` parameter is a string
    representing the path to the Exomiser application properties file. This file contains configuration
    settings for the Exomiser tool. If this parameter is not provided, the function will attempt to
    locate the application properties file automatically based on the Exomiser
    :type exomiser_application_properties: str
    :param exomiser_url: The `exomiser_url` parameter is the URL where the Exomiser database files can
    be downloaded from. It is used to construct the download URLs for the phenotype and assembly files
    :type exomiser_url: str
    :param exomiser_release: The `exomiser_release` parameter is used to specify the version of the
    Exomiser data to download. If it is set to "default", "auto", or "config", the function will attempt
    to retrieve the version from the `exomiser.application.properties` file. If it is
    :type exomiser_release: str
    :param exomiser_phenotype_release: The `exomiser_phenotype_release` parameter is used to specify the
    release version of the Exomiser phenotype database. If not provided, it will default to the value
    specified in the `application.properties` file or the latest available release
    :type exomiser_phenotype_release: str
    :param exomiser_remm_release: The `exomiser_remm_release` parameter is used to specify the version
    of the ReMM (Regulatory Mendelian Mutation) database to download. If the value is set to "default",
    "auto", or "config", it will try to retrieve the version from the `application.properties`
    :type exomiser_remm_release: str
    :param exomiser_remm_url: The `exomiser_remm_url` parameter is the URL where the ReMM (Regulatory
    Mendelian Mutation) database can be downloaded from. It is used in the function to construct the
    download URL for the ReMM database files, defaults to https://kircherlab.bihealth.org/download/ReMM
    :type exomiser_remm_url: str (optional)
    :param exomiser_cadd_release: The `exomiser_cadd_release` parameter is used to specify the version
    of the CADD (Combined Annotation Dependent Depletion) database to download. If the value is set to
    "default", "auto", or "config", it will try to retrieve the version from the `exom
    :type exomiser_cadd_release: str
    :param exomiser_cadd_url: The `exomiser_cadd_url` parameter is the URL where the CADD (Combined
    Annotation Dependent Depletion) database files can be downloaded from. It is used to construct the
    download URLs for the CADD database files, defaults to https://kircherlab.bihealth.org/download/CADD
    :type exomiser_cadd_url: str (optional)
    :param exomiser_cadd_url_snv_file: The parameter `exomiser_cadd_url_snv_file` is the name of the
    file containing the SNV (Single Nucleotide Variant) data for the CADD (Combined Annotation Dependent
    Depletion) database, defaults to whole_genome_SNVs.tsv.gz
    :type exomiser_cadd_url_snv_file: str (optional)
    :param exomiser_cadd_url_indel_file: The parameter `exomiser_cadd_url_indel_file` is the name of the
    INDEL file that will be downloaded from the CADD database, defaults to InDels.tsv.gz
    :type exomiser_cadd_url_indel_file: str (optional)
    :param threads: The `threads` parameter specifies the number of threads to use for parallel
    processing. It determines how many tasks can be executed simultaneously. Increasing the number of
    threads can potentially speed up the execution time of the function, especially if there are
    multiple cores available on the machine
    :type threads: int (optional)
    """

    log.info(f"Download Exomiser {assemblies}")

    # Variables
    transcript_source_default = "refseq"
    exomiser_release_default = "2109"

    # Full Path
    exomiser_folder = full_path(exomiser_folder)
    exomiser_application_properties = full_path(exomiser_application_properties)

    # Create folder if not exists
    if not os.path.exists(exomiser_folder):
        Path(exomiser_folder).mkdir(parents=True, exist_ok=True)

    # application.properties as dict
    exomiser_application_properties_dict = {}

    # Find exomiser_application_properties
    if not exomiser_application_properties:
        exomiser_jar = get_bin(bin="exomiser-cli*.jar", tool="exomiser", bin_type="jar", default_folder=f"{DEFAULT_TOOLS_FOLDER}/exomiser")
        if exomiser_jar and os.path.exists(exomiser_jar):
            exomiser_jar_dirname = os.path.dirname(exomiser_jar)
            exomiser_application_properties = os.path.join(exomiser_jar_dirname, "application.properties")

    if exomiser_application_properties and os.path.exists(exomiser_application_properties):
        configs = Properties()
        with open(exomiser_application_properties, 'rb') as read_prop:
            configs.load(read_prop)
        for item in configs.items():
            exomiser_application_properties_dict[item[0]] = item[1][0]

    log.debug(exomiser_application_properties_dict)

    for assembly in assemblies:

        log.info(f"Download Exomiser ['{assembly}']")

        if not exomiser_release or exomiser_release.lower() in ["default", "auto", "config"]:
            exomiser_release_found = exomiser_application_properties_dict.get(f"exomiser.{assembly}.data-version", exomiser_release_default)
        else:
            exomiser_release_found = exomiser_release

        if not exomiser_phenotype_release or exomiser_phenotype_release.lower() in ["default", "auto", "config"]:
            exomiser_phenotype_release_found = exomiser_application_properties_dict.get(f"exomiser.phenotype.data-version", exomiser_release_found)
        else:
            exomiser_phenotype_release_found = exomiser_phenotype_release


        log.debug(f"exomiser_release: {exomiser_release_found}")
        log.debug(f"exomiser_phenotype_release: {exomiser_phenotype_release_found}")

        # Phenotype
        exomiser_phenotype_filename = f"{exomiser_phenotype_release_found}_phenotype.zip"
        exomiser_phenotype_filename_base = f"{exomiser_phenotype_release_found}_phenotype"
        exomiser_phenotype_file = os.path.join(exomiser_folder, exomiser_phenotype_filename)
        exomiser_phenotype_file_base = os.path.join(exomiser_folder, exomiser_phenotype_filename_base)
        exomiser_download_phenotype_url = os.path.join(exomiser_url, "data", exomiser_phenotype_filename)

        # Download Zip file 
        if not os.path.exists(exomiser_phenotype_file):
            log.info(f"Download Exomiser {assemblies} - Download Phenotype '{exomiser_phenotype_filename}'...")
            download_file(url=exomiser_download_phenotype_url, dest_file_path=exomiser_phenotype_file, threads=threads)
        else:
            log.info(f"Download Exomiser {assemblies} - Database Phenotype '{exomiser_phenotype_filename}' already exists")

        # Extract Zip file
        if not os.path.exists(exomiser_phenotype_file_base):
            log.info(f"Download Exomiser {assemblies} - Extract Phenotype '{exomiser_phenotype_filename}'...")
            extract_file(file_path=exomiser_phenotype_file, path=None, threads=threads)
        else:
            log.info(f"Download Exomiser {assemblies} - Database Phenotype '{exomiser_phenotype_filename}' already extracted")

        # Assembly
        exomiser_assembly_filename = f"{exomiser_release_found}_{assembly}.zip"
        exomiser_assembly_filename_base = f"{exomiser_release_found}_{assembly}"
        exomiser_assembly_file = os.path.join(exomiser_folder, exomiser_assembly_filename)
        #exomiser_assembly_file_base = os.path.join(exomiser_folder, exomiser_assembly_filename_base)
        exomiser_assembly_folder = os.path.join(exomiser_folder, assembly)
        exomiser_download_assembly_url = os.path.join(exomiser_url, "data", exomiser_assembly_filename)

        # Download Zip file 
        if not os.path.exists(exomiser_assembly_file):
            log.info(f"Download Exomiser {assemblies} - Download Data '{exomiser_assembly_filename}'...")
            download_file(url=exomiser_download_assembly_url, dest_file_path=exomiser_assembly_file, threads=threads)
        else:
            log.info(f"Download Exomiser {assemblies} - Database Data '{exomiser_assembly_filename}' already exists")

        # Extract Zip file
        #if not os.path.exists(exomiser_assembly_folder):
        if not os.path.exists(os.path.join(exomiser_assembly_folder,exomiser_assembly_filename_base)):
            log.info(f"Download Exomiser {assemblies} - Extract Data '{exomiser_assembly_filename}'...")
            extract_file(file_path=exomiser_assembly_file, path=exomiser_assembly_folder, threads=threads)
        else:
            log.info(f"Download Exomiser {assemblies} - Database Data '{exomiser_assembly_filename}' already extracted")

        # Link Phenotype to assembly
        phenotype_assembly_link_src = f"../{exomiser_phenotype_filename_base}"
        phenotype_assembly_link_dest = os.path.join(exomiser_assembly_folder, exomiser_phenotype_filename_base)
        if not os.path.exists(phenotype_assembly_link_dest):
            log.info(f"Download Exomiser {assemblies} - Create Data '{os.path.basename(phenotype_assembly_link_dest)}' Phenotype link...")
            
            os.symlink(phenotype_assembly_link_src, phenotype_assembly_link_dest)
        else:
            log.info(f"Download Exomiser {assemblies} - Database Data '{os.path.basename(phenotype_assembly_link_dest)}' Phenotype link already created")

        # Generate properties for assembly
        exomiser_application_properties_assembly_file = os.path.join(exomiser_assembly_folder, "application.properties")
        exomiser_application_properties_assembly = Properties()
        exomiser_application_properties_assembly["exomiser.data-directory"] = exomiser_assembly_folder
        exomiser_application_properties_assembly[f"exomiser.{assembly}.data-version"] = exomiser_release_found
        exomiser_application_properties_assembly[f"exomiser.{assembly}.variant-white-list-path"] = f"{exomiser_release_found}_{assembly}_clinvar_whitelist.tsv.gz"
        exomiser_application_properties_assembly[f"exomiser.phenotype.data-version"] = exomiser_phenotype_release_found
        exomiser_application_properties_assembly[f"exomiser.{assembly}.transcript-source"] = transcript_source_default


        # REMM
        exomiser_remm_download = False

        if exomiser_remm_release:
            # Search release if autodetection
            if exomiser_remm_release.lower() in ["default", "auto", "config"]:
                exomiser_remm_release_found = exomiser_application_properties_dict.get("remm.version", None)
                log.debug(f"exomiser_remm_release found: {exomiser_remm_release}")
            else:
                exomiser_remm_release_found = exomiser_remm_release
            # Force create of path and download
            exomiser_remm_download = True
        else:
            # Check REMM release
            exomiser_remm_release_found = exomiser_application_properties_dict.get("remm.version", None)

        # Check if path is mandatory in application.properties
        if exomiser_remm_release_found and exomiser_application_properties_dict.get(f"exomiser.{assembly}.remm-path", None):
            exomiser_remm_download = True

        exomiser_assembly_remm_folder = os.path.join(exomiser_assembly_folder, "remm")
        exomiser_remm_path = os.path.join(exomiser_assembly_remm_folder, f"ReMM.v{exomiser_remm_release_found}.{assembly}.tsv.gz")
        
        if exomiser_remm_download and not os.path.exists(exomiser_remm_path):
            if not os.path.exists(exomiser_assembly_remm_folder):
                Path(exomiser_assembly_remm_folder).mkdir(parents=True, exist_ok=True)
            exomiser_remm_path_tbi = f"{exomiser_remm_path}.tbi"
            exomiser_remm_path_md5 = exomiser_remm_path.replace(".tsv.gz", ".md5")
            log.info(f"Download Exomiser {assemblies} - Download REMM database '{os.path.basename(exomiser_remm_path)}'...")
            exomiser_download_assembly_remm_url = os.path.join(exomiser_remm_url, f"ReMM.v{exomiser_remm_release_found}.{assembly}.tsv.gz") 
            exomiser_download_assembly_remm_tbi_url = os.path.join(exomiser_remm_url, f"ReMM.v{exomiser_remm_release_found}.{assembly}.tsv.gz.tbi")
            exomiser_download_assembly_remm_md5_url = os.path.join(exomiser_remm_url, f"ReMM.v{exomiser_remm_release_found}.{assembly}.md5")
            download_file(url=exomiser_download_assembly_remm_url, dest_file_path=exomiser_remm_path, threads=threads)
            download_file(url=exomiser_download_assembly_remm_tbi_url, dest_file_path=exomiser_remm_path_tbi, threads=threads)
            download_file(url=exomiser_download_assembly_remm_md5_url, dest_file_path=exomiser_remm_path_md5, threads=threads)
        else:
            log.debug(f"Download Exomiser {assemblies} - Database REMM not downloaded")

        log.debug(f"{exomiser_remm_release_found} and {exomiser_remm_path} and os.path.exists({exomiser_remm_path})")
        if exomiser_remm_release_found and exomiser_remm_path and os.path.exists(exomiser_remm_path):
            log.info(f"Download Exomiser {assemblies} - Database REMM '{os.path.basename(exomiser_remm_path)}' already exists")
            exomiser_application_properties_assembly["remm.version"] = exomiser_remm_release_found
            exomiser_application_properties_assembly[f"exomiser.{assembly}.remm-path"] = "/".join(["${exomiser.data-directory}", "remm", f"ReMM.v{exomiser_remm_release_found}.{assembly}.tsv.gz"])


        # CADD
        exomiser_cadd_download = False
        
        if exomiser_cadd_release:
            if exomiser_cadd_release.lower() in ["default", "auto", "config"]:
                exomiser_cadd_release_found = exomiser_application_properties_dict.get("cadd.version", None)
                log.debug(f"exomiser_cadd_release found: {exomiser_cadd_release}")
            else:
                exomiser_cadd_release_found = exomiser_cadd_release
            # Force create of path and download
            exomiser_cadd_download = True
        else:
            # Check CADD release
            exomiser_cadd_release_found = exomiser_application_properties_dict.get("cadd.version", None)

        log.debug(f"exomiser_cadd_release: {exomiser_cadd_release_found}")

        # Check if path is mandatory in application.properties
        if exomiser_cadd_release_found and exomiser_application_properties_dict.get(f"exomiser.{assembly}.cadd-snv-path", None) and exomiser_application_properties_dict.get(f"exomiser.{assembly}.cadd-in-del-path", None) :
            exomiser_cadd_download = True

        if exomiser_cadd_release_found:
            exomiser_assembly_cadd_folder = os.path.join(exomiser_assembly_folder, "cadd", exomiser_cadd_release_found)
            exomiser_cadd_snv_path = os.path.join(exomiser_assembly_cadd_folder, exomiser_cadd_url_snv_file)
            exomiser_cadd_indel_path = os.path.join(exomiser_assembly_cadd_folder, exomiser_cadd_url_indel_file)
        else:
            exomiser_cadd_snv_path = None
            exomiser_cadd_indel_path = None
            exomiser_cadd_download = False
        
        if exomiser_cadd_download and (not os.path.exists(exomiser_cadd_snv_path) or not os.path.exists(exomiser_cadd_indel_path)):

            if not os.path.exists(exomiser_assembly_cadd_folder):
                Path(exomiser_assembly_cadd_folder).mkdir(parents=True, exist_ok=True)
            
            # Genome build from gencode
            log.debug(f"Download Exomiser {assemblies} - Search for assembly '{assembly}' into gencode...")
            genome_build_switch_to_gencode = genome_build_switch(assembly=assembly)

            # SNV
            if not os.path.exists(exomiser_cadd_snv_path):
                exomiser_cadd_snv_path_tbi = f"{exomiser_cadd_snv_path}.tbi"
                log.info(f"Download Exomiser {assemblies} - Download CADD database '{os.path.basename(exomiser_cadd_snv_path)}'...")
                exomiser_download_assembly_cadd_snv_url = os.path.join(exomiser_cadd_url, f"v{exomiser_cadd_release_found}", genome_build_switch_to_gencode, exomiser_cadd_url_snv_file) 
                exomiser_download_assembly_cadd_snv_tbi_url = os.path.join(exomiser_cadd_url, f"v{exomiser_cadd_release_found}", genome_build_switch_to_gencode, f"{exomiser_cadd_url_snv_file}.tbi")
                download_file(url=exomiser_download_assembly_cadd_snv_url, dest_file_path=exomiser_cadd_snv_path, threads=threads)
                download_file(url=exomiser_download_assembly_cadd_snv_tbi_url, dest_file_path=exomiser_cadd_snv_path_tbi, threads=threads)

            # INDEL
            if not os.path.exists(exomiser_cadd_indel_path):
                exomiser_cadd_indel_path_tbi = f"{exomiser_cadd_indel_path}.tbi"
                log.info(f"Download Exomiser {assemblies} - Download CADD database '{os.path.basename(exomiser_cadd_indel_path)}'...")
                exomiser_download_assembly_cadd_indel_url = os.path.join(exomiser_cadd_url, f"v{exomiser_cadd_release_found}", genome_build_switch_to_gencode, exomiser_cadd_url_indel_file) 
                exomiser_download_assembly_cadd_indel_tbi_url = os.path.join(exomiser_cadd_url, f"v{exomiser_cadd_release_found}", genome_build_switch_to_gencode, f"{exomiser_cadd_url_indel_file}.gz.tbi")
                download_file(url=exomiser_download_assembly_cadd_indel_url, dest_file_path=exomiser_cadd_indel_path, threads=threads)
                download_file(url=exomiser_download_assembly_cadd_indel_tbi_url, dest_file_path=exomiser_cadd_indel_path_tbi, threads=threads)

        else:
            log.debug(f"Download Exomiser {assemblies} - Database CADD not downloaded")

        log.debug(f"{exomiser_cadd_release_found} and {exomiser_cadd_snv_path} and {exomiser_cadd_indel_path}")
        if exomiser_cadd_release_found and exomiser_cadd_snv_path and os.path.exists(exomiser_cadd_snv_path) and exomiser_cadd_indel_path and os.path.exists(exomiser_cadd_indel_path) :
            log.info(f"Download Exomiser {assemblies} - Database CADD '{os.path.basename(exomiser_cadd_snv_path)}' already exists")
            log.info(f"Download Exomiser {assemblies} - Database CADD '{os.path.basename(exomiser_cadd_indel_path)}' already exists")
            exomiser_application_properties_assembly["cadd.version"] = exomiser_cadd_release_found
            exomiser_application_properties_assembly[f"exomiser.{assembly}.cadd-snv-path"] = "/".join(["${exomiser.data-directory}", "cadd", exomiser_cadd_release_found, exomiser_cadd_url_snv_file])
            exomiser_application_properties_assembly[f"exomiser.{assembly}.cadd-in-del-path"] = "/".join(["${exomiser.data-directory}", "cadd", exomiser_cadd_release_found, exomiser_cadd_url_indel_file])

        
        # Create application.properties
        with open(exomiser_application_properties_assembly_file, "wb") as f:
            exomiser_application_properties_assembly.store(f, encoding="utf-8")
        

    return True

def databases_download_dbsnp(assemblies:list, dbsnp_folder:str = DEFAULT_DBSNP_FOLDER, dbsnp_releases:list = ["b156"], dbsnp_release_default:str = None, dbsnp_url:str = DEFAULT_DBSNP_URL, dbsnp_url_files:dict = None, dbsnp_url_files_prefix:str = "GCF_000001405", dbsnp_assemblies_map:dict = {"hg19": "25", "hg38": "40"}, genomes_folder: str = DEFAULT_GENOME_FOLDER, threads:int = 1, dbsnp_vcf:bool = False, dbsnp_parquet:bool = False, dbsnp_parquet_explode_infos:bool = True) -> str:
    """
    The function `databases_download_dbsnp` downloads dbSNP files, generates VCF files, and converts
    them to Parquet format.
    
    :param assemblies: A list of genome assemblies for which to download dbSNP data
    :type assemblies: list
    :param dbsnp_folder: The folder where the dbSNP files will be downloaded and stored
    :type dbsnp_folder: str
    :param dbsnp_releases: List of releases to download. Default:["b156"]
    :type dbsnp_releases: list
    :param dbsnp_release_default: Default release to link in default folder. Default: first release in dbsnp_releases
    :type dbsnp_release_default: str
    :param dbsnp_url: The `dbsnp_url` parameter is a string that represents the base URL where the dbSNP
    files are located. This URL is used to construct the full URL for downloading the dbSNP files
    :type dbsnp_url: str
    :param dbsnp_url_files: The `dbsnp_url_files` parameter is a dictionary that maps assembly names to
    specific dbSNP URL files. It allows you to provide custom dbSNP URL files for specific assemblies
    instead of using the default file naming convention
    :type dbsnp_url_files: dict
    :param dbsnp_url_files_prefix: The `dbsnp_url_files_prefix` parameter is a string that represents the
    prefix of the dbSNP file name for a specific assembly. It is used to construct the full URL of the
    dbSNP file to be downloaded. By default, the value is set to "GCF_000001405"
    :type dbsnp_url_files_prefix: str (optional)
    :param dbsnp_assemblies_map: The `dbsnp_assemblies_map` parameter is a dictionary that maps assembly
    names to their corresponding dbSNP versions. It is used to construct the dbSNP file name based on
    the assembly name. For example, if the assembly is "hg19", the corresponding dbSNP version is "
    :type dbsnp_assemblies_map: dict
    :param genomes_folder: The `genomes_folder` parameter is a string that specifies the folder where the
    genome index files are located. These index files are used for generating the VCF file from the
    downloaded dbSNP file
    :type genomes_folder: str
    :param threads: The `threads` parameter specifies the number of threads to use for downloading and
    processing the dbSNP files, defaults to 1
    :type threads: int (optional)
    :param dbsnp_vcf: A boolean flag indicating whether to generate a VCF file from the downloaded
    dbSNP data. If set to True, the function will generate a VCF file. If set to False, the function
    will not generate a VCF file, defaults to False
    :type dbsnp_vcf: bool (optional)
    :param dbsnpparquet: A boolean flag indicating whether to generate a Parquet file from the
    downloaded dbSNP data. If set to True, a Parquet file will be generated; if set to False, no Parquet
    file will be generated, defaults to False
    :type dbsnp_parquet: bool (optional)
    """

    # Import Database object
    from howard.objects.database import Database

    # Log
    log.info(f"Download dbSNP {assemblies}")
    
    # Database config       
    conn_config = {"threads": threads}

    # Genomes folder
    if not genomes_folder:
        genomes_folder = DEFAULT_GENOME_FOLDER
    log.debug(f"genomes_folder: {genomes_folder}")

    # Full path
    genomes_folder = full_path(genomes_folder)
    dbsnp_folder = full_path(dbsnp_folder)

    # Default release
    dbsnp_release_default_found = None
    if dbsnp_release_default:
        dbsnp_release_default_found = dbsnp_release_default
    elif len(dbsnp_releases) > 0:
        dbsnp_release_default_found = dbsnp_releases[0]
    log.debug(f"dbsnp_release_default_found: {dbsnp_release_default_found}")

    # dbSNP assemblies map
    if isinstance(dbsnp_assemblies_map, str):
        dbsnp_assemblies_map = json.loads(dbsnp_assemblies_map)
    log.debug(f"dbsnp_assemblies_map: {dbsnp_assemblies_map}")

    # dbSNP URL files
    if isinstance(dbsnp_url_files, str):
        dbsnp_url_files = json.loads(dbsnp_url_files)
    log.debug(f"dbsnp_url_files: {dbsnp_url_files}")

    # Files to generate
    log.debug(f"dbsnp_vcf: {dbsnp_vcf}")
    log.debug(f"dbsnp_parquet: {dbsnp_parquet}")

    for assembly in assemblies:

        # Log
        log.info(f"Download dbSNP {[assembly]}")

        # Folder
        dbsnp_folder_assembly = os.path.join(dbsnp_folder, assembly)

        for dbsnp_release in dbsnp_releases:

            # Log
            log.info(f"Download dbSNP {[assembly]} - Release {[dbsnp_release]}")

            # Folder
            dbsnp_folder_assembly_release = os.path.join(dbsnp_folder_assembly, dbsnp_release)

            # Create folder if not exists
            if not os.path.exists(dbsnp_folder_assembly_release):
                Path(dbsnp_folder_assembly_release).mkdir(parents=True, exist_ok=True)

            # Construct URL
            if dbsnp_url_files and dbsnp_url_files.get(assembly, None):
                dbsnp_url_file = dbsnp_url_files.get(assembly)
            elif dbsnp_assemblies_map.get(assembly, None):
                dbsnp_url_file = f"{dbsnp_url_files_prefix}.{dbsnp_assemblies_map.get(assembly)}.gz"
            else:
                dbsnp_url_file = None

            if dbsnp_url_file:
                dbsnp_url_path = os.path.join(dbsnp_url, dbsnp_release, "VCF", dbsnp_url_file)
            else:
                log.error("No dbSNP file from URL provided")
                raise ValueError("No dbSNP file from URL provided")

            # Construct File
            dbsnp_file = os.path.join(dbsnp_folder_assembly_release, os.path.basename(dbsnp_url_path))

            ### Download dbSNP
            if not os.path.exists(dbsnp_file):

                # Log
                log.info(f"Download dbSNP {[assembly]} - Release {[dbsnp_release]} - Download '{dbsnp_url_file}'...")

                # Create folder if not exists
                if not os.path.exists(os.path.dirname(dbsnp_file)):
                    Path(os.path.dirname(dbsnp_file)).mkdir(parents=True, exist_ok=True)

                # Download file
                download_file(url=dbsnp_url_path, dest_file_path=dbsnp_file, threads=threads)

                # Index tbi
                dbsnp_url_path_tbi = f"{dbsnp_url_path}.tbi"
                dbsnp_file_tbi = f"{dbsnp_file}.tbi"
                download_file(url=dbsnp_url_path_tbi, dest_file_path=dbsnp_file_tbi, threads=threads)

            else:

                # Log
                log.info(f"Download dbSNP {[assembly]} - Release {[dbsnp_release]} - Database {os.path.basename(dbsnp_file)} already downloaded")


            ### Generate VCF and/or Parquet

            # Construct VCF File
            dbsnp_vcf_file = os.path.join(dbsnp_folder_assembly_release, "dbsnp.vcf.gz")
            dbsnp_vcf_file_buffer = os.path.join(dbsnp_folder_assembly_release, "dbsnp.buffer.vcf.gz")

            # Construct File
            dbsnp_parquet_file = os.path.join(dbsnp_folder_assembly_release, "dbsnp.parquet")
            dbsnp_parquet_file_buffer = os.path.join(dbsnp_folder_assembly_release, "dbsnp.buffer.parquet")
            dbsnp_parquet_hdr_file = os.path.join(dbsnp_folder_assembly_release, "dbsnp.parquet.hdr")

            # Files to generate or already generated
            dbsnp_files_to_generate = []
            dbsnp_files_already_generated = []

            if os.path.exists(dbsnp_file) and (dbsnp_vcf and not os.path.exists(dbsnp_vcf_file)) or (dbsnp_parquet and not os.path.exists(dbsnp_parquet_file)):

                if (dbsnp_vcf and not os.path.exists(dbsnp_vcf_file)):
                    write_vcf = True
                    dbsnp_files_to_generate.append(os.path.basename(dbsnp_vcf_file))
                else:
                    write_vcf = False
                    if dbsnp_vcf:
                        dbsnp_files_already_generated.append(os.path.basename(dbsnp_vcf_file))

                if (dbsnp_parquet and not os.path.exists(dbsnp_parquet_file)):
                    write_parquet = True
                    dbsnp_files_to_generate.append(os.path.basename(dbsnp_parquet_file))
                else:
                    write_parquet = False
                    if dbsnp_parquet:
                        dbsnp_files_already_generated.append(os.path.basename(dbsnp_parquet_file))

                # Log
                if dbsnp_files_already_generated:
                    log.info(f"Download dbSNP {[assembly]} - Release {[dbsnp_release]} - Files {dbsnp_files_already_generated} already generated")
                if dbsnp_files_to_generate:
                    log.info(f"Download dbSNP {[assembly]} - Release {[dbsnp_release]} - Files {dbsnp_files_to_generate} generation...")

                # Generate Header

                # Genome index
                genome_index = os.path.join(genomes_folder, assembly, f"{assembly}.fa.fai")
                if not os.path.exists(genome_index):
                    log.error("No genome index")
                    raise ValueError("No genome index")

                # Database
                conn_config_db = conn_config
                db = Database(database=dbsnp_file, format="vcf", conn_config=conn_config_db)
                
                # Header
                header_file_tmp = f"{dbsnp_file}.tmp.hdr"
                db_header = db.get_header(header_file=header_file_tmp)
                db_header_file = db.get_header_file(header_file=header_file_tmp)
                db_header_list = db.read_header_file(header_file=db_header_file)

                # Header info/format...
                db_header_list_new = db_header_list[:-1]
                # Header #CHROM line
                db_header_list_chrom = [db_header_list[-1]]
                # Header Contig
                res = db.query(query=f"""
                            SELECT
                                    column0 AS chr,
                                    concat(
                                        '##contig=<ID=',
                                        column0,
                                        ',length=',
                                        column1,
                                        '>\n'
                                    ) AS contig
                            FROM read_csv_auto('{genome_index}')
                            """)
                db_header_list_chrs = list(res.df()['chr'])
                db_header_list_contigs = list(res.df()['contig'])

                # Write new heaer with contigs
                db_header_new = db.get_header_from_list(header_list=db_header_list_new+db_header_list_contigs+db_header_list_chrom)
                vcf.Writer(open(header_file_tmp, 'w'), db_header_new)

                # Header
                header_file = f"{dbsnp_vcf_file}.hdr"
                shutil.copy(src=header_file_tmp, dst=header_file)

                # INFO fields
                query_select_info_fields_array = []
                if write_parquet and dbsnp_parquet_explode_infos:
                    for info in db_header.infos:
                        if db_header.infos[info].type == "Flag":
                            query_select_info_fields_array.append(f"""
                                concat(';', INFO, ';') LIKE '%;{info};%'
                                AS {info}                             
                            """)
                        else:
                            query_select_info_fields_array.append(f"""
                                CASE
                                    WHEN concat(';', INFO) NOT LIKE '%;{info}=%' THEN NULL
                                    ELSE REGEXP_EXTRACT(concat(';', INFO), ';{info}=([^;]*)',1)
                                END AS {info}                             
                            """)
                query_select_info_fields = " , ".join(query_select_info_fields_array)

                # Chunk CSV
                chunksize = 1000000
                skiprows = len(db_header_list)
                names = ["CHROM_OLD", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"]
                chunk_csv = pd.read_csv(dbsnp_file, low_memory=False, chunksize = chunksize, sep="\t", skiprows=skiprows, names=names)
                chunk_i = 0

                # Open CSV filewith pgzip compression
                with pgzip.open(dbsnp_vcf_file_buffer, mode="w", thread=threads) as f:

                    # fetch chunk
                    for df_chunk in chunk_csv:

                        # Chunk i
                        chunk_i += 1
    
                        # Write CSV
                        if write_vcf:
                            # Query
                            query=f"""
                                SELECT 
                                        concat(
                                            'chr',
                                            replace(replace(replace(regexp_extract("CHROM_OLD", 'NC_[0]*([0-9]*)_?', 1), '23', 'X'), '24', 'Y'), '25', 'M')
                                        ) AS '#CHROM',
                                        "POS", "ID", "REF",
                                        UNNEST(string_split("ALT", ',')) AS 'ALT',
                                        "QUAL", "FILTER", "INFO"
                                FROM df_chunk
                                WHERE list_contains({db_header_list_chrs}, "#CHROM")
                                """
                            res = db.query(query=query)
                            # Use polars to parallelize csv write and an infile with pgzip to parallelise compression
                            res.pl().write_csv(f, separator="\t", include_header=False)

                        # Write Parquet
                        if write_parquet:
                            # Query
                            query=f"""
                                SELECT 
                                        concat(
                                            'chr',
                                            replace(replace(replace(regexp_extract("CHROM_OLD", 'NC_[0]*([0-9]*)_?', 1), '23', 'X'), '24', 'Y'), '25', 'M')
                                        ) AS '#CHROM',
                                        "POS", "ID", "REF",
                                        UNNEST(string_split("ALT", ',')) AS 'ALT',
                                        "QUAL", "FILTER", "INFO",
                                        {query_select_info_fields}
                                FROM df_chunk
                                WHERE list_contains({db_header_list_chrs}, "#CHROM")
                                """
                            res = db.query(query=query)
                            # Use pandas an to_parquet with append option ans fastparquet engine (no append with other df like polars or pyarrow)
                            if chunk_i == 1:
                                append = False
                            else:
                                append = True
                            res.df().to_parquet(dbsnp_parquet_file_buffer, engine='fastparquet', append=append)

                        # Log
                        log.debug(f"chunk {chunk_i} - {chunksize*chunk_i} lines processed")

                # Header

                if write_vcf:
                    # Concat and Compress CSV with header to VCF in BGZIP compression type
                    concat_and_compress_files(input_files=[header_file_tmp, dbsnp_vcf_file_buffer], output_file=dbsnp_vcf_file, threads=threads, compression_type="bgzip", sort=False, index=False)

                if write_parquet:
                    # Move Parquet buffer to Parquet
                    os.rename(dbsnp_parquet_file_buffer, dbsnp_parquet_file)
                    # Copy header file
                    shutil.copy(src=header_file, dst=dbsnp_parquet_hdr_file)

                # Clean
                remove_if_exists([header_file_tmp,dbsnp_vcf_file_buffer])

            else:

                if os.path.exists(dbsnp_vcf_file):
                    dbsnp_files_already_generated.append(os.path.basename(dbsnp_vcf_file))
                if os.path.exists(dbsnp_vcf_file):
                    dbsnp_files_already_generated.append(os.path.basename(dbsnp_parquet_file))
                if dbsnp_files_already_generated:
                    log.info(f"Download dbSNP {[assembly]} - Release {[dbsnp_release]} - Files {dbsnp_files_already_generated} already generated")

        # Generate default release link
        if dbsnp_release_default_found:
            dbsnp_release_default_link_src = dbsnp_release_default_found
            dbsnp_release_default_link_dest = os.path.join(dbsnp_folder_assembly, "default")
            if not os.path.exists(dbsnp_release_default_link_dest) or dbsnp_release_default:
                log.info(f"Download dbSNP {[assembly]} - Release {[dbsnp_release_default_found]} - Defined as default")
                if not os.path.exists(dbsnp_release_default_link_dest) or os.path.islink(dbsnp_release_default_link_dest):
                    remove_if_exists(dbsnp_release_default_link_dest)
                    os.symlink(dbsnp_release_default_link_src, dbsnp_release_default_link_dest)
                else:
                    log.warning(f"Download dbSNP {[assembly]} - Default release path is not a symlink. Change not allowed")
            else:
                dbsnp_release_default_already = os.path.basename(os.path.realpath(dbsnp_release_default_link_dest))
                log.info(f"Download dbSNP {[assembly]} - Release {[dbsnp_release_default_already]} - Already defined as default")

    return True




def databases_download_hgmd(assemblies:list, hgmd_file:str, hgmd_folder:str = DEFAULT_ANNOTATIONS_FOLDER, output_basename:str = None, threads:int = None, genomes_folder:str = None, to_parquet:bool = True, to_tsv:bool = True) -> bool:
    """
    The `databases_download_hgmd` function converts an HGMD database file into VCF, Parquet, and TSV
    formats.
    
    :param assemblies: A list of assemblies for which the HGMD database should be downloaded and
    converted. Only one assembly can be specified
    :type assemblies: list
    :param hgmd_file: The `hgmd_file` parameter is a string that represents the path to the HGMD
    database file in VCF format. This file contains the variants and their associated information
    :type hgmd_file: str
    :param hgmd_folder: The `hgmd_folder` parameter is a string that represents the path to the folder
    where the HGMD database files will be stored. If no value is provided, it will use the
    `DEFAULT_ANNOTATIONS_FOLDER` constant as the default value
    :type hgmd_folder: str
    :param output_basename: The `output_basename` parameter is a string that specifies the base name for
    the output files. If not provided, it will be set as the base name of the input HGMD file without
    the assembly information
    :type output_basename: str
    :param threads: The `threads` parameter specifies the number of threads to use for processing the
    HGMD database. It determines the level of parallelism and can help speed up the conversion process
    :type threads: int
    :param genomes_folder: The `genomes_folder` parameter is a string that specifies the folder where
    the genome files are located. If this parameter is not provided, it will default to a constant value
    `DEFAULT_GENOME_FOLDER`
    :type genomes_folder: str
    :param to_parquet: The `to_parquet` parameter is a boolean value that specifies whether the HGMD
    database should be converted to the Parquet format or not. If set to `True`, the database will be
    converted to Parquet format. If set to `False`, the conversion will be skipped, defaults to True
    :type to_parquet: bool (optional)
    :param to_tsv: The `to_tsv` parameter is a boolean value that specifies whether the HGMD database
    should be converted to TSV format or not. If set to `True`, the function will generate a TSV file
    from the HGMD database. If set to `False`, the TSV conversion will be, defaults to True
    :type to_tsv: bool (optional)
    :return: a boolean value indicating whether the HGMD database conversion was successful or not.
    """

    # Check assemblies
    if not assemblies or len(assemblies) > 1 or len(assemblies) == 0:
        log.error("Uniq assembly is mandatory")
        raise ValueError("Uniq assembly is mandatory")

    # Full path
    hgmd_file = full_path(hgmd_file)
    hgmd_folder = full_path(hgmd_folder)
    

    # Check HGMD file
    if not hgmd_file or not os.path.exists(hgmd_file):
        log.error(f"HGMD file DOES NOT exist '{hgmd_file}'")
        raise ValueError(f"HGMD file DOES NOT exist '{hgmd_file}'")

    # Log
    log.info(f"Convert HGMD database {assemblies}")

    # genomes folder
    if not genomes_folder:
        genomes_folder = DEFAULT_GENOME_FOLDER

    # Full path
    genomes_folder = full_path(genomes_folder)

    # Assembly
    assembly = assemblies[0]

    # Create folder if not exists
    hgmd_folder_assembly = os.path.join(hgmd_folder,assembly)
    if not os.path.exists(hgmd_folder_assembly):
        Path(hgmd_folder_assembly).mkdir(parents=True, exist_ok=True)

    # Output basename
    if not output_basename:
        output_basename = os.path.basename(hgmd_file).replace(f'_{assembly}.vcf.gz', '')

    output_vcf = os.path.join(hgmd_folder_assembly, output_basename+".vcf.gz")

    if os.path.exists(output_vcf):

        log.info(f"Convert HGMD database {[assembly]} - File '{os.path.basename(output_vcf)}' already exists")

    else:

        log.info(f"Convert HGMD database {[assembly]} - File '{os.path.basename(hgmd_file)}' processing...")

        with TemporaryDirectory() as tmp_dir:
        
            # Import Database object
            from howard.objects.database import Database
            from howard.objects.variants import Variants

            # Database
            database = Database(database=hgmd_file, assembly=assembly)

            # Create duckDB connexion
            db_copy = duckdb.connect(config={"threads":threads})

            # Set max expression depth for big sub databases (e.g. gnomAD, ALL)
            db_copy.execute("SET max_expression_depth TO 10000")

            # Log
            log.debug(f"""Convert HGMD database {[assembly]} - Check chromosomes and number of variants...""")

            # Check chomosomes and total number of variants
            query_chromosomes = f"""
                        SELECT distinct concat('chr', "#CHROM") AS chromosome, count(*) AS nb_variants
                        FROM read_csv('{hgmd_file}', AUTO_DETECT=TRUE, ALL_VARCHAR=1, SEP='\t')
                        GROUP BY "#CHROM"
                    """
            chromosomes = sorted(list(db_copy.query(query_chromosomes).df()["chromosome"]))
            nb_variants = sum(db_copy.query(query_chromosomes).df()["nb_variants"])
            
            # Log
            log.debug(f"""Convert HGMD database {[assembly]} - Found {len(chromosomes)} chromosomes""")
            log.debug(f"""Convert HGMD database {[assembly]} - Found {nb_variants} variants""")

            log.debug(f"""Convert HGMD database {[assembly]} - Generate files...""")

            # Chromosome sizes
            genomes_sizes = {}
            genomes_sizes_file = os.path.join(genomes_folder, assembly, f"{assembly}.fa.sizes")
            if os.path.exists(genomes_sizes_file):
                with open(genomes_sizes_file, "r") as f:
                    for line in f:
                        genomes_sizes[line.split("\t")[0]] = line.split("\t")[1].strip()

            ### Header

            # Init
            header_list = []
            header_file = os.path.join(tmp_dir, "header.hdr")
            database_header = database.get_header_file()

            # Read header file to dict
            with mgzip.open(database_header, "r") as f:
                for line in f:
                    #log.debug(line)
                    if line.decode().startswith("#"):
                        header_list.append(line.decode().strip())
                    else:
                        break

            # Remove last line with #CHROM
            header_list = header_list[:-1]

            # Append INFO MC
            header_list.append('##INFO=<ID=MC,Number=1,Type=String,Description="HGMD Accession number">')

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

            ### Variants

            # Init
            variants_file = os.path.join(tmp_dir, "variants.vcf")

            # Query COPY TO
            query_variants = f"""
                COPY (
                        SELECT
                            concat('chr', "#CHROM") AS '#CHROM',
                            POS,
                            ID,
                            REF,
                            ALT,
                            QUAL,
                            FILTER,
                            concat(
                                replace(INFO, '"', ''),
                                ';',
                                'MC=',
                                ID
                                ) AS 'INFO'
                        FROM read_csv('{hgmd_file}', AUTO_DETECT=TRUE, ALL_VARCHAR=1, SEP='\t', SAMPLE_SIZE=10000000)
                    )
                TO '{variants_file}' (FORMAT csv, delimiter '\t', HEADER 0)
                    """
            db_copy.query(query_variants)
            
            # Create HGMD VCF
            log.info(f"""Convert HGMD database {[assembly]} - Generate VCF file '{os.path.basename(output_vcf)}'""")
            concat_and_compress_files(input_files=[header_file, variants_file], output_file=output_vcf, sort=True, index=True, compression_type = "bgzip")

            # Create HGMD Parquet
            if to_parquet or to_tsv:
                output_parquet = os.path.join(hgmd_folder_assembly, output_basename+".parquet")
                log.info(f"""Convert HGMD database {[assembly]} - Generate Parquet file '{os.path.basename(output_parquet)}'""")
                hgmd_database_to_parquet = Variants(input=output_vcf, output=output_parquet, config={"threads":threads}, param={"explode_infos": True}, load=False)
                hgmd_database_to_parquet.load_data(sample_size=-1)
                hgmd_database_to_parquet.export_output(export_header=True)

            # Create HGMD TSV
            if to_tsv:
                output_tsv = os.path.join(hgmd_folder_assembly, output_basename+".tsv")
                log.info(f"""Convert HGMD database {[assembly]} - Generate TSV file '{os.path.basename(output_tsv)}'""")
                hgmd_database_to_tsv = Variants(input=output_parquet, output=output_tsv, config={"threads":threads}, param={"explode_infos": True}, load=False)
                hgmd_database_to_tsv.load_data(sample_size=-1)
                hgmd_database_to_tsv.export_output(export_header=True)


    return True


