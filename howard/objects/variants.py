import csv
import gc
import gzip
import io
import multiprocessing
import os
import random
import re
import shlex
import sqlite3
import subprocess
from tempfile import NamedTemporaryFile, TemporaryDirectory
import tempfile
import duckdb
import json
import yaml
import argparse
import Bio.bgzf as bgzf
import pandas as pd
from pyfaidx import Fasta
import numpy as np
import vcf
import logging as log
import fastparquet as fp
from multiprocesspandas import applyparallel

from howard.commons import *
from howard.objects.database import *
from howard.tools.databases import *
from howard.utils import *


class Variants:

    def __init__(self, conn=None, input: str = None, output: str = None, config: dict = {}, param: dict = {}, load: bool = False) -> None:
        """
        The function `__init__` initializes the variables, sets the input, output, config, param, connexion and
        header

        :param conn: the connection to the database
        :param input: the input file
        :param output: the output file
        :param config: a dictionary containing the configuration of the model
        :param param: a dictionary containing the parameters of the model
        """

        # Init variables
        self.init_variables()

        # Input
        self.set_input(input)

        # Config
        self.set_config(config)

        # Param
        self.set_param(param)

        # Output
        self.set_output(output)

        # connexion
        self.set_connexion(conn)

        # Header
        self.set_header()

        # Load data
        if load:
            self.load_data()


    def set_input(self, input: str = None) -> None:
        """
        The function takes a file name as input, splits the file name into a name and an extension, and
        then sets the input_name, input_extension, and input_format attributes of the class

        :param input: The input file
        """
        self.input = input
        
        # Input format
        if input:
            input_name, input_extension = os.path.splitext(self.input)
            self.input_name = input_name
            self.input_extension = input_extension
            self.input_format = self.input_extension.replace(".", "")


    def set_config(self, config: dict) -> None:
        """
        This function takes in a config object and sets it as the config object for the class

        :param config: The configuration object
        """
        self.config = config


    def set_param(self, param: dict) -> None:
        """
        This function takes in a param object and sets it as the param object for the class

        :param param: The paramters object
        """
        self.param = param


    def init_variables(self) -> None:
        """
        This function initializes the variables that will be used in the rest of the class
        """
        self.prefix = "howard"
        self.table_variants = "variants"
        self.dataframe = None

        self.comparison_map = {
            "gt": ">",
            "gte": ">=",
            "lt": "<",
            "lte": "<=",
            "equals": "=",
            "contains": "SIMILAR TO"
        }

        self.code_type_map = {
            "Integer": 0,
            "String": 1,
            "Float": 2,
            "Flag": 3
        }

        self.code_type_map_to_sql = {
            "Integer": "INTEGER",
            "String": "VARCHAR",
            "Float": "FLOAT",
            "Flag": "VARCHAR"
        }

        self.index_additionnal_fields = []


    def get_indexing(self) -> bool:
        """
        It returns the value of the key "indexing" in the dictionary. If the key is not present, it
        returns False.
        :return: The value of the indexing parameter.
        """
        return self.get_param().get("indexing", False)


    def get_connexion_config(self) -> dict:
        """
        The function `get_connexion_config` returns a dictionary containing the configuration for a
        connection, including the number of threads and memory limit.
        :return: a dictionary containing the configuration for the Connexion library.
        """

        # config
        config = self.get_config()

        # Connexion config
        connexion_config = {}
        threads = self.get_threads()

        # Threads
        if threads:
            connexion_config["threads"] = threads

        # Memory
        if config.get("memory", None):
            connexion_config["memory_limit"] = config.get("memory")

        # Temporary directory
        if config.get("tmp", None):
            connexion_config["temp_directory"] = config.get("tmp")
        
        # Access
        if config.get("access", None):
            access = config.get("access")
            if access in ["RO"]:
                access = "READ_ONLY"
            elif access in ["RW"]:
                access = "READ_WRITE"
            connexion_db = self.get_connexion_db()
            if connexion_db in ":memory:":
                access = "READ_WRITE"
            connexion_config["access_mode"] = access

        return connexion_config

    def get_duckdb_settings(self) -> dict:
        """
        The function `get_duckdb_settings` retrieves DuckDB settings from a configuration file or a
        string.
        :return: The function `get_duckdb_settings` returns a dictionary object `duckdb_settings_dict`.
        """

        # config
        config = self.get_config()

        # duckdb settings
        duckdb_settings_dict = {}
        if config.get("duckdb_settings", None):
            duckdb_settings = config.get("duckdb_settings")
            # duckdb setting is a file
            if os.path.exists(duckdb_settings):
                with open(duckdb_settings) as json_file:
                    duckdb_settings_dict = yaml.safe_load(json_file)
            # duckdb settings is a string
            else:
                duckdb_settings_dict = json.loads(duckdb_settings)

        return duckdb_settings_dict


    def set_connexion_db(self) -> str:
        """
        The function `set_connexion_db` returns the appropriate database connection string based on the
        input format and connection type.
        :return: the value of the variable `connexion_db`.
        """

        # Default connexion db
        default_connexion_db = ":memory:"

        # Find connexion db
        if self.get_input_format() in ["db", "duckdb"]:
            connexion_db = self.get_input()
        elif self.get_connexion_type() in ["memory", default_connexion_db, None]:
            connexion_db = default_connexion_db
        elif self.get_connexion_type() in ["tmpfile"]:
            tmp_name = tempfile.mkdtemp(prefix=self.get_prefix(
            ), dir=self.get_tmp_dir(), suffix=".db")
            connexion_db = f"{tmp_name}/tmp.db"
        elif self.get_connexion_type() != "":
            connexion_db = self.get_connexion_type()
        else:
            connexion_db = default_connexion_db
        
        # Set connexion db
        self.connexion_db = connexion_db

        return connexion_db


    def set_connexion(self, conn) -> None:
        """
        It creates a connection to the database

        :param conn: The connection to the database. If not provided, a new connection to an in-memory
        database is created
        """

        # Connexion db
        connexion_db = self.set_connexion_db()

        # Connexion config
        connexion_config = self.get_connexion_config()

        # Connexion format
        connexion_format = self.get_config().get("connexion_format", "duckdb")
        # Set connexion format
        self.connexion_format = connexion_format

        # Connexion
        if not conn:
            if connexion_format in ["duckdb"]:
                conn = duckdb.connect(connexion_db, config=connexion_config)
                # duckDB settings
                duckdb_settings = self.get_duckdb_settings()
                if duckdb_settings:
                    for setting in duckdb_settings:
                        setting_value = duckdb_settings.get(setting)
                        if isinstance(setting_value, str):
                            setting_value = f"'{setting_value}'"
                        conn.execute(f"PRAGMA {setting}={setting_value};")
            elif connexion_format in ["sqlite"]:
                conn = sqlite3.connect(connexion_db)

        # Set connexion
        self.conn = conn

        # Log
        log.debug(f"connexion_format: {connexion_format}")
        log.debug(f"connexion_db: {connexion_db}")
        log.debug(f"connexion config: {connexion_config}")
        log.debug(f"connexion duckdb settings: {self.get_duckdb_settings()}")


    def set_output(self, output: str = None) -> None:
        """
        If the config file has an output key, set the output to the value of that key. Otherwise, set
        the output to the input

        :param output: The name of the output file
        """
        self.output = output

        # Output format
        if self.output:
            output_name, output_extension = os.path.splitext(self.output)
            self.output_name = output_name
            self.output_extension = output_extension
            self.output_format = self.output_extension.replace(".", "")
        else:
            self.output_name = None
            self.output_extension = None
            self.output_format = None


    def set_header(self) -> None:
        """
        It reads the header of a VCF file and stores it as a list of strings and as a VCF object
        """
        input_file = self.get_input()
        default_header_list = [
            '##fileformat=VCFv4.2',
            '#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO'
            ]

        if input_file:

            input_format = self.get_input_format()
            input_compressed = self.get_input_compressed()
            config = self.get_config()
            header_list = default_header_list
            if input_format in ["vcf", "hdr", "tsv", "csv", "psv", "parquet", "db", "duckdb"]:
                # header provided in param
                if config.get("header_file", None):
                    with open(config.get("header_file"), 'rt') as f:
                        header_list = self.read_vcf_header(f)
                # within a vcf file format (header within input file itsself)
                elif input_format in ["vcf", "hdr"]:
                    # within a compressed vcf file format (.vcf.gz)
                    if input_compressed:
                        with bgzf.open(input_file, 'rt') as f:
                            header_list = self.read_vcf_header(f)
                    # within an uncompressed vcf file format (.vcf)
                    else:
                        with open(input_file, 'rt') as f:
                            header_list = self.read_vcf_header(f)
                # header provided in default external file .hdr
                elif os.path.exists((input_file+".hdr")):
                    with open(input_file+".hdr", 'rt') as f:
                        header_list = self.read_vcf_header(f)
                else:
                    try: # Try to get header info fields and file columns

                        with tempfile.TemporaryDirectory() as tmpdir:

                            # Create database
                            db_for_header = Database(database=input_file)
                            
                            # Get header columns for infos fields
                            db_header_from_columns = db_for_header.get_header_from_columns()

                            # Get real columns in the file
                            db_header_columns = db_for_header.get_columns()

                            # Write header file
                            header_file_tmp = os.path.join(tmpdir,'header')
                            f = open(header_file_tmp, 'w')
                            vcf.Writer(f, db_header_from_columns)
                            f.close()
                            
                            # Replace #CHROM line with rel columns
                            header_list = db_for_header.read_header_file(header_file=header_file_tmp)
                            header_list[-1] = "\t".join(db_header_columns)

                    except:

                        log.warning(f"No header for file {input_file}. Set as default VCF header")
                        header_list = default_header_list

            else:  # try for unknown format ?

                log.error(f"Input file format '{input_format}' not available")
                raise ValueError(
                    f"Input file format '{input_format}' not available")

            if not header_list:
                header_list = default_header_list

            # header as list
            self.header_list = header_list

            # header as VCF object
            self.header_vcf = vcf.Reader(io.StringIO("\n".join(header_list)))

        else:

            self.header_list = None
            self.header_vcf = None


    def get_query_to_df(self, query: str = "", limit:int = None) -> pd.DataFrame:
        """
        > The function `get_query_to_df` takes a query as a string and returns a pandas dataframe

        :param query: str = ""
        :type query: str
        :return: A dataframe
        """

        # Connexion format
        connexion_format = self.get_connexion_format()

        # Limit in query
        if limit:
            pd.set_option('display.max_rows', limit)
            if connexion_format in ["duckdb"]:
                df = self.conn.execute(query).fetch_record_batch(limit).read_next_batch().to_pandas()
            elif connexion_format in ["sqlite"]:
                df = next(pd.read_sql_query(query, self.conn, chunksize=limit))
                
        # Full query
        else:
            if connexion_format in ["duckdb"]:
                df = self.conn.execute(query).df()
            elif connexion_format in ["sqlite"]:
                df = pd.read_sql_query(query, self.conn)

        return df


    def get_overview(self) -> None:
        """
        The function prints the input, output, config, and dataframe of the current object
        """
        table_variants_from = self.get_table_variants(clause="from")
        sql_columns = self.get_header_columns_as_sql()
        sql_query_export = f"SELECT {sql_columns} FROM {table_variants_from}"
        df = self.get_query_to_df(sql_query_export)
        log.info("Input:  " + str(self.get_input()) +
                 " [" + str(str(self.get_input_format())) + "]")
        log.info("Output: " + str(self.get_output()) +
                 " [" + str(str(self.get_output_format())) + "]")
        log.info("Config: ")
        for d in str(json.dumps(self.get_config(), indent=4, sort_keys=True)).split("\n"):
            log.info("\t" + str(d))
        log.info("Param: ")
        for d in str(json.dumps(self.get_param(), indent=4, sort_keys=True)).split("\n"):
            log.info("\t" + str(d))
        log.info("Sample list: " + str(self.get_header_sample_list()))
        log.info("Dataframe: ")
        for d in str(df).split("\n"):
            log.info("\t" + str(d))

        # garbage collector
        del df
        gc.collect()

        return None


    def get_stats(self) -> dict:
        """
        The `get_stats` function calculates and returns various statistics of the current object,
        including information about the input file, variants, samples, header fields, quality, and
        SNVs/InDels.
        :return: a dictionary containing various statistics of the current object. The dictionary has
        the following structure:
        """

        # Log
        log.info(f"Stats Calculation...")

        # table varaints
        table_variants_from = self.get_table_variants()

        # stats dict
        stats = {
            "Infos": {}
        }

        ### File
        input_file = self.get_input()
        stats["Infos"]["Input file"] = input_file

        # Header
        header_infos = self.get_header().infos
        header_formats = self.get_header().formats
        header_infos_list = list(header_infos)
        header_formats_list = list(header_formats)


        ### Variants
        
        stats["Variants"] = {}

        # Variants by chr
        sql_query_nb_variant_by_chrom = f"SELECT \"#CHROM\" as CHROM, count(*) as count FROM {table_variants_from} GROUP BY \"#CHROM\""
        df_nb_of_variants_by_chrom = self.get_query_to_df(sql_query_nb_variant_by_chrom)
        nb_of_variants_by_chrom = df_nb_of_variants_by_chrom.sort_values(by=['CHROM'], kind='quicksort')
        
        # Total number of variants
        nb_of_variants = nb_of_variants_by_chrom["count"].sum()

        # Calculate percentage
        nb_of_variants_by_chrom['percent'] = nb_of_variants_by_chrom['count'].apply(lambda x: (x / nb_of_variants))

        stats["Variants"]["Number of variants by chromosome"] = nb_of_variants_by_chrom.to_dict(orient="index")

        stats["Infos"]["Number of variants"] = int(nb_of_variants)


        ### Samples
        nb_of_samples = len(self.get_header_sample_list())
        if "FORMAT" in self.get_header_columns() and "DP" in header_formats_list:
            stats["Infos"]["Number of samples"] = nb_of_samples
        elif nb_of_samples:
            stats["Infos"]["Number of samples"] = "not a VCF format"

        samples = {}
        
        if "GT" in header_formats_list and "FORMAT" in self.get_header_columns():
            i = 0
            for sample in self.get_header_sample_list():
                i += 1
                sql_query_samples = f"""
                    SELECT  '{sample}' as sample,
                            REGEXP_EXTRACT("{sample}", '^([0-9/|.]*)[:]*',1) as genotype,
                            count(REGEXP_EXTRACT("{sample}", '^([0-9/|.]*)[:]*',1)) as count,
                            concat((count(REGEXP_EXTRACT("{sample}", '^([0-9/|.]*)[:]*',1))/{nb_of_variants})) as percentage
                    FROM {table_variants_from}
                    WHERE (
                        regexp_matches("{sample}", '^[0-9]([/|][0-9])+')
                      )
                    GROUP BY genotype
                    """
                if i > 5:
                    break
                sql_query_genotype_df = self.conn.execute(sql_query_samples).df()
                sample_genotype_count = sql_query_genotype_df["count"].sum()
                if len(sql_query_genotype_df):
                    samples[f"{sample} - {sample_genotype_count} variants"] = sql_query_genotype_df.to_dict(orient="index")
            
            stats["Samples"] = samples


        ### INFO and FORMAT fields
        header_types_df = {}
        header_types_list = {"List of INFO fields": header_infos, "List of FORMAT fields": header_formats}
        i = 0
        for header_type in header_types_list:
            
            header_type_infos = header_types_list.get(header_type)
            header_infos_dict = {}

            for info in header_type_infos:

                i += 1
                header_infos_dict[i] = {}

                # ID
                header_infos_dict[i]["id"] = info

                # num
                genotype_map = {
                    None: ".",
                    -1: "A",
                    -2: "G",
                    -3: "R"
                }
                if header_type_infos[info].num in genotype_map.keys():
                    header_infos_dict[i]["Number"] = genotype_map.get(header_type_infos[info].num)
                else:
                    header_infos_dict[i]["Number"] = header_type_infos[info].num

                # type
                if header_type_infos[info].type:
                    header_infos_dict[i]["Type"] = header_type_infos[info].type
                else:
                    header_infos_dict[i]["Type"] = "."

                # desc
                if header_type_infos[info].desc != None:
                    header_infos_dict[i]["Description"] = header_type_infos[info].desc
                else:
                    header_infos_dict[i]["Description"] = ""

            if len(header_infos_dict):
                header_types_df[header_type] = pd.DataFrame.from_dict(header_infos_dict, orient="index").to_dict(orient="index")

        # Stats
        stats["Infos"]["Number of INFO fields"] = len(header_infos_list)
        stats["Infos"]["Number of FORMAT fields"] = len(header_formats_list)
        stats["Header"] = header_types_df


        ### QUAL
        if "QUAL" in self.get_header_columns():
            sql_query_qual = f"""
                    SELECT
                        avg(CAST(QUAL AS INTEGER)) AS Average,
                        min(CAST(QUAL AS INTEGER)) AS Minimum,
                        max(CAST(QUAL AS INTEGER)) AS Maximum,
                        stddev(CAST(QUAL AS INTEGER)) AS StandardDeviation,
                        median(CAST(QUAL AS INTEGER)) AS Median,
                        variance(CAST(QUAL AS INTEGER)) AS Variance
                    FROM {table_variants_from}
                    WHERE QUAL NOT IN ('.')
                    """
            
            qual = self.conn.execute(sql_query_qual).df().to_dict(orient="index")
            stats["Quality"] = {
                "Stats": qual
            }

        ### SNV and InDel

        sql_query_snv = f"""
            
            SELECT Type, count FROM (

                    SELECT
                        'Total' AS Type,
                        count(*) AS count
                    FROM {table_variants_from}

                    UNION

                    SELECT
                        'MNV' AS Type,
                        count(*) AS count
                    FROM {table_variants_from}
                    WHERE len(REF) > 1 AND len(ALT) > 1
                    AND len(REF) = len(ALT)

                    UNION

                    SELECT
                        'InDel' AS Type,
                        count(*) AS count
                    FROM {table_variants_from}
                    WHERE len(REF) > 1 OR len(ALT) > 1
                    AND len(REF) != len(ALT)
                    
                    UNION

                    SELECT
                        'SNV' AS Type,
                        count(*) AS count
                    FROM {table_variants_from}
                    WHERE len(REF) = 1 AND len(ALT) = 1

                )

            ORDER BY count DESC

                """
        snv_indel = self.conn.execute(sql_query_snv).df().to_dict(orient="index")

        sql_query_snv_substitution = f"""
                SELECT
                    concat(REF, '>', ALT) AS 'Substitution',
                    count(*) AS count
                FROM {table_variants_from}
                WHERE len(REF) = 1 AND len(ALT) = 1
                GROUP BY REF, ALT
                ORDER BY count(*) DESC
                """
        snv_substitution = self.conn.execute(sql_query_snv_substitution).df().to_dict(orient="index")
        stats["Variants"]["Counts"] = snv_indel
        stats["Variants"]["Substitutions"] = snv_substitution

        return stats
    

    def stats_to_file(self, file:str = None) -> str:
        """
        The function `stats_to_file` takes a file name as input, retrieves statistics, serializes them
        into a JSON object, and writes the JSON object to the specified file.
        
        :param file: The `file` parameter is a string that represents the file path where the JSON data
        will be written
        :type file: str
        :return: the name of the file that was written to.
        """

        # Get stats
        stats = self.get_stats()

        # Serializing json
        json_object = json.dumps(stats, indent=4)
        
        # Writing to sample.json
        with open(file, "w") as outfile:
            outfile.write(json_object)

        return file


    def print_stats(self, output_file:str = None, json_file:str = None) -> None:
        """
        The `print_stats` function generates a markdown file and prints the statistics contained in a
        JSON file in a formatted manner.
        
        :param output_file: The `output_file` parameter is a string that specifies the path and filename
        of the output file where the stats will be printed in Markdown format. If no `output_file` is
        provided, a temporary directory will be created and the stats will be saved in a file named
        "stats.md" within that
        :type output_file: str
        :param json_file: The `json_file` parameter is a string that represents the path to the JSON
        file where the statistics will be saved. If no value is provided, a temporary directory will be
        created and a default file name "stats.json" will be used
        :type json_file: str
        :return: The function `print_stats` does not return any value. It has a return type annotation
        of `None`.
        """

        with tempfile.TemporaryDirectory() as tmpdir:

            # Files
            if not output_file:
                output_file = os.path.join(tmpdir, "stats.md")
            if not json_file:
                json_file = os.path.join(tmpdir, "stats.json")

            # Create folders
            if not os.path.exists(os.path.dirname(output_file)):
                Path(os.path.dirname(output_file)).mkdir(parents=True, exist_ok=True)
            if not os.path.exists(os.path.dirname(json_file)):
                Path(os.path.dirname(json_file)).mkdir(parents=True, exist_ok=True)

            # Create stats JSON file
            stats_file = self.stats_to_file(file=json_file)

            # Print stats file
            with open(stats_file) as f:
                stats = yaml.safe_load(f)

            # Output 
            output_title = []
            output_index = []
            output = []

            # Title
            output_title.append("# HOWARD Stats")

            # Index
            output_index.append("## Index")

            # Process sections
            for section in stats:
                infos = stats.get(section)
                section_link = "#" + section.lower().replace(" ", "-")
                output.append(f"## {section}")
                output_index.append(f"- [{section}]({section_link})")
                
                if len(infos):
                    for info in infos:
                        try:
                            df = pd.DataFrame.from_dict(infos.get(info), orient='index')
                            is_df = True
                        except:
                            try:
                                df = pd.DataFrame.from_dict(json.loads((infos.get(info))), orient="index")
                                is_df = True
                            except:
                                is_df = False
                        if is_df:
                            output.append(f"### {info}")
                            info_link = "#" + info.lower().replace(" ", "-")
                            output_index.append(f"   - [{info}]({info_link})")
                            output.append(f"{df.to_markdown(index=False)}")
                        else:
                            output.append(f"- {info}: {infos.get(info)}")
                else:
                    output.append(f"NA")

            # Write stats in markdown file
            with open(output_file, 'w') as fp:
                for item in output_title:
                    fp.write("%s\n" % item)
                for item in output_index:
                    fp.write("%s\n" % item)
                for item in output:
                    fp.write("%s\n" % item)

            # Output stats in markdown
            print("")
            print("\n\n".join(output_title))
            print("")
            print("\n\n".join(output))
            print("")

        return None


    def get_input(self) -> str:
        """
        It returns the value of the input variable.
        :return: The input is being returned.
        """
        return self.input


    def get_input_format(self, input_file: str = None) -> str:
        """
        It returns the format of the input variable.
        :return: The format is being returned.
        """
        if not input_file:
            input_file = self.get_input()
        input_format = get_file_format(input_file)
        return input_format


    def get_input_compressed(self, input_file: str = None) -> str:
        """
        It returns the format of the input variable.
        :return: The format is being returned.
        """
        if not input_file:
            input_file = self.get_input()
        input_compressed = get_file_compressed(input_file)
        return input_compressed


    def get_output(self) -> str:
        """
        It returns the output of the neuron.
        :return: The output of the neural network.
        """
        return self.output


    def get_output_format(self, output_file: str = None) -> str:
        """
        It returns the format of the input variable.
        :return: The format is being returned.
        """
        if not output_file:
            output_file = self.get_output()
        output_format = get_file_format(output_file)

        return output_format


    def get_config(self) -> dict:
        """
        It returns the config
        :return: The config variable is being returned.
        """
        return self.config


    def get_param(self) -> dict:
        """
        It returns the param
        :return: The param variable is being returned.
        """
        return self.param


    def get_connexion_db(self) -> str:
        """
        It returns the connexion_db attribute of the object
        :return: The connexion_db is being returned.
        """
        return self.connexion_db


    def get_prefix(self) -> str:
        """
        It returns the prefix of the object.
        :return: The prefix is being returned.
        """
        return self.prefix


    def get_table_variants(self, clause: str = "select") -> str:
        """
        This function returns the table_variants attribute of the object

        :param clause: the type of clause the table will be used. Either "select" or "from" (optional),
        defaults to select (optional)
        :return: The table_variants attribute of the object.
        """

        # Access
        access = self.get_config().get("access", None)

        # Clauses "select", "where", "update"
        if clause in ["select", "where", "update"]:
            table_variants = self.table_variants
        # Clause "from"
        elif clause in ["from"]:
            # For Read Only
            if self.get_input_format() in ["parquet"] and access in ["RO"]:
                input_file = self.get_input()
                table_variants = f"'{input_file}' as variants"
            # For Read Write
            else:
                table_variants = f"{self.table_variants} as variants"
        else:
            table_variants = self.table_variants
        return table_variants


    def get_tmp_dir(self) -> str:
        """
        It returns the value of the tmp_dir key in the config dictionary, or /tmp if the key doesn't
        exist

        :return: The value of the key "tmp_dir" in the config file.
        """
        return self.get_config().get("tmp_dir", "/tmp")


    def get_connexion_type(self) -> str:
        """
        If the connexion type is not in the list of allowed connexion types, raise a ValueError

        :return: The connexion type is being returned.
        """
        return self.get_config().get("connexion_type", "memory")


    def get_connexion(self):
        """
        It returns the connection object

        :return: The connection object.
        """
        return self.conn


    def close_connexion(self) -> None:
        """
        This function closes the connection to the database.
        :return: The connection is being closed.
        """
        return self.conn.close()


    def get_header(self, type: str = "vcf"):
        """
        This function returns the header of the VCF file as a list of strings

        :param type: the type of header you want to get, defaults to vcf (optional)
        :return: The header of the vcf file.
        """
        
        if self.header_vcf:
            if type == "vcf":
                return self.header_vcf
            elif type == "list":
                return self.header_list
        else:
            if type == "vcf":
                header = vcf.Reader(io.StringIO("\n".join(vcf_required)))
                return header
            elif type == "list":
                return vcf_required


    def get_header_length(self, file:str = None) -> int:
        """
        The function `get_header_length` returns the length of the header list, excluding the #CHROM
        line.
        
        :param file: The `file` parameter is an optional argument that specifies the path to a VCF
        header file. If this argument is provided, the function will read the header from the specified
        file and return the length of the header list minus 1 (to exclude the #CHROM line)
        :type file: str
        :return: the length of the header list, excluding the #CHROM line.
        """

        if file:
            return len(self.read_vcf_header_file(file=file)) -1 
        elif self.get_header(type="list"):
            return len(self.get_header(type="list")) -1
        else:
            return 0


    def get_header_columns(self) -> str:
        """
        This function returns the header list of a VCF

        :return: The length of the header list.
        """
        if self.get_header():
            return self.get_header(type="list")[-1]
        else:
            return ""


    def get_header_columns_as_list(self) -> list:
        """
        This function returns the header list of a VCF

        :return: The length of the header list.
        """
        if self.get_header():
            return self.get_header_columns().strip().split("\t")
        else:
            return []


    def get_header_columns_as_sql(self) -> str:
        """
        This function retruns header length (without #CHROM line)

        :return: The length of the header list.
        """
        sql_column_list = []
        for col in self.get_header_columns_as_list():
            sql_column_list.append(f"\"{col}\"")
        return ",".join(sql_column_list)


    def get_header_sample_list(self) -> list:
        """
        This function retruns header length (without #CHROM line)

        :return: The length of the header list.
        """
        return self.header_vcf.samples


    def get_verbose(self) -> bool:
        """
        It returns the value of the "verbose" key in the config dictionary, or False if the key doesn't
        exist

        :return: The value of the key "verbose" in the config dictionary.
        """
        return self.get_config().get("verbose", False)


    def get_connexion_format(self) -> str:
        """
        It returns the connexion format of the object.
        :return: The connexion_format is being returned.
        """
        connexion_format = self.connexion_format
        if connexion_format not in ["duckdb", "sqlite"]:
            log.error(f"Unknown connexion format {connexion_format}")
            raise ValueError(f"Unknown connexion format {connexion_format}")
        else:
            return connexion_format


    def insert_file_to_table(self, file, columns: str, header_len: int = 0, sep: str = '\t', chunksize: int = 1000000) -> None:
        """
        The function reads a file in chunks, and inserts each chunk into a table

        :param file: the file to be loaded
        :param columns: a string of the column names separated by commas
        :param header_len: the number of lines to skip at the beginning of the file, defaults to 0
        (optional)
        :param sep: the separator used in the file, defaults to \t (optional)
        :param chunksize: The number of rows to read in at a time, defaults to 1000000 (optional)
        """

        # Config
        chunksize = self.get_config().get("load", {}).get("chunk", chunksize)
        connexion_format = self.get_connexion_format()

        log.debug("chunksize: "+str(chunksize))

        if chunksize:
            for chunk in pd.read_csv(file, skiprows=header_len, sep=sep, chunksize=chunksize, engine="c"):
                if connexion_format in ["duckdb"]:
                    sql_insert_into = f"INSERT INTO variants ({columns}) SELECT {columns} FROM chunk"
                    self.conn.execute(sql_insert_into)
                elif connexion_format in ["sqlite"]:
                    chunk.to_sql("variants", self.conn,
                                 if_exists='append', index=False)


    def load_data(self, input_file:str = None, drop_variants_table:bool = False, sample_size:int = 20480) -> None:
        """
        The `load_data` function reads a VCF file and inserts it into a table, with options to drop the
        table before loading the data and specify a sample size.
        
        :param input_file: The path to the input file. This is the VCF file that will be loaded into the
        table
        :type input_file: str
        :param drop_variants_table: The `drop_variants_table` parameter is a boolean flag that
        determines whether the variants table should be dropped before loading the data. If set to
        `True`, the variants table will be dropped. If set to `False` (default), the variants table will
        not be dropped, defaults to False
        :type drop_variants_table: bool (optional)
        :param sample_size: The `sample_size` parameter determines the number of rows to be sampled from
        the input file. If it is set to `None`, the default value of 20480 will be used, defaults to
        20480
        :type sample_size: int (optional)
        """

        log.info("Loading...")

        # change input file
        if input_file:
            self.set_input(input_file)
            self.set_header()

        # drop variants table
        if drop_variants_table:
            self.drop_variants_table()

        # get table variants
        table_variants = self.get_table_variants()

        # Access
        access = self.get_config().get("access", None)
        log.debug(f"access: {access}")

        # Input format and compress
        input_format = self.get_input_format()
        input_compressed = self.get_input_compressed()
        log.debug(f"input_format: {input_format}")
        log.debug(f"input_compressed: {input_compressed}")

        # input_compressed_format
        if input_compressed:
            input_compressed_format = "gzip"
        else:
            input_compressed_format = "none"
        log.debug(f"input_compressed_format: {input_compressed_format}")

        # Connexion format
        connexion_format = self.get_connexion_format()

        # Sample size
        if not sample_size:
            sample_size = -1
        log.debug(f"sample_size: {sample_size}")

        # Load data
        log.debug(f"Load Data from {input_format}")
        if input_format in ["vcf", "tsv", "csv", "psv"]:

            # delimiter
            delimiter = file_format_delimiters.get(input_format, "\t")

            # Skip header lines
            skip = self.get_header_length(file=self.input)

            if connexion_format in ["duckdb"]:
                if access in ["RO"]:
                    sql_vcf = f"""
                    CREATE VIEW {table_variants} AS
                        SELECT *
                        FROM read_csv('{self.input}', auto_detect=True, skip={skip}, delim='{delimiter}', sample_size={sample_size}, compression={input_compressed_format})
                    """
                else:
                    sql_vcf = f"""
                    CREATE TABLE {table_variants} AS 
                        SELECT *
                        FROM read_csv('{self.input}', auto_detect=True, skip={skip}, delim='{delimiter}', sample_size={sample_size}, compression={input_compressed_format})
                    """
                self.conn.execute(sql_vcf)

            elif connexion_format in ["sqlite"]:

                # Main structure
                structure = {
                    "#CHROM": "VARCHAR",
                    "POS": "INTEGER",
                    "ID": "VARCHAR",
                    "REF": "VARCHAR",
                    "ALT": "VARCHAR",
                    "QUAL": "VARCHAR",
                    "FILTER": "VARCHAR",
                    "INFO": "VARCHAR",
                }

                # Strcuture with samples
                structure_complete = structure
                if self.get_header_sample_list():
                    structure["FORMAT"] = "VARCHAR"
                    for sample in self.get_header_sample_list():
                        structure_complete[sample] = "VARCHAR"

                # Columns list for create and insert
                sql_create_table_columns = []
                sql_create_table_columns_list = []
                for column in structure_complete:
                    column_type = structure_complete[column]
                    sql_create_table_columns.append(f"\"{column}\" {column_type} default NULL")
                    sql_create_table_columns_list.append(f"\"{column}\"")

                # Create database
                log.debug(f"Create Table {table_variants}")
                sql_create_table_columns_sql = ", ".join(sql_create_table_columns)
                sql_create_table_columns_list_sql = ", ".join(
                    sql_create_table_columns_list)
                sql_create_table = f"CREATE TABLE IF NOT EXISTS {table_variants} ({sql_create_table_columns_sql})"
                self.conn.execute(sql_create_table)

                # chunksize define length of file chunk load file
                chunksize = 100000

                # delimiter
                delimiter = file_format_delimiters.get(input_format, "\t")

                # Load the input file
                with open(self.input, "rt") as input_file:
                    # Use the appropriate file handler based on the input format
                    if input_compressed:
                        input_file = bgzf.open(self.input, "rt")
                    if input_format in ["vcf"]:
                        header_len = self.get_header_length()
                    else:
                        header_len = 0

                    # Insert the file contents into a table
                    self.insert_file_to_table(
                        input_file,
                        columns=sql_create_table_columns_list_sql,
                        header_len=header_len,
                        sep=delimiter,
                        chunksize=chunksize,
                    )

        elif self.input_format in ["parquet"]:

            # Load Parquet
            if connexion_format in ["duckdb"]:

                if os.path.isdir(self.input):
                    list_of_parquet = glob.glob(os.path.join(self.input,"**/*parquet"), recursive=True)
                    if list_of_parquet:
                        list_of_parquet_level_path = "*/" * (list_of_parquet[0].replace(self.input, "").count('/')-1)
                        sql_form = f"read_parquet('{self.input}/{list_of_parquet_level_path}*parquet', hive_partitioning=1)"
                    else:
                        log.error(f"Input file '{self.input}' not a compatible partitionned parquet folder")
                        raise ValueError(f"Input file '{self.input}' not a compatible partitionned parquet folder")
                else:
                    sql_form = f"read_parquet('{self.input}')"

                if access in ["RO"]:
                    sql_parquet = f"CREATE VIEW {table_variants} AS SELECT * FROM {sql_form}"
                else:
                    sql_parquet = f"CREATE TABLE {table_variants} AS SELECT * FROM {sql_form}"

                self.conn.execute(sql_parquet)

            else:
                log.error(f"Input file format '{self.input_format}' not compatilbe with database format '{connexion_format}'")
                raise ValueError(f"Input file format '{self.input_format}' not compatilbe with database format '{connexion_format}'")


        elif self.input_format in ["db", "duckdb"]:

            if connexion_format in ["duckdb"]:
                log.debug(f"Input file format '{self.input_format}' duckDB")
            else:
                log.error(f"Input file format '{self.input_format}' not compatilbe with database format '{connexion_format}'")
                raise ValueError(f"Input file format '{self.input_format}' not compatilbe with database format '{connexion_format}'")

        else:
            log.error(f"Input file format '{self.input_format}' not available")
            raise ValueError(f"Input file format '{self.input_format}' not available")

        # Explode INFOS fields into table fields
        if self.get_explode_infos():
            self.explode_infos(prefix=self.get_explode_infos_prefix(), fields=self.get_explode_infos_fields(), force=True)

        # Create index after insertion
        self.create_indexes()


    def get_explode_infos(self) -> bool:
        """
        The function `get_explode_infos` returns the value of the "explode_infos" parameter, defaulting
        to False if it is not set.
        :return: The method is returning the value of the "explode_infos" parameter, which is a boolean
        value. If the parameter is not present, it will return False.
        """

        return self.get_param().get("explode_infos", False)


    def get_explode_infos_fields(self, explode_infos_fields:str = None, remove_fields_not_in_header:bool = False) -> list:
        """
        The `get_explode_infos_fields` function returns a list of exploded information fields based on
        the input parameter `explode_infos_fields`.
        
        :param explode_infos_fields: The `explode_infos_fields` parameter is a string that specifies the
        fields to be exploded. It can be set to "ALL" to explode all fields, or it can be a
        comma-separated list of field names to explode
        :type explode_infos_fields: str
        :param remove_fields_not_in_header: The parameter `remove_fields_not_in_header` is a boolean
        flag that determines whether to remove fields that are not present in the header. If it is set
        to `True`, any field that is not in the header will be excluded from the list of exploded
        information fields. If it is set to `, defaults to False
        :type remove_fields_not_in_header: bool (optional)
        :return: The function `get_explode_infos_fields` returns a list of exploded information fields.
        If the `explode_infos_fields` parameter is not provided or is set to None, it returns an empty
        list. If the parameter is provided and its value is "ALL", it also returns an empty list.
        Otherwise, it returns a list of exploded information fields after removing any spaces and
        splitting the string by commas.
        """

        # If no fields, get it in param
        if not explode_infos_fields:
            explode_infos_fields = self.get_param().get("explode_infos_fields", None)

        # If no fields, defined as all fields in header using keyword
        if not explode_infos_fields:
            explode_infos_fields = '*'

        # If fields list not empty
        if explode_infos_fields:

            # Input fields list
            if isinstance(explode_infos_fields, str):
                fields_input = explode_infos_fields.split(",")
            elif isinstance(explode_infos_fields, list):
                fields_input = explode_infos_fields
            else:
                fields_input = []

            # Fields list without * keyword
            fields_without_all = fields_input.copy()
            if "*".casefold() in (item.casefold() for item in fields_without_all):
                fields_without_all.remove("*")

            # Fields in header
            fields_in_header = sorted(list(set(self.get_header().infos)))

            # Construct list of fields
            fields_output = []
            for field in fields_input:

                # Strip field
                field = field.strip()

                # format keyword * in regex
                if field.upper() in ["*"]:
                    field = '.*'

                # Find all fields with pattern
                r = re.compile(field)
                fields_search = sorted(list(filter(r.match, fields_in_header)))
                
                # Remove fields input from search
                if fields_search != [field]:
                    fields_search = sorted(list(set(fields_search).difference(fields_input)))

                # If field is not in header (avoid not well formatted header)
                if not fields_search and not remove_fields_not_in_header:
                    fields_search = [field]
                
                # Add found fields
                for new_field in fields_search:
                    # Add field, if not already exists, and if it is in header (if asked)
                    if new_field not in fields_output and (not remove_fields_not_in_header or new_field in fields_in_header) and new_field not in [".*"]:
                        fields_output.append(new_field)
                        
            return fields_output

        else:

            return []
        

    def get_explode_infos_prefix(self, explode_infos_prefix:str = None) -> str:
        """
        The function `get_explode_infos_prefix` returns the value of the `explode_infos_prefix` parameter, or
        the value of `self.get_param().get("explode_infos_prefix", None)` if `explode_infos_prefix` is
        not provided.
        
        :param explode_infos_prefix: The parameter `explode_infos_prefix` is a string that specifies a
        prefix to be used for exploding or expanding information
        :type explode_infos_prefix: str
        :return: the value of the variable `explode_infos_prefix`.
        """

        if not explode_infos_prefix:
            explode_infos_prefix = self.get_param().get("explode_infos_prefix", "")

        return explode_infos_prefix


    def add_column(self, table_name, column_name, column_type, default_value=None, drop:bool = False) -> dict:
        """
        The `add_column` function adds a column to a SQLite or DuckDB table with a default value if it
        doesn't already exist.
        
        :param table_name: The name of the table to which you want to add a column
        :param column_name: The parameter "column_name" is the name of the column that you want to add
        to the table
        :param column_type: The `column_type` parameter specifies the data type of the column that you
        want to add to the table. It should be a string that represents the desired data type, such as
        "INTEGER", "TEXT", "REAL", etc
        :param default_value: The `default_value` parameter is an optional parameter that specifies the
        default value for the newly added column. If a default value is provided, it will be assigned to
        the column for any existing rows that do not have a value for that column
        :param drop: The `drop` parameter is a boolean flag that determines whether to drop the column
        if it already exists in the table. If `drop` is set to `True`, the function will drop the
        existing column before adding the new column. If `drop` is set to `False` (default),, defaults
        to False
        :type drop: bool (optional)
        :return: a boolean value indicating whether the column was successfully added to the table.
        """
        
        # added
        added = False
        dropped = False

        # Check if the column already exists in the table
        query = f""" SELECT * FROM {table_name} LIMIT 0 """
        columns = self.get_query_to_df(query).columns.tolist()
        if column_name in columns:
            log.debug(f"The {column_name} column already exists in the {table_name} table")
            if drop:
                self.drop_column(table_name=table_name, column_name=column_name)
                dropped = True
            else:
                return None
        else:
            log.debug(f"The {column_name} column NOT exists in the {table_name} table")

        # Add column in table
        add_column_query = f""" ALTER TABLE {table_name} ADD COLUMN "{column_name}" {column_type} """
        if default_value is not None:
            add_column_query += f" DEFAULT {default_value}"
        self.execute_query(add_column_query)
        added = not dropped
        log.debug(f"The {column_name} column was successfully added to the {table_name} table")
        
        if added:
            added_column = {
                "table_name": table_name,
                "column_name": column_name,
                "column_type": column_type,
                "default_value": default_value
            }
        else:
            added_column = None

        return added_column


    def drop_column(self, column:dict = None, table_name:str = None, column_name:str = None) -> bool:
        """
        The `drop_column` function drops a specified column from a given table in a database and returns
        True if the column was successfully dropped, and False if the column does not exist in the
        table.
        
        :param column: The `column` parameter is a dictionary that contains information about the column
        you want to drop. It has two keys:
        :type column: dict
        :param table_name: The `table_name` parameter is the name of the table from which you want to
        drop a column
        :type table_name: str
        :param column_name: The `column_name` parameter is the name of the column that you want to drop
        from the table
        :type column_name: str
        :return: a boolean value. It returns True if the column was successfully dropped from the table,
        and False if the column does not exist in the table.
        """

        # Find column infos
        if column:
            table_name = column.get("table_name", None)
            column_name = column.get("column_name", None)

        if not table_name and not column_name:
            return False

        # Removed
        removed = False

        # Check if the column already exists in the table
        query = f""" SELECT * FROM {table_name} LIMIT 0 """
        columns = self.get_query_to_df(query).columns.tolist()
        if column_name in columns:
            log.debug(f"The {column_name} column exists in the {table_name} table")
        else:
            log.debug(f"The {column_name} column NOT exists in the {table_name} table")
            return False
        
        # Add column in table # ALTER TABLE integers DROP k
        add_column_query = f""" ALTER TABLE {table_name} DROP "{column_name}" """
        self.execute_query(add_column_query)
        removed = True
        log.debug(f"The {column_name} column was successfully dropped to the {table_name} table")
        
        return removed
    

    def explode_infos(self, prefix: str = None, create_index: bool = False, fields: list = None, force:bool = False, proccess_all_fields_together:bool = False) -> list:
        """
        The `explode_infos` function takes a VCF file and explodes the INFO fields into individual
        columns, returning a list of added columns.
        
        :param prefix: The `prefix` parameter is a string that is used as a prefix for the exploded INFO
        fields. If the `prefix` is not provided or is set to `None`, the function will use the value of
        `self.get_explode_infos_prefix()` as the prefix
        :type prefix: str
        :param create_index: The `create_index` parameter is a boolean flag that specifies whether to
        create indexes on the exploded INFO fields. If set to `True`, indexes will be created; if set to
        `False`, indexes will not be created. The default value is `False`, defaults to False
        :type create_index: bool (optional)
        :param fields: The `fields` parameter is a list of INFO fields that you want to explode into
        individual columns. If this parameter is not provided, all INFO fields will be exploded
        :type fields: list
        :param force: The `force` parameter is a boolean flag that determines whether to drop and
        recreate the column if it already exists in the table. If `force` is set to `True`, the column
        will be dropped and recreated. If `force` is set to `False`, the column will not be dropped,
        defaults to False
        :type force: bool (optional)
        :param proccess_all_fields_together: The `proccess_all_fields_together` parameter is a boolean
        flag that determines whether to process all the INFO fields together or individually. If set to
        `True`, all the INFO fields will be processed together. If set to `False`, each INFO field will
        be processed individually, defaults to False
        :type proccess_all_fields_together: bool (optional)
        :return: The function `explode_infos` returns a list of added columns.
        """
        
        # drop indexes
        self.drop_indexes()

        # connexion format 
        connexion_format = self.get_connexion_format()

        # Access
        access = self.get_config().get("access", None)

        # Added columns
        added_columns = []

        if access not in ["RO"]:

            # prefix
            if prefix in [None, True] or not isinstance(prefix, str):
                if self.get_explode_infos_prefix() not in [None, True]:
                    prefix = self.get_explode_infos_prefix()
                else:
                    prefix = "INFO/"

            # table variants
            table_variants = self.get_table_variants(clause="select")

            # extra infos
            try:
                extra_infos = self.get_extra_infos()
            except:
                extra_infos = []

            # Header infos
            header_infos = self.get_header().infos

            log.debug(
                f"Explode INFO fields - ADD [{len(header_infos)}] annotations fields")

            sql_info_alter_table_array = []

            # Info fields to check
            fields_list = list(header_infos)
            if fields:
                fields_list += fields
            fields_list = set(fields_list)

            # If no fields
            if not fields:
                fields = []

            # Translate fields if patterns
            fields = self.get_explode_infos_fields(explode_infos_fields=fields)

            for info in fields:

                info_id_sql = prefix+info

                if (info in fields_list or prefix+info in fields_list or info in extra_infos):

                    log.debug(
                        f"Explode INFO fields - ADD '{info}' annotations fields")

                    if info in header_infos:
                        info_type = header_infos[info].type
                        info_num = header_infos[info].num
                    else:
                        info_type = "String"
                        info_num = 0

                    type_sql = self.code_type_map_to_sql.get(info_type, "VARCHAR")
                    if info_num != 1:
                        type_sql = "VARCHAR"

                    # Add field
                    added_column = self.add_column(table_name=table_variants, column_name=info_id_sql, column_type=type_sql, default_value="null", drop=force)
                    
                    if added_column:
                        added_columns.append(added_column)

                    if added_column or force:
                        
                        # add field to index
                        self.index_additionnal_fields.append(info_id_sql)

                        # Update field array
                        if connexion_format in ["duckdb"]:
                            update_info_field = f"""
                            "{info_id_sql}" =
                                CASE
                                    WHEN REGEXP_EXTRACT(concat(';', INFO), ';{info}=([^;]*)',1) IN ('','.') THEN NULL
                                    ELSE REGEXP_EXTRACT(concat(';', INFO), ';{info}=([^;]*)',1)
                                END
                            """
                        elif connexion_format in ["sqlite"]:
                            update_info_field = f"""
                                "{info_id_sql}" =
                                    CASE
                                        WHEN instr(INFO, '{info}=') = 0 THEN NULL
                                        WHEN instr(substr(INFO, instr(INFO, '{info}=')+{len(info)+1}),';') = 0 THEN substr(substr(INFO, instr(INFO, '{info}=')+{len(info)+1}), instr(substr(INFO, instr(INFO, '{info}=')+{len(info)+1}), '=')+1)
                                        ELSE substr(substr(INFO, instr(INFO, '{info}=')+{len(info)+1}), instr(substr(INFO, instr(INFO, '{info}=')+{len(info)+1}), '=')+1, instr(substr(INFO, instr(INFO, '{info}=')+{len(info)+1}),';')-instr(substr(INFO, instr(INFO, '{info}=')+{len(info)+1}), '=')-1)
                                    END
                            """
                            
                        sql_info_alter_table_array.append(update_info_field)

            if sql_info_alter_table_array:

                # By chromosomes
                try:
                    chromosomes_list = list(self.get_query_to_df(
                        f""" SELECT "#CHROM" FROM {table_variants} GROUP BY "#CHROM" """)["#CHROM"])
                except:
                    chromosomes_list = [None]

                for chrom in chromosomes_list:
                    log.debug(
                        f"Explode INFO fields - Chromosome {chrom}...")
                    
                    # Where clause
                    where_clause = ""
                    if chrom and len(chromosomes_list) > 1:
                        where_clause = f""" WHERE "#CHROM" = '{chrom}' """

                    # Update table
                    if proccess_all_fields_together:
                        sql_info_alter_table_array_join = ", ".join(
                            sql_info_alter_table_array)
                        if sql_info_alter_table_array_join:
                            sql_info_alter_table = f"""
                                UPDATE {table_variants}
                                SET {sql_info_alter_table_array_join}
                                {where_clause}
                                """
                            log.debug(f"Explode INFO fields - Explode all {len(sql_info_alter_table_array)} fields...")
                            #log.debug(sql_info_alter_table)
                            self.conn.execute(sql_info_alter_table)
                    else:
                        sql_info_alter_num = 0
                        for sql_info_alter in sql_info_alter_table_array:
                            sql_info_alter_num += 1
                            sql_info_alter_table = f"""
                                UPDATE {table_variants}
                                SET {sql_info_alter}
                                {where_clause}
                                """
                            log.debug(f"Explode INFO fields - Explode field {sql_info_alter_num}/{len(sql_info_alter_table_array)}...")
                            #log.debug(sql_info_alter_table)
                            self.conn.execute(sql_info_alter_table)

        # create indexes
        if create_index:
            self.create_indexes()

        return added_columns


    def create_indexes(self) -> None:
        """
        Create indexes on the table after insertion
        """

        # Access
        access = self.get_config().get("access", None)

        # get table variants
        table_variants = self.get_table_variants("FROM")

        if self.get_indexing() and access not in ["RO"]:
            # Create index
            sql_create_table_index = f'CREATE INDEX IF NOT EXISTS idx_{self.get_table_variants()} ON {table_variants} ("#CHROM", "POS", "REF", "ALT")'
            self.conn.execute(sql_create_table_index)
            sql_create_table_index = f'CREATE INDEX IF NOT EXISTS idx_{self.get_table_variants()}_chrom ON {table_variants} ("#CHROM")'
            self.conn.execute(sql_create_table_index)
            sql_create_table_index = f'CREATE INDEX IF NOT EXISTS idx_{self.get_table_variants()}_pos ON {table_variants} ("POS")'
            self.conn.execute(sql_create_table_index)
            sql_create_table_index = f'CREATE INDEX IF NOT EXISTS idx_{self.get_table_variants()}_ref ON {table_variants} ( "REF")'
            self.conn.execute(sql_create_table_index)
            sql_create_table_index = f'CREATE INDEX IF NOT EXISTS idx_{self.get_table_variants()}_alt ON {table_variants} ("ALT")'
            self.conn.execute(sql_create_table_index)
            for field in self.index_additionnal_fields:
                sql_create_table_index = f""" CREATE INDEX IF NOT EXISTS "idx_{self.get_table_variants()}_{field}" ON {table_variants} ("{field}") """
                self.conn.execute(sql_create_table_index)


    def drop_indexes(self) -> None:
        """
        Create indexes on the table after insertion
        """

        # Access
        access = self.get_config().get("access", None)

        # get table variants
        table_variants = self.get_table_variants("FROM")

        # Get database format
        connexion_format = self.get_connexion_format()

        if access not in ["RO"]:
            if connexion_format in ["duckdb"]:
                sql_list_indexes = f"SELECT index_name FROM duckdb_indexes WHERE table_name='{table_variants}'"
            elif connexion_format in ["sqlite"]:
                sql_list_indexes = f"SELECT name FROM sqlite_master WHERE type='index' AND tbl_name='{table_variants}';"

            list_indexes = self.conn.execute(sql_list_indexes)
            index_names = [row[0] for row in list_indexes.fetchall()]
            for index in index_names:
                sql_drop_table_index = f""" DROP INDEX IF EXISTS "{index}" """
                self.conn.execute(sql_drop_table_index)
                

    def read_vcf_header(self, f) -> list:
        """
        It reads the header of a VCF file and returns a list of the header lines

        :param f: the file object
        :return: The header lines of the VCF file.
        """

        header_list = []
        for line in f:
            header_list.append(line)
            if line.startswith('#CHROM'):
                break
        return header_list


    def read_vcf_header_file(self, file:str = None) -> list:
        """
        The function `read_vcf_header_file` reads the header of a VCF file, either from a compressed or
        uncompressed file.
        
        :param file: The `file` parameter is a string that represents the path to the VCF header file
        that you want to read. It is an optional parameter, so if you don't provide a value, it will
        default to `None`
        :type file: str
        :param compressed: The `compressed` parameter is a boolean flag that indicates whether the VCF
        file is compressed or not. If `compressed` is set to `True`, it means that the VCF file is
        compressed using the BGZF compression format. If `compressed` is set to `False`, it means that,
        defaults to False
        :type compressed: bool (optional)
        :return: a list.
        """

        if self.get_input_compressed(input_file=file):
            with bgzf.open(file, 'rt') as f:
                return self.read_vcf_header(f=f)
        else:
            with open(file, 'rt') as f:
                return self.read_vcf_header(f=f)


    def execute_query(self, query: str):
        """
        It takes a query as an argument, executes it, and returns the results

        :param query: The query to be executed
        :return: The result of the query is being returned.
        """
        if query:
            return self.conn.execute(query)  # .fetchall()
        else:
            return None


    def export_output(self, output_file: str = None, output_header: str = None, export_header: bool = True, query: str = None, parquet_partitions:list = None, chunk_size:int = None, threads:int = None, sort:bool = False, index:bool = False, order_by:str = None) -> bool:
        """
        The `export_output` function exports data from a VCF file to a specified output file in various
        formats, including VCF, CSV, TSV, PSV, and Parquet.
        
        :param output_file: The `output_file` parameter is a string that specifies the name of the
        output file to be generated by the function. This is where the exported data will be saved
        :type output_file: str
        :param output_header: The `output_header` parameter is a string that specifies the name of the
        file where the header of the VCF file will be exported. If this parameter is not provided, the
        header will be exported to a file with the same name as the `output_file` parameter, but with
        the extension "
        :type output_header: str
        :param export_header: The `export_header` parameter is a boolean flag that determines whether
        the header of a VCF file should be exported to a separate file or not. If `export_header` is
        True, the header will be exported to a file. If `export_header` is False, the header will not
        be, defaults to True
        :type export_header: bool (optional)
        :param query: The `query` parameter is an optional SQL query that can be used to filter and
        select specific data from the VCF file before exporting it. If provided, only the data that
        matches the query will be exported
        :type query: str
        :param parquet_partitions: The `parquet_partitions` parameter is a list that specifies the
        columns to be used for partitioning the Parquet file during export. Partitioning is a way to
        organize data in a hierarchical directory structure based on the values of one or more columns.
        This can improve query performance when working with large datasets
        :type parquet_partitions: list
        :param chunk_size: The `chunk_size` parameter specifies the number of
        records in batch when exporting data in Parquet format. This parameter is used for
        partitioning the Parquet file into multiple files.
        :type chunk_size: int
        :param threads: The `threads` parameter is an optional parameter that specifies the number of
        threads to be used during the export process. It determines the level of parallelism and can
        improve the performance of the export operation. If not provided, the function will use the
        default number of threads
        :type threads: int
        :param sort: The `sort` parameter is a boolean flag that determines whether the output file
        should be sorted or not. If `sort` is set to `True`, the output file will be sorted based on the
        genomic coordinates of the variants. By default, the value of `sort` is `False`, defaults to
        False
        :type sort: bool (optional)
        :param index: The `index` parameter is a boolean flag that determines whether an index should be
        created on the output file. If `index` is True, an index will be created. If `index` is False,
        no index will be created. The default value is False, defaults to False
        :type index: bool (optional)
        :param order_by: The `order_by` parameter is a string that specifies the column(s) to use for
        sorting the output file. This parameter is only applicable when exporting data in VCF format
        :type order_by: str
        :return: a boolean value. It checks if the output file exists and returns True if it does, or
        None if it doesn't.
        """

        # Log
        log.info("Exporting...")

        # If no output, get it
        if not output_file:
            output_file = self.get_output()

        # If not threads
        if not threads:
            threads = self.get_threads()

        # Auto header name with extension
        if export_header or output_header:
            if not output_header:
                output_header = f"{output_file}.hdr"
            # Export header
            self.export_header(output_file=output_file)

        # Parquet partition
        if not parquet_partitions:
            parquet_partitions = self.get_param().get("parquet_partitions", None)

        # Chunk size
        if not chunk_size:
            chunk_size = self.get_config().get("chunk_size", None)

        # Order by
        if not order_by:
            order_by = self.get_param().get("order_by", "")

        # Header in output
        header_in_output = self.get_param().get("header_in_output", False)

        # Database
        database_source=self.get_connexion()

        # Connexion format
        connexion_format = self.get_connexion_format()
        
        # Explode infos
        if self.get_explode_infos():
            self.explode_infos(prefix=self.get_explode_infos_prefix(), fields=self.get_explode_infos_fields(), force=False)
        
        # Tmp files to remove
        tmp_to_remove = []

        #if connexion_format in ["sqlite"] or query:
        if connexion_format in ["sqlite"]:

            # Export in Parquet
            random_tmp = ''.join(random.choice(string.ascii_lowercase) for i in range(10))
            database_source = f"""{output_file}.{random_tmp}.database_export.parquet"""
            tmp_to_remove.append(database_source)

            # Table Variants
            table_variants = self.get_table_variants()

            # Create export query
            # if query:
            #     sql_query_export_subquery = f"""
            #         SELECT * FROM ({query})
            #         """
            #if connexion_format in ["sqlite"]:
            sql_query_export_subquery = f"""
                SELECT * FROM {table_variants}
                """

            # Write source file
            fp.write(database_source, self.get_query_to_df(sql_query_export_subquery))

        # Create database
        database = Database(database=database_source, table="variants", header_file=output_header, conn_config=self.get_connexion_config())
        
        # Existing colomns header
        #existing_columns_header = database.get_header_file_columns(output_header)
        existing_columns_header = database.get_header_columns_from_database()

        # Export file
        database.export(output_database=output_file, output_header=output_header, existing_columns_header=existing_columns_header, parquet_partitions=parquet_partitions, chunk_size=chunk_size, threads=threads, sort=sort, index=index, header_in_output=header_in_output, order_by=order_by, query=query)
        
        # Remove
        remove_if_exists(tmp_to_remove)

        return (os.path.exists(output_file) or None) and (os.path.exists(output_file) or None)


    def get_extra_infos(self, table: str = None) -> list:
        """
        > This function returns a list of columns that are in the table but not in the header

        The function is called `get_extra_infos` and it takes two arguments: `self` and `table`. The
        `self` argument is a reference to the object that called the function. The `table` argument is
        the name of the table that we want to get the extra columns from

        :param table: The table to get the extra columns from. If not specified, it will use the
        variants table
        :param format: The format of the output. If it's "sql", it will return a string of the extra
        columns separated by commas. If it's "list", it will return a list of the extra columns
        :return: A list of columns that are in the table but not in the header
        """

        header_columns = []

        if not table:
            table = self.get_table_variants(clause="from")
            header_columns = self.get_header_columns()

        # Check all columns in the database
        query = f""" SELECT * FROM {table} LIMIT 1 """
        log.debug(f"query {query}")
        table_columns = self.get_query_to_df(query).columns.tolist()
        extra_columns = []

        # Construct extra infos (not in header)
        for column in table_columns:
            if column not in header_columns:
                extra_columns.append(column)

        return extra_columns


    def get_extra_infos_sql(self, table: str = None) -> str:
        """
        It returns a string of the extra infos, separated by commas, and each extra info is surrounded
        by double quotes

        :param table: The name of the table to get the extra infos from. If None, the default table is
        used
        :type table: str
        :return: A string of the extra infos
        """

        return ", ".join(['"' + str(elem) + '"' for elem in self.get_extra_infos(table=table)])


    def export_header(self, header_name: str = None, output_file: str = None) -> str:
        """
        It takes a VCF file, and writes the header to a new file

        :param header_name: the name of the header file to be created. If not specified, the header will
        be written to the output file
        :return: The name of the temporary header file.
        """

        if not header_name and not output_file:
            output_file = self.get_output()

        if self.get_header():

            # Get header object
            haeder_obj = self.get_header()

            # Create database
            db_for_header = Database(database=self.get_input())
            
            # Get real columns in the file
            db_header_columns = db_for_header.get_columns()

            with tempfile.TemporaryDirectory() as tmpdir:

                # Write header file
                header_file_tmp = os.path.join(tmpdir,'header')
                f = open(header_file_tmp, 'w')
                vcf.Writer(f, haeder_obj)
                f.close()
                
                # Replace #CHROM line with rel columns
                header_list = db_for_header.read_header_file(header_file=header_file_tmp)
                header_list[-1] = "\t".join(db_header_columns)

            tmp_header_name = output_file + ".hdr"

            f = open(tmp_header_name, 'w')
            for line in header_list:
                f.write(line)
            f.close()

        return tmp_header_name


    def export_variant_vcf(self, vcf_file, file_type: str = "gz", remove_info: bool = False, add_samples: bool = True, list_samples:list = [], compression: int = 1, index: bool = False, threads:int = None) -> None:
        """
        The `export_variant_vcf` function takes a VCF file and a list of samples, and returns a VCF file
        with only the samples in the list.
        
        :param vcf_file: The name of the file to write the VCF data to
        :param file_type: The `file_type` parameter specifies the type of the output file. It can be
        either "vcf" or "gz" (compressed VCF file). By default, it is set to "vcf", defaults to gz
        :type file_type: str (optional)
        :param remove_info: The `remove_info` parameter is a boolean flag that determines whether to
        remove the INFO field from the output VCF file. If set to `True`, the INFO field will be
        removed. If set to `False`, the INFO field will be included in the output file. If you want to
        remove, defaults to False
        :type remove_info: bool (optional)
        :param add_samples: A boolean parameter that determines whether the samples should be added to
        the VCF file or not. If set to True, the samples will be added. If set to False, the samples
        will be removed. The default value is True, defaults to True
        :type add_samples: bool (optional)
        :param list_samples: The `list_samples` parameter is a list of samples that you want to include
        in the output VCF file. By default, all samples will be included. If you provide a list of
        samples, only those samples will be included in the output file
        :type list_samples: list
        :param compression: The `compression` parameter determines the level of compression for the
        output VCF file. It ranges from 1 to 9, with 1 being the fastest and 9 being the most
        compressed. The default value is 1, defaults to 1
        :type compression: int (optional)
        :param index: The `index` parameter is a boolean flag that determines whether or not to create
        an index for the output VCF file. If `index` is set to `True`, the output VCF file will be
        indexed using tabix. If `index` is set to `False`, no index will, defaults to False
        :type index: bool (optional)
        :param threads: The `threads` parameter specifies the number of threads to use for exporting the
        VCF file. It is an optional parameter, so if it is not provided, the code will use the value
        returned by the `get_threads()` method
        :type threads: int
        """

        # Extract VCF
        log.debug("Export VCF...")

        connexion_format = self.get_connexion_format()

        table_variants = self.get_table_variants()
        sql_query_hard = ""
        sql_query_sort = ""
        sql_query_limit = ""

        # Info fields
        if remove_info:
            if type(remove_info) != str:
                remove_info = "."
            info_field = f"""'{remove_info}' as INFO"""
        else:
            info_field = "INFO"
        # samples fields
        if add_samples or list_samples:
            if not list_samples:
                samples_fields = " , FORMAT , " + \
                    " , ".join(self.get_header_sample_list())
            else:
                samples_fields = " , FORMAT , " + \
                    " , ".join(list_samples)
            log.debug(f"samples_fields: {samples_fields}")
        else:
            samples_fields = ""

        # Header (without "#CHROM")
        tmp_header = NamedTemporaryFile(
            prefix=self.get_prefix(), dir=self.get_tmp_dir(), delete=False)
        tmp_header_name = tmp_header.name
        with open(tmp_header_name, 'w') as header_f:
            for head in self.get_header("list"):
                if not head.startswith("#CHROM"):
                    # Clean head for malformed header
                    head_clean = head
                    head_clean = re.subn('##FORMAT=<ID=(.*),Number=(.*),Type=Flag', r'##FORMAT=<ID=\1,Number=\2,Type=String', head_clean,2)[0]
                    # Write header
                    header_f.write(head_clean)

        # Variants
        tmp_variants = NamedTemporaryFile(
            prefix=self.get_prefix(), dir=self.get_tmp_dir(), suffix="", delete=False)
        tmp_variants_name = tmp_variants.name
        select_fields = """ "#CHROM", POS, ID, REF, ALT, QUAL, FILTER """

        sql_query_select = f""" SELECT {select_fields}, {info_field} {samples_fields} FROM {table_variants} WHERE 1 {sql_query_hard} {sql_query_sort} {sql_query_limit} """

        if connexion_format in ["duckdb"]:
            sql_query_export = f"COPY ({sql_query_select}) TO '{tmp_variants_name}' WITH (FORMAT CSV, DELIMITER '\t', HEADER, QUOTE '', COMPRESSION 'gzip')"
            self.conn.execute(sql_query_export)
        elif connexion_format in ["sqlite"]:
            cursor = pd.read_sql(sql_query_select, self.conn)
            cursor.to_csv(tmp_variants_name, sep='\t', compression='gzip', quoting='', index=False)

        # Threads
        if not threads:
            threads = self.get_threads()

        # export format
        if file_type in ["vcf"]:
            compression_type = "none"
        else:
            compression_type = "bgzip"

        concat_and_compress_files(input_files=[tmp_header_name, tmp_variants_name], output_file=vcf_file, compression_type=compression_type, threads=threads, sort=True, index=index, compression_level=compression)


    def run_commands(self, commands: list = [], threads: int = 1) -> None:
        """
        It takes a list of commands and runs them in parallel using the number of threads specified

        :param commands: A list of commands to run
        :param threads: The number of threads to use, defaults to 1 (optional)
        """

        run_parallel_commands(commands, threads)


    def get_threads(self) -> int:
        """
        It returns the number of threads to use for the current job
        :return: The number of threads.
        """

        input_thread = self.config.get("threads", None)
        if not input_thread:
            threads = 1
        elif int(input_thread) <= 0:
            threads = os.cpu_count()
        else:
            threads = int(input_thread)
        return threads


    def update_from_vcf(self, vcf_file: str) -> None:
        """
        > If the database is duckdb, then use the parquet method, otherwise use the sqlite method

        :param vcf_file: the path to the VCF file
        """

        connexion_format = self.get_connexion_format()

        if connexion_format in ["duckdb"]:
            self.update_from_vcf_duckdb(vcf_file)
        elif connexion_format in ["sqlite"]:
            self.update_from_vcf_sqlite(vcf_file)


    def update_from_vcf_duckdb(self, vcf_file: str) -> None:
        """
        It takes a VCF file and updates the INFO column of the variants table in the database with the
        INFO column of the VCF file

        :param vcf_file: the path to the VCF file
        """

        # varaints table
        table_variants = self.get_table_variants()

        # Loading VCF into temporaire table
        skip = self.get_header_length(file=vcf_file)
        vcf_df = pd.read_csv(vcf_file, sep='\t', engine='c', skiprows=skip, header=0, low_memory=False, )
        sql_query_update = f"""
        UPDATE {table_variants} as table_variants
            SET INFO = concat(
                            CASE
                                WHEN INFO NOT IN ('', '.')
                                THEN INFO
                                ELSE ''
                            END,
                            (
                                SELECT 
                                    concat(
                                        CASE
                                            WHEN table_variants.INFO NOT IN ('','.') AND table_parquet.INFO NOT IN ('','.')
                                            THEN ';'
                                            ELSE ''
                                        END
                                        ,
                                        CASE
                                            WHEN table_parquet.INFO NOT IN ('','.')
                                            THEN table_parquet.INFO
                                            ELSE ''
                                        END
                                    )
                                FROM vcf_df as table_parquet
                                        WHERE CAST(table_parquet.\"#CHROM\" AS VARCHAR) = CAST(table_variants.\"#CHROM\" AS VARCHAR)
                                        AND table_parquet.\"POS\" = table_variants.\"POS\"
                                        AND table_parquet.\"ALT\" = table_variants.\"ALT\"
                                        AND table_parquet.\"REF\" = table_variants.\"REF\"
                            )
                        )
            ;
            """
        self.conn.execute(sql_query_update)


    def update_from_vcf_sqlite(self, vcf_file: str) -> None:
        """
        It creates a temporary table in the SQLite database, loads the VCF file into the temporary
        table, then updates the INFO column of the variants table with the INFO column of the temporary
        table

        :param vcf_file: The path to the VCF file you want to update the database with
        """

        # Create a temporary table for the VCF
        table_vcf = 'tmp_vcf'
        sql_create = f"CREATE TEMPORARY TABLE {table_vcf} AS SELECT * FROM variants WHERE 0"
        self.conn.execute(sql_create)

        # Loading VCF into temporaire table
        vcf_df = pd.read_csv(vcf_file, sep='\t', comment='#',
                             header=None, low_memory=False)
        vcf_df.columns = ['#CHROM', 'POS', 'ID',
                          'REF', 'ALT', 'QUAL', 'FILTER', 'INFO']
        vcf_df.to_sql(table_vcf, self.conn, if_exists='append', index=False)

        # Update table 'variants' with VCF data
        # warning: CONCAT as || operator
        sql_query_update = f"""
            UPDATE variants as table_variants
            SET INFO = CASE
                            WHEN INFO NOT IN ('', '.')
                            THEN INFO
                            ELSE ''
                        END ||
                        (
                        SELECT 
                            CASE 
                                WHEN table_variants.INFO NOT IN ('','.') 
                                    AND table_vcf.INFO NOT IN ('','.')  
                                THEN ';' 
                                ELSE '' 
                            END || 
                            CASE 
                                WHEN table_vcf.INFO NOT IN ('','.') 
                                THEN table_vcf.INFO 
                                ELSE '' 
                            END
                        FROM {table_vcf} as table_vcf
                        WHERE table_vcf.\"#CHROM\" = table_variants.\"#CHROM\"
                            AND table_vcf.\"POS\" = table_variants.\"POS\"
                            AND table_vcf.\"ALT\" = table_variants.\"ALT\"
                            AND table_vcf.\"REF\" = table_variants.\"REF\"
                        )
        """
        self.conn.execute(sql_query_update)

        # Drop temporary table
        sql_drop = f"DROP TABLE {table_vcf}"
        self.conn.execute(sql_drop)


    def drop_variants_table(self) -> None:
        """
        > This function drops the variants table
        """

        table_variants = self.get_table_variants()
        sql_table_variants = f"DROP TABLE IF EXISTS {table_variants}"
        self.conn.execute(sql_table_variants)


    def set_variant_id(self, variant_id_column:str = "variant_id", force:bool = None) -> str:
        """
        It adds a column to the variants table called `variant_id` and populates it with a hash of the
        `#CHROM`, `POS`, `REF`, and `ALT` columns
        
        :param variant_id_column: The name of the column to be created in the variants table, defaults
        to variant_id
        :type variant_id_column: str (optional)
        :param force: If True, the variant_id column will be created even if it already exists
        :type force: bool
        :return: The name of the column that contains the variant_id
        """

        # Assembly
        assembly = self.get_param().get("assembly", self.get_config().get("assembly", DEFAULT_ASSEMBLY))

        # INFO/Tag prefix
        prefix = self.get_explode_infos_prefix()

        # Explode INFO/SVTYPE
        self.explode_infos(prefix=prefix, fields=["SVTYPE"])

        # variants table
        table_variants = self.get_table_variants()

        # variant_id column
        if not variant_id_column:
            variant_id_column = "variant_id"

        # Creta variant_id column
        if "variant_id" not in self.get_extra_infos() or force:

            # Create column
            self.add_column(table_name=table_variants, column_name=variant_id_column, column_type="UBIGINT", default_value="0")

            # Update column
            self.conn.execute(
                f"""
                    UPDATE {table_variants}
                    SET "{variant_id_column}" = hash('{assembly}', "#CHROM", "POS", "REF", "ALT", '"{prefix}SVTYPE"')
                """)
        
        # return variant_id column name
        return variant_id_column


    def get_variant_id_column(self, variant_id_column:str = "variant_id", force:bool = None) -> str:
        """
        This function returns the variant_id column name
        
        :param variant_id_column: The name of the column in the dataframe that contains the variant IDs,
        defaults to variant_id
        :type variant_id_column: str (optional)
        :param force: If True, will force the variant_id to be set to the value of variant_id_column. If
        False, will only set the variant_id if it is not already set. If None, will set the variant_id
        if it is not already set, or if it is set
        :type force: bool
        :return: The variant_id column name.
        """
        return self.set_variant_id(variant_id_column=variant_id_column, force=force)
    

    ###
    # Annotation
    ###

    def annotation(self) -> None:
        """
        It annotates the VCF file with the annotations specified in the config file.
        """

        # Config
        config = self.get_config()

        # Param
        param = self.get_param()

        # Param - Assembly
        assembly = param.get("assembly", config.get("assembly", None))
        if not assembly:
            assembly = DEFAULT_ASSEMBLY
            log.warning(f"Default assembly '{assembly}'")

        # annotations databases folders
        annotations_databases = set(
            config.get("folders", {}).get("databases", {}).get("annotations", [DEFAULT_ANNOTATIONS_FOLDER])
            + config.get("folders", {}).get("databases", {}).get("parquet", ["/databases/parquet/current"])
            + config.get("folders", {}).get("databases", {}).get("bcftools", ["/databases/bcftools/current"])
        )

        if param.get("annotations"):

            # Log
            log.info("Annotations - Check annotation parameters")

            if not "annotation" in param:
                param["annotation"] = {}

            # List of annotations parameters
            annotations_list_input = param.get("annotations", {})
            annotations_list = {}

            # Explode annotations if ALL
            for annotation_file in annotations_list_input:

                if annotation_file.upper() == "ALL" or annotation_file.upper().startswith("ALL:"):

                    # check ALL parameters (formats, releases)
                    annotation_file_split = annotation_file.split(":")
                    database_formats = ["parquet"]
                    database_releases = ["current"]
                    if len(annotation_file_split) > 1:
                        database_formats = annotation_file_split[1].split("|")
                    if len(annotation_file_split) > 2:
                        database_releases = annotation_file_split[2].split("|")

                    # Scan for availabled databases
                    log.info(f"Annotations - Check annotation parameters - Scan existing databases - Assembly {[assembly]} - Formats {database_formats} - Releases {database_releases}...")
                    databases_infos_dict = databases_infos(database_folder_releases=database_releases, database_formats=database_formats, assembly=assembly, config=config)
                    log.info(f"Annotations - Check annotation parameters - Scan existing databases - {len(databases_infos_dict)} databases found")

                    # Add found databases in annotation parameters
                    for database_infos in databases_infos_dict.keys():
                        annotations_list[database_infos] = {'INFO': None}

                else:
                    annotations_list[annotation_file] = annotations_list_input[annotation_file]
            
            # Check each databases
            if len(annotations_list):

                log.info(f"Annotations - Check annotation parameters - Check {len(annotations_list)} databases...")

                for annotation_file in annotations_list:

                    # Init
                    annotations = annotations_list.get(
                        annotation_file, None)
                    
                    # Annotation snpEff
                    if annotation_file == "snpeff":
                        log.debug(f"Quick Annotation snpEff")
                        if "snpeff" not in param["annotation"]:
                            param["annotation"]["snpeff"] = {}
                        if "options" not in param["annotation"]["snpeff"]:
                            param["annotation"]["snpeff"]["options"] = ""
                    
                    # Annotation Annovar
                    elif annotation_file.startswith("annovar"):
                        log.debug(f"Quick Annotation Annovar")
                        if "annovar" not in param["annotation"]:
                            param["annotation"]["annovar"] = {}
                        if "annotations" not in param["annotation"]["annovar"]:
                            param["annotation"]["annovar"]["annotations"] = {}
                        annotation_file_split = annotation_file.split(":")
                        if len(annotation_file_split) > 1:
                            annotation_file_annotation = annotation_file_split[1]
                            param["annotation"]["annovar"]["annotations"][annotation_file_annotation] = annotations

                    # Annotation Exomiser
                    elif annotation_file.startswith("exomiser"):
                        log.debug(f"Quick Annotation Exomiser")
                        if "exomiser" not in param["annotation"]:
                            param["annotation"]["exomiser"] = {}
                        annotation_file_split = annotation_file.split(":")
                        if len(annotation_file_split) > 1:
                            annotation_file_options = annotation_file_split[1]

                            annotation_file_options_split = annotation_file_options.split('|')
                            log.debug(annotation_file_options_split)
                            
                            for annotation_file_options_split_option in annotation_file_options_split:
                                log.debug(annotation_file_options_split_option)
                                annotation_file_options_split_option_var_val = annotation_file_options_split_option.split("=")
                                annotation_file_options_split_option_var = annotation_file_options_split_option_var_val[0].strip()
                                annotation_file_options_split_option_val = annotation_file_options_split_option_var_val[1].strip()
                                log.debug(annotation_file_options_split_option_val)
                                log.debug(annotation_file_options_split_option_val)
                                if annotation_file_options_split_option_val:
                                    if not annotation_file_options_split_option_val:
                                        annotation_file_options_split_option_val = None
                                    else:
                                        annotation_file_options_split_option_val = annotation_file_options_split_option_val.replace("+", ",").replace(" ", "")
                                    param["annotation"]["exomiser"][annotation_file_options_split_option_var] = annotation_file_options_split_option_val

                            log.debug(param["annotation"]["exomiser"])

                    # Annotation Parquet or BCFTOOLS
                    else:
                        
                        # Find file
                        annotation_file_found = None

                        if os.path.exists(annotation_file):
                            annotation_file_found = annotation_file

                        else:
                            # Find within assembly folders
                            for annotations_database in annotations_databases:
                                found_files = find_all(annotation_file, os.path.join(annotations_database, assembly))
                                if len(found_files) > 0:
                                    annotation_file_found = found_files[0]
                                    break
                            if not annotation_file_found and not assembly:
                                # Find within folders
                                for annotations_database in annotations_databases:
                                    found_files = find_all(annotation_file, annotations_database)
                                    if len(found_files) > 0:
                                        annotation_file_found = found_files[0]
                                        break
                        log.debug(f"for {annotation_file} annotation_file_found={annotation_file_found}")

                        if annotation_file_found:

                            database = Database(database=annotation_file_found)
                            quick_annotation_format = database.get_format()
                            quick_annotation_is_compressed = database.is_compressed()
                            quick_annotation_is_indexed = os.path.exists(f"{annotation_file_found}.tbi")
                            bcftools_preference = True

                            # Check Annotation Tool
                            annotation_tool = None
                            if bcftools_preference and quick_annotation_format in ["vcf", "bed"] and quick_annotation_is_compressed and quick_annotation_is_indexed:
                                annotation_tool = "bcftools"
                            elif quick_annotation_format in ["vcf", "bed", "tsv", "tsv", "csv", "json", "tbl", "parquet", "duckdb"]:
                                annotation_tool = "parquet"
                            else:
                                log.error(
                                    f"Quick Annotation File {annotation_file_found} - Format {quick_annotation_format} not supported yet")
                                raise ValueError(
                                    f"Quick Annotation File {annotation_file_found} - Format {quick_annotation_format} not supported yet"
                                )
                            
                            log.debug(f"Quick Annotation File {annotation_file} - Annotation tool: {annotation_tool}")
                            
                            # Annotation Tool dispatch
                            if annotation_tool:
                                if annotation_tool not in param["annotation"]:
                                    param["annotation"][annotation_tool] = {}
                                if "annotations" not in param["annotation"][annotation_tool]:
                                    param["annotation"][annotation_tool]["annotations"] = {}
                                param["annotation"][annotation_tool]["annotations"][annotation_file_found] = annotations

                        else:
                            log.error(f"Quick Annotation File {annotation_file} does NOT exist")

                self.set_param(param)

        if param.get("annotation", None):
            log.info("Annotations")
            if param.get("annotation", {}).get("parquet", None):
                log.info("Annotations 'parquet'...")
                self.annotation_parquet()
            if param.get("annotation", {}).get("bcftools", None):
                log.info("Annotations 'bcftools'...")
                self.annotation_bcftools()
            if param.get("annotation", {}).get("annovar", None):
                log.info("Annotations 'annovar'...")
                self.annotation_annovar()
            if param.get("annotation", {}).get("snpeff", None):
                log.info("Annotations 'snpeff'...")
                self.annotation_snpeff()
            if param.get("annotation", {}).get("exomiser", None):
                log.info("Annotations 'exomiser'...")
                self.annotation_exomiser()
            if param.get("annotation", {}).get("varank", None):
                log.info("Annotations 'varank'...")

        # Explode INFOS fields into table fields
        if self.get_explode_infos():
            self.explode_infos(prefix=self.get_explode_infos_prefix(), fields=self.get_explode_infos_fields(), force=True)


    def annotation_bcftools(self, threads: int = None) -> None:
        """
        This function annotate with bcftools

        :param threads: Number of threads to use
        :return: the value of the variable "return_value".
        """

        # DEBUG
        log.debug("Start annotation with bcftools databases")

        # Threads
        if not threads:
            threads = self.get_threads()
        log.debug("Threads: "+str(threads))

        # DEBUG
        delete_tmp = True
        if self.get_config().get("verbosity", "warning") in ["debug"]:
            delete_tmp = False
            log.debug("Delete tmp files/folders: "+str(delete_tmp))

        # Config
        databases_folders = set(
            self.get_config().get("folders", {}).get("databases", {}).get("annotations", ["."])
            + self.get_config().get("folders", {}).get("databases", {}).get("bcftools", ["."])
        )
        log.debug("Databases annotations: " + str(databases_folders))

        # Param
        annotations = self.get_param().get("annotation", {}).get(
            "bcftools", {}).get("annotations", None)
        log.debug("Annotations: " + str(annotations))

        # Assembly
        assembly = self.get_param().get("assembly", self.get_config().get("assembly", DEFAULT_ASSEMBLY))

        # Data
        table_variants = self.get_table_variants()

        # Check if not empty
        log.debug("Check if not empty")
        sql_query_chromosomes = f"""SELECT count(*) as count FROM {table_variants} as table_variants"""
        sql_query_chromosomes_df = self.get_query_to_df(sql_query_chromosomes)
        if not sql_query_chromosomes_df["count"][0]:
            log.info(f"VCF empty")
            return

        # Export in VCF
        log.debug("Create initial file to annotate")
        tmp_vcf = NamedTemporaryFile(prefix=self.get_prefix(
        ), dir=self.get_tmp_dir(), suffix=".vcf.gz", delete=False)
        tmp_vcf_name = tmp_vcf.name

        # VCF header
        vcf_reader = self.get_header()
        log.debug("Initial header: " + str(vcf_reader.infos))

        # Existing annotations
        for vcf_annotation in self.get_header().infos:

            vcf_annotation_line = self.get_header().infos.get(vcf_annotation)
            log.debug(
                f"Existing annotations in VCF: {vcf_annotation} [{vcf_annotation_line}]")

        if annotations:

            tmp_ann_vcf_list = []
            commands = []
            tmp_files = []
            err_files = []

            for annotation in annotations:
                annotation_fields = annotations[annotation]

                if not annotation_fields:
                    annotation_fields = {"INFO": None}

                log.debug(f"Annotation '{annotation}'")
                log.debug(
                    f"Annotation '{annotation}' - fields: {annotation_fields}")

                # Create Database
                database = Database(database=annotation, databases_folders=databases_folders, assembly=assembly)

                # Find files
                db_file = database.get_database()
                db_hdr_file = database.get_header_file()
                db_file_type = database.get_format()
                db_tbi_file = f"{db_file}.tbi"
                db_file_compressed = database.is_compressed()

                # Check if compressed
                if not db_file_compressed:
                    log.error(f"Annotation '{annotation}' - {db_file} NOT compressed file")
                    raise ValueError(f"Annotation '{annotation}' - {db_file} NOT compressed file")

                # Check if indexed
                if not os.path.exists(db_tbi_file):
                    log.error(f"Annotation '{annotation}' - {db_file} NOT indexed file")
                    raise ValueError(f"Annotation '{annotation}' - {db_file} NOT indexed file")

                # Check index - try to create if not exists
                if not os.path.exists(db_file) or not os.path.exists(db_hdr_file):
                    log.error("Annotation failed: database not valid")
                    log.error(f"Annotation annotation file: {db_file}")
                    log.error(f"Annotation annotation header: {db_hdr_file}")
                    log.error(f"Annotation annotation index: {db_tbi_file}")
                    raise ValueError(
                        f"Annotation failed: database not valid - annotation file {db_file} / annotation header {db_hdr_file} / annotation index {db_tbi_file} / annotation compression {db_file_compressed}")
                else:

                    log.debug(
                        f"Annotation '{annotation}' - file: " + str(db_file) + " and " + str(db_hdr_file))

                    # Load header as VCF object
                    db_hdr_vcf = Variants(input=db_hdr_file)
                    db_hdr_vcf_header_infos = db_hdr_vcf.get_header().infos
                    log.debug("Annotation database header: " +
                              str(db_hdr_vcf_header_infos))

                    # For all fields in database
                    if "ALL" in annotation_fields or "INFO" in annotation_fields:
                        annotation_fields = {
                            key: key for key in db_hdr_vcf_header_infos}
                        log.debug(
                            "Annotation database header - All annotations added: " + str(annotation_fields))

                    # Create file for field rename
                    log.debug("Create file for field rename")
                    tmp_rename = NamedTemporaryFile(prefix=self.get_prefix(
                    ), dir=self.get_tmp_dir(), suffix=".rename", delete=False)
                    tmp_rename_name = tmp_rename.name
                    tmp_files.append(tmp_rename_name)

                    # Number of fields
                    nb_annotation_field = 0
                    annotation_list = []

                    for annotation_field in annotation_fields:

                        # field new name, if parametered SKIPPED !!!!!! not managed actually TODO
                        annotation_fields_new_name = annotation_fields.get(
                            annotation_field, annotation_field)
                        if not annotation_fields_new_name:
                            annotation_fields_new_name = annotation_field

                        # Check if field is in DB and if field is not elready in input data
                        if annotation_field in db_hdr_vcf.get_header().infos and annotation_fields_new_name not in self.get_header().infos:

                            log.info(
                                f"Annotation '{annotation}' - '{annotation_field}' -> 'INFO/{annotation_fields_new_name}'")

                            # Add INFO field to header
                            db_hdr_vcf_header_infos_number = db_hdr_vcf_header_infos[
                                annotation_field].num or "."
                            db_hdr_vcf_header_infos_type = db_hdr_vcf_header_infos[
                                annotation_field].type or "String"
                            db_hdr_vcf_header_infos_description = db_hdr_vcf_header_infos[
                                annotation_field].desc or f"{annotation_field} description"
                            db_hdr_vcf_header_infos_source = db_hdr_vcf_header_infos[
                                annotation_field].source or "unknown"
                            db_hdr_vcf_header_infos_version = db_hdr_vcf_header_infos[
                                annotation_field].version or "unknown"

                            vcf_reader.infos[annotation_fields_new_name] = vcf.parser._Info(
                                annotation_fields_new_name,
                                db_hdr_vcf_header_infos_number,
                                db_hdr_vcf_header_infos_type,
                                db_hdr_vcf_header_infos_description,
                                db_hdr_vcf_header_infos_source,
                                db_hdr_vcf_header_infos_version,
                                self.code_type_map[db_hdr_vcf_header_infos_type]
                            )

                            annotation_list.append(annotation_field)

                            nb_annotation_field += 1

                        else:

                            if annotation_field not in db_hdr_vcf.get_header().infos:
                                log.warning(
                                    f"Annotation '{annotation}' - '{annotation_field}' - not available in vcf/bed file")
                            if annotation_fields_new_name in self.get_header().infos:
                                log.warning(
                                    f"Annotation '{annotation}' - '{annotation_fields_new_name}' - already exists (skipped)")

                    log.info(
                        f"Annotation '{annotation}' - {nb_annotation_field} annotations available in vcf/bed file")

                    annotation_infos = ",".join(annotation_list)

                    if annotation_infos != "":

                        # Protect header for bcftools (remove "#CHROM" and variants line)
                        log.debug(
                            "Protect Header file - remove #CHROM line if exists")
                        tmp_header_vcf = NamedTemporaryFile(prefix=self.get_prefix(
                        ), dir=self.get_tmp_dir(), suffix=".hdr", delete=False)
                        tmp_header_vcf_name = tmp_header_vcf.name
                        tmp_files.append(tmp_header_vcf_name)
                        # Command
                        if db_hdr_file.endswith(".gz"):
                            command_extract_header = f"zcat {db_hdr_file} | grep '^##' > {tmp_header_vcf_name}"
                        else:
                            command_extract_header = f"cat {db_hdr_file} | grep '^##' > {tmp_header_vcf_name}"
                        # Run
                        run_parallel_commands([command_extract_header], 1)

                        # Find chomosomes
                        log.debug("Find chromosomes ")
                        sql_query_chromosomes = f"""SELECT table_variants.\"#CHROM\" as CHROM FROM {table_variants} as table_variants GROUP BY table_variants.\"#CHROM\""""
                        sql_query_chromosomes_df = self.get_query_to_df(
                            sql_query_chromosomes)
                        chomosomes_list = list(
                            sql_query_chromosomes_df["CHROM"])

                        log.debug("Chromosomes found: " +
                                  str(list(chomosomes_list)))

                        # Add rename info
                        run_parallel_commands(
                            [f"echo 'INFO/{annotation_field} {annotation_fields_new_name}' >> {tmp_rename_name}"], 1)

                        # BED columns in the annotation file
                        if db_file_type in ["bed"]:
                            annotation_infos = "CHROM,POS,POS," + annotation_infos

                        for chrom in chomosomes_list:

                            # Create BED on initial VCF
                            log.debug(
                                "Create BED on initial VCF: " + str(tmp_vcf_name))
                            tmp_bed = NamedTemporaryFile(prefix=self.get_prefix(
                            ), dir=self.get_tmp_dir(), suffix=".bed", delete=False)
                            tmp_bed_name = tmp_bed.name
                            tmp_files.append(tmp_bed_name)

                            # Detecte regions
                            log.debug(
                                f"Annotation '{annotation}' - Chromosome '{chrom}' - Start detecting regions...")
                            window = 1000000
                            sql_query_intervals_for_bed = f"""
                                SELECT  \"#CHROM\",
                                        CASE WHEN \"POS\"-{window}-1 < 0 THEN 0 ELSE \"POS\"-{window}-1 END,
                                        \"POS\"+{window}
                                FROM {table_variants} as table_variants
                                WHERE table_variants.\"#CHROM\" = '{chrom}'
                            """
                            regions = self.conn.execute(
                                sql_query_intervals_for_bed).fetchall()
                            merged_regions = merge_regions(regions)
                            log.debug(
                                f"Annotation '{annotation}' - Chromosome '{chrom}' - Stop detecting regions...")

                            header = ["#CHROM", "START", "END"]
                            with open(tmp_bed_name, "w") as f:
                                # Write the header with tab delimiter
                                f.write("\t".join(header) + "\n")
                                for d in merged_regions:
                                    # Write each data row with tab delimiter
                                    f.write("\t".join(map(str, d)) + "\n")

                            # Tmp files
                            tmp_annotation_vcf = NamedTemporaryFile(prefix=self.get_prefix(
                            ), dir=self.get_tmp_dir(), suffix=".vcf.gz", delete=False)
                            tmp_annotation_vcf_name = tmp_annotation_vcf.name
                            tmp_files.append(tmp_annotation_vcf_name)
                            tmp_ann_vcf_list.append(
                                f"{tmp_annotation_vcf_name}")
                            tmp_annotation_vcf_name_err = tmp_annotation_vcf_name + ".err"
                            err_files.append(tmp_annotation_vcf_name_err)

                            # Annotate Command
                            log.debug(
                                f"Annotation '{annotation}' - add bcftools command")

                            command_annotate = f"bcftools annotate --regions-file={tmp_bed_name} -a {db_file} -h {tmp_header_vcf_name} -c {annotation_infos} --rename-annots={tmp_rename_name} {tmp_vcf_name} -o {tmp_annotation_vcf_name} -Oz 2>>{tmp_annotation_vcf_name_err} && tabix {tmp_annotation_vcf_name} 2>>{tmp_annotation_vcf_name_err} "

                            commands.append(command_annotate)

            # if some commands
            if commands:

                # Export VCF file
                self.export_variant_vcf(vcf_file=tmp_vcf_name, file_type="gz",
                                        remove_info=True, add_samples=False, compression=1, index=True)

                # Threads
                # calculate threads for annotated commands
                if commands:
                    threads_bcftools_annotate = round(threads / len(commands))
                else:
                    threads_bcftools_annotate = 1

                if not threads_bcftools_annotate:
                    threads_bcftools_annotate = 1

                # Add threads option to bcftools commands
                if threads_bcftools_annotate > 1:
                    commands_threaded = []
                    for command in commands:
                        commands_threaded.append(command.replace(
                            "bcftools annotate ", f"bcftools annotate --threads={threads_bcftools_annotate} "))
                    commands = commands_threaded

                # Command annotation multithreading
                log.debug(f"Annotation - Annotation commands: " + str(commands))
                log.info(f"Annotation - Annotation multithreaded in " +
                         str(len(commands)) + " commands")

                run_parallel_commands(commands, threads)

                # Merge
                tmp_ann_vcf_list_cmd = " ".join(tmp_ann_vcf_list)

                if tmp_ann_vcf_list_cmd:

                    # Tmp file
                    tmp_annotate_vcf = NamedTemporaryFile(prefix=self.get_prefix(
                    ), dir=self.get_tmp_dir(), suffix=".vcf.gz", delete=True)
                    tmp_annotate_vcf_name = tmp_annotate_vcf.name
                    tmp_annotate_vcf_name_err = tmp_annotate_vcf_name + ".err"
                    err_files.append(tmp_annotate_vcf_name_err)

                    # Tmp file remove command
                    tmp_files_remove_command = ""
                    if tmp_files:
                        tmp_files_remove_command = " && rm -f " + \
                            " ".join(tmp_files)

                    # Command merge
                    merge_command = f"bcftools merge --force-samples --threads={threads} {tmp_vcf_name} {tmp_ann_vcf_list_cmd} -o {tmp_annotate_vcf_name} -Oz 2>>{tmp_annotate_vcf_name_err} {tmp_files_remove_command}"
                    log.info(f"Annotation - Annotation merging " +
                             str(len(commands)) + " annotated files")
                    log.debug(f"Annotation - merge command: {merge_command}")
                    run_parallel_commands([merge_command], 1)

                    # Error messages
                    log.info(f"Error/Warning messages:")
                    error_message_command_all = []
                    error_message_command_warning = []
                    error_message_command_err = []
                    for err_file in err_files:
                        with open(err_file, 'r') as f:
                            for line in f:
                                message = line.strip()
                                error_message_command_all.append(message)
                                if line.startswith('[W::'):
                                    error_message_command_warning.append(
                                        message)
                                if line.startswith('[E::'):
                                    error_message_command_err.append(
                                        f"{err_file}: " + message)
                    # log info
                    for message in list(set(error_message_command_err + error_message_command_warning)):
                        log.info(f"   {message}")
                    # debug info
                    for message in list(set(error_message_command_all)):
                        log.debug(f"   {message}")
                    # failed
                    if len(error_message_command_err):
                        log.error("Annotation failed: Error in commands")
                        raise ValueError(
                            "Annotation failed: Error in commands")

                    # Update variants
                    log.info(f"Annotation - Updating...")
                    self.update_from_vcf(tmp_annotate_vcf_name)


    def annotation_exomiser(self, threads: int = None) -> None:
        """
        This function annotate with Exomiser

        This function uses args as parameters, in section "annotation" -> "exomiser", with sections:
        - "analysis" (dict/file):
            Full analysis dictionnary parameters (see Exomiser docs).
            Either a dict, or a file in JSON or YAML format.
            These parameters may change depending on other parameters (e.g. phenotipicFeatures/HPO)
            Default : None
        - "preset" (string):
            Analysis preset (available in config folder).
            Used if no full "analysis" is provided.
            Default: "exome"
        - "phenopacket" (dict/file):
            Samples and phenotipic features parameters (see Exomiser docs).
            Either a dict, or a file in JSON or YAML format.
            Default: None
        - "subject" (dict):
            Sample parameters (see Exomiser docs).
            Example: 
                "subject": 
                    {
                        "id": "ISDBM322017",
                        "sex": "FEMALE"
                    }
            Default: None
        - "sample" (string):
            Sample name to construct "subject" section:
                "subject": 
                    {
                        "id": "<sample>",
                        "sex": "UNKNOWN_SEX"
                    }
            Default: None
        - "phenotypicFeatures" (dict)
            Phenotypic features to construct "subject" section.
            Example: 
                "phenotypicFeatures":
                    [
                        { "type": { "id": "HP:0001159", "label": "Syndactyly" } },
                        { "type": { "id": "HP:0000486", "label": "Strabismus" } }
                    ]
        - "hpo" (list)
            List of HPO ids as phenotypic features.
            Example:
                "hpo": ['0001156', '0001363', '0011304', '0010055']
            Default: []
        - "outputOptions" (dict):
            Output options (see Exomiser docs).
            Default:
                "output_options" =
                    {
                        "outputContributingVariantsOnly": False,
                        "numGenes": 0,
                        "outputFormats": ["TSV_VARIANT", "VCF"]
                    }
        - "transcript_source" (string):
            Transcript source (either "refseq", "ucsc", "ensembl")
            Default: "refseq"
        - "exomiser_to_info" (boolean):
            Add exomiser TSV file columns as INFO fields in VCF.
            Default: False
        - "release" (string):
            Exomise database release.
            If not exists, database release will be downloaded (take a while).
            Default: None (provided by application.properties configuration file)
        - "exomiser_application_properties" (file):
            Exomiser configuration file (see Exomiser docs).
            Useful to automatically download databases (especially for specific genome databases).

        Notes:
        - If no sample in parameters, first sample in VCF will be chosen
        - If no HPO found, "hiPhivePrioritiser" analysis step will be switch off
        
        :param threads: The number of threads to use
        :return: None.
        """

        # DEBUG
        log.debug("Start annotation with Exomiser databases")

        # Threads
        if not threads:
            threads = self.get_threads()
        log.debug("Threads: "+str(threads))

        # Config
        config = self.get_config()
        log.debug("Config: " + str(config))

        # Config - Folders - Databases
        databases_folders = config.get("folders", {}).get(
            "databases", {}).get("exomiser", f"{DEFAULT_DATABASE_FOLDER}/exomiser/current")
        if not os.path.exists(databases_folders):
            log.error(f"Databases annotations: {databases_folders} NOT found")
        log.debug("Databases annotations: " + str(databases_folders))

        # Config - Java
        java_bin = config.get("tools", {}).get("java", {}).get("bin", "java")
        log.debug("Java bin: " + str(java_bin))

        # Config - Exomiser
        exomiser_jar = get_bin(bin="exomiser-cli*.jar", tool="exomiser", bin_type="jar", config=config, default_folder=f"{DEFAULT_TOOLS_FOLDER}/exomiser")
        log.debug("Exomiser bin: " + str(exomiser_jar))

        # Config - Java
        java_bin = config.get("tools", {}).get("java", {}).get("bin", "java")
        
        # Param
        param = self.get_param()
        log.debug("Param: " + str(param))

        # Param - Exomiser
        param_exomiser = param.get("annotation", {}).get("exomiser", {})
        log.debug(f"Param Exomiser: {param_exomiser}")

        # Param - Assembly
        assembly = param.get("assembly", config.get("assembly", DEFAULT_ASSEMBLY))
        log.debug("Assembly: " + str(assembly))

        # Data
        table_variants = self.get_table_variants()

        # Check if not empty
        log.debug("Check if not empty")
        sql_query_chromosomes = f"""SELECT count(*) as count FROM {table_variants} as table_variants"""
        if not self.get_query_to_df(f"{sql_query_chromosomes}")["count"][0]:
            log.info(f"VCF empty")
            return False

        # VCF header
        vcf_reader = self.get_header()
        log.debug("Initial header: " + str(vcf_reader.infos))

        # Samples
        samples = self.get_header_sample_list()
        if not samples:
            log.error("No Samples in VCF")
            return False
        log.debug(f"Samples: {samples}")
        
        # Memory limit
        if config.get("memory", None):
            memory_limit = config.get("memory", "8G")
        else:
            memory_limit = "8G"
        log.debug(f"memory_limit: {memory_limit}")

        # Exomiser java options
        exomiser_java_options = f" -Xmx{memory_limit} -XX:+UseParallelGC -XX:ParallelGCThreads={threads} "
        log.debug(f"Exomiser java options: {exomiser_java_options}")

        # Download Exomiser (if not exists)
        exomiser_release = param_exomiser.get("release", None)
        exomiser_application_properties = param_exomiser.get("exomiser_application_properties", None)
        databases_download_exomiser(assemblies=[assembly], exomiser_folder=databases_folders, exomiser_release=exomiser_release, exomiser_phenotype_release=exomiser_release, exomiser_application_properties=exomiser_application_properties)

        # Force annotation
        force_update_annotation = True

        if "Exomiser" not in self.get_header().infos or force_update_annotation:
            log.debug("Start annotation Exomiser")

            with TemporaryDirectory(dir=self.get_tmp_dir()) as tmp_dir:

                #tmp_dir = "/tmp/exomiser"

                ### ANALYSIS ###
                ################

                # Create analysis.json through analysis dict
                # either analysis in param or by default
                # depending on preset exome/genome)

                # Init analysis dict
                param_exomiser_analysis_dict = {}

                # analysis from param
                param_exomiser_analysis = param_exomiser.get("analysis", {})

                # If analysis in param -> load anlaysis json
                if param_exomiser_analysis:

                    # If param analysis is a file and exists
                    if isinstance(param_exomiser_analysis, str) and os.path.exists(param_exomiser_analysis):
                        # Load analysis file into analysis dict (either yaml or json)
                        with open(param_exomiser_analysis) as json_file:
                            param_exomiser_analysis_dict = yaml.safe_load(json_file)

                    # If param analysis is a dict
                    elif isinstance(param_exomiser_analysis, dict):
                        # Load analysis dict into analysis dict (either yaml or json)
                        param_exomiser_analysis_dict = param_exomiser_analysis

                    # Error analysis type
                    else:
                        log.error(f"Analysis type unknown. Check param file.")
                        raise ValueError(f"Analysis type unknown. Check param file.")

                # Case no input analysis config file/dict
                # Use preset (exome/genome) to open default config file
                if not param_exomiser_analysis_dict:

                    # default preset
                    default_preset = "exome"

                    # Get param preset or default preset
                    param_exomiser_preset = param_exomiser.get("preset", default_preset)

                    # Try to find if preset is a file
                    if os.path.exists(param_exomiser_preset):
                        # Preset file is provided in full path
                        param_exomiser_analysis_default_config_file = param_exomiser_preset
                    elif os.path.exists(os.path.join(folder_config, param_exomiser_preset)):
                        # Preset file is provided a basename in config folder (can be a path with subfolders)
                        param_exomiser_analysis_default_config_file = os.path.join(folder_config, param_exomiser_preset)
                    else:
                        # Construct preset file
                        param_exomiser_analysis_default_config_file = os.path.join(folder_config, f"preset-{param_exomiser_preset}-analysis.json")

                    # If preset file exists
                    if os.path.exists(param_exomiser_analysis_default_config_file):
                        # Load prest file into analysis dict (either yaml or json)
                        with open(param_exomiser_analysis_default_config_file) as json_file:
                            #param_exomiser_analysis_dict[""] = json.load(json_file)
                            param_exomiser_analysis_dict["analysis"] = yaml.safe_load(json_file)

                    # Error preset file
                    else:
                        log.error(f"No analysis preset config file ({param_exomiser_analysis_default_config_file})")
                        raise ValueError(f"No analysis preset config file ({param_exomiser_analysis_default_config_file})")

                # If no analysis dict created
                if not param_exomiser_analysis_dict:
                    log.error(f"No analysis config")
                    raise ValueError(f"No analysis config")

                # Log
                log.debug(f"Pre analysis dict: {param_exomiser_analysis_dict}")


                ### PHENOPACKET ###
                ###################

                # If no PhenoPacket in analysis dict -> check in param
                if "phenopacket" not in param_exomiser_analysis_dict:
                    
                    # If PhenoPacket in param -> load anlaysis json
                    if param_exomiser.get("phenopacket", None):

                        param_exomiser_phenopacket = param_exomiser.get("phenopacket")

                        # If param phenopacket is a file and exists
                        if isinstance(param_exomiser_phenopacket, str) and os.path.exists(param_exomiser_phenopacket):
                            # Load phenopacket file into analysis dict (either yaml or json)
                            with open(param_exomiser_phenopacket) as json_file:
                                param_exomiser_analysis_dict["phenopacket"] = yaml.safe_load(json_file)

                        # If param phenopacket is a dict
                        elif isinstance(param_exomiser_phenopacket, dict):
                            # Load phenopacket dict into analysis dict (either yaml or json)
                            param_exomiser_analysis_dict["phenopacket"] = param_exomiser_phenopacket

                        # Error phenopacket type
                        else:
                            log.error(f"Phenopacket type unknown. Check param file.")
                            raise ValueError(f"Phenopacket type unknown. Check param file.")
                    
                # If no PhenoPacket in analysis dict -> construct from sample and HPO in param
                if "phenopacket" not in param_exomiser_analysis_dict:
                    
                    # Init PhenoPacket
                    param_exomiser_analysis_dict["phenopacket"] = {
                        "id": "analysis",
                        "proband": {}
                    }

                    ### Add subject ###

                    # If subject exists
                    param_exomiser_subject = param_exomiser.get("subject", {})

                    # If subject not exists -> found sample ID
                    if not param_exomiser_subject:

                        # Found sample ID in param
                        sample = param_exomiser.get("sample", None)

                        # Find sample ID (first sample)
                        if not sample:
                            sample_list = self.get_header_sample_list()
                            if len(sample_list) > 0:
                                sample = sample_list[0]
                            else:
                                log.error(f"No sample found")
                                raise ValueError(f"No sample found")

                        # Create subject
                        param_exomiser_subject = {"id": sample,"sex": "UNKNOWN_SEX"}

                    # Add to dict
                    param_exomiser_analysis_dict["phenopacket"]["subject"] = param_exomiser_subject


                    ### Add "phenotypicFeatures" ###

                    # If phenotypicFeatures exists
                    param_exomiser_phenotypicfeatures = param_exomiser.get("phenotypicFeatures", [])

                    # If phenotypicFeatures not exists -> Try to infer from hpo list
                    if not param_exomiser_phenotypicfeatures:

                        # Found HPO in param
                        param_exomiser_hpo = param_exomiser.get("hpo", [])

                        # Split HPO if list in string format separated by comma
                        if isinstance(param_exomiser_hpo, str):
                            param_exomiser_hpo = param_exomiser_hpo.split(',')

                        # Create HPO list
                        for hpo in param_exomiser_hpo:
                            hpo_clean = re.sub("[^0-9]", '', hpo)
                            param_exomiser_phenotypicfeatures.append({
                                "type": {
                                    "id": f"HP:{hpo_clean}",
                                    "label": f"HP:{hpo_clean}"
                                }
                            })

                    # Add to dict
                    param_exomiser_analysis_dict["phenopacket"]["phenotypicFeatures"] = param_exomiser_phenotypicfeatures

                    # If phenotypicFeatures not exists -> Remove hiPhivePrioritiser step
                    if not param_exomiser_phenotypicfeatures:
                        for step in param_exomiser_analysis_dict.get("analysis", {}).get("steps", []):
                            if "hiPhivePrioritiser" in step:
                                param_exomiser_analysis_dict.get("analysis", {}).get("steps", []).remove(step)


                ### Add Input File ###

                # Initial file name and htsFiles
                tmp_vcf_name = os.path.join(tmp_dir, "initial.vcf.gz")
                param_exomiser_analysis_dict["phenopacket"]["htsFiles"] = [
                    {
                        "uri": tmp_vcf_name,
                        "htsFormat": "VCF",
                        "genomeAssembly": assembly
                    }
                ]


                ### Add metaData ###

                # If metaData not in analysis dict
                if "metaData" not in param_exomiser_analysis_dict:
                    param_exomiser_analysis_dict["phenopacket"]["metaData"] = {
                        "created": f"{datetime.datetime.now()}".replace(" ", "T") + "Z",
                        "createdBy": "howard",
                        "phenopacketSchemaVersion": 1
                    }

                ### OutputOptions ###

                # Init output result folder
                output_results = os.path.join(tmp_dir, "results")

                # If no outputOptions in analysis dict
                if "outputOptions" not in param_exomiser_analysis_dict:

                    # default output formats
                    defaut_output_formats = ["TSV_VARIANT", "VCF"]

                    # Get outputOptions in param
                    output_options = param_exomiser.get("outputOptions", None)

                    # If no output_options in param -> check 
                    if not output_options:
                        output_options = {
                            "outputContributingVariantsOnly": False,
                            "numGenes": 0,
                            "outputFormats": defaut_output_formats
                        }

                    # Replace outputDirectory in output options 
                    output_options["outputDirectory"] = output_results
                    output_options["outputFileName"] = "howard"

                    # Add outputOptions in analysis dict
                    param_exomiser_analysis_dict["outputOptions"] = output_options

                else:
                    
                    # Replace output_results and output format (if exists in param)
                    param_exomiser_analysis_dict["outputOptions"]["outputDirectory"] = output_results
                    param_exomiser_analysis_dict["outputOptions"]["outputFormats"] = list(set(param_exomiser_analysis_dict.get("outputOptions",{}).get("outputFormats", []) + ["TSV_VARIANT", "VCF"]))

                # log
                log.debug(f"Pre analysis dict: {param_exomiser_analysis_dict}")


                ### ANALYSIS FILE ###
                #####################

                ### Full JSON analysis config file ###

                exomiser_analysis = os.path.join(tmp_dir, "analysis.json")
                with open(exomiser_analysis, 'w') as fp:
                    json.dump(param_exomiser_analysis_dict, fp, indent=4)

                ### SPLIT analysis and sample config files

                # Splitted analysis dict
                param_exomiser_analysis_dict_for_split = param_exomiser_analysis_dict.copy()

                # Phenopacket JSON file
                exomiser_analysis_phenopacket = os.path.join(tmp_dir, "analysis_phenopacket.json")
                with open(exomiser_analysis_phenopacket, 'w') as fp:
                    json.dump(param_exomiser_analysis_dict_for_split.get("phenopacket"), fp, indent=4)

                # Analysis JSON file without Phenopacket parameters
                param_exomiser_analysis_dict_for_split.pop("phenopacket")
                exomiser_analysis_analysis = os.path.join(tmp_dir, "analysis_analysis.json")
                with open(exomiser_analysis_analysis, 'w') as fp:
                    json.dump(param_exomiser_analysis_dict_for_split, fp, indent=4)
                

                ### INITAL VCF file ###
                #######################

                ### Create list of samples to use and include inti initial VCF file ####

                # Subject (main sample)
                # Get sample ID in analysis dict
                sample_subject = param_exomiser_analysis_dict.get("phenopacket",{}).get("subject",{}).get("id", None)
                sample_proband = param_exomiser_analysis_dict.get("phenopacket",{}).get("proband",{}).get("subject",{}).get("id", None)
                sample = []
                if sample_subject:
                    sample.append(sample_subject)
                if sample_proband:
                    sample.append(sample_proband)

                # Get sample ID within Pedigree
                pedigree_persons_list = param_exomiser_analysis_dict.get("phenopacket",{}).get("pedigree",{}).get("persons",{})

                # Create list with all sample ID in pedigree (if exists)
                pedigree_persons = []
                for person in pedigree_persons_list:
                    pedigree_persons.append(person.get("individualId"))

                # Concat subject sample ID and samples ID in pedigreesamples
                samples = list(set(sample + pedigree_persons))

                # Check if sample list is not empty
                if not samples:
                    log.error(f"No samples found")
                    raise ValueError(f"No samples found")

                # Create VCF with sample (either sample in param or first one by default)
                # Export VCF file
                self.export_variant_vcf(vcf_file=tmp_vcf_name, file_type="gz", remove_info=True, add_samples=True, list_samples=samples, compression=1, index=False)

                
                ### Execute Exomiser ###
                ########################

                # Init command
                exomiser_command = ""

                # Command exomiser options
                exomiser_options = f" --spring.config.location={databases_folders}/{assembly}/application.properties --exomiser.data-directory={databases_folders}/{assembly} "

                # Release
                exomiser_release = param_exomiser.get("release", None)
                if exomiser_release:
                    exomiser_options += f" --exomiser.phenotype.data-version={exomiser_release} "
                    exomiser_options += f" --exomiser.{assembly}.data-version={exomiser_release} "
                    exomiser_options += f" --exomiser.{assembly}.variant-white-list-path={exomiser_release}_{assembly}_clinvar_whitelist.tsv.gz "

                # transcript_source
                transcript_source = param_exomiser.get("transcript_source", None) # ucsc, refseq, ensembl
                if transcript_source:
                    exomiser_options += f" --exomiser.{assembly}.transcript-source={transcript_source} "

                # If analysis contain proband param
                if param_exomiser_analysis_dict.get("phenopacket", {}).get("proband", {}):
                    exomiser_command = f" {exomiser_java_options} -jar {exomiser_jar} --analysis={exomiser_analysis_analysis} --sample={exomiser_analysis_phenopacket} {exomiser_options} "
                
                # If no proband (usually uniq sample)
                else:
                    exomiser_command = f" {exomiser_java_options} -jar {exomiser_jar} --analysis={exomiser_analysis} {exomiser_options}"

                # Log
                log.debug(f"{java_bin} {exomiser_command}")

                # Run command 
                result = subprocess.call([java_bin] + exomiser_command.split(), stdout=subprocess.PIPE)
                if result:
                    log.error("Exomiser command failed")
                    raise ValueError("Exomiser command failed")

                
                ### RESULTS ###
                ###############

                ### Annotate with TSV fields ###

                # Init result tsv file
                exomiser_to_info = param_exomiser.get("exomiser_to_info", False)

                # Init result tsv file
                output_results_tsv = os.path.join(output_results,"howard.variants.tsv")

                # Parse TSV file and explode columns in INFO field
                if exomiser_to_info and os.path.exists(output_results_tsv):

                    # Log
                    log.debug("Exomiser columns to VCF INFO field")

                    # Retrieve columns and types
                    query = f""" SELECT * FROM read_csv('{output_results_tsv}', auto_detect=True, delim='\t', sample_size=-1) LIMIT 0 """
                    output_results_tsv_df = self.get_query_to_df(query)
                    output_results_tsv_columns = output_results_tsv_df.columns.tolist()

                    # Init concat fields for update 
                    sql_query_update_concat_fields = []

                    # Fields to avoid
                    fields_to_avoid = ["CONTIG", "START", "END", "REF", "ALT", "QUAL", "FILTER", "GENOTYPE"]

                    # List all columns to add into header
                    for header_column in output_results_tsv_columns:

                        # If header column is enable
                        if header_column not in fields_to_avoid:

                            # Header info type
                            header_info_type = "String"
                            header_column_df = output_results_tsv_df[header_column]
                            header_column_df_dtype = header_column_df.dtype
                            if header_column_df_dtype == object:
                                if pd.to_numeric(header_column_df, errors='coerce').notnull().all():
                                    header_info_type = "Float"
                            else:
                                header_info_type = "Integer"

                            # Header info
                            characters_to_validate = ['-']
                            pattern = '[' +  ''.join(characters_to_validate) +  ']'
                            header_info_name = re.sub(pattern, '_', f"Exomiser_{header_column}".replace("#",""))
                            header_info_number = "."
                            header_info_description = f"Exomiser {header_column} annotation"
                            header_info_source = "Exomiser"
                            header_info_version = "unknown"
                            header_info_code = CODE_TYPE_MAP[header_info_type]
                            vcf_reader.infos[header_info_name] = vcf.parser._Info(
                                header_info_name,
                                header_info_number,
                                header_info_type,
                                header_info_description,
                                header_info_source,
                                header_info_version,
                                header_info_code
                            )
                            
                            # Add field to add for update to concat fields
                            sql_query_update_concat_fields.append(f"""
                                CASE
                                    WHEN table_parquet."{header_column}" NOT IN ('','.')
                                    THEN concat(
                                        '{header_info_name}=',
                                        table_parquet."{header_column}",
                                        ';'
                                        )

                                    ELSE ''
                                END
                            """)

                    # Update query
                    sql_query_update = f"""
                        UPDATE {table_variants} as table_variants
                            SET INFO = concat(
                                            CASE
                                                WHEN INFO NOT IN ('', '.')
                                                THEN INFO
                                                ELSE ''
                                            END,
                                            CASE
                                                WHEN table_variants.INFO NOT IN ('','.')
                                                THEN ';'
                                                ELSE ''
                                            END,
                                            (
                                            SELECT 
                                                concat(
                                                    {",".join(sql_query_update_concat_fields)}
                                                )
                                            FROM read_csv('{output_results_tsv}', auto_detect=True, delim='\t', sample_size=-1) as table_parquet
                                                    WHERE concat('chr', CAST(table_parquet.\"CONTIG\" AS STRING)) = table_variants.\"#CHROM\"
                                                    AND table_parquet.\"START\" = table_variants.\"POS\"
                                                    AND table_parquet.\"ALT\" = table_variants.\"ALT\"
                                                    AND table_parquet.\"REF\" = table_variants.\"REF\"
                                            )
                                        )
                            ;
                        """
                    
                    # Update
                    self.conn.execute(sql_query_update)


                ### Annotate with VCF INFO field ###
                
                # Init result VCF file
                output_results_vcf = os.path.join(output_results,"howard.vcf.gz")

                # If VCF exists
                if os.path.exists(output_results_vcf):

                    # Log
                    log.debug("Exomiser result VCF update variants")

                    # Find Exomiser INFO field annotation in header
                    with gzip.open(output_results_vcf, 'rt') as f:
                        header_list = self.read_vcf_header(f)
                    exomiser_vcf_header = vcf.Reader(io.StringIO("\n".join(header_list)))
                    
                    # Add annotation INFO field to header
                    vcf_reader.infos["Exomiser"] = exomiser_vcf_header.infos["Exomiser"]
                
                    # Update variants with VCF
                    self.update_from_vcf(output_results_vcf)
        
        return True
     

    def annotation_snpeff(self, threads: int = None) -> None:
        """
        This function annotate with snpEff

        :param threads: The number of threads to use
        :return: the value of the variable "return_value".
        """

        # DEBUG
        log.debug("Start annotation with snpeff databases")

        # Threads
        if not threads:
            threads = self.get_threads()
        log.debug("Threads: "+str(threads))

        # DEBUG
        delete_tmp = True
        if self.get_config().get("verbosity", "warning") in ["debug"]:
            delete_tmp = False
            log.debug("Delete tmp files/folders: "+str(delete_tmp))

        # Config
        config = self.get_config()
        log.debug("Config: " + str(config))

        # Config - Folders - Databases
        databases_folders = config.get("folders", {}).get(
            "databases", {}).get("snpeff", ["."])
        log.debug("Databases annotations: " + str(databases_folders))

        # Config - Java
        java_bin = config.get("tools", {}).get("java", {}).get("bin", "java")

        # Config - snpEff
        snpeff_jar = config.get("tools", {}).get(
            "snpeff", {}).get("jar", "snpEff.jar")
        snpeff_databases = config.get("folders", {}).get(
            "databases", {}).get("snpeff", DEFAULT_SNPEFF_FOLDER)

        # Config - check tools
        if not os.path.exists(java_bin):
            log.warning(
                f"Annotation warning: no java bin '{java_bin}'. Try to find...")
            # Try to find java
            try:
                java_bin = find_all('java', '/usr/bin')[0]
            except:
                log.error(f"Annotation failed: no java bin '{java_bin}'")
                raise ValueError(
                    f"Annotation failed: no java bin '{java_bin}'")
            if not os.path.exists(java_bin):
                log.error(f"Annotation failed: no java bin '{java_bin}'")
                raise ValueError(
                    f"Annotation failed: no java bin '{java_bin}'")
            else:
                log.warning(f"Annotation warning: snpEff jar bin found '{java_bin}'")

        if not os.path.exists(snpeff_jar):
            log.warning(
                f"Annotation warning: no snpEff jar '{snpeff_jar}'. Try to find...")
            # Try to find snpEff.jar
            try:
                snpeff_jar = find_all('snpEff.jar', '/')[0]
            except:
                log.error(f"Annotation failed: no snpEff jar '{snpeff_jar}'")
                raise ValueError(
                    f"Annotation failed: no snpEff jar '{snpeff_jar}'")
            if not os.path.exists(snpeff_jar):
                log.error(f"Annotation failed: no snpEff jar '{snpeff_jar}'")
                raise ValueError(
                    f"Annotation failed: no snpEff jar '{snpeff_jar}'")
            else:
                log.warning(
                    f"Annotation warning: snpEff jar found '{snpeff_jar}'")

        if snpeff_databases is not None and snpeff_databases != "":
            log.debug(f"Create snpEff databases folder")
            if not os.path.exists(snpeff_databases):
                os.makedirs(snpeff_databases)

        # Param
        param = self.get_param()
        log.debug("Param: " + str(param))

        # Param
        options = param.get("annotation", {}).get(
            "snpeff", {}).get("options", None)
        log.debug("Options: " + str(options))

        # Param - Assembly
        assembly = param.get("assembly", config.get("assembly", DEFAULT_ASSEMBLY))

        # Param - Options
        snpeff_options = param.get("annotation", {}).get(
            "snpeff", {}).get("options", "")
        snpeff_stats = param.get("annotation", {}).get(
            "snpeff", {}).get("stats", None)
        snpeff_csvstats = param.get("annotation", {}).get(
            "snpeff", {}).get("csvStats", None)
        if snpeff_stats:
            snpeff_stats = snpeff_stats.replace("OUTPUT", self.get_output())
            snpeff_options += f" -stats {snpeff_stats}"
        if snpeff_csvstats:
            snpeff_csvstats = snpeff_csvstats.replace(
                "OUTPUT", self.get_output())
            snpeff_options += f" -csvStats {snpeff_csvstats}"

        # Data
        table_variants = self.get_table_variants()

        # Check if not empty
        log.debug("Check if not empty")
        sql_query_chromosomes = f"""SELECT count(*) as count FROM {table_variants} as table_variants"""
        #if not self.conn.execute(f"{sql_query_chromosomes}").df()["count"][0]:
        if not self.get_query_to_df(f"{sql_query_chromosomes}")["count"][0]:
            log.info(f"VCF empty")
            return

        # Export in VCF
        log.debug("Create initial file to annotate")
        tmp_vcf = NamedTemporaryFile(prefix=self.get_prefix(
        ), dir=self.get_tmp_dir(), suffix=".vcf.gz", delete=True)
        tmp_vcf_name = tmp_vcf.name

        # VCF header
        vcf_reader = self.get_header()
        log.debug("Initial header: " + str(vcf_reader.infos))

        # Existing annotations
        for vcf_annotation in self.get_header().infos:

            vcf_annotation_line = self.get_header().infos.get(vcf_annotation)
            log.debug(
                f"Existing annotations in VCF: {vcf_annotation} [{vcf_annotation_line}]")

        # Memory limit
        if config.get("memory", None):
            memory_limit = config.get("memory", "8G")
        else:
            memory_limit = "8G"
        log.debug(f"memory_limit: {memory_limit}")

        # snpEff java options
        snpeff_java_options = f" -Xmx{memory_limit} -XX:+UseParallelGC -XX:ParallelGCThreads={threads} "
        log.debug(f"Exomiser java options: {snpeff_java_options}")

        force_update_annotation = True

        if "ANN" not in self.get_header().infos or force_update_annotation:

            # Check snpEff database
            log.debug(f"Check snpEff databases {[assembly]}")
            databases_download_snpeff(folder=snpeff_databases, assemblies=[assembly], config=config)

            # Export VCF file
            self.export_variant_vcf(vcf_file=tmp_vcf_name, file_type="gz",
                                    remove_info=True, add_samples=False, compression=1, index=True)

            # Tmp file
            err_files = []
            tmp_annotate_vcf = NamedTemporaryFile(prefix=self.get_prefix(
            ), dir=self.get_tmp_dir(), suffix=".vcf", delete=False)
            tmp_annotate_vcf_name = tmp_annotate_vcf.name
            tmp_annotate_vcf_name_err = tmp_annotate_vcf_name + ".err"
            err_files.append(tmp_annotate_vcf_name_err)

            # Command
            snpeff_command = f"{java_bin} {snpeff_java_options} -jar {snpeff_jar} {assembly} -dataDir {snpeff_databases} {snpeff_options} {tmp_vcf_name} 1>{tmp_annotate_vcf_name} 2>>{tmp_annotate_vcf_name_err}"
            log.debug(f"Annotation - snpEff command: {snpeff_command}")
            run_parallel_commands([snpeff_command], 1)

            # Error messages
            log.info(f"Error/Warning messages:")
            error_message_command_all = []
            error_message_command_warning = []
            error_message_command_err = []
            for err_file in err_files:
                with open(err_file, 'r') as f:
                    for line in f:
                        message = line.strip()
                        error_message_command_all.append(message)
                        if line.startswith('[W::'):
                            error_message_command_warning.append(message)
                        if line.startswith('[E::'):
                            error_message_command_err.append(
                                f"{err_file}: " + message)
            # log info
            for message in list(set(error_message_command_err + error_message_command_warning)):
                log.info(f"   {message}")
            # debug info
            for message in list(set(error_message_command_all)):
                log.debug(f"   {message}")
            # failed
            if len(error_message_command_err):
                log.error("Annotation failed: Error in commands")
                raise ValueError("Annotation failed: Error in commands")

            # Find annotation in header
            with open(tmp_annotate_vcf_name, 'rt') as f:
                header_list = self.read_vcf_header(f)
            annovar_vcf_header = vcf.Reader(
                io.StringIO("\n".join(header_list)))

            for ann in annovar_vcf_header.infos:
                if ann not in self.get_header().infos:
                    vcf_reader.infos[ann] = annovar_vcf_header.infos.get(ann)

            # Update variants
            log.info(f"Annotation - Updating...")
            self.update_from_vcf(tmp_annotate_vcf_name)

        else:
            if "ANN" in self.get_header().infos:
                log.debug(
                    f"Existing snpEff annotations in VCF")
            if force_update_annotation:
                log.debug(
                    f"Existing snpEff annotations in VCF - annotation forced")


    def annotation_annovar(self, threads: int = None) -> None:
        """
        It takes a VCF file, annotates it with Annovar, and then updates the database with the new
        annotations

        :param threads: number of threads to use
        :return: the value of the variable "return_value".
        """

        # DEBUG
        log.debug("Start annotation with Annovar databases")

        # Threads
        if not threads:
            threads = self.get_threads()
        log.debug("Threads: "+str(threads))

        # Tmp en Err files
        tmp_files = []
        err_files = []

        # DEBUG
        delete_tmp = True
        if self.get_config().get("verbosity", "warning") in ["debug"]:
            delete_tmp = False
            log.debug("Delete tmp files/folders: "+str(delete_tmp))

        # Config
        config = self.get_config()
        log.debug("Config: " + str(config))

        # Config - Folders - Databases
        databases_folders = config.get("folders", {}).get(
            "databases", {}).get("annovar", ["."])
        log.debug("Databases annotations: " + str(databases_folders))

        # Config - annovar
        annovar_bin = config.get("tools", {}).get(
            "annovar", {}).get("bin", "table_annovar.pl")
        annovar_databases = config.get("folders", {}).get(
            "databases", {}).get("annovar", DEFAULT_ANNOVAR_FOLDER)

        if not os.path.exists(annovar_bin):
            log.warning(
                f"Annotation warning: no annovar bin '{annovar_bin}'. Try to find...")
            # Try to find table_annovar.pl
            try:
                annovar_bin = find_all('table_annovar.pl', '/')[0]
            except:
                log.error(f"Annotation failed: no annovar bin '{annovar_bin}'")
                raise ValueError(
                    f"Annotation failed: no annovar bin '{annovar_bin}'")
            if not os.path.exists(annovar_bin):
                log.error(f"Annotation failed: no annovar bin '{annovar_bin}'")
                raise ValueError(
                    f"Annotation failed: no annovar bin '{annovar_bin}'")
            else:
                log.warning(
                    f"Annotation warning: annovar bin found '{annovar_bin}'")

        if annovar_databases != "" and not os.path.exists(annovar_databases):
            os.makedirs(annovar_databases)

        # Param
        param = self.get_param()
        log.debug("Param: " + str(param))

        # Param - options
        options = param.get("annotation", {}).get(
            "annovar", {}).get("options", {})
        log.debug("Options: " + str(options))

        # Param - annotations
        annotations = param.get("annotation", {}).get(
            "annovar", {}).get("annotations", {})
        log.debug("Annotations: " + str(annotations))

        # Param - Assembly
        assembly = param.get("assembly", config.get("assembly", DEFAULT_ASSEMBLY))

        # Annovar database assembly
        annovar_databases_assembly = f"{annovar_databases}/{assembly}"
        if annovar_databases_assembly != "" and not os.path.exists(annovar_databases_assembly):
            os.makedirs(annovar_databases_assembly)

        # Data
        table_variants = self.get_table_variants()

        # Check if not empty
        log.debug("Check if not empty")
        sql_query_chromosomes = f"""SELECT count(*) as count FROM {table_variants} as table_variants"""
        sql_query_chromosomes_df = self.get_query_to_df(sql_query_chromosomes)
        if not sql_query_chromosomes_df["count"][0]:
            log.info(f"VCF empty")
            return

        # VCF header
        vcf_reader = self.get_header()
        log.debug("Initial header: " + str(vcf_reader.infos))

        # Existing annotations
        for vcf_annotation in self.get_header().infos:

            vcf_annotation_line = self.get_header().infos.get(vcf_annotation)
            log.debug(
                f"Existing annotations in VCF: {vcf_annotation} [{vcf_annotation_line}]")

        force_update_annotation = True

        if annotations:

            commands = []
            tmp_annotates_vcf_name_list = []

            # Export in VCF
            log.debug("Create initial file to annotate")
            tmp_vcf = NamedTemporaryFile(prefix=self.get_prefix(
            ), dir=self.get_tmp_dir(), suffix=".vcf.gz", delete=False)
            tmp_vcf_name = tmp_vcf.name
            tmp_files.append(tmp_vcf_name)
            tmp_files.append(tmp_vcf_name+".tbi")

            # Export VCF file
            self.export_variant_vcf(vcf_file=tmp_vcf_name, file_type="gz",
                                    remove_info=".", add_samples=False, compression=1, index=True)

            # Create file for field rename
            log.debug("Create file for field rename")
            tmp_rename = NamedTemporaryFile(prefix=self.get_prefix(
            ), dir=self.get_tmp_dir(), suffix=".rename", delete=False)
            tmp_rename_name = tmp_rename.name
            tmp_files.append(tmp_rename_name)

            # Check Annovar database
            log.debug(f"Check Annovar databases {[assembly]}: {list(annotations.keys())}")
            databases_download_annovar(folder=annovar_databases, files=list(annotations.keys()), assemblies = [assembly])
            
            for annotation in annotations:
                annotation_fields = annotations[annotation]

                if not annotation_fields:
                    annotation_fields = {"INFO": None}

                log.info(f"Annotations Annovar - database '{annotation}'")
                log.debug(
                    f"Annotation '{annotation}' - fields: {annotation_fields}")

                # Tmp file for annovar
                err_files = []
                tmp_annotate_vcf_directory = TemporaryDirectory(prefix=self.get_prefix(
                ), dir=self.get_tmp_dir(), suffix=".annovar")
                tmp_annotate_vcf_prefix = tmp_annotate_vcf_directory.name + "/annovar"
                tmp_annotate_vcf_name_annovar = tmp_annotate_vcf_prefix + \
                    "." + assembly + "_multianno.vcf"
                tmp_annotate_vcf_name_err = tmp_annotate_vcf_directory.name + "/.err"
                err_files.append(tmp_annotate_vcf_name_err)
                tmp_files.append(tmp_annotate_vcf_name_err)

                # Tmp file final vcf annotated by annovar
                tmp_annotate_vcf = NamedTemporaryFile(prefix=self.get_prefix(
                ), dir=self.get_tmp_dir(), suffix=".vcf.gz", delete=False)
                tmp_annotate_vcf_name = tmp_annotate_vcf.name
                tmp_annotates_vcf_name_list.append(tmp_annotate_vcf_name)
                tmp_files.append(tmp_annotate_vcf_name)
                tmp_files.append(tmp_annotate_vcf_name+".tbi")

                # Number of fields
                annotation_list = []
                annotation_renamed_list = []

                for annotation_field in annotation_fields:

                    # field new name, if parametered SKIPPED !!!!!! not managed actually TODO
                    annotation_fields_new_name = annotation_fields.get(
                        annotation_field, annotation_field)
                    if not annotation_fields_new_name:
                        annotation_fields_new_name = annotation_field

                    if force_update_annotation or annotation_fields_new_name not in self.get_header().infos:
                        annotation_list.append(annotation_field)
                        annotation_renamed_list.append(
                            annotation_fields_new_name)
                    else:  # annotation_fields_new_name in self.get_header().infos and not force_update_annotation:
                        log.warning(
                            f"Annotation '{annotation}' - '{annotation_fields_new_name}' - already exists (skipped)")

                    # Add rename info
                    run_parallel_commands(
                        [f"echo 'INFO/{annotation_field} {annotation_fields_new_name}' >> {tmp_rename_name}"], 1)

                # log.debug("fields_to_removed: " + str(fields_to_removed))
                log.debug("annotation_list: " + str(annotation_list))

                # protocol
                protocol = annotation

                # argument
                argument = ""

                # operation
                operation = "f"
                if annotation in ["refGene", "refGeneWithVer"] or annotation.startswith("ensGene"):
                    operation = "g"
                    if options.get("genebase", None):
                        argument = f"""'{options.get("genebase","")}'"""
                elif annotation in ["cytoBand"]:
                    operation = "r"

                # argument option
                argument_option = ""
                if argument != "":
                    argument_option = " --argument " + argument

                # command options
                command_options = f""" --nastring . --vcfinput --polish --dot2underline --thread {threads} """  # --intronhgvs 10
                for option in options:
                    if option not in ["genebase"]:
                        command_options += f""" --{option}={options[option]}"""

                # Command

                # Command - Annovar
                command_annovar = f"""{annovar_bin} {tmp_vcf_name} {annovar_databases_assembly} --buildver {assembly} --outfile {tmp_annotate_vcf_prefix} --remove --protocol {protocol} --operation {operation} {argument_option} {command_options} 2>>{tmp_annotate_vcf_name_err} && mv {tmp_annotate_vcf_name_annovar} {tmp_annotate_vcf_name}.tmp.vcf """
                tmp_files.append(f"{tmp_annotate_vcf_name}.tmp.vcf")

                # Command - start pipe
                command_annovar += f""" && bcftools view --threads={threads} {tmp_annotate_vcf_name}.tmp.vcf 2>>{tmp_annotate_vcf_name_err} """

                # Command - Clean INFO/ANNOVAR_DATE (due to Annovar issue with multiple TAGS!)
                command_annovar += """ | sed "s/ANNOVAR_DATE=[^;\t]*;//gi" """

                # Command - Special characters (refGene annotation)
                command_annovar += """ | sed "s/\\\\\\x3b/,/gi" """

                # Command - Clean empty fields (with value ".")
                command_annovar += ''' | awk -F'\\t' -v OFS='\\t' '{if ($0 ~ /^#/) print; else {split($8,a,";");for(i=1;i<=length(a);i++) {split(a[i],b,"=");if(b[2]!=".") {c[b[1]]=b[2]}}; split($8,d,";");for(i=1;i<=length(d);i++) {split(d[i],e,"=");if(c[e[1]]!="") {if(info!="") {info=info";"}; info=info""e[1]"="c[e[1]]}}; if(info!="") {$8=info} else {$8=""}; delete c; info=""; print}}' '''

                # Command - Extract only needed fields, and remove ANNOVAR fields, and compress and index final file
                annovar_fields_to_keep = [
                    "INFO/ANNOVAR_DATE", "INFO/ALLELE_END"]
                if "ALL" not in annotation_list and "INFO" not in annotation_list:
                    # for ann in annotation_renamed_list:
                    for ann in annotation_list:
                        annovar_fields_to_keep.append(f"^INFO/{ann}")

                command_annovar += f""" | bcftools annotate --threads={threads} -x {",".join(annovar_fields_to_keep)} --rename-annots={tmp_rename_name} -o {tmp_annotate_vcf_name} -Oz 2>>{tmp_annotate_vcf_name_err} """

                # Command - indexing
                command_annovar += f"""  && tabix {tmp_annotate_vcf_name} """

                log.debug(f"Annotation - Annovar command: {command_annovar}")
                run_parallel_commands([command_annovar], 1)

                # Error messages
                log.info(f"Error/Warning messages:")
                error_message_command_all = []
                error_message_command_warning = []
                error_message_command_err = []
                for err_file in err_files:
                    with open(err_file, 'r') as f:
                        for line in f:
                            message = line.strip()
                            error_message_command_all.append(message)
                            if line.startswith('[W::') or line.startswith('WARNING'):
                                error_message_command_warning.append(
                                    message)
                            if line.startswith('[E::') or line.startswith('ERROR'):
                                error_message_command_err.append(
                                    f"{err_file}: " + message)
                # log info
                for message in list(set(error_message_command_err + error_message_command_warning)):
                    log.info(f"   {message}")
                # debug info
                for message in list(set(error_message_command_all)):
                    log.debug(f"   {message}")
                # failed
                if len(error_message_command_err):
                    log.error("Annotation failed: Error in commands")
                    raise ValueError(
                        "Annotation failed: Error in commands")

            if tmp_annotates_vcf_name_list:

                # List of annotated files
                tmp_annotates_vcf_name_to_merge = " ".join(
                    tmp_annotates_vcf_name_list)

                # Tmp file
                tmp_annotate_vcf = NamedTemporaryFile(prefix=self.get_prefix(
                ), dir=self.get_tmp_dir(), suffix=".vcf.gz", delete=False)
                tmp_annotate_vcf_name = tmp_annotate_vcf.name
                tmp_files.append(tmp_annotate_vcf_name)
                tmp_annotate_vcf_name_err = tmp_annotate_vcf_name + ".err"
                err_files.append(tmp_annotate_vcf_name_err)
                tmp_files.append(tmp_annotate_vcf_name_err)

                # Command merge
                merge_command = f"bcftools merge --force-samples --threads={threads} {tmp_vcf_name} {tmp_annotates_vcf_name_to_merge} -o {tmp_annotate_vcf_name} -Oz 2>>{tmp_annotate_vcf_name_err} "
                log.info(f"Annotation Annovar - Annotation merging " +
                         str(len(commands)) + " annotated files")
                log.debug(f"Annotation - merge command: {merge_command}")
                run_parallel_commands([merge_command], 1)

                # Find annotation in header
                with bgzf.open(tmp_annotate_vcf_name, 'rt') as f:
                    header_list = self.read_vcf_header(f)
                annovar_vcf_header = vcf.Reader(
                    io.StringIO("\n".join(header_list)))

                for ann in annovar_vcf_header.infos:
                    if ann not in self.get_header().infos:
                        vcf_reader.infos[ann] = annovar_vcf_header.infos.get(
                            ann)

                # Update variants
                log.info(f"Annotation Annovar - Updating...")
                self.update_from_vcf(tmp_annotate_vcf_name)

            # Clean files
            # Tmp file remove command
            if True:
                tmp_files_remove_command = ""
                if tmp_files:
                    tmp_files_remove_command = " ".join(tmp_files)
                clean_command = f" rm -f {tmp_files_remove_command} "
                log.debug(f"Annotation Annovar - Annotation cleaning ")
                log.debug(f"Annotation - cleaning command: {clean_command}")
                run_parallel_commands([clean_command], 1)

    # NEW def
    def annotation_parquet(self, threads: int = None) -> None:
        """
        It takes a VCF file, and annotates it with a parquet file

        :param threads: number of threads to use for the annotation
        :return: the value of the variable "result".
        """

        # DEBUG
        log.debug("Start annotation with parquet databases")

        # Threads
        if not threads:
            threads = self.get_threads()
        log.debug("Threads: "+str(threads))

        # DEBUG
        delete_tmp = True
        if self.get_config().get("verbosity", "warning") in ["debug"]:
            delete_tmp = False
            log.debug("Delete tmp files/folders: "+str(delete_tmp))

        # Config
        databases_folders = set(
            self.get_config().get("folders", {}).get("databases", {}).get("annotations", ["."])
            + self.get_config().get("folders", {}).get("databases", {}).get("parquet", ["."])
        )
        log.debug("Databases annotations: " + str(databases_folders))

        # Param
        annotations = self.get_param().get("annotation", {}).get(
            "parquet", {}).get("annotations", None)
        log.debug("Annotations: " + str(annotations))

        # Assembly
        assembly = self.get_param().get("assembly", self.get_config().get("assembly", DEFAULT_ASSEMBLY))

        # Data
        table_variants = self.get_table_variants()

        # Check if not empty
        log.debug("Check if not empty")
        sql_query_chromosomes_df = self.get_query_to_df(
            f"""SELECT count(*) as count FROM {table_variants} as table_variants LIMIT 1""")
        if not sql_query_chromosomes_df["count"][0]:
            log.info(f"VCF empty")
            return

        # VCF header
        vcf_reader = self.get_header()
        log.debug("Initial header: " + str(vcf_reader.infos))

        # Nb Variants POS
        log.debug("NB Variants Start")
        nb_variants = self.conn.execute(
            f"SELECT count(*) AS count FROM variants").fetchdf()["count"][0]
        log.debug("NB Variants Stop")

        # Existing annotations
        for vcf_annotation in self.get_header().infos:

            vcf_annotation_line = self.get_header().infos.get(vcf_annotation)
            log.debug(
                f"Existing annotations in VCF: {vcf_annotation} [{vcf_annotation_line}]")

        # prefix
        prefix=self.get_explode_infos_prefix()
        
        # Added columns
        added_columns = []

        # explode infos
        added_columns += self.explode_infos(prefix=prefix)

        # drop indexes
        log.debug(f"Drop indexes...")
        self.drop_indexes()
        
        if annotations:

            for annotation in annotations:

                # Annotation Name
                annotation_name = os.path.basename(annotation)
                            
                # Annotation fields
                annotation_fields = annotations[annotation]
                if not annotation_fields:
                    annotation_fields = {"INFO": None}

                log.debug(f"Annotation '{annotation_name}'")
                log.debug(
                    f"Annotation '{annotation_name}' - fields: {annotation_fields}")

                # Create Database
                database = Database(database=annotation, databases_folders=databases_folders, assembly=assembly)

                # Find files
                parquet_file = database.get_database()
                parquet_hdr_file = database.get_header_file()
                parquet_format = database.get_format()
                parquet_type = database.get_type()
                
                # Check if files exists
                if not parquet_file or not parquet_hdr_file:
                    log.error("Annotation failed: file not found")
                    raise ValueError("Annotation failed: file not found")
                else:
                    # Get parquet connexion
                    parquet_sql_attach = database.get_sql_database_attach(output="query")
                    if parquet_sql_attach:
                        self.conn.execute(parquet_sql_attach)
                    parquet_file_link = database.get_sql_database_link()
                    # Log
                    log.debug(f"Annotation '{annotation_name}' - file: " +
                              str(parquet_file) + " and " + str(parquet_hdr_file))

                    # Database full header columns
                    parquet_hdr_vcf_header_columns = database.get_header_file_columns(parquet_hdr_file)
                    # Log
                    log.debug("Annotation database header columns : " +
                              str(parquet_hdr_vcf_header_columns))

                    # Load header as VCF object
                    parquet_hdr_vcf_header_infos = database.get_header().infos
                    # Log
                    log.debug("Annotation database header: " +
                              str(parquet_hdr_vcf_header_infos))

                    # Get extra infos
                    parquet_columns = database.get_extra_columns()
                    # Log
                    log.debug("Annotation database Columns: " +
                              str(parquet_columns))
                    
                    # Add extra columns if "ALL" in annotation_fields
                    # if "ALL" in annotation_fields:
                    #     allow_add_extra_column = True
                    if "ALL" in annotation_fields and database.get_extra_columns():
                        for extra_column in database.get_extra_columns():
                            if extra_column not in annotation_fields and extra_column.replace("INFO/","") not in parquet_hdr_vcf_header_infos:
                                parquet_hdr_vcf_header_infos[extra_column] = vcf.parser._Info(
                                    extra_column,
                                    ".",
                                    "String",
                                    f"{extra_column} description",
                                    "unknown",
                                    "unknown",
                                    self.code_type_map["String"]
                                )
                    
                    # For all fields in database
                    annotation_fields_ALL = False
                    if "ALL" in annotation_fields or "INFO" in annotation_fields:
                        annotation_fields_ALL = True
                        annotation_fields = {
                            key: key for key in parquet_hdr_vcf_header_infos}
                        
                        log.debug(
                            "Annotation database header - All annotations added: " + str(annotation_fields))

                    # Init 
                    
                    # List of annotation fields to use
                    sql_query_annotation_update_info_sets = []

                    # List of annotation to agregate
                    sql_query_annotation_to_agregate = []

                    # Number of fields
                    nb_annotation_field = 0

                    # Annotation fields processed
                    annotation_fields_processed = []
                    
                    # Columns mapping
                    map_columns = database.map_columns(columns=annotation_fields, prefixes=["INFO/"])

                    # Fetch Anotation fields
                    for annotation_field in annotation_fields:
                        
                        # annotation_field_column
                        annotation_field_column = map_columns.get(annotation_field, "INFO")

                        # field new name, if parametered
                        annotation_fields_new_name = annotation_fields.get(
                            annotation_field, annotation_field)
                        if not annotation_fields_new_name:
                            annotation_fields_new_name = annotation_field

                        # To annotate
                        force_update_annotation = False
                        if annotation_field in parquet_hdr_vcf_header_infos and (force_update_annotation or (annotation_fields_new_name not in self.get_header().infos)):
                            
                            # Add field to annotation to process list
                            annotation_fields_processed.append(
                                annotation_fields_new_name)

                            # Sep between fields in INFO
                            nb_annotation_field += 1
                            if nb_annotation_field > 1:
                                annotation_field_sep = ";"
                            else:
                                annotation_field_sep = ""

                            log.info(
                                f"Annotation '{annotation_name}' - '{annotation_field}' -> 'INFO/{annotation_fields_new_name}'")

                            # Add INFO field to header
                            parquet_hdr_vcf_header_infos_number = parquet_hdr_vcf_header_infos[
                                annotation_field].num or "."
                            parquet_hdr_vcf_header_infos_type = parquet_hdr_vcf_header_infos[
                                annotation_field].type or "String"
                            parquet_hdr_vcf_header_infos_description = parquet_hdr_vcf_header_infos[
                                annotation_field].desc or f"{annotation_field} description"
                            parquet_hdr_vcf_header_infos_source = parquet_hdr_vcf_header_infos[
                                annotation_field].source or "unknown"
                            parquet_hdr_vcf_header_infos_version = parquet_hdr_vcf_header_infos[
                                annotation_field].version or "unknown"

                            vcf_reader.infos[annotation_fields_new_name] = vcf.parser._Info(
                                annotation_fields_new_name,
                                parquet_hdr_vcf_header_infos_number,
                                parquet_hdr_vcf_header_infos_type,
                                parquet_hdr_vcf_header_infos_description,
                                parquet_hdr_vcf_header_infos_source,
                                parquet_hdr_vcf_header_infos_version,
                                self.code_type_map[parquet_hdr_vcf_header_infos_type]
                            )

                            # Annotation/Update query fields
                            # Found in INFO column
                            if annotation_field_column == "INFO" and "INFO" in parquet_hdr_vcf_header_columns:
                                sql_query_annotation_update_info_sets.append(f"""
                                CASE WHEN REGEXP_EXTRACT(concat(';', table_parquet.INFO), ';{annotation_field}=([^;]*)',1) NOT IN ('','.')
                                        THEN concat('{annotation_field_sep}', '{annotation_fields_new_name}=', REGEXP_EXTRACT(concat(';', table_parquet.INFO), ';{annotation_field}=([^;]*)',1))
                                        ELSE ''
                                    END
                                """)
                            # Found in a specific column
                            else:
                                sql_query_annotation_update_info_sets.append(f"""
                                CASE WHEN table_parquet."{annotation_field_column}" NOT IN ('','.')
                                        THEN concat('{annotation_field_sep}', '{annotation_fields_new_name}=', replace(table_parquet."{annotation_field_column}", ';', ','))
                                        ELSE ''
                                    END
                                """)
                                sql_query_annotation_to_agregate.append(f""" string_agg(DISTINCT table_parquet_from."{annotation_field_column}", ',') AS "{annotation_field_column}" """)

                        # Not to annotate
                        else:

                            if force_update_annotation:
                                annotation_message = "forced"
                            else:
                                annotation_message = "skipped"

                            if annotation_field not in parquet_hdr_vcf_header_infos:
                                log.warning(
                                    f"Annotation '{annotation_name}' - '{annotation_field}' [{nb_annotation_field}] - not available in parquet file")
                            if annotation_fields_new_name in self.get_header().infos:
                                log.warning(
                                    f"Annotation '{annotation_name}' - '{annotation_fields_new_name}' [{nb_annotation_field}] - already exists in header ({annotation_message})")

                    # Check if ALL fields have to be annotated. Thus concat all INFO field
                    allow_annotation_full_info = True
                    
                    if parquet_type in ["regions"]:
                        allow_annotation_full_info = False

                    if allow_annotation_full_info and nb_annotation_field == len(annotation_fields) and annotation_fields_ALL and ("INFO" in parquet_hdr_vcf_header_columns and "INFO" in database.get_extra_columns()):
                        log.debug("Column INFO annotation enabled")
                        sql_query_annotation_update_info_sets = []
                        sql_query_annotation_update_info_sets.append(
                            f" table_parquet.INFO ")

                    if sql_query_annotation_update_info_sets:

                        # Annotate
                        log.info(f"Annotation '{annotation_name}' - Annotation...")

                        # Join query annotation update info sets for SQL
                        sql_query_annotation_update_info_sets_sql = ",".join(
                            sql_query_annotation_update_info_sets)

                        # Check chromosomes list (and variant max position)
                        sql_query_chromosomes_max_pos = f""" SELECT table_variants."#CHROM" as CHROM, MAX(table_variants."POS") as MAX_POS, MIN(table_variants."POS")-1 as MIN_POS FROM {table_variants} as table_variants GROUP BY table_variants."#CHROM" """
                        sql_query_chromosomes_max_pos_df = self.conn.execute(
                            sql_query_chromosomes_max_pos).df()

                        # Create dictionnary with chromosomes (and max position)
                        sql_query_chromosomes_max_pos_dictionary = sql_query_chromosomes_max_pos_df.groupby('CHROM').apply(
                            lambda x: {'max_pos': x['MAX_POS'].max(), 'min_pos': x['MIN_POS'].min()}).to_dict()

                        # Affichage du dictionnaire
                        log.debug("Chromosomes max pos found: " +
                                  str(sql_query_chromosomes_max_pos_dictionary))

                        # nb_of_variant_annotated
                        nb_of_query = 0
                        nb_of_variant_annotated = 0
                        query_dict = {}

                        for chrom in sql_query_chromosomes_max_pos_dictionary:

                            # nb_of_variant_annotated_by_chrom
                            nb_of_variant_annotated_by_chrom = 0

                            # Get position of the farthest variant (max position) in the chromosome
                            sql_query_chromosomes_max_pos_dictionary_max_pos = sql_query_chromosomes_max_pos_dictionary.get(
                                chrom, {}).get("max_pos")
                            sql_query_chromosomes_max_pos_dictionary_min_pos = sql_query_chromosomes_max_pos_dictionary.get(
                                chrom, {}).get("min_pos")

                            # Autodetect range of bases to split/chunk
                            log.debug(
                                f"Annotation '{annotation_name}' - Chromosome '{chrom}' - Start Autodetection Intervals...")

                            batch_annotation_databases_step = None
                            batch_annotation_databases_ncuts = 1

                            # Create intervals from 0 to max position variant, with the batch window previously defined
                            sql_query_intervals = split_interval(
                                sql_query_chromosomes_max_pos_dictionary_min_pos, sql_query_chromosomes_max_pos_dictionary_max_pos, step=batch_annotation_databases_step, ncuts=batch_annotation_databases_ncuts)

                            log.debug(
                                f"Annotation '{annotation_name}' - Chromosome '{chrom}' - Stop Autodetection Intervals")

                            # Interval Start/Stop
                            sql_query_interval_start = sql_query_intervals[0]

                            # For each interval
                            for i in sql_query_intervals[1:]:

                                # Interval Start/Stop
                                sql_query_interval_stop = i

                                log.debug(
                                    f"Annotation '{annotation_name}' - Chromosome '{chrom}' - Interval [{sql_query_interval_start}-{sql_query_interval_stop}] ...")

                                log.debug(
                                    f"Annotation '{annotation_name}' - Chromosome '{chrom}' - Interval [{sql_query_interval_start}-{sql_query_interval_stop}] - Start detecting regions...")

                                regions = [
                                    (chrom, sql_query_interval_start, sql_query_interval_stop)]

                                log.debug(
                                    f"Annotation '{annotation_name}' - Chromosome '{chrom}' - Interval [{sql_query_interval_start}-{sql_query_interval_stop}] - Stop detecting regions")

                                # Fusion des régions chevauchantes
                                if regions:

                                    # Number of regions
                                    nb_regions = len(regions)

                                    # create where caluse on regions
                                    clause_where_regions_variants = create_where_clause(
                                        regions, table="table_variants")

                                    log.debug(
                                        f"Annotation '{annotation_name}' - Chromosome '{chrom}' - Interval [{sql_query_interval_start}-{sql_query_interval_stop}] - {nb_regions} regions...")

                                    # Annotation with regions database
                                    if parquet_type in ["regions"]:
                                        sql_query_annotation_from_clause = f"""
                                            FROM (
                                                SELECT 
                                                    '{chrom}' AS \"#CHROM\",
                                                    table_variants_from.\"POS\" AS \"POS\",
                                                    {",".join(sql_query_annotation_to_agregate)}
                                                FROM {table_variants} as table_variants_from
                                                LEFT JOIN {parquet_file_link} as table_parquet_from ON (
                                                    table_parquet_from.\"#CHROM\" in ('{chrom}')
                                                    AND table_variants_from.\"POS\" <= table_parquet_from.\"END\"
                                                    AND (table_variants_from.\"POS\" >= (table_parquet_from.\"START\"+1)
                                                        OR table_variants_from.\"POS\" + (len(table_variants_from.\"REF\")-1) >= (table_parquet_from.\"START\"+1)
                                                        )
                                                )
                                                WHERE table_variants_from.\"#CHROM\" in ('{chrom}')
                                                GROUP BY table_variants_from.\"POS\"
                                                )
                                                as table_parquet
                                        """
                                        
                                        sql_query_annotation_where_clause = """
                                            table_parquet.\"#CHROM\" = table_variants.\"#CHROM\"
                                            AND table_parquet.\"POS\" = table_variants.\"POS\"
                                        """

                                    # Annotation with variants database
                                    else:
                                        sql_query_annotation_from_clause = f"""
                                            FROM {parquet_file_link} as table_parquet 
                                        """
                                        sql_query_annotation_where_clause = f"""
                                            table_parquet.\"#CHROM\" in ('{chrom}')
                                            AND ( {clause_where_regions_variants} )
                                            AND table_parquet.\"#CHROM\" = table_variants.\"#CHROM\"
                                            AND table_parquet.\"POS\" = table_variants.\"POS\"
                                            AND table_parquet.\"ALT\" = table_variants.\"ALT\"
                                            AND table_parquet.\"REF\" = table_variants.\"REF\"
                                        """

                                    # Create update query
                                    sql_query_annotation_chrom_interval_pos = f"""
                                        UPDATE {table_variants} as table_variants
                                            SET INFO = 
                                                concat(
                                                    CASE WHEN table_variants.INFO NOT IN ('','.')
                                                        THEN table_variants.INFO
                                                        ELSE ''
                                                    END
                                                    ,
                                                    CASE WHEN table_variants.INFO NOT IN ('','.')
                                                              AND (
                                                                concat({sql_query_annotation_update_info_sets_sql})
                                                                )
                                                                NOT IN ('','.') 
                                                         THEN ';'
                                                         ELSE ''
                                                    END
                                                    ,
                                                    {sql_query_annotation_update_info_sets_sql}
                                                    )
                                            {sql_query_annotation_from_clause}
                                            WHERE {sql_query_annotation_where_clause}
                                            ;
                                        """
                                    
                                    # Add update query to dict
                                    query_dict[f"{chrom}:{sql_query_interval_start}-{sql_query_interval_stop}"] = sql_query_annotation_chrom_interval_pos

                                    log.debug(
                                        "Create SQL query: " + str(sql_query_annotation_chrom_interval_pos))

                                    # Interval Start/Stop
                                    sql_query_interval_start = sql_query_interval_stop

                            # nb_of_variant_annotated
                            nb_of_variant_annotated += nb_of_variant_annotated_by_chrom

                        nb_of_query = len(query_dict)
                        num_query = 0

                        # SET max_expression_depth TO x
                        self.conn.execute("SET max_expression_depth TO 10000")

                        for query_name in query_dict:
                            query = query_dict[query_name]
                            num_query += 1
                            log.info(
                                f"Annotation '{annotation_name}' - Annotation - Query [{num_query}/{nb_of_query}] {query_name}...")
                            result = self.conn.execute(query)
                            nb_of_variant_annotated_by_query = result.df()[
                                "Count"][0]
                            nb_of_variant_annotated += nb_of_variant_annotated_by_query
                            log.info(
                                f"Annotation '{annotation_name}' - Annotation - Query [{num_query}/{nb_of_query}] {query_name} - {nb_of_variant_annotated_by_query} variants annotated")

                        log.info(
                            f"Annotation '{annotation_name}' - Annotation of {nb_of_variant_annotated} variants out of {nb_variants} (with {nb_of_query} queries)")

                    else:

                        log.info(
                            f"Annotation '{annotation_name}' - No Annotations available")

                    log.debug("Final header: " + str(vcf_reader.infos))


    ###
    # Prioritization
    ###

    def prioritization(self) -> None:
        """
        It takes a VCF file, and adds a bunch of new INFO fields to it, based on the values of other
        INFO fields
        """

        log.info(f"Prioritization... ")

        # Config
        config = self.get_config()
        config_profiles = config.get(
            "prioritization", {}).get("config_profiles", None)

        # Param
        param = self.get_param()
        config_profiles = param.get("prioritization", {}).get(
            "config_profiles", config_profiles)
        profiles = param.get("prioritization", {}).get("profiles", None)
        pzfields = param.get("prioritization", {}).get(
            "pzfields", ["PZFlag", "PZScore"])
        default_profile = param.get("prioritization", {}).get(
            "default_profile", None)
        pzfields_sep = param.get("prioritization", {}).get("pzfields_sep", "_")
        prioritization_score_mode = param.get("prioritization", {}).get(
            "prioritization_score_mode", "HOWARD")

        # Profiles are in files
        if config_profiles and os.path.exists(config_profiles):
            with open(config_profiles) as profiles_file:
                config_profiles = json.load(profiles_file)
        else:
            log.error("NO Profiles configuration")
            raise ValueError(f"NO Profiles configuration")

        # If no profiles provided, all profiles in the config profiles
        if not profiles:
            profiles = list(config_profiles.keys())

        if not default_profile:
            default_profile = profiles[0]

        log.debug("Profiles availables: " + str(list(config_profiles.keys())))
        log.debug("Profiles to check: " + str(list(profiles)))

        # Variables
        table_variants = self.get_table_variants(clause="update")

        # Added columns
        added_columns = []

        # Create list of PZfields
        # List of PZFields
        list_of_pzfields_original = pzfields + \
            [pzfield+pzfields_sep+profile for pzfield in pzfields for profile in profiles]
        list_of_pzfields = []
        log.debug(f"{list_of_pzfields_original}")

        # Remove existing PZfields to use if exists
        for pzfield in list_of_pzfields_original:
            if self.get_header().infos.get(pzfield, None) is None:
                list_of_pzfields.append(pzfield)
                log.debug(
                    f"VCF Input - Header - PZfield '{pzfield}' not in VCF")
            else:
                log.debug(
                    f"VCF Input - Header - PZfield '{pzfield}' already in VCF")

        if list_of_pzfields:

            # Explode Infos fields
            explode_infos_prefix = self.get_explode_infos_prefix()
            added_columns += self.explode_infos(prefix=explode_infos_prefix)
            extra_infos = self.get_extra_infos()


            # PZfields tags description
            PZfields_INFOS = {
                'PZTags': {
                    'ID': 'PZTags',
                    'Number': '.',
                    'Type': 'String',
                    'Description': 'Variant tags based on annotation criteria'
                },
                'PZScore': {
                    'ID': 'PZScore',
                    'Number': 1,
                    'Type': 'Integer',
                    'Description': 'Variant score based on annotation criteria'
                },
                'PZFlag': {
                    'ID': 'PZFlag',
                    'Number': 1,
                    'Type': 'String',
                    'Description': 'Variant flag based on annotation criteria'
                },
                'PZComment': {
                    'ID': 'PZComment',
                    'Number': '.',
                    'Type': 'String',
                    'Description': 'Variant comment based on annotation criteria'
                },
                'PZInfos': {
                    'ID': 'PZInfos',
                    'Number': '.',
                    'Type': 'String',
                    'Description': 'Variant infos based on annotation criteria'
                }
            }

            # Create INFO fields if not exist
            for field in PZfields_INFOS:
                field_ID = PZfields_INFOS[field]["ID"]
                field_description = PZfields_INFOS[field]["Description"]
                if field_ID not in self.get_header().infos and field_ID in pzfields:
                    if field != "PZTags":
                        field_description = PZfields_INFOS[field]["Description"] + \
                            f", profile {default_profile}"
                    self.get_header().infos[field_ID] = vcf.parser._Info(
                        field_ID, PZfields_INFOS[field]["Number"], PZfields_INFOS[field]["Type"], field_description, 'unknown', 'unknown', code_type_map[PZfields_INFOS[field]["Type"]])

            # Create INFO fields if not exist for each profile
            for profile in config_profiles:
                if profile in profiles or profiles == []:
                    for field in PZfields_INFOS:
                        if field != "PZTags":
                            field_ID = PZfields_INFOS[field]["ID"] + \
                                pzfields_sep+profile
                            field_description = PZfields_INFOS[field]["Description"] + \
                                f", profile {profile}"
                            if field_ID not in self.get_header().infos and field in pzfields:
                                self.get_header().infos[field_ID] = vcf.parser._Info(
                                    field_ID, PZfields_INFOS[field]["Number"], PZfields_INFOS[field]["Type"], field_description, 'unknown', 'unknown', code_type_map[PZfields_INFOS[field]["Type"]])

            # Header
            for pzfield in list_of_pzfields:
                if re.match("PZScore.*", pzfield):
                    added_column = self.add_column(table_name=table_variants, column_name=pzfield, column_type="INTEGER", default_value="0")
                elif re.match("PZFlag.*", pzfield):
                    added_column = self.add_column(table_name=table_variants, column_name=pzfield, column_type="BOOLEAN", default_value="1")
                else:
                    added_column = self.add_column(table_name=table_variants, column_name=pzfield, column_type="STRING", default_value="''")
                added_columns.append(added_column)

            # Profiles
            if profiles:

                # foreach profile in configuration file
                for profile in config_profiles:

                    # If profile is asked in param, or ALL are asked (empty profile [])
                    if profile in profiles or profiles == []:
                        log.info(f"Profile '{profile}'")

                        sql_set_info_option = ""
                        
                        sql_set_info = []

                        # PZ fields set

                        # PZScore
                        if f"PZScore{pzfields_sep}{profile}" in list_of_pzfields:
                            sql_set_info.append(
                                f"""
                                    concat(
                                        'PZScore{pzfields_sep}{profile}=',
                                        PZScore{pzfields_sep}{profile}
                                    ) 
                                """)
                            if profile == default_profile and "PZScore" in list_of_pzfields:
                                sql_set_info.append(
                                    f"""
                                        concat(
                                            'PZScore=',
                                            PZScore{pzfields_sep}{profile}
                                        )
                                    """)
                                
                        # PZFlag
                        if f"PZFlag{pzfields_sep}{profile}" in list_of_pzfields:
                            sql_set_info.append(
                                f"""
                                    concat(
                                        'PZFlag{pzfields_sep}{profile}=',
                                        CASE 
                                            WHEN PZFlag{pzfields_sep}{profile}==1
                                            THEN 'PASS'
                                            WHEN PZFlag{pzfields_sep}{profile}==0
                                            THEN 'FILTERED'
                                        END
                                    ) 
                                """)
                            if profile == default_profile and "PZFlag" in list_of_pzfields:
                                sql_set_info.append(
                                    f"""
                                        concat(
                                            'PZFlag=',
                                            CASE 
                                                WHEN PZFlag{pzfields_sep}{profile}==1
                                                THEN 'PASS'
                                                WHEN PZFlag{pzfields_sep}{profile}==0
                                                THEN 'FILTERED'
                                            END
                                        )
                                    """)

                        # PZComment
                        if f"PZComment{pzfields_sep}{profile}" in list_of_pzfields:
                            sql_set_info.append(
                                f"""
                                    CASE
                                        WHEN PZComment{pzfields_sep}{profile} NOT IN ('')
                                        THEN concat('PZComment{pzfields_sep}{profile}=', PZComment{pzfields_sep}{profile})
                                        ELSE ''
                                    END
                                """)
                            if profile == default_profile and "PZComment" in list_of_pzfields:
                                sql_set_info.append(
                                    f"""
                                        CASE
                                            WHEN PZComment{pzfields_sep}{profile} NOT IN ('')
                                            THEN concat('PZComment=', PZComment{pzfields_sep}{profile})
                                            ELSE ''
                                        END
                                    """)
                        
                        # PZInfos
                        if f"PZInfos{pzfields_sep}{profile}" in list_of_pzfields:
                            sql_set_info.append(
                                f"""
                                    CASE
                                        WHEN PZInfos{pzfields_sep}{profile} NOT IN ('')
                                        THEN concat('PZInfos{pzfields_sep}{profile}=', PZInfos{pzfields_sep}{profile})
                                        ELSE ''
                                    END
                                """)
                            if profile == default_profile and "PZInfos" in list_of_pzfields:
                                sql_set_info.append(
                                    f"""
                                        CASE
                                            WHEN PZInfos{pzfields_sep}{profile} NOT IN ('')
                                            THEN concat('PZInfos=', PZInfos{pzfields_sep}{profile})
                                            ELSE ''
                                        END
                                    """)
                                
                        # Merge PZfields
                        sql_set_info_option = ""
                        sql_set_sep = ""
                        for sql_set in sql_set_info:
                            if sql_set_sep:
                                sql_set_info_option += f"""
                                    , concat('{sql_set_sep}', {sql_set})
                                """
                            else:
                                sql_set_info_option += f"""
                                    , {sql_set}
                                """
                            sql_set_sep = ";"


                        sql_queries = []
                        for annotation in config_profiles[profile]:

                            # Check if annotation field is present
                            if not f"{explode_infos_prefix}{annotation}" in extra_infos:
                                log.debug(
                                    f"Annotation '{annotation}' not in data")
                                continue
                            else:
                                log.debug(
                                    f"Annotation '{annotation}' in data")

                            # For each criterions
                            for criterion in config_profiles[profile][annotation]:
                                criterion_type = criterion['type']
                                criterion_value = criterion['value']
                                criterion_score = criterion.get('score', 0)
                                criterion_flag = criterion.get('flag', 'PASS')
                                criterion_flag_bool = (
                                    criterion_flag == "PASS")
                                criterion_comment = ", ".join(criterion.get('comment', [])).replace(
                                    '\'', '\'\'').replace(';', ',').replace('\t', ' ')
                                criterion_infos = str(criterion).replace(
                                    '\'', '\'\'').replace(';', ',').replace('\t', ' ')

                                sql_set = []
                                sql_set_info = []

                                # PZ fields set
                                if f"PZScore{pzfields_sep}{profile}" in list_of_pzfields:
                                    if prioritization_score_mode == "HOWARD":
                                        sql_set.append(
                                            f"PZScore{pzfields_sep}{profile} = PZScore{pzfields_sep}{profile} + {criterion_score}")
                                    elif prioritization_score_mode == "VaRank":
                                        sql_set.append(
                                            f"PZScore{pzfields_sep}{profile} = CASE WHEN {criterion_score}>PZScore{pzfields_sep}{profile} THEN {criterion_score} END")
                                    else:
                                        sql_set.append(
                                            f"PZScore{pzfields_sep}{profile} = PZScore{pzfields_sep}{profile} + {criterion_score}")
                                if f"PZFlag{pzfields_sep}{profile}" in list_of_pzfields:
                                    sql_set.append(
                                        f"PZFlag{pzfields_sep}{profile} = PZFlag{pzfields_sep}{profile} AND {criterion_flag_bool}")
                                if f"PZComment{pzfields_sep}{profile}" in list_of_pzfields:
                                    sql_set.append(
                                        f"""
                                            PZComment{pzfields_sep}{profile} = 
                                                concat(
                                                    PZComment{pzfields_sep}{profile},
                                                    CASE 
                                                        WHEN PZComment{pzfields_sep}{profile}!=''
                                                        THEN ', '
                                                        ELSE ''
                                                    END,
                                                    '{criterion_comment}'
                                                )
                                        """
                                        )
                                if f"PZInfos{pzfields_sep}{profile}" in list_of_pzfields:
                                    sql_set.append(
                                        f"""
                                            PZInfos{pzfields_sep}{profile} = 
                                                concat(
                                                    PZInfos{pzfields_sep}{profile},
                                                    '{criterion_infos}'
                                                )
                                        """
                                        )
                                sql_set_option = ",".join(sql_set)
                                
                                # Criterion and comparison
                                try:
                                    float(criterion_value)
                                    sql_update = f"""
                                        UPDATE {table_variants} \
                                        SET {sql_set_option} \
                                        WHERE "{explode_infos_prefix}{annotation}" NOT IN ('','.') \
                                        AND "{explode_infos_prefix}{annotation}"{comparison_map[criterion_type]}{criterion_value}
                                    """
                                except:
                                    contains_option = ""
                                    if criterion_type == "contains":
                                        contains_option = ".*"
                                    sql_update = f"""
                                    UPDATE {table_variants} \
                                        SET {sql_set_option} \
                                        WHERE "{explode_infos_prefix}{annotation}" SIMILAR TO '{contains_option}{criterion_value}{contains_option}'
                                        """
                                sql_queries.append(sql_update)

                        log.info(
                            f"""Profile '{profile}' - Prioritization... """)

                        if sql_queries:

                            for sql_query in sql_queries:
                                log.debug(f"""Profile '{profile}' - Prioritization query: {sql_query}... """)
                                self.conn.execute(sql_query)

                        log.info(f"""Profile '{profile}' - Update... """)
                        sql_query_update = f"""
                            UPDATE {table_variants}
                            SET INFO =  
                                concat(
                                    CASE
                                        WHEN INFO NOT IN ('','.')
                                        THEN concat(INFO, ';')
                                        ELSE ''
                                    END
                                    {sql_set_info_option}
                                )
                        """
                        self.conn.execute(sql_query_update)

        else:

            log.warning(f"No profiles in parameters")


        # Remove added columns
        for added_column in added_columns:
            self.drop_column(column=added_column)

        # Explode INFOS fields into table fields
        if self.get_explode_infos():
            self.explode_infos(prefix=self.get_explode_infos_prefix(), fields=self.get_explode_infos_fields(), force=True)


    ###
    # HGVS
    ###


    def annotation_hgvs(self, threads:int = None) -> None:
        """
        The `annotation_hgvs` function performs HGVS annotation on a set of variants using genomic
        coordinates and alleles.
        
        :param threads: The `threads` parameter is an optional integer that specifies the number of
        threads to use for parallel processing. If no value is provided, it will default to the number
        of threads obtained from the `get_threads()` method
        :type threads: int
        """

        # Function for each partition of the Dask Dataframe
        def partition_function(partition):
            """
            The function `partition_function` applies the `annotation_hgvs_partition` function to
            each row of a DataFrame called `partition`.
            
            :param partition: The parameter "partition" is a pandas DataFrame that contains the data
            to be processed
            :return: the result of applying the "annotation_hgvs_partition" function to each row of
            the "partition" dataframe along the axis 1.
            """
            return partition.apply(annotation_hgvs_partition, axis=1)

        def annotation_hgvs_partition(row) -> str:
            """
            The function `annotation_hgvs_partition` takes in a row of data and returns a string
            containing a list of HGVS names associated with the given genomic coordinates and alleles.
            
            :param row: A dictionary-like object that contains the values for the following keys:
            :return: a string that contains the HGVS names associated with the given row of data.
            """

            chr = row["CHROM"]
            pos = row["POS"]
            ref = row["REF"]
            alt = row["ALT"]

            # Find list of associated transcripts
            transcripts_list = list(polars_conn.execute(f"""
                SELECT transcript
                FROM refseq_df
                WHERE CHROM='{chr}'
                AND POS={pos} 
            """)["transcript"])


            # Full HGVS annotation in list
            hgvs_full_list = []

            for transcript_name in transcripts_list:

                # Transcript
                transcript = get_transcript(transcripts=transcripts, transcript_name=transcript_name)
                # Exon
                if use_exon:
                    exon=transcript.find_exon_number(pos)
                else:
                    exon = None
                # Protein
                transcript_protein = None
                if use_protein or add_protein or full_format:
                    transcripts_protein = list(polars_conn.execute(f"""
                        SELECT protein
                        FROM refseqlink_df
                        WHERE transcript='{transcript_name}'
                        LIMIT 1
                    """)["protein"])
                    if len(transcripts_protein):
                        transcript_protein = transcripts_protein[0]

                # HGVS name
                hgvs_name = format_hgvs_name(chr, pos, ref, alt, genome=genome, transcript=transcript, transcript_protein=transcript_protein, exon=exon, use_gene=use_gene, use_protein=use_protein, full_format=full_format, use_version=use_version, codon_type=codon_type)
                hgvs_full_list.append(hgvs_name)
                if add_protein and not use_protein and not full_format:
                    hgvs_name = format_hgvs_name(chr, pos, ref, alt, genome=genome, transcript=transcript, transcript_protein=transcript_protein, exon=exon, use_gene=use_gene, use_protein=True, full_format=False, use_version=use_version, codon_type=codon_type)
                    hgvs_full_list.append(hgvs_name)

            # Create liste of HGVS annotations
            hgvs_full = ",".join(hgvs_full_list)

            return hgvs_full

        log.info(f"HGVS Annotation... ")

        # Polars connexion
        polars_conn = pl.SQLContext(register_globals=True, eager_execution=True)

        # Config
        config = self.get_config()
        
        # Databases
        # Genome
        databases_genomes_folders = config.get("folders", {}).get(
            "databases", {}).get("genomes", DEFAULT_GENOME_FOLDER)
        databases_genome = config.get("databases", {}).get("genome", "")
        # refseq database folder
        databases_refseq_folders = config.get("folders", {}).get(
            "databases", {}).get("refseq", DEFAULT_REFSEQ_FOLDER)
        # refseq
        databases_refseq = config.get("databases", {}).get("refSeq", "")
        # refSeqLink
        databases_refseqlink = config.get("databases", {}).get("refSeqLink", "")
        
        # Param
        param = self.get_param()
        
        # HGVS Param
        param_hgvs = param.get("hgvs",{})
        use_exon = param_hgvs.get("use_exon",False)
        use_gene = param_hgvs.get("use_gene",False)
        use_protein = param_hgvs.get("use_protein",False)
        add_protein = param_hgvs.get("add_protein",False)
        full_format = param_hgvs.get("full_format",False)
        use_version = param_hgvs.get("use_version",False)
        codon_type = param_hgvs.get("codon_type","3")

        # refSseq refSeqLink
        databases_refseq = param_hgvs.get("refseq", databases_refseq)
        databases_refseqlink = param_hgvs.get("refseqlink", databases_refseqlink)

        # Assembly
        assembly = param.get("assembly", config.get("assembly", DEFAULT_ASSEMBLY))

        # Genome
        genome_file = None
        if find_genome(databases_genome):
            genome_file = find_genome(databases_genome)
        else:
            genome_file = find_genome(genome_path=databases_genomes_folders, assembly=assembly)
        log.debug("Genome: "+str(genome_file))
        
        # refSseq
        refseq_file = find_file_prefix(input_file=databases_refseq, prefix="ncbiRefSeq", folder=databases_refseq_folders, assembly=assembly)
        log.debug("refSeq: "+str(refseq_file))

        # refSeqLink
        refseqlink_file = find_file_prefix(input_file=databases_refseqlink, prefix="ncbiRefSeqLink", folder=databases_refseq_folders, assembly=assembly)
        log.debug("refSeqLink: "+str(refseqlink_file))

        # Threads
        if not threads:
            threads = self.get_threads()
        log.debug("Threads: "+str(threads))

        # Variables
        table_variants = self.get_table_variants(clause="update")

        # Get variants SNV and InDel only
        query_variants = f"""
            SELECT "#CHROM" AS CHROM, POS, REF, ALT
            FROM {table_variants}
            WHERE REF ~ '^[A-Za-z]+$' AND ALT ~ '^[A-Za-z]+$'
            """
        df_variants = self.get_query_to_df(query_variants)

        # Added columns
        added_columns = []

        # Add hgvs column in variants table
        hgvs_column_name = "hgvs_" + str(random.randrange(1000))
        added_column = self.add_column(table_variants, hgvs_column_name, "STRING", default_value=None)
        added_columns.append(added_column)

        log.debug(f"refSeq loading...")
        # refSeq in duckDB
        refseq_table = get_refseq_table(conn=self.conn, refseq_table="refseq", refseq_file=refseq_file)
        # Loading all refSeq in Dataframe
        refseq_query = f"""
            SELECT df_variants.CHROM, df_variants.POS, {refseq_table}.name AS transcript
            FROM {refseq_table}
            JOIN df_variants ON (
                {refseq_table}.chrom = df_variants.CHROM
                AND {refseq_table}.txStart<=df_variants.POS
                AND {refseq_table}.txEnd>=df_variants.POS
            )
        """
        refseq_df = self.conn.query(refseq_query).pl()


        if refseqlink_file:
            log.debug(f"refSeqLink loading...")
            # refSeqLink in duckDB
            refseqlink_table = get_refseq_table(conn=self.conn, refseq_table="refseqlink", refseq_file=refseqlink_file)
            # Loading all refSeqLink in Dataframe
            protacc_column = "protAcc_with_ver"
            mrnaacc_column = "mrnaAcc_with_ver"
            refseqlink_query = f"""
                SELECT {refseq_table}.chrom, {protacc_column} AS protein, {mrnaacc_column} AS transcript
                FROM {refseqlink_table} 
                JOIN {refseq_table} ON ({refseq_table}.name = {refseqlink_table}.mrnaAcc_with_ver)
                WHERE protAcc_without_ver IS NOT NULL
            """
            # Polars Dataframe
            refseqlink_df = self.conn.query(f"{refseqlink_query}").pl()


        # Read RefSeq transcripts into a python dict/model.
        log.debug(f"Transcripts loading...")
        with tempfile.TemporaryDirectory() as tmpdir:
            transcripts_query = f"""
                COPY (
                    SELECT {refseq_table}.*
                    FROM {refseq_table}
                    JOIN df_variants ON (
                        {refseq_table}.chrom=df_variants.CHROM
                        AND {refseq_table}.txStart<=df_variants.POS
                        AND {refseq_table}.txEnd>=df_variants.POS
                    )
                )
                TO '{tmpdir}/transcript.tsv' (DELIMITER '\t');
            """
            self.conn.query(transcripts_query)
            with open(f'{tmpdir}/transcript.tsv') as infile:
                transcripts = read_transcripts(infile)

        # Polars connexion
        polars_conn = pl.SQLContext(register_globals=True, eager_execution=True)

        log.debug("Genome loading...")
        # Read genome sequence using pyfaidx.
        genome = Fasta(genome_file)

        log.debug("Start annotation HGVS...")

        # Create 
        # a Dask Dataframe from Pandas dataframe with partition as number of threads
        ddf = dd.from_pandas(df_variants, npartitions=threads)
        
        # Use dask.dataframe.apply() to apply function on each partition
        ddf[hgvs_column_name] = ddf.map_partitions(partition_function)

        # Convert Dask DataFrame to Pandas Dataframe
        df = ddf.compute()
        
        # Convert Pandas dataframe to parquet (due to error in cast VARCHAR -> NULL ???)
        with tempfile.TemporaryDirectory() as tmpdir:
            df_parquet = os.path.join(tmpdir,"df.parquet")
            df.to_parquet(df_parquet)

            # Update hgvs column
            update_variant_query = f"""
                UPDATE {table_variants}
                SET "{hgvs_column_name}"=df."{hgvs_column_name}"
                FROM read_parquet('{df_parquet}') as df
                WHERE variants."#CHROM" = df.CHROM
                AND variants.POS = df.POS
                AND variants.REF = df.REF
                AND variants.ALT = df.ALT
                AND df."{hgvs_column_name}" NOT IN ('') AND df."{hgvs_column_name}" NOT NULL
                """
            self.execute_query(update_variant_query)

        # Update INFO column
        sql_query_update = f"""
            UPDATE {table_variants}
            SET INFO = 
                concat(
                    CASE 
                        WHEN INFO NOT IN ('','.')
                        THEN concat(INFO, ';')
                        ELSE ''
                    END,
                    'hgvs=',
                    {hgvs_column_name}
                )
            WHERE "{hgvs_column_name}" NOT IN ('') AND "{hgvs_column_name}" NOT NULL
            """
        self.execute_query(sql_query_update)

        # Remove added columns
        for added_column in added_columns:
            self.drop_column(column=added_column)


    ###
    # Calculation
    ###

    def get_operations_help(self, operations_config_dict:dict = {}, operations_config_file:str = None) -> list:

        # Init
        operations_help = []

        # operations
        operations = self.get_operations_config(operations_config_dict=operations_config_dict, operations_config_file=operations_config_file)
        for op in operations:
            op_name = operations[op].get("name", op).upper()
            op_description = operations[op].get("description", op_name)
            op_available = operations[op].get("available", False)
            if op_available:
                operations_help.append(f"   {op_name}: {op_description}")

        # Sort operations
        operations_help.sort()

        # insert header
        operations_help.insert(0, "Available calculation operations:")

        # Return
        return operations_help


    def get_operations_config(self, operations_config_dict:dict = {}, operations_config_file:str = None) -> dict:

        operations_config_default = {
            "variant_chr_pos_alt_ref":
                {
                    "type": "sql",
                    "name": "variant_chr_pos_alt_ref",
                    "description": "Create a variant ID with chromosome, position, alt and ref",
                    "available": False,
                    "output_column_name": "variant_chr_pos_alt_ref",
                    "output_column_type": "String",
                    "output_column_description": "variant ID with chromosome, position, alt and ref",
                    "operation_query": """ concat("#CHROM", '_', "POS", '_', "REF", '_', "ALT") """,
                    "operation_info": True,
                },
            "VARTYPE":
                {
                    "type": "sql",
                    "name": "VARTYPE",
                    "description": "Variant type (e.g. SNV, INDEL, MNV, BND...)",
                    "available": True,
                    "output_column_name": "VARTYPE",
                    "output_column_type": "String",
                    "output_column_description": "Variant type: SNV if X>Y, MOSAIC if X>Y,Z or X,Y>Z, INDEL if XY>Z or X>YZ",
                    "operation_query": """
                        CASE
                            WHEN "SVTYPE" NOT NULL THEN "SVTYPE"
                            WHEN LENGTH(REF) = 1 AND LENGTH(ALT) = 1 THEN 'SNV'
                            WHEN REF LIKE '%,%' OR ALT LIKE '%,%' THEN 'MOSAIC'
                            WHEN LENGTH(REF) == LENGTH(ALT) AND LENGTH(REF) > 1 THEN 'MNV'
                            WHEN LENGTH(REF) <> LENGTH(ALT) THEN 'INDEL'
                            ELSE 'UNDEFINED'
                        END
                        """,
                    "info_fields": ["SVTYPE"],
                    "operation_info": True,
                },
            "snpeff_hgvs":
                {
                    "type": "python",
                    "name": "snpeff_hgvs",
                    "description": "HGVS nomenclatures from snpEff annotation",
                    "available": True,
                    "function_name": "calculation_extract_snpeff_hgvs",
                    "function_params": []
                },
            "NOMEN":
                {
                    "type": "python",
                    "name": "NOMEN",
                    "description": "NOMEN information (e.g. NOMEN, CNOMEN, PNOMEN...) from HGVS nomenclature field",
                    "available": True,
                    "function_name": "calculation_extract_nomen",
                    "function_params": []
                },
            "FINDBYPIPELINE":
                {
                    "type": "python",
                    "name": "FINDBYPIPELINE",
                    "description": "Number of pipeline that identify the variant (for multi pipeline VCF)",
                    "available": True,
                    "function_name": "calculation_find_by_pipeline",
                    "function_params": ["findbypipeline"]
                },
            "FINDBYSAMPLE":
                {
                    "type": "python",
                    "name": "FINDBYSAMPLE",
                    "description": "Number of sample that have a genotype for the variant (for multi sample VCF)",
                    "available": True,
                    "function_name": "calculation_find_by_pipeline",
                    "function_params": ["findbysample"]
                },
            "GENOTYPECONCORDANCE":
                {
                    "type": "python",
                    "name": "GENOTYPECONCORDANCE",
                    "description": "Concordance of genotype for multi caller VCF",
                    "available": True,
                    "function_name": "calculation_genotype_concordance",
                    "function_params": []
                },
            "BARCODE":
                {
                    "type": "python",
                    "name": "BARCODE",
                    "description": "BARCODE as VaRank tool",
                    "available": True,
                    "function_name": "calculation_barcode",
                    "function_params": []
                },
            "TRIO":
                {
                    "type": "python",
                    "name": "TRIO",
                    "description": "Inheritance for a trio family",
                    "available": True,
                    "function_name": "calculation_trio",
                    "function_params": []
                },
            "VAF":
                {
                    "type": "python",
                    "name": "VAF",
                    "description": "Variant Allele Frequency (VAF) harmonization",
                    "available": True,
                    "function_name": "calculation_vaf_normalization",
                    "function_params": []
                },
            "VAF_stats":
                {
                    "type": "python",
                    "name": "VAF_stats",
                    "description": "Variant Allele Frequency (VAF) statistics",
                    "available": True,
                    "function_name": "calculation_genotype_stats",
                    "function_params": ["VAF"]
                },
            "DP_stats":
                {
                    "type": "python",
                    "name": "DP_stats",
                    "description": "Depth (DP) statistics",
                    "available": True,
                    "function_name": "calculation_genotype_stats",
                    "function_params": ["DP"]
                },
            "variant_id":
                {
                    "type": "python",
                    "name": "variant_id",
                    "description": "Variant ID generated from variant position and type",
                    "available": True,
                    "function_name": "calculation_variant_id",
                    "function_params": []
                },
        }

        # Create with default operations
        operations_config = operations_config_default

        # Replace operations from dict
        for operation_config in operations_config_dict:
            operations_config[operation_config] = operations_config_dict[operation_config]

        # Replace operations from file
        if operations_config_file:
            if os.path.exists(operations_config_file):
                with open(operations_config_file) as operations_config_file_content:
                    operations_config_file_dict = json.load(operations_config_file_content)
                for operation_config in operations_config_file_dict:
                    operations_config[operation_config] = operations_config_file_dict[operation_config]
            else:
                log.error(f"Operations config file '{operations_config_file}' does NOT exist")
                raise ValueError(f"Operations config file '{operations_config_file}' does NOT exist")

        return operations_config


    def calculation(self, operations:dict = None, operations_config_dict:dict = {}, operations_config_file:str = None) -> None:
        """
        It takes a list of operations, and for each operation, it checks if it's a python or sql
        operation, and then calls the appropriate function

        param json example:
            "calculation": {
                "NOMEN": {
                    "options": {
                        "hgvs_field": "hgvs"
                    },
                "middle" : null
            }
        """

        # operations config
        operations_config = self.get_operations_config(operations_config_dict=operations_config_dict, operations_config_file=operations_config_file)

        # Upper keys
        operations_config = {k.upper(): v for k, v in operations_config.items()}
        
        # Operations for calculation
        if not operations:
            operations = self.get_param().get("calculation", {})

        # For each operations
        for operation_name in operations:
            operation_name = operation_name.upper()
            if operation_name in operations_config:
                log.info(f"Calculation '{operation_name}'")
                operation = operations_config[operation_name]
                operation_type = operation.get("type", "sql")
                if operation_type == "python":
                    self.calculation_process_function(operation=operation, operation_name=operation_name)
                elif operation_type == "sql":
                    self.calculation_process_sql(operation=operation, operation_name=operation_name)
                else:
                    log.error(f"Operations config: Type '{operation_type}' NOT available")
                    raise ValueError(f"Operations config: Type '{operation_type}' NOT available")
            else:
                log.error(f"Operations config: Calculation '{operation_name}' NOT available")
                raise ValueError(f"Operations config: Calculation '{operation_name}' NOT available")

        # Explode INFOS fields into table fields
        if self.get_explode_infos():
            self.explode_infos(prefix=self.get_explode_infos_prefix(), fields=self.get_explode_infos_fields(), force=True)


    def calculation_process_sql(self, operation:dict, operation_name:str = "unknown") -> None:
        """
        The `calculation_process_sql` function takes in a mathematical operation as a string and
        performs the operation, updating the specified table with the result.
        
        :param operation: The `operation` parameter is a dictionary that contains information about the
        mathematical operation to be performed. It includes the following keys:
        :type operation: dict
        :param operation_name: The `operation_name` parameter is a string that represents the name of
        the mathematical operation being performed. It is used for logging and error handling purposes,
        defaults to unknown
        :type operation_name: str (optional)
        """

        # table variants
        table_variants = self.get_table_variants(clause="alter")

        # Operation infos
        operation_name = operation.get('name', 'unknown')
        log.debug(f"process sql {operation_name}")
        output_column_name = operation.get('output_column_name',operation_name)
        output_column_type = operation.get('output_column_type', 'String')
        prefix = operation.get('explode_infos_prefix', '')
        output_column_type_sql = code_type_map_to_sql.get(
            output_column_type, "VARCHAR")
        output_column_description = operation.get('output_column_description', f'{operation_name} operation')
        operation_query = operation.get('operation_query', None)
        if isinstance(operation_query, list):
            operation_query = " ".join(operation_query)
        operation_info_fields = operation.get('info_fields', [])
        operation_info_fields_check = operation.get('info_fields_check', False)
        operation_info = operation.get('operation_info', True)

        if operation_query:

            # Info fields check
            operation_info_fields_check_result = True
            if operation_info_fields_check:
                header_infos = self.get_header().infos
                for info_field in operation_info_fields:
                    operation_info_fields_check_result = operation_info_fields_check_result and info_field in header_infos

            # If info fields available 
            if operation_info_fields_check_result:

                # Added_columns
                added_columns = []

                # Create VCF header field
                vcf_reader = self.get_header()
                vcf_reader.infos[output_column_name] = vcf.parser._Info(
                    output_column_name,
                    ".",
                    output_column_type,
                    output_column_description,
                    "howard calculation",
                    "0",
                    self.code_type_map.get(output_column_type)
                )

                # Explode infos if needed
                log.debug(f"calculation_process_sql prefix {prefix}")
                added_columns += self.explode_infos(prefix=prefix, fields=[output_column_name] + operation_info_fields, force=True)

                # Create column
                added_column = self.add_column(table_name=table_variants, column_name=prefix+output_column_name, column_type=output_column_type_sql, default_value="null")
                added_columns.append(added_column)

                # Operation calculation
                try:

                    # Query to update calculation column
                    sql_update = f"""
                        UPDATE {table_variants}
                        SET "{prefix}{output_column_name}" = ({operation_query})
                    """
                    self.conn.execute(sql_update)

                    # Add to INFO
                    if operation_info:
                        sql_update_info = f"""
                            UPDATE {table_variants}
                            SET "INFO" =
                                concat(
                                    CASE
                                        WHEN "INFO" IS NOT NULL
                                        THEN concat("INFO", ';')
                                        ELSE ''
                                    END,
                                    '{output_column_name}=',
                                    "{prefix}{output_column_name}"
                                )
                            WHERE "{prefix}{output_column_name}" IS NOT NULL AND "{prefix}{output_column_name}" NOT IN ('')
                        """
                        self.conn.execute(sql_update_info)

                except:
                    log.error(f"Operations config: Calculation '{operation_name}' query failed")
                    raise ValueError(f"Operations config: Calculation '{operation_name}' query failed")

                # Remove added columns
                for added_column in added_columns:
                    log.debug(f"added_column: {added_column}")
                    self.drop_column(column=added_column)

            else:
                log.error(f"Operations config: Calculation '{operation_name}' DOES NOT contain all mandatory fields {operation_info_fields}")
                raise ValueError(f"Operations config: Calculation '{operation_name}' DOES NOT contain all mandatory fields {operation_info_fields}")

        else:
            log.error(f"Operations config: Calculation '{operation_name}' query NOT defined")
            raise ValueError(f"Operations config: Calculation '{operation_name}' query NOT defined")
        

    def calculation_process_function(self, operation:dict, operation_name:str = "unknown") -> None:
        """
        The `calculation_process_function` takes in an operation dictionary and performs the specified
        function with the given parameters.
        
        :param operation: The `operation` parameter is a dictionary that contains information about the
        operation to be performed. It has the following keys:
        :type operation: dict
        :param operation_name: The `operation_name` parameter is a string that represents the name of
        the operation being performed. It is used for logging purposes, defaults to unknown
        :type operation_name: str (optional)
        """
        
        operation_name = operation['name']
        log.debug(f"process sql {operation_name}")
        function_name = operation['function_name']
        function_params = operation['function_params']
        getattr(self, function_name)(*function_params)


    def calculation_variant_id(self) -> None:
        """
        The function `calculation_variant_id` adds a variant ID annotation to a VCF file header and
        updates the INFO field of a variants table with the variant ID.
        """

        # variant_id annotation field
        variant_id_tag = self.get_variant_id_column()

        # variant_id hgvs tags"
        vcf_infos_tags = {
            variant_id_tag: "howard variant ID annotation",
        }

        # Variants table
        table_variants = self.get_table_variants()
        
        # Header
        vcf_reader = self.get_header()

        # Add variant_id to header
        vcf_reader.infos[variant_id_tag] = vcf.parser._Info(
                    variant_id_tag,
                    ".",
                    "String",
                    vcf_infos_tags.get(
                        variant_id_tag, "howard variant ID annotation"),
                    "howard calculation",
                    "0",
                    self.code_type_map.get("String")
                )
        
        # Update
        sql_update = f"""
            UPDATE {table_variants}
            SET "INFO" = 
                concat(
                    CASE
                        WHEN "INFO" IS NULL OR "INFO" IN ('','.')
                        THEN ''
                        ELSE concat("INFO", ';')
                    END,
                    '{variant_id_tag}=',
                    "{variant_id_tag}"
                )
        """
        self.conn.execute(sql_update)


    def calculation_extract_snpeff_hgvs(self) -> None:
        """
        The function `calculation_extract_snpeff_hgvs` extracts HGVS nomenclatures from the SnpEff
        annotation field in a VCF file and adds them as a new column in the variants table.
        """

        # SnpEff annotation field
        snpeff_ann = "ANN"

        # SnpEff annotation field
        snpeff_hgvs = "snpeff_hgvs"

        # Snpeff hgvs tags
        vcf_infos_tags = {
            snpeff_hgvs: "HGVS nomenclatures from snpEff annotation",
        }

        # Prefix
        prefix = self.get_explode_infos_prefix()
        if prefix:
            prefix = "INFO/"

        # snpEff fields
        speff_ann_infos = prefix+snpeff_ann
        speff_hgvs_infos = prefix+snpeff_hgvs

        # Variants table
        table_variants = self.get_table_variants()
        
        # Header
        vcf_reader = self.get_header()

        # Add columns
        added_columns = []

        # Explode HGVS field in column
        added_columns += self.explode_infos(fields=[snpeff_ann])

        if "ANN" in vcf_reader.infos:
                
            log.debug(vcf_reader.infos["ANN"])

            # Create variant id
            variant_id_column = self.get_variant_id_column()

            # Create dataframe
            dataframe_snpeff_hgvs = self.get_query_to_df(
                f""" SELECT "{variant_id_column}", "{speff_ann_infos}" FROM {table_variants} """)
            
            # Create main NOMEN column
            dataframe_snpeff_hgvs[speff_hgvs_infos] = dataframe_snpeff_hgvs[speff_ann_infos].apply(
                    lambda x: extract_snpeff_hgvs(str(x)))

            # Add snpeff_hgvs to header
            vcf_reader.infos[snpeff_hgvs] = vcf.parser._Info(
                        snpeff_hgvs,
                        ".",
                        "String",
                        vcf_infos_tags.get(
                            snpeff_hgvs, "snpEff hgvs annotations"),
                        "howard calculation",
                        "0",
                        self.code_type_map.get("String")
                    )
            
            # Update
            sql_update = f"""
                UPDATE variants
                SET "INFO" = 
                    concat(
                        CASE
                            WHEN "INFO" IS NULL OR "INFO" IN ('','.')
                            THEN ''
                            ELSE concat("INFO", ';')
                        END,
                        CASE 
                            WHEN dataframe_snpeff_hgvs."{speff_hgvs_infos}" NOT IN ('','.','NaN')
                            AND dataframe_snpeff_hgvs."{speff_hgvs_infos}" NOT NULL
                            THEN concat(
                                    '{snpeff_hgvs}=',
                                    dataframe_snpeff_hgvs."{speff_hgvs_infos}"
                                )
                            ELSE ''
                        END
                    )
                FROM dataframe_snpeff_hgvs
                WHERE {table_variants}."{variant_id_column}" = dataframe_snpeff_hgvs."{variant_id_column}"

            """
            self.conn.execute(sql_update)

            # Delete dataframe
            del dataframe_snpeff_hgvs
            gc.collect()
        
        else:

            log.warning("No snpEff annotation. Please Anotate with snpEff before use this calculation option")

        # Remove added columns
        for added_column in added_columns:
            self.drop_column(column=added_column)


    def calculation_extract_nomen(self) -> None:
        """
        This function extracts the HGVS nomenclature from the calculation/identification of NOMEN.
        """

        # NOMEN field
        field_nomen_dict = "NOMEN_DICT"

        # NOMEN structure
        nomen_dict = {
            "NOMEN": "NOMEN hgvs nomenclature considered as reference hgvs (official transcript, first otherwise)",
            "CNOMEN": "CNOMEN hgvs nomenclature at DNA level related to a transcript (TNOMEN)",
            "RNOMEN": "RNOMEN hgvs nomenclature at RNA level related to a transcript (TNOMEN)",
            "NNOMEN": "NNOMEN hgvs nomenclature for non-coding variant",
            "PNOMEN": "PNOMEN hgvs nomenclature at Protein level related to a transcript (TNOMEN)",
            "TVNOMEN": "TVNOMEN hgvs transcript with version (if any) used (e.g. for CNOMEN and PNOMEN)",
            "TNOMEN": "TNOMEN hgvs transcript used (e.g. for CNOMEN and PNOMEN)",
            "VNOMEN": "VNOMEN hgvs transcript version used (e.g. for CNOMEN and PNOMEN)",
            "ENOMEN": "ENOMEN hgvs exon nomenclature related to a transcript (TNOMEN)",
            "GNOMEN": "GNOMEN hgvs gene nomenclature related to a transcript (TNOMEN)",
        }

        # Param
        param = self.get_param()

        # Prefix
        prefix = self.get_explode_infos_prefix()

        # Header
        vcf_reader = self.get_header()

        # Get HGVS field
        hgvs_field = param.get("calculation", {}).get(
            "NOMEN", {}).get("options", {}).get("hgvs_field", "hgvs")
        
        # Get transcripts
        transcripts_file = param.get("calculation", {}).get(
            "NOMEN", {}).get("options", {}).get("transcripts", None)
        log.debug(f"Transcript file '{transcripts_file}'")
        transcripts = []
        if transcripts_file:
            if os.path.exists(transcripts_file):
                log.debug(f"Transcript file '{transcripts_file}' does exist")
                transcripts_dataframe = pd.read_csv(transcripts_file, sep="\t", header=None, names=['transcript', 'gene'])
                transcripts = transcripts_dataframe.iloc[:, 0].tolist()
                log.debug(f"Transcripts DF")
                log.debug(transcripts_dataframe)
                log.debug(f"Transcripts: {transcripts}")
            else:
                log.error(f"Transcript file '{transcripts_file}' does NOT exist")
                raise ValueError(f"Transcript file '{transcripts_file}' does NOT exist")
            
        # Added columns
        added_columns = []

        # Explode HGVS field in column
        added_columns += self.explode_infos(fields=[hgvs_field])

        # extra infos
        extra_infos = self.get_extra_infos()
        extra_field = prefix+hgvs_field

        if extra_field in extra_infos:

            # Create dataframe
            dataframe_hgvs = self.get_query_to_df(
                f""" SELECT "#CHROM", "POS", "REF", "ALT", "{extra_field}" FROM variants """)

            # Create main NOMEN column
            dataframe_hgvs[field_nomen_dict] = dataframe_hgvs[extra_field].apply(
                lambda x: find_nomen(str(x), transcripts=transcripts))

            # Explode NOMEN Structure and create SQL set for update
            sql_nomen_fields = []
            for nomen_field in nomen_dict:

                # Explode each field into a column
                dataframe_hgvs[nomen_field] = dataframe_hgvs[field_nomen_dict].apply(
                    lambda x: dict(x).get(nomen_field, ""))
                
                # Create VCF header field
                vcf_reader.infos[nomen_field] = vcf.parser._Info(
                    nomen_field,
                    ".",
                    "String",
                    nomen_dict.get(
                        nomen_field, "howard calculation NOMEN"),
                    "howard calculation",
                    "0",
                    self.code_type_map.get("String")
                )
                sql_nomen_fields.append(
                    f"""
                        CASE 
                            WHEN dataframe_hgvs."{nomen_field}" NOT NULL AND dataframe_hgvs."{nomen_field}" NOT IN ('')
                            THEN concat(
                                    ';{nomen_field}=',
                                    dataframe_hgvs."{nomen_field}"
                                )
                            ELSE ''
                        END
                    """)

            # SQL set for update
            sql_nomen_fields_set = ", ".join(sql_nomen_fields)

            # Update
            sql_update = f"""
                UPDATE variants
                SET "INFO" = 
                    concat(
                        CASE
                            WHEN "INFO" IS NULL
                            THEN ''
                            ELSE "INFO"
                        END,
                        {sql_nomen_fields_set}
                    )
                FROM dataframe_hgvs
                WHERE variants."#CHROM" = dataframe_hgvs."#CHROM"
                    AND variants."POS" = dataframe_hgvs."POS" 
                    AND variants."REF" = dataframe_hgvs."REF"
                    AND variants."ALT" = dataframe_hgvs."ALT"
            """
            self.conn.execute(sql_update)

            # Delete dataframe
            del dataframe_hgvs
            gc.collect()

        # Remove added columns
        for added_column in added_columns:
            self.drop_column(column=added_column)


    def calculation_find_by_pipeline(self, tag:str = "findbypipeline") -> None:
        """
        The function `calculation_find_by_pipeline` performs a calculation to find the number of
        pipeline/sample for a variant and updates the variant information in a VCF file.
        
        :param tag: The `tag` parameter is a string that represents the annotation field for the
        "findbypipeline" information in the VCF file. It is used to create the annotation field in the
        VCF header and to update the corresponding field in the variants table, defaults to
        findbypipeline
        :type tag: str (optional)
        """

        # if FORMAT and samples
        if "FORMAT" in self.get_header_columns_as_list() and self.get_header_sample_list():

            # findbypipeline annotation field
            findbypipeline_tag = tag
            
            # VCF infos tags
            vcf_infos_tags = {
                findbypipeline_tag: f"Number of pipeline/sample for a variant ({findbypipeline_tag})",
            }

            # Prefix
            prefix = self.get_explode_infos_prefix()

            # Field
            findbypipeline_infos = prefix+findbypipeline_tag

            # Variants table
            table_variants = self.get_table_variants()

            # Header
            vcf_reader = self.get_header()

            # Create variant id
            variant_id_column = self.get_variant_id_column()

            # variant_id, FORMAT and samples
            samples_fields = f" {variant_id_column}, FORMAT , " + \
                    " , ".join(self.get_header_sample_list())
            
            # Create dataframe
            dataframe_findbypipeline = self.get_query_to_df(
                f""" SELECT {samples_fields} FROM {table_variants} """)
        
            # Create findbypipeline column
            dataframe_findbypipeline[findbypipeline_infos] = dataframe_findbypipeline.apply(lambda row: findbypipeline(row, samples=self.get_header_sample_list()), axis=1)
            
            # Add snpeff_hgvs to header
            vcf_reader.infos[findbypipeline_tag] = vcf.parser._Info(
                        findbypipeline_tag,
                        ".",
                        "String",
                        vcf_infos_tags.get(
                            findbypipeline_tag, "Find in pipeline/sample"),
                        "howard calculation",
                        "0",
                        self.code_type_map.get("String")
                    )

            # Update
            sql_update = f"""
                UPDATE variants
                SET "INFO" = 
                    concat(
                        CASE
                            WHEN "INFO" IS NULL OR "INFO" IN ('','.')
                            THEN ''
                            ELSE concat("INFO", ';')
                        END,
                        CASE 
                            WHEN dataframe_findbypipeline."{findbypipeline_infos}" NOT IN ('','.')
                                AND dataframe_findbypipeline."{findbypipeline_infos}" NOT NULL
                            THEN concat(
                                    '{findbypipeline_tag}=',
                                    dataframe_findbypipeline."{findbypipeline_infos}"
                                )
                            ELSE ''
                        END
                    )
                FROM dataframe_findbypipeline
                WHERE variants."{variant_id_column}" = dataframe_findbypipeline."{variant_id_column}"
            """
            self.conn.execute(sql_update)

            # Delete dataframe
            del dataframe_findbypipeline
            gc.collect()


    def calculation_genotype_concordance(self) -> None:
        """
        The function `calculation_genotype_concordance` calculates the genotype concordance for
        multi-caller VCF files and updates the variant information in the database.
        """

        # if FORMAT and samples
        if "FORMAT" in self.get_header_columns_as_list() and self.get_header_sample_list():

            # genotypeconcordance annotation field
            genotypeconcordance_tag = "genotypeconcordance"
            
            # VCF infos tags
            vcf_infos_tags = {
                genotypeconcordance_tag: "Concordance of genotype for multi caller VCF",
            }

            # Prefix
            prefix = self.get_explode_infos_prefix()

            # Field
            genotypeconcordance_infos = prefix+genotypeconcordance_tag

            # Variants table
            table_variants = self.get_table_variants()

            # Header
            vcf_reader = self.get_header()

            # Create variant id
            variant_id_column = self.get_variant_id_column()

            # variant_id, FORMAT and samples
            samples_fields = f" {variant_id_column}, FORMAT , " + \
                    " , ".join(self.get_header_sample_list())
            
            # Create dataframe
            dataframe_genotypeconcordance = self.get_query_to_df(
                f""" SELECT {samples_fields} FROM {table_variants} """)

        
            # Create genotypeconcordance column
            dataframe_genotypeconcordance[genotypeconcordance_infos] = dataframe_genotypeconcordance.apply(lambda row: genotypeconcordance(row, samples=self.get_header_sample_list()), axis=1)
            
            # Add genotypeconcordance to header
            vcf_reader.infos[genotypeconcordance_tag] = vcf.parser._Info(
                        genotypeconcordance_tag,
                        ".",
                        "String",
                        vcf_infos_tags.get(
                            genotypeconcordance_tag, "snpEff hgvs annotations"),
                        "howard calculation",
                        "0",
                        self.code_type_map.get("String")
                    )

            # Update
            sql_update = f"""
                UPDATE variants
                SET "INFO" = 
                    concat(
                        CASE
                            WHEN "INFO" IS NULL OR "INFO" IN ('','.')
                            THEN ''
                            ELSE concat("INFO", ';')
                        END,
                        CASE
                            WHEN dataframe_genotypeconcordance."{genotypeconcordance_infos}" NOT IN ('','.')
                                AND dataframe_genotypeconcordance."{genotypeconcordance_infos}" NOT NULL
                            THEN concat(
                                    '{genotypeconcordance_tag}=',
                                    dataframe_genotypeconcordance."{genotypeconcordance_infos}"
                                )
                            ELSE ''
                        END
                    )
                FROM dataframe_genotypeconcordance
                WHERE variants."{variant_id_column}" = dataframe_genotypeconcordance."{variant_id_column}"
            """
            self.conn.execute(sql_update)

            # Delete dataframe
            del dataframe_genotypeconcordance
            gc.collect()


    def calculation_barcode(self) -> None:
        """
        The `calculation_barcode` function calculates barcode values for variants in a VCF file and
        updates the INFO field in the file with the calculated barcode values.
        """

        # if FORMAT and samples
        if "FORMAT" in self.get_header_columns_as_list() and self.get_header_sample_list():

            # barcode annotation field
            barcode_tag = "barcode"
            
            # VCF infos tags
            vcf_infos_tags = {
                "barcode": "barcode calculation (VaRank)",
            }

            # Prefix
            prefix = self.get_explode_infos_prefix()

            # Field
            barcode_infos = prefix+barcode_tag

            # Variants table
            table_variants = self.get_table_variants()

            # Header
            vcf_reader = self.get_header()

            # Create variant id
            variant_id_column = self.get_variant_id_column()

            # variant_id, FORMAT and samples
            samples_fields = f" {variant_id_column}, FORMAT , " + \
                    " , ".join(self.get_header_sample_list())
            
            # Create dataframe
            dataframe_barcode = self.get_query_to_df(
                f""" SELECT {samples_fields} FROM {table_variants} """)

            # Create barcode column
            dataframe_barcode[barcode_infos] = dataframe_barcode.apply(lambda row: barcode(row, samples=self.get_header_sample_list()), axis=1)
            
            # Add barcode to header
            vcf_reader.infos[barcode_tag] = vcf.parser._Info(
                        barcode_tag,
                        ".",
                        "String",
                        vcf_infos_tags.get(
                            barcode_tag, "snpEff hgvs annotations"),
                        "howard calculation",
                        "0",
                        self.code_type_map.get("String")
                    )

            # Update
            sql_update = f"""
                UPDATE {table_variants}
                SET "INFO" = 
                    concat(
                        CASE
                            WHEN "INFO" IS NULL OR "INFO" IN ('','.')
                            THEN ''
                            ELSE concat("INFO", ';')
                        END,
                        CASE
                            WHEN dataframe_barcode."{barcode_infos}" NOT IN ('','.')
                            AND dataframe_barcode."{barcode_infos}" NOT NULL
                            THEN concat(
                                    '{barcode_tag}=',
                                    dataframe_barcode."{barcode_infos}"
                                )
                            ELSE ''
                        END
                    )
                FROM dataframe_barcode
                WHERE {table_variants}."{variant_id_column}" = dataframe_barcode."{variant_id_column}"
            """
            self.conn.execute(sql_update)

            # Delete dataframe
            del dataframe_barcode
            gc.collect()


    def calculation_trio(self) -> None:
        """
        The `calculation_trio` function performs trio calculations on a VCF file by adding trio
        information to the INFO field of each variant.
        """

        # if FORMAT and samples
        if "FORMAT" in self.get_header_columns_as_list() and self.get_header_sample_list():

            # trio annotation field
            trio_tag = "trio"
            
            # VCF infos tags
            vcf_infos_tags = {
                "trio": "trio calculation",
            }

            # Param
            param = self.get_param()

            # Prefix
            prefix = self.get_explode_infos_prefix()

            # Trio param
            trio_ped = param.get("calculation", {}).get("TRIO", {})
            if trio_ped:
                trio_samples = [
                    trio_ped.get("father",""),
                    trio_ped.get("mother",""),
                    trio_ped.get("child","")
                ]
            else:
                trio_samples = self.get_header_sample_list()[0:3]

            log.debug(f"Param for trio sample: {trio_ped}")
            log.debug(f"List of trio sample: {trio_samples}")

            # Field
            trio_infos = prefix+trio_tag

            # Variants table
            table_variants = self.get_table_variants()

            # Header
            vcf_reader = self.get_header()

            # Create variant id
            variant_id_column = self.get_variant_id_column()

            # variant_id, FORMAT and samples
            samples_fields = f" {variant_id_column}, FORMAT , " + \
                    " , ".join(self.get_header_sample_list())
            
            # Create dataframe
            dataframe_trio = self.get_query_to_df(
                f""" SELECT {samples_fields} FROM {table_variants} """)

        
            # Create trio column
            dataframe_trio[trio_infos] = dataframe_trio.apply(lambda row: trio(row, samples=trio_samples), axis=1)
            
            # Add trio to header
            vcf_reader.infos[trio_tag] = vcf.parser._Info(
                        trio_tag,
                        ".",
                        "String",
                        vcf_infos_tags.get(
                            trio_tag, "snpEff hgvs annotations"),
                        "howard calculation",
                        "0",
                        self.code_type_map.get("String")
                    )

            # Update
            sql_update = f"""
                UPDATE {table_variants}
                SET "INFO" = 
                    concat(
                        CASE
                            WHEN "INFO" IS NULL OR "INFO" IN ('','.')
                            THEN ''
                            ELSE concat("INFO", ';')
                        END,
                        CASE
                            WHEN dataframe_trio."{trio_infos}" NOT IN ('','.')
                             AND dataframe_trio."{trio_infos}" NOT NULL
                            THEN concat(
                                    '{trio_tag}=',
                                    dataframe_trio."{trio_infos}"
                                )
                            ELSE ''
                        END
                    )
                FROM dataframe_trio
                WHERE {table_variants}."{variant_id_column}" = dataframe_trio."{variant_id_column}"
            """
            self.conn.execute(sql_update)

            # Delete dataframe
            del dataframe_trio
            gc.collect()


    def calculation_vaf_normalization(self) -> None:
        """
        The `calculation_vaf_normalization` function calculates the VAF (Variant Allele Frequency)
        normalization for each sample in a VCF file and updates the FORMAT and INFO fields accordingly.
        :return: The function does not return anything.
        """

        # if FORMAT and samples
        if "FORMAT" in self.get_header_columns_as_list() and self.get_header_sample_list():

            # vaf_normalization annotation field
            vaf_normalization_tag = "VAF"
            
            # VCF infos tags
            vcf_infos_tags = {
                "VAF": "VAF Variant Frequency",
            }

            # Prefix
            prefix = self.get_explode_infos_prefix()

            # Variants table
            table_variants = self.get_table_variants()

            # Header
            vcf_reader = self.get_header()

            # Do not calculate if VAF already exists
            if "VAF" in vcf_reader.formats:
                log.debug("VAF already on genotypes")
                return

            # Create variant id
            variant_id_column = self.get_variant_id_column()

            # variant_id, FORMAT and samples
            samples_fields = f" {variant_id_column}, FORMAT , " + \
                    " , ".join(self.get_header_sample_list())
            
            # Create dataframe
            dataframe_vaf_normalization = self.get_query_to_df(
                f""" SELECT {variant_id_column}, FORMAT, {samples_fields} FROM {table_variants} """)

            vaf_normalization_set = []

            # for each sample vaf_normalization 
            for sample in self.get_header_sample_list():
                dataframe_vaf_normalization[sample] = dataframe_vaf_normalization.apply(lambda row: vaf_normalization(row, sample=sample), axis=1)
                vaf_normalization_set.append(f""" "{sample}" = dataframe_vaf_normalization."{sample}" """)

            # Add VAF to FORMAT
            dataframe_vaf_normalization["FORMAT"] = dataframe_vaf_normalization["FORMAT"].apply(lambda x: str(x)+":VAF")
            vaf_normalization_set.append(f""" "FORMAT" = dataframe_vaf_normalization."FORMAT" """)

            # Add vaf_normalization to header
            vcf_reader.formats[vaf_normalization_tag] = vcf.parser._Format(
                        id=vaf_normalization_tag,
                        num="1",
                        type="Float",
                        desc=vcf_infos_tags.get(
                            vaf_normalization_tag, "VAF Variant Frequency"),
                        type_code=self.code_type_map.get("Float")
                    )

            # Create fields to add in INFO
            sql_vaf_normalization_set = " , ".join(vaf_normalization_set)

            # Update
            sql_update = f"""
                UPDATE {table_variants}
                SET {sql_vaf_normalization_set}
                FROM dataframe_vaf_normalization
                WHERE variants."{variant_id_column}" = dataframe_vaf_normalization."{variant_id_column}"

            """
            self.conn.execute(sql_update)

            # Delete dataframe
            del dataframe_vaf_normalization
            gc.collect()


    def calculation_genotype_stats(self, info:str = "VAF") -> None:
        """
        The `calculation_genotype_stats` function calculates genotype statistics for a given information
        field in a VCF file and updates the INFO column of the variants table with the calculated
        statistics.
        
        :param info: The `info` parameter is a string that represents the type of information for which
        genotype statistics are calculated. It is used to generate various VCF info tags for the
        statistics, such as the number of occurrences, the list of values, the minimum value, the
        maximum value, the mean, the median, defaults to VAF
        :type info: str (optional)
        """

        # if FORMAT and samples
        if "FORMAT" in self.get_header_columns_as_list() and self.get_header_sample_list():

            # vaf_stats annotation field
            vaf_stats_tag = info+"_stats"
            
            # VCF infos tags
            vcf_infos_tags = {
                info+"_stats_nb": f"genotype {info} Statistics - number of {info}",
                info+"_stats_list": f"genotype {info} Statistics - list of {info}",
                info+"_stats_min": f"genotype {info} Statistics - min {info}",
                info+"_stats_max": f"genotype {info} Statistics - max {info}",
                info+"_stats_mean": f"genotype {info} Statistics - mean {info}",
                info+"_stats_mediane": f"genotype {info} Statistics - mediane {info}",
                info+"_stats_stdev": f"genotype {info} Statistics - standard deviation {info}",
            }
        
            # Prefix
            prefix = self.get_explode_infos_prefix()

            # Field
            vaf_stats_infos = prefix+vaf_stats_tag

            # Variants table
            table_variants = self.get_table_variants()

            # Header
            vcf_reader = self.get_header()

            # Create variant id
            variant_id_column = self.get_variant_id_column()

            # variant_id, FORMAT and samples
            samples_fields = f" {variant_id_column}, FORMAT , " + \
                    " , ".join(self.get_header_sample_list())
            
            # Create dataframe
            dataframe_vaf_stats = self.get_query_to_df(
                f""" SELECT {samples_fields} FROM {table_variants} """)

            # Create vaf_stats column
            dataframe_vaf_stats[vaf_stats_infos] = dataframe_vaf_stats.apply(lambda row: genotype_stats(row, samples=self.get_header_sample_list(), info=info), axis=1)
            
            # List of vcf tags
            sql_vaf_stats_fields = []

            # Check all VAF stats infos
            for stat in vcf_infos_tags:

                # Extract stats
                dataframe_vaf_stats[stat] = dataframe_vaf_stats[vaf_stats_infos].apply(lambda x: dict(x).get(stat, ""))

                # Add snpeff_hgvs to header
                vcf_reader.infos[stat] = vcf.parser._Info(
                            stat,
                            ".",
                            "String",
                            vcf_infos_tags.get(
                                stat, "genotype statistics"),
                            "howard calculation",
                            "0",
                            self.code_type_map.get("String")
                        )
                
                if len(sql_vaf_stats_fields):
                    sep = ";"
                else:
                    sep = ""

                # Create fields to add in INFO
                sql_vaf_stats_fields.append(
                    f"""
                        CASE
                            WHEN dataframe_vaf_stats."{stat}" NOT NULL
                            THEN concat(
                                    '{sep}{stat}=',
                                    dataframe_vaf_stats."{stat}"
                                )
                            ELSE ''
                        END
                    """)
            
            # SQL set for update
            sql_vaf_stats_fields_set = ",  ".join(sql_vaf_stats_fields)

            # Update
            sql_update = f"""
                UPDATE variants
                SET "INFO" = 
                    concat(
                        CASE
                            WHEN "INFO" IS NULL OR "INFO" IN ('','.')
                            THEN ''
                            ELSE concat("INFO", ';')
                        END,
                        {sql_vaf_stats_fields_set}
                    )
                FROM dataframe_vaf_stats
                WHERE variants."{variant_id_column}" = dataframe_vaf_stats."{variant_id_column}"

            """
            self.conn.execute(sql_update)

            # Delete dataframe
            del dataframe_vaf_stats
            gc.collect()


