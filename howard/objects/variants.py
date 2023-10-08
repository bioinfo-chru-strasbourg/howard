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


    def set_connexion(self, conn) -> None:
        """
        It creates a connection to the database

        :param conn: The connection to the database. If not provided, a new connection to an in-memory
        database is created
        """

        # config
        config = self.get_config()

        # Connexion config
        connexion_config = {}
        if config.get("threads", None):
            connexion_config["threads"] = config.get("threads")
        if config.get("memory", None):
            connexion_config["memory_limit"] = config.get("memory")
        default_connexion_db = ":memory:"

        # Connexion format
        connexion_format = self.get_config().get("connexion_format", "duckdb")

        # Conn
        connexion_db = default_connexion_db
        if not conn:
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

            if connexion_format in ["duckdb"]:
                conn = duckdb.connect(connexion_db, config=connexion_config)
            elif connexion_format in ["sqlite"]:
                conn = sqlite3.connect(connexion_db)

        # Compression
        # PRAGMA force_compression, expected Auto, Uncompressed, Constant, RLE, Dictionary, PFOR, BitPacking, FSST, Chimp, Patas
        # conn.execute("PRAGMA force_compression='Patas';")

        # Set connexion
        self.connexion_format = connexion_format
        self.connexion_db = connexion_db
        self.conn = conn

        log.debug(f"connexion_format: {connexion_format}")
        log.debug(f"connexion_db: {connexion_db}")
        log.debug(f"connexion config: {connexion_config}")


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


    def get_query_to_df(self, query: str = "") -> pd.DataFrame:
        """
        > The function `get_query_to_df` takes a query as a string and returns a pandas dataframe

        :param query: str = ""
        :type query: str
        :return: A dataframe
        """

        connexion_format = self.get_connexion_format()
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


    def get_stats(self) -> None:
        """
        The function prints statistics of the current object
        """

        # Log
        log.info(f"Stats Calculation...")

        # table varaints
        table_variants_from = self.get_table_variants()

        # Samples
        nb_of_samples = len(self.get_header_sample_list())

        # Variants by chr
        sql_query_nb_variant_by_chrom = f"SELECT \"#CHROM\" as CHROM, count(*) as count FROM {table_variants_from} GROUP BY \"#CHROM\""
        log.debug(f"Query Variants by Chr: {sql_query_nb_variant_by_chrom}")
        df_nb_of_variants_by_chrom = self.get_query_to_df(sql_query_nb_variant_by_chrom)
        nb_of_variants_by_chrom = df_nb_of_variants_by_chrom.sort_values(by=['CHROM'], kind='quicksort')

        # Total number of variants
        nb_of_variants = nb_of_variants_by_chrom["count"].sum()

        # Calculate percentage
        nb_of_variants_by_chrom['percent'] = nb_of_variants_by_chrom['count'].apply(lambda x: (x / nb_of_variants) * 100)

        # Genotypes
        genotypes = {}
        for sample in self.get_header_sample_list():
            sql_query_genotype = f"""
                SELECT  '{sample}' as sample,
                        REGEXP_EXTRACT("{sample}", '^([0-9/|.]*)[:]*',1) as genotype,
                        count(REGEXP_EXTRACT("{sample}", '^([0-9/|.]*)[:]*',1)) as count,
                        concat((count(REGEXP_EXTRACT("{sample}", '^([0-9/|.]*)[:]*',1))*100/{nb_of_variants}), '%') as percentage
                FROM {table_variants_from}
                WHERE 1
                GROUP BY genotype
                """
            genotypes[sample] = self.conn.execute(sql_query_genotype).df()

        # Output
        log.info("Number of Sample(s): " + str(nb_of_samples))
        log.info("Number of Variant(s): " + str(nb_of_variants))
        log.info("Number of Variant(s) by chromosomes:")
        for d in str(nb_of_variants_by_chrom).split("\n"):
            log.info("\t" + str(d))
        log.info("Number of Sample(s): " + str(nb_of_samples))
        log.info(f"Genotypes:")
        for sample in genotypes:
            for d in str(genotypes[sample]).split("\n"):
                log.info("\t" + str(d))

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

        if clause in ["select", "where", "update"]:
            table_variants = self.table_variants
        elif clause in ["from"]:
            if self.get_input_format() in ["parquet"] and access in ["RO"]:
                input_file = self.get_input()
                table_variants = f"'{input_file}' as variants"
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


    def load_data(self, input_file:str = None, drop_variants_table:bool = False) -> None:
        """
        It reads a VCF file and inserts it into a table
        
        :param input_file: The path to the input file
        :type input_file: str
        :param drop_variants_table: If True, the variants table will be dropped before loading the data,
        defaults to False
        :type drop_variants_table: bool (optional)
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

        # Input format and compress
        input_format = self.get_input_format()
        input_compressed = self.get_input_compressed()

        # Connexion format
        connexion_format = self.get_connexion_format()

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
                        FROM read_csv('{self.input}', auto_detect=True, skip={skip}, delim='{delimiter}')
                    """
                else:
                    sql_vcf = f"""
                    CREATE TABLE {table_variants} AS 
                        SELECT *
                        FROM read_csv('{self.input}', auto_detect=True, skip={skip}, delim='{delimiter}')
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
                    #sql_create_table_columns.append(f"\"{column}\" {column_type}")
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
        if self.get_param().get("explode_infos", None) is not None:
            self.explode_infos(
                prefix=self.get_param().get("explode_infos", None))

        # Create index after insertion
        self.create_indexes()


    def add_column(self, table_name, column_name, column_type, default_value=None):
        """
        Adds a column to a SQLite or DuckDB table with a default value if it doesn't already exist.

        table_name: the name of the table
        column_name: the name of the column to add
        column_type: the data type of the column to add
        default_value: the default value of the column to add (optional)
        """

        # Check if the column already exists in the table
        query = f""" SELECT * FROM {table_name} LIMIT 0 """
        columns = self.get_query_to_df(query).columns.tolist()
        if column_name in columns:
            log.debug(f"The {column_name} column already exists in the {table_name} table")
            return
        else:
            log.debug(f"The {column_name} column NOT exists in the {table_name} table")
        
        # Add column in table
        add_column_query = f""" ALTER TABLE {table_name} ADD COLUMN "{column_name}" {column_type} """
        if default_value is not None:
            add_column_query += f" DEFAULT {default_value}"
        self.execute_query(add_column_query)
        log.debug(f"The {column_name} column was successfully added to the {table_name} table")
        
        return


    def explode_infos(self, prefix: str = None, create_index: bool = False, fields: list = None, update: bool = True, force:bool = False) -> None:
        """
        The function takes a VCF file and explodes the INFO fields into individual columns
        """

        # drop indexes
        self.drop_indexes()

        # connexion format 
        connexion_format = self.get_connexion_format()

        # Access
        access = self.get_config().get("access", None)

        if access not in ["RO"]:

            # prefix
            if prefix in [None, True] or type(prefix) != str:
                if self.get_param().get("explode_infos",None) not in [None, True]:
                    prefix = self.get_param().get("explode_infos","INFO/")
                else:
                    prefix = "INFO/"

            # table variants
            table_variants = self.get_table_variants(clause="select")

            # extra infos
            extra_infos = self.get_extra_infos()

            log.debug(
                f"Explode INFO fields - ADD [{len(self.get_header().infos)}] annotations fields")

            sql_info_alter_table_array = []

            # Info fields to check
            fields_list = list(self.get_header().infos)
            if fields:
                fields_list += fields
            fields_list = set(fields_list)

            for info in fields_list:

                info_id_sql = prefix+info

                if (fields is None or info in fields or prefix+info in fields) and (update or info not in extra_infos):

                    log.debug(
                        f"Explode INFO fields - ADD '{info}' annotations fields")

                    if info in self.get_header().infos:
                        info_type = self.get_header().infos[info].type
                        info_num = self.get_header().infos[info].num
                    else:
                        info_type = "String"
                        info_num = 0

                    type_sql = self.code_type_map_to_sql.get(info_type, "VARCHAR")
                    if info_num != 1:
                        type_sql = "VARCHAR"

                    # Add field
                    self.add_column(table_name=table_variants, column_name=info_id_sql, column_type=type_sql, default_value="null")

                    # add field to index
                    self.index_additionnal_fields.append(info_id_sql)

                    # Update field array
                    if connexion_format in ["duckdb"]:
                        update_info_field = f"""
                        "{info_id_sql}" =
                            CASE
                                WHEN REGEXP_EXTRACT(concat(';', INFO), ';{info}=([^;]*)',1) == '' THEN NULL
                                WHEN REGEXP_EXTRACT(concat(';', INFO), ';{info}=([^;]*)',1) == '.' THEN NULL
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

            # By chromosomes
            chromosomes_df = self.get_query_to_df(
                f""" SELECT "#CHROM" FROM {table_variants} GROUP BY "#CHROM" """)

            for chrom in chromosomes_df["#CHROM"]:
                log.debug(
                    f"Explode INFO fields - Chromosome {chrom}...")
                # Update table
                sql_info_alter_table_array_join = ", ".join(
                    sql_info_alter_table_array)
                if sql_info_alter_table_array_join:
                    sql_info_alter_table = f"""
                        UPDATE {table_variants}
                        SET {sql_info_alter_table_array_join}
                        WHERE "#CHROM" = '{chrom}'
                        """
                    log.debug(
                        f"Explode INFO fields - ADD [{len(self.get_header().infos)}]: {sql_info_alter_table}")
                    self.conn.execute(sql_info_alter_table)

        # create indexes
        if create_index:
            self.create_indexes()


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


    def export_output(self, output_file: str = None, output_header: str = None, export_header: bool = True, query: str = None, parquet_partitions:list = None, threads:int = None, sort:bool = False, index:bool = False) -> bool:
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
        :param threads: Number of threads (optional)
        :type threads: int
        :param sort: sort output file, only if VCF format (optional)
        :type sort: bool
        :param index: index output file, only if VCF format (optional)
        :type index: int
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

        # Database
        database_source=self.get_connexion()

        # Connexion format
        connexion_format = self.get_connexion_format()

        # Tmp files to remove
        tmp_to_remove = []

        if connexion_format in ["sqlite"] or query:

            # Export in Parquet
            random_tmp = ''.join(random.choice(string.ascii_lowercase) for i in range(10))
            database_source = f"""{output_file}.{random_tmp}.database_export.parquet"""
            tmp_to_remove.append(database_source)

            # Table Variants
            table_variants = self.get_table_variants()

            # Create export query
            if query:
                sql_query_export_subquery = f"""
                    SELECT * FROM ({query})
                    """
            elif connexion_format in ["sqlite"]:
                sql_query_export_subquery = f"""
                    SELECT * FROM {table_variants}
                    """

            # Write source file
            fp.write(database_source, self.get_query_to_df(sql_query_export_subquery))

        # Create database
        database = Database(database=database_source, table="variants", header_file=output_header)

        # Export file
        database.export(output_database=output_file, parquet_partitions=parquet_partitions, threads=threads, sort=sort, index=index)

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
        query = f""" SELECT * FROM {table} LIMIT 1 """
        table_columns = self.get_query_to_df(query).columns.tolist()
        extra_columns = []
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
                    head_clean = head.replace("Type=Flag", "Type=String")
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
            # import csv
            # with pgzip.open(tmp_variants_name, 'wt') as f:
            #     writer = csv.writer(f, delimiter='\t',
            #                         quotechar='', quoting=csv.QUOTE_NONE)
            #     cursor = self.conn.execute(sql_query_select)
            #     writer.writerow([i[0] for i in cursor.description])
            #     writer.writerows(cursor)
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

        concat_and_compress_files(input_files=[tmp_header_name, tmp_variants_name], output_file=vcf_file, compression_type=compression_type, threads=threads, sort=True, index=index)


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

        # merged VCF
        tmp_merged_vcf = NamedTemporaryFile(prefix=self.get_prefix(
        ), dir=self.get_tmp_dir(), suffix=".parquet", delete=True)
        tmp_merged_vcf_name = tmp_merged_vcf.name

        # Fix error with multithreading
        config = self.get_config()
        config["threads"] = 1

        mergeVCF = Variants(
            None, vcf_file, tmp_merged_vcf_name, config=config)
        mergeVCF.load_data()
        mergeVCF.export_output(export_header=False)
        del mergeVCF
        gc.collect()

        # Original number of variants
        query = "SELECT count(*) AS count FROM variants"
        count_original_variants = self.conn.execute(query).df()["count"][0]

        # Annotated number of variants
        query = f"SELECT count(*) AS count FROM '{tmp_merged_vcf_name}' as table_parquet"
        count_annotated_variants = self.conn.execute(query).df()["count"][0]

        if count_original_variants != count_annotated_variants:
            log.warning(
                f"Update from VCF - Discodance of number of variants between database ({count_original_variants}) and VCF for update ({count_annotated_variants})")

        if count_annotated_variants:
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
                                FROM '{tmp_merged_vcf_name}' as table_parquet
                                        WHERE table_parquet.\"#CHROM\" = table_variants.\"#CHROM\"
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
        prefix = "INFO/"

        # Explode INFO/SVTYPE
        self.explode_infos(prefix=prefix,fields=["SVTYPE"])

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
            if not "annotation" in param:
                param["annotation"] = {}
            for annotation_file in param.get("annotations"):

                annotations = param.get("annotations").get(
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
                        # for annotation_file_ann in annotation_file_annotation.split("+"):
                        #     param["annotation"]["annovar"]["annotations"][annotation_file_ann] = annotations

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
                            # log.debug(f"find all: {annotation_file} {found_files}")
                            #raise ValueError("FIND ALL")
                            if len(found_files) > 0:
                                annotation_file_found = found_files[0]
                                break
                        if not annotation_file_found and not assembly:
                            # Find within folders
                            for annotations_database in annotations_databases:
                                found_files = find_all(annotation_file, annotations_database)
                                # log.debug(f"find all2: {annotation_file} {found_files}")
                                # raise ValueError("FIND ALL")
                                if len(found_files) > 0:
                                    annotation_file_found = found_files[0]
                                    break
                    log.debug(f"for {annotation_file} annotation_file_found={annotation_file_found}")

                    if annotation_file_found:

                        database = Database(database=annotation_file_found)
                        quick_annotation_format = database.get_format()

                        # Check Annotation Tool
                        annotation_tool = None
                        if quick_annotation_format in ["tsv", "tsv", "csv", "json", "tbl", "parquet", "duckdb"]:
                            annotation_tool = "parquet"
                        elif quick_annotation_format in ["vcf", "bed"]:
                            annotation_tool = "bcftools"
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
        if self.get_param().get("explode_infos", None) is not None:
            self.explode_infos(
                prefix=self.get_param().get("explode_infos", None))


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

                # if not db_file or (not db_hdr_file and db_file_format not in ["gz"]):
                if not db_file or not db_hdr_file:
                    log.error("Annotation failed: file not found")
                    log.error(f"Annotation annotation file: {db_file}")
                    log.error(f"Annotation annotation header: {db_hdr_file}")
                    raise ValueError(
                        f"Annotation failed: databases not found - annotation file {db_file} / annotation header {db_hdr_file}")
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

        :param threads: The number of threads to use
        :return: the value of the variable "return_value".
        """

        # DEBUG
        log.debug("Start annotation with Exomiser databases")

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

        # # Param - HPO
        # options_hpo = param_exomiser.get("hpo", None)
        # log.debug("HPO: " + str(options_hpo))

        # options = param_exomiser.get("options", None)
        # log.debug("Options: " + str(options))

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

        # # Export in VCF
        # log.debug("Create initial file to annotate")
        # tmp_vcf = NamedTemporaryFile(prefix=self.get_prefix(), dir=self.get_tmp_dir(), suffix=".vcf.gz", delete=True)
        # tmp_vcf_name = tmp_vcf.name
        # log.debug(f"VCF exported: {tmp_vcf_name}")

        # VCF header
        vcf_reader = self.get_header()
        log.debug("Initial header: " + str(vcf_reader.infos))

        # Samples
        samples = self.get_header_sample_list()
        if not samples:
            log.error("No Samples in VCF")
            return False
        log.debug(f"Samples: {samples}")
        
        # Existing annotations
        # for vcf_annotation in self.get_header().infos:

        #     vcf_annotation_line = self.get_header().infos.get(vcf_annotation)
        #     log.debug(
        #         f"Existing annotations in VCF: {vcf_annotation} [{vcf_annotation_line}]")

        force_update_annotation = True

        if "Exomiser" not in self.get_header().infos or force_update_annotation:
            log.debug("Start annotation Exomiser")

            with TemporaryDirectory(dir=self.get_tmp_dir()) as tmp_dir:

                

                # time java -XX:ParallelGCThreads=8 -Xms2g -Xmx4g -jar exomiser-cli-13.2.0/exomiser-cli-13.2.0.jar  --analysis=exomiser-cli-13.2.0/examples/test.yml --spring.config.location=/databases/exomiser/current/hg19 --exomiser.data-directory=/databases/exomiser/current/hg19 --exomiser.hg19.data-version=2302 --exomiser.phenotype.data-version=2302
                # time java -XX:ParallelGCThreads=8 -Xms2g -Xmx4g -jar exomiser-cli-13.2.0/exomiser-cli-13.2.0.jar  --analysis=exomiser-cli-13.2.0/examples/test.yml --exomiser.data-directory=/databases/exomiser/current/hg19 --exomiser.hg19.data-version=2302 --exomiser.phenotype.data-version=2302
                # time java -XX:ParallelGCThreads=8 -Xms2g -Xmx4g -jar exomiser-cli-13.2.0/exomiser-cli-13.2.0.jar  --sample=exomiser-cli-13.2.0/examples/sample.test.yml --exomiser.data-directory=/databases/exomiser/current/hg19 --exomiser.hg19.data-version=2302 --exomiser.phenotype.data-version=2302
                # time java -XX:ParallelGCThreads=8 -Xms2g -Xmx4g -jar exomiser-cli-13.2.0/exomiser-cli-13.2.0.jar  --sample=exomiser-cli-13.2.0/examples/sample.test.json --exomiser.data-directory=/databases/exomiser/current/hg19 --exomiser.hg19.data-version=2302 --exomiser.phenotype.data-version=2302
                # # java -jar exomiser-cli-13.2.0.jar --sample examples/sample.test.json --preset exome --output-directory=/tmp/exomiser_output --output-format=HTML,JSON,TSV_GENE,TSV_VARIANT,VCF --outputContributingVariantsOnly=False --spring.config.location=/databases/exomiser/current/hg19/application.properties --exomiser.data-directory=/databases/exomiser/current/hg19 --exomiser.hg19.data-version=2302


                ######  examples/sample.test.json

                # {
                #   "id": "manuel",
                #   "subject": { "id": "manuel", "sex": "UNKNOWN_SEX" },
                #   "phenotypicFeatures": [
                #     { "type": { "id": "HP:0001159" } },
                #     { "type": { "id": "HP:0000486" } },
                #     { "type": { "id": "HP:0000327" } },
                #     { "type": { "id": "HP:0000520" } },
                #     { "type": { "id": "HP:0000316" } },
                #     { "type": { "id": "HP:0000244" } }
                #   ],
                #   "htsFiles": [
                #     {
                #       "uri": "examples/Pfeiffer.vcf",
                #       "htsFormat": "VCF",
                #       "genomeAssembly": "hg19"
                #     }
                #   ]
                # }


                ###### examples/output-options.test.yml
                # #
                # ---
                # outputContributingVariantsOnly: false
                # #numGenes options: 0 = all or specify a limit e.g. 500 for the first 500 results
                # numGenes: 0
                # #minExomiserGeneScore: 0.7
                # # Path to the desired output directory. Will default to the 'results' subdirectory of the exomiser install directory
                # outputDirectory: /tmp/exomiser_results
                # # Filename for the output files. Will default to {input-vcf-filename}-exomiser
                # outputFileName: exomiser-output
                # #out-format options: HTML, JSON, TSV_GENE, TSV_VARIANT, VCF (default: HTML)
                # outputFormats: [ HTML, JSON, TSV_GENE, TSV_VARIANT, VCF ]

                # Create sample YML/JSON
                #exomiser_config_sample = {}

                # with open('/databases/exomiser/13.2.0/exomiser-cli-13.2.0/examples/sample.test.yml', 'r') as file:
                #     configuration = yaml.safe_load(file)

                # with open('/databases/exomiser/13.2.0/exomiser-cli-13.2.0/examples/sample.test.json', 'w') as json_file:
                #     json.dump(configuration, json_file)
                    
                # output = json.dumps(json.load(open('/databases/exomiser/13.2.0/exomiser-cli-13.2.0/examples/sample.test.json')), indent=2)
                # print(output)

                # {
                # "id": "manuel",
                # "subject": { "id": "manuel", "sex": "UNKNOWN_SEX" },
                # "phenotypicFeatures": [
                #     { "type": { "id": "HP:0001159" } },
                #     { "type": { "id": "HP:0000486" } },
                #     { "type": { "id": "HP:0000327" } },
                #     { "type": { "id": "HP:0000520" } },
                #     { "type": { "id": "HP:0000316" } },
                #     { "type": { "id": "HP:0000244" } }
                # ],
                # "htsFiles": [
                #     {
                #     "uri": "exomiser-cli-13.2.0/examples/Pfeiffer.vcf",
                #     "htsFormat": "VCF",
                #     "genomeAssembly": "hg19"
                #     }
                # ]
                # }


                # STEPS
                
                # Initial file name
                log.debug("Create initial file to annotate")
                tmp_vcf_name = os.path.join(tmp_dir, "initial.vcf.gz")
                log.debug(f"VCF exported: {tmp_vcf_name}")


                # Create analysis.json (either analysis in param or by default, depending on preset exome/genome)
                log.debug(f"file_folder: {file_folder}")
                log.debug(f"folder_main: {folder_main}")
                if os.path.exists(folder_main):
                    log.debug(f"main folder exists: {folder_main}")
                if os.path.exists(folder_config):
                    log.debug(f"config folder exists: {folder_config}")

                log.debug(f"param_exomiser: {param_exomiser}")

                param_exomiser_analysis_dict = {}
                param_exomiser_analysis = param_exomiser.get("analysis", {})
                #param_exomiser_analysis = param_exomiser.get("analysis", {"test": "truc"})
                #param_exomiser_analysis = param_exomiser.get("analysis", os.path.join(folder_config, "preset-exome-analysis.json"))
                if param_exomiser_analysis:
                    if isinstance(param_exomiser_analysis, str) and os.path.exists(param_exomiser_analysis):
                        log.debug(f"is str/file")
                        with open(param_exomiser_analysis) as json_file:
                            param_exomiser_analysis_dict = json.load(json_file)
                    elif isinstance(param_exomiser_analysis, dict):
                        log.debug(f"is dict")
                        param_exomiser_analysis_dict = param_exomiser_analysis
                    else:
                        log.error(f"analysis type unknown")

                log.debug(param_exomiser_analysis_dict)

                # Case no input analysis config file/dict
                # Use preset (exome/genome) to open default config file
                #param_exomiser_preset = param_exomiser.get("preset", "")
                if not param_exomiser_analysis_dict or "analysis" in param_exomiser_analysis_dict:
                    param_exomiser_preset = param_exomiser.get("preset", "exome")
                    param_exomiser_analysis_default_config_file = os.path.join(folder_config, f"preset-{param_exomiser_preset}-analysis.json")
                    if os.path.exists(param_exomiser_analysis_default_config_file):
                        log.debug(f"is str/file default config ({param_exomiser_analysis_default_config_file})")
                        with open(param_exomiser_analysis_default_config_file) as json_file:
                            param_exomiser_analysis_dict["analysis"] = json.load(json_file)
                    else:
                        log.error(f"No analysis default config file ({param_exomiser_analysis_default_config_file})")

                if not param_exomiser_analysis_dict:
                    log.error(f"No analysis config")

                log.debug(param_exomiser_analysis_dict)


                

                # Add input VCF ???
                    # "genomeAssembly": "hg19",
                    # "vcf": "/databases/exomiser/13.2.0/exomiser-cli-13.2.0/examples/Pfeiffer.test.hom.vcf",
                    # "ped": null,
                    # "proband": null,
                    # "hpoIds": [
                    #     "HP:0001156",
                    #     "HP:0001363",
                    #     "HP:0011304",
                    #     "HP:0010055"
                    # ]
                    #
                    # OR
                    # 
                    # 
                    # "phenopacket": {
                    #     "id": "manuell",
                    #     "subject": {
                    #         "id": "manuel",
                    #         "sex": "MALE"
                    #     },
                    #     "phenotypicFeatures": [
                    #     {
                    #         "type": {
                    #         "id": "HP:0001159",
                    #         "label": "Syndactyly"
                    #         }
                    #     },
                    #     {
                    #         "type": {
                    #         "id": "HP:0000486",
                    #         "label": "Strabismus"
                    #         }
                    #     },
                    #     {
                    #         "type": {
                    #         "id": "HP:0000327",
                    #         "label": "Hypoplasia of the maxilla"
                    #         }
                    #     },
                    #     {
                    #         "type": {
                    #         "id": "HP:0000520",
                    #         "label": "Proptosis"
                    #         }
                    #     },
                    #     {
                    #         "type": {
                    #         "id": "HP:0000316",
                    #         "label": "Hypertelorism"
                    #         }
                    #     },
                    #     {
                    #         "type": {
                    #         "id": "HP:0000244",
                    #         "label": "Brachyturricephaly"
                    #         }
                    #     }
                    #     ],
                    #     "htsFiles": [
                    #     {
                    #         "uri": "/databases/exomiser/13.2.0/exomiser-cli-13.2.0/examples/Pfeiffer.vcf",
                    #         "htsFormat": "VCF",
                    #         "genomeAssembly": "hg19"
                    #     }
                    #     ],
                    #     "metaData": {
                    #     "created": "2019-11-12T13:47:51.948Z",
                    #     "createdBy": "julesj",
                    #     "resources": [
                    #         {
                    #         "id": "hp",
                    #         "name": "human phenotype ontology",
                    #         "url": "http://purl.obolibrary.org/obo/hp.owl",
                    #         "version": "hp/releases/2019-11-08",
                    #         "namespacePrefix": "HP",
                    #         "iriPrefix": "http://purl.obolibrary.org/obo/HP_"
                    #         }
                    #     ],
                    #     "phenopacketSchemaVersion": 1
                    #     }
                    # },

                if "phenopacket" not in param_exomiser_analysis_dict:

                    param_exomiser_analysis_dict["phenopacket"] = {
                        "id": "analysis"
                    }

                    # Add subject
                    param_exomiser_subject = param_exomiser.get("subject", {})

                    if not param_exomiser_subject:
                        
                        # Find sample ID (first sample)
                        sample_list = self.get_header_sample_list()
                        #log.debug(f"sample_list: {sample_list}")
                        if len(sample_list) > 0:
                            sample = sample_list[0]
                        else:
                            log.error(f"No sample found")
                            # RAISE!!!

                        param_exomiser_subject = {"id": sample,"sex": "UNKNOWN_SEX"}

                    param_exomiser_analysis_dict["phenopacket"]["subject"] = param_exomiser_subject

                    # Add "phenotypicFeatures"
                    param_exomiser_phenotypicFeatures = param_exomiser.get("phenotypicFeatures", [])

                    # Try to infer from hpo list
                    if not param_exomiser_phenotypicFeatures:
                        param_exomiser_hpo = param_exomiser.get("hpo", [])

                        # if not param_exomiser_hpo:
                        #     log.debug(f"No HPO found")

                        #log.debug(f"HPO: {param_exomiser_hpo}")
                        for hpo in param_exomiser_hpo:
                            #log.debug(f"hpo: {hpo}")
                            #hpo_clean = hpo.replace(["0", "1"], "L")
                            characters_to_remove = ['h', 'p', 'o', 'H', 'P', 'O', ' ', ':']
                            pattern = '[' +  ''.join(characters_to_remove) +  ']'
                            hpo_clean = re.sub(pattern, '', hpo)
                            #log.debug(f"hpo_clean: {hpo_clean}")
                            param_exomiser_phenotypicFeatures.append({
                                "type": {
                                    "id": f"HP:{hpo_clean}",
                                    "label": "Syndactyly"
                                }
                            })

                    #log.debug(f"param_exomiser_phenotypicFeatures: {param_exomiser_phenotypicFeatures}")

                    if not param_exomiser_phenotypicFeatures:
                        param_exomiser_analysis_dict["phenopacket"]["phenotypicFeatures"] = []
                        #log.error(f"No phenotypicFeatures found")
                        # THEN REMOVE "hiPhivePrioritiser": {}
                        #param_exomiser_analysis_dict["analysis"]
                        for step in param_exomiser_analysis_dict.get("analysis", {}).get("steps", []):
                            #log.debug(f"step: {step}")
                            if "hiPhivePrioritiser" in step:
                                param_exomiser_analysis_dict.get("analysis", {}).get("steps", []).remove(step)

                        #param_exomiser_analysis_dict["analysis"].pop('k4', None)

                    # Add Input File
                    #if "htsFiles" not in param_exomiser_analysis_dict:
                    param_exomiser_analysis_dict["phenopacket"]["htsFiles"] = [
                        {
                            "uri": tmp_vcf_name,
                            "htsFormat": "VCF",
                            "genomeAssembly": assembly
                        }
                    ]


                    if "metaData" not in param_exomiser_analysis_dict:
                        log.debug(f"NO 'metaData'")
                        param_exomiser_analysis_dict["phenopacket"]["metaData"] = {
                            "created": f"{datetime.datetime.now()}".replace(" ", "T") + "Z",
                            "createdBy": "howard",
                            "phenopacketSchemaVersion": 1
                        }
                    # if "metaData" not in param_exomiser_analysis_dict:
                    #     log.debug(f"NO 'metaData'")
                    #     # param_exomiser_analysis_dict["phenopacket"]["metaData"] =  {
                    #     #     "created": f"{datetime.datetime.now()}".replace(" ", "T"),
                    #     #     "createdBy": "howard",
                    #     #     "resources": [
                    #     #         {
                    #     #         "id": "hp",
                    #     #         "name": "human phenotype ontology",
                    #     #         "url": "http://purl.obolibrary.org/obo/hp.owl",
                    #     #         "version": "unknown",
                    #     #         "namespacePrefix": "HP",
                    #     #         "iriPrefix": "http://purl.obolibrary.org/obo/HP_"
                    #     #         }
                    #     #     ],
                    #     #     "phenopacketSchemaVersion": 1
                    #     # }
                    #     param_exomiser_analysis_dict["phenopacket"]["metaData"] =  {
                    #         "created": f"{datetime.datetime.now()}".replace(" ", "T") + "Z",
                    #         "createdBy": "julesj",
                    #         "resources": [
                    #             {
                    #             "id": "hp",
                    #             "name": "human phenotype ontology",
                    #             "url": "http://purl.obolibrary.org/obo/hp.owl",
                    #             "version": "hp/releases/2019-11-08",
                    #             "namespacePrefix": "HP",
                    #             "iriPrefix": "http://purl.obolibrary.org/obo/HP_"
                    #             }
                    #         ],
                    #         "phenopacketSchemaVersion": 1
                    #     }


                    param_exomiser_analysis_dict["phenopacket"]["phenotypicFeatures"] = param_exomiser_phenotypicFeatures


                # Find samples
                # Main sample
                sample = param_exomiser_analysis_dict.get("phenopacket",{}).get("subject",{}).get("id", None)
                if sample is None:
                    sample = []
                else:
                    sample = [sample]
                # Pedigree
                pedigree_persons_list = param_exomiser_analysis_dict.get("phenopacket",{}).get("pedigree",{}).get("persons",{})
                pedigree_persons = []
                for person in pedigree_persons_list:
                    pedigree_persons.append(person.get("individualId"))

                # Concat samples
                samples = sample + pedigree_persons

                # Check if sample
                if not samples:
                    log.error(f"No samples found")
                    # RAISE !!!

                # OutputOptions
                output_results = os.path.join(tmp_dir, "results")
                if "outputOptions" not in param_exomiser_analysis_dict:
                    log.debug(f"NO 'outputOptions'")
                    param_exomiser_analysis_dict["outputOptions"] = {
                        "outputContributingVariantsOnly": False,
                        "numGenes": 0,
                        "outputDirectory": output_results,
                        "outputFileName": "howard",
                        #"outputFormats": ["HTML", "JSON", "TSV_GENE", "TSV_VARIANT", "VCF"]
                        "outputFormats": ["TSV_VARIANT", "VCF"]
                    }

                else:
                    log.debug(f"Found 'outputOptions'")
                    param_exomiser_analysis_dict["outputOptions"]["outputDirectory"] = output_results
                    param_exomiser_analysis_dict["outputOptions"]["outputFormats"] = list(set(param_exomiser_analysis_dict.get("outputOptions",{}).get("outputFormats", []) + ["TSV_VARIANT", "VCF"]))

                    
                # Create VCF with sample (either sample in param or first one by default)

                # Export VCF file
                self.export_variant_vcf(vcf_file=tmp_vcf_name, file_type="gz", remove_info=True, add_samples=True, list_samples=samples, compression=1, index=False)

                
                log.debug(f"param_exomiser_analysis_dict: {param_exomiser_analysis_dict}")

                exomiser_analysis = os.path.join(tmp_dir, "analysis.json")
                with open(exomiser_analysis, 'w') as fp:
                    json.dump(param_exomiser_analysis_dict, fp)

                with open(exomiser_analysis) as json_file:
                    analysis_json = json.load(json_file)
                    log.debug(analysis_json)
                    log.debug(json.dumps(analysis_json, indent=4))


                shutil.copyfile(tmp_vcf_name, "/tmp/initial.vcf.gz")
                shutil.copyfile(exomiser_analysis, "/tmp/analysis.json")


                # java -jar exomiser-cli-13.2.0.jar --analysis /tmp/analysis.json --vcf /tmp/initial.vcf.gz --spring.config.location=/databases/exomiser/current/hg19/application.properties --exomiser.data-directory=/databases/exomiser/current/hg19
                # java -jar exomiser-cli-13.2.0.jar --analysis /tmp/analysis.tmp.json --spring.config.location=/databases/exomiser/current/hg19/application.properties --exomiser.data-directory=/databases/exomiser/current/hg19
                # java -jar exomiser-cli-13.2.0.jar --analysis /tool/tests/tmp/analysis.tmp.json --spring.config.location=/databases/exomiser/current/hg19/application.properties --exomiser.data-directory=/databases/exomiser/current/hg19
                # 

                #exit()
                
                # Execute Exomiser
                # --analysis examples/test.json --spring.config.location=/databases/exomiser/current/hg19/application.properties --exomiser.data-directory=/databases/exomiser/current/hg19
                # databases_folders
                exomiser_command = ""

                #exomiser_java_options = f" -XX:ParallelGCThreads={threads} -Xms2g -Xmx4g "
                exomiser_java_options = f" -XX:ParallelGCThreads={threads} "
                #exomiser_java_options = ""
                


                #exomiser_command = f"{java_bin} -jar {exomiser_jar} --analysis={exomiser_analysis} --spring.config.location={databases_folders}/{assembly}/application.properties --exomiser.data-directory={databases_folders}/{assembly}"
                exomiser_command = f" {exomiser_java_options} -jar {exomiser_jar} --analysis={exomiser_analysis} --spring.config.location={databases_folders}/{assembly}/application.properties --exomiser.data-directory={databases_folders}/{assembly}"

                log.debug(exomiser_command)

                result = subprocess.call([java_bin] + exomiser_command.split(), stdout=subprocess.PIPE)
                #result = subprocess.call([java_bin] + exomiser_command.split())

                log.debug(glob.glob(f'{output_results}/*'))

                log.debug(self.conn.query(f"SELECT \"#CHROM\", POS, REF, ALT, INFO FROM {table_variants} as table_variants").df())


                output_results_tsv = os.path.join(output_results,"howard.variants.tsv")
                if os.path.exists(output_results_tsv):


                    query = f""" SELECT * FROM read_csv('{output_results_tsv}', auto_detect=True, delim='\t') LIMIT {DTYPE_LIMIT_AUTO} """
                    output_results_tsv_df = self.get_query_to_df(query)
                    output_results_tsv_columns = output_results_tsv_df.columns.tolist()
                    log.debug(output_results_tsv_df.columns)
                    log.debug(output_results_tsv_df.dtypes)
                    log.debug(output_results_tsv_columns)

                    sql_query_update_concat_fields = []

                    # List all columns to add into header
                    for header_column in output_results_tsv_columns:

                        if header_column not in ["CONTIG", "START", "END", "REF", "ALT", "QUAL", "FILTER", "GENOTYPE"]:

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
                            #header_info_name = f"Exomiser_{header_column}".replace("#","").replace("-","")
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
                            # vcf.parser._Info
                            log.debug(f"{header_column} >>> {header_info_name}")
                            #vcf_reader.infos[header_info_name] = database_header_to_add

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

                        # ",".join(sql_query_update_concat_fields)

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
                                        FROM read_csv('{output_results_tsv}', auto_detect=True, delim='\t') as table_parquet
                                                WHERE concat('chr', CAST(table_parquet.\"CONTIG\" AS STRING)) = table_variants.\"#CHROM\"
                                                AND table_parquet.\"START\" = table_variants.\"POS\"
                                                AND table_parquet.\"ALT\" = table_variants.\"ALT\"
                                                AND table_parquet.\"REF\" = table_variants.\"REF\"
                                        )
                                    )
                        ;
                        """
                    log.debug(sql_query_update)
                    self.conn.execute(sql_query_update)


                #, auto_detect=True, skip={skip}, delim='{delimiter}'

                log.debug(self.conn.query("SELECT \"#CHROM\", POS, REF, ALT, INFO FROM variants").df())

                output_results_vcf = os.path.join(output_results,"howard.vcf.gz")
                if os.path.exists(output_results_vcf):

                    # Find annotation in header
                    with gzip.open(output_results_vcf, 'rt') as f:
                        header_list = self.read_vcf_header(f)
                    exomiser_vcf_header = vcf.Reader(
                        io.StringIO("\n".join(header_list)))
                    log.debug(exomiser_vcf_header.infos)
                    log.debug(exomiser_vcf_header.infos["Exomiser"])

                    vcf_reader.infos["Exomiser"] = exomiser_vcf_header.infos["Exomiser"]
                #    vcf_reader.infos["Exomiser"] = vcf.parser._Info(
                #         "Exomiser",
                #         ".",
                #         "String",
                #         "A pipe-separated set of values for the proband allele(s) from the record with one per compatible MOI following the format: {RANK|ID|GENE_SYMBOL|ENTREZ_GENE_ID|MOI|P-VALUE|EXOMISER_GENE_COMBINED_SCORE|EXOMISER_GENE_PHENO_SCORE|EXOMISER_GENE_VARIANT_SCORE|EXOMISER_VARIANT_SCORE|CONTRIBUTING_VARIANT|WHITELIST_VARIANT|FUNCTIONAL_CLASS|HGVS|EXOMISER_ACMG_CLASSIFICATION|EXOMISER_ACMG_EVIDENCE|EXOMISER_ACMG_DISEASE_ID|EXOMISER_ACMG_DISEASE_NAME}",
                #         "Exomiser",
                #         "unknown",
                #         CODE_TYPE_MAP["String"]
                #     )
                    # vcf.parser._Info
                    #vcf_reader.infos["Exomiser"] = database_header_to_add

                    log.debug(f"VCF result found: {output_results_vcf}")
                    # Update variants
                    log.info(f"Annotation - Updating...")
                    self.update_from_vcf(output_results_vcf)

                #log.debug(self.conn.query("SELECT INFO FROM variants WHERE INFO LIKE '%Exomiser%'").df())
                log.debug(self.conn.query("SELECT \"#CHROM\", POS, REF, ALT, INFO FROM variants").df())



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
            snpeff_command = f"{java_bin} -Xmx4g -jar {snpeff_jar} {assembly} -dataDir {snpeff_databases} {snpeff_options} {tmp_vcf_name} 1>{tmp_annotate_vcf_name} 2>>{tmp_annotate_vcf_name_err}"
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
        prefix=self.get_param().get("explode_infos", None)
        
        # explode infos
        self.explode_infos(prefix=prefix)

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

                    #if allow_annotation_full_info and nb_annotation_field == len(annotation_fields) and annotation_fields_ALL and "INFO" in parquet_hdr_vcf_header_columns:
                    #if allow_annotation_full_info and nb_annotation_field == len(annotation_fields) and annotation_fields_ALL and ("INFO" in parquet_hdr_vcf_header_columns or "INFO" in database.get_extra_columns()):
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
            explode_infos_prefix = self.get_param().get("explode_infos", "INFO/")
            if explode_infos_prefix == True:
                explode_infos_prefix = "INFO/"
            self.explode_infos(prefix=explode_infos_prefix)
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
                    self.add_column(table_name=table_variants, column_name=pzfield, column_type="INTEGER", default_value="0")
                elif re.match("PZFlag.*", pzfield):
                    self.add_column(table_name=table_variants, column_name=pzfield, column_type="BOOLEAN", default_value="1")
                else:
                    self.add_column(table_name=table_variants, column_name=pzfield, column_type="STRING", default_value="''")

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
                            if not explode_infos_prefix+annotation in extra_infos:
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

        # Explode INFOS fields into table fields
        if self.get_param().get("explode_infos", None) is not None:
            self.explode_infos(
                prefix=self.get_param().get("explode_infos", None))


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

        # Add hgvs column in variants table
        self.add_column(table_variants, "hgvs", "STRING", default_value=None)

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
        ddf["hgvs"] = ddf.map_partitions(partition_function)

        # Convert Dask DataFrame to Pandas Dataframe
        df = ddf.compute()
        
        # Convert Pandas dataframe to parquet (due to error in cast VARCHAR -> NULL ???)
        with tempfile.TemporaryDirectory() as tmpdir:
            df_parquet = os.path.join(tmpdir,"df.parquet")
            df.to_parquet(df_parquet)

            # Update hgvs column
            update_variant_query = f"""
                UPDATE {table_variants}
                SET hgvs=df.hgvs
                FROM read_parquet('{df_parquet}') as df
                WHERE variants."#CHROM" = df.CHROM
                AND variants.POS = df.POS
                AND variants.REF = df.REF
                AND variants.ALT = df.ALT
                AND df.hgvs NOT IN ('') AND df.hgvs NOT NULL
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
                    hgvs
                )
            WHERE hgvs NOT IN ('') AND hgvs NOT NULL
            """
        self.execute_query(sql_query_update)


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
                            WHEN "INFO/SVTYPE" NOT NULL THEN "INFO/SVTYPE"
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
        if self.get_param().get("explode_infos", None) is not None:
            self.explode_infos(
                prefix=self.get_param().get("explode_infos", None))


    def calculation_process_sql(self, operation:dict, operation_name:str = "unknown") -> None:
        """
        This function takes in a string of a mathematical operation and returns the result of that
        operation

        :param operation: The operation to be performed
        """

        # Param
        param = self.get_param()
        prefix = param.get("explode_infos", "INFO/")

        # table variants
        table_variants = self.get_table_variants(clause="alter")

        # Operation infos
        operation_name = operation.get('name', 'unknown')
        log.debug(f"process sql {operation_name}")
        output_column_name = operation.get('output_column_name',operation_name)
        output_column_type = operation.get('output_column_type', 'String')
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

                # Create column
                self.add_column(table_name=table_variants, column_name=prefix+output_column_name, column_type=output_column_type_sql, default_value="null")

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
                self.explode_infos(fields=[output_column_name] + operation_info_fields, force=True)

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

            else:
                log.error(f"Operations config: Calculation '{operation_name}' DOES NOT contain all mandatory fields {operation_info_fields}")
                raise ValueError(f"Operations config: Calculation '{operation_name}' DOES NOT contain all mandatory fields {operation_info_fields}")

        else:
            log.error(f"Operations config: Calculation '{operation_name}' query NOT defined")
            raise ValueError(f"Operations config: Calculation '{operation_name}' query NOT defined")


    def calculation_process_function(self, operation:dict, operation_name:str = "unknown") -> None:
        """
        This function takes in a string, and returns a string

        :param operation: The operation to be performed
        """
        operation_name = operation['name']
        log.debug(f"process sql {operation_name}")
        function_name = operation['function_name']
        function_params = operation['function_params']
        getattr(self, function_name)(*function_params)

    # Operation functions


    def calculation_variant_id(self) -> None:

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

        # SnpEff annotation field
        snpeff_ann = "ANN"

        # SnpEff annotation field
        snpeff_hgvs = "snpeff_hgvs"


        # Snpeff hgvs tags
        vcf_infos_tags = {
            snpeff_hgvs: "HGVS nomenclatures from snpEff annotation",
        }

        # Param
        param = self.get_param()
        prefix = param.get("explode_infos", "INFO/")
        if prefix == True:
            prefix = "INFO/"

        speff_ann_infos = prefix+snpeff_ann
        speff_hgvs_infos = prefix+snpeff_hgvs

        # Variants table
        table_variants = self.get_table_variants()
        
        # Header
        vcf_reader = self.get_header()

        # Explode HGVS field in column
        self.explode_infos(fields=[snpeff_ann])

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
        prefix = param.get("explode_infos", "INFO/")
        if prefix == True:
            prefix = "INFO/"

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
        # Explode HGVS field in column
        self.explode_infos(fields=[hgvs_field])

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
                            WHEN dataframe_hgvs."{nomen_field}" NOT NULL
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


    def calculation_find_by_pipeline(self, tag:str = "findbypipeline") -> None:

        # if FORMAT and samples
        if "FORMAT" in self.get_header_columns_as_list() and self.get_header_sample_list():

            # findbypipeline annotation field
            findbypipeline_tag = tag
            
            # VCF infos tags
            vcf_infos_tags = {
                findbypipeline_tag: f"Number of pipeline/sample for a variant ({findbypipeline_tag})",
            }

            # Param
            param = self.get_param()
            prefix = param.get("explode_infos", "INFO/")
            if prefix == True:
                prefix = "INFO/"

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

        # if FORMAT and samples
        if "FORMAT" in self.get_header_columns_as_list() and self.get_header_sample_list():

            # genotypeconcordance annotation field
            genotypeconcordance_tag = "genotypeconcordance"
            
            # VCF infos tags
            vcf_infos_tags = {
                genotypeconcordance_tag: "Concordance of genotype for multi caller VCF",
            }

            # Param
            param = self.get_param()
            prefix = param.get("explode_infos", "INFO/")
            if prefix == True:
                prefix = "INFO/"

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

        # if FORMAT and samples
        if "FORMAT" in self.get_header_columns_as_list() and self.get_header_sample_list():

            # barcode annotation field
            barcode_tag = "barcode"
            
            # VCF infos tags
            vcf_infos_tags = {
                "barcode": "barcode calculation (VaRank)",
            }

            # Param
            param = self.get_param()
            prefix = param.get("explode_infos", "INFO/")
            if prefix == True:
                prefix = "INFO/"

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
            prefix = param.get("explode_infos", "INFO/")
            if prefix == True:
                prefix = "INFO/"

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

        # if FORMAT and samples
        if "FORMAT" in self.get_header_columns_as_list() and self.get_header_sample_list():

            # vaf_normalization annotation field
            vaf_normalization_tag = "VAF"
            
            # VCF infos tags
            vcf_infos_tags = {
                "VAF": "VAF Variant Frequency",
            }

            # Param
            param = self.get_param()
            prefix = param.get("explode_infos", "INFO/")
            if prefix == True:
                prefix = "INFO/"

            #vaf_normalization_infos = prefix+vaf_normalization_tag

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
        
            # Param
            param = self.get_param()
            prefix = param.get("explode_infos", "INFO/")
            if prefix == True:
                prefix = "INFO/"

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


