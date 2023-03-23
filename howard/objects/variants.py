import csv
import gc
import gzip
import io
import multiprocessing
import os
import re
import shlex
import sqlite3
import subprocess
from tempfile import NamedTemporaryFile, TemporaryDirectory
import tempfile
import duckdb
import json
import argparse
import Bio.bgzf as bgzf
import pandas as pd
import vcf
import logging as log

from howard.commons import *


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

    # INIT section

    def set_input(self, input: str = None) -> None:
        """
        The function takes a file name as input, splits the file name into a name and an extension, and
        then sets the input_name, input_extension, and input_format attributes of the class

        :param input: The input file
        """
        self.input = input
        # Input format
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
        if config.get("memory_limit", None):
            connexion_config["memory_limit"] = config.get("memory_limit")
        default_connexion_db = ":memory:"

        # Connexion format
        connexion_format = self.get_config().get("connexion_format", "duckdb")

        # Conn
        connexion_db = default_connexion_db
        if not conn:
            if self.get_input_format() in ["db", "duckdb"]:
                connexion_db = self.get_input()
                # self.set_output(self.get_input())
            elif self.get_output_format() in ["db", "duckdb"]:
                connexion_db = self.get_output()
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
        config = self.get_config()
        header_list = []
        if self.input_format in ["gz"]:
            with bgzf.open(input_file, 'rt') as f:
                header_list = self.read_vcf_header(f)
        elif self.input_format in ["parquet", "tsv", "csv", "psv", "db", "duckdb"]:
            if config.get("header_file", None):
                with open(config.get("header_file"), 'rt') as f:
                    header_list = self.read_vcf_header(f)
            elif os.path.exists((input_file+".hdr")):
                with open(input_file+".hdr", 'rt') as f:
                    header_list = self.read_vcf_header(f)
            else:
                log.error(f"No header for file {input_file}")
                raise ValueError(f"No header for file {input_file}")
        else: # such as VCF
            with open(input_file, 'rt') as f:
                header_list = self.read_vcf_header(f)

        # header as list
        self.header_list = header_list

        # header as VCF object
        self.header_vcf = vcf.Reader(io.StringIO("\n".join(header_list)))

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

        # table varaints
        table_variants_from = self.get_table_variants()

        # Samples
        nb_of_samples = len(self.get_header_sample_list())

        # Variants
        sql_query_nb_variant = f"SELECT count(*) as count FROM {table_variants_from}"
        nb_of_variants = self.conn.execute(
            sql_query_nb_variant).df()["count"][0]

        # Variants by chr
        sql_query_nb_variant = f"SELECT \"#CHROM\" as CHROM, count(*) as count, (count(*)*100/{nb_of_variants}) || '%' as percentage FROM {table_variants_from} GROUP BY \"#CHROM\""
        nb_of_variants_by_chrom = self.conn.execute(
            sql_query_nb_variant).df().sort_values(by=['CHROM'], kind='quicksort')

        # TS/TV
        sql_query_nb_variant = f"""
            SELECT REF, ALT, count(*) as count FROM {table_variants_from} WHERE len(REF)=1 AND len(ALT)=1 GROUP BY REF, ALT
            ORDER BY REF asc, ALT desc
            """
        tstv = self.conn.execute(sql_query_nb_variant).df()

        # Genotypes
        genotypes = {}
        for sample in self.get_header_sample_list():
            sql_query_genotype = f"""
                SELECT  '{sample}' as sample,
                        REGEXP_EXTRACT("{sample}", '^([0-9/|.]*)[:]*',1) as genotype,
                        count(REGEXP_EXTRACT("{sample}", '^([0-9/|.]*)[:]*',1)) as count,
                        (count(REGEXP_EXTRACT("{sample}", '^([0-9/|.]*)[:]*',1))*100/{nb_of_variants}) || '%' as percentage
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
        log.info("Ti/Ts/Tv Ratio:")
        for d in str(tstv).split("\n"):
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

    def get_input_format(self) -> str:
        """
        It returns the format of the input variable.
        :return: The format is being returned.
        """
        return self.input_format

    def get_output(self) -> str:
        """
        It returns the output of the neuron.
        :return: The output of the neural network.
        """
        return self.output

    def get_output_format(self) -> str:
        """
        It returns the format of the input variable.
        :return: The format is being returned.
        """
        return self.output_format

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
        if type == "vcf":
            return self.header_vcf
        elif type == "list":
            return self.header_list

    def get_header_length(self) -> int:
        """
        This function retruns header length (without #CHROM line)

        :return: The length of the header list.
        """
        return len(self.header_list) - 1

    def get_header_columns(self) -> list:
        """
        This function returns the header list of a VCF

        :return: The length of the header list.
        """
        return self.header_list[-1]

    def get_header_columns_as_sql(self) -> str:
        """
        This function retruns header length (without #CHROM line)

        :return: The length of the header list.
        """
        sql_column_list = []
        for col in self.get_header_columns().strip().split("\t"):
            sql_column_list.append(f"\"{col}\"")
        return ",".join(sql_column_list)

    def get_header_sample_list(self):
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
        return self.connexion_format

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
                else:
                    log.error(
                        f"Failed insert file into table. Unknown connexion format {connexion_format}")
                    raise ValueError(
                        f"Failed insert file into table. Unknown connexion format {connexion_format}")

    def load_data(self) -> None:
        """
        It reads a VCF file and inserts it into a table
        """

        log.info("Loading data...")

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
            sql_create_table_columns.append(f"\"{column}\" {column_type}")
            sql_create_table_columns_list.append(f"\"{column}\"")


        # get table variants
        table_variants = self.get_table_variants()

        # Create database
        log.debug(f"Create Table {table_variants}")
        sql_create_table_columns_sql = ", ".join(sql_create_table_columns)
        sql_create_table_columns_list_sql = ", ".join(
            sql_create_table_columns_list)
        sql_create_table = f"CREATE TABLE IF NOT EXISTS {table_variants} ({sql_create_table_columns_sql})"
        self.conn.execute(sql_create_table)

        # chunksize define length of file chunk load file
        # chunksize = 100000
        # chunksize = 1000000
        chunksize = 100000

        # Access
        access = self.get_config().get("access", None)

        # Load data
        log.debug(f"Load Data from {self.input_format}")
        if self.input_format in ["vcf", "gz", "tsv", "csv", "psv"]:

            # delimiter
            delimiters = {"vcf": "\t", "gz": "\t", "tsv": "\t", "csv": ",", "psv": "|"}
            delimiter = delimiters.get(self.input_format, "\t")

            # Load the input file
            with open(self.input, "rt") as input_file:
                # Use the appropriate file handler based on the input format
                if self.input_format == "gz":
                    input_file = bgzf.open(self.input, "rt")
                if self.input_format == "gz" or self.input_format == "vcf":
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

            if access in ["RO"]:
                # print("NO loading data")
                self.drop_variants_table()
                sql_parquet = f"CREATE VIEW {table_variants} AS SELECT * FROM '{self.input}'"
                #self.conn.execute(sql_view)
            else:
                sql_parquet = f"COPY {table_variants} FROM '{self.input}'"
            
            self.conn.execute(sql_parquet)

        elif self.input_format in ["db", "duckdb"]:
            log.debug(f"Input file format '{self.input_format}' duckDB")

        else:
            log.error(f"Input file format '{self.input_format}' not available")
            raise ValueError(
                f"Input file format '{self.input_format}' not available")

        # Explode INFOS fields into table fields
        if self.get_param().get("explode_infos", None):
            self.explode_infos(
                prefix=self.get_param().get("explode_infos", None))

        # Create index after insertion
        self.create_indexes()


    def explode_infos(self, prefix: str = None, create_index: bool = False) -> None:
        """
        The function takes a VCF file and explodes the INFO fields into individual columns
        """

        # drop indexes
        self.drop_indexes()

        # Access
        access = self.get_config().get("access", None)

        if access not in ["RO"]:

            # prefix
            if not prefix or type(prefix) != str:
                prefix = "INFO/"

            # table variants
            table_variants = self.get_table_variants(clause="select")

            log.debug(
                f"Explode INFO fields - ADD [{len(self.get_header().infos)}] annotations fields")

            sql_info_alter_table_array = []

            for info in self.get_header().infos:
                log.debug(
                    f"Explode INFO fields - ADD {info} annotations fields")

                info_id_sql = prefix+info
                type_sql = self.code_type_map_to_sql.get(
                    self.get_header().infos[info].type, "VARCHAR")
                if self.get_header().infos[info].num != 1:
                    type_sql = "VARCHAR"

                # Add field
                sql_info_alter_table = f"ALTER TABLE {table_variants} ADD COLUMN IF NOT EXISTS \"{info_id_sql}\" {type_sql} DEFAULT null"
                log.debug(
                    f"Explode INFO fields - ADD {info} annotations fields: {sql_info_alter_table}")
                self.conn.execute(sql_info_alter_table)

                # add field to index
                self.index_additionnal_fields.append(info_id_sql)

                # Update field array
                update_info_field = f"\"{info_id_sql}\" = CASE WHEN REGEXP_EXTRACT(INFO, '[^;]*{info}=([^;]*)',1) == '' THEN NULL WHEN REGEXP_EXTRACT(INFO, '{info}=([^;]*)',1) == '.' THEN NULL ELSE REGEXP_EXTRACT(INFO, '{info}=([^;]*)',1) END"
                sql_info_alter_table_array.append(update_info_field)

            # # Update table
            sql_info_alter_table_array_join = ", ".join(
                sql_info_alter_table_array)
            sql_info_alter_table = f"UPDATE {table_variants} SET {sql_info_alter_table_array_join}"
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

        # print("create indexes")

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
                # print(f"create {field}")
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

        if access not in ["RO"]:
            list_indexes = self.conn.execute(
                f"SELECT index_name FROM duckdb_indexes WHERE table_name='{table_variants}'")
            for index in list_indexes.df()["index_name"]:
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
            # print(line)
            header_list.append(line)
            if line.startswith('#CHROM'):
                break
        return header_list

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

    def export_output(self, export_header: bool = True) -> None:
        """
        It takes a VCF file, and outputs a VCF file

        :param export_header: If True, the header will be exported to a file. If False, the header will
        be exported to a temporary file, defaults to True (optional)
        """

        # Export  header
        if export_header:
            header_name = self.export_header()
        else:
            # Header
            tmp_header = NamedTemporaryFile(
                prefix=self.get_prefix(), dir=self.get_tmp_dir())
            tmp_header_name = tmp_header.name
            f = open(tmp_header_name, 'w')
            vcf.Writer(f, self.header_vcf)
            f.close()
            header_name = tmp_header_name

        if self.get_output():

            output_file = self.get_output()
            sql_columns = self.get_header_columns_as_sql()
            table_variants = self.get_table_variants()
            sql_query_hard = ""
            sql_query_sort = ""
            sql_query_limit = ""

            # delimiter
            delimiters = {"vcf": "\t", "gz": "\t", "tsv": "\t", "csv": ",", "psv": "|"}
            delimiter = delimiters.get(self.get_output_format(), "\t")

            # Threads
            threads = self.get_threads()

            log.debug(f"Export file: {output_file}")

            # Extra columns
            sql_extra_columns = ""
            if self.get_param().get("export_extra_infos", None):
                sql_extra_columns = self.get_extra_infos_sql()

            log.debug(f"Export extra columns: {sql_extra_columns}")

            if self.get_output_format() in ["parquet"]:

                # Export parquet
                sql_query_export = f"COPY (SELECT {sql_columns} {sql_extra_columns} FROM {table_variants} WHERE 1 {sql_query_hard} {sql_query_sort} {sql_query_limit}) TO '{output_file}' WITH (FORMAT PARQUET)"
                self.conn.execute(sql_query_export)

            # elif self.get_output_format() in ["db", "duckdb"]:

            #     log.debug("Export in DuckDB. Nothing to do")

            elif self.get_output_format() in ["tsv", "csv", "psv"]:

                # Export TSV/CSV
                sql_query_export = f"COPY (SELECT {sql_columns} {sql_extra_columns} FROM {table_variants} WHERE 1 {sql_query_hard} {sql_query_sort} {sql_query_limit}) TO '{output_file}' WITH (FORMAT CSV, DELIMITER '{delimiter}', HEADER)"
                self.conn.execute(sql_query_export)

            elif self.get_output_format() in ["vcf", "gz"]:

                # Extract VCF
                # Variants
                tmp_variants = NamedTemporaryFile(prefix=self.get_prefix(
                ), dir=self.get_tmp_dir(), suffix=".gz", delete=False)
                tmp_variants_name = tmp_variants.name
                sql_query_export = f"COPY (SELECT {sql_columns} FROM {table_variants} WHERE 1 {sql_query_hard} {sql_query_sort} {sql_query_limit}) TO '{tmp_variants_name}' WITH (FORMAT CSV, DELIMITER '\t', HEADER, QUOTE '')"
                self.conn.execute(sql_query_export)

                # VCF
                command = f"grep '^#CHROM' -v {header_name} > {output_file}.vcf; cat {tmp_variants_name} >> {output_file}.vcf"
                subprocess.run(command, shell=True)

                if self.get_output_format() in ["vcf"]:
                    command = f"mv {output_file}.vcf {output_file}"
                    subprocess.run(command, shell=True)
                elif self.get_output_format() in ["gz"]:
                    bgzip_command = get_bgzip(threads=threads)
                    command = f""" {bgzip_command} {output_file}.vcf > {output_file} && rm {output_file}.vcf """
                    subprocess.run(command, shell=True)


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
        table_columns = self.conn.execute(query).df().columns.tolist()
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

    def export_header(self, header_name: str = None) -> str:
        """
        It takes a VCF file, and writes the header to a new file

        :param header_name: the name of the header file to be created. If not specified, the header will
        be written to the output file
        :return: The name of the temporary header file.
        """

        if not header_name:
            output_file = self.get_output()
        tmp_header_name = output_file + ".hdr"
        f = open(tmp_header_name, 'w')
        vcf.Writer(f, self.get_header())
        f.close()
        return tmp_header_name

    def export_variant_vcf(self, vcf_file, file_type: str = "vcf", remove_info: bool = False, add_samples: bool = True, compression: int = 1, index: bool = False) -> None:
        """
        It takes a VCF file and a list of samples, and returns a VCF file with only the samples in the
        list

        :param vcf_file: the name of the file to write to
        :param file_type: The file type of the output file. Can be "vcf" or "gz", defaults to vcf
        (optional)
        :param remove_info: If you want to remove the INFO field, set this to True. If you want to
        remove a specific INFO field, set this to the name of the INFO field, defaults to False
        (optional)
        :param add_samples: If True, the samples will be added to the VCF. If False, the samples will be
        removed, defaults to True (optional)
        :param compression: 1-9, 1 being the fastest and 9 being the most compressed, defaults to 1
        (optional)
        :param index: If True, the output VCF will be indexed with tabix, defaults to False (optional)
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
        if add_samples:
            self.get_header_sample_list()
            samples_fields = " , FORMAT , " + \
                " , ".join(self.get_header_sample_list())
        else:
            samples_fields = ""

        # Header
        tmp_header = NamedTemporaryFile(
            prefix=self.get_prefix(), dir=self.get_tmp_dir(), delete=False)
        tmp_header_name = tmp_header.name
        f = open(tmp_header_name, 'w')
        vcf.Writer(f, self.header_vcf)
        f.close()

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
            with gzip.open(tmp_variants_name, 'wt') as f:
                writer = csv.writer(f, delimiter='\t',
                                    quotechar='', quoting=csv.QUOTE_NONE)
                cursor = self.conn.execute(sql_query_select)
                writer.writerow([i[0] for i in cursor.description])
                writer.writerows(cursor)

        # Create output
        # Cat header and variants
        command_gzip = ""
        if file_type in ["gz"]:
            command_gzip = f" | bgzip -l {compression} -c "
        if index:
            command_tabix = f" && tabix {vcf_file}"
        else:
            command_tabix = ""

        command = f"grep '^#CHROM' -v {tmp_header_name} {command_gzip} > {vcf_file}; gzip -dc {tmp_variants_name} {command_gzip} >> {vcf_file} {command_tabix}"

        subprocess.run(command, shell=True)

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
        return int(self.config.get("threads", 1))

    def annotation(self) -> None:
        """
        It annotates the VCF file with the annotations specified in the config file.
        """

        param = self.get_param()

        if param.get("annotations"):
            if not "annotation" in param:
                param["annotation"] = {}
            for annotation_file in param.get("annotations"):
                annotations = param.get("annotations").get(
                    annotation_file, None)
                if annotation_file == "snpeff":
                    if "snpeff" not in param["annotation"]:
                        param["annotation"]["snpeff"] = {}
                    if "options" not in param["annotation"]["snpeff"]:
                        param["annotation"]["snpeff"]["options"] = ""
                if annotation_file.startswith("annovar"):
                    if "annovar" not in param["annotation"]:
                        param["annotation"]["annovar"] = {}
                    if "annotations" not in param["annotation"]["annovar"]:
                        param["annotation"]["annovar"]["annotations"] = {}
                    annotation_file_split = annotation_file.split(":")
                    if len(annotation_file_split) > 1:
                        annotation_file_annotation = annotation_file_split[1]
                        param["annotation"]["annovar"]["annotations"][annotation_file_annotation] = annotations
                elif os.path.exists(annotation_file):
                    log.debug(f"Quick Annotation File {annotation_file}")
                    quick_annotation_file = annotation_file
                    quick_annotation_name, quick_annotation_extension = os.path.splitext(
                        annotation_file)
                    quick_annotation_format = quick_annotation_extension.replace(
                        ".", "")
                    # print(quick_annotation_format)
                    if quick_annotation_format in ["gz"]:
                        quick_annotation_format_name, quick_annotation_format_extension = os.path.splitext(
                            quick_annotation_name)
                        quick_annotation_type = quick_annotation_format_extension.replace(
                            ".", "")
                    else:
                        quick_annotation_type = quick_annotation_format
                    # print(quick_annotation_type)
                    format = None
                    if quick_annotation_type in ["parquet", "duckdb"]:
                        format = "parquet"
                    elif quick_annotation_type in ["vcf", "bed"]:
                        format = "bcftools"
                    else:
                        log.error(
                            f"Quick Annotation File {quick_annotation_file} - format {quick_annotation_type} not supported yet")
                        raise ValueError(
                            f"Quick Annotation File {quick_annotation_file} - format {quick_annotation_type} not supported yet"
                        )
                    if format:
                        if format not in param["annotation"]:
                            param["annotation"][format] = {}
                        if "annotations" not in param["annotation"][format]:
                            param["annotation"][format]["annotations"] = {}
                        param["annotation"][format]["annotations"][quick_annotation_file] = annotations
                else:
                    log.error(
                        f"Quick Annotation File {annotation_file} does NOT exist")

            self.set_param(param)

        # print(param)
        # return

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
            if param.get("annotation", {}).get("varank", None):
                log.info("Annotations 'varank'...")

        # Explode INFOS fields into table fields
        if self.get_param().get("explode_infos", None):
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
        databases_folders = self.config.get("folders", {}).get(
            "databases", {}).get("bcftools", ["."])
        log.debug("Databases annotations: " + str(databases_folders))

        # Param
        annotations = self.param.get("annotation", {}).get(
            "bcftools", {}).get("annotations", None)
        log.debug("Annotations: " + str(annotations))

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

                # Find vcf/bed file and header file
                db_file = None
                db_hdr_file = None
                for databases_folder in databases_folders:
                    db_file = None
                    db_hdr_file = None
                    log.debug("Annotation file check: " + annotation + " or " +
                              str(databases_folder+"/"+annotation+".{vcf,bed}"))

                    # VCF .vcf BED .bed
                    if os.path.exists(annotation):
                        db_file = annotation
                    elif os.path.exists(databases_folder+"/"+annotation+".vcf"):
                        db_file = databases_folder+"/"+annotation+".vcf"
                    elif os.path.exists(databases_folder+"/"+annotation+".vcf.gz"):
                        db_file = databases_folder+"/"+annotation+".vcf.gz"
                    # BED .bed
                    if os.path.exists(annotation):
                        db_file = annotation
                    elif os.path.exists(databases_folder+"/"+annotation+".bed"):
                        db_file = databases_folder+"/"+annotation+".bed"
                    elif os.path.exists(databases_folder+"/"+annotation+".bed.gz"):
                        db_file = databases_folder+"/"+annotation+".bed.gz"
                    if not db_file:
                        continue

                    # Header .hdr
                    if os.path.exists(db_file+".hdr"):
                        db_hdr_file = db_file+".hdr"

                    # parquet and hdr found
                    if db_file and db_hdr_file:
                        break

                # Database format and type
                db_file_name, db_file_extension = os.path.splitext(db_file)
                db_file_format = db_file_extension.replace(".", "")
                if db_file_format in ["gz"]:
                    db_file_format_name, db_file_format_extension = os.path.splitext(
                        db_file_name)
                    db_file_type = db_file_format_extension.replace(".", "")
                else:
                    db_file_type = db_file_format

                # try to extract header
                if db_file_type in ["vcf", "bed"] and not db_hdr_file:
                    log.debug(f"Try to extract header of file {db_file}")
                    tmp_extract_header = NamedTemporaryFile(prefix=self.get_prefix(
                    ), dir=self.get_tmp_dir(), suffix=".hdr", delete=False)
                    tmp_extract_header_name = tmp_extract_header.name
                    tmp_files.append(tmp_extract_header_name)
                    command_extract_header = f"bcftools view -h {db_file} > {tmp_extract_header_name} 2>/dev/null"
                    run_parallel_commands([command_extract_header], threads)
                    db_hdr_file = tmp_extract_header_name

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

                        # Protect header for bcftools (remove "#CHROM" line)
                        log.debug(
                            "Protect Header file - remove #CHROM line if exists")
                        tmp_header_vcf = NamedTemporaryFile(prefix=self.get_prefix(
                        ), dir=self.get_tmp_dir(), suffix=".hdr", delete=False)
                        tmp_header_vcf_name = tmp_header_vcf.name
                        tmp_files.append(tmp_header_vcf_name)
                        run_parallel_commands(
                            [f"grep '^#CHROM' -v {db_hdr_file} > {tmp_header_vcf_name}"], 1)

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

                        if True:

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
                                # print(sql_query_intervals_for_bed)
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
                        log.info(message)
                    # debug info
                    for message in list(set(error_message_command_all)):
                        log.debug(message)
                    # failed
                    if len(error_message_command_err):
                        log.error("Annotation failed: Error in commands")
                        raise ValueError(
                            "Annotation failed: Error in commands")

                    # Update variants
                    log.info(f"Annotation - Updating...")
                    self.update_from_vcf(tmp_annotate_vcf_name)

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
            "databases", {}).get("snpeff", "")

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
                log.warning(f"Annotation warning: snpEff jar bin '{java_bin}'")

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

        if not os.path.exists(snpeff_databases):
            log.warning(
                f"Annotation warning: no snpEff database '{snpeff_databases}'. Try to find...")
            # Try to find snpeff database
            try:
                snpeff_databases = os.path.dirname(os.path.dirname(
                    find_all('snpEffectPredictor.bin', '/')[0]))
            except:
                log.error(
                    f"Annotation failed: no snpEff database '{snpeff_databases}'")
                raise ValueError(
                    f"Annotation failed: no snpEff database '{snpeff_databases}'")
            if not os.path.exists(snpeff_databases):
                log.error(
                    f"Annotation failed: no snpEff database '{snpeff_databases}'")
                raise ValueError(
                    f"Annotation failed: no snpEff database '{snpeff_databases}'")
            else:
                log.warning(
                    f"Annotation warning: snpEff database found '{snpeff_databases}'")

        # Param
        param = self.get_param()
        log.debug("Param: " + str(param))

        # Param
        options = param.get("annotation", {}).get(
            "snpeff", {}).get("options", None)
        log.debug("Options: " + str(options))

        # Param - Assembly
        assembly = param.get("assembly", "hg19")

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
        if not self.conn.execute(f"{sql_query_chromosomes}").df()["count"][0]:
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
                log.info(message)
            # debug info
            for message in list(set(error_message_command_all)):
                log.debug(message)
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
            "databases", {}).get("annovar", "")

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

        if not os.path.exists(annovar_databases):
            log.warning(
                f"Annotation warning: no annovar database '{annovar_databases}'. Try to find...")
            # Try to find annovar database
            try:
                annovar_databases = os.path.dirname(
                    find_all('annovar_downdb.log', '/')[0])
            except:
                log.error(
                    f"Annotation failed: no annovar database '{annovar_databases}'")
                raise ValueError(
                    f"Annotation failed: no annovar database '{annovar_databases}'")
            if not os.path.exists(annovar_databases):
                log.error(
                    f"Annotation failed: no annovar database '{annovar_databases}'")
                raise ValueError(
                    f"Annotation failed: no annovar database '{annovar_databases}'")
            else:
                log.warning(
                    f"Annotation warning: annovar database found '{annovar_databases}'")

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
        assembly = param.get("assembly", "hg19")

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

            for annotation in annotations:
                annotation_fields = annotations[annotation]

                if not annotation_fields:
                    annotation_fields = {"INFO": None}

                log.debug(f"Annotation '{annotation}'")
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
                    operation = "gx"
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
                command_annovar = f"""{annovar_bin} {tmp_vcf_name} {annovar_databases} --buildver {assembly} --outfile {tmp_annotate_vcf_prefix} --remove --protocol {protocol} --operation {operation} {argument_option} {command_options} 2>>{tmp_annotate_vcf_name_err} && mv {tmp_annotate_vcf_name_annovar} {tmp_annotate_vcf_name}.tmp.vcf """
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
                if True:
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
                        log.info(message)
                    # debug info
                    for message in list(set(error_message_command_all)):
                        log.debug(message)
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
                log.info(f"Annotation - Annotation merging " +
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
                log.info(f"Annotation - Updating...")
                self.update_from_vcf(tmp_annotate_vcf_name)

            # Clean files
            # Tmp file remove command
            if True:
                tmp_files_remove_command = ""
                if tmp_files:
                    tmp_files_remove_command = " ".join(tmp_files)
                clean_command = f" rm -f {tmp_files_remove_command} "
                log.info(f"Annotation - Annotation cleaning ")
                log.debug(f"Annotation - cleaning command: {clean_command}")
                run_parallel_commands([clean_command], 1)

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
        databases_folders = self.config.get("folders", {}).get(
            "databases", {}).get("parquet", ["."])
        log.debug("Databases annotations: " + str(databases_folders))

        # Param
        annotations = self.param.get("annotation", {}).get(
            "parquet", {}).get("annotations", None)
        log.debug("Annotations: " + str(annotations))

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

        # explode infos
        self.explode_infos(prefix=self.get_param().get("explode_infos", None))

        # drop indexes
        log.debug(f"Drop indexes...")
        self.drop_indexes()

        if annotations:

            for annotation in annotations:
                annotation_fields = annotations[annotation]

                if not annotation_fields:
                    annotation_fields = {"INFO": None}

                log.debug(f"Annotation '{annotation}'")
                log.debug(
                    f"Annotation '{annotation}' - fields: {annotation_fields}")

                # Find parquet file and header file
                parquet_file = None
                parquet_hdr_file = None
                for databases_folder in databases_folders:
                    parquet_file = None
                    parquet_hdr_file = None
                    log.debug("Annotation file check: " + annotation +
                              " or " + str(databases_folder+"/"+annotation+".parquet"))

                    # Parquet .parquet
                    if os.path.exists(annotation):
                        parquet_file = annotation
                    elif os.path.exists(databases_folder+"/"+annotation+".parquet"):
                        parquet_file = databases_folder+"/"+annotation+".parquet"
                    if not parquet_file:
                        continue

                    # Header .hdr
                    if os.path.exists(parquet_file+".hdr"):
                        parquet_hdr_file = parquet_file+".hdr"

                    # parquet and hdr found
                    if parquet_file and parquet_hdr_file:
                        break

                if not parquet_file or not parquet_hdr_file:
                    log.error("Annotation failed: file not found")
                    raise ValueError("Annotation failed: file not found")
                else:

                    parquet_file_link = f"'{parquet_file}'"

                    # Database format and type
                    parquet_file_name, parquet_file_extension = os.path.splitext(
                        parquet_file)
                    parquet_file_basename = os.path.basename(parquet_file)
                    parquet_file_format = parquet_file_extension.replace(
                        ".", "")
                    parquet_file_type = parquet_file_format

                    if parquet_file_format in ["db", "duckdb", "sqlite"]:
                        parquet_file_as_duckdb_name = parquet_file_basename.replace(
                            ".", "_")
                        if parquet_file_format in ["sqlite"]:
                            parquet_file_format_attached_type = ", TYPE SQLITE"
                        else:
                            parquet_file_format_attached_type = ""
                        log.debug(
                            f"Annotation '{annotation}' - attach database : " + str(parquet_file))
                        self.conn.execute(
                            f"ATTACH DATABASE '{parquet_file}' AS {parquet_file_as_duckdb_name} (READ_ONLY{parquet_file_format_attached_type})")
                        parquet_file_link = f"{parquet_file_as_duckdb_name}.variants"
                    elif parquet_file_format in ["parquet"]:
                        parquet_file_link = f"'{parquet_file}'"

                    log.debug(f"Annotation '{annotation}' - file: " +
                              str(parquet_file) + " and " + str(parquet_hdr_file))

                    # Load header as VCF object
                    parquet_hdr_vcf = Variants(input=parquet_hdr_file)
                    parquet_hdr_vcf_header_infos = parquet_hdr_vcf.get_header().infos
                    log.debug("Annotation database header: " +
                              str(parquet_hdr_vcf_header_infos))

                    # get extra infos
                    parquet_columns = self.get_extra_infos(
                        table=parquet_file_link)
                    # print(parquet_columns)

                    # For all fields in database
                    annotation_fields_ALL = False
                    if "ALL" in annotation_fields or "INFO" in annotation_fields:
                        annotation_fields_ALL = True
                        annotation_fields = {
                            key: key for key in parquet_hdr_vcf_header_infos}
                        log.debug(
                            "Annotation database header - All annotations added: " + str(annotation_fields))

                    # List of annotation fields to use
                    sql_query_annotation_update_info_sets = []

                    # Number of fields
                    nb_annotation_field = 0

                    # Annotation fields processed
                    annotation_fields_processed = []

                    for annotation_field in annotation_fields:

                        # annotation_field_column
                        if annotation_field in parquet_columns:
                            annotation_field_column = annotation_field
                        elif "INFO/" + annotation_field in parquet_columns:
                            annotation_field_column = "INFO/" + annotation_field
                        else:
                            annotation_field_column = "INFO"

                        # field new name, if parametered
                        annotation_fields_new_name = annotation_fields.get(
                            annotation_field, annotation_field)
                        if not annotation_fields_new_name:
                            annotation_fields_new_name = annotation_field

                        # check annotation field in data
                        annotation_field_exists_on_variants = 0
                        if annotation_fields_new_name not in self.get_header().infos:
                            sampling_annotation_field_exists_on_variants = 10000
                            sql_query_chromosomes = f"""
                                SELECT 1 AS count
                                FROM (SELECT * FROM {table_variants} as table_variants LIMIT {sampling_annotation_field_exists_on_variants})
                                WHERE ';' || INFO LIKE '%;{annotation_fields_new_name}=%'
                                LIMIT 1
                                """
                            annotation_field_exists_on_variants = len(
                                self.conn.execute(f"{sql_query_chromosomes}").df()["count"])
                            log.debug(f"Annotation field {annotation_fields_new_name} found in variants: " + str(
                                annotation_field_exists_on_variants))

                        # To annotate
                        force_update_annotation = False
                        if annotation_field in parquet_hdr_vcf.get_header().infos and (force_update_annotation or (annotation_fields_new_name not in self.get_header().infos and not annotation_field_exists_on_variants)):

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
                                f"Annotation '{annotation}' - '{annotation_field}' -> 'INFO/{annotation_fields_new_name}'")

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
                            if annotation_field_column == "INFO":
                                sql_query_annotation_update_info_sets.append(f"""
                                || CASE WHEN REGEXP_EXTRACT(';' || table_parquet.INFO, ';{annotation_field}=([^;]*)',1) NOT IN ('','.')
                                        THEN '{annotation_field_sep}' || '{annotation_fields_new_name}=' || REGEXP_EXTRACT(';' || table_parquet.INFO, ';{annotation_field}=([^;]*)',1)
                                        ELSE ''
                                    END
                                """)
                            # Found in a specific column
                            else:
                                sql_query_annotation_update_info_sets.append(f"""
                                || CASE WHEN table_parquet."{annotation_field_column}" NOT IN ('','.')
                                        THEN '{annotation_field_sep}' || '{annotation_fields_new_name}=' || table_parquet."{annotation_field_column}"
                                        ELSE ''
                                    END
                                """)

                        # Not to annotate
                        else:

                            if force_update_annotation:
                                annotation_message = "forced"
                            else:
                                annotation_message = "skipped"

                            if annotation_field not in parquet_hdr_vcf.get_header().infos:
                                log.warning(
                                    f"Annotation '{annotation}' - '{annotation_field}' [{nb_annotation_field}] - not available in parquet file")
                            if annotation_fields_new_name in self.get_header().infos:
                                log.warning(
                                    f"Annotation '{annotation}' - '{annotation_fields_new_name}' [{nb_annotation_field}] - already exists in header ({annotation_message})")
                            if annotation_field_exists_on_variants:
                                log.warning(
                                    f"Annotation '{annotation}' - '{annotation_fields_new_name}' [{nb_annotation_field}] - already exists in variants ({annotation_message})")

                    # Check if ALL fields have to be annotated. Thus concat all INFO field
                    allow_annotation_full_info = True
                    if allow_annotation_full_info and nb_annotation_field == len(annotation_fields) and annotation_fields_ALL:
                        sql_query_annotation_update_info_sets = []
                        sql_query_annotation_update_info_sets.append(
                            f"|| table_parquet.INFO ")

                    if sql_query_annotation_update_info_sets:

                        # Annotate
                        log.info(f"Annotation '{annotation}' - Annotation...")

                        # Join query annotation update info sets for SQL
                        sql_query_annotation_update_info_sets_sql = " ".join(
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
                                f"Annotation '{annotation}' - Chromosome '{chrom}' - Start Autodetection Intervals...")

                            batch_annotation_databases_step = None
                            batch_annotation_databases_ncuts = 1

                            # Create intervals from 0 to max position variant, with the batch window previously defined
                            sql_query_intervals = split_interval(
                                sql_query_chromosomes_max_pos_dictionary_min_pos, sql_query_chromosomes_max_pos_dictionary_max_pos, step=batch_annotation_databases_step, ncuts=batch_annotation_databases_ncuts)

                            log.debug(
                                f"Annotation '{annotation}' - Chromosome '{chrom}' - Stop Autodetection Intervals")

                            # Interval Start/Stop
                            sql_query_interval_start = sql_query_intervals[0]

                            # For each interval
                            for i in sql_query_intervals[1:]:

                                # Interval Start/Stop
                                sql_query_interval_stop = i

                                log.debug(
                                    f"Annotation '{annotation}' - Chromosome '{chrom}' - Interval [{sql_query_interval_start}-{sql_query_interval_stop}] ...")

                                log.debug(
                                    f"Annotation '{annotation}' - Chromosome '{chrom}' - Interval [{sql_query_interval_start}-{sql_query_interval_stop}] - Start detecting regions...")

                                regions = [
                                    (chrom, sql_query_interval_start, sql_query_interval_stop)]

                                log.debug(
                                    f"Annotation '{annotation}' - Chromosome '{chrom}' - Interval [{sql_query_interval_start}-{sql_query_interval_stop}] - Stop detecting regions")

                                # Fusion des régions chevauchantes
                                if regions:

                                    # Number of regions
                                    nb_regions = len(regions)

                                    # create where caluse on regions
                                    clause_where_regions_variants = create_where_clause(
                                        regions, table="table_variants")
                                    clause_where_regions_parquet = create_where_clause(
                                        regions, table="table_parquet")

                                    log.debug(
                                        f"Annotation '{annotation}' - Chromosome '{chrom}' - Interval [{sql_query_interval_start}-{sql_query_interval_stop}] - {nb_regions} regions...")

                                    sql_query_annotation_chrom_interval_pos = f"""
                                        UPDATE {table_variants} as table_variants
                                            SET INFO = CASE WHEN table_variants.INFO NOT IN ('','.') THEN table_variants.INFO ELSE '' END || CASE WHEN table_variants.INFO NOT IN ('','.') AND ('' {sql_query_annotation_update_info_sets_sql}) NOT IN ('','.') THEN ';' ELSE '' END {sql_query_annotation_update_info_sets_sql}
                                            FROM {parquet_file_link} as table_parquet
                                            WHERE ( {clause_where_regions_parquet} )
                                                AND table_parquet.\"#CHROM\" = table_variants.\"#CHROM\"
                                                AND table_parquet.\"POS\" = table_variants.\"POS\"
                                                AND table_parquet.\"ALT\" = table_variants.\"ALT\"
                                                AND table_parquet.\"REF\" = table_variants.\"REF\";
                                                """
                                    query_dict[f"{chrom}:{sql_query_interval_start}-{sql_query_interval_stop}"] = sql_query_annotation_chrom_interval_pos

                                    log.debug(
                                        "Create SQL query: " + str(sql_query_annotation_chrom_interval_pos))

                                    # Interval Start/Stop
                                    sql_query_interval_start = sql_query_interval_stop

                            # nb_of_variant_annotated
                            nb_of_variant_annotated += nb_of_variant_annotated_by_chrom

                        nb_of_query = len(query_dict)
                        num_query = 0
                        for query_name in query_dict:
                            query = query_dict[query_name]
                            num_query += 1
                            log.info(
                                f"Annotation '{annotation}' - Annotation - Query [{num_query}/{nb_of_query}] {query_name}...")
                            result = self.conn.execute(query)
                            nb_of_variant_annotated_by_query = result.df()[
                                "Count"][0]
                            nb_of_variant_annotated += nb_of_variant_annotated_by_query
                            log.info(
                                f"Annotation '{annotation}' - Annotation - Query [{num_query}/{nb_of_query}] {query_name} - {nb_of_variant_annotated_by_query} variants annotated")

                        log.info(
                            f"Annotation '{annotation}' - Annotation of {nb_of_variant_annotated} variants out of {nb_variants} (with {nb_of_query} queries)")

                    else:

                        log.info(
                            f"Annotation '{annotation}' - No Annotations available")

                    log.debug("Final header: " + str(vcf_reader.infos))

    def update_from_vcf(self, vcf_file: str) -> None:
        """
        > If the database is duckdb, then use the parquet method, otherwise use the sqlite method

        :param vcf_file: the path to the VCF file
        """

        connexion_format = self.get_connexion_format()

        if connexion_format in ["duckdb"]:
            self.update_from_vcf_parquet(vcf_file)
        elif connexion_format in ["sqlite"]:
            self.update_from_vcf_sqlite(vcf_file)

    def update_from_vcf_parquet(self, vcf_file: str) -> None:
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
        mergeVCF = Variants(
            None, vcf_file, tmp_merged_vcf_name, config=self.get_config())
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
                SET INFO = (
                            SELECT CASE WHEN table_variants.INFO NOT IN ('','.') THEN table_variants.INFO ELSE '' END || CASE WHEN table_variants.INFO NOT IN ('','.') AND table_parquet.INFO NOT IN ('','.')  THEN ';' ELSE '' END || CASE WHEN table_parquet.INFO NOT IN ('','.') THEN table_parquet.INFO ELSE '' END
                            FROM '{tmp_merged_vcf_name}' as table_parquet
                                    WHERE table_parquet.\"#CHROM\" = table_variants.\"#CHROM\"
                                    AND table_parquet.\"POS\" = table_variants.\"POS\"
                                    AND table_parquet.\"ALT\" = table_variants.\"ALT\"
                                    AND table_parquet.\"REF\" = table_variants.\"REF\"
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
        sql_query_update = f"""
            UPDATE variants as table_variants
            SET INFO = (
                SELECT CASE 
                    WHEN table_variants.INFO NOT IN ('','.') 
                    THEN table_variants.INFO 
                    ELSE '' 
                END || 
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
        sql_table_variants = f"DROP TABLE {table_variants}"
        self.conn.execute(sql_table_variants)

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
        profiles = param.get("prioritization", {}).get("profiles", ["default"])
        pzfields = param.get("prioritization", {}).get(
            "pzfields", ["PZFlag", "PZScore"])
        default_profile = param.get("prioritization", {}).get(
            "default_profile", ["default"])
        pzfields_sep = param.get("prioritization", {}).get("pzfields_sep", "_")
        prioritization_score_mode = param.get("prioritization", {}).get(
            "prioritization_score_mode", "HOWARD")

        # Profiles are in files
        if os.path.exists(config_profiles):
            with open(config_profiles) as profiles_file:
                config_profiles = json.load(profiles_file)

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
                    self.conn.execute(
                        f"ALTER TABLE {table_variants} ADD COLUMN IF NOT EXISTS {pzfield} INTEGER DEFAULT 0")
                elif re.match("PZFlag.*", pzfield):
                    self.conn.execute(
                        f"ALTER TABLE {table_variants} ADD COLUMN IF NOT EXISTS {pzfield} BOOLEAN DEFAULT 1")
                else:
                    self.conn.execute(
                        f"ALTER TABLE {table_variants} ADD COLUMN IF NOT EXISTS {pzfield} STRING DEFAULT ''")

            # Profiles
            if profiles:

                # foreach profile in configuration file
                for profile in config_profiles:

                    # If profile is asked in param, or ALL are asked (empty profile [])
                    if profile in profiles or profiles == []:
                        log.info(f"Profile '{profile}'")

                        sql_queries = []
                        for annotation in config_profiles[profile]:

                            # Check if annotation field is present
                            if not explode_infos_prefix+annotation in extra_infos:
                                log.debug(
                                    f"Annotation '{annotation}' not in data")
                                continue

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
                                    sql_set_info.append(
                                        f" || 'PZScore{pzfields_sep}{profile}=' || PZScore{pzfields_sep}{profile} ")
                                    if profile == default_profile and "PZScore" in list_of_pzfields:
                                        sql_set_info.append(
                                            f" || 'PZScore=' || PZScore{pzfields_sep}{profile} ")
                                if f"PZFlag{pzfields_sep}{profile}" in list_of_pzfields:
                                    sql_set.append(
                                        f"PZFlag{pzfields_sep}{profile} = PZFlag{pzfields_sep}{profile} AND {criterion_flag_bool}")
                                    sql_set_info.append(
                                        f" || 'PZFlag{pzfields_sep}{profile}=' || CASE WHEN PZFlag{pzfields_sep}{profile}==1 THEN 'PASS' WHEN PZFlag{pzfields_sep}{profile}==0 THEN 'FILTERED' END ")
                                    if profile == default_profile and "PZFlag" in list_of_pzfields:
                                        sql_set_info.append(
                                            f" || 'PZFlag=' || PZFlag{pzfields_sep}{profile} ")
                                if f"PZComment{pzfields_sep}{profile}" in list_of_pzfields:
                                    sql_set.append(
                                        f"PZComment{pzfields_sep}{profile} = PZComment{pzfields_sep}{profile} || CASE WHEN PZComment{pzfields_sep}{profile}!='' THEN ', ' ELSE '' END || '{criterion_comment}'")
                                    sql_set_info.append(
                                        f" || 'PZComment{pzfields_sep}{profile}=' || PZComment{pzfields_sep}{profile}  ")
                                    if profile == default_profile and "PZComment" in list_of_pzfields:
                                        sql_set_info.append(
                                            f" || 'PZComment=' || PZComment{pzfields_sep}{profile} ")
                                if f"PZInfos{pzfields_sep}{profile}" in list_of_pzfields:
                                    sql_set.append(
                                        f"PZInfos{pzfields_sep}{profile} = PZInfos{pzfields_sep}{profile} || '{criterion_infos}'")
                                    sql_set_info.append(
                                        f" || 'PZInfos{pzfields_sep}{profile}=' || PZInfos{pzfields_sep}{profile}")
                                    if profile == default_profile and "PZInfos" in list_of_pzfields:
                                        sql_set_info.append(
                                            f" || 'PZInfos=' || PZInfos{pzfields_sep}{profile} ")
                                sql_set_option = ",".join(sql_set)
                                sql_set_info_option = " || ';' ".join(
                                    sql_set_info)

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
                                self.conn.execute(sql_query)

                            log.info(f"""Profile '{profile}' - Update... """)
                            sql_query_update = f"""
                                UPDATE {table_variants}
                                SET INFO = 
                                    CASE WHEN INFO NOT IN ('','.') THEN INFO || ';' ELSE '' END
                                    {sql_set_info_option}

                            """
                            self.conn.execute(sql_query_update)

        else:

            log.warning(f"No profiles in parameters")
