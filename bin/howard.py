import io
import multiprocessing
import os
import re
import subprocess
from tempfile import NamedTemporaryFile
import tempfile
import duckdb
import json
import argparse
import Bio.bgzf as bgzf
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import vcf
import logging as log


### fonctions


### logging

def set_log_level(verbosity):
    configs = {
        "debug": log.DEBUG,
        "info": log.INFO,
        "warning": log.WARNING,
        "error": log.ERROR,
        "critical": log.CRITICAL,
    }
    if verbosity not in configs.keys():
        raise ValueError(
            "Unknown verbosity level:" + verbosity + "\nPlease use any in:" + configs.keys()
        )
    log.basicConfig(
        format="#[%(asctime)s] [%(levelname)s] %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        level=configs[verbosity],
    )

def split_interval(start, end, step=None, ncuts=None):
    if step is None and ncuts is None:
        raise ValueError("Either step or ncuts must be provided")
    if step is not None and ncuts is not None:
        raise ValueError("Only one of step or ncuts must be provided")
    if step is not None:
        return list(range(start, end, step)) + [end]
    if ncuts is not None:
        step = (end - start) / ncuts
        return [start + i*step for i in range(ncuts+1)] + [end]

### Run parallel bash commands

def run_parallel_commands(commands, threads):
    pool = multiprocessing.Pool(threads)
    results = []
    for cmd in commands:
        results.append(pool.apply_async(command, args=(cmd,), error_callback=lambda e: print(e)))
    pool.close()
    pool.join()
    return results


def command(command):
    return subprocess.Popen(command, shell=True).wait()



### Run parallel python functions

def run_parallel_functions(functions, threads):
    pool = multiprocessing.Pool(threads)
    results = []
    for func in functions:
        results.append(pool.apply_async(func, args=(1, "hello"), error_callback=lambda e: print(e)))
    pool.close()
    pool.join()
    return results


# Functions to call for parallel processes


def function_query(obj,query):
    print(query)
    return obj.execute_query(query)




class VCFDataObject:
    def __init__(self, conn=None, input=None, output=None, config={}, param={}):
        """
        The function `__init__` initializes the variables, sets the input, output, config, param, connexion and
        header
        
        :param input: the input file
        :param output: the output file
        :param config: a dictionary containing the configuration of the model
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

    # INIT section

    def set_input(self, input):
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

    def set_config(self, config):
        """
        This function takes in a config object and sets it as the config object for the class
        
        :param config: The configuration object
        """
        self.config = config

    def set_param(self, param):
        """
        This function takes in a param object and sets it as the param object for the class
        
        :param param: The paramters object
        """
        self.param = param

    def init_variables(self):
        """
        This function initializes the variables that will be used in the rest of the class
        """
        self.prefix = "howard"
        self.table_variants = "variants"
        self.dataframe = None

        self.code_type_map = {
            "Integer": 0,
            "String": 1,
            "Float": 2,
            "Flag": 3
        }


    # SET section

    def set_connexion(self, conn):
        """
        It creates a connection to the database
        """
        if not conn:
            conn = duckdb.connect(":memory:")
        self.conn = conn #duckdb.connect(self.get_connexion_db())

    def set_output(self, output=None):
        """
        If the config file has an output key, set the output to the value of that key. Otherwise, set
        the output to the input
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

    def set_header(self):
        """
        It reads the header of a VCF file and stores it as a list of strings and as a VCF object
        """
        input = self.get_input()
        config = self.get_config()
        header_list = []
        if self.input_format in ["gz"]:
            with bgzf.open(input, 'rt') as f:
                header_list = self.read_vcf_header(f)
        elif self.input_format in ["vcf"]:
            with open(input, 'rt') as f:
                header_list = self.read_vcf_header(f)
        # elif self.input_format in ["tsv", "csv", "psv"]:
        #     # check hdr file
        #     with open(input, 'rt') as f:
        #         header_list = self.read_vcf_header(f)
        elif self.input_format in ["parquet", "tsv", "csv", "psv"]:
            #print(self.config)
            #print("TODO: find header in another file in param, or if parquet.hdr exists")
            if config.get("header_file",None):
                #print("header in config")
                with open(config.get("header_file"), 'rt') as f:
                    header_list = self.read_vcf_header(f)
            elif os.path.exists((input+".hdr")):
                with open(input+".hdr", 'rt') as f:
                    header_list = self.read_vcf_header(f)
            else:
                log.error(f"No header for file {input}!!!")
                raise ValueError(f"No header for file {input}!!!")
                #return False
        elif self.input_format in ["db", "duckdb"]:
            print(self.config)
            print("TODO: find in db, table header, list of lines in header")
        else:
            with open(input, 'rt') as f:
                header_list = self.read_vcf_header(f)

        # header as list
        self.header_list = header_list

        # header as VCF object
        self.header_vcf = vcf.Reader(io.StringIO("\n".join(header_list)))

    def set_dataframe(self, dataframe):
        """
        The function takes in a dataframe and sets it as the dataframe attribute of the class
        
        :param dataframe: The dataframe that you want to plot
        """
        self.dataframe = dataframe

    # GET section

    def get_overview(self):
        """
        The function prints the input, output, config, and dataframe of the current object
        """
        table_variants_from = self.get_table_variants(clause="from")
        table_variants_select = self.get_table_variants(clause="select")
        sql_columns = self.get_header_columns_as_sql()
        sql_query_export = f"SELECT {sql_columns} FROM {table_variants_from}"
        df = self.conn.execute(sql_query_export).df()
        log.info("Input:  " + str(self.get_input()) + " [" + str(str(self.get_input_format())) + "]")
        log.info("Output: " + str(self.get_output()) + " [" + str(str(self.get_output_format())) + "]")
        #log.info("Config: " + str(self.get_config()))
        for d in str(json.dumps(self.get_config(), indent=4, sort_keys=True)).split("\n"):
            log.info(str(d))
        #log.info("Config: " + str(json.dump(self.get_config(),indent=4)))
        #log.info("Param:  " + str(self.get_param()))
        for d in str(json.dumps(self.get_param(), indent=4, sort_keys=True)).split("\n"):
            log.info(str(d))
        log.info("Sample list: " + str(self.get_header_sample_list()))
        #log.info("Dataframe\n"+str(df))
        log.info("Dataframe: ")
        for d in str(df).split("\n"):
            log.info(str(d))
        
        #print(str(json.dumps(self.get_config(), indent=4, sort_keys=True)))

    def get_stats(self):
        """
        The function prints statistics of the current object
        """
        table_variants_from = self.get_table_variants(clause="from")
        table_variants_select = self.get_table_variants(clause="select")
        sql_columns = self.get_header_columns_as_sql()

        # Samples
        nb_of_samples = len(self.get_header_sample_list())
        
    
        # Variants
        sql_query_nb_variant = f"SELECT count(*) as count FROM {table_variants_from}"
        nb_of_variants = self.conn.execute(sql_query_nb_variant).df()["count"][0]

        # Variants by chr
        sql_query_nb_variant = f"SELECT \"#CHROM\" as CHROM, count(*) as count, (count(*)*100/{nb_of_variants}) || '%' as percentage FROM {table_variants_from} GROUP BY \"#CHROM\""
        nb_of_variants_by_chrom = self.conn.execute(sql_query_nb_variant).df().sort_values(by=['CHROM'], kind='quicksort')

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
                        REGEXP_EXTRACT({sample}, '^([0-9/|\.]*)[:]*',1) as genotype,
                        count(REGEXP_EXTRACT({sample}, '^([0-9/|\.]*)[:]*',1)) as count,
                        (count(REGEXP_EXTRACT({sample}, '^([0-9/|\.]*)[:]*',1))*100/{nb_of_variants}) || '%' as percentage
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
            log.info(str(d))
        log.info("Ti/Ts/Tv Ratio:")
        for d in str(tstv).split("\n"):
            log.info(str(d))
        log.info("Number of Sample(s): " + str(nb_of_samples))
        log.info(f"Genotypes:")
        for sample in genotypes:
            for d in str(genotypes[sample]).split("\n"):
                log.info(str(d))




    def get_input(self):
        """
        It returns the value of the input variable.
        :return: The input is being returned.
        """
        return self.input

    def get_input_format(self):
        """
        It returns the format of the input variable.
        :return: The format is being returned.
        """
        return self.input_format

    def get_output(self):
        """
        It returns the output of the neuron.
        :return: The output of the neural network.
        """
        return self.output

    def get_output_format(self):
        """
        It returns the format of the input variable.
        :return: The format is being returned.
        """
        return self.output_format

    def get_config(self):
        """
        It returns the config
        :return: The config variable is being returned.
        """
        return self.config

    def get_param(self):
        """
        It returns the param
        :return: The param variable is being returned.
        """
        return self.param

    def get_connexion_db(self):
        """
        It returns the connexion_db attribute of the object
        :return: The connexion_db is being returned.
        """
        self.set_connexion_db()
        return self.connexion_db

    def get_prefix(self):
        """
        It returns the prefix of the object.
        :return: The prefix is being returned.
        """
        return self.prefix

    def get_table_variants(self, clause="select"):
        """
        This function returns the table_variants attribute of the object
        :param clause: the type of clause the table will be used. Either "select" or "from" (optional)
        :return: The table_variants attribute of the object.
        """
        access = self.get_config().get("access",None)
        #access = "RO"
        #print(access)
        #table_variants = f"{self.table_variants}"
        if clause in ["select"]:
            table_variants = self.table_variants
        elif clause in ["from"]:
            if self.get_input_format() in ["parquet"] and access in ["RO"]:
                input = self.get_input()
                table_variants = f"'{input}' as variants"
            else:
                table_variants = f"{self.table_variants} as variants"
        else:
            table_variants = self.table_variants
        #print(table_variants)
        return table_variants

    def get_tmp_dir(self):
        """
        It returns the value of the tmp_dir key in the config dictionary, or /tmp if the key doesn't
        exist
        :return: The value of the key "tmp_dir" in the config file.
        """
        return self.get_config().get("tmp_dir", "/tmp")

    def get_connexion_type(self):
        """
        > If the connexion type is not in the list of allowed connexion types, raise a ValueError
        :return: The connexion type is being returned.
        """
        connexion_type = self.get_config().get("connexion_type", "memory")
        # check connexion type
        if connexion_type not in ["memory", "tmpfile"]:
            log.error(f"Connexion type {connexion_type} failed!!! Either 'memory' or 'tmpfile'")
            raise ValueError(
                f"Connexion type {connexion_type} failed!!! Either 'memory' or 'tmpfile'")
        return connexion_type

    def get_connexion(self):
        """
        It returns the connection object
        :return: The connection object.
        """
        return self.conn

    def get_header(self, type="vcf"):
        """
        This function returns the header of the VCF file as a list of strings
        :param type: the type of header you want to get, defaults to vcf (optional)
        :return: The header of the vcf file.
        """
        if type == "vcf":
            return self.header_vcf
        elif type == "list":
            return self.header_list

    def get_infos_field_as_dict(self, info_field:str = "") -> dict:
        """
        It takes a string of the form `
        
        :param info_field: The INFO field to parse
        :type info_field: str
        :return: A dictionary with the info field as key and the value as value.
        """

        #string = '##INFO=<ID=non_topmed_AF_popmax,Number=.,Type=String,Description="non_topmed_AF_popmax">'

        # Split in parts
        parts = info_field.strip().split('<')[1].split('>')[0].split(',')

        # Create dictionnary
        info_dict = {part.split('=')[0]: part.split('=')[1].replace('"', '') for part in parts}

        # Check other variables?
        # if 'Source' in info_dict:
        #     info_dict['Source'] = info_dict['source'].replace('"', '')
        # if 'Version' in info_dict:
        #     info_dict['Version'] = info_dict['version'].replace('"', '')

        # Return dict
        return info_dict

    def get_header_infos(self):
        """
        This function retruns header length (without #CHROM line)
        :return: The length of the header list.
        """
        header_infos = {}
        for line in self.header_list:
            try:
                found = re.match('##INFO=<ID=([^,]*),.*', line)
                if found:
                    # explode line
                    #print(line)
                    line_dict = self.get_infos_field_as_dict(line)
                    # Add exploded line
                    header_infos[found.groups()[0]] = line_dict
            except AttributeError:
                pass
        
        return header_infos

    def get_header_length(self):
        """
        This function retruns header length (without #CHROM line)
        :return: The length of the header list.
        """
        return len(self.header_list) - 1

    def get_header_columns(self):
        """
        This function retruns header length (without #CHROM line)
        :return: The length of the header list.
        """
        return self.header_list[-1]

    def get_header_columns_as_sql(self):
        """
        This function retruns header length (without #CHROM line)
        :return: The length of the header list.
        """
        sql_column_list = []
        for col in self.get_header_columns().strip().split("\t"):
            sql_column_list.append(f"\"{col}\"")
        sql_column = ",".join(sql_column_list)
        return sql_column

    def get_header_sample_list(self):
        """
        This function retruns header length (without #CHROM line)
        :return: The length of the header list.
        """
        return self.header_vcf.samples

    def get_verbose(self):
        """
        It returns the value of the "verbose" key in the config dictionary, or False if the key doesn't
        exist
        :return: The value of the key "verbose" in the config dictionary.
        """
        return self.get_config().get("verbose", False)

    def get_dataframe(self):
        """
        This function returns the dataframe of the class
        :return: The dataframe is being returned.
        """
        return self.dataframe

    # PROCESS section

    def insert_file_to_table(self, file, columns, header_len=0, sep='\t', chunksize=100000):
        for chunk in pd.read_csv(file, skiprows=header_len, sep=sep, chunksize=chunksize):
            #print(chunk)
            sql_insert_into = f"INSERT INTO variants ({columns}) SELECT {columns} FROM chunk"
            self.conn.execute(sql_insert_into)
            

    def append_to_parquet_table(self, dataframe, filepath=None, writer=None):
        """Method writes/append dataframes in parquet format.

        This method is used to write pandas DataFrame as pyarrow Table in parquet format. If the methods is invoked
        with writer, it appends dataframe to the already written pyarrow table.

        :param dataframe: pd.DataFrame to be written in parquet format.
        :param filepath: target file location for parquet file.
        :param writer: ParquetWriter object to write pyarrow tables in parquet format.
        :return: ParquetWriter object. This can be passed in the subsequenct method calls to append DataFrame
            in the pyarrow Table
        """
        table = pa.Table.from_pandas(dataframe)
        if writer is None:
            writer = pq.ParquetWriter(filepath, table.schema)
        writer.write_table(table=table)
        return writer



    def load_data(self):
        """
        > The function loads the data from the input file into the database
        """

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
            structure["FORMAT"]="VARCHAR"
            for sample in self.get_header_sample_list():
                structure_complete[sample] = "VARCHAR"

        # Columns list for create and insert
        sql_create_table_columns = []
        sql_create_table_columns_list = []
        for column in structure_complete:
            column_type = structure_complete[column]
            sql_create_table_columns.append(f"\"{column}\" {column_type}")
            sql_create_table_columns_list.append(f"\"{column}\"")

        # Create database
        sql_create_table_columns_sql = ", ".join(sql_create_table_columns)
        sql_create_table_columns_list_sql = ", ".join(sql_create_table_columns_list)
        sql_create_table = f"CREATE TABLE IF NOT EXISTS {self.table_variants} ({sql_create_table_columns_sql})"
        self.conn.execute(sql_create_table)

        # create index
        # CREATE INDEX idx_table_parquet ON table_parquet ("#CHROM", "POS", "ALT", "REF");
        sql_create_table_index = f'CREATE INDEX idx_{self.table_variants} ON {self.table_variants} ("#CHROM", "POS", "ALT", "REF")'
        self.conn.execute(sql_create_table_index)

        # chunksize = 100000
        chunksize = 100000

        # Access
        access = self.get_config().get("access",None)

        # tmp parquet
        tmp_parquet = "/tmp/tmp.parquet"

        # Load data
        if self.input_format in ["vcf", "gz", "tsv", "csv", "psv"]:
            header_len = self.get_header_length()

            # delimiter
            delimiter = "\t"
            if self.input_format in ["vcf", "gz"]:
                delimiter = "\t"
            elif self.input_format in ["tsv"]:
                delimiter = "\t"
            elif self.input_format in ["csv"]:
                delimiter = ","
            elif self.input_format in ["psv"]:
                delimiter = "|"
                
            # Load VCF
            if self.input_format in ["gz"]:
                with bgzf.open(self.input, 'rt') as file:
                    self.insert_file_to_table(file, columns=sql_create_table_columns_list_sql, header_len=header_len, sep=delimiter, chunksize=chunksize)
            elif self.input_format in ["vcf"]:
                with open(self.input, 'rt') as file:
                    self.insert_file_to_table(file, columns=sql_create_table_columns_list_sql, header_len=header_len, sep=delimiter, chunksize=chunksize)
            else:
                with open(self.input, 'rt') as file:
                    self.insert_file_to_table(file, columns=sql_create_table_columns_list_sql, sep=delimiter, chunksize=chunksize)

            # By parquet
            # if False:
            #     writer = None
            #     if self.input_format in ["gz"]:
            #         with bgzf.open(self.input, 'rt') as file:
            #             for chunk in pd.read_csv(file, skiprows=header_len, sep='\t', chunksize=chunksize):
            #                 writer = self.append_to_parquet_table(chunk, tmp_parquet, writer)
            #     else:
            #         with open(self.input, 'rt') as file:
            #             for chunk in pd.read_csv(file, skiprows=header_len, sep='\t', chunksize=chunksize):
            #                 writer = self.append_to_parquet_table(chunk, tmp_parquet, writer)
            #     if writer:
            #         writer.close()
            #     if access in ["RO"]:
            #         #print("NO loading data")
            #         self.drop_variants_table()
            #         sql_view = f"CREATE VIEW {self.table_variants} AS SELECT * FROM '{tmp_parquet}'"
            #         self.conn.execute(sql_view)
            #     else:
            #         sql_insert_table = f"COPY {self.table_variants} FROM '{tmp_parquet}'"
            #         self.conn.execute(sql_insert_table)

            #     # DEVEL
            #     df = pd.read_parquet(tmp_parquet)
            #     print(df)

            pass
        elif self.input_format in ["parquet"]:
            # Load Parquet
            
            #access = "RO"
            #print(access)
            if access in ["RO"]:
                #print("NO loading data")
                self.drop_variants_table()
                sql_view = f"CREATE VIEW {self.table_variants} AS SELECT * FROM '{self.input}'"
                self.conn.execute(sql_view)
            else:
                sql_insert_table = f"COPY {self.table_variants} FROM '{self.input}'"
                self.conn.execute(sql_insert_table)
            pass
        elif self.input_format in ["db", "duckdb"]:
            # Load DuckDB
            if access in ["RO"]:
                self.drop_variants_table()
                sql_view = f"CREATE VIEW {self.table_variants} AS SELECT * FROM input_conn.{self.table_variants}"
                self.conn.execute(sql_view)
            else:
                input_conn = duckdb.connect(self.input)
                sql_insert_table = f"INSERT INTO {self.table_variants} SELECT * FROM input_conn.{self.table_variants}"
                self.conn.execute(sql_insert_table)
            pass
        else:
            log.error(f"Input file format '{self.input_format}' not available")
            raise ValueError("Input file format '{self.input_format}' not available")

    def load_data_naive(self):
        """
        > The function loads the data from the input file into the database
        """
        if self.input_format in ["vcf", "gz"]:
            header_len = self.get_header_length()
            # Load VCF
            if self.input_format in ["gz"]:
                with bgzf.open(self.input, 'rt') as file:
                    dataframe = pd.read_csv(
                        file, skiprows=header_len, sep='\t')
            else:
                dataframe = pd.read_csv(
                    self.input, skiprows=header_len, sep='\t')
            self.set_dataframe(dataframe)
            # Variants reading
            self.conn.execute(
                f"CREATE TABLE {self.table_variants} AS SELECT * FROM dataframe")
            pass
        elif self.input_format in ["parquet"]:
            # Load Parquet
            pass
        elif self.input_format in ["db", "duckdb"]:
            # Load DuckDB
            pass
        else:
            log.error(f"Input file format '{self.input_format}' not available")
            raise ValueError(f"Input file format '{self.input_format}' not available")

    def read_vcf_header(self, f):
        """
        It reads the header of a VCF file and returns a list of the header lines
        
        :param f: the file object
        :return: The header lines of the VCF file.
        """
        header_list = []
        for line in f:
            #print(line)
            header_list.append(line)
            if line.startswith('#CHROM'):
                break
        return header_list

    def execute_query(self, query):
        """
        It takes a query as an argument, executes it, and returns the results
        
        :param query: The query to be executed
        :return: The result of the query is being returned.
        """
        if query:
            return self.conn.execute(query) #.fetchall()
        else:
            return None

    def export_output(self):
        #print("Export output")
        if self.get_output():
            #print(f"Export output {self.get_output()} now...")

            output_file = self.get_output()
            sql_column = self.get_header_columns_as_sql()
            table_variants = self.get_table_variants()
            sql_query_hard = ""
            sql_query_sort = ""
            sql_query_limit = ""

            if self.get_output_format() in ["parquet"]:
                
                # delimiter
                delimiter = "\t"
                if self.get_output_format() in ["csv"]:
                    delimiter = ","
                if self.get_output_format() in ["psv"]:
                    delimiter = "|"

                # Header
                tmp_header_name = output_file + ".hdr"
                f = open(tmp_header_name, 'w')
                vcf_writer = vcf.Writer(f, self.header_vcf)
                f.close()

                # Export parquet
                sql_query_export = f"COPY (SELECT {sql_column} FROM {table_variants} WHERE 1 {sql_query_hard} {sql_query_sort} {sql_query_limit}) TO '{output_file}' WITH (FORMAT PARQUET)"
                self.conn.execute(sql_query_export)

            elif self.get_output_format() in ["db", "duckdb"]:

                # Export duckDB
                sql_query_export = f"EXPORT DATABASE '{output_file}'"
                self.conn.execute(sql_query_export)

            elif self.get_output_format() in ["tsv", "csv", "psv"]:

                # delimiter
                delimiter = "\t"
                if self.get_output_format() in ["csv"]:
                    delimiter = ","
                if self.get_output_format() in ["psv"]:
                    delimiter = "|"

                # Header
                tmp_header_name = output_file + ".hdr"
                f = open(tmp_header_name, 'w')
                vcf_writer = vcf.Writer(f, self.header_vcf)
                f.close()

                # Export TSV/CSV
                #sql_query_export = f"EXPORT DATABASE '{output_file}'  (FORMAT CSV, DELIMITER '{delimiter}')"
                sql_query_export = f"COPY (SELECT {sql_column} FROM {table_variants} WHERE 1 {sql_query_hard} {sql_query_sort} {sql_query_limit}) TO '{output_file}' WITH (FORMAT CSV, DELIMITER '{delimiter}', HEADER)"
                self.conn.execute(sql_query_export)

            elif self.get_output_format() in ["vcf", "gz"]:
                # Extract VCF
                #print("#[INFO] VCF Output - Extract VCF...")

                # Header
                tmp_header = NamedTemporaryFile(prefix=self.get_prefix(), dir=self.get_tmp_dir())
                tmp_header_name = tmp_header.name
                f = open(tmp_header_name, 'w')
                vcf_writer = vcf.Writer(f, self.header_vcf)
                f.close()

                # Variants
                tmp_variants = NamedTemporaryFile(prefix=self.get_prefix(), dir=self.get_tmp_dir())
                tmp_variants_name = tmp_variants.name
                sql_query_export = f"COPY (SELECT {sql_column} FROM {table_variants} WHERE 1 {sql_query_hard} {sql_query_sort} {sql_query_limit}) TO '{tmp_variants_name}' WITH (FORMAT CSV, DELIMITER '\t', HEADER, QUOTE '', COMPRESSION 'gzip')"
                #print(sql_query_export)
                self.conn.execute(sql_query_export)

                # Create output
                # Cat header and variants
                command_gzip = ""
                if self.output_format in ["gz"]:
                    command_gzip = f" | bgzip -c "
                command = f"grep '^#CHROM' -v {tmp_header_name} {command_gzip} > {output_file}; bgzip -dc {tmp_variants_name} {command_gzip} >> {output_file}"
                #print(command)
                subprocess.run(command, shell=True)

    def export_variant_vcf(self, vcf_file, type="vcf", remove_info=False):
        
        # Extract VCF
        #print("#[INFO] Export VCF...")

        sql_column = self.get_header_columns_as_sql()
        table_variants = self.get_table_variants()
        sql_query_hard = ""
        sql_query_sort = ""
        sql_query_limit = ""
        if remove_info:
            info_field = "'.' as INFO"
        else:
            info_field = "INFO"

        # Header
        tmp_header = NamedTemporaryFile(prefix=self.get_prefix(), dir=self.get_tmp_dir())
        tmp_header_name = tmp_header.name
        f = open(tmp_header_name, 'w')
        vcf_writer = vcf.Writer(f, self.header_vcf)
        f.close()

        # Variants
        tmp_variants = NamedTemporaryFile(prefix=self.get_prefix(), dir=self.get_tmp_dir())
        tmp_variants_name = tmp_variants.name
        sql_query_export = f"COPY (SELECT \"#CHROM\", POS, ID, REF, ALT, QUAL, FILTER, {info_field} FROM {table_variants} WHERE 1 {sql_query_hard} {sql_query_sort} {sql_query_limit}) TO '{tmp_variants_name}' WITH (FORMAT CSV, DELIMITER '\t', HEADER, QUOTE '', COMPRESSION 'gzip')"
        #print(sql_query_export)
        self.conn.execute(sql_query_export)

        # Create output
        # Cat header and variants
        command_gzip = ""
        if type in ["gz"]:
            command_gzip = f" | bgzip -c "
        command = f"grep '^#CHROM' -v {tmp_header_name} {command_gzip} > {vcf_file}; bgzip -dc {tmp_variants_name} {command_gzip} >> {vcf_file}"
        #print(command)
        subprocess.run(command, shell=True)


    def run_commands(self,commands=[],threads=1):
        # commands = [
        #     "sleep 1",
        #     "sleep 3",
        #     "sleep 1",
        #     "sleep 2"
        # ]
        #print("run parallel commands")
        run_parallel_commands(commands, threads)


    def set_id_null(self, threads=1):
        functions = [function_query(self,"UPDATE variants SET ID='.' WHERE REF='A'"), function_query(self,"UPDATE variants SET ID='.' WHERE REF='C'"), function_query(self,"UPDATE variants SET ID='.' WHERE REF='T'"), function_query(self,"UPDATE variants SET ID='.' WHERE REF='G'")]
        run_parallel_functions(functions, threads)


    def get_threads(self):
        return int(self.config.get("threads",1))

    def annotation_bcftools(self, threads=None):

        log.debug("Start annotation with bcftools and with vcf and bed databases")
        log.debug(self.param.get("annotation",{}).get("bcftools"))

        return

        if not threads:
            threads = self.get_threads()

        delete_tmp=False

        input_vcf = self.get_input()

        # Export in VCF
        tmp_vcf = NamedTemporaryFile(prefix=self.get_prefix(), dir=self.get_tmp_dir(), suffix=".vcf.gz", delete=delete_tmp)
        tmp_vcf_name = tmp_vcf.name
        self.export_variant_vcf(vcf_file=tmp_vcf_name,type="gz")
        run_parallel_commands([f"tabix {tmp_vcf_name}"], 1)

        annotate_dir = "/mnt/STARK/databases/annotations/20230214/hg19/"
        # annotations_list = {
        #     "clinvar_20210123": "INFO",
        #     "cosmic70": "INFO",
        #     "dbscsnv11": "INFO",
        #     "esp6500siv2_aa": "INFO",
        #     "esp6500siv2_all": "INFO",
        #     "esp6500siv2_ea": "INFO",
        #     "exac03": "INFO",
        #     "gene4denovo201907": "INFO",
        #     "gme": "INFO",
        #     "gnomad211_exome": "INFO",
        #     "icgc28": "INFO",
        #     "intervar_20180118": "INFO",
        #     "mcap13": "INFO",
        #     "mitimpact24": "INFO",
        #     "nci60": "INFO",
        #     "popfreq_max_20150413": "INFO",
        #     #"regsnpintron": "INFO",
        # }
        # annotations_list = {
        #     "clinvar_20210123": "INFO",
        #     "cosmic70": "INFO",
        #     "dbscsnv11": "INFO",
        #     "esp6500siv2_aa": "INFO",
        #     "esp6500siv2_all": "INFO",
        #     "esp6500siv2_ea": "INFO",
        #     "exac03": "INFO",
        #     "gene4denovo201907": "INFO",
        #     "gme": "INFO",
        #     "gnomad211_exome": "INFO",
        #     "icgc28": "INFO",
        #     "intervar_20180118": "INFO",
        #     "mcap13": "INFO",
        #     "mitimpact24": "INFO",
        #     "nci60": "INFO",
        # }
        annotations_list = {
            "cosmic70": "INFO",
            "nci60": "INFO",
            "revel": "INFO",
            "clinvar_20210123": "INFO",
            "avsnp150": "INFO",
            }
        # annotations_list = {
        #     "cosmic70": "INFO",
        #     "nci60": "INFO",
        #     }
        # annotations_list = {
        #     "cosmic70": "INFO",
        #     }

        commands = []
        tmp_ann_vcf_list = []
        for ann in annotations_list:
            ann_infos = annotations_list[ann]
            tmp_ann_vcf = NamedTemporaryFile(prefix=self.get_prefix(), dir=self.get_tmp_dir(), suffix=".vcf.gz", delete=delete_tmp)
            tmp_ann_vcf_name = tmp_ann_vcf.name
            tmp_ann_vcf_list.append(f"{tmp_ann_vcf_name}")
            commands.append(f"bcftools annotate -a {annotate_dir}/{ann}.vcf.gz -h {annotate_dir}/{ann}.vcf.hdr -c {ann_infos} {tmp_vcf_name} 2>/dev/null | bgzip -c > {tmp_ann_vcf_name} && tabix  {tmp_ann_vcf_name}")
        run_parallel_commands(commands, threads)

        # check 

        # Merge
        tmp_annotate_vcf = NamedTemporaryFile(prefix=self.get_prefix(), dir=self.get_tmp_dir(), suffix=".vcf.gz", delete=delete_tmp)
        tmp_annotate_vcf_name = tmp_annotate_vcf.name
        tmp_ann_vcf_list_cmd = " ".join(tmp_ann_vcf_list)
        merge_command = f"bcftools merge --no-index --threads={threads} {input_vcf} {tmp_ann_vcf_list_cmd} 2>/dev/null | bgzip --threads={threads} -c > {tmp_annotate_vcf_name}"
        run_parallel_commands([merge_command], 1)

        self.update_from_vcf(tmp_annotate_vcf_name)
        #self.update_from_vcf_brutal(tmp_annotate_vcf_name)

    def annotation_parquet(self, threads=None):

        # DEBUG
        log.debug("Start annotation with parquet databases")
        
        # Threads
        if not threads:
            threads = self.get_threads()
        log.debug("Threads: "+str(threads))

        # DEBUG
        delete_tmp=True
        if self.get_config().get("verbosity","warning") in ["debug"]:
            delete_tmp=False
            log.debug("Delete tmp files/folders: "+str(delete_tmp))

        # Config
        databases_folders = self.config.get("folders",{}).get("databases",{}).get("parquet",["."])
        log.debug("Databases annotations: " + str(databases_folders))

        # Param
        annotations = self.param.get("annotation",{}).get("parquet",{}).get("annotations",None)
        log.debug("Annotations: " + str(annotations))
        
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
        tmp_vcf = NamedTemporaryFile(prefix=self.get_prefix(), dir=self.get_tmp_dir(), suffix=".vcf.gz", delete=delete_tmp)
        tmp_vcf_name = tmp_vcf.name
        self.export_variant_vcf(vcf_file=tmp_vcf_name,type="gz")
        run_parallel_commands([f"tabix {tmp_vcf_name}"], 1)
        log.debug("Create initial file to annotate: " +str(tmp_vcf_name))

        # VCF header
        vcf_reader = self.get_header()
        log.debug("Initial header: " +str(vcf_reader.infos))

        # Existing annotations
        for vcf_annotation in self.get_header_infos():
            vcf_annotation_line = self.get_header_infos().get(vcf_annotation)
            log.debug(f"Existing annotations in VCF: {vcf_annotation} [{vcf_annotation_line}]")

        if annotations:
            for annotation in annotations:
                annotation_fields = annotations[annotation]
                log.debug(f"Annotation with database '{annotation}'")
                log.debug(f"Annotation with database '{annotation}' - fields: {annotation_fields}")

                # Find parquet file and header file
                parquet_file = None
                parquet_hdr_file = None
                for databases_folder in databases_folders:
                    parquet_file = None
                    parquet_hdr_file = None
                    log.debug("Annotation with database file check: " +str(databases_folder+"/"+annotation+".parquet"))
                    # Parquet .parquet
                    if os.path.exists(annotation):
                        parquet_file = annotation
                    elif os.path.exists(databases_folder+"/"+annotation+".parquet"):
                        parquet_file = databases_folder+"/"+annotation+".parquet"
                    # Header .hdr
                    if os.path.exists(parquet_file+".hdr"):
                        parquet_hdr_file = parquet_file+".hdr"
                    elif os.path.exists(databases_folder+"/"+annotation+".parquet.hdr"):
                        parquet_hdr_file = databases_folder+"/"+annotation+".parquet.hdr"
                    if parquet_file and parquet_hdr_file:
                        break

                if not parquet_file or not parquet_hdr_file:
                    log.error("Annotation with database failed: file not found")
                    raise ValueError("Annotation with database failed: file not found")
                else:
                
                    log.debug(f"Annotation with database '{annotation}' - file: " + str(parquet_file) + " and " + str(parquet_hdr_file) )
            
                    # Load header as VCF object
                    parquet_hdr_vcf = VCFDataObject(input=parquet_hdr_file)
                    parquet_hdr_vcf_header_infos = parquet_hdr_vcf.get_header_infos()
                    log.debug("Annotation database header: " +str(parquet_hdr_vcf_header_infos))
                    # Alternative ith PyVCF3
                    # parquet_hdr_vcf_2 = vcf.Reader(io.StringIO(parquet_hdr_file))
                    # log.debug("Annotation database header 2: " +str(parquet_hdr_vcf_2.infos))

                    # For all fields in database
                    if "ALL" in annotation_fields or "INFO" in annotation_fields:
                        annotation_fields = {key: key for key in parquet_hdr_vcf_header_infos}
                        log.debug("Annotation database header - All annotations added: " +str(annotation_fields))

                    # List of annotation fields to use
                    sql_query_annotations_list = []

                    # Number of fields
                    nb_annotation_field = 0

                    for annotation_field in annotation_fields:

                        # field new name, if parametered
                        annotation_fields_new_name = annotation_fields.get(annotation_field,annotation_field)
                        if not annotation_fields_new_name:
                            annotation_fields_new_name = annotation_field

                        if annotation_field in parquet_hdr_vcf.get_header_infos() and annotation_field not in self.get_header_infos():

                            # Sep between fields in INFO
                            nb_annotation_field += 1
                            if nb_annotation_field > 1:
                                annotation_field_sep = ";"
                            else:
                                annotation_field_sep = ""

                            log.info(f"Annotation with database '{annotation}' - '{annotation_field}' -> 'INFO/{annotation_fields_new_name}'")
                            
                            # Annotation query fields
                            sql_query_annotations_list.append(f"|| '{annotation_field_sep}' ||'{annotation_fields_new_name}=' || REGEXP_EXTRACT(';' || table_parquet.INFO, ';{annotation_field}=([^;]*)',1) ")

                            
                            # Add INFO field to header
                            parquet_hdr_vcf_header_infos_number = parquet_hdr_vcf_header_infos[annotation_field].get("Number", ".")
                            parquet_hdr_vcf_header_infos_type = parquet_hdr_vcf_header_infos[annotation_field].get("Type", "String")
                            parquet_hdr_vcf_header_infos_description = parquet_hdr_vcf_header_infos[annotation_field].get("Description", f"{annotation_field} description")
                            parquet_hdr_vcf_header_infos_source = parquet_hdr_vcf_header_infos[annotation_field].get("Source", "unknown")
                            parquet_hdr_vcf_header_infos_version = parquet_hdr_vcf_header_infos[annotation_field].get("Version", "unknown")
                            #print(parquet_hdr_vcf_header_infos_ID)
                            vcf_reader.infos[annotation_fields_new_name] = vcf.parser._Info(
                                        annotation_fields_new_name,
                                        parquet_hdr_vcf_header_infos_number,
                                        parquet_hdr_vcf_header_infos_type,
                                        parquet_hdr_vcf_header_infos_description,
                                        parquet_hdr_vcf_header_infos_source,
                                        parquet_hdr_vcf_header_infos_version,
                                        self.code_type_map[parquet_hdr_vcf_header_infos_type]
                                    )

                        else:
                            if annotation_field in parquet_hdr_vcf.get_header_infos():
                                log.warning(f"Annotation with database '{annotation}' - '{annotation_field}' [{nb_annotation_field}] - already exists (skipped)")
                            elif annotation_field not in self.get_header_infos():
                                log.warning(f"Annotation with database '{annotation}' - '{annotation_field}' [{nb_annotation_field}] - not available in parquet file")
                    

                    # Annotate
                    log.info(f"Annotation with database '{annotation}' - Annotation...")

                    sql_query_annotations_list_sql = " ".join(sql_query_annotations_list)

                    
                    log.debug("Find chromosomes ")
                    sql_query_chromosomes = f"""SELECT table_variants.\"#CHROM\" as CHROM FROM {table_variants} as table_variants GROUP BY table_variants.\"#CHROM\""""
                    results = self.conn.execute(f"{sql_query_chromosomes}").df()["CHROM"]
                    log.debug("Chromosomes found:\n" + str(list(results)))

                    # Check chromosomes list (and variant max position)
                    sql_query_chromosomes_max_pos = f"""SELECT table_variants.\"#CHROM\" as CHROM, MAX(table_variants.\"POS\") as MAX_POS FROM {table_variants} as table_variants GROUP BY table_variants.\"#CHROM\""""
                    sql_query_chromosomes_max_pos_df = self.conn.execute(sql_query_chromosomes_max_pos).df()
                    # Create dictionnary with chromosomes (and max position)
                    sql_query_chromosomes_max_pos_dictionary = sql_query_chromosomes_max_pos_df.groupby('CHROM').apply(lambda x: {'max_pos': x['MAX_POS'].max()}).to_dict()

                    # Affichage du dictionnaire
                    log.debug("Chromosomes max pos found:\n" + str(sql_query_chromosomes_max_pos_dictionary))

                    # Batch parameters
                    param_batch_annotation_databases_window = self.get_param().get("annotation",{}).get("parquet",{}).get("batch",{}).get("window",0)
                    param_batch_annotation_databases_auto = self.get_param().get("annotation",{}).get("parquet",{}).get("batch",{}).get("auto","each_chrom")
                    param_batch_annotation_databases_batch = self.get_param().get("annotation",{}).get("parquet",{}).get("batch",{}).get("batch",1000)

                    # Init
                    batch_annotation_databases_window = param_batch_annotation_databases_window


                    # nb_of_variant_annotated
                    nb_of_variant_annotated = 0

                    for chrom in sql_query_chromosomes_max_pos_dictionary:

                        # nb_of_variant_annotated_by_chrom
                        nb_of_variant_annotated_by_chrom = 0

                        # Autodetect range of bases to split/chunk
                        if not param_batch_annotation_databases_window and (not batch_annotation_databases_window or param_batch_annotation_databases_auto == "each_chrom"):
                            # Query to detect window of "batch" number of variant in the chromosome
                            autodetect_range = f"""SELECT table_parquet.\"POS\" as POS FROM {table_variants} as table_variants
                                INNER JOIN '{parquet_file}' as table_parquet ON
                                table_parquet.\"#CHROM\" = '{chrom}' AND table_variants.\"#CHROM\" = '{chrom}'
                                AND table_parquet.\"#CHROM\" = table_variants.\"#CHROM\"
                                AND table_parquet.\"POS\" = table_variants.\"POS\"
                                AND table_parquet.\"ALT\" = table_variants.\"ALT\"
                                AND table_parquet.\"REF\" = table_variants.\"REF\"
                                LIMIT {param_batch_annotation_databases_batch}
                                """
                            autodetect_range_results = self.conn.execute(f"{autodetect_range}").df()["POS"]

                            # Window is max position, if "batch" variants were found, otherwise Maximum Effort!!!
                            if len(autodetect_range_results) == param_batch_annotation_databases_batch:
                                batch_annotation_databases_window = autodetect_range_results[len(autodetect_range_results)-1]
                            else:
                                batch_annotation_databases_window = 1000000000000000 # maximum effort

                        # Get position of the farthest variant (max position) in the chromosome
                        sql_query_chromosomes_max_pos_dictionary_max_pos = sql_query_chromosomes_max_pos_dictionary.get(chrom,{}).get("max_pos")
                        # Create intervals from 0 to max position variant, with the batch window previously defined
                        sql_query_intervals = split_interval(0, sql_query_chromosomes_max_pos_dictionary_max_pos, step=batch_annotation_databases_window, ncuts=None)
                        #print("sql_query_intervals: " + str(sql_query_intervals))

                        # Interval Start/Stop
                        sql_query_interval_start = sql_query_intervals[0]
                        
                        # For each interval
                        for i in sql_query_intervals[1:]:

                            # Interval Start/Stop
                            sql_query_interval_stop = i
                            
                            log.debug(f"Annotation with database '{annotation}' - Chromosome '{chrom}' - Interval [{sql_query_interval_start}-{sql_query_interval_stop}] ...")
                            
                            # Create query to update
                            sql_query_annotation_chrom_interval_pos = f"""
                                UPDATE {table_variants} as table_variants
                                    SET INFO = CASE WHEN table_variants.INFO NOT IN ('','.') THEN table_variants.INFO ELSE '' END || CASE WHEN table_variants.INFO NOT IN ('','.') THEN ';' ELSE '' END {sql_query_annotations_list_sql}
                                    FROM '{parquet_file}' as table_parquet
                                    WHERE table_parquet.\"#CHROM\" = table_variants.\"#CHROM\"
                                        AND table_parquet.\"POS\" = table_variants.\"POS\"
                                        AND table_parquet.\"ALT\" = table_variants.\"ALT\"
                                        AND table_parquet.\"REF\" = table_variants.\"REF\"
                                        AND table_parquet.\"#CHROM\" = '{chrom}'
                                        AND table_parquet.\"POS\" > {sql_query_interval_start}
                                        AND table_parquet.\"POS\" <= {sql_query_interval_stop}
                                        """
                            log.debug("Create SQL query: " + str(sql_query_annotation_chrom_interval_pos))
                            
                            #Process query
                            result = self.conn.execute(sql_query_annotation_chrom_interval_pos)
                            nb_of_variant_annotated_by_chrom_by_interval = result.df()["Count"][0]
                            log.debug(f"Annotation with database '{annotation}' - Chromosome '{chrom}' - Interval [{sql_query_interval_start}-{sql_query_interval_stop}] - {nb_of_variant_annotated_by_chrom_by_interval} variants annotated")

                            # nb_of_variant_annotated_by_chrom
                            nb_of_variant_annotated_by_chrom += nb_of_variant_annotated_by_chrom_by_interval
                            

                            # Interval Start/Stop
                            sql_query_interval_start = sql_query_interval_stop

                            # Number of variant annotated
                            # if False:
                            #     log.debug("Create SQL query: " + str(sql_query_annotation_chrom_interval_pos))
                            #     result = self.conn.execute(sql_query_annotation_chrom_interval_pos) 
                            #     nb_of_variant_annotated = result.df()["Count"][0]
                            #     log.info(f"Annotation with database '{annotation}' - {nb_of_variant_annotated} variants annotated")
                    
                        # nb_of_variant_annotated
                        log.debug(f"Annotation with database '{annotation}' - Chromosome '{chrom}' - {nb_of_variant_annotated_by_chrom} variants annotated")
                        nb_of_variant_annotated += nb_of_variant_annotated_by_chrom


                    log.info(f"Annotation with database '{annotation}' - Annotation of {nb_of_variant_annotated} variants")

                    log.debug("Final header: " +str(vcf_reader.infos))

                    # DEBUG
                    if self.get_config().get("verbosity","warning") in ["debug"]:
                        #results = self.conn.execute(f"SELECT * FROM {table_variants} WHERE INFO LIKE '%;CLN%'")
                        results = self.conn.execute(f"SELECT * FROM {table_variants} ")
                        log.info("Annotation with database results:\n" + str(results.df()))
                        #print(results.fetchall())


        # # DEBUG
        # results = self.conn.execute(f"SELECT * FROM {table_variants} WHERE INFO LIKE '%nci60%'").df()
        # log.debug("Annotation with database results:\n" + str(results))
        # #print(results.fetchall())


        return

            

            



    def update_from_vcf(self,vcf_file):
        #print("vcf_file: "+vcf_file)
        table_variants = self.get_table_variants()
        tmp_merged_vcf = NamedTemporaryFile(prefix=self.get_prefix(), dir=self.get_tmp_dir(), suffix=".parquet", delete=True)
        tmp_merged_vcf_name = tmp_merged_vcf.name
        mergeVCF = VCFDataObject(None, vcf_file, tmp_merged_vcf_name, config=self.get_config())
        #print("PARQUET")
        #print(tmp_merged_vcf_name)
        mergeVCF.load_data()
        #mergeVCF.get_overview()
        mergeVCF.export_output()
        sql_query_update = f"""
        UPDATE {table_variants} as table_variants
            SET INFO = (SELECT table_parquet.INFO FROM '{tmp_merged_vcf_name}' as table_parquet
                                WHERE table_parquet.\"#CHROM\" = table_variants.\"#CHROM\"
                                AND table_parquet.\"POS\" = table_variants.\"POS\"
                                AND table_parquet.\"ALT\" = table_variants.\"ALT\"
                                AND table_parquet.\"REF\" = table_variants.\"REF\"
                        )
            ;
            """
        #print(sql_query_update)
        self.conn.execute(sql_query_update)
        #print(results.fetchall())
        
        # sql_query_update = f"""SELECT INFO FROM {table_variants};"""
        # print(sql_query_update)
        # results = self.conn.execute(sql_query_update)
        # print(results.fetchall())



    def update_from_vcf_brutal(self,vcf_file):
        table_variants = self.get_table_variants()
        self.drop_variants_table()
        self.set_input(vcf_file)
        self.set_header()
        self.load_data()


    def drop_variants_table(self):
        table_variants = self.get_table_variants()
        sql_table_variants = f"DROP TABLE {table_variants}"
        self.conn.execute(sql_table_variants)



# Main function
def main():
    """
    It loads a VCF file in multiple format (VCF, parquet, DB), and process, query, export data
    """
    
    # Args
    parser = argparse.ArgumentParser(
        description="Load a VCF file in multiple format (VCF, parquet, DB), and process, query, export data")
    parser.add_argument(
        "--input", help="Input file path (format: vcf, vcf.gz, parquet or db) Required", required=True)
    parser.add_argument(
        "--output", help="Output file path (format: vcf, vcf.gz, parquet or db) Required", required=False)
    parser.add_argument(
        "--config", help="Configuration file (format: JSON) (default: {})", default="{}")
    parser.add_argument(
        "--param", help="Parameters file (format: JSON) (default: {})", default="{}")
    parser.add_argument(
        "--query", help="Query (format: SQL) (default: 'SELECT * FROM variants LIMIT 5')", default=None)
    parser.add_argument("--overview", "--overview_header", help="Overview after loading data", action="store_true")
    parser.add_argument("--overview_footer", help="Overview before data processing", action="store_true")
    parser.add_argument("--stats", "--stats_header", help="Statistics after loading data", action="store_true")
    parser.add_argument("--stats_footer", help="Statistics before data processing", action="store_true")
    parser.add_argument("--verbose", help="Verbose", action="store_true")
    parser.add_argument("--debug", help="Debug", action="store_true")
    parser.add_argument("--verbosity", help="Verbosity level: CRITICAL, ERROR, WARNING, INFO, DEBUG or NOTSET", required=False, default="warning")
    args = parser.parse_args()

    # Verbose
    if args.verbose:
        args.verbosity = "info"
    
    # Debug
    if args.debug:
        args.verbosity = "debug"
    
    # Logging
    set_log_level(args.verbosity)


    # Load configuration in JSON format
    if os.path.exists(args.config):
        with open(args.config) as config_file:
            config = json.load(config_file)
    else:
        config = json.loads(args.config)

    config["verbosity"] = args.verbosity

    # Load parameters in JSON format
    if os.path.exists(args.param):
        with open(args.param) as param_file:
            param = json.load(param_file)
    else:
        param = json.loads(args.param)

    
    # Create VCF object
    vcfdata_obj = VCFDataObject(None, args.input, args.output, config, param)

    # Connexion
    #connexion_db = ":memory:"
    
    if vcfdata_obj.get_input_format() in ["db", "duckdb"]:
        connexion_db = vcfdata_obj.get_input()

    if vcfdata_obj.get_connexion_type() in ["memory", None]:
        connexion_db = ":memory:"
    elif vcfdata_obj.get_connexion_type() in ["tmpfile"]:
        tmp_name = tempfile.mkdtemp(prefix=vcfdata_obj.get_prefix(), dir=vcfdata_obj.get_tmp_dir(), suffix=".db")
        connexion_db = f"{tmp_name}/tmp.db"
    elif vcfdata_obj.get_connexion_type() != "":
        connexion_db = vcfdata_obj.get_connexion_type()
    else:
        connexion_db = ":memory:"
    
    #print("connexion_db: "+str(connexion_db))
    connexion_config = {}
    if config.get("threads",None):
        connexion_config["threads"] = config.get("threads")
    if config.get("memory_limit",None):
        connexion_config["memory_limit"] = config.get("memory_limit")
        

    conn = duckdb.connect(connexion_db, config=connexion_config)

    vcfdata_obj.set_connexion(conn)


    # Load data from input file
    log.info("Loading data...")
    vcfdata_obj.load_data()

    # Overview
    if args.overview:
        vcfdata_obj.get_overview()

    # Stats
    if args.stats:
        vcfdata_obj.get_stats()


    # Query
    if args.query or param.get("query",None):
        log.info("Querying...")
        if args.query:
            result = vcfdata_obj.execute_query(args.query)
        elif param.get("query",None):
            result = vcfdata_obj.execute_query(param.get("query",None))
        print(result.df())

    # Annotate
    if param.get("annotation",None):
        log.info("Annotations")
        if param.get("annotation",{}).get("parquet",None):
            log.info("Annotations 'parquet'...")
            vcfdata_obj.annotation_parquet()
        if param.get("annotation",{}).get("bcftools",None) and False:
            log.info("Annotations 'bcftools'...")
            vcfdata_obj.annotation_bcftools()
        if param.get("annotation",{}).get("annovar",None):
            log.info("Annotations 'annovar'...")
        if param.get("annotation",{}).get("snpeff",None):
            log.info("Annotations 'snpeff'...")
        if param.get("annotation",{}).get("varank",None):
            log.info("Annotations 'varank'...")
        # for annotation_method in param.get("annotation",None):
        #     print(annotation_method)
        #vcfdata_obj.annotate()
        #vcfdata_obj.annotate_parquet()
        # # Overview
        # if args.overview:
        #     vcfdata_obj.get_overview()
        # # Stats
        # if args.stats:
        #     vcfdata_obj.get_stats()


    # Output
    if args.output:
        log.info("Exporting...")
        vcfdata_obj.export_output()


    # Overview footer
    if args.overview_footer:
        vcfdata_obj.get_overview()

    # Stats footer
    if args.stats_footer:
        vcfdata_obj.get_overview()


if __name__ == "__main__":
    main()
