import gc
import glob
import io
import math
import os
from pathlib import Path
import random
import re
import sqlite3
import string
from tempfile import TemporaryDirectory
import tempfile
import duckdb  # type: ignore
import json
import yaml  # type: ignore
import Bio.bgzf as bgzf  # type: ignore
import pandas as pd  # type: ignore
import vcf  # type: ignore
import logging as log
import fastparquet as fp  # type: ignore

from howard.functions.commons import (
    DEFAULT_CHUNK_SIZE,
    DEFAULT_ASSEMBLY,
    CODE_TYPE_MAP,
    cast_columns_query,
    clean_annotation_field,
    convert_markdown_to_html,
    convert_markdown_to_pdf,
    detect_column_type,
    duckdb_has_spilled,
    escape_markdown_table_chars,
    extract_memory_in_go,
    full_path,
    get_file_compressed,
    get_file_format,
    get_memory,
    get_random,
    remove_if_exists,
    run_parallel_commands,
    vcf_required,
    file_format_delimiters,
    code_type_map_to_sql,
    sort_contigs,
    choose_update_strategy_safe,
)

from howard.objects.database import Database


# Mixins for Variants class
from .variants_mixin.view import variants_view
from .variants_mixin.transcripts import variants_transcripts
from .variants_mixin.annotation import variants_annotation
from .variants_mixin.calculation import variants_calculation
from .variants_mixin.prioritization import variants_prioritization
from .variants_mixin.tmp import variants_tmp


class Variants(
    variants_view,
    variants_transcripts,
    variants_annotation,
    variants_calculation,
    variants_prioritization,
    variants_tmp,
):

    def __init__(
        self,
        conn=None,
        input: str = None,
        output: str = None,
        config: dict = {},
        param: dict = {},
        load: bool = False,
    ) -> None:
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

        # Samples
        self.set_samples()

        # Temporary files
        self.set_tmp_files()

        # Load data
        if load:
            self.load_data(input)


    def get_access(self, default: str = None) -> str:
        """

        The function `get_access` retrieves the access level from the parameters or configuration,
        returning a default value if not found.
        :param default: The `default` parameter in the `get_access` function is used to specify a default
        access level to return if the access level is not found in either the parameters or the
        configuration. If the access level is not found in both places, the function will return the
        value of the `default` parameter
        :type default: str

        :return: The method `get_access` is returning a string value that represents the access level.
        """

        # Find access in param then in config, otherwise return default
        return self.get_param().get("access", self.get_config().get("access", default))

    def get_contigs(self) -> dict:
        """
        Get contigs from VCF header
        """

        # Get header
        header = self.get_header()

        # Construct dict)
        contigs = {contig: header.contigs[contig] for contig in header.contigs}

        return contigs

    def load_header(
        self,
        header=None,
        table: str = None,
        drop: bool = False,
        view_name: str = "header",
    ) -> str:
        """
        Load header in a table, with INFO, FORMAT, FILTERS, SAMPLES and METADATA

        Args:
            header (vcfobject, optional): VCF object from pyVCF. Defaults to None (header of the Variants object).
            table (str, optional): Table name of the header table. Defaults to None (defined as 'header' later).
            drop (bool, optional): Drop table if exists. Defaults to False.
            view_name (str, optional): Name of the table. Defaults to 'header'.

        Returns:
            str: Name of the table, None otherwise

        """

        def create_header_table(conn):
            """
            Create header table

            Args:
                conn (conn): Database connexion.

            """

            # Columns
            columns = [
                "section VARCHAR",
                "id VARCHAR",
                "number VARCHAR",
                "type VARCHAR",
                "description VARCHAR",
            ]

            # Query create
            query_create = f"""
            CREATE OR REPLACE table {view_name} (
                {', '.join(columns)},
                PRIMARY KEY (section, id)
            );
            """

            # Execute
            conn.execute(query_create)

        def insert_header(conn, vcf_header):
            """
            Insert header into table

            Args:
                conn (conn): Database connexion.

            """

            # Init
            inserts = []

            # Add INFO section
            for info_id, info in vcf_header.infos.items():
                inserts.append(
                    (
                        "INFO",
                        info_id,
                        str(info.num if info.num is not None else "."),
                        info.type if info.type is not None else "",
                        info.desc if info.desc is not None else "",
                    )
                )

            # Add FORMAT section
            for format_id, format in vcf_header.formats.items():
                inserts.append(
                    (
                        "FORMAT",
                        format_id,
                        str(format.num if format.num is not None else "."),
                        format.type if format.type is not None else "",
                        format.desc if format.desc is not None else "",
                    )
                )

            # Add FILTER section
            for filter_id, filter in vcf_header.filters.items():
                inserts.append(
                    (
                        "FILTER",
                        filter_id,
                        "",
                        "",
                        filter.desc if filter.desc is not None else "",
                    )
                )

            # Add Samples
            for sample_id in vcf_header.samples:
                inserts.append(
                    (
                        "SAMPLE",
                        sample_id,
                        "",
                        "",
                        "",
                    )
                )

            # Add Metadata
            for key, value in vcf_header.metadata.items():
                inserts.append(
                    (
                        "METADATA",
                        key,
                        "",
                        "",
                        str(value) if value is not None else "",
                    )
                )

            # Create query of insert with parameters
            query_insert = f"""
            INSERT INTO {view_name} (section, id, number, type, description) VALUES (?, ?, ?, ?, ?);
            """
            conn.executemany(query_insert, inserts)

        # Get header is None
        if header is None:
            header = self.get_header()

        # Header table
        if table is None:
            table = "header"

        # If header is not None
        if header is not None:

            # Connexion
            conn = self.get_connexion()

            # Drop table
            if drop:
                query_drop = f"""
                DROP TABLE IF EXISTS {table}
                """
                conn.execute(query_drop)

            # Create table
            create_header_table(conn)
            insert_header(conn, header)

            return table

        else:

            return None

    def set_samples(self, samples: list = None) -> list:
        """
        The function `set_samples` sets the samples attribute of an object to a provided list or
        retrieves it from a parameter dictionary.

        :param samples: The `set_samples` method is a method of a class that takes a list of samples as
        input and sets the `samples` attribute of the class to the provided list. If no samples are
        provided, it tries to get the samples from the class's parameters using the `get_param` method
        :type samples: list
        :return: The `samples` list is being returned.
        """

        if not samples:
            samples = self.get_param().get("samples", {}).get("list", None)

        self.samples = samples

        return samples

    def get_samples(self) -> list:
        """
        This function returns a list of samples.
        :return: The `get_samples` method is returning the `samples` attribute of the object.
        """

        return self.samples

    def get_samples_check(self) -> bool:
        """
        This function returns the value of the "check" key within the "samples" dictionary retrieved
        from the parameters.
        :return: The method `get_samples_check` is returning the value of the key "check" inside the
        "samples" dictionary, which is nested inside the dictionary returned by the `get_param()`
        method. If the key "check" is not found, it will return `False`.
        """

        return self.get_param().get("samples", {}).get("check", True)

    def set_input(self, input: str = None) -> None:
        """
        The function `set_input` takes a file name as input, extracts the name and extension, and sets
        attributes in the class accordingly.

        :param input: The `set_input` method in the provided code snippet is used to set attributes
        related to the input file. Here's a breakdown of the parameters and their usage in the method:
        :type input: str
        """

        if input and not isinstance(input, str):
            try:
                self.input = input.name
            except:
                log.error(f"Input file '{input} in bad format")
                raise ValueError(f"Input file '{input} in bad format")
        else:
            self.input = input

        # Input format
        if input:
            input_name, input_extension = os.path.splitext(self.input)
            self.input_name = input_name
            self.input_extension = input_extension
            self.input_format = self.input_extension.replace(".", "")

    def set_config(self, config: dict) -> None:
        """
        The set_config function takes a config object and assigns it as the configuration object for the
        class.

        :param config: The `config` parameter in the `set_config` function is a dictionary object that
        contains configuration settings for the class. When you call the `set_config` function with a
        dictionary object as the argument, it will set that dictionary as the configuration object for
        the class
        :type config: dict
        """

        self.config = config

    def set_param(self, param: dict) -> None:
        """
        This function sets a parameter object for the class based on the input dictionary.

        :param param: The `set_param` method you provided takes a dictionary object as input and sets it
        as the `param` attribute of the class instance
        :type param: dict
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
            "contains": "SIMILAR TO",
        }

        self.code_type_map = {"Integer": 0, "String": 1, "Float": 2, "Flag": 3}

        self.code_type_map_to_sql = code_type_map_to_sql

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
        # if config.get("memory", None):
        #     connexion_config["memory_limit"] = config.get("memory")
        if self.get_memory():
            connexion_config["memory_limit"] = self.get_memory()

        # Temporary directory
        if config.get("temp_directory", None):
            connexion_config["temp_directory"] = config.get("temp_directory")
        else:
            config["temp_directory"] = os.path.join(
                self.get_tmp_dir(), f"duckdb_temp_" + get_random(10)
            )
            self.set_config(config)
            connexion_config["temp_directory"] = config["temp_directory"]

        # Access
        access = self.get_access(default=None)
        if access:
            access_mode = access
            if access in ["RO"]:
                access_mode = "READ_ONLY"
            elif access in ["RW"]:
                access_mode = "READ_WRITE"
            connexion_db = self.get_connexion_db()
            if connexion_db in ":memory:":
                access_mode = "READ_WRITE"
            connexion_config["access_mode"] = access_mode

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
            duckdb_settings = full_path(duckdb_settings)
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
            connexion_db = self.get_tmp_dir() + f"/howard.{get_random()}.tmp.db"
        elif self.get_connexion_type() != "":
            connexion_db = self.get_connexion_type()
        else:
            connexion_db = default_connexion_db

        # Set connexion db
        self.connexion_db = connexion_db

        return connexion_db

    def set_connexion(self, conn) -> None:
        """
        The function `set_connexion` creates a connection to a database, with options for different
        database formats and settings.

        :param conn: The `conn` parameter in the `set_connexion` method is the connection to the
        database. If a connection is not provided, a new connection to an in-memory database is created.
        The method then proceeds to set up the connection based on the specified format (e.g., duckdb or
        sqlite
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
                # duckDB settings arrow large buffer size
                conn.execute("SET arrow_large_buffer_size=true")
                # settings = conn.execute("SELECT * FROM duckdb_settings()").df()
                # log.debug(f"DuckDB settings after connexion:\n{settings.to_string()}")
            elif connexion_format in ["sqlite"]:
                conn = sqlite3.connect(connexion_db)

        # Set connexion
        self.conn = conn

        # Log
        log.debug(f"connexion_format: {connexion_format}")
        log.debug(f"connexion_db: {connexion_db}")
        log.debug(f"connexion config: {connexion_config}")
        log.debug(f"connexion duckdb settings: {self.get_duckdb_settings()}")
        log.debug("connexion duckdb settings: arrow_large_buffer_size=true")

    def set_output(self, output: str = None) -> None:
        """
        The `set_output` function in Python sets the output file based on the input or a specified key
        in the config file, extracting the output name, extension, and format.

        :param output: The `output` parameter in the `set_output` method is used to specify the name of
        the output file. If the config file has an 'output' key, the method sets the output to the value
        of that key. If no output is provided, it sets the output to `None`
        :type output: str
        """

        if output and not isinstance(output, str):
            self.output = output.name
        else:
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
            "##fileformat=VCFv4.2",
            "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO",
        ]

        # Full path
        input_file = full_path(input_file)

        if input_file:

            input_format = self.get_input_format()
            input_compressed = self.get_input_compressed()
            config = self.get_config()
            header_list = default_header_list
            if input_format in [
                "vcf",
                "hdr",
                "tsv",
                "csv",
                "psv",
                "parquet",
                "db",
                "duckdb",
                "json",
            ]:
                # header provided in param
                if config.get("header_file", None):
                    with open(config.get("header_file"), "rt") as f:
                        header_list = self.read_vcf_header(f)
                # within a vcf file format (header within input file itsself)
                elif input_format in ["vcf", "hdr"] and not os.path.isdir(input_file):
                    # within a compressed vcf file format (.vcf.gz)
                    if input_compressed:
                        with bgzf.open(input_file, "rt") as f:
                            header_list = self.read_vcf_header(f)
                    # within an uncompressed vcf file format (.vcf)
                    else:
                        with open(input_file, "rt") as f:
                            header_list = self.read_vcf_header(f)
                # header provided in default external file .hdr
                elif os.path.exists((input_file + ".hdr")):
                    with open(input_file + ".hdr", "rt") as f:
                        header_list = self.read_vcf_header(f)
                else:
                    try:  # Try to get header info fields and file columns

                        with tempfile.TemporaryDirectory() as tmpdir:

                            # Create database
                            db_for_header = Database(database=input_file)

                            # Get header columns for infos fields
                            db_header_from_columns = (
                                db_for_header.get_header_from_columns()
                            )

                            # Get real columns in the file
                            db_header_columns = db_for_header.get_columns()

                            # Write header file
                            header_file_tmp = os.path.join(tmpdir, "header")
                            f = open(header_file_tmp, "w")
                            vcf.Writer(f, db_header_from_columns)
                            f.close()

                            # Replace #CHROM line with rel columns
                            header_list = db_for_header.read_header_file(
                                header_file=header_file_tmp
                            )
                            header_list[-1] = "\t".join(db_header_columns)

                    except:

                        log.warning(
                            f"No header for file {input_file}. Set as default VCF header"
                        )
                        header_list = default_header_list

            else:  # try for unknown format ?

                log.error(f"Input file format '{input_format}' not available")
                raise ValueError(f"Input file format '{input_format}' not available")

            if not header_list:
                header_list = default_header_list

            # header as list
            self.header_list = header_list

            # header as VCF object
            self.header_vcf = vcf.Reader(io.StringIO("\n".join(header_list)))

        else:

            self.header_list = None
            self.header_vcf = None

    def get_query_to_df(self, query: str = "", limit: int = None) -> pd.DataFrame:
        """
        The `get_query_to_df` function takes a query as a string and returns the result as a pandas
        DataFrame based on the connection format. It supports both limited and full queries.

        :param query: The `query` parameter in the `get_query_to_df` function is a string that
        represents the SQL query you want to execute. This query will be used to fetch data from a
        database and convert it into a pandas DataFrame
        :type query: str
        :param limit: The `limit` parameter in the `get_query_to_df` function is used to specify the
        maximum number of rows to be returned in the resulting dataframe. If a limit is provided, the
        function will only fetch up to that number of rows from the database query result. If no limit
        is specified,
        :type limit: int
        :return: A pandas DataFrame is being returned by the `get_query_to_df` function.
        """

        # Connexion format
        connexion_format = self.get_connexion_format()

        # Limit in query
        if limit:

            # Panda settings
            pd.set_option("display.max_rows", limit)

            # DuckDB connexion
            if connexion_format in ["duckdb"]:

                result = self.get_connexion().execute(query).fetch_record_batch(limit)
                if result is None:
                    df = result.df()
                else:
                    try:
                        df = result.read_next_batch().to_pandas()
                    except StopIteration:
                        df = self.get_connexion().execute(query).df()[0:limit]

            # SQLite connexion
            elif connexion_format in ["sqlite"]:
                df = next(pd.read_sql_query(query, self.conn, chunksize=limit))

        # Full query without limit
        else:

            # DuckDB connexion
            if connexion_format in ["duckdb"]:
                df = self.get_connexion().execute(query).df()

            # SQLite connexion
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
        log.info(
            "Input:  "
            + str(self.get_input())
            + " ["
            + str(str(self.get_input_format()))
            + "]"
        )
        log.info(
            "Output: "
            + str(self.get_output())
            + " ["
            + str(str(self.get_output_format()))
            + "]"
        )
        log.info("Config: ")
        for d in str(json.dumps(self.get_config(), indent=4, sort_keys=True)).split(
            "\n"
        ):
            log.info("\t" + str(d))
        log.info("Param: ")
        for d in str(json.dumps(self.get_param(), indent=4, sort_keys=True)).split(
            "\n"
        ):
            log.info("\t" + str(d))
        log.info("Sample list: " + str(self.get_header_sample_list()))
        log.info("Dataframe: ")
        for d in str(df).split("\n"):
            log.info("\t" + str(d))

        # garbage collector
        del df
        gc.collect()

        return None

    def get_stats(
        self,
        table: str = None,
        table_view: str = None,
        annotations_stats: bool = False,
        queries: dict = None,
        queries_view: str = None,
    ) -> dict:
        """
        Calculate and return various statistics of the current object, including information about the input file,
        variants, samples, header fields, quality, and SNVs/InDels.

        :param table: The name of the table containing variant data. If not provided, the default table is used.
        :type table: str, optional
        :param table_view: The name of the table view to be used for statistics calculation. If not provided, a new view is created.
        :type table_view: str, optional
        :param annotations_stats: Whether to calculate annotation statistics. Defaults to False.
        :type annotations_stats: bool, optional
        :param queries: The `queries` parameter is a dictionary that contains queries to be executed
        and added to the statistics. The keys of the dictionary are the names of the queries, and the
        values are the SQL queries to be executed.
        :type queries: dict
        :param queries_view: The `queries_view` parameter is a string that represents the name of the
        view to be used for the queries. If no value is provided, a new view will be created.
        :type queries_view: str

        :return: A dictionary containing various statistics of the current object. The dictionary has the following structure:

            - **Infos** (*dict*): General information about the input file and header fields.
                - **Input file** (*str*): The path to the input file.
                - **Header Infos** (*list*): List of INFO fields in the header.
                - **Header Formats** (*list*): List of FORMAT fields in the header.
                - **Number of INFO fields** (*int*): Number of INFO fields in the header.
                - **Number of FORMAT fields** (*int*): Number of FORMAT fields in the header.
                - **Number of samples** (*int*): Number of samples in the dataset.
                - **Number of variants** (*int*): Total number of variants in the dataset.

            - **Variants** (*dict*): Statistics about the variants.
                - **By Chromosome** (*list*): List of dictionaries with chromosome names and variant counts.
                - **By Type** (*list*): List of dictionaries with variant types and counts.
                - **By Quality** (*list*): List of dictionaries with quality scores and counts.
                - **By Filter** (*list*): List of dictionaries with filter values and counts.

            - **Samples** (*dict*): Statistics about the samples.
                - **Variants in samples** (*list*): List of dictionaries with sample names, variant counts, and percentages.

            - **Header** (*dict*): Detailed information about the header fields.
                - **List of INFO fields** (*dict*): Dictionary with detailed information about INFO fields.
                - **List of FORMAT fields** (*dict*): Dictionary with detailed information about FORMAT fields.
                - **List of FILTER fields** (*dict*): Dictionary with detailed information about FILTER fields.

            - **Annotations** (*dict*, optional): Annotation statistics, if `annotations_stats` is True.
                - **Stats** (*dict*): Dictionary with annotation statistics.

            - **Quality** (*dict*): Quality statistics.
                - **QUAL** (*dict*): Dictionary with quality statistics (average, minimum, maximum, standard deviation, median, variance).

        :rtype: dict

        :example:

        .. code-block:: python

            stats = get_stats(table="variants_table", table_view="variants_view", annotations_stats=True)
            print(stats)
        """

        # Log
        log.info(f"Stats Calculation...")

        # table variants
        if table is None:
            table_variants_from = self.get_table_variants()
        else:
            table_variants_from = table

        # table view
        if table_view is None:
            variants_view_stats_name = "variants_view_stats_" + get_random()
        else:
            variants_view_stats_name = table_view

        # Tables to remove
        tables_to_remove = []

        # Percent_round
        percent_round = 2

        # Sample struct column name
        sample_struct_column = "SAMPLES"

        # Info struct column name
        info_struct_column = None

        # Sample struct column format needed
        if annotations_stats:
            info_prefix_column = ""
        else:
            info_prefix_column = None

        # Create view
        variants_view_stats_name = self.create_annotations_view(
            table=table_variants_from,
            view=variants_view_stats_name,
            view_type="table",
            view_mode="full",
            info_prefix_column=info_prefix_column,
            info_struct_column=info_struct_column,
            sample_struct_column=sample_struct_column,
            formats=["GT"],
            fields_needed=[
                "#CHROM",
                "POS",
                "REF",
                "ALT",
                "QUAL",
                "FILTER",
            ]
            + self.get_header_sample_list(),
        )
        tables_to_remove.append(variants_view_stats_name)

        # stats dict
        stats = {"Infos": {}}

        ### File
        input_file = self.get_input()
        stats["Infos"]["Input file"] = input_file

        # Header
        header_infos = self.get_header().infos
        header_formats = self.get_header().formats
        header_infos_list = list(header_infos)
        header_formats_list = list(header_formats)
        header_table = self.load_header()

        # Stat section
        stats["Stats"] = {}

        ### Variants

        # Variants by chr
        sql_query_nb_variant_by_chrom = f'SELECT "#CHROM" as CHROM, count(*) as count FROM {variants_view_stats_name} GROUP BY "#CHROM"'
        df_nb_of_variants_by_chrom = self.get_query_to_df(sql_query_nb_variant_by_chrom)
        nb_of_variants_by_chrom = df_nb_of_variants_by_chrom.sort_values(
            by=["CHROM"], kind="quicksort"
        )

        # Total number of variants
        nb_of_variants = nb_of_variants_by_chrom["count"].sum()

        # Calculate percentage
        nb_of_variants_by_chrom["percent"] = nb_of_variants_by_chrom["count"].apply(
            lambda x: round((x * 100 / nb_of_variants), percent_round)
        )

        # Add to stats dict the number of variants by chromosome and the total number of variants
        stats["Stats"]["Variants by chromosome"] = nb_of_variants_by_chrom.to_dict(
            orient="index"
        )

        # Add to stats dict the total number of variants
        stats["Infos"]["Number of variants"] = int(nb_of_variants)

        ### Samples

        # Init
        samples = {}
        nb_of_samples = 0

        # Check Samples
        if "GT" in header_formats_list and "FORMAT" in self.get_header_columns():
            log.debug(f"Check samples...")

            # Samples stats
            samples_stats = {}

            # Get samples stats by genotype for each sample in the header
            for sample in self.get_header_sample_list():
                sql_query_samples = f"""
                    SELECT 
                        '{sample}' as 'sample',
                        SAMPLES."{sample}".GT as 'genotype',
                        count(SAMPLES."{sample}".GT) as 'count',
                        ROUND((count(SAMPLES."{sample}".GT)*100/{nb_of_variants}), {percent_round}) as 'percent'
                    FROM {variants_view_stats_name}
                    WHERE  SAMPLES."{sample}".GT IS NOT NULL
                    GROUP BY genotype
                    ORDER BY genotype
                """

                # Get samples stats by genotype for each sample in the header
                sql_query_genotype_df = self.get_connexion().execute(sql_query_samples).df()
                non_null_genotypes = sql_query_genotype_df[
                    sql_query_genotype_df["genotype"].str.contains(r"\d")
                ]
                sample_genotype_count = non_null_genotypes["count"].sum()

                # Add to samples dict the samples stats by genotype for each sample in the header
                if len(sql_query_genotype_df):

                    # Number of samples
                    nb_of_samples += 1

                    # Add to samples dict the samples stats by genotype for each sample in the header
                    samples[sample] = sql_query_genotype_df.to_dict(orient="index")

                    # Add to samples stats dict the samples stats by genotype for each sample in the header
                    samples_stats[sample] = {
                        "Sample": f"{sample}",
                        "count": int(sample_genotype_count),
                        "percent": round(
                            (sample_genotype_count * 100 / nb_of_variants),
                            percent_round,
                        ),
                    }

            # Add to stats dict the samples stats by genotype for each sample in the header
            stats["Samples"] = samples
            stats["Infos"]["Number of samples"] = nb_of_samples
            stats["Stats"]["Variants by sample"] = samples_stats

        else:

            samples_stats = {}

        ### INFO and FORMAT fields
        header_types_df = {}
        header_types_list = {
            "INFO": {
                "label": "List of INFO fields",
                "fields": {
                    "id": "INFO",
                    "number": "Number",
                    "type": "Type",
                    "description": "Description",
                },
            },
            "FORMAT": {
                "label": "List of FORMAT fields",
                "fields": {
                    "id": "FORMAT",
                    "number": "Number",
                    "type": "Type",
                    "description": "Description",
                },
            },
            "FILTER": {
                "label": "List of FILTER fields",
                "fields": {
                    "id": "FILTER",
                    "description": "Description",
                },
            },
        }

        # Init
        header_types_df = {}

        # Get header types for INFO and FORMAT fields
        for header_section, header_info in header_types_list.items():
            label = header_info["label"]
            fields = header_info["fields"]

            # Construire la liste des champs à sélectionner
            select_fields = ", ".join(fields.keys())

            # SQL query
            sql_query_header = f"""
                SELECT {select_fields}
                FROM {header_table}
                WHERE section = '{header_section}'
            """
            header_infos_df = self.get_query_to_df(sql_query_header)
            header_infos_dict = {}

            # Add to header_types_df the header types for INFO and FORMAT fields
            for i, row in header_infos_df.iterrows():
                header_infos_dict[i] = {
                    new if new else original: row[original]
                    for original, new in fields.items()
                }

            # Add to header_types_df the header types for INFO and FORMAT fields
            if len(header_infos_dict):

                # Add to header_types_df the header types for INFO and FORMAT fields
                header_types_df[label] = pd.DataFrame.from_dict(
                    header_infos_dict, orient="index"
                ).to_dict(orient="index")

                # Add to stats dict the number of INFO and FORMAT fields
                stats["Infos"][f"Number of {header_section} fields"] = len(
                    header_types_df[label]
                )

        # Add to stats dict the header types for INFO and FORMAT fields
        stats["Header"] = header_types_df

        # Annotations stats
        if annotations_stats:

            # Init
            sql_queries_info = []

            # Get header infos list
            for field in header_infos_list:

                # Create a table with a field by line (only for INFO section), and le number of distinct value on variants table, and the number of variants with a value
                sql_queries_info.append(f"""
                        SELECT
                            '{field}' AS 'Annotation',
                            count(distinct "{field}") as 'Distinct values',
                            count("{field}") as 'Annotated Variants',
                            ROUND((count("{field}") * 100 / {nb_of_variants}), {percent_round}) as 'Percent',
                        FROM
                            {variants_view_stats_name}
                        WHERE
                            "{field}" IS NOT NULL AND TRIM(CAST("{field}" AS VARCHAR)) NOT IN ('','.')
                    """)

            # Join all queries
            sql_query_info = f""" UNION ALL """.join(sql_queries_info)

            # Get info stats
            info_stats = self.get_query_to_df(sql_query_info)

            # Add to stats dict the annotations stats
            stats["Annotations"] = {"Distribution": info_stats.to_dict(orient="index")}

        ### Quality stats
        log.debug(f"Quality stats...")

        ### QUAL
        log.debug(f"Quality stats: QUAL...")

        if "QUAL" in self.get_header_columns():

            # SQL query
            sql_query_qual = f"""
                    SELECT
                        avg(CAST(QUAL AS INTEGER)) AS Average,
                        min(CAST(QUAL AS INTEGER)) AS Minimum,
                        max(CAST(QUAL AS INTEGER)) AS Maximum,
                        stddev(CAST(QUAL AS INTEGER)) AS StandardDeviation,
                        median(CAST(QUAL AS INTEGER)) AS Median,
                        variance(CAST(QUAL AS INTEGER)) AS Variance
                    FROM {variants_view_stats_name}
                    WHERE CAST(QUAL AS VARCHAR) NOT IN ('.')
                    """

            # Get quality stats
            qual_stats = self.get_connexion().execute(sql_query_qual).df().to_dict(orient="index")

        else:

            # Empty quality stats
            qual_stats = {}

        ### FILTER
        log.debug(f"Quality stats: FILTER...")

        if "FILTER" in self.get_header_columns():

            # SQL query
            sql_query_filter = f"""
                WITH split_filter AS (
                    SELECT
                        TRIM(UNNEST(STRING_SPLIT(CASE WHEN TRIM(FILTER) = '' OR FILTER IS NULL THEN '.' ELSE FILTER END, ';'))) AS filter_value
                    FROM
                        {variants_view_stats_name}
                )
                SELECT
                    filter_value,
                    COUNT(*) AS 'count',
                    ROUND((count * 100 / {nb_of_variants}), {percent_round}) AS 'percent'
                FROM
                    split_filter
                GROUP BY
                    filter_value
                ORDER BY
                    count DESC
            """

            # Get filter stats
            filter_stats = (
                self.get_connexion().execute(sql_query_filter).df().to_dict(orient="index")
            )

        else:

            # Empty filter stats
            filter_stats = {}

        ### SNV and InDel

        # SQL query
        sql_query_snv = f"""
            
            SELECT Type, count, ROUND((count * 100 / {nb_of_variants}), {percent_round}) AS 'percent' FROM (

                    SELECT
                        'Total' AS Type,
                        count(*) AS count
                    FROM {variants_view_stats_name}
                    
                    UNION

                    SELECT
                        'SNV' AS Type,
                        count(*) AS count
                    FROM {variants_view_stats_name}
                    WHERE len(REF) = 1 AND len(ALT) = 1

                    UNION

                    SELECT
                        'MNV' AS Type,
                        count(*) AS count
                    FROM {variants_view_stats_name}
                    WHERE len(REF) > 1 AND len(ALT) > 1
                    AND len(REF) = len(ALT)

                    UNION

                    SELECT
                        'InDel' AS Type,
                        count(*) AS count
                    FROM {variants_view_stats_name}
                    WHERE len(REF) > 1 OR len(ALT) > 1
                    AND len(REF) != len(ALT)

                )

            ORDER BY 
            CASE
                WHEN Type = 'Total' THEN 1
                WHEN Type = 'SNV' THEN 2
                WHEN Type = 'MNV' THEN 3
                WHEN Type = 'InDel' THEN 4
            END

                """

        # Get SNV and InDel stats
        snv_indel = self.get_query_to_df(sql_query_snv).to_dict(orient="index")

        # Substitutions
        sql_query_snv_substitution = f"""
                SELECT
                    concat(REF, '>', ALT) AS 'Substitution',
                    count(*) AS count,
                    ROUND((count * 100 / {nb_of_variants}), {percent_round}) AS 'percent'
                FROM {variants_view_stats_name}
                WHERE len(REF) = 1 AND len(ALT) = 1
                GROUP BY REF, ALT
                ORDER BY count(*) DESC
                """
        snv_substitution = self.get_query_to_df(sql_query_snv_substitution).to_dict(
            orient="index"
        )

        # Add to stats dict the SNV and InDel stats
        stats["Stats"]["Variant types"] = snv_indel
        stats["Stats"]["Substitutions"] = snv_substitution
        stats["Stats"]["Quality"] = qual_stats
        stats["Stats"]["Filters"] = filter_stats

        # Queries
        if queries is not None:

            # Create full annotations view
            variants_view_query_name = queries_view
            variants_view_query_name = self.create_annotations_view(
                table=table_variants_from,
                view=variants_view_query_name,
                view_type="view",
                view_mode="explore",
                info_prefix_column="",
                info_struct_column="INFOS",
                sample_struct_column="SAMPLES",
                formats=None,
                fields_needed_all=True,
                drop_view=True,
            )
            if queries_view is not None:
                tables_to_remove.append(variants_view_query_name)

            # Stats queries section
            stats["Queries"] = {}

            # For each query
            for query_infos in queries.items():

                # Query name and query
                query_name = query_infos[0]
                query = query_infos[1]

                # Query cast
                query_cast = cast_columns_query(query=query, conn=self.get_connexion())

                # Query execute
                query_res = self.get_query_to_df(query_cast).to_dict(orient="index")
                stats["Queries"][query_name] = query_res

        # Remove table or view
        self.remove_tables_or_views(tables=tables_to_remove)

        return stats

    def stats_to_file(
        self,
        file: str = None,
        annotations_stats: bool = False,
        queries: dict = None,
        queries_view: str = None,
    ) -> str:
        """
        The function `stats_to_file` takes a file name as input, retrieves statistics, serializes them
        into a JSON object, and writes the JSON object to the specified file.

        :param file: The `file` parameter is a string that represents the file path where the JSON data
        will be written
        :type file: str
        :param annotations_stats: The `annotations_stats` parameter is a boolean that specifies whether
        to calculate annotation statistics. If `annotations_stats` is set to True, annotation statistics
        will be calculated. If `annotations_stats` is set to False, annotation statistics will not be
        calculated. The default value is False.
        :type annotations_stats: bool
        :param queries: The `queries` parameter is a dictionary that contains queries to be executed
        and added to the statistics. The keys of the dictionary are the names of the queries, and the
        values are the SQL queries to be executed.
        :type queries: dict
        :param queries_view: The `queries_view` parameter is a string that represents the name of the
        view to be used for the queries. If no value is provided, a new view will be created.
        :type queries_view: str

        :return: The name of the file that was written to.
        """

        # Get stats
        stats = self.get_stats(
            annotations_stats=annotations_stats,
            queries=queries,
            queries_view=queries_view,
        )

        # Serializing json
        json_object = json.dumps(stats, indent=4)

        # Writing to sample.json
        with open(file, "w") as outfile:
            outfile.write(json_object)

        return file

    def print_stats(
        self,
        stdout: bool = False,
        output_file: str = None,
        json_file: str = None,
        html_file: str = None,
        pdf_file: str = None,
        annotations_stats: bool = False,
        queries: dict = None,
        queries_view: str = None,
    ) -> None:
        """
        The `print_stats` function generates a markdown file and prints the statistics contained in a
        JSON file in a formatted manner.

        :param stdout: The `stdout` parameter is a boolean that specifies whether to print the stats
        directly to the standard output. If `stdout` is set to True, the stats will be printed to the
        standard output. If `stdout` is set to False, the stats will not be printed to the standard
        output. The default value is False.
        :type stdout: bool
        :param output_file: The `output_file` parameter is a string that specifies the path and filename
        of the output file where the stats will be printed in Markdown format. If no `output_file` is
        provided, a temporary directory will be created and the stats will be saved in a file named
        "stats.md" within that
        :type output_file: str
        :param json_file: The `json_file` parameter is a string that represents the path to the JSON
        file where the statistics will be saved. If no value is provided, a temporary directory will be
        created and a default file name "stats.json" will be used
        :type json_file: str
        :param html_file: The `html_file` parameter is a string that specifies the path and filename of
        the output file where the stats will be printed in HTML format. If no `html_file` is provided,
        a temporary directory will be created and the stats will be saved in a file named "stats.html"
        within that
        :type html_file: str
        :param pdf_file: The `pdf_file` parameter is a string that specifies the path and filename of the
        output file where the stats will be printed in PDF format. If no `pdf_file` is provided, a
        temporary directory will be created and the stats will be saved in a file named "stats.pdf"
        within that
        :type pdf_file: str
        :param annotations_stats: Whether to calculate annotation statistics. Defaults to False.
        :type annotations_stats: bool, optional
        :param queries: The `queries` parameter is a dictionary that contains queries to be executed
        and added to the statistics. The keys of the dictionary are the names of the queries, and the
        values are the SQL queries to be executed.
        :type queries: dict
        :param queries_view: The `queries_view` parameter is a string that represents the name of the
        view to be used for the queries. If no value is provided, a new view will be created.
        :type queries_view: str
        of `None`.

        :return: The function `print_stats` does not return any value. It has a return type annotation

        """

        # Full path
        output_file = full_path(output_file)
        json_file = full_path(json_file)

        # Create stats file in temporary directory
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
            stats_file = self.stats_to_file(
                file=json_file,
                annotations_stats=annotations_stats,
                queries=queries,
                queries_view=queries_view,
            )

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
            output_index.append("## Table of context")

            # Process sections
            for section in stats:
                infos = stats.get(section)
                section_link = "#" + section.lower().replace(" ", "-")
                output.append(f"\n")
                output.append(f"## {section}")
                output_index.append(f"- [{section}]({section_link})")

                if len(infos):

                    # For each info
                    for info in infos:

                        # Check if dataframe or not
                        try:
                            df = pd.DataFrame.from_dict(infos.get(info), orient="index")
                            is_df = True
                        except:
                            try:
                                df = pd.DataFrame.from_dict(
                                    json.loads((infos.get(info))), orient="index"
                                )
                                is_df = True
                            except:
                                is_df = False

                        # If dataframe is a dataframe
                        if is_df:
                            df = df.map(escape_markdown_table_chars)
                            output.append(f"### {info}")
                            info_link = "#" + info.lower().replace(" ", "-")
                            output_index.append(f"   - [{info}]({info_link})")
                            output.append(f"{df.to_markdown(index=False)}")

                        # If not a dataframe
                        else:
                            output.append(f"- {info}: {infos.get(info)}")

                else:

                    # If no info
                    output.append(f"NA")

            # Write stats in markdown file
            with open(output_file, "w") as fp:
                for item in output_title:
                    fp.write("%s\n" % item)
                fp.write("\n")
                for item in output_index:
                    fp.write("%s\n" % item)
                fp.write("\n")
                for item in output:
                    fp.write("%s\n" % item)

            # Output stats in markdown
            if stdout:
                print("")
                print("\n\n".join(output_title))
                print("")
                print("\n\n".join(output))
                print("")

            # Generate HTML and PDF files
            if html_file:
                convert_markdown_to_html(output_file, html_file)
            if pdf_file:
                convert_markdown_to_pdf(output_file, pdf_file)

        return None

    def get_input(self) -> str:
        """
        It returns the value of the input variable.
        :return: The input is being returned.
        """
        return self.input

    def get_input_format(self, input_file: str = None) -> str:
        """
        This function returns the format of the input variable, either from the provided input file or
        by prompting for input.

        :param input_file: The `input_file` parameter in the `get_input_format` method is a string that
        represents the file path of the input file. If no `input_file` is provided when calling the
        method, it will default to `None`
        :type input_file: str
        :return: The format of the input variable is being returned.
        """

        if not input_file:
            input_file = self.get_input()
        input_format = get_file_format(input_file)
        return input_format

    def get_input_compressed(self, input_file: str = None) -> str:
        """
        The function `get_input_compressed` returns the format of the input variable after compressing
        it.

        :param input_file: The `input_file` parameter in the `get_input_compressed` method is a string
        that represents the file path of the input file. If no `input_file` is provided when calling the
        method, it will default to `None` and the method will then call `self.get_input()` to
        :type input_file: str
        :return: The function `get_input_compressed` returns the compressed format of the input
        variable.
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
        The function `get_output_format` returns the format of the input variable or the output file if
        provided.

        :param output_file: The `output_file` parameter in the `get_output_format` method is a string
        that represents the file path of the output file. If no `output_file` is provided when calling
        the method, it will default to the output obtained from the `get_output` method of the class
        instance. The
        :type output_file: str
        :return: The format of the input variable is being returned.
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
        access = self.get_access(default=None)

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

    def close_connexion(self) -> str:
        """
        This function closes the connection to the database.
        :return: The connection is being closed.
        """

        log.debug(f"Close connexion...")
        self.get_connexion().close()

        connexion_db = self.get_connexion_db()

        # Remove connexion db file
        if os.path.exists(connexion_db) and connexion_db != self.get_connexion_type():
            log.debug(f"Remove connexion db file: {connexion_db}")
            remove_if_exists([connexion_db])

        log.debug(f"Connexion '{connexion_db}' closed.")

        # Clean temporary files
        tmp_files = self.clean_tmp_files()
        log.debug(f"Temporary files: {len(tmp_files)} files removed.")

        return connexion_db

    def get_header(self, type: str = "vcf", vcf_file: str = None) -> list:
        """
        This function returns the header of the VCF file as a list of strings

        :param type: the type of header you want to get, defaults to vcf (optional)
        :return: The header of the vcf file.
        """

        if vcf_file:
            vcf_file_header = self.read_vcf_header_file(file=vcf_file)
            if type == "vcf":
                header = vcf.Reader(io.StringIO("\n".join(vcf_file_header)))
                return header
            elif type == "list":
                return vcf_file_header
        elif self.header_vcf:
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

    def get_header_infos_list(self) -> list:
        """
        This function retrieves a list of information fields from the header.
        :return: A list of information fields from the header.
        """

        # Init
        infos_list = []

        for field in self.get_header().infos:
            infos_list.append(field)

        return infos_list

    def get_header_length(self, file: str = None) -> int:
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
            return len(self.read_vcf_header_file(file=file)) - 1
        elif self.get_header(type="list"):
            return len(self.get_header(type="list")) - 1
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
            sql_column_list.append(f'"{col}"')
        return ",".join(sql_column_list)

    def get_header_sample_list(
        self, check: bool = False, samples: list = None, samples_force: bool = False, check_format: bool = False
    ) -> list:
        """
        The function `get_header_sample_list` returns a list of samples from a VCF header, with optional
        checking and filtering based on input parameters.

        :param check: The `check` parameter in the `get_header_sample_list` function is a boolean
        parameter that determines whether to check if the samples in the list are properly defined as
        genotype columns. If `check` is set to `True`, the function will verify if each sample in the
        list is defined as a, defaults to False
        :type check: bool (optional)
        :param samples: The `samples` parameter in the `get_header_sample_list` function is a list that
        allows you to specify a subset of samples from the header. If you provide a list of sample
        names, the function will check if each sample is defined in the header. If a sample is not found
        in the
        :type samples: list
        :param samples_force: The `samples_force` parameter in the `get_header_sample_list` function is
        a boolean parameter that determines whether to force the function to return the sample list
        without checking if the samples are genotype columns. If `samples_force` is set to `True`, the
        function will return the sample list without performing, defaults to False
        :type samples_force: bool (optional)
        :param check_format: The `check_format` parameter in the `get_header_sample_list` function is a boolean
        parameter that determines whether to check the format of the genotype columns. If `check_format` is set to `True`,
        the function will verify if each sample's genotype column conforms to the expected format, defaults to False
        :type check_format: bool (optional)
        :return: The function `get_header_sample_list` returns a list of samples based on the input
        parameters and conditions specified in the function.
        """

        # Init
        samples_list = []

        # Determine the initial list of samples based on the provided input and the VCF header
        if samples is None:
            samples_list = self.header_vcf.samples
        else:
            samples_checked = []
            for sample in samples:
                if sample in self.header_vcf.samples:
                    samples_checked.append(sample)
                else:
                    log.warning(f"Sample '{sample}' not defined in header")
            samples_list = samples_checked

            # Force sample list without checking if is_genotype_column
            if samples_force:
                log.warning(f"Samples {samples_list} not checked if genotypes")
                return samples_list

        if check:

            # Check which samples are genotype columns
            samples_checked = self.is_genotype_columns(columns=samples_list, check_format=check_format)

            # Check issues for each sample not conforming
            for sample in samples_list:
                if sample not in samples_checked:
                    log.warning(f"Sample '{sample}' is not a genotype column")
                    try:
                        list_of_non_conforming = Database(database=self.get_input()).is_genotype_column_non_conforming(column=sample)
                        if list_of_non_conforming is False:
                            log.debug(f"Sample '{sample}' does not exist")
                        elif len(list_of_non_conforming) > 0:
                            log.debug(f"Non-conforming genotypes for sample '{sample}':\n{list_of_non_conforming.head()}")
                        else:
                            log.debug(f"Unknown issue with sample '{sample}'. Check FORMAT and genotype consistency.")
                    except Exception as e:
                        log.debug(f"Error checking non-conforming genotype column for sample '{sample}': {e}")


            
            # List of samples that are genotype columns
            samples_list = samples_checked
            #log.debug(f"Samples checked: {samples_list}")

        # Return samples list
        return samples_list

    def sort_contigs(self) -> None:
        """
        This function sort contigs

        :return: None
        """

        # Sort contigs
        header = self.get_header()
        header = sort_contigs(header)

        # Return
        return None

    # def is_genotype_column(self, column: str = None, check_format: bool = False) -> bool:
    #     """
    #     This function checks if a given column is a genotype column in a database.

    #     :param column: The `column` parameter in the `is_genotype_column` method is a string that
    #     represents the column name in a database table. This method checks if the specified column is a
    #     genotype column in the database. If a column name is provided, it calls the `is_genotype_column`
    #     method of
    #     :type column: str
    #     :return: The `is_genotype_column` method is returning a boolean value. If the `column` parameter
    #     is not None, it calls the `is_genotype_column` method of the `Database` class with the specified
    #     column name and returns the result. If the `column` parameter is None, it returns False.
    #     """

    #     if column is not None:
    #         return Database(database=self.get_input()).is_genotype_column(column=column, check_format=check_format)
    #     else:
    #         return False

    def is_genotype_columns(self, columns: list[str] = None, check_format: bool = False) -> list[str]:
        """
        This function checks which columns in a database are genotype columns.

        :param columns: The `columns` parameter is a list of column names in a database table. This method checks which of the specified columns are genotype columns in the database. If no columns are provided, it checks all columns.
        :type columns: list[str]
        :param check_format: Whether to check the FORMAT column for consistency with genotype columns.
        :type check_format: bool
        :return: A list of column names that are genotype columns.
        :rtype: list[str]
        """

        if columns is not None:
            return Database(database=self.get_input()).is_genotype_columns(columns=columns, check_format=check_format)
        else:
            return []

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

    def insert_file_to_table(
        self,
        file,
        columns: str,
        header_len: int = 0,
        sep: str = "\t",
        chunksize: int = 1000000,
    ) -> None:
        """
        The function reads a file in chunks and inserts each chunk into a table based on the specified
        database format.

        :param file: The `file` parameter is the file that you want to load into a table. It should be
        the path to the file on your system
        :param columns: The `columns` parameter in the `insert_file_to_table` function is a string that
        should contain the names of the columns in the table where the data will be inserted. The column
        names should be separated by commas within the string. For example, if you have columns named
        "id", "name
        :type columns: str
        :param header_len: The `header_len` parameter in the `insert_file_to_table` function specifies
        the number of lines to skip at the beginning of the file before reading the actual data. This
        parameter allows you to skip any header information present in the file before processing the
        data, defaults to 0
        :type header_len: int (optional)
        :param sep: The `sep` parameter in the `insert_file_to_table` function is used to specify the
        separator character that is used in the file being read. In this case, the default separator is
        set to `\t`, which represents a tab character. You can change this parameter to a different
        separator character if needed, defaults to \t
        :type sep: str (optional)
        :param chunksize: The `chunksize` parameter specifies the number of rows to read in at a time
        when processing the file in chunks. In the provided code snippet, the default value for
        `chunksize` is set to 1000000. This means that the file will be read in chunks of 1,, defaults
        to 1000000
        :type chunksize: int (optional)
        """

        # Config
        chunksize = self.get_config().get("load", {}).get("chunk", chunksize)
        connexion_format = self.get_connexion_format()

        log.debug("chunksize: " + str(chunksize))

        if chunksize:
            for chunk in pd.read_csv(
                file, skiprows=header_len, sep=sep, chunksize=chunksize, engine="c"
            ):
                if connexion_format in ["duckdb"]:
                    sql_insert_into = (
                        f"INSERT INTO variants ({columns}) SELECT {columns} FROM chunk"
                    )
                    #self.conn.execute(sql_insert_into)
                    self.get_connexion().execute(sql_insert_into)
                elif connexion_format in ["sqlite"]:
                    chunk.to_sql("variants", self.get_connexion(), if_exists="append", index=False)

    def load_data(
        self,
        input_file: str = None,
        drop_variants_table: bool = False,
        sample_size: int = 20480,
    ) -> None:
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
        # access = self.get_config().get("access", None)
        access = self.get_access(default=None)
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

        # DuckDB connexion
        if connexion_format in ["duckdb"]:

            # Database already exists
            if input_format in ["db", "duckdb"]:

                if connexion_format in ["duckdb"]:
                    log.debug(f"Input file format '{input_format}' duckDB")
                else:
                    log.error(
                        f"Input file format '{input_format}' not compatilbe with database format '{connexion_format}'"
                    )
                    raise ValueError(
                        f"Input file format '{input_format}' not compatilbe with database format '{connexion_format}'"
                    )

            # Load from existing database format
            else:

                try:
                    # Create Table or View
                    database = Database(database=self.input)
                    sql_from = database.get_sql_from(sample_size=sample_size)

                    log.debug(f"Load Data into {table_variants}...")
                    if access in ["RO"]:
                        sql_load = (
                            f"CREATE VIEW {table_variants} AS SELECT * FROM {sql_from}"
                        )
                    else:
                        sql_load = (
                            f"CREATE TABLE {table_variants} AS SELECT * FROM {sql_from}"
                        )
                    self.get_connexion().execute(sql_load)
                    log.debug(f"Load Data into {table_variants} - done.")

                except Exception as e:
                    # Format not available
                    msg_err = f"Load Data into {table_variants} - failed to load data: {str(e)}"
                    log.error(msg_err)
                    raise ValueError(msg_err)

        # SQLite connexion
        elif connexion_format in ["sqlite"] and input_format in [
            "vcf",
            "tsv",
            "csv",
            "psv",
        ]:

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
                sql_create_table_columns.append(
                    f'"{column}" {column_type} default NULL'
                )
                sql_create_table_columns_list.append(f'"{column}"')

            # Create database
            log.debug(f"Create Table {table_variants}")
            sql_create_table_columns_sql = ", ".join(sql_create_table_columns)
            sql_create_table_columns_list_sql = ", ".join(sql_create_table_columns_list)
            sql_create_table = f"CREATE TABLE IF NOT EXISTS {table_variants} ({sql_create_table_columns_sql})"
            self.get_connexion().execute(sql_create_table)

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

        else:
            log.error(
                f"Connexion format '{connexion_format}' not available with format '{input_format}'"
            )
            raise ValueError(
                f"Connexion format '{connexion_format}' not available with format '{input_format}'"
            )

        # Add INFO column if not exists
        if access not in ["RO"] and "INFO" not in self.get_header_columns_as_list():
            log.debug("INFO column not found, adding it")
            # Add INFO column
            self.add_column(
                table_name=table_variants,
                column_name="INFO",
                column_type="VARCHAR",
                default_value=None,
            )

        # # Explode INFOS fields into table fields
        # if self.get_explode_infos():
        #     self.explode_infos(
        #         prefix=self.get_explode_infos_prefix(),
        #         fields=self.get_explode_infos_fields(),
        #         force=True,
        #     )

        # Create index after insertion
        self.create_indexes()

    def get_explode_infos(self) -> bool:
        """
        The function `get_explode_infos` returns the value of the "explode_infos" parameter, defaulting
        to False if it is not set.
        :return: The method is returning the value of the "explode_infos" parameter, which is a boolean
        value. If the parameter is not present, it will return False.
        """

        return self.get_param().get("explode", {}).get("explode_infos", False)

    def get_explode_infos_fields(
        self,
        explode_infos_fields: str = None,
        remove_fields_not_in_header: bool = False,
    ) -> list:
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
            explode_infos_fields = (
                self.get_param().get("explode", {}).get("explode_infos_fields", None)
            )

        # If no fields, defined as all fields in header using keyword
        if not explode_infos_fields:
            explode_infos_fields = "*"

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
                    field = ".*"

                # Find all fields with pattern
                # Check if field is in header (to prevent special caracters in field such as '+', e.g. 'GERP++_RS')
                if field in fields_in_header:
                    fields_search = [field]
                else:
                    r = re.compile(rf"^{field}$")
                    fields_search = sorted(list(filter(r.match, fields_in_header)))

                # Remove fields input from search
                if field in fields_search:
                    fields_search = [field]
                elif fields_search != [field]:
                    fields_search = sorted(
                        list(set(fields_search).difference(fields_input))
                    )

                # If field is not in header (avoid not well formatted header)
                if not fields_search and not remove_fields_not_in_header:
                    fields_search = [field]

                # Add found fields
                for new_field in fields_search:
                    # Add field, if not already exists, and if it is in header (if asked)
                    if (
                        new_field not in fields_output
                        and (
                            not remove_fields_not_in_header
                            or new_field in fields_in_header
                        )
                        and new_field not in [".*"]
                    ):
                        fields_output.append(new_field)

            return fields_output

        else:

            return []

    def get_explode_infos_prefix(self, explode_infos_prefix: str = None) -> str:
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
            explode_infos_prefix = (
                self.get_param().get("explode", {}).get("explode_infos_prefix", "")
            )

        return explode_infos_prefix

    def add_column(
        self,
        table_name,
        column_name,
        column_type,
        default_value=None,
        drop: bool = False,
    ) -> dict:
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
        if column_name.upper() in [c.upper() for c in columns]:
            log.debug(
                f"The {column_name} column already exists in the {table_name} table"
            )
            if drop:
                self.drop_column(table_name=table_name, column_name=column_name)
                dropped = True
            else:
                return None
        else:
            log.debug(f"The {column_name} column NOT exists in the {table_name} table")

        # Add column in table
        add_column_query = (
            f""" ALTER TABLE {table_name} ADD COLUMN "{column_name}" {column_type} """
        )
        if default_value is not None:
            add_column_query += f" DEFAULT {default_value}"
        log.debug(f"add_column_query: {add_column_query}")
        self.execute_query(add_column_query)
        added = not dropped
        log.debug(
            f"The {column_name} column was successfully added to the {table_name} table"
        )

        if added:
            added_column = {
                "table_name": table_name,
                "column_name": column_name,
                "column_type": column_type,
                "default_value": default_value,
            }
        else:
            added_column = None

        return added_column

    def drop_column(
        self, column: dict = None, table_name: str = None, column_name: str = None
    ) -> bool:
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
            if isinstance(column, dict):
                table_name = column.get("table_name", None)
                column_name = column.get("column_name", None)
            elif isinstance(column, str):
                table_name = self.get_table_variants()
                column_name = column
            else:
                table_name = None
                column_name = None

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
        log.debug(
            f"The {column_name} column was successfully dropped to the {table_name} table"
        )

        return removed

    def get_batch_split(
        self,
        table: str = None,
        block: int = 1000,
        nb_lines: int = None,
        use_memory: bool = True,
    ) -> int:
        """
        Calculate the batch size for processing data based on the number of rows in the table and available memory.

        Args:
            table (str, optional): The name of the table to evaluate. If None, the default variants table is used.
            block (int, optional): The block size to use for the calculation. Default is 1000.
            nb_lines (int, optional): The number of lines in the table. If None, it will be calculated.
            use_memory (bool, optional): Whether to consider available memory in the calculation to ponderate block size (memory*block). Default is True.

        Returns:
            int: The calculated batch size.
        """

        # Get table variants if no table
        if table is None:
            table = self.get_table_variants()

        # Evaluate split
        log.debug("Evaluate batch size by parameter")

        # Count numbber of variants in table variants
        if nb_lines is None:
            nb_lines = self.get_connexion().execute(f"""
                        SELECT count(1)
                        FROM {table}
                    """).fetchone()[0]

        # Check memory
        if not use_memory:
            memory = 1
        else:
            memory = extract_memory_in_go(
                get_memory(self.get_config(), self.get_param())
            )

        # Avaluate block size using block size (e.g. 1000 viarants) and memory
        block_size = block * memory

        # Calculate batch
        batch = round(nb_lines / block_size) + 1

        # Return
        return batch

    def explode_infos(
        self,
        prefix: str = None,
        create_index: bool = False,
        fields: list = None,
        fields_just_add: list = [],
        fields_not_exists: bool = True,
        detect_type_list: bool = True,
        force: bool = False,
        proccess_all_fields_together: bool = False,
        fields_forced_as_varchar: bool = False,
        table: str = None,
        table_source: str = None,
        table_dest: str = None,
        table_key: list = None,
    ) -> list:
        """
        Explode the INFO fields of a VCF file into individual columns in a specified table.

        Args:
            prefix (str, optional): A prefix for the exploded INFO fields. If not provided, the function
                will use the value of `self.get_explode_infos_prefix()`.
            create_index (bool, optional): Whether to create indexes on the exploded INFO fields. Defaults to False.
            fields (list, optional): A list of INFO fields to explode into individual columns. If not provided,
                all INFO fields will be exploded.
            fields_just_add (list, optional): A list of INFO fields to add as individual columns without exploding values.
            fields_not_exists (bool, optional): Whether to add fields that do not exist in the table. Defaults to True.
            detect_type_list (bool, optional): Whether to detect if the field is a list type. Defaults to True.
            force (bool, optional): Whether to drop and recreate a column if it already exists in the table. Defaults to False.
            proccess_all_fields_together (bool, optional): Whether to process all INFO fields together or individually.
                Defaults to False.
            fields_forced_as_varchar (bool, optional): Whether to force all fields to be treated as VARCHAR. Defaults to False.
            table (str, optional): The name of the table where the exploded INFO fields will be added as individual columns.
            table_source (str, optional): The name of the source table containing the INFO fields.
            table_dest (str, optional): The name of the destination table where the exploded INFO fields will be added.
            table_key (list, optional): A list of keys to use for identifying rows in the table.

        Returns:
            list: A list of added columns.
        """

        # drop indexes
        self.drop_indexes()

        # connexion format
        connexion_format = self.get_connexion_format()
        if connexion_format in ["sqlite"]:
            msg_err = (
                f"Connexion format '{connexion_format}' not available for explode infos"
            )
            log.error(msg_err)
            raise ValueError(msg_err)

        # Access
        # access = self.get_config().get("access", None)
        access = self.get_access(default=None)

        # Added columns
        added_columns = []

        if access not in ["RO"]:

            # Translate fields if patterns
            fields = self.get_explode_infos_fields(explode_infos_fields=fields)

            if fields is None or len(fields) == 0:
                return []

            # prefix
            if prefix in [None, True] or not isinstance(prefix, str):
                if self.get_explode_infos_prefix() not in [None, True]:
                    prefix = self.get_explode_infos_prefix()
                else:
                    prefix = "INFO/"

            # table variants
            if table is None:
                table = self.get_table_variants(clause="select")

            # table source
            if table_source is None:
                table_source = table

            # table dest
            if table_dest is None:
                table_dest = table

            # table key
            if table_key is None:
                table_key = ["#CHROM", "POS", "REF", "ALT"]

            # Check source table columns
            try:
                table_source_struct = self.get_columns(table=table_source)
            except:
                table_source_struct = []
            try:
                table_dest_struct = self.get_columns(table=table_dest)
            except:
                table_dest_struct = []

            if "INFO" not in table_source_struct:
                msg_err = f"Column 'INFO' not found in table '{table_source}'"
                log.warning(msg_err)
                # return None
                # raise ValueError(msg_err)

            # Header infos
            header_infos = self.get_header().infos

            log.debug(
                f"Explode INFO fields - [{len(header_infos)}] annotations fields in header"
            )

            # Create view with all fields
            view_source = "view_source_" + str(random.randint(10000, 100000))
            view_source = self.create_annotations_view(
                table=table_source,
                fields=fields,
                view=view_source,
                view_type="view",
                view_mode="explore",
                info_prefix_column=prefix,
                fields_needed=table_key,
                fields_not_exists=fields_not_exists,
                fields_forced_as_varchar=fields_forced_as_varchar,
                detect_type_list=detect_type_list,
            )

            # Describe view source
            describe_query = f"DESCRIBE {view_source}"
            res = self.execute_query(describe_query)
            description_dict = {row[0]: {"type": row[1]} for row in res.fetchall()}

            # View source structure
            view_source_struct = self.get_columns(table=view_source)

            # Set fields
            sql_info_alter_table_array = []

            for info in fields:

                info_id_sql = prefix + info

                if info_id_sql in table_dest_struct:
                    log.debug(f"Field '{info_id_sql}' already exists in table")

                if (
                    info_id_sql in view_source_struct
                    and "INFO" in table_source_struct
                    and (info_id_sql not in table_dest_struct or force)
                ):

                    if "INFO" not in table_source_struct:
                        msg_err = f"Column 'INFO' not found in table '{table_source}' - Column 'INFO' needed!!!"
                        log.error(msg_err)
                        raise ValueError(msg_err)

                    if info_id_sql in table_dest_struct and force:
                        log.debug(
                            f"Explode INFO fields - Force '{info}' annotations fields update from 'INFO' column"
                        )

                    log.debug(f"Explode INFO fields - ADD '{info}' annotations fields")

                    # Get field type
                    type_sql = description_dict.get(info_id_sql, {})["type"]

                    # Add field
                    added_column = self.add_column(
                        table_name=table_dest,
                        column_name=info_id_sql,
                        column_type=type_sql,
                        default_value="null",
                        drop=force,
                    )

                    # Added column
                    if added_column:
                        added_columns.append(added_column)
                        log.debug(
                            f"Explode INFO fields - ADD '{info}' annotations fields - added"
                        )
                    else:
                        log.debug(
                            f"Explode INFO fields - ADD '{info}' annotations fields - not added"
                        )

                    # if added_column or force: #fileds_just_add
                    if (added_column or force) and not info in fields_just_add:

                        # add field to index
                        self.index_additionnal_fields.append(info_id_sql)

                        update_info_field = f"""
                            "{info_id_sql}" = {view_source}."{info_id_sql}"
                            """

                        # Set field append
                        sql_info_alter_table_array.append(update_info_field)

            if sql_info_alter_table_array:

                # Where clause join
                where_clause_join = f"""
                    {" AND ".join([f'"{table_dest}"."{key}" = "{view_source}"."{key}"' for key in table_key])}
                """

                # Evaluate block size
                batch_split = self.get_batch_split()

                # Insert by batch
                for batch_index in range(batch_split):

                    log.debug(
                        f"Explode INFO fields - Process batch [{batch_index+1}/{batch_split}]..."
                    )

                    where_clause = where_clause_join

                    # where clause
                    if batch_split > 1:
                        where_clause += (
                            f" AND ({table_dest}.POS % {batch_split}) = {batch_index} "
                        )
                    else:
                        where_clause += ""

                    # Update table
                    if proccess_all_fields_together:
                        sql_info_alter_table_array_join = ", ".join(
                            sql_info_alter_table_array
                        )
                        if sql_info_alter_table_array_join:
                            sql_info_alter_table = f"""
                                UPDATE {table_dest}
                                SET {sql_info_alter_table_array_join}
                                FROM {view_source}
                                WHERE {where_clause}
                                """
                            log.debug(
                                f"Explode INFO fields - Explode all {len(sql_info_alter_table_array)} fields..."
                            )
                            # log.debug(sql_info_alter_table)
                            self.get_connexion().execute(sql_info_alter_table)
                    else:
                        sql_info_alter_num = 0
                        for sql_info_alter in sql_info_alter_table_array:
                            sql_info_alter_num += 1
                            sql_info_alter_table = f"""
                                UPDATE {table_dest}
                                SET {sql_info_alter}
                                FROM {view_source}
                                WHERE {where_clause}
                                """
                            log.debug(
                                f"Explode INFO fields - Explode field {sql_info_alter_num}/{len(sql_info_alter_table_array)}..."
                            )
                            # log.debug(sql_info_alter_table)
                            self.get_connexion().execute(sql_info_alter_table)

            # Remove view_source
            self.remove_tables_or_views(tables=[view_source])

        # create indexes
        if create_index:
            self.create_indexes()

        return added_columns

    def create_indexes(self) -> None:
        """
        Create indexes on the table after insertion
        """

        # Access
        # access = self.get_config().get("access", None)
        access = self.get_access(default=None)

        # get table variants
        table_variants = self.get_table_variants("FROM")

        if self.get_indexing() and access not in ["RO"]:
            # Create index
            sql_create_table_index = f'CREATE INDEX IF NOT EXISTS idx_{self.get_table_variants()} ON {table_variants} ("#CHROM", "POS", "REF", "ALT")'
            self.get_connexion().execute(sql_create_table_index)
            sql_create_table_index = f'CREATE INDEX IF NOT EXISTS idx_{self.get_table_variants()}_chrom ON {table_variants} ("#CHROM")'
            self.get_connexion().execute(sql_create_table_index)
            sql_create_table_index = f'CREATE INDEX IF NOT EXISTS idx_{self.get_table_variants()}_pos ON {table_variants} ("POS")'
            self.get_connexion().execute(sql_create_table_index)
            sql_create_table_index = f'CREATE INDEX IF NOT EXISTS idx_{self.get_table_variants()}_ref ON {table_variants} ( "REF")'
            self.get_connexion().execute(sql_create_table_index)
            sql_create_table_index = f'CREATE INDEX IF NOT EXISTS idx_{self.get_table_variants()}_alt ON {table_variants} ("ALT")'
            self.get_connexion().execute(sql_create_table_index)
            for field in self.index_additionnal_fields:
                sql_create_table_index = f""" CREATE INDEX IF NOT EXISTS "idx_{self.get_table_variants()}_{field}" ON {table_variants} ("{field}") """
                self.get_connexion().execute(sql_create_table_index)

    def drop_indexes(self) -> None:
        """
        Create indexes on the table after insertion
        """

        # Access
        # access = self.get_config().get("access", None)
        access = self.get_access(default=None)

        # get table variants
        table_variants = self.get_table_variants("FROM")

        # Get database format
        connexion_format = self.get_connexion_format()

        if access not in ["RO"]:
            if connexion_format in ["duckdb"]:
                sql_list_indexes = f"SELECT index_name FROM duckdb_indexes WHERE table_name='{table_variants}'"
            elif connexion_format in ["sqlite"]:
                sql_list_indexes = f"SELECT name FROM sqlite_master WHERE type='index' AND tbl_name='{table_variants}';"

            list_indexes = self.get_connexion().execute(sql_list_indexes)
            index_names = [row[0] for row in list_indexes.fetchall()]
            for index in index_names:
                sql_drop_table_index = f""" DROP INDEX IF EXISTS "{index}" """
                self.get_connexion().execute(sql_drop_table_index)

    def read_vcf_header(self, f) -> list:
        """
        It reads the header of a VCF file and returns a list of the header lines

        :param f: the file object
        :return: The header lines of the VCF file.
        """

        header_list = []
        for line in f:
            header_list.append(line)
            if line.startswith("#CHROM"):
                break
        return header_list

    def read_vcf_header_file(self, file: str = None) -> list:
        """
        The `read_vcf_header_file` function reads the header of a VCF file, handling both compressed and
        uncompressed files.

        :param file: The `file` parameter is a string that represents the path to the VCF header file
        that you want to read. It is an optional parameter, so if you don't provide a value, it will
        default to `None`
        :type file: str
        :return: The function `read_vcf_header_file` returns a list.
        """

        if self.get_input_compressed(input_file=file):
            with bgzf.open(file, "rt") as f:
                return self.read_vcf_header(f=f)
        else:
            with open(file, "rt") as f:
                return self.read_vcf_header(f=f)

    def execute_query(self, query: str):
        """
        It takes a query as an argument, executes it, and returns the results

        :param query: The query to be executed
        :return: The result of the query is being returned.
        """
        if query:
            return self.get_connexion().execute(query)  # .fetchall()
        else:
            return None

    def export_output(
        self,
        output_file: str | None = None,
        output_header: str | None = None,
        export_header: bool = True,
        explode_infos: bool = True,
        explode_infos_prefix: list = None,
        explode_infos_fields: list = None,
        header_in_output: bool = None,
        query: str | None = None,
        parquet_partitions: list | None = None,
        chunk_size: int | None = None,
        threads: int | None = None,
        sort: bool = False,
        index: bool = False,
        order_by: str | None = None,
        fields_to_rename: dict | None = None,
        force_cast_as_flat: bool = False,
        **kwargs
    ) -> bool:
        """
        The `export_output` function exports data from a VCF file to various formats, including VCF,
        CSV, TSV, PSV, and Parquet, with options for customization such as filtering, sorting, and
        partitioning.

        :param output_file: The `output_file` parameter is a string that specifies the name of the
        output file where the exported data will be saved
        :type output_file: str | None
        :param output_header: The `output_header` parameter is a string that specifies the name of the
        file where the header of the VCF file will be exported. If this parameter is not provided, the
        header will be exported to a file with the same name as the `output_file` parameter, but with
        the extension "
        :type output_header: str | None
        :param export_header: The `export_header` parameter is a boolean flag that determines whether
        the header of a VCF file should be exported to a separate file or not. If `export_header` is
        True, the header will be exported to a file. If `export_header` is False, the header will not
        be, defaults to True
        :type export_header: bool (optional)
        :param explode_infos: The `explode_infos` parameter is a boolean flag that determines whether
        the INFO fields in the VCF file should be exploded into individual columns in the output file.
        If `explode_infos` is set to True, the INFO fields will be exploded. If `explode_infos` is set
        to False, the INFO fields will not be exploded. By default, the INFO fields are exploded
        :type explode_infos: bool (optional)
        :param explode_infos_prefix: The `explode_infos_prefix` parameter is a string that specifies a
        prefix to be added to the names of the exploded INFO fields in the output file. This allows for better organization and identification of the INFO fields in the exported data. If not provided, a default prefix will be used
        :type explode_infos_prefix: str | None
        :param explode_infos_fields: The `explode_infos_fields` parameter is a list of specific INFO fields that you want to explode into individual columns in the output file. If this parameter is provided, only the specified INFO fields will be exploded. If it is not provided or is set to None
        , all INFO fields will be exploded by default. This allows for selective extraction of relevant information from the VCF file during the export process
        :type explode_infos_fields: list | None
        :param query: The `query` parameter in the `export_output` function is an optional SQL query
        that can be used to filter and select specific data from the VCF file before exporting it. If
        provided, only the data that matches the query will be exported. This allows you to customize
        the exported data based on
        :type query: str | None
        :param header_in_output: The `header_in_output` parameter is a boolean flag that determines
        whether the header should be included in the output file. If `header_in_output` is set to `True`,
        the header will be included in the output file. If `header_in_output` is set to `False`, the
        header will not be included in the output file. By default, the header is included in the output
        file
        :type header_in_output: bool (optional)
        :param parquet_partitions: The `parquet_partitions` parameter is a list that specifies the
        columns to be used for partitioning the Parquet file during export. Partitioning is a way to
        organize data in a hierarchical directory structure based on the values of one or more columns.
        This can improve query performance when working with large datasets
        :type parquet_partitions: list | None
        :param chunk_size: The `chunk_size` parameter specifies the number of records in a batch when
        exporting data in Parquet format. This parameter is used for partitioning the Parquet file into
        multiple files. It helps in optimizing the export process by breaking down the data into
        manageable chunks for processing and storage
        :type chunk_size: int | None
        :param threads: The `threads` parameter in the `export_output` function specifies the number of
        threads to be used during the export process. It determines the level of parallelism and can
        improve the performance of the export operation. If this parameter is not provided, the function
        will use the default number of threads
        :type threads: int | None
        :param sort: The `sort` parameter in the `export_output` function is a boolean flag that
        determines whether the output file should be sorted based on genomic coordinates of the
        variants. If `sort` is set to `True`, the output file will be sorted. If `sort` is set to
        `False`,, defaults to False
        :type sort: bool (optional)
        :param index: The `index` parameter in the `export_output` function is a boolean flag that
        determines whether an index should be created on the output file. If `index` is set to `True`,
        an index will be created on the output file. If `index` is set to `False`, no, defaults to False
        :type index: bool (optional)
        :param order_by: The `order_by` parameter in the `export_output` function is a string that
        specifies the column(s) to use for sorting the output file. This parameter is only applicable
        when exporting data in VCF format. It allows you to specify the column(s) based on which the
        output file should be
        :type order_by: str | None
        :param fields_to_rename: The `fields_to_rename` parameter is a dictionary that specifies the
        mapping of field names to be renamed during the export process. This parameter allows you to
        customize the output field names before exporting the data. Each key-value pair in the
        dictionary represents the original field name as the key and the new field name
        :type fields_to_rename: dict | None
        :param force_cast_as_flat: Only for Parquet format. The `force_cast_as_flat` parameter is a boolean
        flag that determines whether to force the export of nested or complex data structures as flat
        structures. If `force_cast_as_flat` is set to `True`, the function will flatten any nested
        structures in the data before exporting it. If `force_cast_as_flat` is set to `False`, the
        function will preserve the original structure of the data during export. By default, it is set
        to False
        :type force_cast_as_flat: bool (optional)
        :return: The `export_output` function returns a boolean value. It checks if the output file
        exists and returns True if it does, or None if it doesn't.
        """

        # Log
        log.info("Exporting...")

        # Full path
        output_file = full_path(output_file)
        output_header = full_path(output_header)

        # Config
        config = self.get_config()

        # Param
        param = self.get_param()

        # Tmp files to remove
        tmp_to_remove = []

        # If no output, get it
        if not output_file:
            output_file = self.get_output()

        # If not threads
        if not threads:
            threads = self.get_threads()

        # Rename fields
        if not fields_to_rename:
            fields_to_rename = param.get("export", {}).get("fields_to_rename", None)
        self.rename_info_fields(fields_to_rename=fields_to_rename)

        # Force cast as flat
        force_cast_as_flat = param.get("export", {}).get(
            "force_cast_as_flat", force_cast_as_flat
        )

        # Auto header name with extension
        if export_header or output_header:
            if not output_header:
                output_header = f"{output_file}.hdr"
            # Export header
            self.export_header(output_file=output_file, query=query)

        # Switch off export header if VCF output
        output_file_type = get_file_format(output_file)
        if output_file_type in ["vcf"]:
            export_header = False
            tmp_to_remove.append(output_header)

        # Chunk size
        if not chunk_size:
            chunk_size = config.get("chunk_size", None)

        # Parquet partition
        if not parquet_partitions:
            parquet_partitions = param.get("export", {}).get("parquet_partitions", None)
        if parquet_partitions and isinstance(parquet_partitions, str):
            parquet_partitions = parquet_partitions.split(",")

        # Order by
        if not order_by:
            order_by = param.get("export", {}).get("order_by", "")

        # Header in output
        if header_in_output is None:
            header_in_output = param.get("export", {}).get("include_header", False)

        # Database
        database_source = self.get_connexion()

        # Connexion format
        connexion_format = self.get_connexion_format()

        # Explode infos
        if self.get_explode_infos() and explode_infos:
            self.explode_infos(
                prefix=explode_infos_prefix or self.get_explode_infos_prefix(),
                fields=explode_infos_fields or self.get_explode_infos_fields(),
                force=False,
                fields_forced_as_varchar=True,
            )

        # if connexion_format in ["sqlite"] or query:
        if connexion_format in ["sqlite"]:

            # Export in Parquet
            random_tmp = "".join(
                random.choice(string.ascii_lowercase) for i in range(10)
            )
            database_source = f"""{output_file}.{random_tmp}.database_export.parquet"""
            tmp_to_remove.append(database_source)

            # Table Variants
            table_variants = self.get_table_variants()

            # Create export query
            sql_query_export_subquery = f"""
                SELECT * FROM {table_variants}
                """

            # Write source file
            fp.write(database_source, self.get_query_to_df(sql_query_export_subquery))

        # Create database
        database = Database(
            database=database_source,
            table="variants",
            header_file=output_header,
            conn_config=self.get_connexion_config(),
        )

        # Existing colomns header
        existing_columns_header = database.get_header_columns_from_database(query=query)

        # Sample list
        if output_file_type in ["vcf"]:
            get_samples = self.get_samples()
            get_samples_check = self.get_samples_check()
            samples_force = get_samples is not None
            sample_list = self.get_header_sample_list(
                check=get_samples_check,
                samples=get_samples,
                samples_force=samples_force,
            )
        else:
            sample_list = None

        # Export file
        database.export(
            output_database=output_file,
            output_header=output_header,
            existing_columns_header=existing_columns_header,
            parquet_partitions=parquet_partitions,
            chunk_size=chunk_size,
            threads=threads,
            sort=sort,
            index=index,
            header_in_output=header_in_output,
            order_by=order_by,
            query=query,
            export_header=export_header,
            sample_list=sample_list,
            force_cast_as_flat=force_cast_as_flat,
        )

        # Remove
        remove_if_exists(tmp_to_remove)

        return (os.path.exists(output_file) or None) and (
            os.path.exists(output_file) or None
        )

    def get_columns(self, table: str = None) -> list:
        """
        The `get_columns` function returns a list of columns in a specified table. If the `table`
        parameter is not provided when calling the function, it will default to using the variants table.

        Args:
            table (str, optional): The name of the table from which you want to retrieve the columns. If not provided,
                it will default to using the variants table.

        Returns:
            list: A list of columns in the specified table.
        """

        if not table:
            table = self.get_table_variants()

        # Use PRAGMA table_info for SQLite or DESCRIBE for other databases
        connexion_format = self.get_connexion_format()
        if connexion_format == "sqlite":
            query = f"PRAGMA table_info({table})"
            columns_info = self.get_query_to_df(query)
            columns = columns_info["name"].tolist()
        else:
            query = f"DESCRIBE {table}"
            columns_info = self.get_query_to_df(query)
            columns = columns_info["column_name"].tolist()

        return columns

    def get_extra_infos(self, table: str = None) -> list:
        """
        The `get_extra_infos` function returns a list of columns that are in a specified table but not
        in the header.

        :param table: The `table` parameter in the `get_extra_infos` function is used to specify the
        name of the table from which you want to retrieve the extra columns that are not present in the
        header. If the `table` parameter is not provided when calling the function, it will default to
        using the variants
        :type table: str
        :return: A list of columns that are in the specified table but not in the header of the table.
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

        return ", ".join(
            ['"' + str(elem) + '"' for elem in self.get_extra_infos(table=table)]
        )

    def export_header(
        self,
        header_name: str = None,
        output_file: str = None,
        output_file_ext: str = ".hdr",
        clean_header: bool = True,
        clean_info_flag: bool = False,
        remove_chrom_line: bool = False,
        query: str | None = None,
    ) -> str:
        """
        The `export_header` function takes a VCF file, extracts the header, modifies it according to
        specified options, and writes it to a new file.

        :param header_name: The `header_name` parameter is the name of the header file to be created. If
        this parameter is not specified, the header will be written to the output file
        :type header_name: str
        :param output_file: The `output_file` parameter in the `export_header` function is used to
        specify the name of the output file where the header will be written. If this parameter is not
        provided, the header will be written to a temporary file
        :type output_file: str
        :param output_file_ext: The `output_file_ext` parameter in the `export_header` function is a
        string that represents the extension of the output header file. By default, it is set to ".hdr"
        if not specified by the user. This extension will be appended to the `output_file` name to
        create the final, defaults to .hdr
        :type output_file_ext: str (optional)
        :param clean_header: The `clean_header` parameter in the `export_header` function is a boolean
        flag that determines whether the header should be cleaned or not. When `clean_header` is set to
        `True`, the function will clean the header by modifying certain lines based on a specific
        pattern. If `clean_header`, defaults to True
        :type clean_header: bool (optional)
        :param clean_info_flag: The `clean_info_flag` parameter in the `export_header` function is a boolean
        flag that determines whether the header should be cleaned for INFO/tags that are 'Flag' type.
        When `clean_info_flag` is set to `True`, the function will replace INFO/tags 'Type' as 'String'.
        Default to False
        :type clean_info_flag: bool (optional)
        :param remove_chrom_line: The `remove_chrom_line` parameter in the `export_header` function is a
        boolean flag that determines whether the #CHROM line should be removed from the header before
        writing it to the output file. If set to `True`, the #CHROM line will be removed; if set to `,
        defaults to False
        :type remove_chrom_line: bool (optional)
        :param query: The `query` parameter in the `export_header` function is an optional SQL query
        string that can be used to filter the columns in the header. If provided, the function will
        retrieve only the columns that match the query. If not provided, all columns in the header will
        be included in the output
        :type query: str | None
        :return: The function `export_header` returns the name of the temporary header file that is
        created.
        """

        if not header_name and not output_file:
            output_file = self.get_output()

        if self.get_header():

            # Get header object
            header_obj = self.get_header()

            # Create database
            db_for_header = Database(database=self.get_input())

            # Get real columns in the file
            db_header_columns = db_for_header.get_columns(sql_query=query)

            with tempfile.TemporaryDirectory() as tmpdir:

                # Write header file
                header_file_tmp = os.path.join(tmpdir, "header")
                f = open(header_file_tmp, "w")
                vcf.Writer(f, header_obj)
                f.close()

                # Replace #CHROM line with rel columns
                header_list = db_for_header.read_header_file(
                    header_file=header_file_tmp
                )
                header_list[-1] = "\t".join(db_header_columns)

                # Remove CHROM line
                if remove_chrom_line:
                    header_list.pop()

                # Clean header
                if clean_header:
                    header_list_clean = []
                    for head in header_list:
                        # Clean head for malformed header
                        head_clean = head
                        if clean_info_flag:
                            head_clean = re.subn(
                                "##FORMAT=<ID=(.*),Number=(.*),Type=Flag",
                                r"##FORMAT=<ID=\1,Number=\2,Type=String",
                                head_clean,
                                2,
                            )[0]
                        # Write header
                        header_list_clean.append(head_clean)
                    header_list = header_list_clean

            tmp_header_name = output_file + output_file_ext

            f = open(tmp_header_name, "w")
            for line in header_list:
                f.write(line)
            f.close()

        return tmp_header_name

    def export_variant_vcf(
        self,
        vcf_file,
        remove_info: bool = False,
        add_samples: bool = True,
        list_samples: list = [],
        where_clause: str = "",
        chrom_mapping_sql: str = None,
        sort: bool = True,
        index: bool = False,
        threads: int | None = None,
    ) -> bool | None:
        """
        The `export_variant_vcf` function exports a VCF file with specified samples, allowing options to
        remove INFO field, add samples, and control compression and indexing.

        :param vcf_file: The `vcf_file` parameter is the name of the file where the VCF data will be
        written to. It is the output file that will contain the filtered VCF data based on the specified
        parameters
        :param remove_info: The `remove_info` parameter in the `export_variant_vcf` function is a
        boolean flag that determines whether to remove the INFO field from the output VCF file. If set
        to `True`, the INFO field will be removed. If set to `False`, the INFO field will be included
        in, defaults to False
        :type remove_info: bool (optional)
        :param add_samples: The `add_samples` parameter is a boolean parameter that determines whether
        the samples should be added to the VCF file or not. If set to True, the samples will be added.
        If set to False, the samples will be removed. The default value is True, defaults to True
        :type add_samples: bool (optional)
        :param list_samples: The `list_samples` parameter is a list of samples that you want to include
        in the output VCF file. By default, all samples will be included. If you provide a list of
        samples, only those samples will be included in the output file
        :type list_samples: list
        :param where_clause: The `where_clause` parameter in the `export_variant_vcf` function is a
        string that represents a SQL WHERE clause. It is used to filter the variants that will be
        exported to the VCF file. The `where_clause` allows you to specify conditions that the variants
        must meet in order to be included in the output VCF file. If no `where_clause` is provided, all
        variants will be exported
        :type where_clause: str
        :param chrom_mapping_sql: The `chrom_mapping_sql` parameter in the `export_variant_vcf` function is a SQL string that is used to map chromosome names in the `#CHROM` column of the VCF file. If no `chrom_mapping_sql` is provided, the `#CHROM` column will remain unchanged. This allows for customization of chromosome naming conventions in the output VCF file
        :type chrom_mapping_sql: str | None
        :param sort: The `sort` parameter in the `export_variant_vcf` function is a boolean flag that
        determines whether the output VCF file should be sorted based on genomic coordinates of the variants.
        If `sort` is set to `True`, the output VCF file will be sorted. If `sort` is set to `False`, the output VCF file
        will not be sorted. By default, the output VCF file is sorted
        :type sort: bool (optional)
        :param index: The `index` parameter in the `export_variant_vcf` function is a boolean flag that
        determines whether or not to create an index for the output VCF file. If `index` is set to
        `True`, the output VCF file will be indexed using tabix. If `index`, defaults to False
        :type index: bool (optional)
        :param threads: The `threads` parameter in the `export_variant_vcf` function specifies the
        number of threads to use for exporting the VCF file. It determines how many parallel threads
        will be used during the export process. More threads can potentially speed up the export process
        by utilizing multiple cores of the processor. If
        :type threads: int | None
        :return: The `export_variant_vcf` function returns the result of calling the `export_output`
        method with various parameters including the output file, query, threads, sort flag, and index
        flag. The `export_output` method is responsible for exporting the VCF data based on the
        specified parameters and configurations provided in the `export_variant_vcf` function.
        """

        # Config
        config = self.get_config()

        # Extract VCF
        log.debug("Export VCF...")

        # Table variants
        table_variants = self.get_table_variants()

        # Threads
        if not threads:
            threads = self.get_threads()

        # Info fields
        if remove_info:
            if not isinstance(remove_info, str):
                remove_info = "."
            info_field = f"""'{remove_info}' as INFO"""
        else:
            info_field = "INFO"

        # Samples fields
        if add_samples:
            if not list_samples:
                list_samples = self.get_header_sample_list()
            if list_samples:
                samples_fields = " , FORMAT , " + " , ".join(
                    [f""" "{sample}" """ for sample in list_samples]
                )
            else:
                samples_fields = ""
            log.debug(f"samples_fields: {samples_fields}")
        else:
            samples_fields = ""

        # Where clause
        if where_clause is None:
            where_clause = ""

        # Columns
        existing_columns = self.get_columns(table=table_variants)
        columns_default_values = {
            "#CHROM": "'chr'",
            "POS": "0",
            "ID": "'.'",
            "REF": "'N'",
            "ALT": "'N'",
            "QUAL": "'0'",
            "FILTER": "'PASS'",
        }

        select_fields_list = []
        for column in ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER"]:
            if column not in existing_columns:
                select_fields_list.append(
                    f"{columns_default_values.get(column, '')} AS '{column}'"
                )
            elif column == "#CHROM" and chrom_mapping_sql:
                chrom_expr = chrom_mapping_sql
                select_fields_list.append(f'{chrom_expr} AS "#CHROM"')
            else:
                select_fields_list.append(f'"{column}"')
        select_fields = ", ".join(select_fields_list)

        # Query
        sql_query_select = f""" SELECT {select_fields}, {info_field} {samples_fields} FROM {table_variants} {where_clause} """
        # log.debug(f"sql_query_select={sql_query_select}")

        return self.export_output(
            output_file=vcf_file,
            output_header=None,
            export_header=True,
            explode_infos=False,
            query=sql_query_select,
            parquet_partitions=None,
            chunk_size=config.get("chunk_size", None),
            threads=threads,
            sort=sort,
            index=index,
            order_by=None,
        )

    def run_commands(self, commands: list = [], threads: int = 1) -> None:
        """
        It takes a list of commands and runs them in parallel using the number of threads specified

        :param commands: A list of commands to run
        :param threads: The number of threads to use, defaults to 1 (optional)
        """

        run_parallel_commands(commands, threads)

    def get_threads(self, default: int = 1) -> int:
        """
        This function returns the number of threads to use for a job, with a default value of 1 if not
        specified.

        :param default: The `default` parameter in the `get_threads` method is used to specify the
        default number of threads to use if no specific value is provided. If no value is provided for
        the `threads` parameter in the configuration or input parameters, the `default` value will be
        used, defaults to 1
        :type default: int (optional)
        :return: the number of threads to use for the current job.
        """

        # Config
        config = self.get_config()

        # Param
        param = self.get_param()

        # Input threads
        input_thread = param.get("threads", config.get("threads", None))

        # Check threads
        if not input_thread:
            threads = default
        elif int(input_thread) <= 0:
            threads = os.cpu_count()
        else:
            threads = int(input_thread)
        return threads

    def get_memory_system(self, type: str = "available", unit: str = "G") -> str:
        """
        This function retrieves the system memory in the system and returns it as a string.

        :param type: The `type` parameter in the `get_memory_system` function specifies the type of
        memory information to retrieve. It can take one of the following values: "total" to retrieve the
        total memory, "used" to retrieve the used memory, "percent" to retrieve the percentage of used memory,
        or "available" to retrieve the available memory. The default value is "available"
        :type type: str (optional)
        :param unit: The `unit` parameter in the `get_memory_system` function specifies the unit of
        measurement for the memory value. It can take one of the following values: "K" for kilobytes,
        "M" for megabytes, "G" for gigabytes, or "T" for terabytes. The default value is "G"
        :type unit: str (optional)

        :return: The function `get_memory_system` returns a string representation of the available
        memory in the system.
        """

        import psutil  # type: ignore

        # Check system memory
        mem = psutil.virtual_memory()

        # Get memory type
        if type == "total":
            memory = str(mem.total)
        elif type == "used":
            memory = str(mem.used)
        elif type == "percent":
            memory = str(mem.percent)
        elif type == "available":
            memory = str(mem.available)
        else:
            memory = str(mem.total)

        # Convert unit
        unit_powers = {"K": 1, "M": 2, "G": 3, "T": 4}
        if unit in unit_powers:
            power = unit_powers[unit]
            memory = str(int(memory) // (1024**power))

        # Return memory
        return f"{memory}{unit}"

    def get_memory(self, default: str = None, available: bool = False) -> str:
        """
        This function retrieves the memory value from parameters or configuration with a default value
        if not found.

        :param default: The `get_memory` function takes in a default value as a string parameter. This
        default value is used as a fallback in case the `memory` parameter is not provided in the
        `param` dictionary or the `config` dictionary. If `memory` is not found in either dictionary,
        the function
        :type default: str
        :param available: A boolean parameter that determines whether to limit the memory to the
        available system memory or not. If set to True, the function will return the minimum value
        between the input memory and the available system memory. If set to False, the function will
        return the input memory or the default value, defaults to False
        :type available: bool (optional)
        :return: The `get_memory` function returns a string value representing the memory parameter. If
        the `input_memory` is provided in the parameters, it will return that value. Otherwise, it will
        return the default value provided as an argument to the function.
        """

        # Config
        config = self.get_config()

        # Param
        param = self.get_param()

        # Input threads
        input_memory = param.get("memory", config.get("memory", None))

        # Avalable memory
        if available:
            available_memory = self.get_memory_system(type="available", unit="G")
            if available_memory is not None and input_memory is not None:
                input_memory = min(input_memory, available_memory)
            else:
                input_memory = available_memory

        # Check threads
        if input_memory:
            memory = input_memory
        else:
            memory = default

        return memory

    def update_from_vcf(
        self,
        vcf_file: str,
        update_existing_fields: bool = False,
        remove_vcf_file: bool = True,
        upper_case: bool = True,
        update_header: bool = False,
        annotation_header_fields_override: dict | None = None,
        chrom_mapping_sql: str = None
    ) -> None:
        """
        > If the database is duckdb, then use the parquet method, otherwise use the sqlite method

        :param vcf_file: the path to the VCF file you want to update the database with
        :type vcf_file: str
        :param update_existing_fields: If True, existing fields in the INFO column will be updated
        with the values from the VCF file. If False, only new fields will be added, defaults to False
        :type update_existing_fields: bool (optional)
        :param remove_vcf_file: If True, the VCF file will be removed after the update is complete,
        defaults to True
        :type remove_vcf_file: bool (optional)
        :param upper_case: If True, the ALT and REF fields will be compared in uppercase, defaults to True
        :type upper_case: bool (optional)
        :param update_header: If True, the header of the VCF file will be updated
        :type update_header: bool (optional)
        :param annotation_header_fields_override: A dictionary that allows you to override specific
        fields in the annotation header when updating the database from a VCF file. The keys of
        the dictionary represent the field names in the annotation header, and the values represent
        the new values that you want to assign to those fields. This parameter is optional and can
        be used to customize the annotation header during the update process. If not provided, the
        default values from the VCF file will be used for the annotation header fields
        :type annotation_header_fields_override: dict | None (optional)
        :param chrom_mapping_sql: SQL string used to map chromosome names in the `#CHROM` column of the VCF file. If not provided, the `#CHROM` column will remain unchanged.
        :type chrom_mapping_sql: str | None (optional)

        :return: None
        """

        # Header

        if update_header:

            # VCF header
            vcf_reader = self.get_header()

            # Find annotation in header
            vcf_file_header = self.get_header(type="vcf", vcf_file=vcf_file)

            # Add annotation to header if not exist
            for ann in vcf_file_header.infos:
                #log.debug(f"Check annotation '{ann}' in header...")
                if ann not in self.get_header().infos or update_existing_fields:
                    ann_info = vcf_file_header.infos.get(ann)
                    vcf_reader.infos[ann] = self.build_info_with_header_override(
                        field_name=ann,
                        field_number=ann_info.num,
                        field_type=ann_info.type,
                        field_description=ann_info.desc,
                        field_source=ann_info.source,
                        field_version=ann_info.version,
                        header_fields_override=annotation_header_fields_override,
                    )

        # Update content of vcf file in database

        connexion_format = self.get_connexion_format()

        if connexion_format in ["duckdb"]:
            self.update_from_vcf_duckdb(
                vcf_file,
                update_existing_fields=update_existing_fields,
                remove_vcf_file=remove_vcf_file,
                upper_case=upper_case,
                chrom_mapping_sql=chrom_mapping_sql
            )
        elif connexion_format in ["sqlite"]:
            self.update_from_vcf_sqlite(vcf_file)

        if remove_vcf_file:
            remove_if_exists([vcf_file])

    def update_from_vcf_duckdb(
        self,
        vcf_file: str,
        update_existing_fields: bool = False,
        remove_vcf_file: bool = True,
        upper_case: bool = True,
        chrom_mapping_sql: str | None = None,
    ) -> None:
        """
        It takes a VCF file and updates the INFO column of the variants table in the database with the
        INFO column of the VCF file

        :param vcf_file: The path to the VCF file you want to update the database with
        :type vcf_file: str
        :param update_existing_fields: If True, existing fields in the INFO column will be updated
        with the values from the VCF file. If False, only new fields will be added, defaults to False
        :type update_existing_fields: bool (optional)
        :param remove_vcf_file: If True, the VCF file will be removed after the update is complete,
        defaults to True
        :type remove_vcf_file: bool (optional)
        :param upper_case: If True, the ALT and REF fields will be compared in uppercase, defaults to True
        :type upper_case: bool (optional)
        :param chrom_mapping_sql: SQL string used to map chromosome names in the `#CHROM` column of the VCF file. If not provided, the `#CHROM` column will remain unchanged.
        :type chrom_mapping_sql: str | None (optional)

        :return: None
        """

        # variants table
        table_variants = self.get_table_variants()

        # Connexion
        conn = self.get_connexion()

        log.info(f"Update variants table from file '{os.path.basename(vcf_file)}'...")

        with TemporaryDirectory(dir=self.get_tmp_dir()) as tmp_dir:

            log.debug(f"Create parquet files from VCF '{vcf_file}'...")

            # Create parquet from VCF
            vcf_file_parquet_path = os.path.join(tmp_dir, "vcf_file.parquet")
            vcf_file_parquet = Variants(
                input=vcf_file, load=True, config={"access": "RO"}
            )

            log.debug(f"Variants input format '{vcf_file_parquet.get_input_format()}'")

            if vcf_file_parquet.get_input_format() == "parquet":

                # list of parquet files
                vcf_file_parquet_path = vcf_file

            else:

                # Export parquet parameters
                chunk_size = self.get_config().get("chunk_size", None)
                threads = self.get_threads()

                # Export parquet files
                log.debug("Export VCF to partitioned parquet...")
                vcf_file_parquet.export_output(
                    output_file=vcf_file_parquet_path,
                    chunk_size=chunk_size,
                    threads=threads,
                    export_header=True,
                )
                log.debug(f"Parquet generated: {vcf_file_parquet_path}")

                if remove_vcf_file:
                    remove_if_exists([vcf_file])

            # Update if fields exist
            if update_existing_fields:
                # list of header columns
                header_columns = self.get_header().infos.keys()
                header_columns_vcf_file_parquet = (
                    vcf_file_parquet.get_header().infos.keys()
                )

                # columns that exist in both
                common_columns = list(
                    set(header_columns).intersection(
                        set(header_columns_vcf_file_parquet)
                    )
                )

                # Remove common columns
                if len(common_columns) > 0:
                    log.debug(f"Common columns to update/remove: {common_columns}")
                    self.rename_info_fields(
                        fields_to_rename=dict.fromkeys(common_columns, None)
                    )
                else:
                    log.debug("No common columns to update/remove")

            # Upper case function for ALT and REF
            if upper_case:
                upper_func = "upper"
            else:
                upper_func = ""

            # Chromosome mapping (from_tool -> internal)
            if chrom_mapping_sql:
                chrom_expr = chrom_mapping_sql
                chrom_select = f'{chrom_expr} AS "#CHROM"'
            else:
                chrom_select = '"#CHROM"'

            # Create table/view from parquet files
            table_source_name = "table_parquet_" + get_random(10)
            sql_query_update = f"""
                CREATE VIEW {table_source_name}
                AS (
                    SELECT {chrom_select}, POS, {upper_func}(REF) as REF, {upper_func}(ALT) as ALT, INFO
                    FROM read_parquet('{vcf_file_parquet_path}')
                    WHERE INFO NOT IN ('','.')
                )
                ;
            """
            # log.debug(f"sql_query_update: {sql_query_update}")
            conn.execute(sql_query_update)

            # Update INFO fields with update_table function
            source = {
                "table": table_source_name,
                "join_keys": ["#CHROM", "POS", "REF", "ALT"],
                "columns": {
                    "INFO": {
                        "columns": ["INFO"],
                        "mode": "append",
                        "separator": ";",
                    }
                },
            }
            self.update_table(
                dest_table=table_variants,
                sources=[source],
                physical_order=True,
                force_strategy=None,
                upper_case=upper_case,
            )

            return None

    def update_from_vcf_sqlite(self, vcf_file: str) -> None:
        """
        It creates a temporary table in the SQLite database, loads the VCF file into the temporary
        table, then updates the INFO column of the variants table with the INFO column of the temporary
        table

        :param vcf_file: The path to the VCF file you want to update the database with
        """

        # Create a temporary table for the VCF
        table_vcf = "tmp_vcf"
        sql_create = (
            f"CREATE TEMPORARY TABLE {table_vcf} AS SELECT * FROM variants WHERE 0"
        )
        self.get_connexion().execute(sql_create)

        # Loading VCF into temporaire table
        vcf_df = pd.read_csv(
            vcf_file, sep="\t", comment="#", header=None, low_memory=False
        )
        vcf_df.columns = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"]
        vcf_df.to_sql(table_vcf, self.conn, if_exists="append", index=False)

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
        self.get_connexion().execute(sql_query_update)

        # Drop temporary table
        sql_drop = f"DROP TABLE {table_vcf}"
        self.get_connexion().execute(sql_drop)

    def drop_variants_table(self) -> None:
        """
        > This function drops the variants table
        """

        table_variants = self.get_table_variants()
        sql_table_variants = f"DROP TABLE IF EXISTS {table_variants}"
        self.get_connexion().execute(sql_table_variants)

    def set_variant_id(
        self, variant_id_column: str = "variant_id", force: bool = None
    ) -> str:
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
        assembly = self.get_param().get(
            "assembly", self.get_config().get("assembly", DEFAULT_ASSEMBLY)
        )

        # INFO/Tag prefix
        prefix = self.get_explode_infos_prefix()

        # Explode INFO/SVTYPE
        added_columns = self.explode_infos(prefix=prefix, fields=["SVTYPE"])

        # variants table
        table_variants = self.get_table_variants()

        # variant_id column
        if not variant_id_column:
            variant_id_column = "variant_id"

        # Creta variant_id column
        if "variant_id" not in self.get_extra_infos() or force:

            # Create column
            self.add_column(
                table_name=table_variants,
                column_name=variant_id_column,
                column_type="UBIGINT",
                default_value="0",
            )

            # Update column
            self.get_connexion().execute(f"""
                    UPDATE {table_variants}
                    SET "{variant_id_column}" = hash('{assembly}', "#CHROM", "POS", "REF", "ALT", '"{prefix}SVTYPE"')
                """)

        # Remove added columns
        for added_column in added_columns:
            self.drop_column(column=added_column)

        # return variant_id column name
        return variant_id_column

    def get_variant_id_column(
        self, variant_id_column: str = "variant_id", force: bool = None
    ) -> str:
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

    def get_columns_type(self, table: str) -> dict:
        """
        Get columns type of a table.
        :param table: The `table` parameter is a string that represents the name of the table
        for which you want to retrieve the column types.
        :type table: str
        :return: The function `get_columns_type` returns a dictionary where the keys are the
        column names of the specified table and the values are the corresponding data types of
        those columns.
        """

        # Get columns info
        query = f"""
            PRAGMA table_info('{table}')
        """
        df_columns = self.get_query_to_df(query)

        # Construct columns type dict
        columns_type = {}
        for _, row in df_columns.iterrows():
            columns_type[row["name"]] = row["type"]

        return columns_type

    ##################
    # Prioritization #
    ##################

    def get_config_default(self, name: str) -> dict:
        """
        The function `get_config_default` returns a dictionary containing default configurations for
        various calculations and prioritizations.

        :param name: The `get_config_default` function returns a dictionary containing default
        configurations for different calculations and prioritizations. The `name` parameter is used to
        specify which specific configuration to retrieve from the dictionary
        :type name: str
        :return: The function `get_config_default` returns a dictionary containing default configuration
        settings for different calculations and prioritizations. The specific configuration settings are
        retrieved based on the input `name` parameter provided to the function. If the `name` parameter
        matches a key in the `config_default` dictionary, the corresponding configuration settings are
        returned. If there is no match, an empty dictionary is returned.
        """

        config_default = {
            "calculations": self.get_config_calculations_default(),
            "prioritizations": self.get_config_prioritizations_default(),
            #     {
            #     "default": {
            #         "ANN": [
            #             {
            #                 "type": "contains",
            #                 "value": "HIGH",
            #                 "score": 5,
            #                 "flag": "PASS",
            #                 "comment": [
            #                     "The variant is assumed to have high (disruptive) impact in the protein, probably causing protein truncation, loss of function or triggering nonsense mediated decay"
            #                 ],
            #             },
            #             {
            #                 "type": "contains",
            #                 "value": "MODERATE",
            #                 "score": 3,
            #                 "flag": "PASS",
            #                 "comment": [
            #                     "A non-disruptive variant that might change protein effectiveness"
            #                 ],
            #             },
            #             {
            #                 "type": "contains",
            #                 "value": "LOW",
            #                 "score": 0,
            #                 "flag": "FILTERED",
            #                 "comment": [
            #                     "Assumed to be mostly harmless or unlikely to change protein behavior"
            #                 ],
            #             },
            #             {
            #                 "type": "contains",
            #                 "value": "MODIFIER",
            #                 "score": 0,
            #                 "flag": "FILTERED",
            #                 "comment": [
            #                     "Usually non-coding variants or variants affecting non-coding genes, where predictions are difficult or there is no evidence of impact"
            #                 ],
            #             },
            #         ],
            #     }
            # },
        }

        return config_default.get(name, None)

    def get_config_json(
        self, name: str, config_dict: dict = {}, config_file: str = None
    ) -> dict:
        """
        The function `get_config_json` retrieves a configuration JSON object with prioritizations from
        default values, a dictionary, and a file.

        :param name: The `name` parameter in the `get_config_json` function is a string that represents
        the name of the configuration. It is used to identify and retrieve the configuration settings
        for a specific component or module
        :type name: str
        :param config_dict: The `config_dict` parameter in the `get_config_json` function is a
        dictionary that allows you to provide additional configuration settings or overrides. When you
        call the `get_config_json` function, you can pass a dictionary containing key-value pairs where
        the key is the configuration setting you want to override or
        :type config_dict: dict
        :param config_file: The `config_file` parameter in the `get_config_json` function is used to
        specify the path to a configuration file that contains additional settings. If provided, the
        function will read the contents of this file and update the configuration dictionary with the
        values found in the file, overriding any existing values with the
        :type config_file: str
        :return: The function `get_config_json` returns a dictionary containing the configuration
        settings.
        """

        # Create with default prioritizations
        config_default = self.get_config_default(name=name)
        configuration = config_default
        # log.debug(f"configuration={configuration}")

        # Replace prioritizations from dict
        for config in config_dict:
            configuration[config] = config_dict[config]

        # Replace prioritizations from file
        config_file = full_path(config_file)
        if config_file:
            if os.path.exists(config_file):
                with open(config_file) as config_file_content:
                    config_file_dict = yaml.safe_load(config_file_content)
                for config in config_file_dict:
                    configuration[config] = config_file_dict[config]
            else:
                msg_error = f"Config '{name}' file '{config_file}' does NOT exist"
                log.error(msg_error)
                raise ValueError(msg_error)

        return configuration

    ################
    # Update table #
    ################

    def update_table_strategy(
        self,
        dest_table: str,
        sources: list,
        mode: str = "append",
        separator: str = ";",
        physical_order: bool = True,
        cleanup: bool = False,
        strategy: str = "ctas",
        chunking: bool = True,
        chunk_size: int = None,
        upper_case: bool = False,
        chromosomes: list = None,
    ) -> None:
        """
        Update dest_table using multiple sources via CTAS.

        :param dest_table: The `dest_table` parameter is a string that represents the name of the
        destination table that you want to update.
        :type dest_table: str
        :param sources: The `sources` parameter is a list of dictionaries, where each dictionary
        represents a source table and the columns to be updated in the destination table. Each
        dictionary should have the following keys:
            - "table": The name of the source table.
            - "join_keys": A list of column names that will be used to join the source table with the
              destination table.
            - "columns": A dictionary where the keys are the destination column names and the values
              are dictionaries with the following keys:
                - "columns": A list of column names from the source table that will be used to
                  update the destination column.
                - "mode": The mode of updating the column, either "append" or "replace".
                - "separator": The separator to use when appending values (only applicable for "append"
                  mode).
        sources example:
            sources = [
                {   # Default source to append INFO column with a table source as a view
                    "table": "calculation_view_name",
                    "join_keys": ["#CHROM", "POS", "REF", "ALT"],
                    "columns": {
                        "INFO": {"columns": ["INFO"], "mode": "append", "separator": ";"}
                    },
                },
                {   # Calculation view source to update INFO, and create AF and ALT columns
                    "table": "calculation_view_name",
                    "join_keys": ["#CHROM", "POS", "REF", "ALT"],
                    "columns": {
                        "INFO": {"columns": ["INFO"], "mode": "append", "separator": ";"},
                        "AF": {"columns": ["REF"], "mode": "replace"},
                        "ALT": {"columns": ["REF"], "mode": "replace"},
                    },
                },
                {   # Clinvar source to add in a new column CLINVAR with concatenation of 2 columns CLNSIG and CLNID
                    "table": "clinvar_table",
                    "join_keys": ["#CHROM", "POS", "REF", "ALT"],
                    "columns": {
                        "CLINVAR": {"columns": ["CLNSIG", "CLNID"], "mode": "append", "separator": "|"},
                    },
                },
            ]
        :type sources: list
        :param mode: The `mode` parameter determines how the columns in the destination table will
        be updated. It can take two values: "append" or "replace". The default value is "append".
            - If `mode` is set to "append", the values from the source tables will be concatenated to the existing
              values in the destination table, separated by the specified `separator`.
            - If `mode` is set to "replace", the values from the source tables will replace the existing values in the
              destination table.
        :type mode: str
        :param separator: The `separator` parameter is a string that specifies the character or
        sequence of characters used to separate values when appending them together. The default value is
        a semicolon (";"). This parameter is only applicable when the `mode` parameter is set to "append".
        :type separator: str
        :param physical_order: The `physical_order` parameter is a boolean value that determines
        whether the resulting table should have a physical order based on the row number. If set to
        `True`, a `_rowid` column will be added to the resulting table, which will contain the row number for each row. If
        set to `False`, the `_rowid` column will not be included. The default value is `True`.
        :type physical_order: bool
        :param cleanup: The `cleanup` parameter is a boolean value that determines whether to drop
        the temporary table created during the update process. If set to `True`, the temporary table
        will be dropped after the update is complete. If set to `False`, the temporary table will be
        retained. The default value is `False`.
        :type cleanup: bool
        :param strategy: The `strategy` parameter determines the method used to update the
        destination table. It can take two values: "ctas" or "merge". The default value is "ctas".
            - If `strategy` is set to "ctas", the update will be performed using a "Create Table As Select" (CTAS) approach. This
              involves creating a new temporary table with the updated data and then replacing the original table with the
              temporary table.
            - If `strategy` is set to "update", the update will be performed using an "UPDATE" statement. This involves directly modifying the existing
              rows in the destination table based on the data from the source tables.
        :type strategy: str
        :param chunking: The `chunking` parameter is a boolean value that determines whether to
        process the update operation in chunks or not. If set to `True`, the update will be
        performed in smaller chunks of data, which can help improve performance and reduce memory
        usage. If set to `False`, the update will be performed on the entire dataset at once. The
        default value is `True`.
        :type chunking: bool
        :param chunk_size: The `chunk_size` parameter is an optional integer that specifies the
        size of the chunks to be processed during the update operation. If not provided, it will
        default to a value from the configuration file or a predefined constant `DEFAULT_CHUNK_SIZE`.
        :type chunk_size: int or None
        :param upper_case: The `upper_case` parameter is a boolean value that determines whether to
        convert the values of the join keys to uppercase before performing the join operation. If set
        to `True`, the join keys will be converted to uppercase using the `upper` function. If set to
        `False`, the join keys will be used as they are. The default value is `False`.
        :type upper_case: bool
        :param chromosomes: The `chromosomes` parameter is a list of chromosome names (strings)
        that can be used to filter the data during the update process. If provided, only the rows
        corresponding to the specified chromosomes will be considered for the update. If not
        provided, all chromosomes will be included in the update process.
        :type chromosomes: list or None

        :return: The function `update_table` does not return anything. It performs an update
        operation on a destination table using multiple source tables and their specified columns.

        """

        # LOG
        log.debug(f"Update table '{dest_table}' using strategy '{strategy}'...")

        conn = self.get_connexion()

        # Chunk size
        if chunk_size is None:
            chunk_size = self.get_config().get("chunk_size", DEFAULT_CHUNK_SIZE)

        # chunking desactivated because of some specific tables (e.g., operations on multiple lines) ???
        # chunking = False
        # Due to creation of entire table before update, chunking should be ok
        # chunking = True

        # Upper case function for ALT and REF
        if upper_case:
            upper_func = "upper"
        else:
            upper_func = ""

        # Default configuration
        # source columns
        default_source_columns = {"INFO": {"columns": ["INFO"]}}
        # key columns
        default_join_keys = ["#CHROM", "POS", "REF", "ALT"]

        # --- Existing dest columns ---
        dest_cols = self.get_columns(dest_table)

        # --- Compute required output columns ---
        required_dest_cols = list(dest_cols)
        for src in sources:
            for dest_col in src.get("columns", default_source_columns):
                if dest_col not in required_dest_cols:
                    required_dest_cols.append(dest_col)

        column_exprs = []
        join_clauses = []
        where_clauses = {}
        where_column_exprs = []
        where_column_set = []
        list_of_sources_table = {}

        # Build JOINs
        for src in sources:

            join_keys = src.get("join_keys", default_join_keys)

            # Find type of each columns in join keys in des_table table
            join_keys_type = {}
            dest_table_columns_type = self.get_columns_type(dest_table)
            for join_key in join_keys:
                join_keys_type[join_key] = dest_table_columns_type.get(
                    join_key, "VARCHAR"
                )

            # Table source
            src_table = src.get("table", None)

            if src_table:

                # List of souce table
                list_of_sources_table[src_table] = True

                # Join clause
                join_clause = f"""
                    LEFT JOIN {src.get("table")}
                    ON {' AND '.join([
                        f'd."{k}" = {src_table}."{k}"' if join_keys_type.get(k, "").upper() not in ["VARCHAR", "TEXT"]
                        else f'{upper_func}(d."{k}") = {upper_func}({src_table}."{k}")'
                        for k in join_keys
                        ])}
                """

                # -- CTAS strategy --

                join_clauses.append(join_clause)

                # -- Update strategy --

                # where clause
                where_clause = (
                    f""" {' AND '.join([f'd."{k}" = n."{k}"' for k in join_keys])} """
                )

                # Store join clauses and keys for each where clause

                # if where clause not exist, create it
                if where_clause not in where_clauses:
                    where_clauses[where_clause] = {
                        "join_clause": [],
                        "join_keys": [],
                    }

                # Append join clause and keys
                where_clauses[where_clause]["join_clause"].append(join_clause)
                where_clauses[where_clause]["join_keys"] += join_keys

        # Rowid for deterministic physical ordering and chunking
        rowid_expr = ", row_number() OVER () AS _rowid"

        # Helper to normalize empty values
        def normalize(val):
            return f"NULLIF(NULLIF({val}, ''), '.')"

        # Build column expressions
        for col in required_dest_cols:

            # --- Sources contributing to this column ---
            update_sources = [
                src
                for src in sources
                if col in src.get("columns", default_source_columns)
            ]

            # --- No update source ---
            if not update_sources:
                if col in dest_cols:
                    column_exprs.append(f'd."{col}"')
                else:
                    column_exprs.append(f'NULL AS "{col}"')
                continue

            # --- Update column ---
            source_values = [
                (
                    normalize(f'{src.get("table")}."{col_name}"')
                    if src.get("table", None)
                    else None
                )
                for src in update_sources
                for col_name in src.get("columns", default_source_columns)[col][
                    "columns"
                ]
            ]

            # --- Remove None values (when no table defined) ---
            source_values = [v for v in source_values if v is not None]

            # --- Determine mode for this column ---
            columns_mode = {
                col: src.get("columns", {}).get(col, {}).get("mode", mode)
                for src in update_sources
            }

            # --- Determine separator for this column ---
            columns_separator = {
                col: src.get("columns", {}).get(col, {}).get("separator", separator)
                for src in update_sources
            }

            # ---------------------------
            # MODE REPLACE
            # ---------------------------

            if columns_mode[col] == "replace":

                if col in dest_cols:
                    replace_candidates = source_values + [f'd."{col}"']
                else:
                    replace_candidates = source_values

                column_exprs.append(
                    f"COALESCE({', '.join(replace_candidates)}) AS \"{col}\""
                )

                continue

            # ---------------------------
            # MODE APPEND (concat + clean)
            # ---------------------------

            pieces = []

            # existing column first (if exists)
            if col in dest_cols:
                pieces.append(normalize(f'd."{col}"'))

            # then all source values
            pieces.extend(source_values)

            # raw concat
            raw_concat = f"""
                concat_ws('{columns_separator[col]}',
                    {", ".join(pieces)}
                )
            """

            # cleanup of semicolons
            if cleanup:
                clean = f"""
                    TRIM(
                        REGEXP_REPLACE(
                            {raw_concat},
                            '{columns_separator[col]}{{2,}}',
                            '{columns_separator[col]}'
                        ),
                        '{columns_separator[col]}'
                    ) AS "{col}"
                """
            else:
                clean = f"""
                            {raw_concat} AS "{col}"
                """

            column_exprs.append(clean)

            where_column_exprs.append(clean)
            where_column_set.append(col)

        # ---------------------------
        # UPDATE
        # ---------------------------

        if strategy == "update":

            # For each where clause, create temp table and update dest_table
            for where_clause in where_clauses.keys():

                # Build join clause
                where_clause_join_clauses = where_clauses.get(where_clause, {}).get(
                    "join_clause", ""
                )

                # Build join keys
                where_clause_join_keys = [
                    f'd."{k}"'
                    for k in set(
                        where_clauses.get(where_clause, {}).get("join_keys", "")
                    )
                ]

                # Create temp table with required columns
                new_dest_table = f"tmp_new_{dest_table}_{get_random(10)}"

                sql = f"""
                    CREATE VIEW {new_dest_table} AS
                        (
                            SELECT
                                {", ".join(where_clause_join_keys)},
                                {", ".join(where_column_exprs)},
                                row_number() OVER () AS _rowid
                            FROM {dest_table} d
                            {" ".join(where_clause_join_clauses)}
                        );
                """

                # log.debug(f"SQL: {sql}")

                # Execute Update Creation View
                log.debug("Execute Update Creation View...")
                conn.execute(sql)

                # Update dest_table with new_dest_table
                log.debug(f"Updating table {dest_table}...")

                if chunking:

                    # Update table {dest_table} with new table
                    # split update by chunk (chunk_size) on _rowid column to avoid transaction too large
                    max_rowid = conn.execute(
                        f"SELECT count(1) AS max_rowid FROM {dest_table}"
                    ).df()["max_rowid"][0]
                    # Handle NaN case or None
                    if max_rowid is None or math.isnan(max_rowid):
                        max_rowid = 0

                    # Range of rowid
                    range_rowid = range(1, int(max_rowid) + 1, chunk_size)

                    # Process chunks
                    chunk_i = 0
                    if max_rowid >= 0:

                        for chunk_start in range_rowid:

                            # Chunk info
                            chunk_i += 1
                            chunk_end = min(
                                chunk_start + chunk_size - 1, int(max_rowid)
                            )
                            log.debug(
                                f"Updating table {dest_table} - rows between {chunk_start} and {chunk_end} [{chunk_i}/{len(range_rowid)}][{chunk_i/len(range_rowid)*100:.2f}%]..."
                            )

                            # Update query without chromosome chunking
                            if chromosomes is None:
                                update_query = f"""
                                    UPDATE {dest_table} AS d
                                    SET
                                        {", ".join([f'"{col}" = n."{col}"' for col in where_column_set])}
                                        FROM {new_dest_table} AS n
                                    WHERE {" ".join(set(where_clauses.keys()))}
                                    -- AND n."#CHROM" = 'chr1'
                                    AND n._rowid BETWEEN {chunk_start} AND {chunk_end}
                                """
                                # log.debug(f"update_query:\n{update_query}")
                                conn.execute(update_query)

                            # Update query with chromosome chunking
                            else:

                                # if chromosomes not provided, as an empty list, get list of chromosomes from dest_table
                                if len(chromosomes) == 0:

                                    # List of chromosomes in dest_table with query
                                    query_chromosomes = f"""
                                        SELECT DISTINCT "#CHROM" AS chrom FROM {dest_table}
                                        ORDER BY TRY_CAST(regexp_extract("#CHROM", '(d+)') AS INTEGER) NULLS LAST, "#CHROM"
                                    """
                                    df_chromosomes = conn.execute(
                                        query_chromosomes
                                    ).df()
                                    chromosomes = df_chromosomes["chrom"].tolist()

                                # Process each chromosome
                                for chrom in chromosomes:

                                    log.debug(
                                        f"Updating table {dest_table} - rows between {chunk_start} and {chunk_end} [{chunk_i}/{len(range_rowid)}][{chunk_i/len(range_rowid)*100:.2f}%] - chromosome {chrom}..."
                                    )

                                    update_query = f"""
                                        UPDATE {dest_table} AS d
                                        SET
                                            {", ".join([f'"{col}" = n."{col}"' for col in where_column_set])}
                                            FROM {new_dest_table} AS n
                                        WHERE {" ".join(set(where_clauses.keys()))}
                                        AND n."#CHROM" = '{chrom}'
                                        AND n._rowid BETWEEN {chunk_start} AND {chunk_end}
                                    """
                                    # log.debug(f"update_query:\n{update_query}")
                                    conn.execute(update_query)

                else:

                    # Update table {dest_table} with new table
                    update_query = f"""
                        UPDATE {dest_table} AS d
                        SET
                            {", ".join([f'"{col}" = n."{col}"' for col in where_column_set])}
                            FROM {new_dest_table} AS n
                            WHERE {" ".join(set(where_clauses.keys()))}
                        """
                    # log.debug(f"update_query:\n{update_query}")
                    conn.execute(update_query)

                # Cleanup
                log.debug("Cleanup temporary view...")
                self.remove_tables_or_views(tables=[new_dest_table])

        # ---------------------------
        # CTAS
        # ---------------------------

        elif strategy == "ctas":

            # Log
            log.debug("Execute CTAS...")

            # Spilling
            # If DuckDB has spilled to disk during previous operations, use Parquet intermediate files for CTAS
            duckdb_temp_directory = self.get_connexion_config().get(
                "temp_directory", ".tmp"
            )
            duckdb_spilled = duckdb_has_spilled(duckdb_temp_directory)
            log.debug(f"DuckDB spilled: {duckdb_spilled}")

            if duckdb_spilled:
                log.debug(
                    "DuckDB has spilled to disk during previous operations. "
                    "CTAS operation will use Parquet intermediate files to avoid out of memory errors. "
                    "This may slow down the operation. "
                    "Consider increasing DuckDB memory limit or adding more RAM to the system."
                )
                ctas_parquet_folder = os.path.join(
                    self.get_tmp_dir(),  # Use Howard temp directory
                    # duckdb_temp_directory, # Use duckDB temp directory to ensure same disk
                    f"ctas_parquet_{get_random(10)}.partition.parquet",
                )
                os.makedirs(ctas_parquet_folder, exist_ok=True)
                log.debug(f"CTAS Parquet folder: {ctas_parquet_folder}")

            if chunking:

                # Create new dest table with required columns OK
                new_dest_table = f"tmp_new_{dest_table}_{get_random(10)}"

                schema_sql = f"""
                    CREATE TABLE {new_dest_table} AS
                    SELECT * FROM {dest_table} WHERE 1=0
                """
                # log.debug(f"Schema SQL for update_table:\n{schema_sql}")
                conn.execute(schema_sql)

                # Update table {dest_table} with new table
                # split update by chunk (chunk_size) on _rowid column to avoid transaction too large
                max_rowid = conn.execute(
                    f"SELECT count(1) AS max_rowid FROM {dest_table}"
                ).df()["max_rowid"][0]
                # Handle NaN case or None
                if max_rowid is None or math.isnan(max_rowid):
                    max_rowid = 0

                # SQL template for update with CTAS strategy
                sql = """
                    INSERT INTO {new_dest_table}
                    WITH d AS (
                        SELECT * {rowid_expr}
                        FROM {dest_table}
                        {dest_where_chrom}
                        QUALIFY _rowid BETWEEN {start} AND {end}
                    ),
                    joined AS (
                        SELECT
                            {join_colunls_exprs}, _rowid
                        FROM d
                        {join_clauses}
                        {join_where_chrom}
                    ),
                    dedup AS (
                        SELECT *
                        FROM joined
                        QUALIFY row_number() OVER (
                            PARTITION BY _rowid
                            ORDER BY "INFO" DESC NULLS LAST
                        ) = 1
                    )
                    SELECT
                        {dest_cols}
                    FROM dedup
                    
                    {order_by}
                """

                # Spilled CTAS SQL template
                sql_parquet = """
                    COPY (
                    WITH d AS (
                        SELECT * {rowid_expr}
                        FROM {dest_table}
                        {dest_where_chrom}
                        QUALIFY _rowid BETWEEN {start} AND {end}
                    ),
                    joined AS (
                        SELECT
                            {join_colunls_exprs}, _rowid
                        FROM d
                        {join_clauses}
                        {join_where_chrom}
                    ),
                    dedup AS (
                        SELECT *
                        FROM joined
                        QUALIFY row_number() OVER (
                            PARTITION BY _rowid
                            ORDER BY "INFO" DESC NULLS LAST
                        ) = 1
                    )
                    SELECT
                        {dest_cols}
                    FROM dedup
                    {order_by}
                    ) TO '{ctas_parquet_folder}/part_{chunk_i}.parquet' (FORMAT 'parquet')
                """

                # Range of rowid
                range_rowid = range(1, int(max_rowid) + 1, chunk_size)

                # Process chunks
                chunk_i = 0
                for chunk_start in range_rowid:

                    # Chunk info
                    chunk_i += 1
                    chunk_end = chunk_start + chunk_size - 1
                    log.debug(
                        f"Execute CTAS on table {dest_table} - [{chunk_i}/{len(range_rowid)}][{chunk_i/len(range_rowid)*100:.2f}%] - rows between {chunk_start} and {chunk_end}..."
                    )

                    # Update query without chromosome chunking
                    if chromosomes is None:

                        if duckdb_spilled:

                            # Prepare CTAS SQL to Parquet file
                            sql_query = sql_parquet.format(
                                rowid_expr=rowid_expr,
                                dest_table=dest_table,
                                start=chunk_start,
                                end=chunk_end,
                                join_colunls_exprs=", ".join(column_exprs),
                                join_clauses=" ".join(join_clauses),
                                dest_cols=", ".join(f'"{col}"' for col in dest_cols),
                                order_by="ORDER BY _rowid" if physical_order else "",
                                dest_where_chrom="",
                                join_where_chrom="",
                                ctas_parquet_folder=ctas_parquet_folder,
                                chunk_i=chunk_i,
                            )

                        else:

                            # Prepare CTAS SQL to table
                            sql_query = sql.format(
                                new_dest_table=new_dest_table,
                                rowid_expr=rowid_expr,
                                dest_table=dest_table,
                                start=chunk_start,
                                end=chunk_end,
                                join_colunls_exprs=", ".join(column_exprs),
                                join_clauses=" ".join(join_clauses),
                                dest_cols=", ".join(f'"{col}"' for col in dest_cols),
                                order_by="ORDER BY _rowid" if physical_order else "",
                                dest_where_chrom="",
                                join_where_chrom="",
                            )

                        # log.debug(f"CTAS SQL for update_table:\n{sql_query}")
                        conn.execute(sql_query)

                    else:

                        # if chromosomes not provided, as an empty list, get list of chromosomes from dest_table
                        if len(chromosomes) == 0:

                            # List of chromosomes in dest_table with query
                            query_chromosomes = f"""
                                SELECT DISTINCT "#CHROM" AS chrom FROM {dest_table}
                                ORDER BY TRY_CAST(regexp_extract("#CHROM", '(d+)') AS INTEGER) NULLS LAST, "#CHROM"
                            """
                            df_chromosomes = conn.execute(query_chromosomes).df()
                            chromosomes = df_chromosomes["chrom"].tolist()

                        # Process each chromosome
                        for chrom in chromosomes:

                            if duckdb_spilled:

                                sql_query = sql_parquet.format(
                                    rowid_expr=rowid_expr,
                                    dest_table=dest_table,
                                    start=chunk_start,
                                    end=chunk_end,
                                    join_colunls_exprs=", ".join(column_exprs),
                                    join_clauses=" ".join(join_clauses),
                                    dest_cols=", ".join(
                                        f'"{col}"' for col in dest_cols
                                    ),
                                    order_by=(
                                        "ORDER BY _rowid" if physical_order else ""
                                    ),
                                    dest_where_chrom=f"""WHERE "#CHROM" LIKE '{chrom}'""",
                                    join_where_chrom=f"""WHERE d."#CHROM" LIKE '{chrom}'""",
                                    ctas_parquet_folder=ctas_parquet_folder,
                                    chunk_i=chunk_i,
                                )

                            else:

                                sql_query = sql.format(
                                    new_dest_table=new_dest_table,
                                    rowid_expr=rowid_expr,
                                    dest_table=dest_table,
                                    start=chunk_start,
                                    end=chunk_end,
                                    join_colunls_exprs=", ".join(column_exprs),
                                    join_clauses=" ".join(join_clauses),
                                    dest_cols=", ".join(
                                        f'"{col}"' for col in dest_cols
                                    ),
                                    order_by=(
                                        "ORDER BY _rowid" if physical_order else ""
                                    ),
                                    dest_where_chrom=f"""WHERE "#CHROM" LIKE '{chrom}'""",
                                    join_where_chrom=f"""WHERE d."#CHROM" LIKE '{chrom}'""",
                                )

                            # Log
                            # log.debug(f"CTAS SQL for update_table:\n{sql_query}")

                            # Execute CTAS
                            log.debug(
                                f"Execute CTAS on table {dest_table} - [{chunk_i}/{len(range_rowid)}][{chunk_i/len(range_rowid)*100:.2f}%] - rows between {chunk_start} and {chunk_end} - chromosome {chrom}..."
                            )
                            conn.execute(sql_query)

            else:

                if duckdb_spilled:

                    # CTAS to Parquet file
                    sql = f"""
                    COPY (
                        WITH d AS (
                            SELECT * {rowid_expr}
                            FROM {dest_table}
                        )
                        SELECT
                            {", ".join(column_exprs)}
                        FROM d
                        {" ".join(join_clauses)}
                        {"ORDER BY d._rowid" if physical_order else ""}
                    ) TO '{ctas_parquet_folder}/part_1.parquet' (FORMAT 'parquet')
                    """

                else:

                    # Create new dest table with required columns OK
                    new_dest_table = f"tmp_new_{dest_table}_{get_random(10)}"
                    sql = f"""
                        CREATE TABLE {new_dest_table} AS
                        WITH d AS (
                            SELECT * {rowid_expr}
                            FROM {dest_table}
                        )
                        SELECT
                            {", ".join(column_exprs)}
                        FROM d
                        {" ".join(join_clauses)}
                        {"ORDER BY d._rowid" if physical_order else ""}
                    """

                # log.debug(f"CTAS SQL for update_table:\n{sql}")

                # Execute CTAS
                log.debug(f"Execute CTAS on table {dest_table}...")
                conn.execute(sql)

            if duckdb_spilled:

                # Create view with all parquet files as new_dest_table
                log.debug(
                    f"Execute CTAS on table {dest_table} - Create view from Parquet files and replace table..."
                )

                # Find all parquet files
                glob_pattern = os.path.join(ctas_parquet_folder, "**", "part_*.parquet")
                parquet_files = sorted(
                    glob.glob(glob_pattern, recursive=True),
                    key=lambda p: (os.path.getmtime(p)),
                )

                # Remove tables and create table from parquet files
                conn.execute(f"DROP TABLE IF EXISTS {new_dest_table}")
                conn.execute(f"DROP TABLE IF EXISTS {dest_table}")
                sql_parquet_create_view = f"CREATE TABLE {dest_table} AS SELECT * FROM read_parquet({parquet_files})"
                # log.debug(f"SQL Parquet Create View:\n{sql_parquet_create_view}")
                conn.execute(sql_parquet_create_view)
                remove_if_exists([ctas_parquet_folder])

            else:

                # Replace dest_table with new_dest_table
                log.debug(f"Execute CTAS on table {dest_table} - Replace table...")
                conn.execute(f"DROP TABLE {dest_table}")
                conn.execute(f"ALTER TABLE {new_dest_table} RENAME TO {dest_table}")

        # ---------------------------
        # FAST
        # ---------------------------

        elif strategy == "fast":

            log.debug("Execute FAST strategy...")

            # Create new dest table name
            new_dest_table = f"tmp_new_{dest_table}_{get_random(10)}"

            # Rename dest_table to new_dest_table
            try:
                conn.execute(f"ALTER TABLE {dest_table} RENAME TO {new_dest_table}")
            except Exception as e:
                log.debug(
                    f"Failed to rename table {dest_table} to {new_dest_table}: {e}"
                )
                log.debug(f"trying to rename view {dest_table} to {new_dest_table}...")
                try:
                    conn.execute(f"ALTER VIEW {dest_table} RENAME TO {new_dest_table}")
                except Exception as e:
                    log.error(
                        f"Failed to rename view {dest_table} to {new_dest_table}: {e}"
                    )
                    raise e

            # Create view with updated data
            sql = f"""
                    CREATE VIEW {dest_table} AS
                    WITH d AS (
                        SELECT * {rowid_expr}
                        FROM {new_dest_table}
                    ),
                    joined AS (
                        SELECT
                            {", ".join(column_exprs)}, _rowid
                        FROM d
                        {" ".join(join_clauses)}
                    ),
                    dedup AS (
                        SELECT *
                        FROM joined
                        QUALIFY row_number() OVER (
                            PARTITION BY _rowid
                            -- ORDER BY "INFO" DESC NULLS LAST
                        ) = 1
                    )
                    SELECT
                        {", ".join(f'"{col}"' for col in dest_cols)}
                    FROM dedup

                    -- {"ORDER BY dedup._rowid" if physical_order else ""}
                """

            # log.debug(f"SQL for FAST strategy:\n{sql}")
            conn.execute(sql)

        else:
            log.error(f"Strategy '{strategy}' NOT available")
            raise ValueError(f"Strategy '{strategy}' NOT available")

        # Remove source tables if cleanup
        if strategy != "fast":
            log.debug("Cleanup source tables/views...")
            self.remove_tables_or_views(tables=list(list_of_sources_table.keys()))

    def update_table(
        self,
        dest_table: str,
        sources: list,
        mode: str = "append",
        separator: str = ";",
        physical_order: bool = True,
        cleanup: bool = False,
        force_strategy: str = None,
        chunk_size: int = None,
        upper_case: bool = False,
        samples: int = 100000,
        chromosomes: list = None,
        only_strategy: bool = False,
    ) -> None:
        """
        Update dest_table using multiple sources via CTAS or hybrid UPDATE.
        Heuristic chooses between UPDATE and CTAS based on number of columns to update.

        :param dest_table: destination table to update
        :type dest_table: str
        :param sources: list of source dicts
            sources example:
            sources = [
                {   # Default source to append INFO column with a table source as a view
                    "table": "calculation_view_name",
                    "join_keys": ["#CHROM", "POS", "REF", "ALT"],
                    "columns": {
                        "INFO": {"columns": ["INFO"], "mode": "append", "separator": ";"}
                    },
                },
                {   # Calculation view source to update INFO, and create AF and ALT columns
                    "table": "calculation_view_name",
                    "join_keys": ["#CHROM", "POS", "REF", "ALT"],
                    "columns": {
                        "INFO": {"columns": ["INFO"], "mode": "append", "separator": ";"},
                        "AF": {"columns": ["REF"], "mode": "replace"},
                        "ALT": {"columns": ["REF"], "mode": "replace"},
                    },
                },
                {   # Clinvar source to add in a new column CLINVAR with concatenation of 2 columns CLNSIG and CLNID
                    "table": "clinvar_table",
                    "join_keys": ["#CHROM", "POS", "REF", "ALT"],
                    "columns": {
                        "CLINVAR": {"columns": ["CLNSIG", "CLNID"], "mode": "append", "separator": "|"},
                    },
                },
            ]
        :type sources: list
        :param mode: "append" or "replace"
        :type mode: str
        :param mode: mode of update ("append" or "replace")
        :type mode: str
        :param separator: separator for concatenation
        :type separator: str
        :param physical_order: keep physical order of dest_table
        :type physical_order: bool
        :param cleanup: clean concatenated values (remove duplicates, start/end separators)
        :type cleanup: bool
        :param force_use_ctas: force use of CTAS (True), UPDATE (False) or heuristic (None)
        :type force_use_ctas: bool or None
        :param ctas_threshold: threshold to use CTAS based on number of columns to update / total columns ratio
        :type ctas_threshold: float
        :param chunk_size: threshold to use CTAS based on dest table size (number of rows)
        :type chunk_size: int
        :param upper_case: upper case join keys
        :type upper_case: bool
        :param samples: number of samples to use for statistics
        :type samples: int
        :param chromosomes: list of chromosomes to chunk update
        :type chromosomes: list or None
        :param only_strategy: only return chosen strategy without performing update
        :type only_strategy: bool

        :return: None
        :rtype: None
        """

        conn = self.get_connexion()

        # Chunk size
        if chunk_size is None:
            chunk_size = self.get_config().get("chunk_size", DEFAULT_CHUNK_SIZE)

        # Default configuration
        default_source_columns = {"INFO": {"columns": ["INFO"]}}
        default_join_keys = ["#CHROM", "POS", "REF", "ALT"]

        # --- Existing dest columns and update columns ---
        update_cols = set()
        for src in sources:
            update_cols.update(src.get("columns", default_source_columns).keys())

        # --- Available memory and threads ---
        memory = int(self.get_memory("1G", available=True).replace("G", ""))
        threads = self.get_threads("1")

        # Heuristics to choose between CTAS or UPDATE
        log.debug(f"force_strategy={force_strategy}")
        if force_strategy is None:

            # --- Choose safe strategy -----------
            strategy, reasoning = choose_update_strategy_safe(
                dest_total_rows=None,
                dest_total_cols=None,
                avg_row_size_bytes=None,
                update_cols_count=len(update_cols),
                update_row_ratio=None,  # unknown yet
                ram_available_gb=memory,
                chunk_size=chunk_size,
                conn=conn,
                dest_table=dest_table,
                sources=sources,
                default_join_keys=default_join_keys,
                samples=samples,
                threads=threads,
            )

        else:
            strategy = force_strategy
            reasoning = [f"Strategy forced to '{strategy}'"]

        # Log chosen strategy
        log.debug(f"Chosen strategy: {strategy}")
        for reason in reasoning:
            log.debug(reason)

        # Return strategy if only_strategy
        if only_strategy:
            log.info(f"Chosen strategy: {strategy}")
            return strategy

        # Perform update with chosen strategy
        self.update_table_strategy(
            dest_table=dest_table,
            sources=sources,
            mode=mode,
            separator=separator,
            physical_order=physical_order,
            cleanup=cleanup,
            strategy=strategy.lower(),
            chunk_size=chunk_size,
            upper_case=upper_case,
            chromosomes=chromosomes,
        )

        return None

    ##############################
    # Function format annotation #
    ##############################

    def annotation_format_to_table(
        self,
        annotation_field: str = "ANN",
        annotation_id: str = "Feature_ID",
        view_name: str = "transcripts",
        view_type: str = "view",
        column_rename: dict = {},
        column_clean: bool = False,
        column_case: str = None,
        column_split: str = "&",
    ) -> str:
        """
        Converts annotation data from a VCF file into a structured table format, ensuring unique values
        and creating a temporary table for further processing or analysis.

        Args:
            annotation_field (str, optional): The field in the VCF file that contains the annotation information
                for each variant. Defaults to "ANN".
            annotation_id (str, optional): The identifier for the annotation feature, used as a column name
                in the resulting table or view. Defaults to "Feature_ID".
            view_name (str, optional): The name of the temporary table that will be created to store the transformed
                annotation data. Defaults to "transcripts".
            view_type (str, optional): The type of the view to be created. Defaults to "view".
            column_rename (dict, optional): A dictionary to specify custom renaming for columns. By providing key-value
                pairs in this dictionary, you can rename specific columns in the resulting table or view.
            column_clean (bool, optional): A flag to determine whether the annotation field should undergo a cleaning process.
                If set to True, the function will clean the annotation field before further processing. Defaults to False.
            column_case (str, optional): Specifies the case transformation to be applied to the column names extracted from
                the annotation data. It allows you to set the case of the column names to either lowercase or uppercase.
            column_split (str, optional): The separator to split field values. Defaults to "&". Set to None to disable splitting.

        Returns:
            str: The name of the view created, which is stored in the variable `view_name`.
        """

        # annotation_id original name
        annotation_id_original = annotation_id

        # Transcript annotation
        if column_rename:
            annotation_id = column_rename.get(annotation_id, annotation_id)

        if column_clean:
            annotation_id = clean_annotation_field(annotation_id)

        # Prefix
        prefix = self.get_explode_infos_prefix()
        if prefix:
            prefix = "INFO/"

        # Variants table
        table_variants = self.get_table_variants()

        # Header
        vcf_reader = self.get_header()

        # Add columns
        added_columns = []
        added_columns_type = {}

        # If annotation_field exists
        if annotation_field in vcf_reader.infos:

            # Extract ANN header
            ann_description = vcf_reader.infos[annotation_field].desc
            # pattern = r"'(.+?)'"
            # pattern = r"Format: (.+?)$"
            # pattern = r"'(.+?)'|Format: (.+?)"
            pattern = r"(?:'|Format:\s*)(.+?)(?:'|$)"
            match = re.search(pattern, ann_description)
            if match:
                ann_header_match = match.group(1).split("|")
                ann_header = []
                ann_header_desc = {}
                for i in range(len(ann_header_match)):
                    ann_header_match[i] = ann_header_match[i].strip()
                    ann_header_info = "".join(
                        char for char in ann_header_match[i] if char.isalnum()
                    )
                    ann_header.append(ann_header_info)
                    ann_header_desc[ann_header_info] = ann_header_match[i]
                if not ann_header_desc:
                    raise ValueError("Invalid header description format")
            else:
                raise ValueError("Invalid header description format")

            # annotation field pattern
            annotation_field_pattern = rf"(^|;)({annotation_field})=([^;]*)?"

            annotation_fields_for_format = []
            for i, header in enumerate(ann_header_desc.values()):
                if header in [annotation_id_original]:
                    annotation_fields_for_format.append(
                        f"SPLIT_PART(annotation, '|', {i+1}) AS '{header}'"
                    )
                else:
                    annotation_fields_for_format.append(
                        f"string_agg(SPLIT_PART(annotation, '|', {i+1}), ',') AS '{header}'"
                    )

            query = f""" 
                WITH exploded_annotations AS (
                    SELECT
                        "#CHROM", POS, REF, ALT,
                        UNNEST(
                            STRING_SPLIT(
                                regexp_extract("INFO", '{annotation_field_pattern}', 3),
                                ','
                            )
                        ) AS annotation
                    FROM {table_variants}
                ),
                split_annotations AS (
                    SELECT
                        "#CHROM", POS, REF, ALT,
                        {", ".join(annotation_fields_for_format)}
                    FROM exploded_annotations
                    GROUP BY "#CHROM", POS, REF, ALT, "{annotation_id_original}"
                )
                SELECT * FROM split_annotations
                LIMIT 10000
                """
            dataframe_annotation_format = self.get_query_to_df(query=query)

            # Init
            query_list_keys = []
            key_i = 0

            for key in dataframe_annotation_format.keys():

                if key in ann_header_desc.values():

                    # Key
                    key_i += 1
                    key_clean = key

                    # key rename
                    if column_rename:
                        key_clean = column_rename.get(key_clean, key_clean)

                    # key clean
                    if column_clean:
                        key_clean = clean_annotation_field(key_clean)

                    # Key case
                    if column_case:
                        if column_case.lower() in ["lower"]:
                            key_clean = key_clean.lower()
                        elif column_case.lower() in ["upper"]:
                            key_clean = key_clean.upper()

                    # Detect column type
                    column_type = detect_column_type(dataframe_annotation_format[key])
                    added_columns_type[key] = column_type
                    log.debug(f"Field '{key}' type detected: {column_type}")

                    # Append key to list
                    if column_split is not None:
                        query_list_keys.append(
                            f""" TRY_CAST(replace(NULLIF(SPLIT_PART(annotation, '|', {key_i}), ''), '{column_split}', ',') AS {column_type}) AS '{prefix}{key_clean}' """
                        )
                    else:
                        query_list_keys.append(
                            f""" TRY_CAST(NULLIF(SPLIT_PART(annotation, '|', {key_i}), '') AS {column_type}) AS '{prefix}{key_clean}' """
                        )

            # Create temporary table
            query_create_view = f"""
                CREATE {view_type} {view_name} AS (
                    WITH exploded_annotations AS (
                        SELECT
                            "#CHROM", POS, REF, ALT,
                            UNNEST(
                                STRING_SPLIT(
                                    regexp_extract("INFO", '{annotation_field_pattern}', 3),
                                    ','
                                )
                            ) AS annotation
                        FROM {table_variants}
                    ),
                    split_annotations AS (
                        SELECT
                            "#CHROM", POS, REF, ALT,
                            {", ".join(query_list_keys)},
                        FROM exploded_annotations
                    )
                    SELECT *, {annotation_id} AS 'transcript' FROM split_annotations
                )
            """
            # log.debug(f"Create view FORMAT:{query_create_view}")
            self.execute_query(query=query_create_view)

        else:

            # Return None
            view_name = None

        return view_name, added_columns, added_columns_type

    ############################
    # Rename and remove fields #
    ############################

    def rename_info_fields(
        self, fields_to_rename: dict = None, table: str = None
    ) -> dict:
        """
        The `rename_info_fields` function renames specified fields in a VCF file header and updates
        corresponding INFO fields in the variants table.

        :param fields_to_rename: The `fields_to_rename` parameter is a dictionary that contains the
        mapping of fields to be renamed in a VCF (Variant Call Format) file. The keys in the dictionary
        represent the original field names that need to be renamed, and the corresponding values
        represent the new names to which the fields should be
        :type fields_to_rename: dict
        :param table: The `table` parameter in the `rename_info_fields` function represents the name of
        the table in which the variants data is stored. This table contains information about genetic
        variants, and the function updates the corresponding INFO fields in this table when renaming
        specified fields in the VCF file header
        :type table: str
        :return: The `rename_info_fields` function returns a dictionary `fields_processed` that contains
        the original field names as keys and their corresponding new names (or None if the field was
        removed) as values after renaming or removing specified fields in a VCF file header and updating
        corresponding INFO fields in the variants table.
        """

        # Config
        config = self.get_config()

        # Access
        # access = config.get("access")
        access = self.get_access(default=None)

        # Init
        fields_processed = {
            "renamed": {},
            "removed": {},
            "not_processed": {},
            "not_found": {},
        }

        if table is None:
            table = self.get_table_variants()

        # Clasue case
        clause_case = []

        # Init
        fields_to_process = {}

        # For each field to rename or remove, one by one
        if fields_to_rename is not None:
            for field_to_rename, field_renamed in fields_to_rename.items():

                # If no field to process
                if field_to_rename != "":

                    # rename empty is remove
                    if field_renamed == "":
                        field_renamed = None

                    # Check if already to process
                    check_field_to_rename_found = False
                    for (
                        field_to_process_to_rename,
                        field_to_process_renamed,
                    ) in fields_to_process.items():
                        if field_to_rename == field_to_process_renamed:
                            # Replace filed to process
                            fields_to_process[field_to_process_to_rename] = (
                                field_renamed
                            )
                            check_field_to_rename_found = True
                    if not check_field_to_rename_found:
                        # Add field to process
                        fields_to_process[field_to_rename] = field_renamed

        if len(fields_to_process) and access not in ["RO"]:

            # Log
            log.info("Rename or remove fields...")

            # Header
            header = self.get_header()

            # For each field
            for field_to_rename, field_renamed in fields_to_process.items():

                # If to rename or remove
                if (
                    field_to_rename != field_renamed
                    and field_to_rename not in ["", None]
                    and field_renamed not in [""]
                ):

                    # Rename
                    if field_renamed is not None:

                        # Case clause
                        clause_case.append(
                            f""" WHEN k = '{field_to_rename}' THEN '{field_renamed}'  """
                        )

                        # Fields processed
                        fields_processed["renamed"][field_to_rename] = field_renamed

                        # Log
                        log.debug(
                            f"Rename or remove fields - field '{field_to_rename}' renamed as '{field_renamed}'"
                        )

                    # Remove
                    else:

                        # Case clause
                        clause_case.append(
                            f""" WHEN k = '{field_to_rename}' THEN NULL  """
                        )

                        # Fields processed
                        fields_processed["removed"][field_to_rename] = field_renamed

                        # Log
                        log.debug(
                            f"Rename or remove fields - field '{field_to_rename}' removed"
                        )

                    # Header
                    if field_to_rename in header.infos:

                        # Rename header if to rename
                        if field_renamed is not None:
                            header.infos[field_renamed] = vcf.parser._Info(
                                field_renamed,
                                header.infos[field_to_rename].num,
                                header.infos[field_to_rename].type,
                                header.infos[field_to_rename].desc,
                                header.infos[field_to_rename].source,
                                header.infos[field_to_rename].version,
                                header.infos[field_to_rename].type_code,
                            )

                        # Remove header, if rename or remove
                        del header.infos[field_to_rename]

                    else:

                        # Log
                        log.warning(
                            f"Rename or remove fields - field '{field_to_rename}' not in header"
                        )

                        # Fields processed
                        fields_processed["not_found"][field_to_rename] = field_renamed

                else:

                    # Fields processed
                    fields_processed["not_pocessed"][field_to_rename] = field_renamed

            # Process
            if len(clause_case):

                # Update query
                query_update = f"""
                    UPDATE {table}
                    SET INFO = renamed_table.INFO                       -- update INFO
                    FROM (
                        SELECT
                            "#CHROM", POS, REF, ALT,                    -- variant id
                            IFNULL(string_agg(kv, ';'), '.') AS INFO    -- INFO
                        FROM (
                            SELECT
                                "#CHROM", POS, REF, ALT,                -- variant id
                                CASE
                                    WHEN k IS NOT NULL                  -- key not null
                                    THEN
                                        CASE 
                                            WHEN v IS NOT NULL          -- value not null
                                            THEN concat(k, '=', v)      -- key-value: either String, Integer, Float
                                            ELSE k                      -- Flag
                                        END
                                    ELSE NULL                           -- remove
                                END AS kv
                            FROM (
                                SELECT "#CHROM", POS, REF, ALT,         -- variant id
                                    CASE
                                        {" ".join(clause_case)}         -- rename or remove
                                        ELSE k                          -- no change
                                    END AS k,                           -- key
                                    v                                   -- value
                                FROM (
                                    SELECT
                                        "#CHROM", POS, REF, ALT,        -- variant id
                                        string_split(kv, '=')[1] AS k,  -- key
                                        string_split(kv, '=')[2] AS v   -- value
                                    FROM (
                                        SELECT "#CHROM", POS, REF, ALT, unnest(string_split(INFO, ';')) AS kv
                                        FROM variants
                                        )
                                    )
                                )
                            )
                        GROUP BY "#CHROM", POS, REF, ALT
                    ) AS renamed_table
                    WHERE {table}."#CHROM" = renamed_table."#CHROM"     -- join
                    AND {table}."POS" = renamed_table."POS"
                    AND {table}."REF" = renamed_table."REF"
                    AND {table}."ALT" = renamed_table."ALT"
                """
                # log.debug(f"query_update={query_update}")
                self.execute_query(query_update)

        return fields_processed

    def recreate_info_fields(
        self, fields_to_rename: dict = None, table: str = None
    ) -> dict:
        """
        The `recreate_info_fields` function renames specified fields in a VCF file header and updates
        corresponding INFO fields in the variants table.

        :param fields_to_rename: The `fields_to_rename` parameter is a dictionary that contains the
        mapping of fields to be renamed in a VCF (Variant Call Format) file. The keys in the dictionary
        represent the original field names that need to be renamed, and the corresponding values
        represent the new names to which the fields should be renamed. Default {}
        :type fields_to_rename: dict
        :param table: The `table` parameter in the `recreate_info_fields` function represents the name of
        the table in which the variants data is stored. This table contains information about genetic
        variants, and the function updates the corresponding INFO fields in this table when renaming
        specified fields in the VCF file header. Default Variants table 'variants'.
        :type table: str
        :return: The `recreate_info_fields` function returns a dictionary `fields_renamed` that contains
        the original field names as keys and their corresponding new names (or None if the field was
        removed) as values after renaming or removing specified fields in a VCF file header and updating
        corresponding INFO fields in the variants table.
        """

        # Init
        config = self.get_config()

        # Access
        # access = config.get("access")
        access = self.get_access(default=None)

        # Table
        if table is None:
            table = self.get_table_variants()

        # Fields to rename or remove
        if fields_to_rename is None:
            fields_to_rename = {}

        # Fields on header
        fields_to_process = {k: k for k in self.get_header().infos.keys()}

        # For each field to rename or remove, one by one
        for field_to_rename, field_renamed in fields_to_rename.items():
            log.debug(f"{field_to_rename}, {field_renamed}")
            if field_to_rename != "" and field_renamed != "":
                log.debug(f"{field_to_rename}, {field_renamed} OK")
                # Check if already to process
                check_field_to_rename_found = False
                for (
                    field_to_process_to_rename,
                    field_to_process_renamed,
                ) in fields_to_process.items():
                    if field_to_rename == field_to_process_renamed:
                        fields_to_process[field_to_process_to_rename] = field_renamed
                        check_field_to_rename_found = True
                if not check_field_to_rename_found:
                    fields_to_process[field_to_rename] = field_renamed

        # Init
        fields_processed = {"removed": {}, "renamed": {}, "not_found": {}}

        # if fields_to_rename is not None and access not in ["RO"]:
        if fields_to_process is not None and access not in ["RO"]:

            log.info(f"Recreate INFO with {len(fields_to_process)} fields...")

            # Create view
            annotation_view_name = "annotation_view_for_recreate_infos_" + str(
                get_random()
            )

            # Header
            header = self.get_header()

            # Update query select clauses
            query_view_select_clause = []

            for field_to_rename, field_renamed in fields_to_process.items():

                # Action
                action = None

                # Header

                # Field to remove
                if field_renamed is None:

                    if field_renamed in header.infos:

                        # Remove in header
                        del header.infos[field_renamed]

                    # Action
                    action = "removed"

                # Field to rename
                elif field_to_rename in header.infos:

                    # Rename header
                    if field_renamed is not None and field_renamed != field_to_rename:
                        header.infos[field_renamed] = vcf.parser._Info(
                            field_renamed,
                            header.infos[field_to_rename].num,
                            header.infos[field_to_rename].type,
                            header.infos[field_to_rename].desc,
                            header.infos[field_to_rename].source,
                            header.infos[field_to_rename].version,
                            header.infos[field_to_rename].type_code,
                        )
                        del header.infos[field_to_rename]

                    # Update query
                    if header.infos[field_renamed].type == "Flag":
                        query_view_select_clause.append(
                            f"""
                                CASE
                                    WHEN "{field_renamed}" IS NOT NULL AND TRY_CAST("{field_renamed}" AS FLOAT) != 0
                                    THEN '{field_renamed};'
                                END
                            """
                        )
                    else:
                        query_view_select_clause.append(
                            f"""
                                CASE
                                    WHEN "{field_renamed}" IS NOT NULL
                                    THEN concat('{field_renamed}=',"{field_renamed}",';')
                                END
                            """
                        )

                    # Action
                    action = "renamed"

                else:

                    log.warning(
                        f"Rename or remove fields - field '{field_to_rename}' not in header"
                    )

                    # Action
                    action = "not_found"

                # List of renamed or removed fields
                if action:
                    fields_processed[action][field_to_rename] = field_renamed

            if len(query_view_select_clause):

                # Create view
                annotation_view_name = self.create_annotations_view(
                    table=table,
                    view=annotation_view_name,
                    view_type="view",
                    view_mode="full",
                    info_prefix_column="",
                    # info_struct_column="INFOS",
                    detect_type_list=False,
                    fields=fields_processed["renamed"].keys(),
                    fields_not_exists=True,
                    fields_forced_as_varchar=True,
                    fields_needed_all=False,
                    fields_to_rename=fields_processed["renamed"],
                    drop_view=True,
                )

                # Log
                log.info(
                    f"Recreate INFO with {len(query_view_select_clause)} found fields..."
                )

                # Query
                query = f"""
                    UPDATE {table}
                    SET
                        INFO = regexp_replace(
                                concat({", ".join(query_view_select_clause)}),
                                ';$',
                                ''
                            )
                    FROM {annotation_view_name}
                    WHERE {table}."#CHROM" = {annotation_view_name}."#CHROM"
                      AND {table}."POS" = {annotation_view_name}."POS"
                      AND {table}."REF" = {annotation_view_name}."REF"
                      AND {table}."ALT" = {annotation_view_name}."ALT"
                """
                # log.debug(f"query={query}")

                # Excecute query
                self.execute_query(query=query)

        return fields_processed

    def calculation_rename_info_fields(
        self,
        section: str = "calculation",
        fields_to_rename: dict = None,
        table: str = None,
        operation_name: str = "RENAME_INFO_FIELDS",
        **kwargs
    ) -> None:
        """
        The `calculation_rename_info_fields` function retrieves parameters from a dictionary, updates
        fields to rename and table if provided, and then calls another function to rename the fields.

        :param fields_to_rename: `fields_to_rename` is a dictionary that contains the fields to be
        renamed in a table. Each key-value pair in the dictionary represents the original field name as
        the key and the new field name as the value
        :type fields_to_rename: dict
        :param table: The `table` parameter in the `calculation_rename_info_fields` method is used to
        specify the name of the table for which the fields are to be renamed. It is a string type
        parameter
        :type table: str
        :param operation_name: The `operation_name` parameter in the `calculation_rename_info_fields`
        method is a string that specifies the name of the operation being performed. In this context, it
        is used as a default value for the operation name if not explicitly provided when calling the
        function, defaults to RENAME_INFO_FIELDS
        :type operation_name: str (optional)
        """

        param_fields_to_rename = None
        param_table = None
        operation_params, _ = self.get_operation_params(
            section=section, operation_params=kwargs, operation_name=operation_name
        )

        # param_fields_to_rename
        if param_fields_to_rename is None:
            param_fields_to_rename = (
                operation_params.get("fields_to_rename")
                or fields_to_rename
                or {}
            )

        # param_table
        if param_table is None:
            param_table = operation_params.get("table") or table or None

        # Init fields_to_rename
        if fields_to_rename is None:
            fields_to_rename = param_fields_to_rename

        # Init table
        if table is None:
            table = param_table

        renamed_fields = self.rename_info_fields(
            fields_to_rename=fields_to_rename, table=table
        )

        log.debug(f"renamed_fields:{renamed_fields}")

    def calculation_recreate_info_fields(
        self,
        section: str = "calculation",
        fields_to_rename: dict = None,
        table: str = None,
        operation_name: str = "RECREATE_INFO_FIELDS",
        **kwargs
    ) -> None:
        """
        The `calculation_recreate_info_fields` function retrieves parameters from a dictionary, recreate
        INFO fields with rename and table if provided, and then calls another function to rename the fields.

        :param fields_to_rename: `fields_to_rename` is a dictionary that contains the fields to be
        renamed in a table. Each key-value pair in the dictionary represents the original field name as
        the key and the new field name as the value
        :type fields_to_rename: dict
        :param table: The `table` parameter in the `calculation_recreate_info_fields` method is used to
        specify the name of the table for which the fields are to be renamed. It is a string type
        parameter
        :type table: str
        :param operation_name: The `operation_name` parameter in the `calculation_recreate_info_fields`
        method is a string that specifies the name of the operation being performed. In this context, it
        is used as a default value for the operation name if not explicitly provided when calling the
        function, defaults to RECREATE_INFO_FIELDS
        :type operation_name: str (optional)
        """

        param_fields_to_rename = None
        param_table = None
        operation_params, _ = self.get_operation_params(
            section=section, operation_params=kwargs, operation_name=operation_name
        )

        ### Parameters for genotype stats calculation

        # param_fields_to_rename
        if param_fields_to_rename is None:
            param_fields_to_rename = (
                operation_params.get("fields_to_rename")
                or fields_to_rename
                or {}
            )

        # param_table
        if param_table is None:
            param_table = operation_params.get("table") or table or None

        # Init fields_to_rename
        if fields_to_rename is None:
            fields_to_rename = param_fields_to_rename

        # Init table
        if table is None:
            table = param_table

        renamed_fields = self.recreate_info_fields(
            fields_to_rename=fields_to_rename, table=table
        )

        log.debug(f"renamed_fields:{renamed_fields}")

    def get_annotation_header_fields_override(self, annotations: dict) -> dict:
        """
            Return header fields override dict (Number/Type/Description), collected
            ONLY from each annotation's own options block (per-annotation level).
            No tool-level/global header_fields config is considered.

            Config location (within each annotation entry, sibling of
            "annotation_fields"):
                "annotation": {
                    "annovar": {
                        "annotations": {
                            "ALL.sites.2015_08": {
                                "annotation_fields": {
                                    "ALL.sites.2015_08": "1000genomesALL"
                                },
                                "options": {
                                    "header_fields": {
                                        "1000genomesALL": {
                                            "number": ".",
                                            "type": "Float",
                                            "description": "1000genomesALL by Annovar"
                                        }
                                    }
                                }
                            }
                        }
                    }
                }

            :param annotations: annotations dict (e.g. param.annotation.<tool>.annotations),
                where each annotation entry may define its own "options.header_fields"
            :return: dict field_name -> {"number": ..., "type": ..., "description": ...}
        """
        header_fields_override = {}
        for annotation in annotations:
            annotation_options = annotations[annotation].get("options", {})
            annotation_header_fields = annotation_options.get("header_fields", {})
            if annotation_header_fields:
                header_fields_override.update(annotation_header_fields)
        return header_fields_override
    

    def build_info_with_header_override(
        self,
        field_name: str,
        header_fields_override: dict = None,
        **kwargs
    ):
        """
        Build a vcf.parser._Info object for INFO field 'field_name', applying
        Number/Type/Description override if defined in header_fields_override
        (see get_annotation_header_fields_override).

        kwargs can include:
            - field_number: The number of values expected for the INFO field. It can be an integer or a string (e.g., "1", "A", "G", "."). If not provided, it will default to "." (unknown).
            - field_type: The data type of the INFO field. It can be one of the following: "Integer", "Float", "Flag", "Character", or "String". If not provided, it will default to "String".
            - field_description: A description of the INFO field. If not provided, it will default to "unknown".
            - field_source: The source of the INFO field. If not provided, it will default to None.
            - field_version: The version of the INFO field. If not provided, it will default to None.

        :param field_name: The `field_name` parameter is a string that represents the name of the INFO field
        for which the vcf.parser._Info object is being built. This field name is used to look up any overrides in the `header_fields_override` dictionary, which may contain custom definitions for the field's number, type, description, source, and version. If an override is found for the specified field name, it will be applied when constructing the vcf.parser._Info object
        :type field_name: str
        :param header_fields_override: The `header_fields_override` parameter is a dictionary that contains
        overrides for the header fields of a VCF (Variant Call Format) file. It allows you to specify custom values for the "number", "type", "description", "source", and "version" attributes of the INFO fields in the VCF header. If an override is provided for a specific field name, it will be used instead of the default values when building the vcf.parser._Info object for that field. If no override is provided, the default values will be used
        :type header_fields_override: dict

        :param kwargs: Additional keyword arguments that can be used to override specific attributes of the INFO field.
        :type kwargs: dict
        

        :return: vcf.parser._Info object
        """

        # Init
        if header_fields_override is None:
            header_fields_override = {}

        # Get field override for this field_name, if any
        field_override = header_fields_override.get(field_name, {})

        # field_override as lowercase keys (e.g. "number", "type", "description")
        field_override = {k.lower(): v for k, v in field_override.items()}

        # Check override type
        override_type = field_override.get("type")
        if override_type and override_type not in CODE_TYPE_MAP:
            msg_err = (
                f"Header field '{field_name}' override - type '{override_type}' "
                f"not valid (should be one of {list(CODE_TYPE_MAP.keys())})"
            )
            log.error(msg_err)
            raise ValueError(msg_err)
            
        # Apply override if defined, otherwise use kwargs or defaults
        field_number = field_override.get("number") or kwargs.get("field_number") or "."
        field_type = field_override.get("type") or kwargs.get("field_type") or "String"
        field_description = field_override.get("description") or kwargs.get("field_description") or "unknown"
        field_source = field_override.get("source") or kwargs.get("field_source") or None
        field_version = field_override.get("version") or kwargs.get("field_version") or None

        # Return vcf.parser._Info object with applied overrides
        return vcf.parser._Info(
            field_name,
            field_number,
            field_type,
            field_description,
            field_source,
            field_version,
            self.code_type_map[field_type],
        )
    