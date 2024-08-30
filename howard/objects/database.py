import hashlib
import duckdb
import sqlite3
import vcf
import glob
import os
import io
import pgzip
import shutil

import polars as pl
import pandas as pd
import Bio.bgzf as bgzf
import pyarrow.parquet as pq
import pyarrow as pa
import logging as log

from tempfile import TemporaryDirectory
from typing import Optional, Union

from howard.functions.commons import (
    load_duckdb_extension,
    get_file_format,
    full_path,
    get_file_compressed,
    get_compression_type,
    remove_if_exists,
    concat_and_compress_files,
    DTYPE_LIMIT_AUTO,
    CODE_TYPE_MAP,
)

SEP_TYPE = {
    "vcf": "\t",
    "tsv": "\t",
    "csv": ",",
    "tbl": "|",
    "bed": "\t",
}

DATABASE_TYPE_NEEDED_COLUMNS = {
    "variants": {
        "#CHROM": ["#CHROM", "CHROM", "CHR", "CHROMOSOME"],
        "POS": ["POS"],
        "REF": ["REF"],
        "ALT": ["ALT"],
    },
    "regions": {
        "#CHROM": ["#CHROM", "CHROM", "CHR", "CHROMOSOME"],
        "START": ["START", "POSITIONSTART", "POS"],
        "END": ["END", "POSITIONEND", "POS"],
    },
    "vcf": {
        "#CHROM": ["#CHROM", "CHROM", "CHR", "CHROMOSOME"],
        "POS": ["POS", "POSITION"],
        "ID": ["ID", "IDENTIFIER"],
        "REF": ["REF", "REFERENCE"],
        "ALT": ["ALT", "ALTERNATIVE"],
        "QUAL": ["QUAL", "QUALITY"],
        "FILTER": ["FILTER"],
        "INFO": ["INFO"],
    },
    "bed": {
        "#CHROM": ["#CHROM", "CHROM", "CHR", "CHROMOSOME"],
        "START": ["START", "POSITIONSTART", "POS"],
        "END": ["END", "POSITIONEND", "POS"],
    },
}

DEFAULT_VCF_HEADER = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"]

DEFAULT_VCF_HEADER_DUCKDB_TYPES = {
    "#CHROM": "STRING",
    "POS": "INT",
    "START": "INT",
    "END": "INT",
    "ID": "VARCHAR",
    "REF": "VARCHAR",
    "ALT": "VARCHAR",
    # "QUAL": "INT",
    "FILTER": "VARCHAR",
    "INFO": "VARCHAR",
}


DEFAULT_HEADER_LIST = ["##fileformat=VCFv4.2", "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO"]

FILE_FORMAT_DELIMITERS = {"vcf": "\t", "tsv": "\t", "csv": ",", "tbl": "|", "bed": "\t"}

DUCKDB_EXTENSION_TO_LOAD = ["sqlite_scanner"]


class Database:

    def __init__(
        self,
        database: str = None,
        format: str = None,
        header: str = None,
        header_file: str = None,
        databases_folders: list = None,
        assembly: str = None,
        conn=None,
        conn_config: dict = {},
        table: str = None,
    ) -> None:
        """
        This is an initialization function for a class that sets up a database and header file for use
        in a DuckDB connection.

        :param database: A string representing the name of the database to be used. If None, the default
        database will be used
        :type database: str
        :param format: The `format` parameter is not described in the docstring, so it is unclear what
        it represents
        :type format: str
        :param header: The `header` parameter is a string that represents the name of the header file
        that contains the column names for the database. It is used in conjunction with the `database`
        parameter to set the header for the database. If the `header` parameter is not provided, the
        header will be set to
        :type header: str
        :param header_file: The `header_file` parameter is a string that represents the file path to the
        header file that contains the column names for the database. It is used in the `set_header()`
        method to set the header attribute of the class
        :type header_file: str
        :param databases_folders: A list of folders where the database files are located. This parameter
        is used in the `set_database()` method to search for the database file in the specified folders.
        If the database file is not found in any of the folders, an error is raised
        :type databases_folders: list
        :param assembly: A string representing the name of the assembly to be used. It is used in
        conjunction with the `set_assembly()` method to set the assembly for the DuckDB connection. If
        the `assembly` parameter is not provided, the default assembly will be used
        :type assembly: str
        :param conn: An optional parameter that represents an existing DuckDBPyConnection object. If
        provided, the class will use this connection instead of creating a new one. If not provided, a
        new connection will be created
        :param conn_config: An optional parameter for DuckDBPyConnection object config (see duckdb.connect)
        :param table: The `table` parameter is a string representing the name of the table in the
        database that will be used in the DuckDB connection. It is used in the `set_table()` method to
        set the table attribute of the class. If the `table` parameter is not provided, the default
        table will
        :type table: str
        """

        # Init
        self.database = database
        self.format = format
        self.header = header
        self.header_file = header_file
        self.databases_folders = databases_folders
        self.assembly = assembly
        self.table = table

        # DuckDB connexion
        create_view = False
        if conn:
            self.conn = conn
        elif type(database) == duckdb.DuckDBPyConnection:
            self.conn = database
        else:
            self.conn = duckdb.connect(config=conn_config)
            create_view = True

        # Install duckDB extensions
        load_duckdb_extension(self.conn, DUCKDB_EXTENSION_TO_LOAD)

        # Check attributes
        self.set_format(format=format)
        self.set_assembly(assembly=assembly)
        self.set_databases_folders(databases_folders=databases_folders)
        self.set_header_file(header_file=header_file)
        self.set_database(
            database=database,
            databases_folders=self.get_database_folders(),
            format=self.get_format(),
            assembly=self.get_assembly(),
        )
        self.set_header(
            database=self.get_database(), header=header, header_file=header_file
        )
        if create_view:
            self.create_view(database=database)

    def set_database(
        self,
        database: str,
        databases_folders: list = None,
        format: str = None,
        assembly: str = None,
    ) -> str:
        """
        This function sets the database attribute of an object to a specified database if it exists or
        can be found in a list of folders.

        :param database: A string representing the name of the database to be set
        :type database: str
        :param databases_folders: The `databases_folders` parameter is a list of folders/directories
        where the `find_database` method will search for the specified database. If the database is
        found in any of these folders, it will be set as the current database. If `databases_folders` is
        not provided, the
        :type databases_folders: list
        :param format: The `format` parameter is an optional string representing the format of the
        database to be searched for. If provided, the `find_database` method will search for the
        database only in the specified format. If not provided, the method will search for the database
        in all formats
        :type format: str
        :param assembly: The `assembly` parameter is an optional string representing the name of the
        assembly to which the database belongs. If provided, the `find_database` method will search for
        the database only in the specified assembly. If not provided, the method will search for the
        database in all assemblies
        :type assembly: str
        """

        if database is None:
            self.database = None
        elif type(database) == duckdb.DuckDBPyConnection:
            self.database = database
        elif self.exists(database=database):
            self.database = database
        elif self.find_database(
            database=database,
            databases_folders=databases_folders,
            database_format=format,
            assembly=assembly,
        ):
            self.database = self.find_database(
                database=database,
                databases_folders=databases_folders,
                database_format=format,
                assembly=assembly,
            )

    def set_databases_folders(self, databases_folders: list = None) -> None:
        """
        This function sets the list of folders where databases are located as an attribute of an object.

        :param databases_folders: databases_folders is a list parameter that contains the paths to the
        folders where the databases are stored. The default value of the parameter is a list with a
        single element, which is the current directory (".")
        :type databases_folders: list
        """

        if databases_folders is None:
            databases_folders = ["."]

        self.databases_folders = databases_folders

    def get_database_folders(self) -> list:
        """
        This function returns a list of database folders.
        :return: The method `get_database_folders` is returning a list of database folders. The specific
        list being returned is stored in the instance variable `databases_folders`.
        """

        return self.databases_folders

    def read_header_file(self, header_file: str = None) -> list:
        """
        This function reads the header of a VCF file and returns a list of the header lines.

        :param header_file: The path to the VCF file header that needs to be read
        :type header_file: str
        :return: a list of header lines of a VCF file.
        """

        # Full path
        header_file = full_path(header_file)

        if not header_file:
            return []

        elif not os.path.exists(header_file):
            return []

        elif os.path.isdir(header_file):
            return []

        else:

            header_file_compressed = get_file_compressed(header_file)
            header_compression_type = get_compression_type(header_file)

            if header_file_compressed:
                if header_compression_type in ["bgzip"]:
                    with bgzf.open(header_file, "rt") as header_lines:
                        header_list = []
                        for line in header_lines:
                            if not line.startswith("#"):
                                break
                            header_list.append(line)
                        return header_list
                else:
                    with pgzip.open(header_file, "rt") as header_lines:
                        header_list = []
                        for line in header_lines:
                            if not line.startswith("#"):
                                break
                            header_list.append(line)
                        return header_list
            else:
                with open(header_file, "rt") as header_lines:
                    header_list = []
                    for line in header_lines:
                        if not line.startswith("#"):
                            break
                        header_list.append(line)
                    return header_list

    def get_header_length(self, header_file: str = None) -> int:
        """
        The `get_header_length` function returns the length of a header file, excluding the first line.

        :param header_file: The `header_file` parameter is a string that represents the file path or
        name of the header file. It is an optional parameter, which means it can be omitted when calling
        the `get_header_length` method
        :type header_file: str
        :return: an integer, which represents the length of the header file.
        """

        # No header is 0 header line
        if not header_file:
            header_list = []
        # header_file is a list of header line
        elif isinstance(header_file, list):
            header_list = header_file
        # header_file is a file path
        elif isinstance(header_file, str):
            header_list = self.read_header_file(header_file)
        # Type not available
        else:
            header_list = []

        nb_line = 0
        for line in header_list:
            if not (
                line.startswith("##") or line.startswith("# ") or line.startswith("#\n")
            ):
                break
            nb_line += 1
        return nb_line

    def get_header_file_columns(self, header_file: str = None) -> list:
        """
        The function `get_header_columns` returns the header list of a VCF file.

        :param header_file: The `header_file` parameter is a string that represents the file path of the
        header file. It is an optional parameter and its default value is `None`
        :type header_file: str
        :return: a list of header columns.
        """

        if not header_file or not self.read_header_file(header_file=header_file):
            return []
        else:
            return (
                self.read_header_file(header_file=header_file)[-1].strip().split("\t")
            )

    def get_header_from_list(self, header_list: list = None) -> vcf:
        """
        The function `get_header_from_list` returns a `vcf.Reader` object with a header generated from a
        given list or a default list.

        :param header_list: The `header_list` parameter is a list of strings representing the header
        lines of a VCF (Variant Call Format) file. It is an optional parameter, meaning it can be
        provided as an argument to the function, but if no argument is provided, a default list of
        header lines will be used
        :type header_list: list
        :return: a `vcf.Reader` object.
        """

        if not header_list:
            header_list = DEFAULT_HEADER_LIST

        # Clean lines
        header_list = [line.replace("\n", "") for line in header_list]

        try:
            return vcf.Reader(io.StringIO("\n".join(header_list)))
        except Exception as inst:
            if str(inst).strip():
                log.error(inst)
                raise ValueError(inst)
            else:
                log.warning("Warning in VCF header")
                log.debug(f"header_list={header_list}")
                return None

    def get_header_from_file(self, header_file: str) -> vcf:
        """
        This function returns a VCF header either from a default list or from a file.

        :param header_file: A string representing the file path of a VCF header file. If this parameter
        is not provided or is an empty string, the function will use a default header list
        :type header_file: str
        :return: a VCF object, which is obtained by calling the `get_header_from_list` method with the
        `header_list` as an argument. The `header_list` is either the default header list or the list
        obtained by reading a header file using the `read_header_file` method.
        """

        if not header_file:
            header_list = DEFAULT_HEADER_LIST
        else:
            header_list = self.read_header_file(header_file=header_file)

        return self.get_header_from_list(header_list)

    def find_header_file(self, database: str = None) -> Optional[str]:
        """
        This function finds the header file for a given database in various formats.

        :param database: The `database` parameter is a string that represents the path to a database
        file. If this parameter is not provided, the `get_database()` method is called to retrieve the
        path to the database file
        :type database: str
        :return: the path to the header file for a given database. If the header is in a separate file,
        it returns the path to that file. If the header is within the database file itself, it returns
        the path to the database file. If the database or its format cannot be determined, it returns
        None.
        """

        if not database:
            database = self.get_database()

        if not database:
            return None

        # database format
        database_format = self.get_format(database=database)

        # extra header file
        database_header_file = None

        # header in extra header file
        database = full_path(database)
        if os.path.exists(f"{database}.hdr"):
            database_header_file = f"{database}.hdr"

        # header within file
        elif database_format in ["vcf", "tsv", "csv", "tbl", "bed"]:
            database_header_file = database

        return database_header_file

    def get_header(
        self,
        database: str = None,
        header_file: str = None,
        header_list: list = None,
        sql_query: str = None,
    ) -> vcf:
        """
        The `get_header` function in Python returns the header of a VCF file from a file, a list, or the
        object itself based on specified conditions.

        :param database: The `database` parameter in the `get_header` function represents a string that
        specifies the database from which the header information should be retrieved or used. It is used
        in various parts of the function to determine how to construct the header of the VCF file
        :type database: str
        :param header_file: The `header_file` parameter in the `get_header` function is a string
        representing the path to a file containing the header information for a VCF file. This parameter
        allows you to specify a file from which the function will read the header information
        :type header_file: str
        :param header_list: The `header_list` parameter in the `get_header` function is a list
        containing the header lines of a VCF file. If provided, the function will construct the header
        from this list using the `get_header_from_list` method. If `header_list` is not provided, the
        function will
        :type header_list: list
        :param sql_query: The `sql_query` parameter in the `get_header` function is a string
        representing an SQL query that can be used to retrieve data from a database. This parameter is
        used in the function to help construct the header of a VCF file based on the query results or
        other conditions specified in the function
        :type sql_query: str
        :return: The `get_header` function returns the header of a VCF file based on different
        conditions:
        """

        if self.header:
            # Construct header from a given header object
            return self.header
        elif header_file:
            # Construct header from a given header file
            return self.get_header_from_file(header_file)
        elif header_list:
            # Construct header from a given header list
            return self.get_header_from_list(header_list)
        elif self.find_header_file(database=database):
            # Construct header from header file found
            return self.get_header_from_file(self.find_header_file(database=database))
        elif self.get_extra_columns(database=database, sql_query=sql_query):
            # Construct header from a list of columns
            return self.get_header_from_columns(
                database=database,
                header_columns=self.get_extra_columns(
                    database=database, sql_query=sql_query
                ),
                sql_query=sql_query,
            )
        # elif "INFO" in self.get_columns(database=database):
        #     # Construct header from annotation in INFO column (in case of VCF database without header, e.g. in TSV)
        #     # TODO
        #     return None
        else:
            # No header
            return None

    def get_header_from_columns(
        self, database: str = None, header_columns: list = [], sql_query: str = None
    ) -> object:
        """
        The function `get_header_from_columns` generates a VCF header based on database columns and adds
        custom annotations to it.

        :param database: The `database` parameter is a string that represents the name of a database. It
        is an optional parameter, and if not provided, the `get_database()` method is called to retrieve
        the default database. This parameter specifies the database from which the columns will be used
        to generate the VCF header
        :type database: str
        :param header_columns: The `header_columns` parameter is a list of column names that will be
        used to generate header information for a VCF file. If no `header_columns` are provided, the
        function will attempt to automatically detect the columns to use based on the database being
        used
        :type header_columns: list
        :param sql_query: The `sql_query` parameter in the `get_header_from_columns` function is used to
        specify a SQL query that will be executed to retrieve column information from the database. This
        query can be customized to fetch specific columns or data based on the requirements of the VCF
        header generation process. If provided,
        :type sql_query: str
        :return: The function `get_header_from_columns` returns a VCF header object that includes
        information about the columns in a database and their data types. The header object is created
        based on the input parameters, including the database name and a list of header columns.
        """

        if not database:
            database = self.get_database()

        database_header = vcf.Reader(io.StringIO("\n".join(DEFAULT_HEADER_LIST)))

        if not header_columns:
            header_columns = self.get_extra_columns(
                database=database, sql_query=sql_query
            )

        database_basename = self.get_database_basename(database=database) or "unknown"

        # Columns query for auto detection of stype

        # Attach if need
        if self.get_sql_database_attach(database=database, output="attach"):
            self.query(
                query=self.get_sql_database_attach(database=database, output="attach"),
            )

        # database columns
        database_query_columns = None
        if sql_query:
            database_query_columns = self.query(query=sql_query)
        if not database_query_columns:
            database_query_columns_sql = f""" SELECT * FROM {self.get_sql_database_link(database=database)} LIMIT {DTYPE_LIMIT_AUTO} """
            database_query_columns = self.query(query=database_query_columns_sql)

        # Remove specific VCF column if is a VCF type
        if (
            self.get_type_from_columns(
                database_columns=self.get_columns(
                    database=database, sql_query=sql_query
                ),
                check_database_type="vcf",
            )
            == "vcf"
        ):
            header_columns = [x for x in header_columns if x not in DEFAULT_VCF_HEADER]

        # List all columns to add into header
        for header_column in header_columns:

            # Header info type
            header_info_type = "String"
            header_column_df = database_query_columns.df()[header_column]
            header_column_df_dtype = header_column_df.dtype
            if header_column_df_dtype == object:
                if pd.to_numeric(header_column_df, errors="coerce").notnull().all():
                    header_info_type = "Float"
            else:
                header_info_type = "Integer"

            # Header info
            header_info_name = header_column
            header_info_number = "."
            header_info_description = f"{header_column} annotation"
            header_info_source = database_basename
            header_info_version = "unknown"
            database_header.infos[header_column] = vcf.parser._Info(
                header_info_name,
                header_info_number,
                header_info_type,
                header_info_description,
                header_info_source,
                header_info_version,
                CODE_TYPE_MAP[header_info_type],
            )

        # Detach if need
        if self.get_sql_database_attach(database=database, output="detach"):
            self.query(
                query=self.get_sql_database_attach(database=database, output="detach"),
            )

        return database_header

    def query(self, query: str = None) -> object:
        """
        This is a Python function that takes in a database and query string as parameters and returns
        the result of the query on the database.

        :param query: The query parameter is a string that represents the SQL query that needs to be
        executed on the database. It can be any valid SQL statement such as SELECT, INSERT, UPDATE,
        DELETE, etc
        :type query: str
        :return: If a query is provided, the method returns the result of the query executed on the
        database. If no query is provided, the method returns None.
        """

        if query:
            return self.get_conn().query(query)
        else:
            return None

    def set_header(
        self, database: str = None, header: vcf = None, header_file: str = None
    ) -> None:
        """
        This function sets the header of a database based on a provided header file or the database
        format.

        :param database: A string representing the name or path of a database file. If not provided, the
        method will attempt to get the database name from the object's attributes
        :type database: str
        :param header: `header` is a variable of type `vcf` (presumably representing a VCF header) that
        can be provided as an argument to the `set_header` method to set the header attribute of the
        object. If `header` is provided, the `header_file` parameter is ignored
        :type header: vcf
        :param header_file: A string representing the file path of a header file. If provided, the
        function will use this header file to set the header attribute of the object
        :type header_file: str
        """

        # Full path
        database = full_path(database)
        header_file = full_path(header_file)

        if header:

            self.header = header
            self.header_file = None

        else:

            if header_file and os.path.exists(header_file):

                # header provided
                self.header = self.get_header(header_file=header_file)
                self.header_file = header_file

            else:

                # default no haeder file
                self.header_file = None

                if not database:
                    database = self.get_database()

                if database:

                    # database format
                    database_format = self.get_format(database=database)

                    # extra header file
                    database_header_file = f"{database}.hdr"

                    # header in extra header file
                    if os.path.exists(database_header_file):
                        self.header = self.get_header(header_file=database_header_file)
                        self.header_file = database_header_file

                    # header within file
                    elif database_format in ["vcf", "tsv", "csv", "tbl", "bed"]:
                        self.header = self.get_header(header_file=database)
                        if self.header:
                            self.header_file = self.get_database()

                    # Not header
                    else:
                        self.header = self.get_header()

                else:

                    self.header = None

    def set_header_file(self, header_file: str = None) -> None:
        """
        This function sets the header file attribute of an object to the value passed as an argument.

        :param header_file: The parameter `header_file` is a string that represents the name or path of
        a header file. This method sets the `header_file` attribute of an object to the value passed as
        an argument. If no argument is passed, the `header_file` attribute remains unchanged
        :type header_file: str
        """

        self.header_file = header_file

    def get_header_columns_from_database(self, database: str = None) -> Optional[list]:
        """
        The function `get_header_columns_from_database` retrieves the column names from a specified
        database table.

        :param database: The `database` parameter is a string that represents the name of the database
        from which you want to retrieve the header columns. If no database is provided, the method will
        use the `get_database()` method to retrieve the default database
        :type database: str
        :return: a list of column names from the specified database.
        """

        if not database:
            database = self.get_database()

        sql_from = self.get_sql_from(database=database)
        if sql_from:
            sql_query = f"SELECT * FROM {sql_from} LIMIT 0"
            try:
                columns_list = list(self.conn.query(sql_query).columns)
            except ValueError:
                columns_list = None
            return columns_list
        else:
            return None

    def get_header_file(
        self,
        header_file: str = None,
        remove_header_line: bool = False,
        replace_header_line: list = None,
        force: bool = False,
        sql_query: str = None,
    ) -> str:
        """
        The function `get_header_file` generates a VCF header file based on specified parameters or a
        default header if needed.

        :param header_file: The `header_file` parameter is a string representing the file path and name
        of the header file. If set to `None`, the default header file path and name will be used
        :type header_file: str
        :param remove_header_line: The `remove_header_line` parameter is a boolean parameter that
        determines whether to remove the `#CHROM` line from the header file. If set to True, the line
        will be removed; otherwise, it will remain in the header file. By default, this parameter is set
        to False, meaning, defaults to False
        :type remove_header_line: bool (optional)
        :param replace_header_line: The `replace_header_line` parameter is a list of columns that can be
        used to replace the header line in the generated header file. For example, if you provide
        `['#CHROM', 'POS', 'ID']`, these columns will be used as the header line in the generated file
        instead
        :type replace_header_line: list
        :param force: The `force` parameter in the `get_header_file` function is a boolean parameter
        that determines whether to force the generation of a header file even if a header file already
        exists. If `force` is set to `True`, the function will replace the existing header file with a
        new one. If, defaults to False
        :type force: bool (optional)
        :param sql_query: The `sql_query` parameter in the `get_header_file` function is used to specify
        a SQL query that can be used to retrieve header information from a database. This query can be
        passed to the function to customize the header generation process based on the query results
        :type sql_query: str
        :return: The function `get_header_file` returns a string which is the name of the header file
        that was generated or None if no header file was generated.
        """

        # Full path
        header_file = full_path(header_file)

        if not header_file:
            header_file = self.header_file

        if header_file and os.path.exists(header_file) and not force:
            return header_file

        with TemporaryDirectory() as tmp_dir:
            header_file_tmp = os.path.join(tmp_dir, "header.hdr")

            if header_file:
                if header_file != self.get_database():
                    if self.get_header(sql_query=sql_query):
                        header = self.get_header(sql_query=sql_query)
                    else:
                        header = self.get_header(
                            header_list=DEFAULT_HEADER_LIST, sql_query=sql_query
                        )
                    # Generate header file
                    if not os.path.exists(header_file_tmp):
                        f = open(header_file_tmp, "w")
                        vcf.Writer(f, header)
                        f.close()
            else:
                # No header generated
                header_file = None

            if (
                header_file
                and header_file != self.get_database()
                and (not os.path.exists(header_file) or force)
            ):
                with open(header_file_tmp, "r") as file:
                    lines = file.readlines()
                with open(header_file, "w") as file:
                    for line in lines:
                        if line.startswith("#CHROM"):
                            if not remove_header_line:
                                if replace_header_line:
                                    file.write("\t".join(replace_header_line) + "\n")
                                else:
                                    file.write(line)
                        else:
                            file.write(line)

        return header_file

    def set_assembly(self, assembly: str = None) -> None:
        """
        This is a function that sets the assembly attribute of an object to a given string value.

        :param assembly: The assembly parameter is a string that represents the name or type of assembly
        that the object belongs to. This method sets the assembly attribute of the object to the value
        passed in as the assembly parameter. If no value is passed in, the assembly attribute remains
        unchanged. The method returns the updated value of the
        :type assembly: str
        """

        self.assembly = assembly

    def get_assembly(self) -> str:
        """
        This function returns the assembly attribute of an object if it exists, otherwise it returns
        None.
        :return: If `self.assembly` is not `None`, then it returns the value of `self.assembly`.
        Otherwise, it returns `None`.
        """

        return self.assembly

    def find_database(
        self,
        database: str = None,
        databases_folders: list = None,
        database_format: str = None,
        assembly: str = None,
    ) -> Optional[str]:
        """
        This function finds a database file in a specified folder or the current directory.

        :param database: The name of the database to be searched for. If not provided, it will call the
        `get_database()` method to get the name of the database. It is a string type parameter
        :type database: str
        :param databases_folders: A list of folders where the function should look for the database
        file. If this parameter is not provided, the function will look for the database file in the
        current directory
        :type databases_folders: list
        :param format: The file format of the database file. It is an optional parameter and if not
        provided, the function will call the `get_format()` method to get the format
        :type format: str
        :param assembly: `assembly` is an optional parameter that represents the name of a subfolder
        where the function should look for the database file. If provided, the function will search for
        the database file in the specified subfolder within each of the `databases_folders`. If not
        provided, the function will only search for
        :type assembly: str
        :return: a string that represents the path to the database file. If the database is not found or
        if no database is specified, it returns None.
        """

        if not database:
            database = self.get_database()

        if not database_format:
            database_format = self.get_format()

        if not assembly:
            assembly = self.get_assembly()

        if not database:
            return None

        elif self.exists(database=database):
            return database

        else:

            if not databases_folders:
                databases_folders = ["."]

            database_file = None

            for format_extension in [
                "",
                f".{database_format}",
                f".{database_format}.gz",
            ]:

                # find in subfolder assemby
                if assembly:
                    for databases_folder in databases_folders:

                        databases_folder = full_path(databases_folder)
                        database_file_check = f"{databases_folder}/{assembly}/{database}{format_extension}"

                        # Log
                        log.debug("Annotation file check: " + database_file_check)

                        # In folder

                        if os.path.exists(database_file_check):
                            database_file = database_file_check

                        # database found
                        if database_file:
                            break

                # find within folder
                if not database_file:

                    # find in folders
                    for databases_folder in databases_folders:

                        databases_folder = full_path(databases_folder)
                        database_file_check = (
                            f"{databases_folder}/{database}{format_extension}"
                        )

                        # Log
                        log.debug("Annotation file check: " + database_file_check)

                        # In folder
                        if os.path.exists(database_file_check):
                            database_file = database_file_check

                        # database found
                        if database_file:
                            break

                # database found
                if database_file:
                    break

            return database_file

    def get_database(self) -> str:
        """
        This function returns the database name as a string.
        :return: The `get_database` method is returning the `database` attribute of the object. The
        return type is a string.
        """

        return self.database

    def get_database_basename(self, database: str = None) -> Optional[str]:
        """
        This function returns the basename of a database file.

        :param database: The parameter `database` is a string that represents the name of a database. If
        it is not provided, the method will use the `get_database()` method to retrieve the current
        database
        :type database: str
        :return: a string which is the basename of the database file. If the database parameter is not
        provided, it gets the current database using the `get_database()` method. If the database
        exists, it returns the basename of the database file using the `os.path.basename()` method. If
        the database does not exist, it returns `None`.
        """

        if not database:
            database = self.get_database()

        if type(database) == duckdb.DuckDBPyConnection:
            return None
        elif database:
            return os.path.basename(database)
        else:
            return None

    def get_database_dirname(self, database: str = None) -> Optional[str]:
        """
        This function returns the directory name of a given database or the current database if none is
        specified.

        :param database: The parameter `database` is a string that represents the path to a database
        file. If it is not provided, the method will call `self.get_database()` to retrieve the path to
        the default database
        :type database: str
        :return: a string that represents the directory name of the specified database file. If no
        database file is specified, it will use the default database file and return its directory name.
        If there is no database file, it will return None.
        """

        if not database:
            database = self.get_database()

        if database:
            return os.path.dirname(database)
        else:
            return None

    def exists(self, database: str = None) -> bool:
        """
        This function checks if a database exists in the specified path or in the default path.

        :param database: The `database` parameter is a string that represents the name or path of a
        database file. If it is not provided, the method will call the `get_database()` method to
        retrieve the default database name/path
        :type database: str
        :return: a boolean value indicating whether the specified database exists or not. If the
        `database` parameter is not provided, it gets the current database using the `get_database()`
        method and checks if it exists using the `os.path.exists()` function.
        """

        # Full path
        database = full_path(database)

        if not database:
            database = self.get_database()

        return database and (
            type(database) == duckdb.DuckDBPyConnection or os.path.exists(database)
        )

    def set_format(self, format: str = None) -> str:
        """
        This is a method in a Python class that sets a format attribute to a specified string.

        :param format: The format parameter is a string that specifies the desired format for the data.
        It is an optional parameter, meaning that if it is not provided, the format attribute of the
        object will not be changed. The function returns a string indicating the current format of the
        object
        :type format: str
        """

        self.format = format

    def get_format(self, database: str = None) -> str:
        """
        This Python function returns the file format of a given database or the current database if none
        is provided.
        Format database:
            - parquet
            - duckdb
            - sqlite
            - vcf
            - csv

        :param database: The `database` parameter is a string that represents the type of database. It
        is an optional parameter and if not provided, the function will call the `get_database()` method
        to retrieve the database type
        :type database: str
        :return: a string that represents the type of database. The type of database can be one of the
        following: "parquet", "duckdb", "sqlite", "vcf", or "csv". The specific type of database is determined by
        the input parameter `database`, which is either passed as an argument to the function or
        obtained by calling the `get_database()` method. The `get_file_format
        """

        if self.format:
            return self.format

        if not database:
            database = self.get_database()

        if type(database) == duckdb.DuckDBPyConnection:
            return "conn"
        else:
            return get_file_format(database)

    def get_type(self, database: str = None, sql_query: str = None) -> Optional[str]:
        """
        The `get_type` function determines the type of a database (variants VCF-like or regions
        BED-like) based on its columns and format.

        :param database: The `database` parameter in the `get_type` function is a string representing
        the name of a database. If this parameter is not provided when calling the function, it will
        attempt to retrieve the database name using the `get_database()` method. This parameter is used
        to specify the database for which you
        :type database: str
        :param sql_query: The `sql_query` parameter in the `get_type` function is used to pass an SQL
        query as a string. This query can be used to filter or manipulate the data before determining
        the type of the database based on its columns. If provided, the function will use this SQL query
        to fetch the
        :type sql_query: str
        :return: The `get_type` function returns a string that represents the type of the database,
        which can be either "variants" (VCF-like) or "regions" (BED-like). If the database is not found
        or does not exist, the function returns None.
        """

        if not database:
            database = self.get_database()

        if database and self.exists(database):
            database_columns = self.get_columns(
                database, table=self.get_database_table(database), sql_query=sql_query
            )
            database_type = self.get_type_from_columns(database_columns)
            if database_type:
                return database_type
            else:
                database_format = self.get_format()
                if database_format in ["vcf"]:
                    return "vcf"
                elif database_format in ["bed"]:
                    return "regions"
                else:
                    return None
        else:
            return None

    def get_database_tables(self, database: str = None) -> Union[str, list, None]:
        """
        This function retrieves a list of tables in a specified database using the DuckDB format.

        :param database: The name of the database for which you want to retrieve the list of tables. If
        no database name is provided, it will use the default database
        :type database: str
        :return: a list of tables in the specified database, or None if the database does not exist or
        the format is not supported.
        """

        if not database:
            database = self.get_database()
        if database and self.exists(database):
            database_format = self.get_format(database)
            if database_format in ["conn"]:
                database_conn = database
                database_tables = list(database_conn.query("SHOW TABLES;").df()["name"])
                return database_tables
            elif database_format in ["duckdb"]:
                database_conn = duckdb.connect(database)
                database_tables = list(database_conn.query("SHOW TABLES;").df()["name"])
                database_conn.close()
                return database_tables
            elif database_format in ["sqlite"]:
                database_conn = sqlite3.connect(database)
                sql_query = "SELECT name FROM sqlite_master WHERE type='table'"
                database_tables = list(
                    pd.read_sql_query(sql_query, database_conn)["name"]
                )
                return database_tables
            else:
                return None
        else:
            return None

    def get_database_table(self, database: str = None) -> Optional[str]:
        """
        This function returns the name of a table in a specified database if it exists and is in a
        supported format.

        :param database: The name of the database to retrieve the table from. If None, it will use the
        default database
        :type database: str
        :return: a string that represents the name of a table in a database, or None if no suitable
        table is found.
        """

        if self.table:
            return self.table

        if not database:
            database = self.get_database()

        if database and self.exists(database):

            database_format = self.get_format(database)
            if database_format in ["duckdb", "sqlite", "conn"]:
                for database_table in self.get_database_tables(database=database):
                    database_columns = self.get_columns(database, table=database_table)
                    database_type = self.get_type_from_columns(database_columns)
                    if database_type:

                        return database_table
                return None
            else:
                return None
        else:
            return None

    def get_type_from_columns(
        self, database_columns: list = [], check_database_type: str = None
    ) -> Optional[str]:
        """
        This function returns the type of a database based on the provided list of columns.

        :param database_columns: a list of column names in a database table
        :type database_columns: list
        :param check_database_type: A database type to check for. If not provided, it defaults
        to all database types defined in the constant `DATABASE_TYPE_NEEDED_COLUMNS`
        :type check_database_type: str
        :return: a string that represents the type of database based on the provided list of columns. If
        the needed columns for a specific database type are not found in the provided list, the function
        returns None.
        """

        if check_database_type:
            check_database_type = [check_database_type]
        else:
            check_database_type = DATABASE_TYPE_NEEDED_COLUMNS.keys()

        for database_type in check_database_type:
            needed_columns = self.get_needed_columns(database_columns, database_type)
            all_needed_columns_found = True
            for col in needed_columns:
                all_needed_columns_found = all_needed_columns_found and (
                    needed_columns[col] is not None
                )
            if all_needed_columns_found:
                return database_type

        return None

    def get_needed_columns(
        self, database_columns: list = [], database_type: str = None
    ) -> dict:
        """
        This function takes a list of database columns and a type, and returns a dictionary of needed
        columns and their corresponding values found in the database columns.

        :param database_columns: A list of column names in a database table
        :type database_columns: list
        :param type: The type of database being used. It is used to determine which columns are needed
        for the specific database type
        :type type: str
        :return: a dictionary containing the columns that are needed for a specific database type, along
        with their corresponding column names in the actual database. The function takes in a list of
        database columns and a database type as input, and uses the `DATABASE_TYPE_NEEDED_COLUMNS`
        dictionary to determine which columns are needed for the specified database type. It then
        searches through the list of database columns to find the
        """

        needed_columns = DATABASE_TYPE_NEEDED_COLUMNS.get(database_type)
        variants_columns_found = {}

        if needed_columns:
            for needed_col in needed_columns:
                variants_columns_found[needed_col] = None
                for possible_col in needed_columns[needed_col]:
                    if database_columns:
                        for existing_column in database_columns:
                            if possible_col.upper() == existing_column.upper():
                                variants_columns_found[needed_col] = existing_column
                                break

        return variants_columns_found

    def get_sql_from(
        self, database: str = None, header_file: str = None, sample_size: int = 20480
    ) -> str:
        """
        This function returns a SQL query string based on the input database format.

        :param database: The parameter "database" is a string that represents the name or path of the
        database that the function will read from. If no value is provided for this parameter, the
        function will call the "get_database()" method to retrieve the default database
        :type database: str
        :return: a string that represents a SQL query to read data from a database file. The specific
        SQL query returned depends on the format of the database file, which is determined by the
        `get_format()` method. The SQL query returned will be in the form of a function call to one of
        the following functions: `read_parquet()`, `read_csv()`, `read_json()`,
        """

        if not database:
            database = self.get_database()

        database_format = self.get_format(database=database)

        sql_from = None

        # Connexion
        if type(database) == duckdb.DuckDBPyConnection:
            sql_from = self.get_database_table(database=database)

        # Parquet
        elif database_format in ["parquet"]:

            # Check for partition
            if os.path.isdir(database):
                list_of_parquet = glob.glob(
                    os.path.join(database, "**/*parquet"), recursive=True
                )
                if list_of_parquet:
                    list_of_parquet_level_path = "*/" * (
                        list_of_parquet[0].replace(database, "").count("/") - 1
                    )
                    sql_from = f"read_parquet('{database}/{list_of_parquet_level_path}*parquet', hive_partitioning=1)"
                else:
                    log.error(
                        f"Input file '{database}' not a compatible partitionned parquet folder"
                    )
                    raise ValueError(
                        f"Input file '{database}' not a compatible partitionned parquet folder"
                    )

            # No partition
            else:
                sql_from = f"read_parquet('{database}')"

        # CSV
        elif database_format in ["vcf", "tsv", "csv", "tbl", "bed"]:

            # Delimiter
            delimiter = SEP_TYPE.get(database_format, "\t")

            # Check infos from database
            table_columns = self.get_table_columns_from_file(
                database=database, header_file=header_file
            )
            header_length = self.get_header_length(header_file=database)
            file_compressed = get_file_compressed(database)
            if file_compressed:
                database_compression = "gzip"
            else:
                database_compression = "none"
            hive_partitioning = 0
            database_ref = database

            # Check for partition
            if os.path.isdir(database):
                list_of_parquet = glob.glob(
                    os.path.join(database, "**/*csv"), recursive=True
                )
                if list_of_parquet:
                    database_compression = "gzip"
                    hive_partitioning = 1
                    list_of_parquet_level_path = "*/" * (
                        list_of_parquet[0].replace(database, "").count("/") - 1
                    )
                    database_ref = f"{database}/{list_of_parquet_level_path}*csv"
                else:
                    log.error(
                        f"Input file '{database}' not a compatible partitionned parquet folder"
                    )
                    raise ValueError(
                        f"Input file '{database}' not a compatible partitionned parquet folder"
                    )

            # Query number of columns detected by duckdb
            query_nb_columns_detected_by_duckdb = f"""
                    SELECT *
                    FROM read_csv('{database_ref}', auto_detect=True, hive_partitioning={hive_partitioning}, compression='{database_compression}', skip={header_length}, delim='{delimiter}', sample_size={sample_size})
                    LIMIT 0
                """
            try:
                nb_columns_detected_by_duckdb = len(
                    self.conn.query(query_nb_columns_detected_by_duckdb).columns
                )
            except ValueError:
                nb_columns_detected_by_duckdb = 0

            # Check table columns
            if not table_columns or (
                nb_columns_detected_by_duckdb != len(table_columns)
            ):
                # Check columns from header
                table_columns = self.get_table_columns_from_format(database=database)

            # If table columns
            if table_columns:

                # Detect and force column type
                table_columns_types_list = []
                for table_column in table_columns:
                    if table_column in DEFAULT_VCF_HEADER_DUCKDB_TYPES:
                        table_columns_types_list.append(
                            f"'{table_column}': {DEFAULT_VCF_HEADER_DUCKDB_TYPES.get(table_column)}"
                        )

                # Create duckdb read_csv types option
                if table_columns_types_list:
                    table_columns_types_list_join_option = (
                        ", types={" + ", ".join(table_columns_types_list) + "}"
                    )
                else:
                    table_columns_types_list_join_option = ""

                # SQL form
                sql_from = f"""read_csv('{database_ref}', names={table_columns}{table_columns_types_list_join_option}, auto_detect=True, compression='{database_compression}', skip={header_length}, delim='{delimiter}', hive_partitioning={hive_partitioning}, sample_size={sample_size})"""

            else:
                sql_from = f"read_csv('{database_ref}', auto_detect=True, compression='{database_compression}', skip={header_length}, delim='{delimiter}', hive_partitioning={hive_partitioning}, sample_size={sample_size})"

        # JSON
        elif database_format in ["json"]:
            sql_from = (
                f"read_json('{database}', auto_detect=True, sample_size={sample_size})"
            )

        # DuckDB
        elif database_format in ["duckdb"]:
            sql_from = f"'{database}'"

        # SQLite
        elif database_format in ["sqlite"]:
            database_table = self.get_database_table(database=database)
            sql_from = f"(SELECT * FROM sqlite_scan('{database}', '{database_table}'))"

        return sql_from

    def get_sql_database_attach(
        self,
        database: str = None,
        output: str = "query",
    ) -> Optional[str]:
        """
        This function returns a SQL query to attach or detach a database based on the specified format
        and output.

        :param database: The name of the database to attach. If not provided, it will try to get the
        default database from the connection
        :type database: str
        :param output: The "output" parameter is a string that specifies the desired output of the
        function. It can take on the following values:, defaults to query
        :type output: str (optional)
        :return: a string that represents a SQL query to attach a database to a DuckDB or SQLite
        database engine. The specific output depends on the value of the `output` parameter, which can
        be "query" (default), "attach", "detach", or "name". If `output` is "query" or "attach", the
        function returns a SQL query to attach the specified database.
        """

        if not database:
            database = self.get_database()

        if not database:
            return None

        database_format = self.get_format(database=database)

        database_attach = None

        if database_format in ["duckdb", "sqlite"]:
            database_name = "database_" + hashlib.sha1(database.encode()).hexdigest()
            database_options = []
            if database_format in ["sqlite"]:
                database_options.append("TYPE SQLITE")
            if database_options:
                database_options_sql = "(" + ",".join(database_options) + ")"
            else:
                database_options_sql = ""

            if output in ["query", "attach"]:
                database_attach = f"ATTACH DATABASE '{database}' AS {database_name} {database_options_sql}"
            elif output in ["detach"]:
                database_attach = f"DETACH {database_name}"
            elif output == "name":
                database_attach = database_name

        return database_attach

    def get_sql_database_link(self, database: str = None) -> str:
        """
        This function returns a SQL database link based on the provided database name or the default
        database.

        :param database: The `database` parameter is a string that represents the name of the database.
        If it is not provided, the method will call the `get_database()` method to retrieve the default
        database
        :type database: str
        :return: a SQL database link as a string. If a database name is provided as an argument, it will
        use that database to construct the link. Otherwise, it will use the default database obtained
        from `self.get_database()`. The link is constructed using the `sql_from` and `sql_table`
        obtained from other methods, and the final link is returned as a string. If the
        """

        if not database:
            database = self.get_database()

        sql_database_link = None

        sql_from = self.get_sql_from(database=database)

        if sql_from:

            database_attach_name = self.get_sql_database_attach(
                database=database, output="name"
            )
            if database_attach_name:
                sql_table = self.get_database_table(database=database)
                sql_from = f"{database_attach_name}.{sql_table}"

            sql_database_link = f"(SELECT * FROM {sql_from})"

        return sql_database_link

    def create_view(
        self, database: str = None, view_name: str = "variants"
    ) -> Optional[str]:
        """
        The `create_view` function creates a view in a specified database or the default database, using
        a SQL database link.

        :param database: The `database` parameter is a string that represents the name of the database.
        If no value is provided, it will use the value returned by the `get_database()` method
        :type database: str
        :param view_name: The `view_name` parameter is a string that specifies the name of the view that
        will be created in the database, defaults to variants
        :type view_name: str (optional)
        :return: the name of the created view.
        """

        if not database:
            database = self.get_database()

        # database link
        sql_database_link = self.get_sql_database_link(database=database)

        # database format
        database_format = self.get_format(database=database)

        # Create view
        if sql_database_link and database_format not in ["duckdb"]:
            try:
                self.get_conn().execute(
                    f"CREATE OR REPLACE VIEW {view_name} AS {sql_database_link}"
                )
                self.view_name = view_name
                return view_name
            except ValueError:
                return None
        else:
            self.view_name = None
            return None

    def get_view(self, database: str = None, create_view: str = None) -> Optional[str]:
        """
        The `get_view` function returns the name of a view in a database, or creates a new view if
        specified.

        :param database: The `database` parameter is a string that represents the name of the database.
        It is an optional parameter and if not provided, the method `get_database()` is called to
        retrieve the database name
        :type database: str
        :param create_view: The `create_view` parameter is a string that represents the name of the view
        that you want to create. If this parameter is provided, the `get_view` method will call the
        `create_view` method and pass the `database` and `view_name` parameters to it
        :type create_view: str
        :return: The method `get_view` returns a string.
        """

        if not database:
            database = self.get_database()

        if create_view:
            return self.create_view(database=database, view_name=create_view)
        elif self.view_name:
            return self.view_name
        else:
            return None

    def is_compressed(self, database: str = None) -> bool:
        """
        This Python function checks if a given file is compressed and returns the format of the
        compression.

        :param database: The `database` parameter is a string that represents the path or name of the
        input file that needs to be checked for compression. If no value is provided for `database`, the
        method calls `get_database()` to retrieve the default database file
        :type database: str
        :return: The function `is_compressed` returns a boolean value indicating whether the input file
        is compressed or not. The function calls another function `get_file_compressed` to determine the
        compression format of the file.
        """

        if not database:
            database = self.get_database()

        if type(database) == duckdb.DuckDBPyConnection:
            return False
        else:
            return get_file_compressed(database)

    def get_header_infos_list(self, database: str = None) -> list:
        """
        This function returns a list of header information for a given database or the current database
        if none is specified.

        :param database: The `database` parameter is a string that represents the name of the database
        from which the header information is to be retrieved. If no database name is provided, the
        method will use the default database name obtained from the `get_database()` method
        :type database: str
        :return: A list of header information from a database, or an empty list if the database header
        is not available.
        """

        if not database:
            database = self.get_database()

        # Database header
        database_header = self.get_header(database=database)

        # Init
        database_header_infos_list = []

        if database_header:
            database_header_infos_list = list(database_header.infos)

        return database_header_infos_list

    def find_column(
        self,
        database: str = None,
        table: str = None,
        column: str = "INFO",
        prefixes: list = ["INFO/"],
        database_columns: list = None,
    ) -> str:
        """
        The `find_column` function searches for a specific column in a database table, with the option
        to search for a column with a specific prefix or within the INFO column header.

        :param database: The name of the database to search for the column in. If not provided, it will
        use the current database that the code is connected to
        :type database: str
        :param table: The "table" parameter is the name of the table in the database where the column is
        located
        :type table: str
        :param column: The "column" parameter is a string that represents the name of the column to
        search for in the database table. By default, it is set to "INFO", but you can change it to
        search for a specific column name, defaults to INFO
        :type column: str (optional)
        :param prefixes: The `prefixes` parameter is a list of strings that are used to search for a
        column with a specific prefix in the database. For example, if the prefixes list contains "DP/",
        the function will search for a column named "DP/INFO" in addition to the default "INFO" column
        :type prefixes: list
        :param database_columns: The `database_columns` parameter is a list that contains the names of
        all the columns in a specific database table. It is used to check if a specific column exists in
        the database. If the `database_columns` parameter is not provided, the function will call the
        `get_columns` method to retrieve
        :type database_columns: list
        :return: a string that represents the name of the column found in the database, based on the
        input parameters. If the column is found, it returns the column name. If the column is not
        found, it returns None.
        """

        if not database:
            database = self.get_database()

        # Database columns
        if not database_columns:
            database_columns = self.get_columns(database=database, table=table)

        # Init
        column_found = None

        # Column exists
        if column in database_columns:
            column_found = column

        # Column with prefix
        elif prefixes:
            for prefix in prefixes:
                if prefix + column in database_columns:
                    column_found = prefix + column
                    break

        # Column in INFO column (test if in header)
        if not column_found and "INFO" in database_columns:
            database_header_infos = self.get_header_infos_list(database=database)
            if column in database_header_infos:
                column_found = "INFO"

        return column_found

    def map_columns(
        self,
        database: str = None,
        table: str = None,
        columns: list = [],
        prefixes: list = ["INFO/"],
    ) -> dict:
        """
        The `map_columns` function maps input columns to their corresponding columns in a specified
        database table, using specified prefixes to filter the columns.

        :param database: The name of the database to search for columns in. If no database is specified,
        the method will use the default database set in the connection
        :type database: str
        :param table: The `table` parameter is the name of the table in the database that you want to
        map the columns for
        :type table: str
        :param columns: A list of column names that you want to map to their corresponding column names
        in the database
        :type columns: list
        :param prefixes: The `prefixes` parameter is a list of strings that are used to filter the
        columns that are searched for. Only columns that start with one of the prefixes in the list will
        be considered. In the code above, the default value for `prefixes` is `["INFO/"]`
        :type prefixes: list
        :return: a dictionary that maps the input columns to their corresponding columns found in the
        specified database and table, with the specified prefixes.
        """

        if not database:
            database = self.get_database()

        # Init
        columns_mapping = {}

        # database_columns
        database_columns = self.get_columns(database=database, table=table)

        for column in columns:
            column_found = self.find_column(
                database=database,
                table=table,
                column=column,
                prefixes=prefixes,
                database_columns=database_columns,
            )
            columns_mapping[column] = column_found

        return columns_mapping

    def get_columns(
        self,
        database: str = None,
        table: str = None,
        header_file: str = None,
        sql_query: str = None,
    ) -> list:
        """
        The function `get_columns` retrieves a list of column names from a specified database and table
        using SQL queries.

        :param database: The `database` parameter in the `get_columns` function is used to specify the
        name of the database from which you want to retrieve the column names. If this parameter is not
        provided, the function will default to using the current database
        :type database: str
        :param table: The `table` parameter in the `get_columns` function represents the name of the
        table in the database for which you want to retrieve the column names. If this parameter is not
        provided, the function will attempt to get the table name from the specified database. If the
        table parameter is not specified and
        :type table: str
        :param header_file: The `header_file` parameter in the `get_columns` function is used to specify
        the file containing the header information for the data source. This information is often used
        in cases where the column names are not explicitly defined in the database schema or where the
        data is stored in a file format that requires additional
        :type header_file: str
        :param sql_query: The `sql_query` parameter in the `get_columns` function is used to specify a
        custom SQL query to retrieve column names from the database table. If a `sql_query` is provided,
        the function will execute that query to get the column names and return them as a list
        :type sql_query: str
        :return: The function `get_columns` returns a list of column names for a given database and
        table. If a SQL query is provided, it executes the query and returns the column names from the
        result. If no database is specified, it uses the current database. It then checks the database
        format and connects to the database accordingly to retrieve the column names using a SQL query.
        If the table parameter is not provided
        """

        if not database:
            database = self.get_database()

        if not table:
            table = self.get_database_table(database=database)

        if sql_query:
            columns_list = list(database.query(sql_query).columns)
            return columns_list

        try:
            if database and self.exists(database):
                database_format = self.get_format(database)
                if database_format in ["conn"]:
                    if table:
                        database_conn = database
                        sql_query = f"SELECT * FROM {table} LIMIT 0"
                        columns_list = list(database_conn.query(sql_query).columns)
                        return columns_list
                elif database_format in ["duckdb"]:
                    if table:
                        database_conn = duckdb.connect(database)
                        sql_query = f"SELECT * FROM {table} LIMIT 0"
                        columns_list = list(database_conn.query(sql_query).columns)
                        database_conn.close()
                        return columns_list
                elif database_format in ["sqlite"]:
                    if table:
                        database_conn = sqlite3.connect(database)
                        sql_query = f"SELECT * FROM {table} LIMIT 0"
                        columns_list = list(
                            pd.read_sql_query(sql_query, database_conn).columns
                        )
                        return columns_list
                elif database_format in [
                    "parquet",
                    "vcf",
                    "tsv",
                    "csv",
                    "tbl",
                    "bed",
                    "json",
                ]:
                    sql_from = self.get_sql_from(
                        database=database, header_file=header_file
                    )
                    sql_query = f"SELECT * FROM {sql_from} LIMIT 0"
                    return list(self.conn.query(sql_query).columns)
        except ValueError:
            return []

        return []

    def get_table_columns_from_format(self, database: str = None) -> list:
        """
        The function `get_table_columns_from_format` returns a list of table columns based on the
        specified database format.

        :param database: The `database` parameter is a string that represents the name of the database.
        It is an optional parameter, which means it has a default value of `None`. If no value is
        provided for the `database` parameter, the `get_database()` method is called to retrieve the
        current database name
        :type database: str
        :return: a list of table columns.
        """

        table_columns = None

        if not database:
            database = self.get_database()

        database_format = self.get_format(database)

        needed_columns = DATABASE_TYPE_NEEDED_COLUMNS.get(database_format, None)
        if needed_columns:
            table_columns = list(needed_columns.keys())
        else:
            table_columns = []

        return table_columns

    def get_table_columns_from_file(
        self,
        database: str = None,
        header_file: str = None,
        header_file_find: bool = True,
    ) -> list:
        """
        The function `get_table_columns_from_file` retrieves the column names from a database or header
        file.

        :param database: The `database` parameter is a string that represents the name or path of the
        database file. If this parameter is not provided, the `get_database()` method is called to
        retrieve the database name or path
        :type database: str
        :param header_file: The `header_file` parameter is a string that represents the file path or
        name of the header file. This file contains the header information for a table, which typically
        includes the names of the columns in the table
        :type header_file: str
        :param header_file_find: Allow header file find if not provided
        :type header_file_find: bool
        :return: a list of table columns.
        """

        table_columns = None

        if not database:
            database = self.get_database()

        if not header_file and header_file_find:
            header_file = self.get_header_file()

        if not header_file and header_file_find:
            header_file = self.find_header_file(database)

        database_format = self.get_format(database=database)
        delimiter = SEP_TYPE.get(database_format, "\t")

        # Try from database file
        try:
            table_header = self.read_header_file(database)
        except ValueError:
            table_header = None

        if table_header:
            try:
                table_columns = (
                    table_header[self.get_header_length(header_file=table_header)]
                    .strip()
                    .split(delimiter)
                )
            except IndexError:
                table_columns = None
        else:
            table_columns = None

        if not table_columns:
            # Try from header file
            try:
                table_header = self.read_header_file(header_file)
            except ValueError:
                table_header = None

            if table_header:
                try:
                    table_columns = (
                        table_header[self.get_header_length(header_file=table_header)]
                        .strip()
                        .split(delimiter)
                    )
                except IndexError:
                    table_columns = None
            else:
                table_columns = None

        return table_columns

    def get_annotations(self, database: str = None) -> object:
        """
        This function returns the annotations of a database or the default database if none is
        specified.

        :param database: The parameter `database` is a string that represents the name of the database
        to retrieve annotations from. If no database name is provided, the method will use the default
        database name obtained from the `get_database()` method
        :type database: str
        :return: The function `get_annotations` returns the `infos` attribute of the header of a
        database. If the `database` parameter is not provided, it gets the current database using the
        `get_database` method. If there is no header, it returns `None`.
        """

        if not database:
            database = self.get_database()

        if self.get_header(database=database):
            return self.get_header(database=database).infos
        else:
            return None

    def get_extra_columns(
        self, database: str = None, database_type: str = None, sql_query: str = None
    ) -> list:
        """
        This Python function returns a list of extra columns in a database table that are not needed
        based on the database type and existing columns.

        :param database: A string representing the name of the database to retrieve columns from. If
        None is provided, the default database will be used
        :type database: str
        :param database_type: The `database_type` parameter in the `get_extra_columns` function
        represents the type of the database for which you want to retrieve the list of extra columns. It
        is used to determine which columns are needed based on the database type and the existing
        columns in the specified database table
        :type database_type: str
        :param sql_query: The `sql_query` parameter in the `get_extra_columns` function is used to pass
        an SQL query that can be used to retrieve specific columns from the database. This query can be
        customized to filter columns based on certain conditions or criteria before analyzing them to
        determine the extra columns that are not needed
        :type sql_query: str
        :return: A list of extra columns in a database table that are not needed based on the database
        type and existing columns.
        """

        if not database:
            database = self.get_database()

        if not database:
            return []

        existing_columns = self.get_columns(
            database=database,
            table=self.get_database_table(database),
            sql_query=sql_query,
        )
        if not database_type:
            database_type = self.get_type(database=database, sql_query=sql_query)
        needed_columns = self.get_needed_columns(
            database_columns=existing_columns, database_type=database_type
        )

        extra_columns = existing_columns.copy()

        for needed_col in needed_columns:
            if needed_columns.get(needed_col) in extra_columns:
                extra_columns.remove(needed_columns.get(needed_col))

        return extra_columns

    def is_vcf(self, database: str = None, sql_query: str = None) -> bool:
        """
        The `is_vcf` function checks if a given database is of type "vcf" by examining its columns and
        their types.

        :param database: The `database` parameter in the `is_vcf` function is a string that represents
        the name of the database that the function will use to check if the file is a VCF (Variant Call
        Format) file. If the `database` parameter is not provided when calling the function, it will
        :type database: str
        :param sql_query: The `sql_query` parameter in the `is_vcf` function is used to pass an SQL
        query string that can be used to filter the columns retrieved from the database. This query can
        be used to narrow down the columns that are considered when checking if the database is of type
        "vcf"
        :type sql_query: str
        :return: The function `is_vcf` returns a boolean value indicating whether the database type is
        "vcf" or not.
        """

        if not database:
            database = self.get_database()

        if not database:
            return False

        database_columns = self.get_columns(
            database=database,
            table=self.get_database_table(database),
            sql_query=sql_query,
        )

        # Assume VCF is 8 needed columns, either only or with extra FORMAT column (assume other are samples)
        return self.get_type_from_columns(
            database_columns=database_columns, check_database_type="vcf"
        ) == "vcf" and ("FORMAT" in database_columns or len(database_columns) == 8)

    def get_conn(self):
        """
        The function returns the connection object.
        :return: The method is returning the value of the instance variable `self.conn`.
        """

        return self.conn

    def is_genotype_column(
        self,
        column: str,
        database: str = None,
        downsampling: int = 1000,
        check_format: bool = True,
    ) -> bool:
        """
        The `is_genotype_column` function in Python checks if a specified column in a database contains
        genotype data based on a regular expression pattern.

        :param column: The `column` parameter is a string that represents the name of a column in a
        database table. It is used to specify the column for which you want to check if it contains
        genotype information based on a regular expression pattern
        :type column: str
        :param database: The `database` parameter in the `is_genotype_column` method is used to specify
        the name of the database from which the data will be queried. If a database is provided, the
        method will query the specified database to check if the given column contains genotype
        information. If no database is provided,
        :type database: str
        :param downsampling: The `downsampling` parameter in the `is_genotype_column` method is an
        integer value that determines the number of rows to be sampled from the database table when
        checking for genotype information in the specified column. This parameter is used to limit the
        number of rows to be processed in order to improve performance, defaults to 1000
        :type downsampling: int (optional)
        :param check_format: The `check_format` parameter in the `is_genotype_column` method is a
        boolean flag that determines whether the function should check the format of the data before
        proceeding with the genotype column analysis. If `check_format` is set to `True`, the function
        will verify if the specified column exists in, defaults to True
        :type check_format: bool (optional)
        :return: The `is_genotype_column` method returns a boolean value. If the specified column in a
        database table contains genotype information, it returns `True`; otherwise, it returns `False`.
        """

        # Table variants
        table_variants_from = self.get_sql_database_link(database=database)

        # Check if format column is present
        if check_format:
            query = f"""
                SELECT * FROM {table_variants_from}
                LIMIT 0
            """
            df = self.query(query=query)
            if "FORMAT" not in df.columns or column not in df.columns:
                return False
            query_format = f"""
                AND (len(string_split(CAST("FORMAT" AS VARCHAR), ':')) = len(string_split(CAST("{column}" AS VARCHAR), ':')) OR "{column}" == './.')
            """
        else:
            query_format = ""

        # Query number of samples
        query_downsampling = f"""
            SELECT "{column}", FORMAT
            FROM {table_variants_from}
            LIMIT {downsampling}
        """
        df_downsampling = self.query(query=query_downsampling)

        # Query to check genotype
        query_genotype = f"""
            SELECT  *
            FROM df_downsampling
            WHERE (
                regexp_matches(CAST("{column}" AS VARCHAR), '^[0-9.]([/|][0-9.])+')
                {query_format}
                )
        """
        df_genotype = self.query(query=query_genotype)

        # return
        return len(df_genotype) == len(df_downsampling)

    def export(
        self,
        output_database: str,
        output_header: str = None,
        header_in_output: bool = True,
        database: str = None,
        table: str = "variants",
        parquet_partitions: list = None,
        threads: int = 1,
        sort: bool = False,
        index: bool = False,
        existing_columns_header: list = [],
        order_by: str = "",
        query: str = None,
        compression_type: str = None,
        chunk_size: int = 1000000,
        export_mode: str = "pyarrow",
        compresslevel: int = 6,
        export_header: bool = True,
        sample_list: list = None,
    ) -> bool:
        """
        The `export` function exports data from a database to a specified output format, compresses it
        if necessary, and returns a boolean value indicating whether the export was successful or not.

        :param output_database: The `output_database` parameter is a string that represents the path and
        filename of the output file to be exported. It specifies where the exported data will be saved
        :type output_database: str
        :param output_header: The `output_header` parameter is an optional string that represents the
        header of the output file. If provided, it specifies the header that will be included in the
        output file. If not provided, the header will be automatically detected based on the output file
        format
        :type output_header: str
        :param header_in_output: The `header_in_output` parameter is a boolean value that determines
        whether the header should be included in the output file. If set to `True`, the header will be
        included in the output file. If set to `False`, the header will not be included in the output
        file. By default,, defaults to True
        :type header_in_output: bool (optional)
        :param database: The `database` parameter is the name of the database from which you want to
        export data. If this parameter is not provided, the function will use the `get_database()`
        method to retrieve the current database
        :type database: str
        :param table: The `table` parameter specifies the name of the table in the database from which
        the data will be exported. By default, if not specified, it is set to "variants", defaults to
        variants
        :type table: str (optional)
        :param parquet_partitions: The `parquet_partitions` parameter is a list that specifies the
        partition columns for the Parquet output format. Each element in the list represents a partition
        column. The partitions are used to organize the data in the Parquet file based on the values of
        the specified columns
        :type parquet_partitions: list
        :param threads: The `threads` parameter in the `export` function is an optional integer that
        specifies the number of threads to use for exporting the data. It determines the level of
        parallelism during the export process. By default, it is set to 1, defaults to 1
        :type threads: int (optional)
        :param sort: The `sort` parameter in the `export` function is a boolean value that specifies
        whether the output file should be sorted based on the genomic coordinates of the variants. If
        `sort` is set to `True`, the output file will be sorted. If `sort` is set to `False`,, defaults
        to False
        :type sort: bool (optional)
        :param index: The `index` parameter is a boolean value that specifies whether to index the
        output file. If `index` is set to `True`, the output file will be indexed. If `index` is set to
        `False` or not provided, the output file will not be indexed. By default,, defaults to False
        :type index: bool (optional)
        :param existing_columns_header: The `existing_columns_header` parameter is a list that
        represents the existing columns in the header of the output file. It is used to determine the
        columns that should be included in the output file. If this parameter is not provided, the
        function will automatically detect the header columns based on the output file format
        :type existing_columns_header: list
        :param order_by: The `order_by` parameter in the `export` function is a string that specifies
        the columns by which the output file should be ordered. You can specify multiple columns
        separated by commas. Each column can be followed by the keyword "ASC" (ascending) or "DESC"
        (descending) to specify
        :type order_by: str
        :param query: The `query` parameter in the `export` function represents a SQL query that
        specifies the data to be exported from the database. If provided, the function will export the
        result of this query. If the `query` parameter is not provided, the function will generate a
        query to export the data from
        :type query: str
        :param compression_type: The `compression_type` parameter in the `export` function specifies the
        type of compression to be applied to the output file. By default, the compression type is set to
        "bgzip". This parameter allows you to choose the compression algorithm for the output file, such
        as "gzip", "bgzip
        :type compression_type: str
        :param chunk_size: The `chunk_size` parameter in the `export` function specifies the size of
        each chunk or batch of data that will be processed during the export operation. It determines
        how many records or lines of data will be included in each chunk that is processed at a time,
        defaults to 1000000
        :type chunk_size: int (optional)
        :param export_mode: The `export_mode` parameter in the `export` function specifies the mode of
        export, which can be either "pyarrow" or "duckdb", defaults to pyarrow
        :type export_mode: str (optional)
        :param compresslevel: The `compresslevel` parameter in the `export` function represents the
        level of compression for gzip. By default, it is set to 6. This parameter allows you to specify
        the compression level when using gzip compression for the output file. The compression level can
        range from 0 (no compression), defaults to 6
        :type compresslevel: int (optional)
        :param export_header: The `export_header` parameter is a boolean flag that determines whether
        the header of a VCF file should be exported to a separate file or not. If `export_header` is
        True, the header will be exported to a file. If `export_header` is False, the header will not
        be, defaults to True
        :type export_header: bool (optional)
        :param sample_list: The `sample_list` parameter in the `export` function is a list that
        specifies the samples to be included in the exported data. If provided, the samples listed in
        this parameter will be included in the output file. If not provided, the function will determine
        the samples to include based on the data
        :type sample_list: list
        :return: The `export` function returns a boolean value indicating whether the export was
        successful or not.
        """

        # Full path
        output_database = full_path(output_database)
        output_header = full_path(output_header)
        database = full_path(database)

        # Database
        if not database:
            database = self.get_database()
        if not database:
            return False

        # Chunk size
        if not chunk_size:
            chunk_size = 1000000
        else:
            chunk_size = int(chunk_size)

        # Export mode
        # Either "pyarrow" (default) or "duckdb"
        if not export_mode:
            export_mode = "pyarrow"

        # Compression level
        if not compresslevel:
            compresslevel = 6

        # Remove output if exists
        remove_if_exists(output_database)

        # Tmp
        tmp_folder = os.path.dirname(output_database)
        if not tmp_folder:
            tmp_folder = "."

        with TemporaryDirectory(
            dir=tmp_folder, prefix="howard_database_export_"
        ) as tmp_dir:

            # tmp files
            tmp_files = []

            # query_set
            query_set = ""

            # Header columns
            if not existing_columns_header and output_header:
                existing_columns_header = self.get_header_file_columns(output_header)

            # Auto-detect output type and compression and delimiter
            output_type = get_file_format(output_database)
            compressed = self.is_compressed(database=output_database)
            delimiter = FILE_FORMAT_DELIMITERS.get(output_type, "\t")

            # database type
            if output_type in ["vcf"]:
                database_type = "vcf"
            elif output_type in ["bed"]:
                database_type = "regions"
            else:
                database_type = self.get_type(database=database, sql_query=query)

            # database object
            # If database is string, then create database conn
            if isinstance(database, str):
                database_conn = Database(database=database).get_conn()
            else:
                database_conn = database

            # Existing columns
            existing_columns = self.get_columns(
                database=database_conn,
                table=self.get_database_table(database=database),
                sql_query=query,
            )

            # Extra columns
            extra_columns = self.get_extra_columns(
                database=database_conn, database_type=output_type, sql_query=query
            )

            # Needed columns
            needed_columns = self.get_needed_columns(
                database_columns=existing_columns, database_type=database_type
            )

            # Order by
            order_by_list = []
            if order_by:
                # Split order by options
                order_by_split = order_by.split(",")
                for order_by_option in order_by_split:
                    # Split order by option
                    order_by_option_split = order_by_option.strip().split(" ")
                    order_by_option_split_column = order_by_option_split[0]
                    if len(order_by_option_split) > 1:
                        order_by_option_split_order = order_by_option_split[1]
                    else:
                        order_by_option_split_order = "ASC"
                    # Chek if column exists
                    if (
                        order_by_option_split_column.replace('"', "").strip()
                        in existing_columns
                    ):
                        order_by_list.append(
                            f"{order_by_option_split_column} {order_by_option_split_order}"
                        )

            # Clean order by
            order_by_clean = ", ".join(order_by_list)

            # Query values
            default_empty_value = ""
            query_export_format = None
            include_header = False
            post_process = False
            order_by_sql = ""
            post_process_just_move = False

            # export options
            export_options = {}

            # VCF
            if output_type in ["vcf"]:
                if not self.is_vcf(database=database_conn, sql_query=query):
                    extra_columns = []
                else:
                    extra_columns = existing_columns_header

                # Check VCF format with extra columns
                extra_columns_clean = []

                # Force samples list in parameter
                if sample_list:
                    if "FORMAT" in extra_columns:
                        extra_columns = ["FORMAT"] + sample_list
                    else:
                        extra_columns = sample_list

                # Check columns
                else:
                    for extra_column in extra_columns:
                        if extra_column not in needed_columns and (
                            extra_column == "FORMAT"
                            or (
                                "FORMAT" in extra_columns_clean
                                and self.is_genotype_column(
                                    database=database, column=extra_column
                                )
                            )
                        ):
                            extra_columns_clean.append(extra_column)
                    extra_columns = extra_columns_clean

                default_empty_value = "."
                query_export_format = f"FORMAT CSV, DELIMITER '{delimiter}', HEADER, QUOTE '', COMPRESSION 'gzip'"
                include_header = True
                post_process = True
                # Export options
                if not compression_type:
                    compression_type = "bgzip"
                export_options = {
                    "format": "CSV",
                    "delimiter": delimiter,
                    "header": True,
                    "quote": None,
                    "compression": compression_type,
                }
                compresslevel = 1

            # TSV/CSV/TBL
            elif output_type in ["tsv", "csv", "tbl"]:
                if output_type in ["csv", "tbl"]:
                    quote = '"'
                else:
                    quote = ""
                query_export_format = f"FORMAT CSV, DELIMITER '{delimiter}', HEADER, QUOTE '{quote}', COMPRESSION 'gzip'"
                if delimiter in ["\t"]:
                    include_header = header_in_output and True
                post_process = True
                if order_by_clean:
                    order_by_sql = f"ORDER BY {order_by_clean}"
                # Export options
                if not compression_type:
                    compression_type = "gzip"
                export_options = {
                    "format": "CSV",
                    "delimiter": delimiter,
                    "header": True,
                    "quote": quote,
                    "compression": compression_type,
                }

            # JSON
            elif output_type in ["json"]:
                query_export_format = "FORMAT JSON, ARRAY TRUE"
                include_header = False
                post_process = True
                if order_by_clean:
                    order_by_sql = f"ORDER BY {order_by_clean}"
                # Export options
                if not compression_type:
                    compression_type = "gzip"
                export_options = {
                    "format": "JSON",
                    "array": True,
                    "compression": compression_type,
                }

            # Parquet
            elif output_type in ["parquet"]:
                query_export_format = "FORMAT PARQUET"
                # Export options
                export_options = {
                    "format": "PARQUET",
                }
                include_header = False
                post_process = True

            # BED
            elif output_type in ["bed"]:
                query_export_format = f"FORMAT CSV, DELIMITER '{delimiter}', HEADER"
                include_header = True
                post_process = True
                if order_by_clean:
                    order_by_sql = f"ORDER BY {order_by_clean}"
                # Export options
                if not compression_type:
                    compression_type = "gzip"
                export_options = {
                    "format": "CSV",
                    "delimiter": delimiter,
                    "header": True,
                    "quote": None,
                    "compression": compression_type,
                }

            # duckDB
            elif output_type in ["duckdb"]:

                # Needed column
                needed_columns = []

                # Export database as Parquet
                database_export_parquet_file = os.path.join(tmp_dir, "output.parquet")
                self.export(
                    database=database, output_database=database_export_parquet_file
                )

                # Create database and connexion
                output_database_conn = duckdb.connect(output_database)

                # Create table in database connexion with Parquet file
                query_copy = f""" 
                    CREATE TABLE {table}
                    AS {self.get_sql_database_link(database=database_export_parquet_file)}
                    """
                output_database_conn.execute(query_copy)

                # Close connexion
                output_database_conn.close()

                # remove tmp
                remove_if_exists([database_export_parquet_file])

                return os.path.exists(output_database)

            # Partition
            if parquet_partitions:
                parquet_partitions_clean = []
                parquet_partitions_array = []
                for parquet_partition in parquet_partitions:
                    parquet_partitions_array.append(
                        parquet_partition.translate({'"': None, "'": None, " ": None})
                    )
                    parquet_partitions_clean.append(
                        '"'
                        + parquet_partition.translate({'"': None, "'": None, " ": None})
                        + '"'
                    )
                parquet_partitions_by = ",".join(parquet_partitions_clean)
                query_export_format += (
                    f", PARTITION_BY ({parquet_partitions_by}), OVERWRITE_OR_IGNORE"
                )
                export_options["partition_by"] = parquet_partitions_array
                if export_options.get("format", None) == "CSV":
                    export_mode = "duckdb"
                    post_process_just_move = True

            # Construct query columns
            query_columns = []

            # Add Needed columns
            for needed_column in needed_columns:
                if needed_columns[needed_column]:
                    query_column_name = needed_columns[needed_column]
                    query_column = f""" "{needed_columns[needed_column]}" """
                else:
                    query_column_name = default_empty_value
                    query_column = f""" '{default_empty_value}' """
                query_column_as = f""" "{needed_column}" """
                if query_column_name == needed_column:
                    query_columns.append(f""" {query_column} """)
                else:
                    query_columns.append(f""" {query_column} AS {query_column_as} """)

            # Add Extra columns
            for extra_column in extra_columns:
                if extra_column not in needed_columns:
                    # query_columns.append(f""" "{extra_column}" AS "{extra_column}" """)
                    query_columns.append(f""" "{extra_column}" """)

            # Query export columns
            query_export_columns = f""" {",".join(query_columns)} """

            if query_columns:

                # Compressed tmp file
                query_output_database_tmp = os.path.join(tmp_dir, "output")

                # Query
                # If no query, generate query of the database
                if not query:
                    query = f"""
                        SELECT {query_export_columns}
                        FROM {self.get_sql_database_link(database=database)}
                        {order_by_sql}
                    """

                # Test empty query
                df = self.conn.execute(query).fetch_record_batch(1)
                query_empty = True
                for d in df:
                    query_empty = False
                    break
                if query_empty:
                    log.error("Export failed: Empty")
                    raise ValueError("Export failed: Empty")

                # Schema names
                schema_names = None

                # Export mode pyarrow
                if export_mode == "pyarrow":

                    # Compression mode
                    # If compress required and compression type as gzip or bgzip
                    # For bgzip compression, recompression will be done, and compression with gzip for tmp file done (level 1)
                    # to reduce tmp file size
                    compression_mode_gzip = compressed and (
                        export_options.get("compression", None) in ["gzip", "bgzip"]
                    )

                    # File stream mode (str or bytes)
                    f_mode = ""
                    if compression_mode_gzip:
                        f_mode = "b"

                    if include_header:

                        # Open stream files (uncompressed and compressed)
                        with open(query_output_database_tmp, mode="w") as f, pgzip.open(
                            query_output_database_tmp,
                            mode="w",
                            thread=threads,
                            compresslevel=compresslevel,
                        ) as f_gz:

                            # Switch to compressed stream file
                            if compression_mode_gzip:
                                f = f_gz

                            # Generate header tmp file
                            query_output_header_tmp = os.path.join(tmp_dir, "header")
                            self.get_header_file(
                                header_file=query_output_header_tmp,
                                remove_header_line=True,
                                sql_query=query,
                            )

                            # Write header to tmp file
                            with open(
                                query_output_header_tmp, "r" + f_mode
                            ) as output_header_tmp:
                                f.write(output_header_tmp.read())

                    # JSON format - Add special "[" character at the beginning of the file
                    if export_options.get("format") in ["JSON"]:

                        # Open stream files (uncompressed and compressed)
                        with open(query_output_database_tmp, mode="a") as f, pgzip.open(
                            query_output_database_tmp,
                            mode="a",
                            thread=threads,
                            compresslevel=compresslevel,
                        ) as f_gz:

                            # Switch to compressed stream file
                            if compression_mode_gzip:
                                f = f_gz
                                f.write(b"[\n")
                            else:
                                f.write("[\n")

                    # Open stream files (uncompressed and compressed) for chunk
                    with open(query_output_database_tmp, mode="a") as f, pgzip.open(
                        query_output_database_tmp,
                        mode="a",
                        thread=threads,
                        compresslevel=compresslevel,
                    ) as f_gz:

                        # Switch to compressed stream file
                        if compression_mode_gzip:
                            f = f_gz

                        # Chunk query with batch of dataframes of chunk_size
                        df = self.conn.execute(query).fetch_record_batch(chunk_size)

                        # id of chunk
                        i = 0

                        # For each chunk dataframe
                        for d in df:

                            # Schema names
                            schema_names = d.schema.names

                            # id of chunk
                            i += 1

                            # Log - number of records
                            log.debug(f"Chunk {i}: records process...")

                            # Check process for first chunk
                            if i == 1:

                                # If include header in file
                                header = export_options.get("header", True)

                                # Parquet output format
                                # Either a folder or a writer
                                if export_options.get("format") in ["PARQUET"]:

                                    # For Parquet with multiple file - folder
                                    if export_options.get(
                                        "partition_by", None
                                    ) or export_options.get("per_thread_output", False):
                                        query_output_database_tmp = (
                                            f"{query_output_database_tmp}.parquet"
                                        )

                                    # For Parquet as a unique file - writer
                                    else:
                                        writer = pq.ParquetWriter(
                                            query_output_database_tmp, d.schema
                                        )

                            else:

                                # Switch of header in file for not first chunk
                                header = False

                            # CSV format
                            if export_options.get("format") in ["CSV"]:

                                # With quote option
                                if export_options.get("quote", None):

                                    # Polars write dataframe
                                    pl.from_arrow(d).write_csv(
                                        file=f,
                                        separator=export_options.get("delimiter", ""),
                                        include_header=header,
                                        quote_char=export_options.get("quote", '"'),
                                    )

                                # Without quote option
                                else:

                                    # Polars write dataframe
                                    pl.from_arrow(d).write_csv(
                                        file=f,
                                        separator=export_options.get("delimiter", ""),
                                        include_header=header,
                                        quote_style="never",
                                    )

                            # JSON format
                            elif export_options.get("format") in ["JSON"]:

                                # Compressed mode gzip
                                if compression_mode_gzip:

                                    # Add comma at the beginning of dataframe (if not the first one) in bytes mode
                                    if i > 1:
                                        f.write(b",\n")

                                    # Write dataframe in bytes mode
                                    f.write(
                                        str.encode(
                                            pl.from_arrow(d)
                                            .write_ndjson()
                                            .replace("\n{", ",\n{")
                                            .replace("[", "")
                                            .replace("]", "")
                                        )
                                    )

                                # Not compressed mode gzip (string mode)
                                else:

                                    # Add comma at the beginning of dataframe (if not the first one) in string mode
                                    if i > 1:
                                        f.write(",\n")

                                    # Write dataframe in string mode
                                    f.write(
                                        pl.from_arrow(d)
                                        .write_ndjson()
                                        .replace("\n{", ",\n{")
                                        .replace("[", "")
                                        .replace("]", "")
                                    )

                            # Parquet format
                            elif export_options.get("format") in ["PARQUET"]:

                                # Partition by fields
                                partition_by = export_options.get("partition_by", None)

                                if partition_by:

                                    # For No partition but split parquet files into a folder
                                    if "None" in partition_by:
                                        partition_by = None

                                    # Pyarrow write
                                    pq.write_to_dataset(
                                        pa.Table.from_batches([d]),
                                        query_output_database_tmp,
                                        partition_cols=partition_by,
                                        use_threads=threads,
                                        existing_data_behavior="overwrite_or_ignore",
                                    )

                                # Parquet in unique file
                                else:
                                    writer.write_batch(d)

                    # Close Parquet writer
                    if (
                        export_options.get("format") in ["PARQUET"]
                        and not export_options.get("partition_by", None)
                        and not export_options.get("per_thread_output", None)
                    ):
                        writer.close()

                    # JSON format - Add special "]" character at the end of the file
                    if export_options.get("format") in ["JSON"]:

                        # Open stream files (uncompressed and compressed)
                        with open(query_output_database_tmp, mode="a") as f, pgzip.open(
                            query_output_database_tmp,
                            mode="a",
                            thread=threads,
                            compresslevel=compresslevel,
                        ) as f_gz:

                            # Switch to compressed stream file
                            if compression_mode_gzip:
                                f = f_gz
                                f.write(b"]\n")
                            else:
                                f.write("]\n")

                # Export mode duckdb
                elif export_mode == "duckdb":

                    # Create COPY TO query
                    query_copy = f""" 
                        {query_set}
                        COPY (
                            {query}
                            )
                        TO '{query_output_database_tmp}'
                        WITH ({query_export_format})
                        """

                    # Export with duckdb
                    self.query(query=query_copy)

                # Export mode unknown
                else:
                    log.error(f"Export mode '{export_mode}' unknown")
                    raise ValueError(f"Export mode '{export_mode}' unknown")

                # Post process
                if post_process:

                    # Log - number of records
                    log.debug("Post processing...")

                    # Input files
                    input_files = []

                    # Export mode duckdb and include header
                    if export_mode == "duckdb" and include_header:

                        # create tmp header file
                        query_output_header_tmp = os.path.join(tmp_dir, "header")
                        tmp_files.append(query_output_header_tmp)
                        self.get_header_file(
                            header_file=query_output_header_tmp, remove_header_line=True
                        )

                        # Add tmp header file for concat and compress
                        input_files.append(query_output_header_tmp)

                    # Add variants file
                    input_files.append(query_output_database_tmp)

                    # Output with concat and compress
                    if not post_process_just_move and (
                        export_mode == "duckdb"
                        or (
                            compressed
                            and export_options.get("compression", None) == "bgzip"
                        )
                        or sort
                        or index
                    ):

                        # Compression type
                        if not compressed:
                            compression_type = "none"
                        else:
                            compression_type = export_options.get(
                                "compression", "bgzip"
                            )

                        # Concat and compress
                        concat_and_compress_files(
                            input_files=input_files,
                            output_file=output_database,
                            compression_type=compression_type,
                            threads=threads,
                            sort=sort,
                            index=index,
                        )

                    # Output already generated file (either compressed in gzip or not compressed, with included header if needed)
                    else:

                        # Move tmp file
                        shutil.move(query_output_database_tmp, output_database)

                # Generate associated header file
                if output_header and export_header:

                    # Log - Generate header
                    log.debug("Generate header...")

                    # Create database
                    database_for_header = Database(database=output_database)

                    # Remove header if exists
                    remove_if_exists([output_header])

                    # Find columns in database
                    if schema_names:
                        header_columns_from_database = schema_names
                    else:
                        header_columns_from_database = (
                            database_for_header.get_header_columns_from_database(
                                database=output_database
                            )
                        )

                    # Generate header file
                    database_for_header.get_header_file(
                        header_file=output_header,
                        replace_header_line=header_columns_from_database,
                        force=True,
                    )

            # Clean tmp files (deprecated)
            remove_if_exists(tmp_files)

            # Return if file exists
            return os.path.exists(output_database) and self.get_type(output_database)
