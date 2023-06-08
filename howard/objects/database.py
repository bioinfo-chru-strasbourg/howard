
# import pyarrow as pa
# import pyarrow.csv as csv
import fileinput
import hashlib
import random
from shutil import copyfileobj
import string
import polars as pl
import pandas as pd
import duckdb
import sqlite3
import vcf
#from Bio import bgzf
import Bio.bgzf as bgzf

from howard.commons import *
from howard.tools.databases import *


SEP_TYPE = {
    "vcf" : "\t",
    "tsv" : "\t",
    "csv" : ",",
    "tbl" : "|",
    "bed" : "\t",
}

DATABASE_TYPE_NEEDED_COLUMNS = {
    "variants":
        {
            "#CHROM": ["#CHROM", "CHROM", "CHR", "CHROMOSOME"],
            "POS": ["POS"],
            "REF": ["REF"],
            "ALT": ["ALT"],
        },
    "regions":
        {
            "#CHROM": ["#CHROM", "CHROM", "CHR", "CHROMOSOME"],
            "START": ["START", "POSITIONSTART", "POS"],
            "END": ["END", "POSITIONEND", "POS"],
        },
    "vcf":
        {
            "#CHROM": ["#CHROM", "CHROM", "CHR", "CHROMOSOME"],
            "POS": ["POS", "POSITION"],
            "ID": ["ID", "IDENTIFIER"],
            "REF": ["REF", "REFERENCE"],
            "ALT": ["ALT", "ALTERNATIVE"],
            "QUAL": ["QUAL", "QUALITY"],
            "FILTER": ["FILTER"],
            "INFO": ["INFO"],
        },
}

DEFAULT_VCF_HEADER = [
    "#CHROM",
    "POS",
    "ID",
    "REF",
    "ALT",
    "QUAL",
    "FILTER",
    "INFO"
]

DEFAULT_HEADER_LIST = [
            '##fileformat=VCFv4.2',
            '#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO'
        ]

CODE_TYPE_MAP = {
            "Integer": 0,
            "String": 1,
            "Float": 2,
            "Flag": 3
        }

FILE_FORMAT_DELIMITERS = {
    "vcf": "\t",
    "tsv": "\t",
    "csv": ",",
    "tbl": "|",
    "bed": "\t"
}

DTYPE_LIMIT_AUTO = 10000


class Database:

    def __init__(self, database:str = None, format:str = None, header:str = None, header_file:str = None, databases_folders:list = None, assembly:str = None, conn = None, table:str = None) -> None:
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
        if conn:
            self.conn = conn
        elif type(database) == duckdb.DuckDBPyConnection:
            self.conn = database
        else:
            self.conn = duckdb.connect()

        # Install sqlite scanner
        self.conn.query("INSTALL sqlite_scanner; LOAD sqlite_scanner; ")

        # Check attributes
        self.set_format(format=format)
        self.set_assembly(assembly=assembly)
        self.set_databases_folders(databases_folders=databases_folders)
        self.set_header_file(header_file=header_file)
        self.set_database(database=database, databases_folders=self.get_database_folders(), format=self.get_format(), assembly=self.get_assembly())
        self.set_header(database=self.get_database(), header=header, header_file=header_file)


    def set_database(self, database:str, databases_folders:list = None, format:str = None, assembly:str = None) -> str:
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
        elif self.find_database(database=database, databases_folders=databases_folders, format=format, assembly=assembly):
            self.database = self.find_database(database=database, databases_folders=databases_folders, format=format, assembly=assembly)


    def set_databases_folders(self, databases_folders:list = ["."]) -> None:
        """
        This function sets the list of folders where databases are located as an attribute of an object.
        
        :param databases_folders: databases_folders is a list parameter that contains the paths to the
        folders where the databases are stored. The default value of the parameter is a list with a
        single element, which is the current directory (".")
        :type databases_folders: list
        """

        self.databases_folders = databases_folders


    def get_database_folders(self) -> list:
        """
        This function returns a list of database folders.
        :return: The method `get_database_folders` is returning a list of database folders. The specific
        list being returned is stored in the instance variable `databases_folders`.
        """

        return self.databases_folders


    def read_header_file(self, header_file:str = None) -> list:
        """
        This function reads the header of a VCF file and returns a list of the header lines.
        
        :param header_file: The path to the VCF file header that needs to be read
        :type header_file: str
        :return: a list of header lines of a VCF file.
        """
        
        if not header_file:
            return []
        
        elif not os.path.exists(header_file):
            return []

        else:

            header_file_compressed = get_file_compressed(header_file)

            if header_file_compressed:
                with bgzf.open(header_file, 'rt') as header_lines:
                    header_list = []
                    for line in header_lines:
                        if not line.startswith('#'):
                            break
                        header_list.append(line)
                    return header_list
            else:
                with open(header_file, 'rt') as header_lines:
                    header_list = []
                    for line in header_lines:
                        if not line.startswith('#'):
                            break
                        header_list.append(line)
                    return header_list
                
    
    def get_header_from_list(self, header_list:list = None) -> vcf:
        """
        This function returns a vcf.Reader object with a header generated from a given list or a default
        list.
        
        :param header_list: A list of strings representing the header lines of a VCF file. If this
        parameter is not provided, the function will use a default list of header lines
        :type header_list: list
        :return: a `vcf.Reader` object that reads the VCF header information from a list of strings. The
        list of strings can be provided as an argument to the function, and if no argument is provided,
        a default list is used.
        """

        if not header_list:
            header_list = DEFAULT_HEADER_LIST
            
        return vcf.Reader(io.StringIO("\n".join(header_list)))
    

    def get_header_from_file(self, header_file:str) -> vcf:
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
        

    def find_header_file(self, database:str = None) -> str:
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
        if os.path.exists(f"{database}.hdr"):
            database_header_file = f"{database}.hdr"

        # header within file
        elif database_format in ["vcf", "tsv", "csv", "tbl", "bed"]:
            database_header_file = database

        return database_header_file
    

    def get_header(self, database:str = None, header_file:str = None, header_list:list = None) -> vcf:
        """
        This function returns the header of a VCF file either from a file, a list, or from the object
        itself.
        
        :param header_file: a string representing the path to a file containing the header information
        for a VCF file
        :type header_file: str
        :param header_list: A list containing the header lines of a VCF file
        :type header_list: list
        :return: The method `get_header` returns an object of type `vcf`. However, if none of the
        conditions are met (i.e. `self.header` is not set, `header_file` is not provided, and
        `header_list` is not provided), then `None` is returned.
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
        elif self.get_extra_columns(database=database):
            # Construct header from a list of columns
            return self.get_header_from_columns(database=database, header_columns=self.get_extra_columns(database=database))
        # elif "INFO" in self.get_columns(database=database):
        #     # Construct header from annotation in INFO column (in case of VCF database without header, e.g. in TSV)
        #     # TODO
        #     return None
        else:
            # No header
            return None

    def get_header_from_columns(self, database:str = None, header_columns:list = []) -> object:
        """
        This function generates a VCF header based on the columns in a database and adds custom
        annotations to it.
        
        :param database: The `database` parameter is a string that represents the name of a database. It
        is an optional parameter and if not provided, the `get_database()` method is called to retrieve
        the default database
        :type database: str
        :param header_columns: A list of column names that will be used to generate header information
        for a VCF file. If no header_columns are provided, the function will attempt to automatically
        detect the columns to use based on the database being used
        :type header_columns: list
        :return: a VCF header object that includes information about the columns in a database and their
        data types. The header object is created based on the input parameters, including the database
        name and a list of header columns.
        """
        
        if not database:
            database = self.get_database()

        database_header = vcf.Reader(io.StringIO("\n".join(DEFAULT_HEADER_LIST)))
        
        if not header_columns:
            header_columns = self.get_extra_columns(database=database)

        database_basename = self.get_database_basename(database=database) or "unknown"

        # Columns query for auto detection of stype
        
        # Attach if need
        if self.get_sql_database_attach(database=database, output="attach"):
            self.query(database=database, query=self.get_sql_database_attach(database=database, output="attach"))

        # database columns
        database_query_columns_sql = f""" SELECT * FROM {self.get_sql_database_link(database=database)} LIMIT {DTYPE_LIMIT_AUTO} """
        database_query_columns = self.query(database=database, query=database_query_columns_sql)

        # Remove specific VCF column if is a VCF type
        if self.get_type_from_columns(database_columns=self.get_columns(database=database), check_database_type="vcf") == "vcf":
            header_columns = [x for x in header_columns if x not in DEFAULT_VCF_HEADER]

        # List all columns to add into header
        for header_column in header_columns:

            # Header info type
            header_info_type = "String"
            header_column_df = database_query_columns.df()[header_column]
            header_column_df_dtype = header_column_df.dtype
            if header_column_df_dtype == object:
                if pd.to_numeric(header_column_df, errors='coerce').notnull().all():
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
                CODE_TYPE_MAP[header_info_type]
            )

        # Detach if need
        if self.get_sql_database_attach(database=database, output="detach"):
            self.query(database=database, query=self.get_sql_database_attach(database=database, output="detach"))

        return database_header


    def query(self, database:str = None, query:str = None) -> object:
        """
        This is a Python function that takes in a database and query string as parameters and returns
        the result of the query on the database.
        
        :param database: A string representing the name of the database to query. If no database is
        provided, the method will attempt to retrieve the default database from the connection object
        :type database: str
        :param query: The query parameter is a string that represents the SQL query that needs to be
        executed on the database. It can be any valid SQL statement such as SELECT, INSERT, UPDATE,
        DELETE, etc
        :type query: str
        :return: If a query is provided, the method returns the result of the query executed on the
        database. If no query is provided, the method returns None.
        """
        
        if not database:
            database = self.get_database()

        if query:
            return self.conn.query(query)
        else:
            return None


    def set_header(self, database:str = None, header:vcf = None, header_file:str = None) -> None:
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
        
        if header:

            self.header = header
            self.header_file = None
        
        else:

            if header_file and os.path.exists(header_file):

                # header provided
                self.header = self.get_header(header_file=header_file)
                self.header_file = header_file
            
            else:

                # default no heder file
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


    def set_header_file(self, header_file:str = None) -> None:
        """
        This function sets the header file attribute of an object to the value passed as an argument.
        
        :param header_file: The parameter `header_file` is a string that represents the name or path of
        a header file. This method sets the `header_file` attribute of an object to the value passed as
        an argument. If no argument is passed, the `header_file` attribute remains unchanged
        :type header_file: str
        """

        self.header_file = header_file


    def get_header_file(self, header_file:str = None, remove_header_line:bool = False) -> str:
        """
        This function generates a VCF header file if it does not exist or generates a default header
        file if the provided header file does not match the database.
        
        :param header_file: A string representing the file path and name of the header file. If None,
        the default header file path and name will be used
        :type header_file: str
        :param remove_header_line: A boolean parameter that determines whether to remove the "#CHROM"
        line from the header file. If set to True, the line will be removed, defaults to False
        :type remove_header_line: bool (optional)
        :return: a string which is the name of the header file that was generated or None if no header
        file was generated.
        """
        
        if not header_file:
            header_file = self.header_file
        
        if header_file:
            if header_file != self.get_database():
                if self.get_header():
                    header = self.get_header()
                else:
                    header = self.get_header(header_list=DEFAULT_HEADER_LIST)
                # Generate header file
                f = open(header_file, 'w')
                vcf.Writer(f, header)
                f.close()
        else:
            # No header generated
            header_file = None

        if header_file and remove_header_line:
            os.system(f"sed -i '/^#CHROM/d' {header_file}")
        
        return header_file


    def set_assembly(self, assembly:str = None) -> None:
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
        

    def find_database(self, database:str = None, databases_folders:list = None, format:str = None, assembly:str = None) -> str:
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

        if not format:
            format = self.get_format()
        
        if not assembly:
            assembly = self.get_assembly()

        if not database:
            return None
        
        elif self.exists(database=database):
            return database
        
        else:

            if not databases_folders:
                databases_folders = ['.']

            database_file = None

            for format_extension in ["", f".{format}", f".{format}.gz"]:

                # find in folders
                for databases_folder in databases_folders:
                    
                    # Log
                    log.debug("Annotation file check: " + str(databases_folder+"/"+database+format_extension))

                    # In folder
                    if os.path.exists(databases_folder+"/"+database+format_extension):
                        database_file = databases_folder+"/"+database+format_extension

                    # database found
                    if database_file:
                        break

                # find in subfolder assemby
                if not database_file and assembly:

                    for databases_folder in databases_folders:
                        
                        # Log
                        log.debug("Annotation file check: " + str(databases_folder+"/"+assembly+"/"+database+format_extension))

                        # In folder
                        if os.path.exists(databases_folder+"/"+assembly+"/"+database+format_extension):
                            database_file = databases_folder+"/"+assembly+"/"+database+format_extension

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


    def get_database_basename(self, database:str = None) -> str:
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


    def get_database_dirname(self, database:str = None) -> str:
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


    def exists(self, database:str = None) -> bool:
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
        
        if not database:
            database = self.get_database()
        
        return database and (type(database) == duckdb.DuckDBPyConnection or os.path.exists(database)) 


    def set_format(self, format:str = None) -> str:
        """
        This is a method in a Python class that sets a format attribute to a specified string.
        
        :param format: The format parameter is a string that specifies the desired format for the data.
        It is an optional parameter, meaning that if it is not provided, the format attribute of the
        object will not be changed. The function returns a string indicating the current format of the
        object
        :type format: str
        """

        self.format = format


    def get_format(self, database:str = None) -> str:
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


    def get_type(self, database:str = None) -> str:
        """
        This function gets the type of a database (variants VCF-like or regions BED-like) based on its columns.
        
        :param database: A string representing the name of a database. If None, the function will
        attempt to retrieve the database name using the get_database() method
        :type database: str
        :return: a string that represents the type of database, which can be either "variants"
        (VCF-like) or "regions" (BED-like). If the database is not found or does not exist, the function
        returns None.
        """
        
        if not database:
            database = self.get_database()
        
        if database and self.exists(database):
            database_columns = self.get_columns(database, table=self.get_database_table(database))
            return self.get_type_from_columns(database_columns)
        else:
            return None
    

    def get_database_tables(self, database:str = None) -> str:
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
            format = self.get_format(database)
            if format in ["conn"]:
                database_conn = database
                database_tables = list(database_conn.query("SHOW TABLES;").df()["name"])
                return database_tables
            elif format in ["duckdb"]:
                database_conn = duckdb.connect(database)
                database_tables = list(database_conn.query("SHOW TABLES;").df()["name"])
                database_conn.close()
                return database_tables
            elif format in ["sqlite"]:
                database_conn = sqlite3.connect(database)
                sql_query = f"SELECT name FROM sqlite_master WHERE type='table'"
                database_tables = list(pd.read_sql_query(sql_query, database_conn)["name"])
                return database_tables
            else:
                return None
        else:
            return None


    def get_database_table(self, database:str = None) -> str:
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


    def get_type_from_columns(self, database_columns:list = [], check_database_type:str = None) -> str:
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
                all_needed_columns_found = all_needed_columns_found and (needed_columns[col] is not None)
            if all_needed_columns_found:
                return database_type

        return None
    

    def get_needed_columns(self, database_columns:list = [], database_type:str = None) -> dict:
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


    def get_sql_from(self, database:str = None) -> str:
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

        database_format = self.get_format(database)

        sql_form = None

        if type(database) == duckdb.DuckDBPyConnection:
            sql_form = self.get_database_table(database)
        elif database_format in ["parquet"]:
            sql_form = f"read_parquet('{database}')"
        elif database_format in ["vcf","tsv", "csv", "tbl", "bed"]:
            delimiter = SEP_TYPE.get(database_format,"\t")
            sql_form = f"read_csv('{database}', auto_detect=True, delim='{delimiter}')"
        elif database_format in ["json"]:
            sql_form = f"read_json('{database}', auto_detect=True)"
        elif database_format in ["duckdb"]:
            sql_form = f"'{database}'"
        elif database_format in ["sqlite"]:
            database_table = self.get_database_table(database=database)
            sql_form = f"(SELECT * FROM sqlite_scan('{database}', '{database_table}'))"

        return sql_form
    

    def get_sql_database_attach(self, database:str = None, output:str = "query", ) -> str:
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
    

    def get_sql_database_link(self, database:str = None) -> str:
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

            database_attach_name = self.get_sql_database_attach(database=database, output="name")
            if database_attach_name:
                sql_table = self.get_database_table(database=database)
                sql_from = f"{database_attach_name}.{sql_table}"
                
            sql_database_link = f"(SELECT * FROM {sql_from})"
        
        return sql_database_link


    def is_compressed(self, database:str = None) -> bool:
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
    

    def get_header_infos_list(self, database:str = None) -> list:
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
    

    def find_column(self, database:str = None, table:str = None, column:str = "INFO", prefixes:list = ["INFO/"]) -> str:
        """
        This function finds a specific column in a database table, with the option to search for a
        column with a specific prefix or within the INFO column header.
        
        :param database: The name of the database to search for the column in. If not provided, it will
        use the current database that the code is connected to
        :type database: str
        :param table: The name of the table in the database where the column is located
        :type table: str
        :param column: The default value for the "column" parameter is "INFO", but it can be changed to
        search for a specific column name, defaults to INFO
        :type column: str (optional)
        :param prefixes: The prefixes parameter is a list of strings that are used to search for a
        column with a specific prefix in the database. For example, if the prefixes list contains "DP/",
        the function will search for a column named "DP/INFO" in addition to the default "INFO" column
        :type prefixes: list
        :return: a string that represents the name of the column found in the database, based on the
        input parameters. If the column is found, it returns the column name. If the column is not
        found, it returns None.
        """

        if not database:
            database = self.get_database()

        # Database columns
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
        if not column_found and "INFO" in database_columns :
            database_header_infos = self.get_header_infos_list(database=database)
            if column in database_header_infos:
                column_found = "INFO"

        return column_found
    

    def map_columns(self, database:str = None, table:str = None, columns:list = [], prefixes:list = ["INFO/"]) -> dict:
        """
        This function maps columns in a database table to their corresponding columns with specified
        prefixes.
        
        :param database: The name of the database to search for columns in. If no database is specified,
        the method will use the default database set in the connection
        :type database: str
        :param table: The name of the table in the database that you want to map the columns for
        :type table: str
        :param columns: A list of column names that you want to map to their corresponding column names
        in the database
        :type columns: list
        :param prefixes: The `prefixes` parameter is a list of strings that are used to filter the
        columns that are searched for. Only columns that start with one of the prefixes in the list will
        be considered. In the code above, the default value for `prefixes` is `["INFO/"]`, which
        :type prefixes: list
        :return: a dictionary that maps the input columns to their corresponding columns found in the
        specified database and table, with the specified prefixes.
        """

        if not database:
            database = self.get_database()

        # Init
        columns_mapping = {}

        for column in columns:
            column_found = self.find_column(database=database, table=table, column=column, prefixes=prefixes)
            columns_mapping[column] = column_found

        return columns_mapping
    

    def get_columns(self, database:str = None, table:str = None) -> list:
        """
        This function retrieves a list of columns from a specified database and table using SQL queries.
        
        :param database: The name of the database to get columns from. If None, it will use the current
        database
        :type database: str
        :param table: The name of the table in the database for which you want to retrieve the column
        names
        :type table: str
        :return: a list of column names for a given database and table. If no database is specified, it
        gets the current database. If the database exists and its format is known (e.g. "duckdb"), it
        connects to the database and retrieves the column names using a SQL query. Otherwise, it
        retrieves the column names using a SQL query on the current connection. If the table parameter
        is
        """
        
        if not database:
            database = self.get_database()

        if not table:
            table = self.get_database_table(database=database)

        if database and self.exists(database):
            database_format = self.get_format(database)
            if database_format in ["conn"]:
                if table:
                    database_conn = database
                    sql_query = f"SELECT * FROM {table} LIMIT 0"
                    columns_list = list(database_conn.query(sql_query).df().columns)
                    return columns_list
            elif database_format in ["duckdb"]:
                if table:
                    database_conn = duckdb.connect(database)
                    sql_query = f"SELECT * FROM {table} LIMIT 0"
                    columns_list = list(database_conn.query(sql_query).df().columns)
                    database_conn.close()
                    return columns_list
            elif database_format in ["sqlite"]:
                if table:
                    database_conn = sqlite3.connect(database)
                    sql_query = f"SELECT * FROM {table} LIMIT 0"
                    columns_list = list(pd.read_sql_query(sql_query, database_conn).columns)
                    return columns_list
            elif database_format in ["parquet", "vcf", "tsv", "csv", "tbl", "bed", "json"]:
                sql_from = self.get_sql_from(database)
                sql_query = f"SELECT * FROM {sql_from} LIMIT 0"
                return list(self.conn.query(sql_query).df().columns)

        return []


    def get_annotations(self, database:str = None) -> object:
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


    def get_extra_columns(self, database:str = None, database_type:str = None) -> list:
        """
        This function returns a list of extra columns in a database table that are not needed.
        
        :param database: A string representing the name of the database to retrieve columns from. If
        None, the default database will be used
        :type database: str
        :return: a list of extra columns that exist in the database table but are not needed based on
        the database type and the existing columns.
        """

        if not database:
            database = self.get_database()

        if not database:
            return []

        existing_columns = self.get_columns(database=database, table = self.get_database_table(database))
        if not database_type:
            database_type = self.get_type(database=database)
        needed_columns = self.get_needed_columns(database_columns=existing_columns, database_type=database_type)

        extra_columns = existing_columns.copy()

        for needed_col in needed_columns:
            if needed_columns.get(needed_col) in extra_columns:
                extra_columns.remove(needed_columns.get(needed_col))

        return extra_columns

       
    def is_vcf(self, database:str = None) -> bool:
        """
        This function checks if a given database is of type "vcf" by getting its columns and checking
        their types.
        
        :param database: The parameter `database` is a string that represents the name of the database
        that the function will use to check if the file is a VCF (Variant Call Format) file. If
        `database` is not provided, the function will use the default database
        :type database: str
        :return: a boolean value indicating whether the database type is "vcf" or not.
        """

        if not database:
            database = self.get_database()

        if not database:
            return False

        database_columns = self.get_columns(database=database, table=self.get_database_table(database))

        # Assume VCF is 8 needed columns, either only or with extra FORMAT column (assume other are samples)
        return self.get_type_from_columns(database_columns=database_columns, check_database_type="vcf") == "vcf" and ("FORMAT" in database_columns or len(database_columns) == 8)


    def get_conn(self):
        """
        The function returns the connection object.
        :return: The method is returning the value of the instance variable `self.conn`.
        """

        return self.conn
    

    def export(self, output_database:str, output_header:str = None, database:str = None, table:str = "variants") -> bool:
        """
        This function exports data from a database to a specified output format and compresses it if
        necessary.
        
        :param output_database: The path and filename of the output file to be exported
        :type output_database: str
        :param output_header: The parameter `output_header` is an optional string that represents the
        header of the output file. If it is not provided, the header will be automatically detected
        based on the output file format
        :type output_header: str
        :param database: The name of the database to export
        :type database: str
        :return: a boolean value indicating whether the export was successful or not.
        """

        if not database:
            database = self.get_database()

        if not database:
            return False
        
        # Header
        if output_header:
            self.get_header_file(header_file=output_header)

        # Auto-detect output type and compression and delimiter
        output_type = get_file_format(output_database)
        compressed = self.is_compressed(database=output_database)
        delimiter = FILE_FORMAT_DELIMITERS.get(output_type, "\t")
        
        # database type
        database_type = self.get_type(database=database)

        # Existing columns
        existing_columns = self.get_columns(database=database, table = self.get_database_table(database))

        # Extra columns
        extra_columns = self.get_extra_columns(database=database, database_type=output_type)
        
        # random
        random_tmp = ''.join(random.choice(string.ascii_lowercase) for i in range(10))

        # Query values
        default_empty_value = ""
        needed_columns = []
        query_export_format = None
        include_header = False

        # VCF
        if output_type in ["vcf"]:
            needed_columns = self.get_needed_columns(database_columns=existing_columns, database_type="vcf")
            if not self.is_vcf(database=database):
                extra_columns = []
            default_empty_value = "."
            query_export_format = f"FORMAT CSV, DELIMITER '{delimiter}', HEADER, QUOTE ''"
            include_header = True

        # TSV/CSV/TBL
        elif output_type in ["tsv", "csv", "tbl"]:
            needed_columns = self.get_needed_columns(database_columns=existing_columns, database_type=database_type)
            query_export_format = f"FORMAT CSV, DELIMITER '{delimiter}', HEADER"
            if delimiter in ["\t"]:
                include_header = True

        # JSON
        elif output_type in ["json"]:
            needed_columns = self.get_needed_columns(database_columns=existing_columns, database_type=database_type)
            query_export_format = f"FORMAT JSON, ARRAY TRUE"
            include_header = False

        # Parquet
        elif output_type in ["parquet"]:
            needed_columns = self.get_needed_columns(database_columns=existing_columns, database_type=database_type)
            query_export_format = f"FORMAT PARQUET"
            include_header = False

        # BED
        elif output_type in ["bed"]:
            needed_columns = self.get_needed_columns(database_columns=existing_columns, database_type="regions")
            query_export_format = f"FORMAT CSV, DELIMITER '{delimiter}', HEADER"
            include_header = True

        # duckDB
        elif output_type in ["duckdb"]:

            # Export database as Parquet
            database_export_parquet_file = f"""{output_database}.{random_tmp}.database_export.parquet"""
            self.export(database=database, output_database=database_export_parquet_file)
            
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

        # else:
        #     log.debug("Not available")

        # Construct query columns
        query_columns = []

        # Add Needed columns
        for needed_column in needed_columns:
            if needed_columns[needed_column]:
                query_column = f""" "{needed_columns[needed_column]}" """
            else:
                query_column = f""" '{default_empty_value}' """
            query_column_as = f""" "{needed_column}" """
            query_columns.append(f""" {query_column} AS {query_column_as} """)

        # Add Extra columns
        for extra_column in extra_columns:
            if extra_column not in needed_columns:
                query_columns.append(f""" "{extra_column}" AS "{extra_column}" """)

        # Query export columns
        query_export_columns = f""" {",".join(query_columns)} """

        if query_columns:

            # Compressed tmp file
            query_output_database_tmp = ""
            if not compressed and not include_header:
                query_output_database_tmp = output_database
            else:
                query_output_database_tmp = f"""{output_database}.{random_tmp}"""
            
            query_copy = f""" 
                COPY (
                    SELECT {query_export_columns}
                    FROM {self.get_sql_database_link(database=database)}
                    )
                TO '{query_output_database_tmp}'
                WITH ({query_export_format})
                """
            
            # Export
            self.query(database=database, query=query_copy)
            
            # Include header
            if include_header:
                # New tmp file
                query_output_database_header_tmp = f"""{query_output_database_tmp}.{random_tmp}"""
                # create tmp header file
                query_output_header_tmp = f"""{query_output_database_tmp}.header.{random_tmp}"""
                self.get_header_file(header_file=query_output_header_tmp, remove_header_line=True)
                # Concat header and database
                concat_file(input_files=[query_output_header_tmp, query_output_database_tmp], output_file=query_output_database_header_tmp)
                # move file
                shutil.move(query_output_database_header_tmp, query_output_database_tmp)
                # remove tmp
                remove_if_exists([query_output_header_tmp])

            # Compress
            if compressed:
                compress_file(input_file=query_output_database_tmp, output_file=output_database)
                # remove tmp
                remove_if_exists([query_output_database_tmp])
            else:
                shutil.move(query_output_database_tmp, output_database)

        return os.path.exists(output_database) and self.get_type(output_database)

