
# import pyarrow as pa
# import pyarrow.csv as csv
import polars as pl
import pandas as pd
import duckdb
import vcf

from howard.commons import *
from howard.tools.databases import *


SEP_TYPE = {
    "vcf" : "\t",
    "tsv" : "\t",
    "csv" : ",",
    "psv" : "|",
    "bed" : "\t",
}

DATABASE_TYPE_NEEDED_COLUMNS = {
    "variants":
        {
            "chromosome": ["#CHROM", "CHROM", "CHR", "CHROMOSOME"],
            "position": ["POS"],
            "reference": ["REF"],
            "alternative": ["ALT"],
        },
    "regions":
        {
            "chromosome": ["#CHROM", "CHROM", "CHR", "CHROMOSOME"],
            "start": ["START", "POSITIONSTART"],
            "end": ["END", "POSITIONEND"],
        },
    "vcf":
        {
            "chromosome": ["#CHROM", "CHROM", "CHR", "CHROMOSOME"],
            "position": ["POS"],
            "reference": ["REF"],
            "alternative": ["ALT"],
            "info": ["INFO"],
        },
}

class Database:

    def __init__(self, database:str = None, header:str = None, databases_folders:list = None, header_file:str = None) -> None:
        """
        This is the initialization function for a class that sets up a database and header file for use
        in a DuckDB connection.
        
        :param database: A string representing the name of the database to be used. If None, the default
        database will be used
        :type database: str
        :param header: The header parameter is a string that represents the name of the header file that
        contains the column names for the database. It is used in conjunction with the database
        parameter to set the header for the database. If the header parameter is not provided, the
        header will be set to None
        :type header: str
        :param databases_folders: The `databases_folders` parameter is a list of folders where the
        database files are located. This parameter is used in the `set_database()` method to search for
        the database file in the specified folders. If the database file is not found in any of the
        folders, an error is raised
        :type databases_folders: list
        :param header_file: The header_file parameter is a string that represents the file path to the
        header file that contains the column names for the database. This parameter is used in the
        set_header method to set the header attribute of the class
        :type header_file: str
        """

        self.database = None
        self.header = None

        self.set_database(database=database, databases_folders=databases_folders)
        self.set_header(database=database, header_file=header_file)
        self.conn = duckdb.connect()


    def set_database(self, database:str, databases_folders:list = None) -> str:
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
        """

        if database is None:
            self.database = None
        elif self.exists(database=database):
            self.database = database
        elif self.find_database(database=database, databases_folders=databases_folders):
            self.database = self.find_database(database=database, databases_folders=databases_folders)


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

        default_header_list = [
            '##fileformat=VCFv4.2',
            '#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO'
        ]

        if not header_list:
            header_list = default_header_list
            
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

        default_header_list = [
            '##fileformat=VCFv4.2',
            '#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO'
            ]

        if not header_file:
            header_list = default_header_list
        else:
            header_list = self.read_header_file(header_file=header_file)

        return self.get_header_from_list(header_list)
        

    def get_header(self, header_file:str = None, header_list:list = None) -> vcf:
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
            return self.header
        elif header_file:
            return self.get_header_from_file(header_file)
        elif header_list:
            return self.get_header_from_list(header_list)
        else:
            return None


    def set_header(self, database:str = None, header_file:str = None) -> None:
        """
        This function sets the header of a database based on a provided header file or the database
        format.
        
        :param database: A string representing the name or path of a database file. If not provided, the
        method will attempt to get the database name from the object's attributes
        :type database: str
        :param header_file: A string representing the file path of a header file. If provided, the
        function will use this header file to set the header attribute of the object. If not provided,
        the function will try to determine the header from the database file
        :type header_file: str
        """

        if header_file and os.path.exists(header_file):

            # header provided
            self.header = self.get_header(header_file=header_file)
        
        else:

            if not database:
                database = self.get_database()

            if database:

                # database format
                database_format = self.get_format(database=database)

                # extra header file
                database_header_file = database + ".hdr"

                # header in extra header file
                if os.path.exists(database_header_file):
                    self.header = self.get_header(header_file=database_header_file)
                # header within file
                elif database_format in ["vcf", "tsv", "csv", "psv", "bed"]:
                    self.header = self.get_header(header_file=database)
                # Not header
                else:
                    self.header = self.get_header()

            else:

                self.header = None


    def find_database(self, database:str = None, databases_folders:list = None) -> str:
        """
        This function finds a database file in a specified folder or the current directory.
        
        :param database: The name of the database to be searched for. If not provided, it will call the
        `get_database()` method to get the name of the database
        :type database: str
        :param databases_folders: A list of folders where the function should look for the database
        file. If this parameter is not provided, the function will look for the database file in the
        current directory
        :type databases_folders: list
        :return: a string that represents the path to the database file. If the database is not found or
        if no database is specified, it returns None.
        """
        
        if not database:
            database = self.get_database()

        if not database:
            return None

        elif self.exists(database=database):
            return database
        
        else:

            if not databases_folders:
                databases_folders = ['.']

            database_file = None
            for databases_folder in databases_folders:
                
                # Log
                log.debug("Annotation file check: " + str(databases_folder+"/"+database))

                # In folder
                if os.path.exists(databases_folder+"/"+database):
                    database_file = databases_folder+"/"+database

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
        
        return database and os.path.exists(database)


    def get_format(self, database:str = None) -> str:
        """
        This Python function returns the file format of a given database or the current database if none
        is provided.
        Format database:
            - parquet
            - duckdb
            - vcf
            - csv
        
        :param database: The `database` parameter is a string that represents the type of database. It
        is an optional parameter and if not provided, the function will call the `get_database()` method
        to retrieve the database type
        :type database: str
        :return: a string that represents the type of database. The type of database can be one of the
        following: "parquet", "duckdb", "vcf", or "csv". The specific type of database is determined by
        the input parameter `database`, which is either passed as an argument to the function or
        obtained by calling the `get_database()` method. The `get_file_format
        """

        if not database:
            database = self.get_database()
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
            if format in ["duckdb"]:
                database_conn = duckdb.connect(database)
                database_tables = list(database_conn.query("SHOW TABLES;").df()["name"])
                database_conn.close()
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
        
        if not database:
            database = self.get_database()

        if database and self.exists(database):
            
            database_format = self.get_format(database)
            if database_format in ["duckdb"]:
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
        :param check_database_type: A list of database types to check for. If not provided, it defaults
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


        for needed_col in needed_columns:
            variants_columns_found[needed_col] = None
            for possible_col in needed_columns[needed_col]:
                if database_columns:
                    for existing_column in database_columns:
                        if possible_col.upper() == existing_column.upper():
                            variants_columns_found[needed_col] = existing_column
                            break

        return variants_columns_found


    def get_sql_from(self, database:str = None):
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
        if database_format in ["parquet"]:
            sql_form = f"read_parquet('{database}')"
        elif database_format in ["vcf", "tsv", "csv", "psv", "bed"]:
            delimiter = SEP_TYPE.get(database_format,"\t")
            sql_form = f"read_csv('{database}', auto_detect=True, delim='{delimiter}')"
        elif database_format in ["json"]:
            sql_form = f"read_json('{database}', auto_detect=True)"
        elif database_format in ["duckdb"]:
            sql_form = f"'{database}'"

        return sql_form
    

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

        return get_file_compressed(database)
    

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

        if database and self.exists(database):
            database_format = self.get_format(database)
            if database_format in ["duckdb"]:
                if table:
                    database_conn = duckdb.connect(database)
                    sql_query = f"SELECT * FROM {table} LIMIT 0"
                    columns_list = list(database_conn.query(sql_query).df().columns)
                    database_conn.close()
                    return columns_list
                else:
                    return []
            elif database_format in ["parquet", "vcf", "tsv", "csv", "psv", "bed", "json"]:
                sql_from = self.get_sql_from(database)
                sql_query = f"SELECT * FROM {sql_from} LIMIT 0"
                return list(self.conn.query(sql_query).df().columns)
            else:
                return []
        return []


    def get_extra_columns(self, database:str = None) -> list:
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
        return self.get_type_from_columns(database_columns=database_columns, check_database_type="vcf") == "vcf"





