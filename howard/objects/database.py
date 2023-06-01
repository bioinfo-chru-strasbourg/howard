
# import pyarrow as pa
# import pyarrow.csv as csv
import polars as pl
import pandas as pd
import duckdb

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

    def __init__(self, database:str = None) -> None:
        self.set_database(database)
        self.conn = duckdb.connect()


    def set_database(self, database:str = None) -> str:
        self.database = database
    
    
    def get_database(self) -> str:
        return self.database


    def exists(self, database:str = None) -> bool:
        if not database:
            database = self.get_database()
        return database and os.path.exists(database)


    def get_format(self, database:str = None) -> str:
        """
        This function get type of database:
            - parquet
            - duckdb
            - vcf
            - csv
        """
        if not database:
            database = self.get_database()
        return get_file_format(database)


    def get_type(self, database:str = None) -> str:
        """
        This function get type of database:
            - variants (VCF-like)
            - regions (BED-like)
        """
        if not database:
            database = self.get_database()
        # print("")
        # print("database: " + str(database))
        if database and self.exists(database):
            
            # format = self.get_format(database)
            # if format in ["duckdb"]:
            #     database_columns = self.get_columns(database, table=self.get_database_table(database))
            # else:
            #     database_columns = self.get_columns(database)
            
            database_columns = self.get_columns(database, table=self.get_database_table(database))

            #print(self.get_extra_columns(database=database))
            
            return self.get_type_from_columns(database_columns)

        else:
            return None
    

    def get_database_tables(self, database:str = None) -> str:
        """

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
            sql_form = f"read_json('{database}')"
        elif database_format in ["duckdb"]:
            sql_form = f"'{database}'"

        return sql_form
    

    def is_compressed(self, database:str = None) -> bool:
        """
        It returns if the input file is compressed.
        :return: The format is being returned.
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
            if database_format == "unknown":
                return []
            elif database_format in ["duckdb"]:
                if table:
                    database_conn = duckdb.connect(database)
                    sql_query = f"SELECT * FROM {table} LIMIT 0"
                    columns_list = list(database_conn.query(sql_query).df().columns)
                    database_conn.close()
                    return columns_list
                else:
                    return []
            else:
                sql_from = self.get_sql_from(database)
                sql_query = f"SELECT * FROM {sql_from} LIMIT 0"
                return list(self.conn.query(sql_query).df().columns)
        
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





