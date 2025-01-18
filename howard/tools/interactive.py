#!/usr/bin/env python

import logging as log
import re
import readline
from datetime import datetime
from tabulate import tabulate  # type: ignore
import pandas as pd  # type: ignore
import duckdb  # type: ignore
import tempfile
import signal
import contextlib
import io
from termcolor import colored  # type: ignore
from colorama import init  # type: ignore

from howard.functions.commons import prompt_color, prompt_mesage, prompt_line_color


def launch_interactive_terminal(
    args=None, variants=None, tmp=None, display_format="dataframe"
):
    """
    Launch an interactive SQL terminal with DuckDB
    """

    # Variants object
    if variants is None:
        log.warning(f"Variants not loaded")
        return None

    # Set up logging
    if tmp is None:
        tmp = tempfile.mkdtemp()

    # Get the DuckDB connection
    conn = variants.get_connexion()

    # Query limit
    query_limit = 1000

    try:
        # Execute query to test connexion
        conn.execute("SELECT 1")
    except Exception as e:
        msg_err = "Variants connexion failed"
        log.warning(msg_err)
        log.warning("Create empty connexion")
        conn = duckdb.connect(":memory:")

    # Create header table
    header_name = "header"
    log.debug(f"Loading table '{header_name}' as 'table'")
    log.info(f"Loading table '{header_name}'...")
    try:
        variants.load_header(drop=True, view_name=header_name)
    except:
        log.debug("View 'header' can not be loaded")

    # Create variants table
    view_name = "variants_view"
    view_type = "table"
    view_mode = "full"
    log.debug(f"Loading table '{view_name}' as '{view_type}' (mode '{view_mode}')")
    log.info(f"Loading table '{view_name}'...")
    try:
        variants.create_annotations_view(
            view=view_name,
            view_type=view_type,
            view_mode=view_mode,
            info_prefix_column="",
            fields_needed_all=True,
            info_struct_column="INFOS",
            sample_struct_column="SAMPLES",
            detect_type_list=True,
            drop_view=True,
        )
    except:
        log.warning(f"View '{view_name}' can not be created")

    # Print welcome message
    log.info("Interactive DuckDB SQL terminal")
    log.info("- 'exit' to quit.")
    log.info("- 'help' for a list of commands")

    # Configure readline for history
    histfile = f"{tmp}/.howard_interactive_history"
    try:
        readline.read_history_file(histfile)
    except FileNotFoundError:
        pass

    # Special commands
    special_commands = [
        "help",
        "tables",
        "history",
        "display",
        "limit",
        "python",
        "exit",
        "quit",
    ]

    # Configure readline for auto-completion
    keywords_sql = keywords = [
        "ADD",
        "ALL",
        "ALTER",
        "AND",
        "ANY",
        "AS",
        "ASC",
        "BETWEEN",
        "BIGINT",
        "BINARY",
        "BLOB",
        "BOOLEAN",
        "BY",
        "CASE",
        "CAST",
        "CHAR",
        "CHARACTER",
        "CHECK",
        "COLUMN",
        "CONSTRAINT",
        "CREATE",
        "CROSS",
        "CURRENT_DATE",
        "CURRENT_TIME",
        "CURRENT_TIMESTAMP",
        "DATABASE",
        "DATE",
        "DECIMAL",
        "DEFAULT",
        "DELETE",
        "DESC",
        "DESCRIBE",
        "DISTINCT",
        "DOUBLE",
        "DROP",
        "ELSE",
        "END",
        "ESCAPE",
        "EXISTS",
        "EXPLAIN",
        "FALSE",
        "FLOAT",
        "FOR",
        "FOREIGN",
        "FROM",
        "FULL",
        "GROUP",
        "HAVING",
        "IF",
        "ILIKE",
        "IN",
        "INNER",
        "INSERT",
        "INTEGER",
        "INTERSECT",
        "INTERVAL",
        "INTO",
        "IS",
        "JOIN",
        "LEFT",
        "LIKE",
        "LIMIT",
        "NATURAL",
        "NOT",
        "NULL",
        "NUMERIC",
        "ON",
        "OR",
        "ORDER",
        "OUTER",
        "PRIMARY",
        "REFERENCES",
        "RIGHT",
        "SELECT",
        "SET",
        "SHOW",
        "SMALLINT",
        "TABLE",
        "THEN",
        "TIME",
        "TIMESTAMP",
        "TRUE",
        "UNION",
        "UNIQUE",
        "UPDATE",
        "USING",
        "VALUES",
        "VARCHAR",
        "VIEW",
        "WHEN",
        "WHERE",
        "WITH",
    ]

    # Create keywords
    keywords = [f"{k.upper()} " for k in keywords_sql + special_commands]

    # tables and columns names
    tables_columns_names = []

    def get_table_and_column_names(reload: bool = False):
        """
        Get a list of table and column names from the database

        Returns:
            list: A list of table and column names

        """
        nonlocal tables_columns_names

        if not len(tables_columns_names) or reload:

            # Check tables
            tables = conn.execute("SHOW TABLES").fetchall()
            table_names = [table[0] for table in tables]

            # Create table.column keywords
            column_names = []
            for table in table_names:
                columns = conn.execute(f"DESCRIBE {table}").fetchall()
                for col in columns:
                    column_names.extend(
                        [
                            f"{table}.{col[0]}",
                        ]
                    )
            tables_columns_names = table_names + column_names

        return tables_columns_names

    def completer(text, state):
        """
        Auto-completion function for readline

        Args:
            text (str): The current input text
            state (int): The current state

        Returns:
            str: The next auto-completion option

        """
        global tables_columns_names

        # Get auto-completion options
        options = [
            i
            for i in keywords + get_table_and_column_names()
            if i.lower().startswith(text.lower())
        ]

        # If no options are found, return None
        if not options:
            return None

        # Find the longest common prefix
        def longest_common_prefix(strs):

            # Helper function to find the longest common prefix of a list of strings
            if not strs:
                return ""

            # Find the shortest string
            shortest = min(strs, key=len)

            # Compare the shortest string with the other strings
            for i, char in enumerate(shortest):
                for other in strs:
                    if other[i] != char:
                        return shortest[:i]
            return shortest

        # Get the longest common prefix of the options
        common_prefix = longest_common_prefix(options)

        # If the common prefix is longer than the input text, return it
        if state == 0:
            # If the common prefix is longer than the input text, return it
            if len(common_prefix) > len(text):
                return common_prefix
            else:
                return options[state]

        # If the common prefix is longer than the input text, return it
        elif state < len(options):
            return options[state]

        # If no more options are available, return None
        else:
            return None

    # Set up readline
    readline.set_completer(completer)
    readline.parse_and_bind("tab: complete")
    readline.set_auto_history(False)

    # Set up signal handler
    def signal_handler(sig, frame):
        print("\nQuery input cancelled. Type 'exit' to quit.")
        raise KeyboardInterrupt

    signal.signal(signal.SIGINT, signal_handler)

    def print_prompt(line: int = 0) -> str:
        init(autoreset=True)
        current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        prompt_log = colored(
            prompt_mesage.format(current_time, title="SQL"),
            color=prompt_color,
            attrs=[],
        )
        if line:
            prompt_line = " ++> "
        else:
            prompt_line = f" >>> "
        prompt_line = colored(prompt_line, color=prompt_line_color)

        colored_prompt = prompt_log + prompt_line

        return colored_prompt

    while True:
        try:
            # Get current date and time
            query_lines = []
            line = input(print_prompt(line=0))

            while True:

                if (
                    line.strip() == ""
                    or line.strip().endswith(";")
                    or line.strip().lower().startswith(tuple(special_commands))
                ):
                    query_lines.append(line.strip())
                    break
                query_lines.append(line.strip())
                line = input(print_prompt(line=1))
            query = " ".join(query_lines)

            # Add to history
            readline.add_history(query)

            # Clean query for special commands
            query_clean = query.replace(";", "").strip().lower()

            # Check for special commands
            ######

            # Check for special commands - exit or quit
            if query_clean in ["exit", "quit"]:
                break

            # Check for special commands - display
            elif query_clean.startswith("display"):
                display_format_input = query_clean.split(" ")[-1]
                if display_format_input not in [
                    "dataframe",
                    "markdown",
                    "simple",
                    "tabulate",
                ]:
                    print(
                        "Unknown display format. Please use 'dataframe', 'markdown', 'simple' or 'tabulate'."
                    )
                    continue
                display_format = display_format_input
                print(f"Display format set to: {display_format}")
                continue

            # Check for special commands - history
            elif query_clean.startswith("history"):
                # Get the history index, if any
                history_idx = query_clean.split(" ")[-1]
                # If no index is provided, show the last command
                history_as_query = None
                try:
                    # Try to convert the history index to an integer
                    history_idx = int(history_idx)
                    # Check if the index is valid
                    if (
                        history_idx <= 0
                        or history_idx > readline.get_current_history_length()
                    ):
                        print("Invalid history index")
                        continue
                    # Get the command at the specified index
                    history_as_query = readline.get_history_item(history_idx)
                    # Print the command
                    print(f"Query: {history_as_query}")
                except ValueError:
                    # If the history index is not an integer, check for special keywords
                    if history_idx.lower() in ["history", "all"]:
                        # Print all history
                        for i in range(readline.get_current_history_length()):
                            print(f"{i + 1}: {readline.get_history_item(i + 1)}")
                    else:
                        print("Invalid history command")
                        continue
                # If no history index is provided, show the last command
                if history_as_query is None:
                    continue
                else:
                    # Set the query to the history command
                    query = history_as_query

            # Check for special commands - help
            elif query_clean.startswith("help"):
                print("Available commands:")
                print(
                    "  history - Show command history (add a number to relaod a specific command)"
                )
                print(
                    "  display - Change display mode for query results (either 'tabulate' or 'dataframe')"
                )
                print(
                    "  limit - Change query limit for query results (e.g. 'limit 100')"
                )
                print("  tables - List all tables")
                print("  Ctrl+C - Cancel query input or execution")
                print("  exit, quit - Exit the terminal")
                continue

            # Check for special commands - tables
            elif query_clean == "tables":
                result = conn.execute("SHOW TABLES").fetchall()
                print(tabulate(result, headers=["Tables"], tablefmt="grid"))
                continue

            # Change query limit
            elif query_clean.startswith("limit"):
                query_limit_input = query_clean.split(" ")[-1]
                try:
                    query_limit_input = int(query_limit_input)
                except:
                    print(
                        f"Unknown limit value '{display_format_input}'. Please use integer value (e.g. 'limit {query_limit}')."
                    )
                    continue
                query_limit = query_limit_input
                print(f"Query limit set to: {query_limit}")
                continue

            # Python
            elif query_clean.startswith("python") and False:
                # print(f"query={query}")
                try:
                    exec_code = query[len("python") :].strip()
                    # Capture the output of the exec command
                    output = io.StringIO()
                    with contextlib.redirect_stdout(output), contextlib.redirect_stderr(
                        output
                    ):
                        try:
                            result = eval(exec_code)
                            if result is not None:
                                print(result)
                        except SyntaxError:
                            exec(exec_code)
                    result = output.getvalue()
                    print(result.strip())
                except Exception as e:
                    print(f"Error executing code: {e}")
                continue

            # Execute the query
            #####

            # Check if the query is a SELECT statement or a DESCRIBE statement
            if (
                query.strip().lower().startswith("select")
                or query.strip().lower().startswith("describe")
                or query.strip().lower().startswith("show")
            ):

                # Try limit query
                query_limited_check = False
                try:
                    query_limited = f"""
                    SELECT *
                    FROM
                        ({re.sub(";$", "", query.strip())})
                    LIMIT {query_limit}
                        """
                    log.debug(f"Try query limitation to {query_limit} lines...")
                    result = conn.execute(query_limited)
                    query_limited_check = True

                # Query without limitation
                except:
                    log.debug(f"Query limitation failed. Try without limitation...")
                    result = conn.execute(query)

                # Fetch rows
                rows = result.fetchmany(query_limit)

                # Fetch rows plus
                rows_plus = result.fetchmany(1)

                if rows_plus or (query_limited_check and len(rows) == query_limit):
                    msg_query_limit = f"Only {query_limit} lines shown (use 'limit' command to change)"
                else:
                    msg_query_limit = None
                # Get column names
                columns = [desc[0] for desc in result.description]

                # Print the results based on the display format
                if display_format in ["dataframe"]:
                    df = pd.DataFrame(rows, columns=columns)
                    print(df)
                elif display_format in ["markdown", "tabulate", "simple"]:
                    df = pd.DataFrame(rows, columns=columns)
                    if display_format in ["tabulate"]:
                        tablefmt = "grid"
                    elif display_format in ["simple"]:
                        tablefmt = "simple"
                    else:
                        tablefmt = "pipe"
                    print(df.to_markdown(tablefmt=tablefmt))
                else:
                    print(
                        "Unknown display format. Please use 'tabulate' or 'dataframe'."
                    )

                # Limit message
                if msg_query_limit:
                    log.warning(msg_query_limit)

            else:

                # Fetch All query
                result = conn.execute(query)
                rows = result.fetchall()

                # Print the number of affected rows for other types of queries
                if result.rowcount == -1:
                    print("Query executed successfully.")
                else:
                    print(f"{result.rowcount} rows affected.")

                # Reload tables and columns names for completer
                get_table_and_column_names(reload=True)

        except KeyboardInterrupt:
            # Handle the signal and return to prompt
            continue
        except Exception as e:
            # Print the error message
            print(f"Error: {e}")

    # Save history
    readline.write_history_file(histfile)
