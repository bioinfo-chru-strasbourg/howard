#!/usr/bin/env python

import logging as log
import readline
from datetime import datetime
from tabulate import tabulate  # type: ignore
import pandas as pd  # type: ignore
import duckdb  # type: ignore
import tempfile
import signal
import contextlib
import io


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

    try:
        # Execute query to test connexion
        conn.execute("SELECT 1")
    except Exception as e:
        msg_err = "Variants connexion failed"
        log.warning(msg_err)
        log.warning("Create empty connexion")
        conn = duckdb.connect(":memory:")

    # Create variants view
    try:
        variants.create_annotations_view(
            view="variants_view",
            info_prefix_column="",
            fields_needed_all=True,
            info_struct_column="INFOS",
            detect_type_list=True,
        )
    except:
        log.debug("View can not be created")

    variants.load_header(drop=True)
    try:
        variants.load_header(drop=True)
    except:
        log.debug("View can not reloaded")

    # Print welcome message
    log.info(f"Interactive DuckDB SQL terminal. Type 'exit' to quit.")
    log.info("Type 'help' for a list of commands.")

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
    keywords = [f"{k} " for k in keywords_sql] + special_commands

    def get_table_and_column_names():
        """
        Get a list of table and column names from the database

        Returns:
            list: A list of table and column names

        """

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
        return table_names + column_names

    def completer(text, state):
        """
        Auto-completion function for readline

        Args:
            text (str): The current input text
            state (int): The current state

        Returns:
            str: The next auto-completion option

        """

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
        print("\nQuery input cancelled. Type 'exit' or 'Ctrl+C' again to quit.")
        raise KeyboardInterrupt

    signal.signal(signal.SIGINT, signal_handler)

    while True:
        try:
            # Get current date and time
            current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            query_lines = []
            line = input(f"#[{current_time}] [SQL]> ")
            while True:

                if (
                    line.strip() == ""
                    or line.strip().endswith(";")
                    or line.strip().lower().startswith(tuple(special_commands))
                ):
                    query_lines.append(line.strip())
                    break
                query_lines.append(line.strip())
                line = input(f"#[{current_time}] [SQL]| ")
            query = " ".join(query_lines)

            # Add to history
            readline.add_history(query)

            # Check for special commands

            # Check for special commands - exit or quit
            if query.lower() in ["exit", "quit"]:
                break

            # Check for special commands - display
            elif query.lower().startswith("display"):
                display_format_input = query.lower().split(" ")[-1]
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
            elif query.lower().startswith("history"):
                # Get the history index, if any
                history_idx = query.lower().split(" ")[-1]
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
            elif query.lower() == "help":
                print("Available commands:")
                print(
                    "  history - Show command history (add a number to relaod a specific command)"
                )
                print(
                    "  display - Change display mode for query results (either 'tabulate' or 'dataframe')"
                )
                print("  tables - List all tables")
                print("  Ctrl+C - Cancel query input or execution")
                print("  exit, quit - Exit the terminal")
                continue

            # Check for special commands - tables
            elif query.lower() == "tables":
                result = conn.execute("SHOW TABLES").fetchall()
                print(tabulate(result, headers=["Tables"], tablefmt="grid"))
                continue

            elif query.lower().startswith("python"):
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
            result = conn.execute(query)

            # Fetch all results
            rows = result.fetchall()

            # Get column names
            columns = [desc[0] for desc in result.description]

            # Check if the query is a SELECT statement or a DESCRIBE statement
            if query.strip().lower().startswith(
                "select"
            ) or query.strip().lower().startswith("describe"):
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
            else:
                # Print the number of affected rows for other types of queries
                if result.rowcount == -1:
                    print("Query executed successfully.")
                else:
                    print(f"{result.rowcount} rows affected.")

        except KeyboardInterrupt:
            # Handle the signal and return to prompt
            continue
        except Exception as e:
            # Print the error message
            print(f"Error: {e}")

    # Save history
    readline.write_history_file(histfile)
