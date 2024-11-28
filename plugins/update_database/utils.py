import gzip
import shutil
import os
import hashlib
import re
from functools import lru_cache
from datetime import datetime
from pathlib import Path
import time
import logging as log


def recursive_chmod(directory, mode):
    directory = Path(directory)  # Convert to a Path object
    # Change permission for the directory and all subdirectories and files
    for path in directory.rglob("*"):
        path.chmod(mode)
    # Finally, change permission for the root directory itself
    directory.chmod(mode)


def now():
    current_date = datetime.now()
    return current_date.strftime("%Y%m%d")


def metaheader_rows(fields, id, number, type, description):
    """
    ##INFO=<ID=STRAND,Number=1,Type=String,Description="Gene strand">
    """
    keys = ["ID", "Number", "Type", "Description"]
    values = list(map(str, [id, number, type, '"' + description + '"']))
    return (
        "##"
        + fields
        + "=<"
        + ",".join(["=".join(val) for val in list(zip(keys, values))])
        + ">"
    )


def extract_gz_file(input_path: str, output_path: str) -> str:
    """
    unzip .gz file
    """
    with gzip.open(input_path, "rb") as f_in:
        with open(output_path, "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)
    print(f"Extracted '{input_path}' to '{output_path}'.")
    os.remove(input_path)
    return output_path


def read_md5_file(file: str):
    return open(file).readline().split()[0]


def get_md5(file: str):
    """
    get md5 of a file
    """
    return hashlib.md5(open(file, "rb").read()).hexdigest()


@lru_cache
def get_compiled_pattern(pattern):
    return re.compile(pattern)


def find_files(path: str, prefix=None, suffix=None) -> list:
    """
    Find files in a specific folder with the given prefix, suffix, or both efficiently.

    :param path: The path of the folder where to search for files.
    :param prefix: The prefix to filter files (optional).
    :param suffix: The suffix to filter files (optional).
    :return: A list of file names that match the given prefix and/or suffix with full path.
    """
    if not os.path.exists(path):
        raise FileNotFoundError(f"The folder {path} does not exist.")

    if not os.path.isdir(path):
        raise NotADirectoryError(f"The path {path} is not a directory.")

    # List comprehension for filtering files based on prefix and/or suffix
    matching_files = [
        os.path.join(path, file_name)
        for file_name in os.listdir(path)
        if (not prefix or file_name.startswith(prefix))
        and (not suffix or file_name.endswith(suffix))
    ]

    return matching_files


def timeit(func):
    """
    Decorator that measures the execution time of a function.

    :param func: Function to measure.
    :return: Wrapped function with timing.
    """

    def wrapper(*args, **kwargs):
        start_time = time.time()
        result = func(*args, **kwargs)
        end_time = time.time()
        execution_time = end_time - start_time
        log.debug(
            f"Function '{func.__name__}' executed in {execution_time:.2f} seconds."
        )
        return result

    return wrapper
