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
import json
import subprocess


def count_row_file(file):
    print("Checking number of rows")
    result = subprocess.run(
        ["bash", "-c", f"zcat {file} | wc -l"],
        check=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    line_count = int(result.stdout.strip())
    return line_count


def read_json(configfile: str) -> dict:
    """
    From json file to python dict
    :param configfile: path of json file
    :return: python dict
    """
    with open(configfile) as js:
        data = json.load(js)
        return data


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


def metaheader_rows(
    fields, description, id=None, number=None, type=None, quoting=False
) -> str:
    """
    From metaheader information to VCF header row
    ##INFO=<ID=STRAND,Number=1,Type=String,Description="Gene strand">

    :param fields: either INFO, FORMAT, FILTER, ALT from vcf
    :param description: Description from metaheader file, could also be the assembly name   for contig
    :param id: The ID of the corresponding field, which is an annotation in the INFO field  of the vcf
    :param number: Number of value for the field, could be 0 for flag, 1, A (match number   of allele) or "." for a list
    :param type: Type of Value Float, Integer, String.. conf vcf specs
    :param quoting: double quote is mandatory for the description metaheader, if there is   no double quote in description string set this param
    :return: processed row of metaheader
    """
    if quoting and description is not None:
        description = '"' + description + '"'
    if fields in ["INFO", "FORMAT"]:
        keys = ["ID", "Number", "Type", "Description"]
        values = list(map(str, [id, number, type, description]))
        return (
            "##"
            + fields
            + "=<"
            + ",".join(["=".join(val) for val in list(zip(keys, values))])
            + ">"
        )
    elif fields in ["ALT", "FILTER"]:
        keys = ["ID", "Description"]
        values = list(map(str, [id, description]))
        return (
            "##"
            + fields
            + "=<"
            + ",".join(["=".join(val) for val in list(zip(keys, values))])
            + ">"
        )
    elif fields == "contig":
        keys = ["ID", "assembly", "length"]
        values = list(map(str, [id, description, number]))
        return (
            "##"
            + fields
            + "=<"
            + ",".join(["=".join(val) for val in list(zip(keys, values))])
            + ">"
        )
    else:
        return "##" + fields + "=" + description


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
