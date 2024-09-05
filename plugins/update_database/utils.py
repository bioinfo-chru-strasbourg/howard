import gzip
import shutil
import os
import hashlib
import re
from functools import lru_cache


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


def get_md5(file: str):
    """
    get md5 of a file
    """
    return hashlib.md5(open(file, "rb").read()).hexdigest()


@lru_cache
def get_compiled_pattern(pattern):
    return re.compile(pattern)
