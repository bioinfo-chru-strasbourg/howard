import gzip
import shutil
import os
import hashlib


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
