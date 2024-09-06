"""
Tests

Usage:
pytest tests/

Coverage:
coverage run -m pytest tests/test_databases_from_annovar.py -x -v --log-cli-level=DEBUG --capture=tee-sys
coverage report --include=howard/* -m
"""

import argparse
import logging as log
import os
import re
from tempfile import TemporaryDirectory

from howard.functions.commons import (
    remove_if_exists,
    find_genome,
    extract_file,
    identical,
)
from howard.tools.tools import from_annovar
from test_needed import tests_config, tests_databases_folder


tests_folder = os.path.dirname(__file__)


def test_from_annovar_dbnsfp():

    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # Init files
        input_vcf = os.path.join(tests_databases_folder, "others", "hg19_dbnsfp42a.txt")
        output_vcf_expected = os.path.join(
            tests_databases_folder, "others", "hg19_dbnsfp42a.vcf"
        )
        output_vcf = os.path.join(tmp_dir, "hg19_dbnsfp42a.vcf.gz")
        output_parquet = os.path.join(tmp_dir, "hg19_dbnsfp42a.parquet")
        genomes_folder = tests_config["folders"]["databases"]["genomes"]
        genome = genomes_folder + "/hg19//hg19.fa"
        annovar_code = "nci60"
        config = {}

        # remove
        remove_if_exists([output_vcf])
        remove_if_exists([output_vcf + ".hdr"])
        remove_if_exists([output_parquet])
        remove_if_exists([output_parquet + ".hdr"])

        # Find genome
        genome = find_genome(genome)

        # prepare arguments for the query function
        args = argparse.Namespace(
            input_annovar=input_vcf,
            output_annovar=output_vcf,
            genome=genome,
            annovar_code=annovar_code,
            annovar_to_parquet=output_parquet,
            threads="2",
            annovar_reduce_memory="disable",
            annovar_multi_variant="disable",
            config=config,
        )

        # Process
        try:
            from_annovar(args)
            assert True
        except:
            assert False

        # Files exists
        assert os.path.exists(output_vcf)
        assert os.path.exists(output_vcf + ".hdr")
        assert os.path.exists(output_parquet)
        assert os.path.exists(output_parquet + ".hdr")

        # Files content
        output_vcf_uncompressed = re.sub(r"\.gz", "", output_vcf)
        extract_file(file_path=output_vcf, path=output_vcf_uncompressed)
        log.debug(f"[{output_vcf_expected}, {output_vcf_uncompressed}]")
        assert identical([output_vcf_expected, output_vcf_uncompressed])


def test_from_annovar_with_multiheader():

    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # Init files
        input_vcf = os.path.join(
            tests_databases_folder, "others", "hg19_nci60_with_multiheader.txt"
        )
        output_vcf_expected = os.path.join(
            tests_databases_folder, "others", "hg19_nci60_with_multiheader.vcf"
        )
        output_vcf = os.path.join(tmp_dir, "hg19_nci60_with_multiheader.vcf.gz")
        output_parquet = os.path.join(tmp_dir, "hg19_nci60_with_multiheader.parquet")
        genomes_folder = tests_config["folders"]["databases"]["genomes"]
        genome = genomes_folder + "/hg19//hg19.fa"
        annovar_code = "nci60"
        config = {}

        # remove
        remove_if_exists([output_vcf])
        remove_if_exists([output_vcf + ".hdr"])
        remove_if_exists([output_parquet])
        remove_if_exists([output_parquet + ".hdr"])

        # Find genome
        genome = find_genome(genome)

        # prepare arguments for the query function
        args = argparse.Namespace(
            input_annovar=input_vcf,
            output_annovar=output_vcf,
            genome=genome,
            annovar_code=annovar_code,
            annovar_to_parquet=output_parquet,
            threads="2",
            annovar_reduce_memory="disable",
            annovar_multi_variant="disable",
            config=config,
        )

        # Process
        try:
            from_annovar(args)
            assert True
        except:
            assert False

        # Files exists
        assert os.path.exists(output_vcf)
        assert os.path.exists(output_vcf + ".hdr")
        assert os.path.exists(output_parquet)
        assert os.path.exists(output_parquet + ".hdr")

        # Files content
        output_vcf_uncompressed = re.sub(r"\.gz", "", output_vcf)
        extract_file(file_path=output_vcf, path=output_vcf_uncompressed)
        log.debug(f"[{output_vcf_expected}, {output_vcf_uncompressed}]")
        assert identical([output_vcf_expected, output_vcf_uncompressed])


def test_from_annovar_with_header():

    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # Init files
        input_vcf = os.path.join(
            tests_databases_folder, "others", "hg19_nci60_with_header.txt"
        )
        output_vcf_expected = os.path.join(
            tests_databases_folder, "others", "hg19_nci60_with_header.vcf"
        )
        output_vcf = os.path.join(tmp_dir, "hg19_nci60_with_header.vcf.gz")
        output_parquet = os.path.join(tmp_dir, "hg19_nci60_with_header.parquet")
        genomes_folder = tests_config["folders"]["databases"]["genomes"]
        genome = genomes_folder + "/hg19//hg19.fa"
        annovar_code = "nci60"
        config = {}

        # remove
        remove_if_exists([output_vcf])
        remove_if_exists([output_vcf + ".hdr"])
        remove_if_exists([output_parquet])
        remove_if_exists([output_parquet + ".hdr"])

        # Find genome
        genome = find_genome(genome)

        # prepare arguments for the query function
        args = argparse.Namespace(
            input_annovar=input_vcf,
            output_annovar=output_vcf,
            genome=genome,
            annovar_code=annovar_code,
            annovar_to_parquet=output_parquet,
            threads="2",
            annovar_reduce_memory="disable",
            annovar_multi_variant="disable",
            config=config,
        )

        # Process
        try:
            from_annovar(args)
            assert True
        except:
            assert False

        # Files exists
        assert os.path.exists(output_vcf)
        assert os.path.exists(output_vcf + ".hdr")
        assert os.path.exists(output_parquet)
        assert os.path.exists(output_parquet + ".hdr")

        # Files content
        output_vcf_uncompressed = re.sub(r"\.gz", "", output_vcf)
        extract_file(file_path=output_vcf, path=output_vcf_uncompressed)
        log.debug(f"[{output_vcf_expected}, {output_vcf_uncompressed}]")
        assert identical([output_vcf_expected, output_vcf_uncompressed])


def test_from_annovar():

    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # Init files
        input_vcf = os.path.join(tests_databases_folder, "others", "hg19_nci60.txt")
        output_vcf_expected = os.path.join(
            tests_databases_folder, "others", "hg19_nci60.vcf"
        )
        output_vcf = os.path.join(tmp_dir, "hg19_nci60.vcf.gz")
        output_parquet = os.path.join(tmp_dir, "hg19_nci60.parquet")
        genomes_folder = tests_config["folders"]["databases"]["genomes"]
        genome = genomes_folder + "/hg19//hg19.fa"
        annovar_code = "nci60"
        config = {}

        # remove
        remove_if_exists([output_vcf])
        remove_if_exists([output_vcf + ".hdr"])
        remove_if_exists([output_parquet])
        remove_if_exists([output_parquet + ".hdr"])

        # Find genome
        genome = find_genome(genome)

        # prepare arguments for the query function
        args = argparse.Namespace(
            input_annovar=input_vcf,
            output_annovar=output_vcf,
            genome=genome,
            annovar_code=annovar_code,
            annovar_to_parquet=output_parquet,
            threads="2",
            annovar_reduce_memory="disable",
            annovar_multi_variant="disable",
            config=config,
        )

        # Process
        try:
            from_annovar(args)
            assert True
        except:
            assert False

        # Files exists
        assert os.path.exists(output_vcf)
        assert os.path.exists(output_vcf + ".hdr")
        assert os.path.exists(output_parquet)
        assert os.path.exists(output_parquet + ".hdr")

        # Files content
        output_vcf_uncompressed = re.sub(r"\.gz", "", output_vcf)
        extract_file(file_path=output_vcf, path=output_vcf_uncompressed)
        log.debug(f"[{output_vcf_expected}, {output_vcf_uncompressed}]")
        assert identical([output_vcf_expected, output_vcf_uncompressed])


def test_from_annovar_reduce_memory():

    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # Init files
        input_vcf = os.path.join(tests_databases_folder, "others", "hg19_nci60.txt")
        output_vcf_expected = os.path.join(
            tests_databases_folder, "others", "hg19_nci60.vcf"
        )
        output_vcf = os.path.join(tmp_dir, "hg19_nci60.vcf.gz")
        output_parquet = os.path.join(tmp_dir, "hg19_nci60.parquet")
        genomes_folder = tests_config["folders"]["databases"]["genomes"]
        genome = genomes_folder + "/hg19//hg19.fa"
        annovar_code = "nci60"
        config = {}

        # remove
        remove_if_exists([output_vcf])
        remove_if_exists([output_vcf + ".hdr"])
        remove_if_exists([output_parquet])
        remove_if_exists([output_parquet + ".hdr"])

        # Find genome
        genome = find_genome(genome)

        # prepare arguments for the query function
        args = argparse.Namespace(
            input_annovar=input_vcf,
            output_annovar=output_vcf,
            genome=genome,
            annovar_code=annovar_code,
            annovar_to_parquet=output_parquet,
            threads="2",
            annovar_reduce_memory="enable",
            annovar_multi_variant="disable",
            config=config,
        )

        # Process
        try:
            from_annovar(args)
            assert True
        except:
            assert False

        # Files exists
        assert os.path.exists(output_vcf)
        assert os.path.exists(output_vcf + ".hdr")
        assert os.path.exists(output_parquet)
        assert os.path.exists(output_parquet + ".hdr")

        # Files content
        output_vcf_uncompressed = re.sub(r"\.gz", "", output_vcf)
        extract_file(file_path=output_vcf, path=output_vcf_uncompressed)
        log.debug(f"[{output_vcf_expected}, {output_vcf_uncompressed}]")
        assert identical([output_vcf_expected, output_vcf_uncompressed])


def test_from_annovar_multi_variant():

    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # Init files
        input_vcf = os.path.join(tests_databases_folder, "others", "hg19_nci60.txt")
        output_vcf_expected = os.path.join(
            tests_databases_folder, "others", "hg19_nci60.vcf"
        )
        output_vcf = os.path.join(tmp_dir, "hg19_nci60.vcf.gz")
        output_parquet = os.path.join(tmp_dir, "hg19_nci60.parquet")
        genomes_folder = tests_config["folders"]["databases"]["genomes"]
        genome = genomes_folder + "/hg19//hg19.fa"
        annovar_code = "nci60"
        config = {}

        # remove
        remove_if_exists([output_vcf])
        remove_if_exists([output_vcf + ".hdr"])
        remove_if_exists([output_parquet])
        remove_if_exists([output_parquet + ".hdr"])

        # Find genome
        genome = find_genome(genome)

        # prepare arguments for the query function
        args = argparse.Namespace(
            input_annovar=input_vcf,
            output_annovar=output_vcf,
            genome=genome,
            annovar_code=annovar_code,
            annovar_to_parquet=output_parquet,
            threads="2",
            annovar_reduce_memory="disable",
            annovar_multi_variant="enable",
            config=config,
        )

        # Process
        try:
            from_annovar(args)
            assert True
        except:
            assert False

        # Files exists
        assert os.path.exists(output_vcf)
        assert os.path.exists(output_vcf + ".hdr")
        assert os.path.exists(output_parquet)
        assert os.path.exists(output_parquet + ".hdr")

        # Files content
        output_vcf_uncompressed = re.sub(r"\.gz", "", output_vcf)
        extract_file(file_path=output_vcf, path=output_vcf_uncompressed)
        log.debug(f"[{output_vcf_expected}, {output_vcf_uncompressed}]")
        assert identical([output_vcf_expected, output_vcf_uncompressed])


def test_from_annovar_reduce_memory_multi_variant():

    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # Init files
        input_vcf = os.path.join(tests_databases_folder, "others", "hg19_nci60.txt")
        output_vcf_expected = os.path.join(
            tests_databases_folder, "others", "hg19_nci60.vcf"
        )
        output_vcf = os.path.join(tmp_dir, "hg19_nci60.vcf.gz")
        output_parquet = os.path.join(tmp_dir, "hg19_nci60.parquet")
        genomes_folder = tests_config["folders"]["databases"]["genomes"]
        genome = genomes_folder + "/hg19//hg19.fa"
        annovar_code = "nci60"
        config = {}

        # remove
        remove_if_exists([output_vcf])
        remove_if_exists([output_vcf + ".hdr"])
        remove_if_exists([output_parquet])
        remove_if_exists([output_parquet + ".hdr"])

        # Find genome
        genome = find_genome(genome)

        # prepare arguments for the query function
        args = argparse.Namespace(
            input_annovar=input_vcf,
            output_annovar=output_vcf,
            genome=genome,
            annovar_code=annovar_code,
            annovar_to_parquet=output_parquet,
            threads="2",
            annovar_reduce_memory="enable",
            annovar_multi_variant="enable",
            config=config,
        )

        # Process
        try:
            from_annovar(args)
            assert True
        except:
            assert False

        # Files exists
        assert os.path.exists(output_vcf)
        assert os.path.exists(output_vcf + ".hdr")
        assert os.path.exists(output_parquet)
        assert os.path.exists(output_parquet + ".hdr")

        # Files content
        output_vcf_uncompressed = re.sub(r"\.gz", "", output_vcf)
        extract_file(file_path=output_vcf, path=output_vcf_uncompressed)
        log.debug(f"[{output_vcf_expected}, {output_vcf_uncompressed}]")
        assert identical([output_vcf_expected, output_vcf_uncompressed])
