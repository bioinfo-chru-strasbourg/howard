# -*- coding: utf-8 -*-
"""
Tests

Usage:
pytest tests/

Coverage:
coverage run -m pytest . -x -v
coverage report --include=howard/* -m
"""

import logging as log
import os
import sys
from tempfile import TemporaryDirectory
import duckdb
import re
import Bio.bgzf as bgzf
import gzip
import pytest

from howard.commons import *
from howard.objects.variants import Variants
from howard.tools.databases import *


# Main tests folder
tests_folder = os.path.dirname(__file__)
tests_data_folder = tests_folder + "/data"
tests_annotations_folder = tests_folder + "/data/annotations"

# Tools folder
tests_tools = "/tools"

# Test config
tests_config = {
  "threads": 1,
  "memory": None,
  "verbosity": "warning",
  "folders": {
    "databases": {
      "root": "",
      "parquet": [f"{tests_folder}/data/annotations"],
      "bcftools": [f"{tests_folder}/data/annotations"],
      "annovar": f"{tests_folder}/data/annotations/annovar",
      "snpeff": f"{tests_folder}/data/annotations/snpeff",
      "varank": f"{tests_folder}/data/annotations/varank",
      "genomes": f"{tests_folder}/data/annotations/genomes",
      "refGene": f"{tests_folder}/data/annotations/refseq",
      "refSeqLink": f"{tests_folder}/data/annotations/refseq",
    }
  },
  "tools": {
    "bcftools": {"bin": "bcftools"},
    "bgzip": {"bin": "bgzip"},
    "snpeff": {"jar": f"{tests_tools}/snpeff/current/bin/snpEff.jar"},
    "java": {"bin": "/usr/bin/java"},
    "annovar": {"bin": f"{tests_tools}/annovar/current/bin/table_annovar.pl"}
  }
}



def download_needed_databases():
    """
    The function `download_needed_databases` downloads needed databases for a specified assembly and
    stores them in the appropriate folders.
    """

    # Construct config dict
    config = tests_config.copy()

    # assembly
    assembly = "hg19"

    # refSeq
    refseq_files = ['ncbiRefSeq.txt', 'ncbiRefSeqLink.txt']
    refseq_folder = config["folders"]["databases"]["refGene"]
    databases_download_refseq(assemblies=[assembly], refseq_folder=refseq_folder, refseq_files=refseq_files)

    # Genome
    genome_folder = config["folders"]["databases"]["genomes"]
    provider = "UCSC"
    contig_regex = None
    threads = 1
    databases_download_genomes(assemblies=[assembly], genome_folder=genome_folder, provider=provider, contig_regex=contig_regex, threads=threads)

