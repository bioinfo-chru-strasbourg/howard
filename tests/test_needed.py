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

from howard.functions.commons import *
from howard.objects.variants import Variants
from howard.functions.databases import *
from howard.tools.tools import main_folder, arguments_dict


# Main tests folder
tests_folder = os.path.dirname(__file__)
tests_data_folder = tests_folder + "/data"
tests_databases_folder = tests_folder + "/databases"
tests_databases_release = "current"
tests_annotations_folder = tests_databases_folder + f"/annotations/{tests_databases_release}/hg19"

# Tools folder
#tests_tools = "~/howard/tools"
tests_tools = full_path("~/howard/tools")

setup_cfg = f'{main_folder}/../setup.cfg'


# Test config
tests_config = {
  "threads": 2,
  "memory": None,
  "folders": {
    "databases": {
      "root": "",
      "annotations": [f"{tests_databases_folder}/annotations/{tests_databases_release}"],
      "annovar": f"{tests_databases_folder}/annovar/{tests_databases_release}",
      "snpeff": f"{tests_databases_folder}/snpeff/{tests_databases_release}",
      "varank": f"{tests_databases_folder}/varank/{tests_databases_release}",
      "genomes": f"{tests_databases_folder}/genomes/{tests_databases_release}",
      "refseq": f"{tests_databases_folder}/refseq/{tests_databases_release}"
    }
  },
  "tools": {
    "bcftools": "bcftools",
    "bgzip": "bgzip",
    "java": "java",
    "snpeff": f"{tests_tools}/snpeff/current/bin/snpEff.jar",
    "annovar": f"{tests_tools}/annovar/current/bin/table_annovar.pl"
  },
  "verbosity": "debug"
}


# Annotation databases
database_files = {
    "parquet" : tests_annotations_folder + "/nci60.parquet",
    "partition_parquet" : tests_annotations_folder + "/nci60.partition.parquet",
    "parquet_without_header" : tests_annotations_folder + "/nci60.without_header.parquet",
    "duckdb" : tests_annotations_folder + "/nci60.duckdb",
    "duckdb_no_annotation_table" : tests_annotations_folder + "/nci60.no_annotation_table.duckdb",
    "sqlite" : tests_annotations_folder + "/nci60.sqlite",
    "vcf" : tests_annotations_folder + "/nci60.vcf",
    "vcf_gz" : tests_annotations_folder + "/nci60.vcf.gz",
    "partition_vcf_gz" : tests_annotations_folder + "/nci60.partition.vcf.gz",
    "vcf_without_header" : tests_annotations_folder + "/nci60.without_header.vcf",
    "vcf_gz_without_header" : tests_annotations_folder + "/nci60.without_header.vcf.gz",
    "tsv" : tests_annotations_folder + "/nci60.tsv",
    "tsv_alternative_columns" : tests_annotations_folder + "/nci60.alternative_columns.tsv",
    "tsv_failed_columns" : tests_annotations_folder + "/nci60.failed_columns.tsv",
    "tsv_lower_columns" : tests_annotations_folder + "/nci60.lower_columns.tsv",
    "tsv_without_header" : tests_annotations_folder + "/nci60.without_header.tsv",
    "tsv_variants" : tests_annotations_folder + "/nci60.variants.tsv",
    "tsv_gz" : tests_annotations_folder + "/nci60.tsv.gz",
    "csv" : tests_annotations_folder + "/nci60.csv",
    "csv_gz" : tests_annotations_folder + "/nci60.csv.gz",
    "tbl" : tests_annotations_folder + "/nci60.tbl",
    "tbl_gz" : tests_annotations_folder + "/nci60.tbl.gz",
    "json" : tests_annotations_folder + "/nci60.json",
    "json_gz" : tests_annotations_folder + "/nci60.json.gz",
    "bed" : tests_annotations_folder + "/annotation_regions.bed",
    "bed_gz" : tests_annotations_folder + "/annotation_regions.bed.gz",
    "refgene" : tests_annotations_folder + "/refGene.bed",
    "refgene_without_header" : tests_annotations_folder + "/refGene.without_header.bed",
    "refgene_gz" : tests_annotations_folder + "/refGene.bed.gz",
    "example_vcf" : tests_data_folder + "/example.vcf",
    "example_vcf_gz" : tests_data_folder + "/example.vcf.gz",
    "example_vcf_gzip" : tests_data_folder + "/example.vcf.gzip",
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
    refseq_folder = config["folders"]["databases"]["refseq"]
    databases_download_refseq(assemblies=[assembly], refseq_folder=refseq_folder, refseq_files=refseq_files)

    # Genome
    genome_folder = config["folders"]["databases"]["genomes"]
    provider = "UCSC"
    contig_regex = None
    threads = 1
    databases_download_genomes(assemblies=[assembly], genomes_folder=genome_folder, provider=provider, contig_regex=contig_regex, threads=threads)

    # Annovar
    annovar_folder = config["folders"]["databases"]["annovar"]
    annovar_files = ["nci60"]
    databases_download_annovar(folder=annovar_folder, files=annovar_files, assemblies=[assembly])

