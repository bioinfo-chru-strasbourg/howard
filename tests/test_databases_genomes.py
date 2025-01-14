# -*- coding: utf-8 -*-
"""
Tests

Usage:
pytest tests/

Coverage:
coverage run -m pytest tests/test_databases.py -x -vv --log-cli-level=DEBUG --capture=tee-sys
coverage report --include=howard/* -m
"""

import os
from tempfile import TemporaryDirectory

from howard.functions.commons import DEFAULT_GENOME_FOLDER
from howard.functions.databases import databases_download_genomes

# from howard.functions.commons import *
# from howard.tools.databases import *

# from test_needed import *


def test_databases_download_genomes():
    """
    The function tests the databases_download_genomes function by checking if genomes are downloaded correctly for
    different assemblies and contig filters.
    """

    import genomepy  # type: ignore

    # Init
    assemblies_config = {
        "sacCer3": {
            "assembly": "sacCer3",
            "contigs": [
                "chrM",
                "chrXI",
                "chrII",
                "chrXVI",
                "chrIII",
                "chrVI",
                "chrV",
                "chrXII",
                "chrVIII",
                "chrXV",
                "chrIV",
                "chrI",
                "chrXIII",
                "chrX",
                "chrIX",
                "chrVII",
                "chrXIV",
            ],
        },
        "sacCer2": {
            "assembly": "sacCer2",
            "contigs": [
                "chrM",
                "2micron",
                "chrXI",
                "chrII",
                "chrXVI",
                "chrIII",
                "chrVI",
                "chrV",
                "chrXII",
                "chrVIII",
                "chrXV",
                "chrIV",
                "chrI",
                "chrXIII",
                "chrX",
                "chrIX",
                "chrVII",
                "chrXIV",
            ],
        },
    }
    threads = 2

    # Uniq assembly not folder provided
    with TemporaryDirectory() as tmpdir:

        assemblies = ["sacCer3"]

        genomes_folder = None
        provider = "UCSC"
        contig_regex = None
        try:
            genome = databases_download_genomes(
                assemblies=assemblies,
                genomes_folder=genomes_folder,
                provider=provider,
                contig_regex=contig_regex,
                threads=threads,
            )
            for assembly in assemblies:
                genome = genomepy.Genome(assembly, genomes_dir=DEFAULT_GENOME_FOLDER)
                assert os.path.exists(genome.genome_file)
                assert (
                    list(genome.keys()).sort()
                    == assemblies_config.get(assembly).get("contigs", []).sort()
                )
        except:
            assert False

    # Uniq assembly
    with TemporaryDirectory() as tmpdir:

        assemblies = ["sacCer3"]

        genomes_folder = tmpdir
        provider = "UCSC"
        contig_regex = None
        try:
            genome = databases_download_genomes(
                assemblies=assemblies,
                genomes_folder=genomes_folder,
                provider=provider,
                contig_regex=contig_regex,
                threads=threads,
            )
            for assembly in assemblies:
                genome = genomepy.Genome(assembly, genomes_dir=genomes_folder)
                assert os.path.exists(genome.genome_file)
                assert (
                    list(genome.keys()).sort()
                    == assemblies_config.get(assembly).get("contigs", []).sort()
                )
        except:
            assert False

    # Multiple assemblies
    with TemporaryDirectory() as tmpdir:

        assemblies = ["sacCer2", "sacCer3"]

        genomes_folder = tmpdir
        provider = "UCSC"
        contig_regex = None
        try:
            genome = databases_download_genomes(
                assemblies=assemblies,
                genomes_folder=genomes_folder,
                provider=provider,
                contig_regex=contig_regex,
                threads=threads,
            )
            for assembly in assemblies:
                genome = genomepy.Genome(assembly, genomes_dir=genomes_folder)
                assert os.path.exists(genome.genome_file)
                assert (
                    list(genome.keys()).sort()
                    == assemblies_config.get(assembly).get("contigs", []).sort()
                )
        except:
            assert False

    # Filtered assembl
    with TemporaryDirectory() as tmpdir:

        assemblies = ["sacCer3"]

        genomes_folder = tmpdir
        provider = "UCSC"
        contig_regex = "^>chrX.*$"
        try:
            genome = databases_download_genomes(
                assemblies=assemblies,
                genomes_folder=genomes_folder,
                provider=provider,
                contig_regex=contig_regex,
                threads=threads,
            )
            for assembly in assemblies:
                genome = genomepy.Genome(assembly, genomes_dir=genomes_folder)
                assert os.path.exists(genome.genome_file)
                assert (
                    list(genome.keys()).sort()
                    == [
                        "chrXI",
                        "chrXVI",
                        "chrXII",
                        "chrXV",
                        "chrXIII",
                        "chrX",
                        "chrXIV",
                    ].sort()
                )
        except:
            assert False
