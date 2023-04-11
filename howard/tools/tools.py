#!/usr/bin/env python

import io
import multiprocessing
import os
import re
import subprocess
from tempfile import NamedTemporaryFile
import tempfile
import duckdb
import json
import argparse
import Bio.bgzf as bgzf
import pandas as pd
import vcf
import logging as log
import sys

# Import Commons
from howard.commons import *

# Import tools
from howard.tools.process import *
from howard.tools.annotation import *
from howard.tools.calculation import *
from howard.tools.prioritization import *
from howard.tools.query import *
from howard.tools.databases import *



# Arguments dict
arguments = {

        # Process & other
        "input": {
            "metavar": "FILE",
            "help": "Input file path\nFormat: BCF, VCF, TSV, CSV, PSV, Parquet or duckDB\nFiles can be compressesd (e.g. vcf.gz, tsv.gz)",
            "required": False
        },
        "output": {
            "metavar": "FILE",
            "help": "Output file path\nFormat: BCF, VCF, TSV, CSV, PSV, Parquet or duckDB\nFiles can be compressesd (e.g. vcf.gz, tsv.gz)",
            "required": False
        },
        "param": {
            "metavar": "JSON",
            "help": "Parameters file\nFormat: JSON\nDefault: {}",
            "default": "{}"
        },
        "query": {
            "metavar": "QUERY",
            "help": "Query in SQL format\nFormat: SQL\nExample: 'SELECT * FROM variants LIMIT 5'",
            "default": None
        },
        "output_query": {
            "metavar": "FILE",
            "help": "Output Query file\nformat: VCF, TSV, Parquet...",
            "default": None
        },
        "annotations": {
            "metavar": "ANNOTATION",
            "help": """Annotation with databases files, or with tools\n"""
                    """Format: list of files in Parquet, VCF, BED, or keywords\n"""
                    """For a Parquet/VCF/BED file, use file path (e.g. '/path/to/file.parquet')\n"""
                    """For snpeff annotation, use keyword 'snpeff'\n"""
                    """For Annovar annotation, use keyword 'annovar' with annovar code (e.g. 'annovar:refGene', 'annovar:cosmic70')\n"""
                    ,
            "nargs": "+",
            "default": None
        },
        "calculations": {
            "metavar": "OPERATION",
            "help": """Calculations on genetic variants information and genotype information\n"""
                    """List of available calculations (see doc for more information):\n"""
                    """ VARTYPE """
                    """ snpeff_hgvs """
                    """ FINDBYPIPELINE """
                    """ GENOTYPECONCORDANCE """
                    """ BARCODE """
                    """ TRIO """
                    """ VAF """
                    """ VAF_STATS """
                    """ DP_STATS """
                    ,
            "nargs": "+",
            "default": None
        },
        "prioritizations": {
            "metavar": "JSON",
            "help": "Prioritization file in JSON format (defines profiles, see doc)",
            "default": None
        },
        "profiles": {
            "metavar": "PROFILES",
            "help": """Prioritization profiles to use (based on file in JSON)\n"""
                    """default: all profiles available""",
            "nargs": "+",
            "default": None
        },
        "default_profile": {
            "metavar": "PROFILE",
            "help": """Prioritization profile by default (see doc)\n"""
                    """default: First profile in JSON file""",
            "default": None
        },
        "pzfields": {
            "metavar": "PZFIELD",
            "help": """Prioritization fields to provide (see doc)\n"""
                    """available: PZScore, PZFlag, PZTags, PZComment, PZInfos\n"""
                    """default: PZScore,PZFlag""",
            "nargs": "+",
            "default": ["PZScore","PZFlag"]
        },
        "prioritization_score_mode": {
            "metavar": "MODE",
            "help": """Prioritization Score mode (see doc)\n"""
                    """available: HOWARD (increment score), VaRank (max score)\n"""
                    """default: HOWARD""",
            "default": '"HOWARD'
        },
        "explode_infos": {
            "help": """Explode VCF INFO/Tag into 'variants' table columns\n"""
                    """default: False""",
            "action": "store_true",
            "default": False
        },
        "overview": {
            "help": "Overview after loading data",
            "action": "store_true"
        },
        "overview_header": {
            "help": "Overview after loading data",
            "action": "store_true"
        },
        "overview_footer": {
            "help": "Overview before data processing",
            "action": "store_true"
        },
        "stats": {
            "help": "Statistics after loading data",
            "action": "store_true"
        },
        "stats_header": {
            "help": "Statistics after loading data",
            "action": "store_true"
        },
        "stats_footer": {
            "help": "Statistics before data processing",
            "action": "store_true"
        },
        "assembly": {
            "metavar": "ASSEMBLY",
            "help": """Assembly to download\n"""
                    """Default: 'hg19'""",
            "nargs": "+",
            "required": False,
            "default": ["hg19"]
        },

        # Databases
        "download-annovar": {
            "metavar": "FOLDER",
            "help": "Download Annovar databases within Annovar folder",
            "required": False
        },
        "download-annovar-files": {
            "metavar": "FILE",
            "help": """Download Annovar databases for a list of Annovar file code (see Annovar Doc)\n"""
                    """Default: All available files\n"""
                    """Example: refGene,gnomad211_exome,cosmic70,clinvar_202*,nci60\n"""
                    """Note: refGene will be at leaset downloaded\n"""
                    """Note2: Only file that not exists or with a different size will be downloaded""",
            "nargs": '+',
            "required": False
        },
        "download-annovar-url": {
            "metavar": "URL",
            "help": """Download Annovar databases URL (see Annovar Doc)\n"""
                    """Default: 'http://www.openbioinformatics.org/annovar/download/'""",
            "required": False,
            "default": "http://www.openbioinformatics.org/annovar/download/"
        },
        "download-snpeff": {
            "metavar": "FOLDER",
            "help": """Download snpEff databases within snpEff folder""",
            "required": False
        },

        # Shared
        "config": {
            "metavar": "JSON",
            "help": "Configuration file\nFormat: JSON\nDefault: {}",
            "required": False,
            "default": "{}"
        },
        "threads": {
            "metavar": "INTEGER",
            "help": "Number of threads (replace config)",
            "required": False,
            "default": None
        },
        "verbosity": {
            "metavar": "LEVEL",
            "help": "Verbosity level (CRITICAL, ERROR, WARNING, INFO, DEBUG or NOTSET)\nDefault: INFO",
            "required": False,
            "default": "warning"
        },
        "quiet": {
            "help": argparse.SUPPRESS,
            "action": "store_true"
        },
        "verbose": {
            "help": argparse.SUPPRESS,
            "action": "store_true"
        },
        "debug": {
            "help": argparse.SUPPRESS,
            "action": "store_true"
        },

    }

   # analysis_parser.add_argument(
    #     "--overview", "--overview_header",
    #     help="Overview after loading data",
    #     action="store_true"
    # )
    # analysis_parser.add_argument(
    #     "--overview_footer",
    #     help="Overview before data processing",
    #     action="store_true"
    # )
    # analysis_parser.add_argument(
    #     "--stats", "--stats_header",
    #     help="Statistics after loading data",
    #     action="store_true"
    # )
    # analysis_parser.add_argument(
    #     "--stats_footer",
    #     help="Statistics before data processing",
    #     action="store_true"
    # )


# Shared arguments
shared_arguments = ["config", "threads", "verbosity", "quiet", "verbose", "debug"]

# Command dict
commands_arguments = {
    "process": {
        "function" : "process",
        "description":  """howard process tool manage genetic variations to:\n"""
                        """- annotates genetic variants with multiple annotation databases/files and tools\n"""
                        """- calculates and normalizes annotations\n"""
                        """- prioritizes variants with profiles (list of citeria) to calculate scores and flags\n"""
                        """- translates into various formats\n"""
                        """- query genetic variants and annotations\n"""
                        """- generates variants statistics""",
        "help":         """Full genetic variations process: annotation, calculation, prioritization, format, query, filter...""",
        "epilog":       """Usage examples:\n"""
                        """   howard process --input=tests/data/example.vcf.gz --output=/tmp/example.annotated.vcf.gz --param=config/param.json\n"""
                        """   howard process --input=tests/data/example.vcf.gz --annotations=snpeff --calculations=snpeff_hgvs,NOMEN --prioritizations=config/prioritization_profiles.json --query='SELECT "INFO/NOMEN" FROM variants'""", 

        "groups": {
            "main": {
                "input": True,
                "output": True,
                "param": False
            },
            "Quick Processes": {
                "annotations": False,
                "calculations": False,
                "prioritizations": False,
                "query": False
            }
        }
    },
    "annotation": {
        "function" : "annotation",
        "description":  """Annotation is mainly based on a build-in Parquet annotation method, and tools such as BCFTOOLS, ANNOVAR and snpEff. It using available databases (see ANNOVAR and snpEff) and homemade databases. Format of databases are: parquet, duckdb, vcf, bed, ANNOVAR and snpEff (ANNOVAR and snpEff databases are automatically downloaded if needed). """,
        "help":         """Annotation of genetic variations using databases/files and tools.""",
        "epilog":       """Usage examples:\n"""
                        """   howard annotation --input=tests/data/example.vcf.gz --output=/tmp/example.annotated.vcf.gz --annotations=tests/data/annotations/dbnsfp42a.parquet,tests/data/annotations/gnomad211_exome.vcf.gz,tests/data/annotations/refGene.bed.gz,snpeff,annovar:refGene""", 
        "groups": {
            "main": {
                "input": True,
                "output": True,
                "annotations": True
            }
        }
    },
    "calculation": {
        "function" : "calculation",
        "description":  """Calculation processes variants information to calculate new information, such as: harmonizes allele frequency (VAF), extracts Nomen (transcript, cNomen, pNomen...) from HGVS fields with an optional list of personalized transcripts, generates VaRank format barcode.""",
        "help":         """Calculation operations on genetic variations and genotype information.\n""",
        "epilog":       """Usage examples:\n"""
                        """   howard calculation --input=tests/data/example.ann.vcf.gz --output=/tmp/example.calculated.vcf.gz --calculations=snpeff_hgvs,NOMEN,VARTYPE""", 
        "groups": {
            "main": {
                "input": True,
                "output": True,
                "calculations": True
            }
        }
    },
    "prioritization": {
        "function" : "prioritization",
        "description":  """Prioritization algorithm uses profiles to flag variants (as passed or filtered), calculate a prioritization score, and automatically generate a comment for each variants (example: 'polymorphism identified in dbSNP. associated to Lung Cancer. Found in ClinVar database').Prioritization profiles are defined in a configuration file. A profile is defined as a list of annotation/value, using wildcards and comparison options (contains, lower than, greater than, equal...). Annotations fields may be quality values (usually from callers, such as 'GQ', 'DP') or other annotations fields provided by annotations tools, such as HOWARD itself (example: COSMIC, Clinvar, 1000genomes, PolyPhen, SIFT). Multiple profiles can be used simultaneously, which is useful to define multiple validation/prioritization levels (example: 'standard', 'stringent', 'rare variants', 'low allele frequency').\n"""
                        """and profiles score/flag profile""",
        "help": "Prioritization of genetic variations based on annotations criteria (profiles)",
        "epilog": """Usage examples:\n"""
                        """   howard prioritization --input=tests/data/example.ann.vcf.gz --output=/tmp/example.prioritized.vcf.gz --prioritizations=config/prioritization_profiles.json --profiles=default,GERMLINE""", 
        "groups": {
            "main": {
                "input": True,
                "output": True,
                "prioritizations": True,
                "profiles": False,
                "default_profile": False,
                "pzfields": False,
                "prioritization_score_mode": False
            }
        }
    },
    "query": {
        "function" : "query",
        "description":  """Query genetic variations in SQL format. """
                        """Data can be loaded into 'variants' table from various format (e.g. VCF, TSV, Parquet...). """
                        """SQL query can also use external data within th request, such as a parquet file """,
        "help": "Query genetic variations in SQL format",
        "epilog": """Usage examples:\n"""
                        """   howard query --input=tests/data/example.vcf.gz --query="SELECT * FROM variants WHERE REF = 'A' AND POS < 100000" n"""
                        """   howard query --input=tests/data/example.vcf.gz --explode_infos --query='SELECT "#CHROM", POS, REF, ALT, "INFO/DP", "INFO/CLNSIG", sample2, sample3 FROM variants WHERE "INFO/DP" >=50 OR "INFO/CLNSIG" NOT NULL ORDER BY "INFO/DP" DESC' \n"""
                        """   howard query --query="SELECT * FROM 'tests/data/annotations/dbnsfp42a.parquet' WHERE \\"INFO/Interpro_domain\\" NOT NULL ORDER BY \\"INFO/SiPhy_29way_logOdds_rankscore\\" DESC" \n"""
                        """   howard query --query="SELECT \\"#CHROM\\", POS, '' AS ID, REF, ALT, '' AS QUAL, '' AS FILTER, STRING_AGG(INFO, ';') AS INFO FROM 'tests/data/annotations/*.parquet' GROUP BY \\"#CHROM\\", POS, REF, ALT" --output=/tmp/full_annotation.parquet """
                        , 
        "groups": {
            "main": {
                "input": False,
                "output": False,
                "query": True,
                "explode_infos": False
            }
        }
    },
    "databases": {
        "function" : "databases",
        "description": """Download databases and needed files for howard and associated tools""",
        "help": """Download databases and needed files for howard and associated tools""",
        "epilog": "", 
        "groups": {
            "main": {
                "assembly": False,
            },
            "snpEff": {
                "download-snpeff": False
            },
            "Annovar": {
                "download-annovar": False,
                "download-annovar-files": False,
                "download-annovar-url": False
            }
        }
    }
}


# get argument
def get_argument(arguments:dict = {}, arg:str = "", required:bool = False) -> dict:
    if arg in arguments:
        arg_infos = arguments.get(arg,{})
        if required != None:
            arg_infos["required"] = required
        return arg_infos
    else:
        return {}

