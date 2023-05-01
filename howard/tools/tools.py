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
from howard.tools.stats import *
from howard.tools.convert import *
from howard.tools.databases import *
from howard.tools.from_annovar import *



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
        "export_infos": {
            "help": """Export VCF INFO/Tag into columns file\n"""
                    """Available only for non VCF format file (e.g. TSV, Parquet...) \n"""
                    """default: False""",
            "action": "store_true",
            "default": False
        },
        "export_infos_prefix": {
            "help": """VCF INFO/Tag prefix for exported columns\n"""
                    """default: 'INFO/'""",
            "metavar": "PREFIX",
            "default": "INFO/"
        },
        "multi_variant": {
            "help": """Variant with multiple annotation lines\n"""
                    """Values: 'auto' (auto-detection), 'enable', 'disable'\n"""
                    """default: 'auto'""",
            "metavar": "BOOL",
            "default": "auto"
        },
        "reduce_memory": {
            "help": """Reduce memory option\n"""
                    """Values: 'auto' (auto-detection), 'enable', 'disable'\n"""
                    """default: 'auto'""",
            "metavar": "BOOL",
            "default": "auto"
        },
        "show_calculations": {
            "help": """Show available calculation operations\n""",
            "action": "store_true"
        },
        "hgvs_field": {
            "help": """HGVS INFO/tag containing a list o HGVS annotations\n"""
                    """default: 'hgvs'""",
            "metavar": "TAG",
            "default": "hgvs"
        },
        "transcripts": {
            "help": """Transcripts file in TSV format\n"""
                    """Format: Transcript in first column, optional Gene in second column \n"""
                    """default: None""",
            "metavar": "FILE",
            "default": None
        },
        "trio_pedigree": {
            "help": """Pedigree Trio for trio inheritance calculation\n"""
                    """Format: JSON file or dict (e.g. 'trio.ped.json', '{"father":"sample1", "mother":"sample2", "child":"sample3"}') \n"""
                    """default: None""",
            "metavar": "JSON",
            "default": None
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
            "help": """Genome Assembly\n"""
                    """Default: 'hg19'""",
            "nargs": "+",
            "required": False,
            "default": ["hg19"]
        },
        "genome": {
            "metavar": "GENOME",
            "help": """Genome file in fasta format\n"""
                    """Default: 'hg19.fa'""",
            "required": False,
            "default": "hg19.fa"
        },
        "to_parquet": {
            "metavar": "FILE",
            "help": """Parquet file conversion\n""",
            "required": False,
            "default": None
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

        # From Annovar
        "annovar-code": {
            "metavar": "CODE",
            "help": """Annovar code, or database name. Usefull to name databases columns""",
            "required": False,
            "default": None
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
            "help": "Number of threads to use\nFormat: INTEGER\nDefault: 1",
            "required": False,
            "default": 1
        },
        "memory": {
            "metavar": "FLOAT[kMG]",
            "help": "Memory to use\nFormat: FLOAT[kMG]\nDefault: 8G",
            "required": False,
            "default": "8G"
        },
        "verbosity": {
            "metavar": "LEVEL",
            "help": "Verbosity level (CRITICAL, ERROR, WARNING, INFO, DEBUG or NOTSET)\nDefault: INFO",
            "required": False,
            "default": "info"
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
shared_arguments = ["config", "threads", "memory", "verbosity", "quiet", "verbose", "debug"]

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
                        """   howard process --input=tests/data/example.vcf.gz --output=/tmp/example.annotated.vcf.gz --param=config/param.json \n"""
                        """   howard process --input=tests/data/example.vcf.gz --annotations=snpeff --calculations=snpeff_hgvs,NOMEN --prioritizations=config/prioritization_profiles.json --query='SELECT "INFO/NOMEN" FROM variants' \n""", 

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
        "description":  """Annotation is mainly based on a build-in Parquet annotation method, and tools such as BCFTOOLS, Annovar and snpEff. It uses available databases (see Annovar and snpEff) and homemade databases. Format of databases are: parquet, duckdb, vcf, bed, Annovar and snpEff (Annovar and snpEff databases are automatically downloaded, see howard databases tool). """,
        "help":         """Annotation of genetic variations using databases/files and tools.""",
        "epilog":       """Usage examples:\n"""
                        """   howard annotation --input=tests/data/example.vcf.gz --output=/tmp/example.howard.vcf.gz --annotations=tests/data/annotations/avsnp150.parquet,tests/data/annotations/dbnsfp42a.parquet,tests/data/annotations/gnomad211_genome.parquet \n"""
                        """   howard annotation --input=tests/data/example.vcf.gz --output=/tmp/example.howard.tsv --annotations=annovar:refGene,snpeff,tests/data/annotations/clinvar_20210123.parquet \n""", 
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
        "description":  """Calculation processes variants information to generate new information, such as: identify variation type (VarType), harmonizes allele frequency (VAF) and calculate sttistics (VAF_stats), extracts Nomen (transcript, cNomen, pNomen...) from an HGVS field (e.g. snpEff, Annovar) with an optional list of personalized transcripts, generates VaRank format barcode, identify trio inheritance.""",
        "help":         """Calculation operations on genetic variations and genotype information.\n""",
        "epilog":       """Usage examples:\n"""
                        """   howard calculation --input=tests/data/example.full.vcf --output=/tmp/example.calculation.tsv --calculations=vartype \n"""
                        """   howard calculation --input=tests/data/example.ann.vcf.gz --output=/tmp/example.calculated.tsv --calculations=snpeff_hgvs,NOMEN --hgvs_field=snpeff_hgvs --transcripts=tests/data/transcripts.tsv \n"""
                        ,
        "groups": {
            "main": {
                "input": True,
                "output": True,
                "calculations": True,
                "show_calculations": False
            },
            "NOMEN calculation": {
                "hgvs_field": False,
                "transcripts": False
            },
            "TRIO calculation": {
                "trio_pedigree": False
            }
        }
    },
    "prioritization": {
        "function" : "prioritization",
        "description":  """Prioritization algorithm uses profiles to flag variants (as passed or filtered), calculate a prioritization score, and automatically generate a comment for each variants (example: 'polymorphism identified in dbSNP. associated to Lung Cancer. Found in ClinVar database'). Prioritization profiles are defined in a configuration file in JSON format. A profile is defined as a list of annotation/value, using wildcards and comparison options (contains, lower than, greater than, equal...). Annotations fields may be quality values (usually from callers, such as 'DP') or other annotations fields provided by annotations tools, such as HOWARD itself (example: COSMIC, Clinvar, 1000genomes, PolyPhen, SIFT). Multiple profiles can be used simultaneously, which is useful to define multiple validation/prioritization levels (example: 'standard', 'stringent', 'rare variants', 'low allele frequency').\n""",
        "help": "Prioritization of genetic variations based on annotations criteria (profiles).",
        "epilog": """Usage examples:\n"""
                        """   howard prioritization --input=tests/data/example.vcf.gz --output=/tmp/example.prioritized.vcf.gz --prioritizations=config/prioritization_profiles.json --profiles=default,GERMLINE \n""", 
        "groups": {
            "main": {
                "input": True,
                "output": True,
                "prioritizations": True
            },
            "Prioritization": {
                "profiles": False,
                "default_profile": False,
                "pzfields": False,
                "prioritization_score_mode": False
            }
        }
    },
    "query": {
        "function" : "query",
        "description":  """Query genetic variations in SQL format. Data can be loaded into 'variants' table from various formats (e.g. VCF, TSV, Parquet...). Using --explode_infos allow query on INFO/tag annotations. SQL query can also use external data within the request, such as a Parquet file(s).  """,
        "help": "Query genetic variations in SQL format.",
        "epilog": """Usage examples:\n"""
                        """   howard query --input=tests/data/example.vcf.gz --query="SELECT * FROM variants WHERE REF = 'A' AND POS < 100000" \n"""
                        """   howard query --input=tests/data/example.vcf.gz --explode_infos --query='SELECT "#CHROM", POS, REF, ALT, "INFO/DP", "INFO/CLNSIG", sample2, sample3 FROM variants WHERE "INFO/DP" >= 50 OR "INFO/CLNSIG" NOT NULL ORDER BY "INFO/DP" DESC' \n"""
                        """   howard query --query="SELECT * FROM 'tests/data/annotations/dbnsfp42a.parquet' WHERE \\"INFO/Interpro_domain\\" NOT NULL ORDER BY \\"INFO/SiPhy_29way_logOdds_rankscore\\" DESC" \n"""
                        """   howard query --query="SELECT \\"#CHROM\\" AS \\"#CHROM\\", POS AS POS, '' AS ID, REF AS REF, ALT AS ALT, '' AS QUAL, '' AS FILTER, STRING_AGG(INFO, ';') AS INFO FROM 'tests/data/annotations/*.parquet' GROUP BY \\"#CHROM\\", POS, REF, ALT" --output=/tmp/full_annotation.tsv \n"""
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
    "stats": {
        "function" : "stats",
        "description":  """Statistics on genetic variations, such as: number of variants, number of samples, statistics by chromosome, genotypes by samples...""",
        "help": "Statistics on genetic variations.",
        "epilog": """Usage examples:\n"""
                        """   howard stats --input=tests/data/example.vcf.gz """, 
        "groups": {
            "main": {
                "input": True
            }
        }
    },
    "convert": {
        "function" : "convert",
        "description":  """Convert genetic variations file to another format. Multiple format are available, such as usual and official VCF and BCF format, but also other formats such as TSV, CSV, PSV and Parquet/duckDB. These formats need a header '.hdr' file to take advantage of the power of howard (especially through INFO/tag definition), and using howard convert tool automatically generate header file fo futher use. """,
        "help": "Convert genetic variations file to another format.",
        "epilog": """Usage examples:\n"""
                        """   howard convert --input=tests/data/example.vcf.gz --output=/tmp/example.tsv --export_infos"""
                        """   howard convert --input=tests/data/example.vcf.gz --output=/tmp/example.parquet""", 
        "groups": {
            "main": {
                "input": True,
                "output": True,
                "export_infos": False,
                "export_infos_prefix": False,
            }
        }
    },
    "databases": {
        "function" : "databases",
        "description": """Download databases and needed files for howard and associated tools""",
        "help": """Download databases and needed files for howard and associated tools""",
        "epilog": """Usage examples:\n"""
                    """   howard databases --assembly=hg19 --download-annovar=/databases/annovar/current --download-annovar-files=refGene,gnomad_exome,dbnsfp42a,cosmic70,clinvar_202*,nci60 --download-snpeff=/databases/annovar/current """, 
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
    },
    "from_annovar": {
        "function" : "from_annovar",
        "description": """(beta) Formatting Annovar database file to other format (VCF and Parquet). Exported Parquet file includes INFO/tags columns as VCF INFO columns had been exploded""",
        "help": """(beta) Formatting Annovar database file to other format (VCF and Parquet)""",
        "epilog": """Usage examples:\n"""
                    """   howard from_annovar --input=/databases/annovar/current/hg19_nci60.txt --output=/databases/annotations/current/hg19/nci60.vcf.gz --to_parquet=/databases/annotations/current/hg19/nci60.parquet --annovar-code=nci60 --genome=/databases/genomes/current/hg19.fa --config=/tool/config/config.json --threads=8 """, 
        "groups": {
            "main": {
                "input": True,
                "output": True,
                "genome": True,
            },
            "Annovar": {
                "annovar-code": False,
            },
            "Parquet": {
                "to_parquet": False,
                #"export_infos": False,
                #"export_infos_prefix": False,
            },
            "Modes": {
                "reduce_memory": False,
                "multi_variant": False
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

