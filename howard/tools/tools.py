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
from howard.tools.hgvs import *
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
            "metavar": "input",
            "help": "Input file path\nFormat: BCF, VCF, TSV, CSV, PSV, Parquet or duckDB\nFiles can be compressesd (e.g. vcf.gz, tsv.gz)",
            "required": False,
            "type": argparse.FileType('r'),
            "gooey": {
                "widget": "FileChooser"
            }
        },
        "output": {
            "metavar": "output",
            "help": "Output file path\nFormat: BCF, VCF, TSV, CSV, PSV, Parquet or duckDB\nFiles can be compressesd (e.g. vcf.gz, tsv.gz)",
            "required": False,
            "type": argparse.FileType('w'),
            "gooey": {
                "widget": "FileSaver"
            }
        },
        "param": {
            "metavar": "param",
            "help": "Parameters file or JSON\nFormat: JSON\nDefault: {}",
            "default": "{}",
            "gooey": {
                "widget": "FileChooser",
                "options": {
                    'initial_value': ''  
                }
            }
        },
        "query": {
            "metavar": "query",
            "help": "Query in SQL format\nFormat: SQL\nExample: 'SELECT * FROM variants LIMIT 5'",
            "default": None,
            "gooey": {
                "widget": "Textarea",
                "options": {
                    'initial_value': 'SELECT * FROM variants'  
                }
            }
        },
        "output_query": {
            "metavar": "output",
            "help": "Output Query file\nformat: VCF, TSV, Parquet...",
            "default": None,
            "gooey": {
                "widget": "FileSaver"
            }
        },
        "annotations": {
            "metavar": "annotations",
            "help": """Annotation with databases files, or with tools\n"""
                    """Format: list of files in Parquet, VCF, BED, or keywords\n"""
                    """For a Parquet/VCF/BED file, use file path (e.g. '/path/to/file.parquet')\n"""
                    """For add all availalbe databases, use 'ALL' keyword:\n"""
                    """   - Use 'ALL:<types>:<releases>'\n"""
                    """   - e.g. 'ALL', 'ALL:parquet:current', 'ALL:parquet,vcf:devel'\n"""
                    """For snpeff annotation, use keyword 'snpeff'\n"""
                    """For Annovar annotation, use keyword 'annovar' with annovar code (e.g. 'annovar:refGene', 'annovar:cosmic70')\n"""
                    ,
            "default": None,
            "gooey": {
                "widget": "MultiFileChooser"
            }
        },
        "calculations": {
            "metavar": "operations",
            "help": """Calculations on genetic variants information and genotype information\n"""
                    """Example: 'VARTYPE,barcode'\n"""
                    """List of available calculations (unsensitive case, see doc for more information):\n"""
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
            "default": None
        },
        "prioritizations": {
            "metavar": "prioritisations",
            "help": "Prioritization file in JSON format (defines profiles, see doc).",
            "default": None,
            "gooey": {
                "widget": "FileChooser"
            }
        },
        "profiles": {
            "metavar": "profiles",
            "help": """Prioritization profiles to use (based on file in JSON).\n"""
                    """default: all profiles available""",
            "default": None
        },
        "default_profile": {
            "metavar": "default profile",
            "help": """Prioritization profile by default (see doc)\n"""
                    """default: First profile in JSON file""",
            "default": None
        },
        "pzfields": {
            "metavar": "pzfields",
            "help": """Prioritization fields to provide (see doc).\n"""
                    """available: PZScore, PZFlag, PZTags, PZComment, PZInfos\n"""
                    """default: PZScore,PZFlag""",
            "default": "PZScore,PZFlag"
        },
        "prioritization_score_mode": {
            "metavar": "prioritization score mode",
            "help": """Prioritization Score mode (see doc).\n"""
                    """available: HOWARD (increment score), VaRank (max score)\n"""
                    """default: HOWARD""",
            "default": 'HOWARD',
            "choices": ["HOWARD", "VaRank"],
            "gooey": {
                "widget": "Dropdown",
                "options": {}
            }
        },

        # Query print options
        "query_limit": {
            "metavar": "query limit",
            "help": """Limit of number of row for query (only for print result, not output).\n"""
                    """default: 10""",
            "default": 10
        },
        "query_print_mode": {
            "metavar": "print mode",
            "help": """Print mode of query result (only for print result, not output).\n"""
                    """Either None (native), 'markdown' or 'tabulate'.\n"""
                    """default: None""",
            "default": None
        },

        # Explode infos
        "explode_infos": {
            "help": """Explode VCF INFO/Tag into 'variants' table columns.\n"""
                    """default: False""",
            "action": "store_true",
            "default": False
        },
        "explode_infos_prefix": {
            "metavar": "explode infos prefix",
            "help": """Explode VCF INFO/Tag with a specific prefix.\n"""
                    """default: ''""",
            "default": ""
        },
        "explode_infos_fields": {
            "metavar": "explode infos list",
            "help": """Explode VCF INFO/Tag specific fields/tags.\n"""
                    """Keyword '*' specify all available fields, except those already specified.\n"""
                    """Pattern (regex) can be used: '.*_score' for fields named with '_score' at the end.\n"""
                    """Examples:\n"""
                    """   - 'HGVS,SIFT,Clinvar' (list of fields)\n"""
                    """   - 'HGVS,*,Clinvar' (list of fields with all other fields at the end)\n"""
                    """   - 'HGVS,.*_score,Clinvar' (list of 2 fields with all scores in the middle)\n"""
                    """   - 'HGVS,.*_score,*' (1 field, scores, all other fields)\n"""
                    """   - 'HGVS,*,.*_score' (1 field and all other fields,\n"""
                    """                        scores included in other fields)\n"""
                    """default: '*'""",
            "default": "*"
        },

        # Include header
        "include_header": {
            "help": """Include header (in VCF format) in output file.\n"""
                    """Only for compatible formats (tab-delimiter format as TSV or BED).\n"""
                    """default: False""",
            "action": "store_true",
            "default": False
        },

        # Sort By
        "order_by": {
            "metavar": "order by",
            "help": """List of columns to sort the result-set in ascending or descending order.\n"""
                    """Use SQL format, and keywords ASC (ascending) and DESC (descending).\n"""
                    """If a column is not available, order will not be considered.\n"""
                    """Order is enable only for compatible format (e.g. TSV, CSV, JSON).\n"""
                    """Examples:\n"""
                    """   - 'ACMG_score DESC'\n"""
                    """   - 'PZFlag DESC, PZScore DESC'\n"""
                    """default: ''""",
            "default": ""
        },

        # Parquet partition
        "parquet_partitions": {
            "metavar": "parquet partitions",
            "help": """Parquet partitioning using huve (only for Parquet export format).\n"""
                    """This option is is faster parallel writing, but memory consuming.\n"""
                    """Use 'None' (string) for NO partition but split parquet files into a folder\n"""
                    """examples: '#CHROM', '#CHROM,REF', 'None'\n"""
                    """default: None""",
            "default": None
        },
        
        # From annovar
        "multi_variant": {
            "metavar": "multi variant",
            "help": """Variant with multiple annotation lines\n"""
                    """Values: 'auto' (auto-detection), 'enable', 'disable'\n"""
                    """default: 'auto'""",
            "default": "auto"
        },
        "reduce_memory": {
            "metavar": "reduce memory",
            "help": """Reduce memory option\n"""
                    """Values: 'auto' (auto-detection), 'enable', 'disable'\n"""
                    """default: 'auto'""",
            "default": "auto"
        },

        # Calculation
        "calculation_config": {
            "metavar": "calculation config",
            "help": """Calculation config file\n"""
                    """Format: JSON""",
            "default": None,
            "gooey": {
                "widget": "FileChooser"
            }
        },
        "show_calculations": {
            "help": """Show available calculation operations""",
            "action": "store_true"
        },
        "hgvs_field": {
            "metavar": "HGVS field",
            "help": """HGVS INFO/tag containing a list o HGVS annotations\n"""
                    """default: 'hgvs'""",
            "default": "hgvs"
        },
        "transcripts": {
            "metavar": "transcripts",
            "help": """Transcripts file in TSV format\n"""
                    """Format: Transcript in first column, optional Gene in second column \n"""
                    """default: None""",
            "default": None,
            "gooey": {
                "widget": "FileChooser"
            }
        },
        "trio_pedigree": {
            "metavar": "trio pedigree",
            "help": """Pedigree Trio for trio inheritance calculation\n"""
                    """Format: JSON file or dict (e.g. 'trio.ped.json', '{"father":"sample1", "mother":"sample2", "child":"sample3"}') \n"""
                    """default: None""",
            "default": None,
            "gooey": {
                "widget": "FileChooser"
            }
        },

        # Other
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

        # Stats
        "stats_md": {
            "metavar": "stats markdown",
            "help": """Stats Output file in MarkDown format\n""",
            "required": False,
            "gooey": {
                "widget": "FileSaver"
            }
        },
        "stats_json": {
            "metavar": "stats json",
            "help": """Stats Output file in JSON format\n""",
            "required": False,
            "gooey": {
                "widget": "FileSaver"
            }
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
            "metavar": "assembly",
            "help": """Genome Assembly\n"""
                    """Default: 'hg19'""",
            "required": False,
            "default": "hg19"
        },
        "genome": {
            "metavar": "genome",
            "help": """Genome file in fasta format\n"""
                    """Default: 'hg19.fa'""",
            "required": False,
            "default": "hg19.fa"
        },
        "to_parquet": {
            "metavar": "to parquet",
            "help": """Parquet file conversion\n""",
            "required": False,
            "default": None,
            "gooey": {
                "widget": "FileSaver"
            }
        },

        # HGVS
        "use_gene": {
            "help": """Use Gene information to generate HGVS annotation\n"""
                    """Example: 'NM_152232(TAS1R2):c.231T>C'""",
            "action": "store_true"
        },
        "use_exon": {
            "help": """Use Exon information to generate HGVS annotation\n"""
                    """Only if 'use_gene' is not enabled\n"""
                    """Example: 'NM_152232(exon2):c.231T>C'""",
            "action": "store_true"
        },
        "use_protein": {
            "help": """Use Protein level to generate HGVS annotation\n"""
                    """Can be used with 'use_exon' or 'use_gene'\n"""
                    """Example: 'NP_689418:p.Cys77Arg'""",
            "action": "store_true"
        },
        "add_protein": {
            "help": """Add Protein level to DNA HGVS annotation\n"""
                    """Example: 'NM_152232:c.231T>C,NP_689418:p.Cys77Arg'""",
            "action": "store_true"
        },
        "full_format": {
            "help": """Generates HGVS annotation in a full format (non-standard)\n"""
                    """Full format use all information to generates an exhaustive annotation.\n"""
                    """Use specifically 'use_exon' to add exon information.\n"""
                    """Example: 'TAS1R2:NM_152232:NP_689418:c.231T>C:p.Cys77Arg'\n"""
                    """         'TAS1R2:NM_152232:NP_689418:exon2:c.231T>C:p.Cys77Arg'""",
            "action": "store_true"
        },
        "use_version": {
            "help": """Generates HGVS annotation with transcript version\n"""
                    """Example: without version 'NM_152232:c.231T>C'\n"""
                    """         with version 'NM_152232.1:c.231T>C'""",
            "action": "store_true"
        },
        "codon_type": {
            "metavar": "Codon type",
            "help": """Amino Acide Codon format type to use to generate HGVS annotation\n"""
                    """Available (default '3'):\n"""
                    """   '1': codon in 1 caracter (e.g. 'C', 'R')\n"""
                    """   '3': codon in 3 caracter (e.g. 'Cys', 'Arg')\n"""
                    """   'FULL': codon in full name (e.g. 'Cysteine', 'Arginine')\n""",
            "required": False,
            "default": "3",
            "choices": ["1", "3", "FULL"],
            "gooey": {
                "widget": "Dropdown",
                "options": {}
            }
        },
        "refgene": {
            "metavar": "refGene",
            "help": """refGene annotation file""",
            "required": False,
            "default": "",
            "gooey": {
                "widget": "FileChooser"
            }
        },
        "refseqlink": {
            "metavar": "refSeqLink",
            "help": """refSeqLink annotation file""",
            "required": False,
            "default": "",
            "gooey": {
                "widget": "FileChooser"
            }
        },

        # Databases

        # Genome
        "download-genomes": {
            "metavar": "genomes",
            "help": "Download Genomes within folder",
            "required": False,
            "gooey": {
                "widget": "DirChooser"
            }
        },
        "download-genomes-provider": {
            "metavar": "genomes provider",
            "help": """Download Genome from an external provider\n"""
                    """Available: GENCODE, Ensembl, UCSC, NCBI\n"""
                    """Default: UCSC\n""",
            "required": False,
            "default": "UCSC"
        },
        "download-genomes-contig-regex": {
            "metavar": "genomes contig regex",
            "help": """Regular expression to select specific chromosome \n"""
                    """Default: None\n"""
                    """Example: 'chr[0-9XYM]+$'\n""",
            "required": False,
            "default": None
        },

        # Annovar
        "download-annovar": {
            "metavar": "Annovar",
            "help": "Download Annovar databases within Annovar folder",
            "required": False,
            "gooey": {
                "widget": "DirChooser"
            }
        },
        "download-annovar-files": {
            "metavar": "Annovar code",
            "help": """Download Annovar databases for a list of Annovar file code (see Annovar Doc)\n"""
                    """Default: All available files\n"""
                    """Example: refGene,gnomad211_exome,cosmic70,clinvar_202*,nci60\n"""
                    """Note: refGene will be at leaset downloaded\n"""
                    """Note2: Only file that not exists or with a different size will be downloaded""",
            "required": False,
            "default": None
        },
        "download-annovar-url": {
            "metavar": "Annovar url",
            "help": """Download Annovar databases URL (see Annovar Doc)\n"""
                    """Default: 'http://www.openbioinformatics.org/annovar/download'""",
            "required": False,
            "default": "http://www.openbioinformatics.org/annovar/download"
        },

        # snpEff
        "download-snpeff": {
            "metavar": "snpEff",
            "help": """Download snpEff databases within snpEff folder""",
            "required": False,
            "gooey": {
                "widget": "DirChooser"
            }
        },

        # refSeq
        "download-refseq": {
            "metavar": "refSeq",
            "help": """Download refSeq databases within refSeq folder""",
            "required": False,
            "gooey": {
                "widget": "DirChooser"
            }
        },
        "download-refseq-url": {
            "metavar": "refSeq url",
            "help": """Download refSeq databases URL (see refSeq WebSite)\n"""
                    """Default: 'http://hgdownload.soe.ucsc.edu/goldenPath'""",
            "required": False,
            "default": "http://hgdownload.soe.ucsc.edu/goldenPath"
        },
        "download-refseq-prefix": {
            "metavar": "refSeq prefix",
            "help": """Check existing refSeq files in refSeq folder\n"""
                    """Default: 'ncbiRefSeq'""",
            "required": False,
            "default": "ncbiRefSeq"
        },
        "download-refseq-files": {
            "metavar": "refSeq files",
            "help": """List of refSeq files to download\n"""
                    """Default: 'ncbiRefSeq.txt,ncbiRefSeqLink.txt'""",
            "required": False,
            "default": "ncbiRefSeq.txt,ncbiRefSeqLink.txt"
        },
        "download-refseq-format-file": {
            "metavar": "refSeq format file",
            "help": """Name of refSeq file to format in BED format\n"""
                    """Exemple: 'ncbiRefSeq.txt'\n"""
                    """Default: None""",
            "required": False,
            "default": None
        },
        "download-refseq-include-utr5": {
            "help": """Formating BED refSeq file including 5'UTR""",
            "action": "store_true"
        },
        "download-refseq-include-utr3": {
            "help": """Formating BED refSeq file including 3'UTR""",
            "action": "store_true"
        },
        "download-refseq-include-chrM": {
            "help": """Formating BED refSeq file including Mitochondiral chromosome 'chrM' or 'chrMT'""",
            "action": "store_true"
        },
        "download-refseq-include-non-canonical-chr": {
            "help": """Formating BED refSeq file including non canonical chromosomes""",
            "action": "store_true"
        },
        "download-refseq-include-non-coding-transcripts": {
            "help": """Formating BED refSeq file including non coding transcripts""",
            "action": "store_true"
        },
        "download-refseq-include-transcript-version": {
            "help": """Formating BED refSeq file including transcript version""",
            "action": "store_true"
        },
        
        # dbNSFP
        "download-dbnsfp": {
            "metavar": "dbNSFP",
            "help": "Download dbNSFP databases within dbNSFP folder",
            "required": False,
            "gooey": {
                "widget": "DirChooser"
            }
        },
        "download-dbnsfp-url": {
            "metavar": "dbNSFP url",
            "help": """Download dbNSFP databases URL (see dbNSFP website)\n"""
                    """Default: 'https://dbnsfp.s3.amazonaws.com'""",
            "required": False,
            "default": "https://dbnsfp.s3.amazonaws.com"
        },
        "download-dbnsfp-release": {
            "metavar": "dnNSFP release",
            "help": """Release of dbNSFP to download (see dbNSFP website)\n"""
                    """Default: '4.4a'""",
            "required": False,
            "default": "4.4a"
        },
        "download-dbnsfp-parquet-size": {
            "metavar": "dbNSFP parquet size",
            "help": """Maximum size (Mb) of data files in Parquet folder.\n"""
                    """Parquet folder are partitioned (hive) by chromosome (sub-folder),\n"""
                    """which contain N data files.\n"""
                    """Default: 100""",
            "required": False,
            "default": 100
        },
        "download-dbnsfp-subdatabases": {
            "help": """Generate dbNSFP sub-databases\n"""
                    """dbNSFP provides multiple databases which are split onto multiple columns.\n"""
                    """This option create a Parquet folder for each sub-database (based on columns names).""",
            "action": "store_true"
        },
        "download-dbnsfp-parquet": {
            "help": """Generate a Parquet file for each Parquet folder.""",
            "action": "store_true"
        },
        "download-dbnsfp-vcf": {
            "help": """Generate a VCF file for each Parquet folder.\n"""
                    """Note: Need genome (see --download-genome)""",
            "action": "store_true"
        },
        "download-dbnsfp-no-files-all": {
            "help": """Not generate database Parquet/VCF file for the entire database ('ALL').\n"""
                    """Only sub-databases files will be generated.\n"""
                    """(see '--download-dbnsfp-subdatabases')""",
            "action": "store_true"
        },
        "download-dbnsfp-add-info": {
            "help": """Add INFO column (VCF format) in Parquet folder and file.\n"""
                    """Useful for speed up full annotation (all available columns).\n"""
                    """Increase memory and space during generation of files.""",
            "action": "store_true"
        },
        "download-dbnsfp-row-group-size": {
            "metavar": "dnNSFP row grooup size",
            "help": """minimum number of rows in a parquet row group (see duckDB doc).\n"""
                    """Lower can reduce memory usage and slightly increase space during generation,\n"""
                    """speed up highly selective queries, slow down whole file queries (e.g. aggregations)\n"""
                    """Default: 100000""",
            "required": False,
            "default": 100000
        },

        # AlphaMissense
        "download-alphamissense": {
            "metavar": "AlphaMissense",
            "help": "Download AlphaMissense databases within Annotations folder",
            "required": False,
            "gooey": {
                "widget": "DirChooser"
            }
        },
        "download-alphamissense-url": {
            "metavar": "AlphaMissense url",
            "help": """Download AlphaMissense databases URL (see AlphaMissense website)\n"""
                    """Default: 'https://storage.googleapis.com/dm_alphamissense'""",
            "required": False,
            "default": "https://storage.googleapis.com/dm_alphamissense"
        },

        # Exomiser
        "download-exomiser": {
            "metavar": "Exomiser",
            "help": """Download Exomiser databases\n"""
                    """Folder where the Exomiser databases will be downloaded and stored.\n"""
                    """If the folder does not exist, it will be created.""",
            "required": False,
            "gooey": {
                "widget": "DirChooser"
            }
        },
        "download-exomiser-application-properties": {
            "metavar": "Exomiser aplication properties",
            "help": """Exomiser Application Properties configuration file (see Exomiser website)\n"""
                    """This file contains configuration settings for the Exomiser tool.\n"""
                    """If this parameter is not provided, the function will attempt to locate\n"""
                    """the application properties file automatically based on the Exomiser.\n"""
                    """Configuration information will be used to download expected releases (if no other parameters)\n"""
                    """CADD and REMM will be downloaded only if 'path' are provided\n""",
            "required": False,
            "default": None,
            "gooey": {
                "widget": "FileChooser"
            }
        },
        "download-exomiser-url": {
            "metavar": "Exomiser url",
            "help": """URL where Exomiser database files can be downloaded from.\n"""
                    """Default: 'http://data.monarchinitiative.org/exomiser'""",
            "required": False,
            "default": "http://data.monarchinitiative.org/exomiser"
        },
        "download-exomiser-release": {
            "metavar": "Exomiser release",
            "help": """Release of Exomiser data to download.\n"""
                    """If "default", "auto", or "config", retrieve from Application Properties file.\n"""
                    """Default: None""",
            "required": False,
            "default": None
        },
        "download-exomiser-phenotype-release": {
            "metavar": "Exomiser phenoptye release",
            "help": """Release of Exomiser phenotype to download.\n"""
                    """If not provided, retrieve from Application Properties file or Exomiser data release\n"""
                    """Default: None""",
            "required": False,
            "default": None
        },
        "download-exomiser-remm-release": {
            "metavar": "Exomiser remm release",
            "help": """Release of ReMM (Regulatory Mendelian Mutation) database to download.\n"""
                    """If "default", "auto", or "config", retrieve from Application Properties file.\n"""
                    """Default: None""",
            "required": False,
            "default": None
        },
        "download-exomiser-remm-url": {
            "metavar": "Exomiser remm url",
            "help": """URL where ReMM (Regulatory Mendelian Mutation) database files can be downloaded from.\n"""
                    """Default: 'https://kircherlab.bihealth.org/download/ReMM'""",
            "required": False,
            "default": "https://kircherlab.bihealth.org/download/ReMM"
        },
        "download-exomiser-cadd-release": {
            "metavar": "Exomiser cadd release",
            "help": """Release of CADD (Combined Annotation Dependent Depletion) database to download.\n"""
                    """If "default", "auto", or "config", retrieve from Application Properties file.\n"""
                    """Default: None""",
            "required": False,
            "default": None
        },
        "download-exomiser-cadd-url": {
            "metavar": "Exomiser cadd url",
            "help": """URL where CADD (Combined Annotation Dependent Depletion) database files can be downloaded from.\n"""
                    """Default: 'https://kircherlab.bihealth.org/download/CADD'""",
            "required": False,
            "default": "https://kircherlab.bihealth.org/download/CADD"
        },
        "download-exomiser-cadd-url-snv-file": {
            "metavar": "Exomiser url snv",
            "help": """Name of the file containing the SNV (Single Nucleotide Variant) data\n"""
                    """for the CADD (Combined Annotation Dependent Depletion) database.\n"""
                    """Default: 'whole_genome_SNVs.tsv.gz'""",
            "required": False,
            "default": "whole_genome_SNVs.tsv.gz"
        },
        "download-exomiser-cadd-url-indel-file": {
            "metavar": "Exomiser cadd url indel",
            "help": """Name of the file containing the INDEL (Insertion-Deletion) data\n"""
                    """for the CADD (Combined Annotation Dependent Depletion) database.\n"""
                    """Default: 'InDels.tsv.gz'""",
            "required": False,
            "default": "InDels.tsv.gz"
        },

        # dbSNP
        "download-dbsnp": {
            "metavar": "dnSNP",
            "help": """Download dbSNP databases\n"""
                    """Folder where the dbSNP databases will be downloaded and stored.\n"""
                    """If the folder does not exist, it will be created.""",
            "required": False,
            "gooey": {
                "widget": "DirChooser"
            }
        },
        "download-dbsnp-releases": {
            "metavar": "dnSNP releases",
            "help": """Release of dbSNP to download\n"""
                    """Example: 'b152,b156'"""
                    """Default: 'b156'""",
            "required": False,
            "default": 'b156'
        },
        "download-dbsnp-release-default": {
            "metavar": "dnSNP release default",
            "help": """Default Release of dbSNP ('default' symlink)\n"""
                    """If None, first release to download will be assigned as dafault\n"""
                    """only if it does not exists\n"""
                    """Example: 'b156'"""
                    """Default: None (first releases by default)""",
            "required": False,
            "default": None
        },
        "download-dbsnp-url": {
            "metavar": "dbSNP url",
            "help": """URL where dbSNP database files can be downloaded from.\n"""
                    """Default: 'https://ftp.ncbi.nih.gov/snp/archive'""",
            "required": False,
            "default": "https://ftp.ncbi.nih.gov/snp/archive"
        },
        "download-dbsnp-url-files": {
            "metavar": "dbSNP url files",
            "help": """Dictionary that maps assembly names to specific dbSNP URL files.\n"""
                    """It allows you to provide custom dbSNP URL files for specific assemblies\n"""
                    """instead of using the default file naming convention\n"""
                    """Default: None""",
            "required": False,
            "default": None
        },
        "download-dbsnp-url-files-prefix": {
            "metavar": "dbSNP url files prefix",
            "help": """String that represents the prefix of the dbSNP file name for a specific assembly.\n"""
                    """It is used to construct the full URL of the dbSNP file to be downloaded.\n"""
                    """Default: 'GCF_000001405'""",
            "required": False,
            "default": "GCF_000001405"
        },
        "download-dbsnp-assemblies-map": {
            "metavar": "dbSNP assemblies map",
            "help": """dictionary that maps assembly names to their corresponding dbSNP versions.\n"""
                    """It is used to construct the dbSNP file name based on the assembly name.\n"""
                    """Default: {"hg19": "25", "hg38": "40"}""",
            "required": False,
            "default": {"hg19": "25", "hg38": "40"},
            "gooey": {
                "options": {
                    'initial_value': '{"hg19": "25", "hg38": "40"}'  
                }
            }
        },
        "download-dbsnp-vcf": {
            "help": """Generate well-formatted VCF from downloaded file:\n"""
                    """- Add and filter contigs associated to assembly\n"""
                    """- Normalize by splitting multiallelics """
                    """- Need genome (see --download-genome)""",
            "action": "store_true"
        },
        "download-dbsnp-parquet": {
            "help": """Generate Parquet file from VCF\n""",
            "action": "store_true"
        },

        # HGMD
        "convert-hgmd": {
            "metavar": "HGMD",
            "help": """Convert HGMD databases\n"""
                    """Folder where the HGMD databases will be stored.\n"""
                    """Fields in VCF, Parquet and TSV will be generated.\n"""
                    """If the folder does not exist, it will be created.""",
            "required": False,
            "gooey": {
                "widget": "DirChooser"
            }
        },
        "convert-hgmd-file": {
            "metavar": "HGMD file",
            "help": """File from HGMD\n"""
                    """Name format 'HGMD_Pro_<release>_<assembly>.vcf.gz'.""",
            "required": False,
            "gooey": {
                "widget": "FileChooser"
            }
        },
        "convert-hgmd-basename": {
            "metavar": "HGMD basename",
            "help": """File output basename\n"""
                    """Generated files will be prefixed by basename.\n"""
                    """Example: 'HGMD_Pro_MY_RELEASE'\n"""
                    """Default: Use input file name without '.vcf.gz'""",
            "required": False
        },

        # Databases parameters
        "generate-param": {
            "metavar": "param",
            "help": """Parameter file (JSON) with all databases found.\n"""
                    """Databases folders scanned are defined in config file.\n"""
                    """Structure of databases follow this structure (see doc):\n"""
                    """   .../<database>/<release>/<assembly>/*.[parquet|vcf.gz|...]""",
            "required": False,
            "gooey": {
                "widget": "FileSaver"
            }
        },
        "generate-param-description": {
            "metavar": "param description",
            "help": """Description file (JSON) with all databases found.\n"""
                    """Contains all databases with description of format, assembly, fields...""",
            "required": False,
            "gooey": {
                "widget": "FileSaver"
            }
        },
        "generate-param-releases": {
            "metavar": "param release",
            "help": """List of database folder releases to check\n"""
                    """Examples: 'current', 'latest'\n"""
                    """Default: 'current'""",
            "required": False,
            "default": "current"
        },
        "generate-param-formats": {
            "metavar": "param formats",
            "help": """List of database formats to check (e.g. parquet, vcf, bed, tsv...)\n"""
                    """Examples: 'parquet', 'parquet,vcf,bed,tsv'\n"""
                    """Default: 'parquet'""",
            "required": False,
            "default": "parquet"
        },
        "generate-param-bcftools": {
            "help": """Generate parameter file with BCFTools annotation for allowed formats\n"""
                    """Allowed formats with BCFTools: 'vcf', 'bed'""",
            "action": "store_true"
        },

        # From Annovar
        "annovar-code": {
            "metavar": "Annovar code",
            "help": """Annovar code, or database name. Usefull to name databases columns""",
            "required": False,
            "default": None
        },

        # Common
        "genomes-folder": {
            "metavar": "genomes",
            "help": """Folder containing genomes\n"""
                    f"""Default: {DEFAULT_GENOME_FOLDER}""",
            "required": False,
            "default": f"{DEFAULT_GENOME_FOLDER}",
            "gooey": {
                "widget": "DirChooser"
            }
        },

        # Shared
        "config": {
            "metavar": "config",
            "help": """Configuration file\n"""
                    """Default: {}""",
            "required": False,
            "default": "{}",
            "gooey": {
                "widget": "FileChooser",
                "options": {
                    'initial_value': ''  
                }
            }
        },
        "assembly": {
            "metavar": "assembly",
            "help": """Default assembly\n"""
                    """Default: 'hg19'""",
            "required": False,
            "default": "hg19"
        },
        "threads": {
            "metavar": "threads",
            "help": """Number of threads to use\n"""
                    """Use -1 to detect number of CPU/cores\n"""
                    """Default: -1""",
            "required": False,
            "default": -1
        },
        "memory": {
            "metavar": "memory",
            "help": """Memory to use (FLOAT[kMG])\n"""
                    """Default: None (80%% of RAM)""",
            "required": False,
            "default": None
        },
        "chunk_size": {
            "metavar": "chunk size",
            "help": """Number of records in batch to export output file.\n"""
                    """The lower the chunk size, the less memory consumption.\n"""
                    """For Parquet partitioning, files size will depend on the chunk size.\n"""
                    """default: 1000000""",
            "default": 1000000
        },
        "tmp": {
            "metavar": "tmp",
            "help": """Temporary folder.\n"""
                    """Especially for duckDB, default '.tmp' (see doc).\n"""
                    """default: None""",
            "default": None,
            "gooey": {
                "widget": "DirChooser"
            }
        },
        "duckdb_settings": {
            "metavar": "duckDB settings",
            "help": """DuckDB settings (see duckDB doc) as JSON (string or file).\n"""
                    """These settings have priority (see options 'threads', 'tmp'...).\n"""
                    """Examples: '{"TimeZone": "GMT", "temp_directory": "/tmp/duckdb", "threads": 8}'\n"""
                    """default: None""",
            "default": None,
            "gooey": {
                "widget": "FileChooser"
            }
        },
        "verbosity": {
            "metavar": "verbosity",
            "help": """Verbosity level\n"""
                    """Available: CRITICAL, ERROR, WARNING, INFO, DEBUG or NOTSET\n"""
                    """Default: INFO""",
            "required": False,
            "choices": ['CRITICAL', 'ERROR', 'WARNING', 'INFO', 'DEBUG', 'NOTSET'],
            "default": "INFO",
            "gooey": {
                "widget": "Dropdown",
                "options": {}
            }
        },
        "log": {
            "metavar": "log",
            "help": """Logs file\n"""
                    """Example: 'my.log'\n"""
                    """Default: None""",
            "default": None,
            "gooey": {
                "widget": "FileChooser"
            }
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


# Shared arguments
shared_arguments = ["config", "threads", "memory", "chunk_size", "tmp", "duckdb_settings", "verbosity", "log", "quiet", "verbose", "debug"]

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
                        """   howard process --input=tests/data/example.vcf.gz --annotations='snpeff' --calculations='snpeff_hgvs' --prioritizations=config/prioritization_profiles.json --explode_infos --output=/tmp/example.annotated.tsv --query='SELECT "#CHROM", POS, ALT, REF, snpeff_hgvs FROM variants' \n""", 

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
                "query": False,
                "explode_infos": False,
                "explode_infos_prefix": False,
                "explode_infos_fields": False,
                "include_header": False
            }
        }
    },
    "annotation": {
        "function" : "annotation",
        "description":  """Annotation is mainly based on a build-in Parquet annotation method, and tools such as BCFTOOLS, Annovar and snpEff. It uses available databases (see Annovar and snpEff) and homemade databases. Format of databases are: parquet, duckdb, vcf, bed, Annovar and snpEff (Annovar and snpEff databases are automatically downloaded, see howard databases tool). """,
        "help":         """Annotation of genetic variations using databases/files and tools.""",
        "epilog":       """Usage examples:\n"""
                        """   howard annotation --input=tests/data/example.vcf.gz --output=/tmp/example.howard.vcf.gz --annotations='tests/databases/annotations/hg19/avsnp150.parquet,tests/databases/annotations/hg19/dbnsfp42a.parquet,tests/databases/annotations/hg19/gnomad211_genome.parquet' \n"""
                        """   howard annotation --input=tests/data/example.vcf.gz --output=/tmp/example.howard.tsv --assembly=hg19 --annotations='annovar:refGene,annovar:cosmic70,snpeff,tests/databases/annotations/hg19/clinvar_20210123.parquet' \n"""
                        """   howard annotation --input=tests/data/example.vcf.gz --output=/tmp/example.howard.tsv --assembly=hg19 --annotations='ALL:parquet' \n""", 
        "groups": {
            "main": {
                "input": True,
                "output": True,
                "annotations": True,
                "assembly": False
            }
        }
    },
    "calculation": {
        "function" : "calculation",
        "description":  """Calculation processes variants information to generate new information, such as: identify variation type (VarType), harmonizes allele frequency (VAF) and calculate sttistics (VAF_stats), extracts Nomen (transcript, cNomen, pNomen...) from an HGVS field (e.g. snpEff, Annovar) with an optional list of personalized transcripts, generates VaRank format barcode, identify trio inheritance.""",
        "help":         """Calculation operations on genetic variations and genotype information.\n""",
        "epilog":       """Usage examples:\n"""
                        """   howard calculation --input=tests/data/example.full.vcf --output=/tmp/example.calculation.tsv --calculations='vartype' \n"""
                        """   howard calculation --input=tests/data/example.ann.vcf.gz --output=/tmp/example.calculated.tsv --calculations='snpeff_hgvs,NOMEN' --hgvs_field=snpeff_hgvs --transcripts=tests/data/transcripts.tsv \n"""
                        """   howard calculation --show_calculations \n"""
                        ,
        "groups": {
            "main": {
                "input": False,
                "output": False,
                "calculations": False,
                "calculation_config": False,
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
    "hgvs": {
        "function" : "hgvs",
        "description":  """HGVS annotation using HUGO HGVS internation Sequence Variant Nomenclature (http://varnomen.hgvs.org/). Annotation refere to refGene and genome to generate HGVS nomenclature for all available transcripts. This annotation add 'hgvs' field into VCF INFO column of a VCF file.""",
        "help":         """HGVS annotation.\n""",
        "epilog":       """Usage examples:\n"""
                        """   howard hgvs --input=tests/data/example.full.vcf --output=/tmp/example.hgvs.vcf \n"""
                        ,
        "groups": {
            "main": {
                "input": True,
                "output": False,
                "assembly": False
            },
            "HGVS": {
                "use_gene": False,
                "use_exon": False,
                "use_protein": False,
                "add_protein": False,
                "full_format": False,
                "codon_type": False
            },
            "Databases": {
                "refgene": False,
                "refseqlink": False,
                "genomes-folder": False
            }
        }
    },
    "prioritization": {
        "function" : "prioritization",
        "description":  """Prioritization algorithm uses profiles to flag variants (as passed or filtered), calculate a prioritization score, and automatically generate a comment for each variants (example: 'polymorphism identified in dbSNP. associated to Lung Cancer. Found in ClinVar database'). Prioritization profiles are defined in a configuration file in JSON format. A profile is defined as a list of annotation/value, using wildcards and comparison options (contains, lower than, greater than, equal...). Annotations fields may be quality values (usually from callers, such as 'DP') or other annotations fields provided by annotations tools, such as HOWARD itself (example: COSMIC, Clinvar, 1000genomes, PolyPhen, SIFT). Multiple profiles can be used simultaneously, which is useful to define multiple validation/prioritization levels (example: 'standard', 'stringent', 'rare variants', 'low allele frequency').\n""",
        "help": "Prioritization of genetic variations based on annotations criteria (profiles).",
        "epilog": """Usage examples:\n"""
                        """   howard prioritization --input=tests/data/example.vcf.gz --output=/tmp/example.prioritized.vcf.gz --prioritizations=config/prioritization_profiles.json --profiles='default,GERMLINE' \n""", 
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
                        """   howard query --input=tests/data/example.vcf.gz --explode_infos --query='SELECT "#CHROM", POS, REF, ALT, DP, CLNSIG, sample2, sample3 FROM variants WHERE DP >= 50 OR CLNSIG NOT NULL ORDER BY DP DESC' \n"""
                        """   howard query --query="SELECT \\\"#CHROM\\\", POS, REF, ALT, \\\"INFO/Interpro_domain\\\" FROM 'tests/databases/annotations/hg19/dbnsfp42a.parquet' WHERE \\\"INFO/Interpro_domain\\\" NOT NULL ORDER BY \\\"INFO/SiPhy_29way_logOdds_rankscore\\\" DESC LIMIT 10" \n"""
                        """   howard query --explode_infos --explode_infos_prefix='INFO/' --query="SELECT \\\"#CHROM\\\", POS, REF, ALT, STRING_AGG(INFO, ';') AS INFO FROM 'tests/databases/annotations/hg19/*.parquet' GROUP BY \\\"#CHROM\\\", POS, REF, ALT" --output=/tmp/full_annotation.tsv  && head -n2 /tmp/full_annotation.tsv \n"""
                        , 
        "groups": {
            "main": {
                "input": False,
                "output": False,
                "query": True
            },
            "Explode infos": {
                "explode_infos": False,
                "explode_infos_prefix": False,
                "explode_infos_fields": False,
            },
            "Query": {
                "query_limit": False,
                "query_print_mode": False,
            },
            "Output": {
                "include_header": False
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
                "input": True,
                "stats_md": False,
                "stats_json": False
            }
        }
    },
    "convert": {
        "function" : "convert",
        "description":  """Convert genetic variations file to another format. Multiple format are available, such as usual and official VCF and BCF format, but also other formats such as TSV, CSV, PSV and Parquet/duckDB. These formats need a header '.hdr' file to take advantage of the power of howard (especially through INFO/tag definition), and using howard convert tool automatically generate header file fo futher use. """,
        "help": "Convert genetic variations file to another format.",
        "epilog": """Usage examples:\n"""
                        """   howard convert --input=tests/data/example.vcf.gz --output=/tmp/example.tsv \n"""
                        """   howard convert --input=tests/data/example.vcf.gz --output=/tmp/example.parquet \n"""
                        """   howard convert --input=tests/data/example.vcf.gz --output=/tmp/example.tsv --explode_infos --explode_infos_fields='CLNSIG,SIFT,DP' --order_by='CLNSIG DESC, DP DESC' \n"""
                        """   howard convert --input=tests/data/example.vcf.gz --output=/tmp/example.tsv --explode_infos --explode_infos_prefix='INFO/' --explode_infos_fields='CLNSIG,SIFT,DP,*' --order_by='"INFO/CLNSIG" DESC, "INFO/DP" DESC' --include_header """,
        "groups": {
            "main": {
                "input": True,
                "output": True,
                "explode_infos": False,
                "explode_infos_prefix": False,
                "explode_infos_fields": False,
                "order_by": False,
                "include_header": False,
                "parquet_partitions": False
            }
        }
    },
    "databases": {
        "function" : "databases",
        "description": """Download databases and needed files for howard and associated tools""",
        "help": """Download databases and needed files for howard and associated tools""",
        "epilog": """Usage examples:\n"""
                    """   howard databases --assembly=hg19 --download-genomes=/databases/genomes/current --download-genomes-provider=UCSC --download-genomes-contig-regex='chr[0-9XYM]+$' --download-annovar=/databases/annovar/current --download-annovar-files='refGene,cosmic70,nci60' --download-snpeff=/databases/snpeff/current --download-refseq=/databases/refseq/current --download-refseq-format-file='ncbiRefSeq.txt' --download-dbnsfp=/databases/dbnsfp/current --download-dbnsfp-release='4.4a' --download-dbnsfp-subdatabases --download-alphamissense=/databases/alphamissense/current --download-exomiser=/databases/exomiser/current --download-dbsnp=/databases/dbsnp/current --download-dbsnp-vcf --threads=8\n"""
                    """   howard databases --generate-param=/tmp/param.json --generate-param-description=/tmp/test.description.json --generate-param-formats=parquet """
                    """\n"""
                    """Notes:\n"""
                    """   - Downloading databases can take a while, depending on network, threads and memory\n"""
                    """   - Proxy: Beware of network and proxy configuration\n"""
                    """   - dbNSFP download: More threads, more memory usage (8 threads ~ 16Gb, 24 threads ~ 32Gb)\n"""
                    ,
        "groups": {
            "main": {
                "assembly": False,
                "genomes-folder": False
            },
            "Genomes": {
                "download-genomes": False,
                "download-genomes-provider": False,
                "download-genomes-contig-regex": False,
            },
            "snpEff": {
                "download-snpeff": False
            },
            "Annovar": {
                "download-annovar": False,
                "download-annovar-files": False,
                "download-annovar-url": False
            },
            "refSeq": {
                "download-refseq": False,
                "download-refseq-url": False,
                "download-refseq-prefix": False,
                "download-refseq-files": False,
                "download-refseq-format-file": False,
                "download-refseq-include-utr5": False,
                "download-refseq-include-utr3": False,
                "download-refseq-include-chrM": False,
                "download-refseq-include-non-canonical-chr": False,
                "download-refseq-include-non-coding-transcripts": False,
                "download-refseq-include-transcript-version": False,
            },
            "dbNSFP": {
                "download-dbnsfp": False,
                "download-dbnsfp-url": False,
                "download-dbnsfp-release": False,
                "download-dbnsfp-parquet-size": False,
                "download-dbnsfp-subdatabases": False,
                "download-dbnsfp-parquet": False,
                "download-dbnsfp-vcf": False,
                "download-dbnsfp-no-files-all": False,
                "download-dbnsfp-add-info": False,
                "download-dbnsfp-row-group-size": False
            },
            "AlphaMissense": {
                "download-alphamissense": False,
                "download-alphamissense-url": False
            },
            "Exomiser": {
                "download-exomiser": False,
                "download-exomiser-application-properties": False,
                "download-exomiser-url": False,
                "download-exomiser-release": False,
                "download-exomiser-phenotype-release": False,
                "download-exomiser-remm-release": False,
                "download-exomiser-remm-url": False,
                "download-exomiser-cadd-release": False,
                "download-exomiser-cadd-url": False,
                "download-exomiser-cadd-url-snv-file": False,
                "download-exomiser-cadd-url-indel-file": False
            },
            "dbSNP": {
                "download-dbsnp": False,
                "download-dbsnp-releases": False,
                "download-dbsnp-release-default": False,
                "download-dbsnp-url": False,
                "download-dbsnp-url-files": False,
                "download-dbsnp-url-files-prefix": False,
                "download-dbsnp-assemblies-map": False,
                "download-dbsnp-vcf": False,
                "download-dbsnp-parquet": False
            },
            "HGMD": {
                "convert-hgmd": False,
                "convert-hgmd-file": False,
                "convert-hgmd-basename": False
            },
            "Parameters file": {
                "generate-param": False,
                "generate-param-description": False,
                "generate-param-releases": False,
                "generate-param-formats": False,
                "generate-param-bcftools": False
            }

        }


    },
    "from_annovar": {
        "function" : "from_annovar",
        "description": """(beta) Formatting Annovar database file to other format (VCF and Parquet). Exported Parquet file includes INFO/tags columns as VCF INFO columns had been exploded""",
        "help": """(beta) Formatting Annovar database file to other format (VCF and Parquet)""",
        "epilog": """Usage examples:\n"""
                    """   howard from_annovar --input=tests/databases/others/hg19_nci60.txt --output=/tmp/nci60.from_annovar.vcf.gz --to_parquet=/tmp/nci60.from_annovar.parquet --annovar-code=nci60 --genome=/databases/genomes/current/hg19.fa --config=/tool/config/config.json --threads=8 """, 
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
            },
            "Modes": {
                "reduce_memory": False,
                "multi_variant": False
            }
        }
    }
}


# get argument
def get_argument(arguments:dict = {}, arg:str = "", required:bool = False, remove_infos:list = ["gooey"], add_metavar:bool = False) -> dict:
    """
    The `get_argument` function retrieves information about a specific argument from a dictionary, and
    can also set its "required" status.
    
    :param arguments: A dictionary containing information about the arguments passed to a function or
    method
    :type arguments: dict
    :param arg: The `arg` parameter is a string that represents the name of the argument that you want
    to retrieve information for
    :type arg: str
    :param required: The `required` parameter is a boolean value that determines whether the argument is
    required or not. If set to True, the function will return an empty dictionary if the argument is not
    found in the `arguments` dictionary. If set to False (default), the function will still return an
    empty dictionary if, defaults to False
    :type required: bool (optional)
    :param remove_infos: The `remove_infos` parameter is a list that contains the names of specific
    information that you want to remove from the argument dictionary. In the code, it is used to remove
    specific argument information such as "gooey" from the `arg_infos` dictionary
    :type remove_infos: list
    :return: a dictionary containing information about a specific argument, specified by the `arg`
    parameter. If the argument is found in the `arguments` dictionary, the function returns a dictionary
    containing the information about that argument. If the argument is not found, an empty dictionary is
    returned. The `required` parameter is used to specify whether the argument is required or not, and
    this information is added to
    """

    if arg in arguments:
        arg_infos = arguments.get(arg,{}).copy()
        for arg_info in remove_infos:
            arg_infos.pop(arg_info, None)
        if required != None:
            arg_infos["required"] = required
        if add_metavar and "metavar" not in arg_infos:
            arg_infos["metavar"] = arg.replace("_", " ")
        return arg_infos
    else:
        return {}


# get_argument_gooey
def get_argument_gooey(arg:str):
    """
    The function `get_argument_gooey` takes an argument and returns the corresponding widget and options
    for the Gooey library in Python.
    
    :param arg: The `arg` parameter is a string that represents the name of the argument you want to
    retrieve information for
    :type arg: str
    :return: The function `get_argument_gooey` returns two values: `widget` and `options`.
    """

    # Init
    argument = get_argument(arguments=arguments, arg=arg, remove_infos=[])
    argument_type = argument.get("type", None)
    gooey_argument = argument.get("gooey", {})

    # Widget
    widget = None
    if gooey_argument.get("widget", None):
        widget = gooey_argument.get("widget")
    else:
        if str(argument_type) == "FileType('r')":
            widget = "FileChooser"
        elif str(argument_type) == "FileType('w')":
            widget = "FileSaver"
    
    # options
    options = gooey_argument.get("options", {})

    # Return
    return widget, options