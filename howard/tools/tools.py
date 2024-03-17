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
import importlib

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
from howard.tools.help import *
from howard.tools.from_annovar import *


# Import gui only if gooey and wx is installed
try:
    check_gooey = importlib.util.find_spec("gooey")
    check_wx = importlib.util.find_spec("wx")
    tool_gui_enable = check_gooey and check_wx
except ImportError:
    tool_gui_enable = False

if tool_gui_enable:
    from howard.tools.gui import *


class PathType(object):

    def __init__(self, exists=True, type='file', dash_ok=True):
        '''exists:
                True: a path that does exist
                False: a path that does not exist, in a valid parent directory
                None: don't care
           type: file, dir, symlink, None, or a function returning True for valid paths
                None: don't care
           dash_ok: whether to allow "-" as stdin/stdout'''

        self.__name__ = "Path"

        assert exists in (True, False, None)
        assert type in ('file','dir','symlink',None) or hasattr(type,'__call__')

        self._exists = exists
        self._type = type
        self._dash_ok = dash_ok

    def __call__(self, string):

        # Full path if not a JSON string
        try:
            json.loads(string)
        except:
            string = full_path(string)

        if string=='-':
            # the special argument "-" means sys.std{in,out}
            if self._type == 'dir':
                raise ValueError('standard input/output (-) not allowed as directory path')
            elif self._type == 'symlink':
                raise ValueError('standard input/output (-) not allowed as symlink path')
            elif not self._dash_ok:
                raise ValueError('standard input/output (-) not allowed')
        else:
            e = os.path.exists(string)
            if self._exists==True:
                if not e:
                    raise ValueError("path does not exist: '%s'" % string)

                if self._type is None:
                    pass
                elif self._type=='file':
                    if not os.path.isfile(string):
                        raise ValueError("path is not a file: '%s'" % string)
                elif self._type=='symlink':
                    if not os.path.symlink(string):
                        raise ValueError("path is not a symlink: '%s'" % string)
                elif self._type=='dir':
                    if not os.path.isdir(string):
                        raise ValueError("path is not a directory: '%s'" % string)
                elif not self._type(string):
                    raise ValueError("path not valid: '%s'" % string)
            else:
                if self._exists==False and e:
                    raise ValueError("path exists: '%s'" % string)

        return string

# Arguments dict
arguments = {

        # Process & other
        "input": {
            "metavar": "input",
            "help": """Input file path.\n"""
                    """Format file must be either VCF, Parquet, TSV, CSV, PSV or duckDB.\n"""
                    """Files can be compressesd (e.g. vcf.gz, tsv.gz).\n""",
            "required": False,
            "default": None,
            "type": PathType(exists=True, type=None),
            "gooey": {
                "widget": "FileChooser",
                "options": {
                    "wildcard":
                        "Parquet file (*.parquet)|*.parquet|"
                        "All files (*)|*"
                }
            }
        },
        "output": {
            "metavar": "output",
            "help": """Output file path.\n"""
                    """Format file must be either VCF, Parquet, TSV, CSV, PSV or duckDB.\n"""
                    """Files can be compressesd (e.g. vcf.gz, tsv.gz).\n""",
            "required": False,
            "default": None,
            "type": PathType(exists=None, type=None),
            "gooey": {
                "widget": "FileSaver"
            }
        },
        "param": {
            "metavar": "param",
            "help": """Parameters JSON file or JSON string.\n""",
            "default": "{}",
            "type": PathType(exists=None, type=None),
            "gooey": {
                "widget": "FileChooser",
                "options": {
                    'initial_value': '',
                    "wildcard":
                        "JSON file (*.json)|*.json|"
                        "All files (*)|*"
                }
            }
        },
        "query": {
            "metavar": "query",
            "help": """Query in SQL format\n"""
                    """(e.g. 'SELECT * FROM variants LIMIT 50').\n""",
            "default": None,
            "type": str,
            "gooey": {
                "widget": "Textarea",
                "options": {
                    'initial_value': 'SELECT * FROM variants'  
                }
            },
            "extra": {
                "param_section": "query"
            }
        },
        "output_query": {
            "metavar": "output",
            "help": """Output Query file.\n"""
                    """Format file must be either VCF, Parquet, TSV, CSV, PSV or duckDB.\n""",
            "default": None,
            "type": PathType(exists=None, type=None),
            "gooey": {
                "widget": "FileSaver",
                "options": {
                    "wildcard":
                        "All files (*)|*",
                }
            }
        },

        # Annotations
        "annotations": {
            "metavar": "annotations",
            "help": """Annotation with databases files, or with tools,\n"""
                    """as a list of files in Parquet, VCF, BED, or keywords\n"""
                    """ (e.g. '/path/to/file.parquet,bcftools:file2.vcf.gz,annovar:refGene,snpeff').\n"""
                    """For a Parquet/VCF/BED file, use file paths\n"""
                    """ (e.g. '/path/to/file.parquet', 'file1.parquet,file2.vcf.gz').\n"""
                    """For BCFTools anotation, use keyword 'bcftools' with file paths\n"""
                    """ (e.g. 'bcftools:/path/to/file.vcf.gz:/path/to/file.bed.gz').\n"""
                    """For Annovar annotation, use keyword 'annovar' with annovar code\n"""
                    """ (e.g. 'annovar:refGene', 'annovar:cosmic70').\n"""
                    """For snpeff annotation, use keyword 'snpeff'.\n"""
                    """For Exomiser annotation, use keyword 'exomiser' with options as key=value\n"""
                    """ (e.g. 'exomiser:preset=exome:transcript_source=refseq').\n"""
                    """For add all availalbe databases files, use 'ALL' keyword,\n"""
                    """ with filters on type and release\n"""
                    """ (e.g. 'ALL', 'ALL:parquet:current', 'ALL:parquet,vcf:current,devel').\n""",
            "default": None,
            "type": str,
            "gooey": {
                "widget": "MultiFileChooser",
                "options": {
                    "default_dir": DEFAULT_ANNOTATIONS_FOLDER,
                    "message": "Database files"
                }
            },
            "extra": {
                "format": "DB[,DB]*[,bcftools:DB[:DB]*][,annovar:KEY[:KEY]*][,snpeff][,exomiser[:var=val]*]",
                "examples": {
                    "Parquet method annotation with 2 Parquet files":
                       "\"annotations\": \"/path/to/database1.parquet,/path/to/database2.parquet\"",
                    "Parquet method annotation with multiple file formats":
                       "\"annotations\": \"/path/to/database1.parquet,/path/to/database2.vcf.gz,/path/to/database2.bed.gz\"",
                    "Parquet method annotation with available Parquet databases in current release (check databases in production)":
                       "\"annotations\": \"ALL:parquet:current\"",
                    "Parquet method annotation with available Parquet databases in latest release (check databases before production)":
                       "\"annotations\": \"ALL:parquet:latest\"",
                    "Annovation with BCFTools":
                       "\"annotations\": \"bcftools:/path/to/database2.vcf.gz:/path/to/database2.bed.gz\"",
                    "Annovation with Annovar (refGene with hgvs and Cosmic)":
                       "\"annotations\": \"annovar:refGene:cosmic70\"",
                    "Annovation with snpEff (default options)":
                       "\"annotations\": \"snpeff\"",
                    "Annovation with Exomiser with options":
                       "\"annotations\": \"exomiser:preset=exome:hpo=0001156+0001363+0011304+0010055:transcript_source=refseq:release=2109\"",
                    "Multiple tools annotations (Parquet method, BCFTools, Annovar, snpEff and Exomiser)":
                       "\"annotations\": \"/path/to/database1.parquet,bcftools:/path/to/database2.vcf.gz,annovar:refGene:cosmic70,snpeff,exomiser:preset=exome:transcript_source=refseq\""
                }
            }
        },
        "annotations_update": {
            "help": """Update option for annotation (Only for Parquet annotation).\n"""
                    """If True, annotation fields will be removed and re-annotated.\n"""
                    """These options will be applied to all annotation databases.\n""",
            "action": "store_true",
            "default": False,
            "gooey": {
                "widget": "BlockCheckbox",
                "options": {
                    'checkbox_label': "Update annotation method"
                }
            },
            "extra": {
                "param_section": "annotation:options"
            }
        },
        "annotations_append": {
            "help": """Append option for annotation (Only for Parquet annotation).\n"""
                    """If True, annotation fields will be annotated only if not annotation exists for the variant.\n"""
                    """These options will be applied to all annotation databases.\n""",
            "action": "store_true",
            "default": False,
            "gooey": {
                "widget": "BlockCheckbox",
                "options": {
                    'checkbox_label': "Append annotation method"
                }
            },
            "extra": {
                "param_section": "annotation:options"
            }
        },

        # Calculations
        "calculations": {
            "metavar": "operations",
            "help": """Quick calculations on genetic variants information and genotype information,\n"""
                    """as a list of operations (e.g. 'VARTYPE,variant_id').\n"""
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
                    """\n""",
            "default": None,
            "type": str
        },

        # Prioritizations
        "prioritizations": {
            "metavar": "prioritisations",
            "help": "Prioritization file in JSON format (defines profiles, see doc).\n",
            "default": None,
            "type": PathType(exists=True, type='file'),
            "gooey": {
                "widget": "FileChooser",
                "options": {
                    "wildcard":
                        "JSON file (*.json)|*.json|"
                        "All files (*)|*"
                }
            }
        },
        "profiles": {
            "metavar": "profiles",
            "help": """List of prioritization profiles to process (based on Prioritization JSON file),\n"""
                    """such as 'default', 'rare variants', 'low allele frequency', 'GERMLINE'.\n"""
                    """By default, all profiles available will be processed.\n""",
            "default": None,
            "type": str
        },
        "default_profile": {
            "metavar": "default profile",
            "help": """Prioritization profile by default (see doc).\n"""
                    """Default is the first profile in the list of prioritization profiles.\n""",
            "default": None,
            "type": str
        },
        "pzfields": {
            "metavar": "pzfields",
            "help": """Prioritization fields to provide (see doc).\n"""
                    """Available: PZScore, PZFlag, PZTags, PZComment, PZInfos\n""",
            "default": "PZScore,PZFlag",
            "type": str
        },
        "prioritization_score_mode": {
            "metavar": "prioritization score mode",
            "help": """Prioritization Score mode (see doc).\n"""
                    """Available: HOWARD (increment score), VaRank (max score)\n""",
            "default": 'HOWARD',
            "type": str,
            "choices": ["HOWARD", "VaRank"],
            "gooey": {
                "widget": "Dropdown",
                "options": {}
            }
        },

        # Query print options
        "query_limit": {
            "metavar": "query limit",
            "help": """Limit of number of row for query (only for print result, not output).\n""",
            "default": 10,
            "type": int,
            "gooey": {
                "widget": "IntegerField",
                "options": {
                    'min': 1, 
                    'max': 10000,
                    'increment': 10
                }
            }
        },
        "query_print_mode": {
            "metavar": "print mode",
            "help": """Print mode of query result (only for print result, not output).\n"""
                    """Either None (native), 'markdown' or 'tabulate'.\n""",
            "choices": [None, "markdown", "tabulate"],
            "default": None,
            "type": str,
            "gooey": {
                "widget": "Dropdown",
                "options": {}
            }
        },

        # Explode infos
        "explode_infos": {
            "help": """Explode VCF INFO/Tag into 'variants' table columns.\n""",
            "action": "store_true",
            "default": False
        },
        "explode_infos_prefix": {
            "metavar": "explode infos prefix",
            "help": """Explode VCF INFO/Tag with a specific prefix.\n""",
            "default": "",
            "type": str
        },
        "explode_infos_fields": {
            "metavar": "explode infos list",
            "help": """Explode VCF INFO/Tag specific fields/tags.\n"""
                    """Keyword `*` specify all available fields, except those already specified.\n"""
                    """Pattern (regex) can be used, such as `.*_score` for fields named with '_score' at the end.\n"""
                    """Examples:\n"""
                    """- 'HGVS,SIFT,Clinvar' (list of fields)\n"""
                    """- 'HGVS,*,Clinvar' (list of fields with all other fields at the end)\n"""
                    """- 'HGVS,.*_score,Clinvar' (list of 2 fields with all scores in the middle)\n"""
                    """- 'HGVS,.*_score,*' (1 field, scores, all other fields)\n"""
                    """- 'HGVS,*,.*_score' (1 field, all other fields, all scores)\n""",
            "default": "*",
            "type": str
        },

        # Include header
        "include_header": {
            "help": """Include header (in VCF format) in output file.\n"""
                    """Only for compatible formats (tab-delimiter format as TSV or BED).\n""",
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
                    """Examples: 'ACMG_score DESC', 'PZFlag DESC, PZScore DESC'.\n""",
            "default": "",
            "type": str,
            "extra": {
                "examples": {
                    "Order by ACMG score in descending order":
                        """ACMG_score DESC""",
                    "Order by PZFlag and PZScore in descending order":
                        """PZFlag DESC, PZScore DESC"""
                }
            }
        },

        # Parquet partition
        "parquet_partitions": {
            "metavar": "parquet partitions",
            "help": """Parquet partitioning using hive (available for any format).\n"""
                    """This option is faster parallel writing, but memory consuming.\n"""
                    """Use 'None' (string) for NO partition but split parquet files into a folder.\n"""
                    """Examples: '#CHROM', '#CHROM,REF', 'None'.\n""",
            "default": None,
            "type": str
        },
        
        # From annovar
        "multi_variant": {
            "metavar": "multi variant",
            "help": """Variant with multiple annotation lines.\n"""
                    """Either 'auto' (auto-detection), 'enable' or 'disable'.\n""",
            "default": "auto",
            "type": str,
            "choices": ["auto", "enable", "disable"],
            "gooey": {
                "widget": "Dropdown",
                "options": {}
            }
        },
        "reduce_memory": {
            "metavar": "reduce memory",
            "help": """Reduce memory option,\n"""
                    """either 'auto' (auto-detection), 'enable' or 'disable'.\n""",
            "default": "auto",
            "type": str,
            "choices": ["auto", "enable", "disable"],
            "gooey": {
                "widget": "Dropdown",
                "options": {}
            }
        },

        # Calculation
        "calculation_config": {
            "metavar": "calculation config",
            "help": """Calculation configuration JSON file.\n""",
            "default": None,
            "type": PathType(exists=True, type='file'),
            "gooey": {
                "widget": "FileChooser",
                "options": {
                    "wildcard":
                        "JSON file (*.json)|*.json|"
                        "All files (*)|*"
                }
            }
        },
        "show_calculations": {
            "help": """Show available calculation operations.\n""",
            "action": "store_true",
            "default": False
        },
        "hgvs_field": {
            "metavar": "HGVS field",
            "help": """HGVS INFO/tag containing a list o HGVS annotations.\n""",
            "default": "hgvs",
            "type": str,
            "extra": {
                "param_section": "calculation:NOMEN:options"
            }
        },
        "transcripts": {
            "metavar": "transcripts",
            "help": """Transcripts TSV file,\n"""
                    """with Transcript in first column, optional Gene in second column.\n""",
            "default": None,
            "type": PathType(exists=True, type='file'),
            "gooey": {
                "widget": "FileChooser",
                "options": {
                    "wildcard":
                        "TSV file (*.tsv)|*.tsv|"
                        "All files (*)|*"
                }
            },
            "extra": {
                "param_section": "calculation:NOMEN:options"
            }
        },
        "trio_pedigree": {
            "metavar": "trio pedigree",
            "help": """Pedigree Trio for trio inheritance calculation.\n"""
                    """either a JSON file or JSON string"""
                    """(e.g. '{"father": "sample1", "mother": "sample2", "child": "sample3"}').\n""",
            "default": None,
            #"type": PathType(exists=True, type='file'),
            "gooey": {
                "widget": "FileChooser",
                "options": {
                    "wildcard":
                        "JSON file (*.json)|*.json|"
                        "All files (*)|*"
                }
            },
            "extra": {
                "param_section": "calculation:TRIO"
            }
        },

        # Stats
        "stats_md": {
            "metavar": "stats markdown",
            "help": """Stats Output file in MarkDown format.\n""",
            "required": False,
            "default": None,
            "type": PathType(exists=None, type='file'),
            "gooey": {
                "widget": "FileSaver",
                "options": {
                    "wildcard":
                        "Markdown file (*.md)|*.md"
                }
            },
            "extra": {
                "examples": {
                    "Export statistics in Markdown format":
                        """"stats_md": '/tmp/stats.md'"""
                }
            }
        },
        "stats_json": {
            "metavar": "stats json",
            "help": """Stats Output file in JSON format.\n""",
            "required": False,
            "default": None,
            "type": PathType(exists=None, type='file'),
            "gooey": {
                "widget": "FileSaver",
                "options": {
                    "wildcard":
                        "JSON file (*.json)|*.json"
                }
            },
            "extra": {
                "examples": {
                    "Export statistics in JSON format":
                        """"stats_json": '/tmp/stats.json'"""
                }
            }
        },

        # Assembly and Genome
        "assembly": {
            "metavar": "assembly",
            "help": """Genome Assembly (e.g. 'hg19', 'hg38').\n""",
            "required": False,
            "default": DEFAULT_ASSEMBLY,
            "type": str,
            "extra": {
                "examples": {
                    "Default assembly for all analysis tools":
                        """"assembly": 'hg19'""",
                    "List of assemblies for databases download tool":
                        """"assembly": 'hg19,hg38'"""
                }
            }
        },
        "genome": {
            "metavar": "genome",
            "help": """Genome file in fasta format (e.g. 'hg19.fa', 'hg38.fa').\n""",
            "required": False,
            "default": "~/howard/databases/genomes/current/hg19/hg19.fa",
            "type": PathType(exists=True, type='file'),
            "gooey": {
                "widget": "FileChooser",
                "options": {
                    "wildcard":
                        "All files (*)|*"
                }
            }
        },
        
        # HGVS
        "hgvs_options": {
            "metavar": "HGVS options",
            "help": """Quick HGVS annotation options.\n"""
                    """This option will skip all other hgvs options.\n"""
                    """Examples:\n"""
                    """- 'default' (for default options)\n"""
                    """- 'full_format' (for full format HGVS annotation)\n"""
                    """- 'use_gene=True:add_protein=true:codon_type=FULL'\n""",
            "required": False,
            "default": None,
            "type": str
        },
        "use_gene": {
            "help": """Use Gene information to generate HGVS annotation\n"""
                    """(e.g. 'NM_152232(TAS1R2):c.231T>C')""",
            "action": "store_true",
            "default": False
        },
        "use_exon": {
            "help": """Use Exon information to generate HGVS annotation\n"""
                    """(e.g. 'NM_152232(exon2):c.231T>C').\n"""
                    """Only if 'use_gene' is not enabled.\n""",
            "action": "store_true",
            "default": False
        },
        "use_protein": {
            "help": """Use Protein level to generate HGVS annotation\n"""
                    """(e.g. 'NP_689418:p.Cys77Arg').\n"""
                    """Can be used with 'use_exon' or 'use_gene'.\n""",
            "action": "store_true",
            "default": False
        },
        "add_protein": {
            "help": """Add Protein level to DNA HGVS annotation """
                    """(e.g 'NM_152232:c.231T>C,NP_689418:p.Cys77Arg').\n""",
            "action": "store_true",
            "default": False
        },
        "full_format": {
            "help": """Generates HGVS annotation in a full format\n"""
                    """by using all information to generates an exhaustive annotation\n"""
                    """(non-standard, e.g. 'TAS1R2:NM_152232:NP_689418:c.231T>C:p.Cys77Arg').\n"""
                    """Use 'use_exon' to add exon information\n"""
                    """(e.g 'TAS1R2:NM_152232:NP_689418:exon2:c.231T>C:p.Cys77Arg').\n""",
            "action": "store_true",
            "default": False
        },
        "use_version": {
            "help": """Generates HGVS annotation with transcript version\n"""
                    """(e.g. 'NM_152232.1:c.231T>C').\n""",
            "action": "store_true",
            "default": False
        },
        "codon_type": {
            "metavar": "Codon type",
            "help": """Amino Acide Codon format type to use to generate HGVS annotation.\n"""
                    """Available:\n"""
                    """- '1': codon in 1 character (e.g. 'C', 'R')\n"""
                    """- '3': codon in 3 character (e.g. 'Cys', 'Arg')\n"""
                    """-'FULL': codon in full name (e.g. 'Cysteine', 'Arginine')\n""",
            "required": False,
            "default": "3",
            "type": str,
            "choices": ["1", "3", "FULL"],
            "gooey": {
                "widget": "Dropdown",
                "options": {}
            }
        },
        "refgene": {
            "metavar": "refGene",
            "help": """Path to refGene annotation file.\n""",
            "required": False,
            "default": None,
            "type": PathType(exists=True, type='file'),
            "gooey": {
                "widget": "FileChooser",
                "options": {
                    "wildcard":
                        "All files (*)|*",
                    "default_dir": DEFAULT_REFSEQ_FOLDER,
                    "default_file": "ncbiRefSeq.txt",
                    "message": "Path to refGene annotation file"
                }
            }
        },
        "refseqlink": {
            "metavar": "refSeqLink",
            "help": """Path to refSeqLink annotation file.\n""",
            "required": False,
            "default": None,
            "type": PathType(exists=True, type='file'),
            "gooey": {
                "widget": "FileChooser",
                "options": {
                    "wildcard":
                        "All files (*)|*",
                    "default_dir": DEFAULT_REFSEQ_FOLDER,
                    "default_file": "ncbiRefSeq.txt",
                    "message": "Path to refGeneLink annotation file"
                }
            }
        },
        "refseq-folder": {
            "metavar": "refseq folder",
            "help": """Folder containing refSeq files.\n""",
            "required": False,
            "default": DEFAULT_REFSEQ_FOLDER,
            "type": PathType(exists=True, type='dir'),
            "gooey": {
                "widget": "DirChooser",
                "options": {
                    "default_dir": DEFAULT_REFSEQ_FOLDER,
                    "message": "Path to refGenefolder"
                }
            }
        },
        
        # Databases

        # Genome
        "download-genomes": {
            "metavar": "genomes",
            "help": """Path to genomes folder\n"""
                    """with Fasta files, indexes,\n"""
                    """and all files generated by pygenome module.\n"""
                    f"""(e.g. '{DEFAULT_GENOME_FOLDER}').\n""",
            "required": False,
            "default": None,
            "type": PathType(exists=None, type='dir'),
            "gooey": {
                "widget": "DirChooser",
                "options": {
                    "default_dir": DEFAULT_DATABASE_FOLDER,
                    "message": "Path to genomes folder"
                }
            }
        },
        "download-genomes-provider": {
            "metavar": "genomes provider",
            "help": """Download Genome from an external provider.\n"""
                    """Available: GENCODE, Ensembl, UCSC, NCBI.\n""",
            "required": False,
            "default": "UCSC",
            "type": str,
            "choices": ["GENCODE", "Ensembl", "UCSC", "NCBI"],
            "gooey": {
                "widget": "Dropdown",
                "options": {}
            }
        },
        "download-genomes-contig-regex": {
            "metavar": "genomes contig regex",
            "help": """Regular expression to select specific chromosome\n"""
                    """(e.g 'chr[0-9XYM]+$').\n""",
            "required": False,
            "default": None,
            "type": str
        },

        # Annovar
        "download-annovar": {
            "metavar": "Annovar",
            "help": """Path to Annovar databases\n"""
                    f"""(e.g. '{DEFAULT_ANNOVAR_FOLDER}').\n""",
            "required": False,
            "type": PathType(exists=None, type='dir'),
            "default": None,
            "gooey": {
                "widget": "DirChooser",
                "options": {
                    "default_dir": DEFAULT_DATABASE_FOLDER,
                    "message": "Path to Annovar databases folder"
                }
            }
        },
        "download-annovar-files": {
            "metavar": "Annovar code",
            "help": """Download Annovar databases for a list of Annovar file code (see Annovar Doc).\n"""
                    """Use None to donwload all available files,\n"""
                    """or Annovar keyword (e.g. 'refGene', 'cosmic70', 'clinvar_202*').\n"""
                    """Note that refGene will at least be downloaded,\n"""
                    """and only files that not already exist or changed will be downloaded.\n""",
            "required": False,
            "default": None,
            "type": str
        },
        "download-annovar-url": {
            "metavar": "Annovar url",
            "help": """Annovar databases URL (see Annovar Doc).\n""",
            "required": False,
            "default": DEFAULT_ANNOVAR_URL,
            "type": str
        },

        # snpEff
        "download-snpeff": {
            "metavar": "snpEff",
            "help": """Download snpEff databases within snpEff folder""",
            "required": False,
            "default": None,
            "type": PathType(exists=None, type='dir'),
            "gooey": {
                "widget": "DirChooser",
                "options": {
                    "default_dir": DEFAULT_DATABASE_FOLDER,
                    "message": "Path to snpEff databases folder"
                }
            }
        },

        # refSeq
        "download-refseq": {
            "metavar": "refSeq",
            "help": """Path to refSeq databases\n"""
                    f"""(e.g. '{DEFAULT_REFSEQ_FOLDER}').\n""",
            "required": False,
            "default": None,
            "type": PathType(exists=None, type='dir'),
            "gooey": {
                "widget": "DirChooser",
                "options": {
                    "default_dir": DEFAULT_DATABASE_FOLDER,
                    "message": "Path to refGene files folder"
                }
            }
        },
        "download-refseq-url": {
            "metavar": "refSeq url",
            "help": """refSeq databases URL (see refSeq WebSite)\n"""
                    f"""(e.g. '{DEFAULT_REFSEQ_URL}')•/n""",
            "required": False,
            "default": DEFAULT_REFSEQ_URL,
            "type": str
        },
        "download-refseq-prefix": {
            "metavar": "refSeq prefix",
            "help": """Check existing refSeq files in refSeq folder.\n""",
            "required": False,
            "default": "ncbiRefSeq",
            "type": str
        },
        "download-refseq-files": {
            "metavar": "refSeq files",
            "help": """List of refSeq files to download.\n""",
            "required": False,
            "default": "ncbiRefSeq.txt,ncbiRefSeqLink.txt",
            "type": str
        },
        "download-refseq-format-file": {
            "metavar": "refSeq format file",
            "help": """Name of refSeq file to convert in BED format\n"""
                    """(e.g. 'ncbiRefSeq.txt').\n"""
                    """Process only if not None.\n""",
            "required": False,
            "default": None,
            "type": str
        },
        "download-refseq-include-utr5": {
            "help": """Formating BED refSeq file including 5'UTR.\n""",
            "action": "store_true",
            "default": False
        },
        "download-refseq-include-utr3": {
            "help": """Formating BED refSeq file including 3'UTR.\n""",
            "action": "store_true",
            "default": False
        },
        "download-refseq-include-chrM": {
            "help": """Formating BED refSeq file including Mitochondiral chromosome 'chrM' or 'chrMT'.\n""",
            "action": "store_true",
            "default": False
        },
        "download-refseq-include-non-canonical-chr": {
            "help": """Formating BED refSeq file including non canonical chromosomes.\n""",
            "action": "store_true",
            "default": False
        },
        "download-refseq-include-non-coding-transcripts": {
            "help": """Formating BED refSeq file including non coding transcripts.\n""",
            "action": "store_true",
            "default": False
        },
        "download-refseq-include-transcript-version": {
            "help": """Formating BED refSeq file including transcript version.\n""",
            "action": "store_true",
            "default": False
        },
        
        # dbNSFP
        "download-dbnsfp": {
            "metavar": "dbNSFP",
            "help": """Download dbNSFP databases within dbNSFP folder"""
                    f"""(e.g. '{DEFAULT_DATABASE_FOLDER}').\n""",
            "required": False,
            "default": None,
            "type": PathType(exists=None, type='dir'),
            "gooey": {
                "widget": "DirChooser",
                "options": {
                    "default_dir": DEFAULT_DATABASE_FOLDER,
                    "message": "Path to dbNSFP databases folder"
                }
            }
        },
        "download-dbnsfp-url": {
            "metavar": "dbNSFP url",
            "help": """Download dbNSFP databases URL (see dbNSFP website)\n"""
                    f"""(e.g. {DEFAULT_DBNSFP_URL}').\n""",
            "required": False,
            "default": DEFAULT_DBNSFP_URL,
            "type": str
        },
        "download-dbnsfp-release": {
            "metavar": "dnNSFP release",
            "help": """Release of dbNSFP to download (see dbNSFP website)\n"""
                    """(e.g. '4.4a').\n""",
            "required": False,
            "default": "4.4a"
        },
        "download-dbnsfp-parquet-size": {
            "metavar": "dbNSFP parquet size",
            "help": """Maximum size (Mb) of data files in Parquet folder.\n"""
                    """Parquet folder are partitioned (hive) by chromosome (sub-folder),\n"""
                    """which contain N data files.\n""",
            "required": False,
            "default": 100,
            "type": int,
            "gooey": {
                "widget": "IntegerField",
                "options": {
                    'min': 1,
                    'max': 100000,
                    'increment': 10
                }
            }
        },
        "download-dbnsfp-subdatabases": {
            "help": """Generate dbNSFP sub-databases.\n"""
                    """dbNSFP provides multiple databases which are split onto multiple columns.\n"""
                    """This option create a Parquet folder for each sub-database (based on columns names).\n""",
            "action": "store_true",
            "default": False
        },
        "download-dbnsfp-parquet": {
            "help": """Generate a Parquet file for each Parquet folder.\n""",
            "action": "store_true",
            "default": False
        },
        "download-dbnsfp-vcf": {
            "help": """Generate a VCF file for each Parquet folder.\n"""
                    """Need genome FASTA file (see --download-genome).\n""",
            "action": "store_true",
            "default": False
        },
        "download-dbnsfp-no-files-all": {
            "help": """Not generate database Parquet/VCF file for the entire database ('ALL').\n"""
                    """Only sub-databases files will be generated.\n"""
                    """(see '--download-dbnsfp-subdatabases').\n""",
            "action": "store_true",
            "default": False
        },
        "download-dbnsfp-add-info": {
            "help": """Add INFO column (VCF format) in Parquet folder and file.\n"""
                    """Useful for speed up full annotation (all available columns).\n"""
                    """Increase memory and space during generation of files.\n""",
            "action": "store_true",
            "default": False
        },
        "download-dbnsfp-row-group-size": {
            "metavar": "dnNSFP row grooup size",
            "help": """Minimum number of rows in a parquet row group (see duckDB doc).\n"""
                    """Lower can reduce memory usage and slightly increase space during generation,\n"""
                    """speed up highly selective queries, slow down whole file queries (e.g. aggregations).\n""",
            "required": False,
            "default": 100000,
            "type": int,
            "gooey": {
                "widget": "IntegerField",
                "options": {
                    'min': 1,
                    'max': 100000000000,
                    'increment': 10000
                }
            }
        },

        # AlphaMissense
        "download-alphamissense": {
            "metavar": "AlphaMissense",
            "help": "Path to AlphaMissense databases",
            "required": False,
            "default": None,
            "type": PathType(exists=None, type='dir'),
            "gooey": {
                "widget": "DirChooser",
                "options": {
                    "default_dir": DEFAULT_DATABASE_FOLDER,
                    "message": "Path to Alphamissense databases folder"
                }
            }
        },
        "download-alphamissense-url": {
            "metavar": "AlphaMissense url",
            "help": """Download AlphaMissense databases URL (see AlphaMissense website)\n"""
                    f"""(e.g. '{DEFAULT_ALPHAMISSENSE_URL}').\n""",
            "required": False,
            "default": DEFAULT_ALPHAMISSENSE_URL,
            "type": str
        },

        # Exomiser
        "download-exomiser": {
            "metavar": "Exomiser",
            "help": """Path to Exomiser databases\n"""
                    f"""(e.g. {DEFAULT_EXOMISER_FOLDER}).\n""",
            "required": False,
            "default": None,
            "type": PathType(exists=None, type='dir'),
            "gooey": {
                "widget": "DirChooser",
                "options": {
                    "default_dir": DEFAULT_DATABASE_FOLDER,
                    "message": "Path to Exomiser databases folder"
                }
            }
        },
        "download-exomiser-application-properties": {
            "metavar": "Exomiser application properties",
            "help": """Exomiser Application Properties configuration file (see Exomiser website).\n"""
                    """This file contains configuration settings for the Exomiser tool.\n"""
                    """If this parameter is not provided, the function will attempt to locate\n"""
                    """the application properties file automatically based on the Exomiser.\n"""
                    """Configuration information will be used to download expected releases (if no other parameters).\n"""
                    """CADD and REMM will be downloaded only if 'path' are provided.\n""",
            "required": False,
            "default": None,
            "type": PathType(exists=True, type='file'),
            "gooey": {
                "widget": "FileChooser",
                "options": {
                    "wildcard":
                        "All files (*)|*",
                    "options": {
                    "default_dir": DEFAULT_EXOMISER_FOLDER,
                    "message": "Path to Exomiser application properties file"
                }
                }
            }
        },
        "download-exomiser-url": {
            "metavar": "Exomiser url",
            "help": """URL where Exomiser database files can be downloaded from\n"""
                    f"""(e.g. '{DEFAULT_EXOMISER_URL}').\n""",
            "required": False,
            "default": DEFAULT_EXOMISER_URL,
            "type": str
        },
        "download-exomiser-release": {
            "metavar": "Exomiser release",
            "help": """Release of Exomiser data to download.\n"""
                    """If "default", "auto", or "config", retrieve from Application Properties file.\n"""
                    """If not provided (None), from Application Properties file (Exomiser data-version) \n"""
                    """or default '2109'.\n""",
            "required": False,
            "default": None,
            "type": str
        },
        "download-exomiser-phenotype-release": {
            "metavar": "Exomiser phenoptye release",
            "help": """Release of Exomiser phenotype to download.\n"""
                    """If not provided (None), from Application Properties file (Exomiser Phenotype data-version)\n"""
                    """or Exomiser release.\n""",
            "required": False,
            "default": None,
            "type": str
        },
        "download-exomiser-remm-release": {
            "metavar": "Exomiser remm release",
            "help": """Release of ReMM (Regulatory Mendelian Mutation) database to download.\n"""
                    """If "default", "auto", or "config", retrieve from Application Properties file.\n""",
            "required": False,
            "default": None,
            "type": str
        },
        "download-exomiser-remm-url": {
            "metavar": "Exomiser remm url",
            "help": """URL where ReMM (Regulatory Mendelian Mutation) database files can be downloaded from\n"""
                    f"""(e.g. '{DEFAULT_EXOMISER_REMM_URL}').\n""",
            "required": False,
            "default": DEFAULT_EXOMISER_REMM_URL,
            "type": str
        },
        "download-exomiser-cadd-release": {
            "metavar": "Exomiser cadd release",
            "help": """Release of CADD (Combined Annotation Dependent Depletion) database to download.\n"""
                    """If "default", "auto", or "config", retrieve from Application Properties file.\n""",
            "required": False,
            "default": None,
            "type": str
        },
        "download-exomiser-cadd-url": {
            "metavar": "Exomiser cadd url",
            "help": """URL where CADD (Combined Annotation Dependent Depletion) database files can be downloaded from\n"""
                    f"""(e.g. '{DEFAULT_EXOMISER_CADD_URL}').\n""",
            "required": False,
            "default": DEFAULT_EXOMISER_CADD_URL,
            "type": str
        },
        "download-exomiser-cadd-url-snv-file": {
            "metavar": "Exomiser url snv file",
            "help": """Name of the file containing the SNV (Single Nucleotide Variant) data\n"""
                    """for the CADD (Combined Annotation Dependent Depletion) database.\n""",
            "required": False,
            "default": "whole_genome_SNVs.tsv.gz",
            "type": str
        },
        "download-exomiser-cadd-url-indel-file": {
            "metavar": "Exomiser cadd url indel",
            "help": """Name of the file containing the INDEL (Insertion-Deletion) data\n"""
                    """for the CADD (Combined Annotation Dependent Depletion) database.\n""",
            "required": False,
            "default": "InDels.tsv.gz",
            "type": str
        },

        # dbSNP
        "download-dbsnp": {
            "metavar": "dnSNP",
            "help": """Path to dbSNP databases\n"""
                    f"""(e.g. '{DEFAULT_DBSNP_FOLDER}').\n""",
            "required": False,
            "default": None,
            "type": PathType(exists=None, type='dir'),
            "gooey": {
                "widget": "DirChooser",
                "options": {
                    "default_dir": DEFAULT_DATABASE_FOLDER,
                    "message": "Path to dbSNP databases folder"
                }
            }
        },
        "download-dbsnp-releases": {
            "metavar": "dnSNP releases",
            "help": """Release of dbSNP to download\n"""
                    """(e.g. 'b152', 'b152,b156').\n""",
            "required": False,
            "default": 'b156',
            "type": str
        },
        "download-dbsnp-release-default": {
            "metavar": "dnSNP release default",
            "help": """Default Release of dbSNP ('default' symlink)\n"""
                    """(e.g. 'b156').\n"""
                    """If None, first release to download will be assigned as default\n"""
                    """only if it does not exists.\n""",
            "required": False,
            "default": None,
            "type": str
        },
        "download-dbsnp-url": {
            "metavar": "dbSNP url",
            "help": """URL where dbSNP database files can be downloaded from.\n"""
                    f"""(e.g. '{DEFAULT_DBSNP_URL}').\n""",
            "required": False,
            "default": DEFAULT_DBSNP_URL,
            "type": str
        },
        "download-dbsnp-url-files": {
            "metavar": "dbSNP url files",
            "help": """Dictionary that maps assembly names to specific dbSNP URL files.\n"""
                    """It allows you to provide custom dbSNP URL files for specific assemblies\n"""
                    """instead of using the default file naming convention.\n""",
            "required": False,
            "default": None,
            "type": str
        },
        "download-dbsnp-url-files-prefix": {
            "metavar": "dbSNP url files prefix",
            "help": """String that represents the prefix of the dbSNP file name for a specific assembly.\n"""
                    """It is used to construct the full URL of the dbSNP file to be downloaded.\n""",
            "required": False,
            "default": "GCF_000001405",
            "type": str
        },
        "download-dbsnp-assemblies-map": {
            "metavar": "dbSNP assemblies map",
            "help": """dictionary that maps assembly names to their corresponding dbSNP versions.\n"""
                    """It is used to construct the dbSNP file name based on the assembly name.\n""",
            "required": False,
            "default": {"hg19": "25", "hg38": "40"},
            "type": str,
            "gooey": {
                "options": {
                    'initial_value': '{"hg19": "25", "hg38": "40"}'  
                }
            }
        },
        "download-dbsnp-vcf": {
            "help": """Generate well-formatted VCF from downloaded file:\n"""
                    """- Add and filter contigs associated to assembly\n"""
                    """- Normalize by splitting multiallelics\n"""
                    """- Need genome (see --download-genome)\n""",
            "action": "store_true",
            "default": False
        },
        "download-dbsnp-parquet": {
            "help": """Generate Parquet file from VCF.\n""",
            "action": "store_true",
            "default": False
        },

        # HGMD
        "convert-hgmd": {
            "metavar": "HGMD",
            "help": """Convert HGMD databases.\n"""
                    """Folder where the HGMD databases will be stored.\n"""
                    """Fields in VCF, Parquet and TSV will be generated.\n"""
                    """If the folder does not exist, it will be created.\n""",
            "required": False,
            "default": None,
            "type": PathType(exists=None, type='dir'),
            "gooey": {
                "widget": "DirChooser"
            }
        },
        "convert-hgmd-file": {
            "metavar": "HGMD file",
            "help": """File from HGMD.\n"""
                    """Name format 'HGMD_Pro_<release>_<assembly>.vcf.gz'.\n""",
            "required": False,
            "default": None,
            "type": PathType(exists=True, type='file'),
            "gooey": {
                "widget": "FileChooser"
            }
        },
        "convert-hgmd-basename": {
            "metavar": "HGMD basename",
            "help": """File output basename.\n"""
                    """Generated files will be prefixed by basename\n"""
                    """(e.g. 'HGMD_Pro_MY_RELEASE')\n"""
                    """By default (None), input file name without '.vcf.gz'.\n""",
            "required": False,
            "default": None,
            "type": str
        },

        # Databases parameters
        "generate-param": {
            "metavar": "param",
            "help": """Parameter file (JSON) with all databases found.\n"""
                    """Databases folders scanned are defined in config file.\n"""
                    """Structure of databases follow this structure (see doc):\n"""
                    """.../<database>/<release>/<assembly>/*.[parquet|vcf.gz|...]\n""",
            "required": False,
            "default": None,
            "type": PathType(exists=None, type=None),
            "gooey": {
                "widget": "FileSaver",
                "options": {
                    "wildcard":
                        "JSON file (*.json)|*.json"
                }
            }
        },
        "generate-param-description": {
            "metavar": "param description",
            "help": """Description file (JSON) with all databases found.\n"""
                    """Contains all databases with description of format, assembly, fields...\n""",
            "required": False,
            "default": None,
            "type": PathType(exists=None, type=None),
            "gooey": {
                "widget": "FileSaver",
                "options": {
                    "wildcard":
                        "JSON file (*.json)|*.json"
                }
            }
        },
        "generate-param-releases": {
            "metavar": "param release",
            "help": """List of database folder releases to check\n"""
                    """(e.g. 'current', 'latest').\n""",
            "required": False,
            "default": "current",
            "type": str
        },
        "generate-param-formats": {
            "metavar": "param formats",
            "help": """List of database formats to check\n"""
                    """(e.g. 'parquet', 'parquet,vcf,bed,tsv').\n""",
            "required": False,
            "default": "parquet",
            "type": str
        },
        "generate-param-bcftools": {
            "help": """Generate parameter JSON file with BCFTools annotation for allowed formats\n"""
                    """(i.e. 'vcf', 'bed').\n""",
            "action": "store_true",
            "default": False
        },

        # From Annovar
        "annovar-code": {
            "metavar": "Annovar code",
            "help": """Annovar code, or database name.\n"""
                    """Usefull to name databases columns.\n""",
            "required": False,
            "default": None,
            "type": str
        },
        "to_parquet": {
            "metavar": "to parquet",
            "help": """Parquet file conversion.\n""",
            "required": False,
            "default": None,
            "type": PathType(exists=None, type=None),
            "gooey": {
                "widget": "FileSaver",
                "options": {
                    "wildcard":
                        "HTML file (*.parquet)|*.parquet",
                }
            }
        },

        # Help
        "help_md": {
            "metavar": "help markdown",
            "help": """Help Output file in MarkDown format.\n""",
            "required": False,
            "default": None,
            "type": PathType(exists=None, type=None),
            "gooey": {
                "widget": "FileSaver",
                "options": {
                    "wildcard":
                        "HTML file (*.md)|*.md",
                }
            }
        },
        "help_html": {
            "metavar": "help html",
            "help": """Help Output file in HTML format.\n""",
            "required": False,
            "default": None,
            "type": PathType(exists=None, type=None),
            "gooey": {
                "widget": "FileSaver",
                "options": {
                    "wildcard":
                        "HTML file (*.html)|*.html",
                }
            }
        },
        "help_json_input": {
            "metavar": "help JSON input",
            "help": """Help input file in JSON format.\n""",
            "required": False,
            "default": None,
            "type": PathType(exists=True, type='file'),
            "gooey": {
                "widget": "FileChooser",
                "options": {
                    "wildcard":
                        "JSON file (*.json)|*.json|"
                        "All files (*)|*",
                }
            }
        },
        "code_type": {
            "metavar": "example code type",
            "help": """Help example code type for input JSON format\n"""
                    """(e.g. 'json', 'bash').\n""",
            "required": False,
            "default": "",
            "type": str
        },
        "help_json_input_title": {
            "metavar": "help JSON input title",
            "help": """Help JSON input title.\n""",
            "required": False,
            "default": "Help",
            "type": str
        },

        # Common
        "genomes-folder": {
            "metavar": "genomes",
            "help": """Folder containing genomes.\n"""
                    f"""(e.g. '{DEFAULT_GENOME_FOLDER}'""",
            "required": False,
            "default": DEFAULT_GENOME_FOLDER,
            "type": PathType(exists=True, type='dir'),
            "gooey": {
                "widget": "DirChooser",
                "options": {
                    "default_dir": DEFAULT_GENOME_FOLDER,
                    "message": "Path to genomes databases folder"
                }
            }
        },

        # Shared
        "config": {
            "metavar": "config",
            "help": """Configuration JSON file or JSON string.\n""",
            "required": False,
            "default": "{}",
            "type": str,
            "gooey": {
                "widget": "FileChooser",
                "options": {
                    'initial_value': '{}'  
                }
            }
        },
        "threads": {
            "metavar": "threads",
            "help": """Specify the number of threads to use for processing HOWARD.\n"""
                    """It determines the level of parallelism,\n"""
                    """either on python scripts, duckdb engine and external tools.\n"""
                    """It and can help speed up the process/tool.\n"""
                    """Use -1 to use all available CPU/cores.\n"""
                    """Either non valid value is 1 CPU/core.\n""",
            "required": False,
            "type": int,
            "default": -1,
            "gooey": {
                "widget": "IntegerField",
                "options": {
                    'min': -1,
                    'max': 1000,
                    'increment': 1
                }
            },
            "extra": {
                "examples": {
                    "# Automatically detect all available CPU/cores":
                        "\"threads\": -1",
                "# Define 8 CPU/cores":
                  "\"threads\": 8"
                }
            }
        },
        "memory": {
            "metavar": "memory",
            "help": """Specify the memory to use in format FLOAT[kMG]\n"""
                    """(e.g. '8G', '12.42G', '1024M').\n"""
                    """It determines the amount of memory for duckDB engine and external tools\n"""
                    """(especially for JAR programs).\n"""
                    """It can help to prevent 'out of memory' failures.\n"""
                    """By default (None) is 80%% of RAM (for duckDB).\n""",
            "required": False,
            "type": str,
            "default": None,
            "extra": {
                "format": "FLOAT[kMG]",
                "examples": {
                    "# Automatically detect all available CPU/cores":
                        "\"threads\": -1",
                "# Define 8 CPU/cores":
                  "\"threads\": 8"
                }
            }
            
        },
        "chunk_size": {
            "metavar": "chunk size",
            "help": """Number of records in batch to export output file.\n"""
                    """The lower the chunk size, the less memory consumption.\n"""
                    """For Parquet partitioning, files size will depend on the chunk size.\n""",
            "required": False,
            "default": 1000000,
            "type": int,
            "gooey": {
                "widget": "IntegerField",
                "options": {
                    'min': 1,
                    'max': 100000000000,
                    'increment': 10000
                }
            },
            "extra": {
                "examples": {
                    "Chunk size of 1.000.000 by default":
                        "\"chunk_size\": 1000000",
                    "Smaller chunk size to reduce Parquet file size and memory usage":
                        "\"chunk_size\": 100000"
                }
            }
        },
        "tmp": {
            "metavar": "Temporary folder",
            "help": """Temporary folder (e.g. '/tmp').\n"""
                    """By default, '.tmp' for duckDB (see doc),"""
                    """external tools and python scripts.\n""",
            "required": False,
            "default": None,
            "type": PathType(exists=True, type='dir'),
            "gooey": {
                "widget": "DirChooser"
            },
            "extra": {
                "examples": {
                    "# System temporary folder":
                        "\"tmp\": \"/tmp\"",
                    "# HOWARD work directory":
                        "\"tmp\": \"~/howard/tmp\"",
                    "# Current work directory":
                        "\"tmp\": \".tmp\""
                }
            }
        },
        "duckdb_settings": {
            "metavar": "duckDB settings",
            "help": """DuckDB settings (see duckDB doc) as JSON (string or file).\n"""
                    """These settings have priority (see options 'threads', 'tmp'...).\n"""
                    """Examples: '{"TimeZone": "GMT", "temp_directory": "/tmp/duckdb", "threads": 8}'.\n""",
            "required": False,
            "default": None,
            "type": PathType(exists=True, type='file'),
            "gooey": {
                "widget": "FileChooser",
                "options": {
                    "wildcard":
                        "JSON file (*.json)|*.json|"
                        "All files (*)|*",
                }
            },
            "extra": {
                "examples": {
                    "DuckDB settings JSON file":
                        "\"duckdb_settings\": \"/path/to/duckdb_config.json\"",
                    "JSON string for Time zone, temporary directory and threads for duckDB":
                        """\"duckdb_settings\": {\n"""
                        """   \"TimeZone\": \"GMT\",\n"""
                        """   \"temp_directory\": \"/tmp/duckdb\",\n"""
                        """   \"threads\": 8\n"""
                        """}"""
                }
            }
        },
        "verbosity": {
            "metavar": "verbosity",
            "help": """Verbosity level\n"""
                    """Available: CRITICAL, ERROR, WARNING, INFO, DEBUG or NOTSET\n"""
                    """- DEBUG: Detailed information, typically of interest only when diagnosing problems.\n"""
                    """- INFO: Confirmation that things are working as expected.\n"""
                    """- WARNING: An indication that something unexpected happened.\n"""
                    """- ERROR: Due to a more serious problem.\n"""
                    """- CRITICAL: A serious error.\n"""
                    """- NOTSET: All messages.\n""",
            "required": False,
            "choices": ['CRITICAL', 'ERROR', 'WARNING', 'INFO', 'DEBUG', 'NOTSET'],
            "default": "INFO",
            "type": str,
            "gooey": {
                "widget": "Dropdown",
                "options": {}
            },
            "extra": {
                "examples": {
                    "Default verbosity":
                        "\"verbosity\": \"INFO\"",
                    "ERROR level (quiet mode)":
                        "\"verbosity\": \"ERROR\"",
                    "For debug":
                        "\"verbosity\": \"DEBUG\""
                }
            }
        },
        "access": {
            "metavar": "access mode",
            "help": """Access mode to variants file or database.\n"""
                    """Either 'RW' for Read and Write, or 'RO' for Read Only.\n""",
            "default": 'RW',
            "type": str,
            "choices": ["RW", "RO"],
            "gooey": {
                "widget": "Dropdown",
                "options": {}
            },
            "extra": {
                "examples": {
                    "Read and Write mode":
                        "\"access\": \"RW\"",
                    "Read only mode":
                        "\"access\": \"RO\""
                }
            }
        },
        "log": {
            "metavar": "log",
            "help": """Logs file\n"""
                    """(e.g. 'my.log').\n""",
            "required": False,
            "default": None,
            "type": PathType(exists=None, type='file'),
            "gooey": {
                "widget": "FileSaver"
            },
            "extra": {
                "examples": {
                    "Relative path to log file":
                        "\"log\": \"my.log\"",
                    "# HOWARD work directory":
                        "\"log\": \"~/howard/log\"",
                    "Full path to log file":
                        "\"log\": \"/tmp/my.log\""
                }
            }
        },
        "quiet": {
            "help": argparse.SUPPRESS,
            "action": "store_true",
            "default": False
        },
        "verbose": {
            "help": argparse.SUPPRESS,
            "action": "store_true",
            "default": False
        },
        "debug": {
            "help": argparse.SUPPRESS,
            "action": "store_true",
            "default": False
        },

        # Only for HELP

    }


# Shared arguments
shared_arguments = ["config", "threads", "memory", "chunk_size", "tmp", "duckdb_settings", "verbosity", "log", "quiet", "verbose", "debug"]

# Command dict
commands_arguments = {
    "query": {
        "function" : "query",
        "description":  """Query genetic variations in SQL format. Data can be loaded into 'variants' table from various formats (e.g. VCF, TSV, Parquet...). Using --explode_infos allow query on INFO/tag annotations. SQL query can also use external data within the request, such as a Parquet file(s).  """,
        "help": "Query genetic variations file in SQL format.",
        "epilog": """Usage examples:\n"""
                        """   howard query --input=tests/data/example.vcf.gz --query="SELECT * FROM variants WHERE REF = 'A' AND POS < 100000" \n"""
                        """   howard query --input=tests/data/example.vcf.gz --explode_infos --query='SELECT "#CHROM", POS, REF, ALT, DP, CLNSIG, sample2, sample3 FROM variants WHERE DP >= 50 OR CLNSIG NOT NULL ORDER BY DP DESC' \n"""
                        """   howard query --query="SELECT \\\"#CHROM\\\", POS, REF, ALT, \\\"INFO/Interpro_domain\\\" FROM 'tests/databases/annotations/current/hg19/dbnsfp42a.parquet' WHERE \\\"INFO/Interpro_domain\\\" NOT NULL ORDER BY \\\"INFO/SiPhy_29way_logOdds_rankscore\\\" DESC LIMIT 10" \n"""
                        """   howard query --explode_infos --explode_infos_prefix='INFO/' --query="SELECT \\\"#CHROM\\\", POS, REF, ALT, STRING_AGG(INFO, ';') AS INFO FROM 'tests/databases/annotations/current/hg19/*.parquet' GROUP BY \\\"#CHROM\\\", POS, REF, ALT" --output=/tmp/full_annotation.tsv  && head -n2 /tmp/full_annotation.tsv \n"""
                        """   howard query --input=tests/data/example.vcf.gz --param=config/param.json"""
                        , 
        "groups": {
            "main": {
                "input": False,
                "output": False,
                "param": False,
                "query": False
            },
            "Explode": {
                "explode_infos": False,
                "explode_infos_prefix": False,
                "explode_infos_fields": False,
            },
            "Query": {
                "query_limit": False,
                "query_print_mode": False
            },
            "Export": {
                "include_header": False,
                "parquet_partition": False
            }
        }
    },
    "stats": {
        "function" : "stats",
        "description":  """Statistics on genetic variations, such as: number of variants, number of samples, statistics by chromosome, genotypes by samples...""",
        "help": "Statistics on genetic variations file.",
        "epilog": """Usage examples:\n"""
                        """   howard stats --input=tests/data/example.vcf.gz \n"""
                        """   howard stats --input=tests/data/example.vcf.gz --stats_md=/tmp/stats.md \n"""
                        """   howard stats --input=tests/data/example.vcf.gz --param=config/param.json """,
        "groups": {
            "main": {
                "input": True,
                "param": False
            },
            "Stats": {
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
                        """   howard convert --input=tests/data/example.vcf.gz --output=/tmp/example.tsv --explode_infos --explode_infos_prefix='INFO/' --explode_infos_fields='CLNSIG,SIFT,DP,*' --order_by='"INFO/CLNSIG" DESC, "INFO/DP" DESC' --include_header \n"""
                        """   howard convert --input=tests/data/example.vcf.gz --output=/tmp/example.tsv --param=config/param.json """,
        "groups": {
            "main": {
                "input": True,
                "output": True,
                "param": False
            },
            "Explode": {
                "explode_infos": False,
                "explode_infos_prefix": False,
                "explode_infos_fields": False,
            },
            "Export": {
                "include_header": False,
                "order_by": False,
                "parquet_partition": False
            }
        }
    },
    "hgvs": {
        "function" : "hgvs",
        "description":  """HGVS annotation using HUGO HGVS internation Sequence Variant Nomenclature (http://varnomen.hgvs.org/). Annotation refere to refGene and genome to generate HGVS nomenclature for all available transcripts. This annotation add 'hgvs' field into VCF INFO column of a VCF file.""",
        "help":         """HGVS annotation (HUGO internation nomenclature) using refGene, genome and transcripts list.\n""",
        "epilog":       """Usage examples:\n"""
                        """   howard hgvs --input=tests/data/example.full.vcf --output=/tmp/example.hgvs.vcf \n"""
                        """   howard hgvs --input=tests/data/example.full.vcf --output=/tmp/example.hgvs.tsv --param=config/param.json \n"""
                        """   howard hgvs --input=tests/data/example.full.vcf --output=/tmp/example.hgvs.vcf --full_format --use_exon \n"""
                        """    \n"""
                        ,
        "groups": {
            "main": {
                "input": True,
                "output": True,
                "param": False,
                "hgvs_options": False,
                "assembly": False
            },
            "HGVS": {
                "use_gene": False,
                "use_exon": False,
                "use_protein": False,
                "add_protein": False,
                "full_format": False,
                "codon_type": False,
                "refgene": False,
                "refseqlink": False
            },
            "Databases": {
                "refseq-folder": False,
                "genomes-folder": False
            }
        }
    },
    "annotation": {
        "function" : "annotation",
        "description":  """Annotation is mainly based on a build-in Parquet annotation method, and tools such as BCFTOOLS, Annovar and snpEff. It uses available databases (see Annovar and snpEff) and homemade databases. Format of databases are: parquet, duckdb, vcf, bed, Annovar and snpEff (Annovar and snpEff databases are automatically downloaded, see howard databases tool). """,
        "help":         """Annotation of genetic variations file using databases/files and tools.""",
        "epilog":       """Usage examples:\n"""
                        """   howard annotation --input=tests/data/example.vcf.gz --output=/tmp/example.howard.vcf.gz --annotations='tests/databases/annotations/current/hg19/avsnp150.parquet,tests/databases/annotations/current/hg19/dbnsfp42a.parquet,tests/databases/annotations/current/hg19/gnomad211_genome.parquet' \n"""
                        """   howard annotation --input=tests/data/example.vcf.gz --output=/tmp/example.howard.tsv --assembly=hg19 --annotations='annovar:refGene,annovar:cosmic70,snpeff,tests/databases/annotations/current/hg19/clinvar_20210123.parquet' \n"""
                        """   howard annotation --input=tests/data/example.vcf.gz --output=/tmp/example.howard.tsv --assembly=hg19 --annotations='ALL:parquet' \n"""
                        """   howard annotation --input=tests/data/example.vcf.gz --output=/tmp/example.howard.tsv --param=config/param.json """, 
        "groups": {
            "main": {
                "input": True,
                "output": True,
                "param": False,
                "annotations": False,
                "assembly": False,
            },
            "annotation": {
                "annotations_update": False,
                "annotations_append": False,
                
            }
        }
    },
    "calculation": {
        "function" : "calculation",
        "description":  """Calculation processes variants information to generate new information, such as: identify variation type (VarType), harmonizes allele frequency (VAF) and calculate sttistics (VAF_stats), extracts Nomen (transcript, cNomen, pNomen...) from an HGVS field (e.g. snpEff, Annovar) with an optional list of personalized transcripts, generates VaRank format barcode, identify trio inheritance.""",
        "help":         """Calculation operations on genetic variations file and genotype information.\n""",
        "epilog":       """Usage examples:\n"""
                        """   howard calculation --input=tests/data/example.full.vcf --output=/tmp/example.calculation.tsv --calculations='vartype' \n"""
                        """   howard calculation --input=tests/data/example.ann.vcf.gz --output=/tmp/example.calculated.tsv --calculations='snpeff_hgvs,NOMEN' --hgvs_field=snpeff_hgvs --transcripts=tests/data/transcripts.tsv \n"""
                        """   howard calculation --input=tests/data/example.ann.vcf.gz --output=/tmp/example.ann.tsv --param=config/param.json \n"""
                        """   howard calculation --show_calculations \n"""
                        ,
        "groups": {
            "main": {
                "input": False,
                "output": False,
                "param": False,
                "calculations": False,
                "calculation_config": False,
                "show_calculations": False
            },
            "NOMEN": {
                "hgvs_field": False,
                "transcripts": False
            },
            "TRIO": {
                "trio_pedigree": False
            }
        }
    },
    "prioritization": {
        "function" : "prioritization",
        "description":  """Prioritization algorithm uses profiles to flag variants (as passed or filtered), calculate a prioritization score, and automatically generate a comment for each variants (example: 'polymorphism identified in dbSNP. associated to Lung Cancer. Found in ClinVar database'). Prioritization profiles are defined in a configuration file in JSON format. A profile is defined as a list of annotation/value, using wildcards and comparison options (contains, lower than, greater than, equal...). Annotations fields may be quality values (usually from callers, such as 'DP') or other annotations fields provided by annotations tools, such as HOWARD itself (example: COSMIC, Clinvar, 1000genomes, PolyPhen, SIFT). Multiple profiles can be used simultaneously, which is useful to define multiple validation/prioritization levels (example: 'standard', 'stringent', 'rare variants', 'low allele frequency').\n""",
        "help": "Prioritization of genetic variations based on annotations criteria (profiles).",
        "epilog": """Usage examples:\n"""
                        """   howard prioritization --input=tests/data/example.vcf.gz --output=/tmp/example.prioritized.vcf.gz --prioritizations=config/prioritization_profiles.json --profiles='default,GERMLINE' \n"""
                        """   howard prioritization --input=tests/data/example.vcf.gz --output=/tmp/example.prioritized.tsv --param=config/param.json """, 
        "groups": {
            "main": {
                "input": True,
                "output": True,
                "prioritizations": False,
                "param": False
            },
            "Prioritization": {
                "profiles": False,
                "default_profile": False,
                "pzfields": False,
                "prioritization_score_mode": False
            }
        }
    },
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
                        """   howard process --input=tests/data/example.vcf.gz --annotations='snpeff' --calculations='snpeff_hgvs' --prioritizations=config/prioritization_profiles.json --explode_infos --output=/tmp/example.annotated.tsv --query='SELECT "#CHROM", POS, ALT, REF, snpeff_hgvs FROM variants' \n"""
                        """   howard process --input=tests/data/example.vcf.gz --hgvs='full_format,use_exon' --explode_infos --output=/tmp/example.annotated.tsv --query='SELECT "#CHROM", POS, ALT, REF, hgvs FROM variants' \n"""
                        """   howard process --input=tests/data/example.vcf.gz --output=/tmp/example.howard.vcf.gz --hgvs='full_format,use_exon' --annotations='tests/databases/annotations/current/hg19/avsnp150.parquet,tests/databases/annotations/current/hg19/dbnsfp42a.parquet,tests/databases/annotations/current/hg19/gnomad211_genome.parquet' --calculation='NOMEN' --explode_infos --query='SELECT NOMEN, REVEL_score, SIFT_score, AF AS 'gnomad_AF', ClinPred_score, ClinPred_pred FROM variants' \n"""
                        , 

        "groups": {
            "main": {
                "input": True,
                "output": True,
                "param": False
            },
            "Quick Processes": {
                "hgvs_options": False,
                "annotations": False,
                "calculations": False,
                "prioritizations": False
            },
            "Query": {
                "query": False,
                "query_limit": False,
                "query_print_mode": False,
            },
            "Export": {
                "include_header": False,
                "order_by": False,
                "parquet_partitions": False
            },
            "Explode": {
                "explode_infos": False,
                "explode_infos_prefix": False,
                "explode_infos_fields": False
            }
        }
    },
    "databases": {
        "function" : "databases",
        "description": """Download databases and needed files for howard and associated tools""",
        "help": """Download databases and needed files for howard and associated tools""",
        "epilog": """Usage examples:\n"""
                    """   howard databases --assembly=hg19 --download-genomes=~/howard/databases/genomes/current --download-genomes-provider=UCSC --download-genomes-contig-regex='chr[0-9XYM]+$' \n"""
                    """   howard databases --assembly=hg19 --download-annovar=~/howard/databases/annovar/current --download-annovar-files='refGene,cosmic70,nci60' \n"""
                    """   howard databases --assembly=hg19 --download-snpeff=~/howard/databases/snpeff/current \n"""
                    """   howard databases --assembly=hg19 --download-refseq=~/howard/databases/refseq/current --download-refseq-format-file='ncbiRefSeq.txt' \n"""
                    """   howard databases --assembly=hg19 --download-dbnsfp=~/howard/databases/dbnsfp/current --download-dbnsfp-release='4.4a' --download-dbnsfp-subdatabases \n"""
                    """   howard databases --assembly=hg19 --download-alphamissense=~/howard/databases/alphamissense/current \n"""
                    """   howard databases --assembly=hg19 --download-exomiser=~/howard/databases/exomiser/current \n"""
                    """   howard databases --assembly=hg19 --download-dbsnp=~/howard/databases/dbsnp/current --download-dbsnp-vcf \n"""
                    """   cd ~/howard/databases && howard databases --assembly=hg19 --download-genomes=genomes/current --download-genomes-provider=UCSC --download-genomes-contig-regex='chr[0-9XYM]+$' --download-annovar=annovar/current --download-annovar-files='refGene,cosmic70,nci60' --download-snpeff=snpeff/current --download-refseq=refseq/current --download-refseq-format-file='ncbiRefSeq.txt' --download-dbnsfp=dbnsfp/current --download-dbnsfp-release='4.4a' --download-dbnsfp-subdatabases --download-alphamissense=alphamissense/current --download-exomiser=exomiser/current --download-dbsnp=dbsnp/current --download-dbsnp-vcf --threads=8 \n"""
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
                    """   howard from_annovar --input=tests/databases/others/hg19_nci60.txt --output=/tmp/nci60.from_annovar.vcf.gz --to_parquet=/tmp/nci60.from_annovar.parquet --annovar-code=nci60 --genome=~/howard/databases/genomes/current/hg19.fa --threads=8 """, 
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
    },
    "gui": {
        "function" : "gui",
        "description": """Graphical User Interface tools""",
        "help": """Graphical User Interface tools""",
        "epilog": """Usage examples:\n"""
                    """   howard gui """,
        "groups": {
            
        }
    },
    "help": {
        "function" : "help",
        "description": """Help tools""",
        "help": """Help tools""",
        "epilog": """Usage examples:\n"""
                    """   howard help --help_md=docs/help.md --help_html=docs/help.html\n"""
                    """   howard help --help_json_input=docs/help.config.json --help_json_input_title='HOWARD Configuration' --help_md=docs/help.config.md --help_html=docs/help.config.html --code_type='json'\n"""
                    """   howard help --help_json_input=docs/help.param.json --help_json_input_title='HOWARD Parameters' --help_md=docs/help.param.md --help_html=docs/help.param.html --code_type='json'""",
        "groups": {
            "main": {
                "help_md": False,
                "help_html": False,
                "help_json_input": False,
                "help_json_input_title": False,
                "code_type": False
            }
        }
    }
}

arguments_dict = {
    "arguments": arguments,
    "commands_arguments": commands_arguments,
    "shared_arguments": shared_arguments 
}
