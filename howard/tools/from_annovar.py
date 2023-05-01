#!/usr/bin/env python

import csv
import gc
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
from functools import partial
import itertools

import pyarrow as pa
import pyarrow.parquet as pq
import pyarrow.csv as csv


from howard.objects.variants import Variants
from howard.objects.annotation import Annotation
from howard.commons import *


# Usage
# time howard from_annovar --input=/databases/annovar/current/hg19_nci60.txt --output=/databases/annotations/current/hg19/nci60.vcf.gz --to_parquet=/databases/annotations/current/hg19/nci60.parquet --annovar-code=nci60 --genome=/databases/genomes/current/hg19.fa --config=/tool/config/config.json --threads=8 
# time howard from_annovar --input=/databases/annovar/current/hg19_cosmic70.txt --output=/databases/annotations/current/hg19/cosmic.vcf.gz --to_parquet=/databases/annotations/current/hg19/cosmic.parquet --annovar-code=cosmic --genome=/databases/genomes/current/hg19.fa --config=/tool/config/config.json --threads=8 
# time howard from_annovar --input=/databases/annovar/current/hg19_clinvar_20221231.txt --output=/databases/annotations/current/hg19/clinvar.vcf.gz --to_parquet=/databases/annotations/current/hg19/clinvar.parquet --annovar-code=clinvar --genome=/databases/genomes/current/hg19.fa --config=/tool/config/config.json --threads=8 

# time howard from_annovar --input=/databases/annovar/current/hg19_gnomad211_exome.txt --output=/databases/annotations/current/hg19/gnomad_exome.vcf.gz --to_parquet=/databases/annotations/current/hg19/gnomad_exome.parquet --annovar-code=gnomad_exome --genome=/databases/genomes/current/hg19.fa --config=/tool/config/config.json --threads=12 --memory=8G

# time howard from_annovar --input=/databases/annovar/current/hg19_dbnsfp42a.txt --output=/databases/annotations/current/hg19/dbnsfp.vcf.gz --to_parquet=/databases/annotations/current/hg19/dbnsfp.parquet --annovar-code=dbnsfp --genome=/databases/genomes/current/hg19.fa --config=/tool/config/config.json --threads=12 --memory=40G --multi_variant=disable --reduce_memory=enable
# time howard from_annovar --input=/databases/annovar/current/hg19_avsnp150.txt --output=/databases/annotations/current/hg19/avsnp150.vcf.gz --to_parquet=/databases/annotations/current/hg19/avsnp150.parquet --annovar-code=avsnp150 --genome=/databases/genomes/current/hg19.fa --config=/tool/config/config.json --threads=12 --memory=40G 


# time howard from_annovar --input=/tmp/nci60.txt --output=/tmp/nci60.vcf.gz --to_parquet=/tmp/nci60.parquet --annovar-code=nci60 --genome=/databases/genomes/current/hg19.fa --config=/tool/config/config.json --threads=8
# time howard from_annovar --input=/tmp/clinvar.txt --output=/tmp/clinvar.vcf.gz --to_parquet=/tmp/clinvar.parquet --annovar-code=clinvar --genome=/databases/genomes/current/hg19.fa --config=/tool/config/config.json --threads=8


# Dictionnaire des types
TYPES = {
    "int": "Integer",
    "int64": "Integer",
    "float": "Float",
    "float64": "Float",
    "object": "String"
}



def from_annovar(args) -> None:

    log.info("Start")

    # Input
    input_file = args.input

    # Output
    output_file = args.output
    
    # Genome
    genome_file = args.genome
    
    # Annovar Code
    annovar_code = args.annovar_code
    
    # To Parquet
    output_file_parquet = args.to_parquet

    # # Export Infos
    # export_infos = args.export_infos

    # # Export Infos Prefix
    # export_infos_prefix = args.export_infos_prefix

    # Multi Variant
    multi_variant = args.multi_variant

    # Reduce memory
    reduce_memory = args.reduce_memory

    # config
    config = args.config

    # Threads
    threads = config.get("threads", None)

    # Threads
    memory = config.get("memory", None)

    # BCFTools
    bcftools = config.get("bcftools", "bcftools")

    # Check parameters

    # Input
    if not os.path.exists(input_file):
        log.error(f"No input file '{input_file}'")
        raise ValueError(f"No input file '{input_file}'")

    # Output
    output_dirname = os.path.dirname(output_file)
    output_file_name, output_file_ext = os.path.splitext(os.path.basename(output_file))
    if output_file_ext not in [".gz"]:
        log.error(f"Output file '{output_file}' without compress '.gz' extension")
        raise ValueError(f"Output file '{output_file}' without compress '.gz' extension")
    if not os.path.exists(output_dirname):
        try:
            os.makedirs(output_dirname, exist_ok=True)
        except:
            log.error(f"Fail create output folder '{output_dirname}'")
            raise ValueError(f"Fail create output folder '{output_dirname}'")

    # To Parquet
    if output_file_parquet:
        output_parquet_dirname = os.path.dirname(output_file_parquet)
        output_file_parquet_name, output_file_parquet_ext = os.path.splitext(os.path.basename(output_file_parquet))
        if output_file_parquet_ext not in [".parquet"]:
            log.error(f"Output file '{output_file_parquet}' without '.parquet' extension")
            raise ValueError(f"Output file '{output_file_parquet}' without '.parquet' extension")
        if not os.path.exists(output_parquet_dirname):
            try:
                os.makedirs(output_parquet_dirname, exist_ok=True)
            except:
                log.error(f"Fail create output folder '{output_parquet_dirname}'")
                raise ValueError(f"Fail create output folder '{output_parquet_dirname}'")

    # Genome
    if not os.path.exists(genome_file):
        log.error(f"No genome file '{genome_file}'")
        raise ValueError(f"No genome file '{genome_file}'")
    if not os.path.exists(genome_file+".fai"):
        log.error(f"No genome index file '{genome_file}.fai'")
        raise ValueError(f"No genome index file '{genome_file}.fai'")
    
    # Annovar code
    if not annovar_code:
        annovar_code = os.path.basename(output_file).replace('.vcf.gz','').replace('.','_')

    # Log
    log.info(f"Input Annovar database: {input_file}")
    log.info(f"Output VCF database: {output_file}")
    if output_file_parquet:
        log.info(f"Output Parquet database: {output_file_parquet}")
    log.info(f"Database name: {annovar_code}")

    # Make output folder
    if not os.path.exists(output_dirname):
        os.makedirs(output_dirname)

    # first VCF

    # Annovar to VCF
    log.info(f"Annovar to VCF and Parquet...")
    annovar_to_vcf(input_file=input_file, output_file=output_file, output_file_parquet=output_file_parquet, annotations=None, header_file=None, database_name=annovar_code, bcftools=bcftools, genome=genome_file, threads=threads, maxmem=memory, remove_annotations=["ID"], reduce_memory=reduce_memory, multi_variant=multi_variant)

    # Header VCF hdr
    log.info(f"VCF Extract header hdr for VCF...")
    command = f"""{bcftools} view -h {output_file} 1>{output_file}.hdr """
    log.debug("bcftools command: " + command)
    subprocess.run(command, shell=True)

    # VCF to Parquet
    if output_file_parquet:
        # File already generated
        # Header VCF hdr
        log.info(f"VCF Extract header hdr for Parquet...")
        command = f"""{bcftools} view -h {output_file} 1>{output_file_parquet}.hdr """
        log.debug("bcftools command: " + command)
        subprocess.run(command, shell=True)

    log.info("End")


   
def annovar_to_vcf(input_file:str, output_file:str, output_file_parquet:str = None, annotations:str = None, header_file:str = None, database_name:str = None, bcftools:str = "bcftools", genome:str = "hg19.fa", threads:int = None, maxmem:str = None, remove_annotations:list = [], reduce_memory:bool = None, multi_variant:bool = None) -> None: 
    """
    This function converts an ANNOVAR file to a VCF file and optionally to a Parquet file, with various
    options for annotations, headers, databases, and memory usage.
    
    :param input_file: The path to the input file in ANNOVAR format that needs to be converted to VCF
    format
    :type input_file: str
    :param output_file: The name of the output VCF file that will be generated by the function
    :type output_file: str
    :param output_file_parquet: output_file_parquet is an optional parameter that specifies the name of
    the output file in Parquet format. If this parameter is not provided, the output will not be saved
    in Parquet format
    :type output_file_parquet: str
    :param annotations: This parameter is used to specify the location of the ANNOVAR annotation
    database files. If not provided, ANNOVAR will use the default location
    :type annotations: str
    :param header_file: The path to a file containing the header information for the VCF output. This
    can be used to customize the output format of the VCF file. If not provided, a default header will
    be used
    :type header_file: str
    :param database_name: The name of the ANNOVAR database used for annotation
    :type database_name: str
    :param bcftools: The path to the bcftools executable, defaults to bcftools
    :type bcftools: str (optional)
    :param genome: The genome parameter specifies the reference genome file to be used for the
    conversion from annovar format to VCF format, defaults to hg19.fa
    :type genome: str (optional)
    :param threads: The number of threads to use for processing. This can speed up the process if your
    computer has multiple cores
    :type threads: int
    :param maxmem: The maximum amount of memory that can be used by the program. It is usually specified
    in units of bytes, kilobytes, megabytes, or gigabytes. For example, "2G" means 2 gigabytes of memory
    :type maxmem: str
    :param remove_annotations: `remove_annotations` is a list of annotations to be removed from the
    output VCF file. These annotations will not be included in the final VCF file
    :type remove_annotations: list
    :param reduce_memory: A boolean parameter that determines whether to reduce memory usage during the
    conversion process. If set to True, the function will attempt to reduce memory usage by using a more
    memory-efficient algorithm, but this may result in slower performance. If set to False, the function
    will use a faster algorithm that may consume more, defaults to False
    :type reduce_memory: bool (optional)
    :param multivariant: A boolean parameter that determines if input file contains multiple annotations
    for each variant (position ref alt). If set to False, the function will attempt to reduce memory usage
    a specific query without 'group by', for a more memory-efficient algorithm. If set to True, the function
    will use a query using 'group by', which may consume more memory. I set to None, the function will
    auto-detemine the parameter value with a sample of variants. Defaults to None (auto)
    :type multivariant: bool (optional)
    """
    
    # Check input file
    if not os.path.exists(input_file):
        raise ValueError("No input file")
        
    # Log
    log.debug("input_file: " + input_file)

    # threads
    threads_connexion = threads
    if not threads:
        threads = 1
    
    # maxmem
    maxmem_connexion = maxmem
    if not maxmem:
        maxmem = "40G"

    # Determine header columns
    header = []
    header_ok = True
    if header_file:
        with open(header_file, 'r') as f:
            header = f.readline().strip().split('\t')
        f.close()
    else:
        
        with open(input_file, 'r') as f:
            for line in f:
                if str(line).startswith("##"):                
                    continue
                elif len(line.strip().split('\t')) < 5:
                    header_ok = False
                elif str(line).startswith("#"):
                    log.debug("Found header")
                    header = line.strip().split('\t')
                    if "#CHROM" not in header:
                        chrom_list_for_header = ["#Chr"]
                        for chr in chrom_list_for_header:
                            if chr in header:
                                header[header.index(chr)] = "#CHROM"
                                log.debug(f"found '{chr}' as chromosome column")
                        if "#CHROM" not in header:
                            log.debug(f"default chromosome column is 0")
                            header[0] = "#CHROM"
                    if "POS" not in header:
                        pos_list_for_header = ["Start"]
                        for pos in pos_list_for_header:
                            if pos in header:
                                header[header.index(pos)] = "POS"
                                log.debug(f"found '{pos}' as position column")
                        if "POS" not in header:
                            log.debug(f"default position column is 1")
                            header[1] = "POS"
                    if "REF" not in header:
                        ref_list_for_header = ["Ref"]
                        for ref in ref_list_for_header:
                            if ref in header:
                                header[header.index(ref)] = "REF"
                                log.debug(f"found '{ref}' as reference column")
                        if "REF" not in header:
                            log.debug(f"default reference column is 3")
                            header[3] = "REF"
                    if "ALT" not in header:
                        alt_list_for_header = ["Alt"]
                        for alt in alt_list_for_header:
                            if alt in header:
                                header[header.index(alt)] = "ALT"
                                log.debug(f"found '{alt}' as alternative column")
                        if "ALT" not in header:
                            log.debug(f"default alternative column is 4")
                            header[4] = "ALT"
                    break
                else:
                    log.debug("NO Header")
                    header = line.strip().split('\t')
                    if len(header) >= 5:
                        log.debug(header)
                        header[0] = "#CHROM"
                        header[1] = "POS"
                        header[2] = "ID"
                        header[3] = "REF"
                        header[4] = "ALT"
                        for h in enumerate(header[5:]):
                            if database_name:
                                prefix = database_name
                            else:
                                prefix = "column"
                            info_tag_name = str(h[0]+1).replace('-', '_').replace('+', '_')
                            if h[0]:
                                column_name = prefix + "_" + info_tag_name
                            else:
                                column_name = prefix
                            if column_name[0].isdigit():
                                column_name = "A" + column_name
                            header[h[0]+5] = column_name
                    else:
                        header_ok = False
                    break
        f.close()

    if not header_ok:
        raise ValueError("Error in header")

    # protect info tag from unauthorized characters
    header_fixed = []
    for h in header:
        h = h.replace('-', '_').replace('+', '').replace('.', '_')
        header_fixed.append(h)
    header = header_fixed

    # determine nb header line
    nb_header_line = 0
    with open(input_file, 'r') as f:
        for line in f:
            if str(line).startswith("##"):
                nb_header_line += 1
                continue
            elif str(line).startswith("#"):
                nb_header_line += 1
            else:
                break
    f.close()
    
    # Check dtype and force if needed
    dtype = {}
    for h in header:
        column_type=None
        if h in ["#CHROM", "REF", "ALT"]:
            column_type = str
            dtype[h] = column_type
        elif h in ["POS", "START", "END"]:
            column_type = int
            dtype[h] = column_type

    # Check format VCF readable
    #nrows = 100000
    log.debug("Check input file struct...")
    nrows_sampling = 1000000
    try:
        df = pd.read_csv(input_file, delimiter='\t', nrows=nrows_sampling, na_values=['.'], header=nb_header_line, names=header, dtype=dtype)
    except:
        log.error("format not supported")
        raise ValueError("format not supported")
    
    # Check multi variant (multiple annotation on same variant CHROM, POS REF, ALT)
    query_multi_variant = """
        SELECT "#CHROM", "POS", "REF", "ALT", count(*) AS count
        FROM df
        GROUP BY "#CHROM", "POS", "REF", "ALT"
        ORDER BY "count" DESC
    """
    query_multi_variant_connexion = duckdb.connect(":memory:")
    res_multi_variant = query_multi_variant_connexion.query(query_multi_variant)
                   
    # Check Multi Variant
    if multi_variant is None or multi_variant.lower() == "auto":
        if res_multi_variant.df()["count"][0] > 1:
            multi_variant = True
            log.debug("""Multi Variant mode enabled: because some variant with multiple annotation lines""")
        else:
            multi_variant = False
    elif multi_variant.lower().startswith("enable"):
        multi_variant = True
    elif multi_variant.lower().startswith("disable"):
        multi_variant = False
    else:
        multi_variant = True

    # Log
    if multi_variant:
        log.info("Multi Variant mode enabled")
    else:
        log.info("Multi Variant mode disabled")

    # Check reduce memory
    if reduce_memory is None or reduce_memory.lower() == "auto":
        log.debug("""Reduce memory mode as None (Auto)""")
        if multi_variant:
            reduce_memory = True
            log.debug("""Reduce memory mode enabled: because Multi Variant mode enabled""")
        else:
            # Check number of variants
            if len(df) == nrows_sampling:
                reduce_memory = True
                log.debug(f"""Reduce memory mode enabled: because more than {nrows_sampling} variants""")
            else:
                reduce_memory = False
    elif reduce_memory.lower().startswith("enable"):
        reduce_memory = True
    elif reduce_memory.lower().startswith("disable"):
        reduce_memory = False
    else:
        reduce_memory = True

    # Log
    if reduce_memory:
        log.info("""Reduce memory mode enabled""")
    else:
        log.info("""Reduce memory mode disabled""")
    
    # Create connexion
    if reduce_memory:
        #log.info("""Reduce memory mode enabled""")
        connexion_type = f"{output_file}.tmp.duckdb"
        remove_if_exists([connexion_type])
    else:
        #log.info("""Reduce memory mode disabled""")
        connexion_type = ":memory:"
    duckdb_config = {}
    if threads_connexion:
        duckdb_config["threads"] = threads_connexion
    if maxmem_connexion:
        duckdb_config["memory_limit"] = maxmem_connexion

    # Connexion
    conn = duckdb.connect(connexion_type, config=duckdb_config)
    delimiter = '\t'

    # Recheck column type les types de chaque colonne
    for col in df.columns:
        if df[col].dtype == object:
            if pd.to_numeric(df[col], errors='coerce').notnull().all():
                df[col] = pd.to_numeric(df[col], errors='coerce').astype(float)
        else:
            df[col] = pd.to_numeric(df[col], errors='coerce')
    column_types = df.dtypes

    # Dictionnary to mapper data types from pandas to SQL
    sql_types_forces = {
        "#CHROM": "STRING",
        "POS": "INTEGER",
        "REF": "STRING",
        "ALT": "STRING",
    }

    # replace if dtype forced
    dtype_final = {}
    for column in dict(column_types):
        if sql_types_forces.get(column,None):
            dtype_final[column] = sql_types_forces.get(column,None)
        else:
            dtype_final[column] = "STRING"

    # column type to VCF type
    vcf_types = {col: TYPES[str(column_types[col])] for col in column_types.index}

    # check chromosomes

    # Indexes
    input_file_index = input_file + ".idx"
    genome_index = genome + ".fai"
    
    # Find chromosomes list in input file
    log.debug("Check original chromosomes")

    input_file_chromosomes = []

    if os.path.exists(input_file_index):

        # Log
        log.debug("Check original chromosomes - from idx file...")

        # Read index of input file
        df = pd.read_csv(input_file_index, delimiter='\t', header=1,  names = ["#CHROM", "col1", "col2", "col3"], dtype = {"#CHROM": str})

        # Unique values
        input_file_chromosomes = df['#CHROM'].unique()

    else:
        
        # Log
        log.debug("Check original chromosomes - from input file...")

        # Set config value
        csv.field_size_limit(sys.maxsize // 10)
        
        # Find chromosomes in corresponding column
        unique_values = set()
        chromosome_index = header.index("#CHROM")
        with open(input_file, newline='') as csvfile:
            reader = csv.reader(csvfile, delimiter="\t")
            for row in reader:
                unique_values.add(row[chromosome_index])

        # Unique values
        input_file_chromosomes = unique_values
    
    # Init lists
    chrom_list_fixed = []
    chrom_list_original = []

    if os.path.exists(genome_index):

        # Check chromosomes in genome
        df_genome = pd.read_csv(genome_index, delimiter='\t', header=None,  names = ["CHROM", "length", "col2", "col3", "col4"], dtype = {"CHROM": str, "Start": int, "End": int})
        genome_chrom_length = {}
        for item in df_genome.iterrows():
            genome_chrom_length[item[1].CHROM] = item[1].length

        # Uniquify
        genome_chrom = df_genome['CHROM'].unique()

        # Fix chromosomes
        for chrom in input_file_chromosomes:
            chrom_fixed = chrom
            chrom_list_original.append(chrom)
            if not str(chrom).startswith("chr"):
                chrom_fixed = "chr" + str(chrom_fixed)
            if chrom_fixed == "chrMT":
                chrom_fixed = "chrM"
            if chrom_fixed == "chr23":
                chrom_fixed = "chrX"
            if chrom_fixed == "chr24":
                chrom_fixed = "chrY"
            if chrom_fixed in genome_chrom:
                chrom_list_fixed.append(chrom_fixed)

        chrom_list_map = {k: v for k, v in zip(chrom_list_original, chrom_list_fixed)}

        # Log
        log.debug("List of original chromosomes: " + str(chrom_list_original))
        log.debug("List of fixed chromosomes: " + str(chrom_list_fixed))

    # Write VCF header

    # open file input in read only, open (temporary) vcf output in write mode
    with bgzf.open(output_file+".tmp.translation.header.vcf.gz", 'wt') as vcf_file:

        # header fileformat
        vcf_file.write('##fileformat=VCFv4.3\n')

        # header INFO TAGS
        for col, type_ in vcf_types.items():
            if col not in ["#CHROM", "POS", "REF", "ALT"] + remove_annotations:
                info_tag = col
                info_tag.replace("-", "_").replace("+", "_")
                vcf_file.write(f"""##INFO=<ID={info_tag},Number=.,Type={type_},Description="{info_tag} annotation">\n""")

        # header contig
        for chrom in chrom_list_fixed:
            log.debug(f"""##contig=<ID={chrom},length={genome_chrom_length[chrom]}>""")
            vcf_file.write(f"""##contig=<ID={chrom},length={genome_chrom_length[chrom]}>\n""")

        # header #CHROM line
        vcf_file.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n')

    # Determine annotations if not annotations in input
    annotations_auto = []
    if not annotations:
        for h in header:
            if h not in ["#CHROM", "POS", "REF", "ALT", "-", ""] and not str(h).startswith("-"):
                #log.debug(h)
                annotations_auto.append(h)
        annotations = annotations_auto

    for remove_annotation in remove_annotations:
        if remove_annotation in annotations:
            annotations.remove(remove_annotation)

    # Check dtypes
    dtype_duckdb = []
    for dtype in dtype_final:
        dtype_duckdb.append(dtype_final[dtype])

    # Columns struct
    columns_struct = {k: v for k, v in zip(header, dtype_duckdb)}

    # Log
    log.debug(f"""duckDB Connexion: {connexion_type}""")
    log.debug(f"""duckDB Config: {duckdb_config}""")

    # Create view of input file
    if reduce_memory and multi_variant:

        # Log
        log.debug(f"Create View From TSV to Parquet")
            
        # Create Parquet file from TSV
        tsv_to_parquet(input_file, f'{output_file}.tmp.annovar.parquet', delim=delimiter, columns=columns_struct, quote=None, nullstr='.', skip=nb_header_line)

        # Query Create view
        query = f""" CREATE VIEW annovar AS SELECT * FROM '{output_file}.tmp.annovar.parquet' """

    else:

        # Log
        log.debug("Create View From TSV")

        # Create view of input file
        query = f""" CREATE VIEW annovar AS SELECT * FROM read_csv_auto('{input_file}', delim='{delimiter}', columns={columns_struct}, names={header}, dtypes={dtype_duckdb}, quote=None, nullstr='.', parallel=True, skip={nb_header_line}) """
    
    conn.execute(query)

    # Check existing columns
    log.debug("Check existing columns")
    query = """ SELECT * FROM annovar LIMIT 0  """
    columns = conn.execute(query).df().columns.tolist()

    # Prepare queries

    # Prepare queries - main columns
    main_columns = """
        CASE WHEN "#CHROM" LIKE 'chr%' THEN '' ELSE 'chr' END ||
        CASE WHEN "#CHROM" = 'MT' THEN 'M' WHEN "#CHROM" = '23' THEN 'X' WHEN "#CHROM" = '24' THEN 'Y' ELSE "#CHROM" END
        AS "#CHROM",
        CASE WHEN "REF" in ('-','') OR "ALT" in ('-','') THEN "POS"-1 ELSE "POS" END
        AS "POS",
        '' AS "ID",
        CASE WHEN "REF" in ('-','') THEN 'N' WHEN "ALT" in ('-','') THEN 'N' || "REF" ELSE "REF" END
        AS "REF",
        CASE WHEN "REF" in ('-','') THEN 'N' || "ALT" WHEN "ALT" in ('-','') THEN 'N' ELSE "ALT" END
        AS "ALT",
        '' AS "QUAL",
        'PASS' AS "FILTER",
    """

    # Prepare queries - other columns
    any_value_list = []
    nb_annotation_column = 0
    for column in columns:
        if column not in ["#CHROM", "POS", "REF", "ALT", "ID"]:
            nb_annotation_column += 1
            if multi_variant:
                any_value_list.append(f""" 
                    CASE WHEN STRING_AGG(DISTINCT "{column}") IS NULL THEN '' ELSE '{column}=' || REPLACE(STRING_AGG(DISTINCT "{column}"),';',',') || ';' END
                """)
            else:
                any_value_list.append(f""" 
                    CASE WHEN "{column}" IS NULL THEN '' ELSE '{column}=' || REPLACE("{column}",';',',') || ';' END
                """)

    # Join to create INFO column
    any_value_sql = " || ".join(any_value_list)

    # Remove last caracter ';'
    any_value_sql = f""" SUBSTR({any_value_sql}, 1, LENGTH({any_value_sql}) - 1) """

    # Create parquet table/view
    log.info("Formatting VCF and Parquet...")

    # Multi Variant Group By
    if multi_variant:
        query_group_by = """ GROUP BY "#CHROM", "POS", "REF", "ALT" """
    else:
        query_group_by = ""

    if reduce_memory and multi_variant:

        # Check list of chromosome
        log.debug("Check chromosomes...")
        if not chrom_list_original:
            query_chrom_list = """ SELECT distinct("#CHROM") FROM annovar """
            chrom_list = list(conn.query(query_chrom_list).df()["#CHROM"])
        else:
            chrom_list = chrom_list_original

        # Window calculation
        window_base = 100000000
        if multi_variant:
            window = round( window_base / nb_annotation_column )
        else:
            window = window_base

        # Insert formatted variants
        for chrom in chrom_list:

            log.info(f"Formatting VCF and Parquet - Chromosome '{chrom}'...")

            # max pos
            log.debug(f"Formatting VCF and Parquet - Chromosome '{chrom}' - Check range...")

            min_pos = 0
            max_pos = genome_chrom_length.get(chrom_list_map.get(chrom,0))

            start = min_pos - 1
            stop = start + window

            while True:

                # Log
                log.debug(f"Formatting VCF and Parquet - Chromosome '{chrom}' - Range {start}-{stop}...")

                # If out of bound
                if start > max_pos or max_pos == 0:
                    break

                # Parquet File
                #log.info(f"Formatting VCF and Parquet - Chromosome '{chrom}' - Exporting...")
                log.debug(f"Create parquet file for chromosome '{chrom}' range {start}-{stop}...")
                query_select_parquet = f"""
                    COPY
                        (
                            SELECT *
                            FROM annovar
                            WHERE "#CHROM" = '{chrom}'
                            AND POS > {start} AND POS <= {stop}
                        )
                    TO '{output_file}.tmp.chrom.parquet' WITH (FORMAT PARQUET)
                """
                res = conn.execute(query_select_parquet)
                start = start + window
                stop = start + window

                # Formatted Parquet file
                #log.info(f"Formatting VCF and Parquet - Chromosome '{chrom}' - Formatting...")
                log.debug(f"Create formatted parquet file for chromosome '{chrom}' range {start}-{stop}...")
                query_select_parquet = f"""
                    COPY
                        (
                            SELECT
                                {main_columns}
                                {any_value_sql} AS "INFO"
                            FROM '{output_file}.tmp.chrom.parquet'
                            {query_group_by}
                        )
                    TO '{output_file}.tmp.chrom.{chrom}.start.{start}.stop.{stop}.parquet' WITH (FORMAT PARQUET)
                """
                res = conn.execute(query_select_parquet)
                count_insert = res.df()["Count"][0] 
                log.debug(f"Create formatted parquet file for chromosome '{chrom}' range {start}-{stop}... {count_insert} variants")
                del res
                gc.collect()

        # Query creation table
        query_parquet = f""" CREATE TABLE parquet AS (SELECT * FROM parquet_scan('{output_file}.tmp.chrom.*.start.*.stop.*.parquet')) """
        res = conn.execute(query_parquet)
        del res
        gc.collect()


    else:

        # Query to select columns
        query_select_parquet = f"""
            SELECT
                {main_columns}
                {any_value_sql} AS "INFO"
            FROM annovar
            {query_group_by} 
        """

        # Log
        log.debug("CREATE view parquet...")

        # Query creation view
        query_parquet = f""" CREATE VIEW parquet AS ({query_select_parquet}) """
        res = conn.execute(query_parquet)
        del res
        gc.collect()


    # Export VCF
    log.info("Exporting VCF...")

    # Check source
    source = """ SELECT "#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO" FROM parquet """

    # Query to copy file
    query = f"""
        COPY ({source})
        TO '{output_file}.tmp.translation.variants.vcf.gz' WITH (FORMAT CSV, COMPRESSION GZIP, DELIMITER '\t', QUOTE '')
    """
    
    # Copy file
    res = conn.execute(query)
    del res
    gc.collect()

    # Close connexion
    log.debug("Close connexion...")
    conn.close()

    # BCFTools to reheader sort normalize
    log.info("Normalizing VCF...")

    # Log
    log.debug("VCF Sorting and Normalization...")

    # Command
    command = f"""zcat {output_file}.tmp.translation.header.vcf.gz {output_file}.tmp.translation.variants.vcf.gz 2>{output_file}.tmp.err | {bcftools} sort --max-mem={maxmem} 2>{output_file}.tmp.err | {bcftools} norm --threads={threads} --check-ref s -f {genome} -Oz -o {output_file} 2>{output_file}.tmp.err """

    # Log
    log.debug("bcftools command: " + command)

    # Run
    subprocess.run(command, shell=True)


    # Export Parquet
    if output_file_parquet:

        # Log
        log.info("Exporting Parquet...")

        # Check source
        parquet_info_explode(input_file=output_file, output_file=output_file_parquet, threads=threads, memory=maxmem_connexion, reduce_memory=reduce_memory)


    # Clean
    clean_command = f""" rm -rf {output_file}.tmp.* {output_file_parquet}.tmp.* """
    subprocess.run(clean_command, shell=True)



def parquet_info_explode(input_file:str, output_file:str, threads:int = None, memory:str = None, reduce_memory:bool = False) -> None:
    """
    This function takes a parquet file, splits it by chromosome, explodes the INFO column, and then
    merges the exploded files back together.
    
    :param input_file: The path to the input file, which can be either a TSV or VCF file
    :type input_file: str
    :param output_file: The name of the output file in Parquet format after exploding the input file
    :type output_file: str
    :param threads: The number of threads to use for processing the parquet file, defaults to None (all)
    :type threads: int (optional)
    :param memory: The among of memory to use for processing the parquet file, defaults to None (all)
    :type memory: str (optional)
    :param reduce_memory: The `reduce_memory` parameter is a boolean flag that determines whether or not
    to use memory reduction techniques during the execution of the function. If set to `True`, the
    function will attempt to reduce memory usage during the execution, which may result in slower
    performance but lower memory usage. If set to `, defaults to False
    :type reduce_memory: bool (optional)
    """

    # Config
    config_ro = {
        "access": "RO"
    }
    config_threads = {}
    param_explode = {
        "explode_infos": True,
        "export_extra_infos": True
    }

    # Threads
    if threads:
        config_ro["threads"] = threads
        config_threads["threads"] = threads

    # Memory
    if memory:
        config_ro["memory_limit"] = memory
        config_threads["memory_limit"] = memory

    if reduce_memory:

        # List of exploded parquet files
        list_of_exploded_files = []

        # Create 
        parquet_input = Variants(input=input_file, config=config_ro)
        parquet_chromosomes = parquet_input.get_query_to_df(query=f""" SELECT "#CHROM" FROM read_csv('{input_file}',AUTO_DETECT=TRUE) WHERE "#CHROM" NOT NULL GROUP BY "#CHROM" """)

        for chrom in parquet_chromosomes["#CHROM"]:
            
            # Split Chrom
            query_chrom = f""" SELECT * FROM read_csv('{input_file}',AUTO_DETECT=TRUE) WHERE "#CHROM"='{chrom}'  """
            query_chrom_output = f"{output_file}.{chrom}.parquet"
            query_chrom_explode_output = f"{output_file}.{chrom}.explode.parquet"

            # Log
            log.info(f"Explode infos chromosome {chrom}")

            # Extract
            log.debug(f"Extract chromosome {chrom}: {query_chrom_explode_output}")
            remove_if_exists([query_chrom_output,query_chrom_output+".hdr"])
            parquet_input.export_output(output_file=query_chrom_output, query=query_chrom, export_header=True)
            
            # Explode
            log.debug(f"Explode infos for chromosome {chrom}: {query_chrom_explode_output}")
            list_of_exploded_files.append(query_chrom_explode_output)
            parquet_explode = Variants(input=query_chrom_output, output=query_chrom_explode_output, config=config_threads, param=param_explode)
            parquet_explode.load_data()
            remove_if_exists([query_chrom_explode_output,query_chrom_explode_output+".hdr"])
            parquet_explode.export_output()
            remove_if_exists([query_chrom_output,query_chrom_output+".hdr"])

        # list_of_exploded_files
        log.info(f"Merge explode infos")
        log.debug(f"Merge explode infos files: {list_of_exploded_files}")
        query_explode = f""" SELECT * FROM read_parquet({list_of_exploded_files}) """
        query_explode_output = output_file
        remove_if_exists([query_explode_output,query_explode_output+".hdr"])
        parquet_explode.export_output(output_file=query_explode_output, query=query_explode, export_header=True)

        # Clear
        for file in list_of_exploded_files:
            remove_if_exists([file,file+".hdr"])

    else:

        # Create 
        parquet_explode = Variants(input=input_file, output=output_file, config=config_threads, param=param_explode)
        parquet_explode.load_data()
        parquet_explode.export_output()

    # Check
    # df = parquet_input.get_query_to_df(f"SELECT * FROM '{query_explode_output}' LIMIT 10 ")
    # log.debug(df)


def tsv_to_parquet(tsv:str, parquet:str, delim:str = None, columns:dict = None, quote:str = None, nullstr:str = None, skip:int = None) -> None:
    """
    The function converts a TSV file to a Parquet file with customizable options.
    
    :param tsv: The path to the TSV file that needs to be converted to Parquet format
    :type tsv: str
    :param parquet: `parquet` is the file path and name of the output Parquet file that will be created
    by the function
    :type parquet: str
    :param delim: The delimiter used in the TSV file to separate columns. If not specified, the default
    delimiter (tab) will be used
    :type delim: str
    :param columns: The `columns` parameter is a dictionary that maps column names to their data types.
    It is used to specify the schema of the resulting Parquet file. For example, if the input TSV file
    has columns "name", "age", and "salary", and we want "name" to be
    :type columns: dict
    :param quote: The `quote` parameter is an optional parameter that specifies the character used to
    quote fields in the TSV file. If not specified, the default quote character is double quotes (")
    :type quote: str
    :param nullstr: The `nullstr` parameter is used to specify the string that represents null values in
    the input TSV file. This parameter is used to correctly interpret and convert null values in the TSV
    file to null values in the resulting Parquet file. For example, if the null value in the TSV
    :type nullstr: str
    :param skip: The `skip` parameter is an optional integer parameter that specifies the number of rows
    to skip at the beginning of the TSV file. This is useful if the TSV file has a header row that
    should not be included in the resulting Parquet file. If `skip` is not specified, no
    :type skip: int
    """

    # Create Schema dict
    columns_pyarrow_dict = []
    for col in columns:
        if columns[col] in ["STRING"]:
            columns_pyarrow_dict.append((col,pa.string()))
        elif columns[col] in ["INTEGER"]:
            columns_pyarrow_dict.append((col,pa.int64()))
        elif columns[col] in ["FLOAT"]:
            columns_pyarrow_dict.append((col,pa.float64()))

    # Create Schema
    schema = pa.schema(columns_pyarrow_dict)

    # CSV options
    convert_options = csv.ConvertOptions()
    read_options = csv.ReadOptions()
    parse_options = csv.ParseOptions()

    # Parameters to VCS options
    convert_options.auto_dict_encode = False
    if delim:
        parse_options.delimiter = delim
    if columns:
        convert_options.column_types = schema
        read_options.column_names = list(columns.keys())
    if quote:
        parse_options.quote_char = quote
    if nullstr:
        convert_options.null_values = nullstr
        convert_options.strings_can_be_null = True
    if skip:
        read_options.skip_rows = skip

    # Read CSV and Write to Parquet
    writer = None
    with csv.open_csv(tsv, convert_options=convert_options, read_options=read_options, parse_options=parse_options) as reader:
        for next_chunk in reader:
            if next_chunk is None:
                break
            if writer is None:
                writer = pq.ParquetWriter(parquet, next_chunk.schema)
            next_table = pa.Table.from_batches([next_chunk])
            writer.write_table(next_table)
    writer.close()

