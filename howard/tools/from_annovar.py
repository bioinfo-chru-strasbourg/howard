#!/usr/bin/env python

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


from howard.objects.variants import Variants
from howard.objects.annotation import Annotation
from howard.commons import *


# Usage
# time howard from_annovar --input=/databases/annovar/current/hg19_nci60.txt --output=/databases/annotations/current/hg19/nci60.vcf.gz --to_parquet=/databases/annotations/current/hg19/nci60.parquet --annovar-code=nci60 --genome=/databases/genomes/current/hg19.fa --config=/tool/config/config.json --threads=8 
# time howard from_annovar --input=/databases/annovar/current/hg19_cosmic70.txt --output=/databases/annotations/current/hg19/cosmic.vcf.gz --to_parquet=/databases/annotations/current/hg19/cosmic.parquet --annovar-code=cosmic --genome=/databases/genomes/current/hg19.fa --config=/tool/config/config.json --threads=8 
# time howard from_annovar --input=/databases/annovar/current/hg19_clinvar_20221231.txt --output=/databases/annotations/current/hg19/clinvar.vcf.gz --to_parquet=/databases/annotations/current/hg19/clinvar.parquet --annovar-code=clinvar --genome=/databases/genomes/current/hg19.fa --config=/tool/config/config.json --threads=8 
# time howard from_annovar --input=/databases/annovar/current/hg19_gnomad211_exome.txt --output=/databases/annotations/current/hg19/gnomad_exome.vcf.gz --to_parquet=/databases/annotations/current/hg19/gnomad_exome.parquet --annovar-code=gnomad_exome --genome=/databases/genomes/current/hg19.fa --config=/tool/config/config.json --threads=12 --reduce_memory

# time howard from_annovar --input=/databases/annovar/current/hg19_dbnsfp42a.txt --output=/databases/annotations/current/hg19/dbnsfp.vcf.gz --to_parquet=/databases/annotations/current/hg19/dbnsfp.parquet --annovar-code=dbnsfp --genome=/databases/genomes/current/hg19.fa --config=/tool/config/config.json --threads=12  --reduce_memory
# time howard from_annovar --input=/databases/annovar/current/hg19_avsnp150.txt --output=/databases/annotations/current/hg19/avsnp150.vcf.gz --to_parquet=/databases/annotations/current/hg19/avsnp150.parquet --annovar-code=avsnp150 --genome=/databases/genomes/current/hg19.fa --config=/tool/config/config.json --threads=12 --reduce_memory

# 
# gnomad211_exome
# dbnsfp42a
# avsnp150


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

    # Export Infos Prefix
    reduce_memory = args.reduce_memory

    # config
    config = args.config

    # Threads
    threads = int(config.get("threads", 1))

    # Threads
    bcftools = config.get("bcftools", "bcftools")

    # Check parameters

    # Input
    if not os.path.exists(input_file):
        log.error(f"No input file '{input_file}'")
        raise ValueError(f"No input file '{input_file}'")
    if not os.path.exists(input_file+".idx"):
        log.error(f"No input index file '{input_file}.idx'")
        raise ValueError(f"No input index file '{input_file}.idx'")

    # Output
    output_dirname = os.path.dirname(output_file)
    output_file_name, output_file_ext = os.path.splitext(os.path.basename(output_file))
    if output_file_ext not in [".gz"]:
        log.error(f"Output file '{output_file}' without compress extension")
        raise ValueError(f"Output file '{output_file}' without compress extension")
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
            log.error(f"Output file '{output_file_parquet}' without compress extension")
            raise ValueError(f"Output file '{output_file_parquet}' without compress extension")
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
    log.info(f"Database name: {annovar_code}")

    # Make output folder
    if not os.path.exists(output_dirname):
        os.makedirs(output_dirname)

    # first VCF

    # Annovar to VCF
    log.info(f"Annovar to VCF and Parquet...")
    annovar_to_vcf(input_file=input_file, output_file=output_file, output_file_parquet=output_file_parquet, annotations=None, header_file=None, database_name=annovar_code, bcftools=bcftools, genome=genome_file, threads=threads, maxmem="40G", remove_annotations=["ID"], reduce_memory=reduce_memory)

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


def annovar_to_vcf(input_file, output_file, output_file_parquet=None, annotations=None, header_file=None, database_name=None, bcftools="bcftools", genome="hg19.fa", threads=1, maxmem="40G", remove_annotations:list = [], reduce_memory=False): #, chrom_col='#CHROM', pos_col='POS', ref_col='REF', alt_col='ALT'):
    """
    This function reads a tab-separated file, converts it to a VCF file, and then normalizes it, with
    options to include specific columns in the INFO field, use a header file, and specify a reference
    genome.
    
    :param input_file: The path to the input file that needs to be converted to VCF format and
    normalized
    :param output_file: The name of the output file that will be generated by the function
    :param output_file_parquet: `output_file_parquet` is an optional parameter that specifies the name
    of the output file in Parquet format. If this parameter is not provided, the function will not
    export the data in Parquet format
    :param annotations: annotations is a list of columns to be included in the INFO field of the VCF. If
    not provided, all columns will be included
    :param header_file: The path to the header file. If provided, it will be used to determine the
    columns of the input file. If not provided, the function will try to determine the columns
    automatically
    :param database_name: The name of the database you want to use for the INFO tags in the VCF file. If
    not provided, a default name will be used
    :param bcftools: bcftools is a software tool used for manipulating and analyzing VCF (Variant Call
    Format) files. It is used in this function to sort and normalize the VCF file, defaults to bcftools
    (optional)
    :param genome: The reference genome to use for the VCF conversion and normalization. It defaults to
    hg19.fa, defaults to hg19.fa (optional)
    :param threads: The number of threads to use for sorting during the conversion process, defaults to
    1 (optional)
    :param maxmem: maxmem is the maximum amount of memory to use for sorting during the VCF sorting and
    normalization step. It is a string that specifies the amount of memory, such as "40G" for 40
    gigabytes, defaults to 40G (optional)
    :param remove_annotations: The list of annotations to be removed from the INFO field of the VCF
    :type remove_annotations: list
    :return: The function does not explicitly return anything, but it does write a VCF file and
    optionally a Parquet file to disk. The function also logs various messages to the console.
    """
    
    # Check input file
    if not os.path.exists(input_file):
        raise ValueError("No input file")
        
    # Log
    log.debug("input_file: " + input_file)

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
    try:
        df = pd.read_csv(input_file, delimiter='\t', nrows=100000, na_values=['.'], header=nb_header_line, names=header, dtype=dtype)
    except:
        log.error("format not supported")
        raise ValueError("format not supported")
    
    # Recheck column type les types de chaque colonne
    for col in df.columns:
        if df[col].dtype == object:
            if pd.to_numeric(df[col], errors='coerce').notnull().all():
                df[col] = pd.to_numeric(df[col], errors='coerce').astype(float)
        else:
            df[col] = pd.to_numeric(df[col], errors='coerce')
    column_types = df.dtypes

    # Dictionnary to mapper data types from pandas to SQL
    sql_types = {
        'int64': 'INTEGER',
        'float64': 'FLOAT',
        'object': 'STRING',
    }

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
            dtype_final[column] = sql_types.get(str(column_types[column]),'STRING')

    # column type to VCF type
    vcf_types = {col: TYPES[str(column_types[col])] for col in column_types.index}

    # check chromosomes
    input_file_index = input_file + ".idx"
    genome_index = genome + ".fai"
        
    chrom_list_fixed = []
    if os.path.exists(input_file_index) and os.path.exists(genome_index):

        # Check chromosomes in genome
        df_genome = pd.read_csv(genome_index, delimiter='\t', header=None,  names = ["CHROM", "length", "col2", "col3", "col4"], dtype = {"CHROM": str, "Start": int, "End": int})
        genome_chrom_length = {}
        for item in df_genome.iterrows():
            genome_chrom_length[item[1].CHROM] = item[1].length

        # Uniquify
        genome_chrom = df_genome['CHROM'].unique()

        # Read index of input file
        df = pd.read_csv(input_file_index, delimiter='\t', header=1,  names = ["#CHROM", "col1", "col2", "col3"], dtype = {"#CHROM": str})
        # Fix chromosomes
        for chrom in df['#CHROM'].unique():
            chrom_fixed = chrom
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
        log.debug("List of chromosomes: " + str(chrom_list_fixed))
        log.debug(df_genome)

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

    # Create view
    conn = duckdb.connect(":memory:", config={"threads": threads, "memory_limit": maxmem})
    delimiter = '\t'

    query = f""" CREATE VIEW annovar AS SELECT * FROM read_csv_auto('{input_file}', auto_detect=True, delim='{delimiter}', names={header}, dtypes={dtype_duckdb}, quote=None, nullstr='.', parallel=True) """
    conn.execute(query)

    # Check existing columns
    query = """ SELECT * FROM annovar LIMIT 0  """
    columns = conn.execute(query).df().columns.tolist()

    # Prepare queries
    any_value_list = []
    explode_info_list = []
    for column in columns:
        if column not in ["#CHROM", "POS", "REF", "ALT", "ID"]:
            any_value_list.append(f""" 
                CASE WHEN STRING_AGG(DISTINCT "{column}") IS NULL THEN '' ELSE '{column}=' || REPLACE(STRING_AGG(DISTINCT "{column}"),';',',') || ';' END
            """)
            explode_info_list.append(f"""
                CASE WHEN STRING_AGG(DISTINCT "{column}") IS NULL THEN NULL ELSE STRING_AGG(DISTINCT "{column}") END
                AS "INFO/{column}"
                """)

    # Join to create INFO column
    any_value_sql = " || ".join(any_value_list)
    # Remove last caracter ';'
    any_value_sql = f""" SUBSTR({any_value_sql}, 1, LENGTH({any_value_sql}) - 1) """
    # Join INFO columns
    explode_info_sql = " , ".join(explode_info_list)


    # Partitioning checking
    #input_file_size = os.path.getsize(input_file)
    #if input_file_size > 1024*1024*1024 or reduce_memory:
    partitioning_by_chromosomes = True
    if reduce_memory:
        partitioning = True
        #partitioning_by_chromosomes = True # "chrom"
    else:
        partitioning = False
        #partitioning_by_chromosomes = False # "chrom"

    # Query to select columns
    query_select_parquet = f"""
        SELECT
            CASE WHEN "#CHROM" LIKE 'chr%' THEN '' ELSE 'chr' END ||
            CASE WHEN "#CHROM" = 'MT' THEN 'M' WHEN "#CHROM" = '23' THEN 'X' WHEN "#CHROM" = '24' THEN 'Y' ELSE "#CHROM" END
            AS '#CHROM',
            CASE WHEN "REF" in ('-','') OR "ALT" in ('-','') THEN "POS"-1 ELSE "POS" END
            AS "POS",
            '' AS "ID",
            CASE WHEN "REF" in ('-','') THEN 'N' WHEN "ALT" in ('-','') THEN 'N' || "REF" ELSE "REF" END
            AS "REF",
            CASE WHEN "REF" in ('-','') THEN 'N' || "ALT" WHEN "ALT" in ('-','') THEN 'N' ELSE "ALT" END
            AS "ALT",
            '' AS "QUAL",
            'PASS' AS "FILTER",
            {any_value_sql} AS "INFO",
            {explode_info_sql}
        FROM annovar
        GROUP BY "#CHROM", "POS", "REF", "ALT" 
    """

    # Partitioning 
    if partitioning:

        # Log
        log.info("Partitioning...")

        if partitioning_by_chromosomes:

            # Log
            log.debug("Partitioning by chromosomes...")

            # Check list of chromosome
            query_chrom_list = """ SELECT distinct("#CHROM") FROM annovar """
            chrom_list = list(conn.query(query_chrom_list).df()["#CHROM"])

            # Query to select columns with filter by chromosome
            query_select_parquet_chrom = f"""
                SELECT
                    CASE WHEN "#CHROM" LIKE 'chr%' THEN '' ELSE 'chr' END ||
                    CASE WHEN "#CHROM" = 'MT' THEN 'M' WHEN "#CHROM" = '23' THEN 'X' WHEN "#CHROM" = '24' THEN 'Y' ELSE "#CHROM" END
                    AS '#CHROM',
                    CASE WHEN "REF" in ('-','') OR "ALT" in ('-','') THEN "POS"-1 ELSE "POS" END
                    AS "POS",
                    '' AS "ID",
                    CASE WHEN "REF" in ('-','') THEN 'N' WHEN "ALT" in ('-','') THEN 'N' || "REF" ELSE "REF" END
                    AS "REF",
                    CASE WHEN "REF" in ('-','') THEN 'N' || "ALT" WHEN "ALT" in ('-','') THEN 'N' ELSE "ALT" END
                    AS "ALT",
                    '' AS "QUAL",
                    'PASS' AS "FILTER",
                    {any_value_sql} AS "INFO",
                    {explode_info_sql}
                FROM annovar
                WHERE "#CHROM" = ?
                GROUP BY "#CHROM", "POS", "REF", "ALT" 
            """

            # Query to copy into files
            query_parquet_partition_chrom = f"""
                COPY ({query_select_parquet_chrom})
                TO '{output_file}.tmp.partition' WITH (FORMAT PARQUET, PARTITION_BY ("#CHROM"), OVERWRITE_OR_IGNORE)
                """
            
            # Copy files
            for chrom in chrom_list:
                #log.info(f"Partitioning {chrom}...")
                res = conn.execute(query_parquet_partition_chrom, [chrom])
                del res
                gc.collect()

        # else:

        #     # Log
        #     log.debug("Partitioning by chromosomes in parallel...")

        #     # Query to copy into files
        #     query_parquet_partition = f"""
        #         COPY ({query_select_parquet})
        #         TO '{output_file}.tmp.partition' WITH (FORMAT PARQUET, PARTITION_BY ("#CHROM"), OVERWRITE_OR_IGNORE)
        #     """

        #     # Copy files
        #     res = conn.query(query_parquet_partition)
        #     del res
        #     gc.collect()

    
    # Export Parquet
    if output_file_parquet:

        # Log
        log.info("Exporting Parquet...")

        # Check source
        if partitioning:
            source = f""" SELECT * FROM parquet_scan('{output_file}.tmp.partition/*/*.parquet', hive_partitioning=1) """
        else:
            source = query_select_parquet

        # Query to copy file
        query_parquet = f"""
            COPY ({source})
            TO '{output_file_parquet}' WITH (FORMAT PARQUET)
        """

        # Copy file
        res = conn.query(query_parquet)
        del res
        gc.collect()


    # Export VCF
    log.info("Exporting VCF...")

    # Check source
    if partitioning:
        source = f""" SELECT "#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO" FROM parquet_scan('{output_file}.tmp.partition/*/*.parquet', hive_partitioning=1) """
    elif output_file_parquet:
        source = f""" SELECT "#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO" FROM '{output_file_parquet}' """
    else:
        source = f""" SELECT "#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO" FROM ({query_select_parquet}) """

    # Query to copy file
    query = f"""
        COPY ({source})
        TO '{output_file}.tmp.translation.variants.vcf.gz' WITH (FORMAT CSV, COMPRESSION GZIP, DELIMITER '\t', QUOTE '')
    """
    
    # Copy file
    res = conn.query(query)
    del res
    gc.collect()

    # Cllose connexion
    conn.close()

    # BCFTools to reheader sort normalize
    if partitioning:
        log.info("VCF Sorting and Normalization by chromosomes...")
        command = f"zcat {output_file}.tmp.translation.header.vcf.gz > {output_file}.tmp.split.vcf; "
        for chrom in chrom_list_fixed:
            command += f"zcat {output_file}.tmp.translation.header.vcf.gz > {output_file}.tmp.translation.{chrom}.vcf; "
            command += f"zcat {output_file}.tmp.translation.variants.vcf.gz | grep -P '{chrom}\t' >> {output_file}.tmp.translation.{chrom}.vcf; "
            command += f"{bcftools} sort --max-mem={maxmem} {output_file}.tmp.translation.{chrom}.vcf 2>{output_file}.tmp.err | {bcftools} norm --threads={threads} --check-ref s -f {genome} 2>{output_file}.tmp.err | {bcftools} view -H 2>{output_file}.tmp.err >> {output_file}.tmp.split.vcf 2>{output_file}.tmp.err;  "
            
        command += f"{bcftools} view {output_file}.tmp.split.vcf -Oz -o {output_file} ; "
        log.debug("bcftools command: " + command)
        subprocess.run(command, shell=True)

    else:
        log.info("VCF Sorting and Normalization...")
        command = f"""zcat {output_file}.tmp.translation.header.vcf.gz {output_file}.tmp.translation.variants.vcf.gz 2>{output_file}.tmp.err | {bcftools} sort --max-mem={maxmem} 2>{output_file}.tmp.err | {bcftools} norm --threads={threads} --check-ref s -f {genome} -Oz -o {output_file} 2>{output_file}.tmp.err """
        log.debug("bcftools command: " + command)
        subprocess.run(command, shell=True)

    # Clean
    clean_command = f""" rm -rf {output_file}.tmp.* {output_file_parquet}.tmp.* """
    subprocess.run(clean_command, shell=True)

