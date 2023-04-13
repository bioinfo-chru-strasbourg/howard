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
from functools import partial
import itertools


from howard.objects.variants import Variants
from howard.objects.annotation import Annotation
from howard.commons import *


# Usage
# time howard from_annovar --input=/databases/annovar/current/hg19_nci60.txt --output=/databases/annotations/current/hg19/nci60.vcf.gz --to_parquet=/databases/annotations/current/hg19/nci60.parquet --annovar-code=nci60 --genome=/databases/genomes/current/hg19.fa --config=/tool/config/config.json --threads=8 
# time howard from_annovar --input=/databases/annovar/current/hg19_cosmic70.txt --output=/databases/annotations/current/hg19/cosmic.vcf.gz --to_parquet=/databases/annotations/current/hg19/cosmic.parquet --annovar-code=cosmic --genome=/databases/genomes/current/hg19.fa --config=/tool/config/config.json --threads=8 
# time howard from_annovar --input=/databases/annovar/current/hg19_clinvar_20221231.txt --output=/databases/annotations/current/hg19/clinvar.vcf.gz --to_parquet=/databases/annotations/current/hg19/clinvar.parquet --annovar-code=clinvar --genome=/databases/genomes/current/hg19.fa --config=/tool/config/config.json --threads=8 
# time howard from_annovar --input=/databases/annovar/current/hg19_gnomad211_exome.txt --output=/databases/annotations/current/hg19/gnomad_exome.vcf.gz --to_parquet=/databases/annotations/current/hg19/gnomad_exome.parquet --annovar-code=gnomad_exome --genome=/databases/genomes/current/hg19.fa --config=/tool/config/config.json --threads=12 

# time howard from_annovar --input=/databases/annovar/current/hg19_dbnsfp42a.txt --output=/databases/annotations/current/hg19/dbnsfp.vcf.gz --to_parquet=/databases/annotations/current/hg19/dbnsfp.parquet --annovar-code=dbnsfp --genome=/databases/genomes/current/hg19.fa --config=/tool/config/config.json --threads=12
# time howard from_annovar --input=/databases/annovar/current/hg19_avsnp150.txt --output=/databases/annotations/current/hg19/avsnp150.vcf.gz --to_parquet=/databases/annotations/current/hg19/avsnp150.parquet --annovar-code=avsnp150 --genome=/databases/genomes/current/hg19.fa --config=/tool/config/config.json --threads=12

# 
# gnomad211_exome
# dbnsfp42a
# avsnp150


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

    
    log.info(f"Input Annovar database: {input_file}")
    log.info(f"Output VCF database: {output_file}")
    log.info(f"Database name: {annovar_code}")

    # Make output folder
    if not os.path.exists(output_dirname):
        os.makedirs(output_dirname)

    # first VCF


    # Annovar to VCF
    log.info(f"Annovar to VCF...")
    output_file_first_vcf = output_file + ".tmp.first.vcf.gz"
    log.debug(output_file_first_vcf)
    annovar_to_vcf(input_file=input_file, output_file=output_file_first_vcf, annotations=None, header_file=None, database_name=annovar_code, bcftools=bcftools, genome=genome_file, threads=threads, maxmem="40G", remove_annotations=["ID"])

    # Merging VCF file
    log.info(f"VCF merging...")
    merge_vcf(input_file=output_file_first_vcf, output_file=output_file)

    # VCF to Parquet
    if output_file_parquet:
        log.info(f"VCF to Parquet...")
        output_file_first_parquet = output_file + ".tmp.first.parquet"
        variants = Variants(input=output_file, output=output_file_first_parquet, load=False)
        log.debug(f"VCF to Parquet loading...")
        variants.load_data()
        log.debug(f"VCF to Parquet exporting...")
        variants.export_output()

        log.info(f"Parquet infos exploding...")
        parquet_info_explode(input_file=output_file_first_parquet, output_file=output_file_parquet, threads=threads)

    

    clean_command = f""" rm -f {output_file_first_vcf} {output_file_first_parquet} {output_file_first_parquet}.hdr """
    log.debug("clean command: " + clean_command)
    subprocess.run(clean_command, shell=True)

    log.info("End")




REMOVE_DIGIT_TABLE = str.maketrans('', '', '0123456789')


# Dictionnaire des types
TYPES = {
    "int": "Integer",
    "int64": "Integer",
    "float": "Float",
    "float64": "Float",
    "object": "String"
}


def replace_dot(val):
    """
    If the value is a dot, return np.nan, otherwise return the value
    
    :param val: The value to be replaced
    :return: the value of the variable val.
    """
    if str(val) == '.':
        return np.nan
    else:
        return val


def replace_char(x):
    """
    > If the input is a string, replace all semicolons, spaces, and equal signs with underscores
    
    The first line of the function is a docstring.  This is a special string that is used to document
    the function.  It is the first thing that is read by the `help` function.  The docstring is enclosed
    in triple quotes.  The first line of the docstring is a one sentence summary of the function.  The
    next line is a blank line.  The next line is a more detailed description of the function.  The last
    line is a description of the return value
    
    :param x: The string to be translated
    :return: the string with the characters ' ;=' replaced with '___'
    """
    if isinstance(x, str):
        translation_table = str.maketrans(' ;=', '___')
        return x.translate(translation_table)
    else:
        return x


def process_chunk(chunk, info_idxs, annotations):
    """
    It takes a chunk of the dataframe, replaces the characters in the chunk, and then reads the lines in
    the chunk
    
    :param chunk: a chunk of the dataframe
    :param info_idxs: a list of the indexes of the columns that contain the information we want to
    extract
    :param annotations: a list of strings that are the names of the columns that you want to extract
    from the file
    :return: A dataframe with the same number of rows as the original dataframe, but with the columns
    replaced by the read_line function.
    """
    chunk_processed_replace_char = chunk.applymap(replace_char)
    chunk_processed_read_line = read_line(chunk_processed_replace_char, info_idxs, annotations)
    return chunk_processed_read_line
    

def read_line(chunk_data, info_idxs, annotations):
    """
    The function takes in a chunk of data, a list of indices for the annotations, and a list of the
    annotations themselves. It then iterates through each row of the chunk, and creates a string that
    contains the chromosome, position, reference allele, alternate allele, and the annotations
    
    :param chunk_data: the dataframe of the chunk
    :param info_idxs: the index of the column in the dataframe that contains the info you want to add to
    the VCF
    :param annotations: the names of the columns in the VCF file that you want to keep
    :return: A list of strings, each string is a line of the VCF file.
    """

    info_list = []

    for index, row in chunk_data.iterrows():

        chrom = str(row["#CHROM"])
        pos = row["POS"]
        ref = row["REF"]
        alt = row["ALT"]

        if not chrom.startswith("chr"):
            chrom = "chr" + chrom
        if chrom == "chrMT":
            chrom = "chrM"
        if chrom == "chr23":
            chrom = "chrX"
        if chrom == "chr24":
            chrom = "chrY"

        ref = str(ref).translate(REMOVE_DIGIT_TABLE)
        alt = str(alt).translate(REMOVE_DIGIT_TABLE)

        if ref in ["-", ""]:
            pos -= 1
            ref = "N"
            alt = "N" + alt
        if alt in ["-", ""]:
            pos -= 1
            ref = "N" + ref
            alt = "N"         

        info = ';'.join([f'{a}={row[info_idxs[i]]}' for i, a in enumerate(annotations) if row[info_idxs[i]] not in ['', '.']])
        
        info_list.append(f'{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\t.\t{info}')

    return info_list


def annovar_to_vcf(input_file, output_file, annotations=None, header_file=None, database_name=None, bcftools="bcftools", genome="hg19.fa", threads=1, maxmem="40G", remove_annotations:list = []): #, chrom_col='#CHROM', pos_col='POS', ref_col='REF', alt_col='ALT'):
    """
    It reads a tab-separated file, converts it to a VCF file, and then normalizes it
    
    :param input_file: the path to the input file
    :param output_file: the name of the output file
    :param annotations: a list of columns to be included in the INFO field of the VCF. If not provided,
    all columns will be included
    :param header_file: the path to the header file. If you don't have one, you can leave it blank
    :param database_name: The name of the database you want to use for the INFO tags
    :param bcftools: the path to the bcftools executable, defaults to bcftools (optional)
    :param genome: the reference genome to use, defaults to hg19.fa (optional)
    :param threads: number of threads to use, defaults to 1 (optional)
    :param maxmem: The maximum amount of memory to use for sorting, defaults to 40G (optional)
    :return: the value of the variable "output_file"
    """
    
    if not os.path.exists(input_file):
        raise ValueError("No input file")
        

    log.debug("")
    log.debug("")
    log.debug("#########################")
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
    auto_format = "VCF"
    try:
        df = pd.read_csv(input_file, delimiter='\t', nrows=100000, na_values=['.'], header=nb_header_line, names=header, dtype=dtype)
        auto_format = "VCF"
    except:
        auto_format = "BED"
        log.error("format not supported")
        raise ValueError("format not supported")
        return

    
    # Recheck column type les types de chaque colonne
    for col in df.columns:
        if df[col].dtype == object:
            if pd.to_numeric(df[col], errors='coerce').notnull().all():
                df[col] = pd.to_numeric(df[col], errors='coerce').astype(float)
        else:
            df[col] = pd.to_numeric(df[col], errors='coerce')
    column_types = df.dtypes

    # replace if dtype forced
    dtype_final = {}
    for column in dict(column_types):
        if dtype.get(column, None):
            dtype_final[column] = dtype[column]
        else:
            dtype_final[column] = column_types[column]


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

    # read file
    log.info("Start reading file")

    # open file input in read only, open (temporary) vcf output in write mode
    with open(input_file, 'r') as f, bgzf.open(output_file+".tmp.translation.vcf.gz", 'wt') as vcf_file:

        # header fileformat
        vcf_file.write('##fileformat=VCFv4.3\n')

        # header INFO TAGS
        for col, type_ in vcf_types.items():
            if col not in ["#CHROM", "POS", "REF", "ALT"] + remove_annotations:
                info_tag = col
                #if col.startswith():
                info_tag.replace("-", "_").replace("+", "_")
                #log.debug(f"""##INFO=<ID={info_tag},Number=.,Type={type_},Description="{col} annotation">""")
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

        # create index of annotations
        info_idxs = [header.index(a) for a in annotations]

        # Read Process

        # chunk parameters
        chunk_idx = 0
        chunksize = 1000000
        #chunksize = 100000

        nb_lines = 0

        # chunk file into blocks (prevent memory usage)
        for chunk in pd.read_csv(f, skiprows=0, sep='\t', chunksize=chunksize, engine="c", header=nb_header_line, names=header, dtype=dtype_final, na_values=['.'], low_memory=True):

            nb_lines += len(chunk)

            # log
            #log.info(f"Reading {len(chunk)} (total {nb_lines})")
            log.info(f"Reading {nb_lines} lines...")

            # Chunk into subblocks for each threads (multithreading)
            chunk_size = len(chunk) // (threads * 1)
            chunk_groups = [chunk.iloc[i:i+chunk_size,:] for i in range(0, len(chunk), chunk_size)]

            # Run process_chunk on each subblock
            with multiprocessing.Pool(processes=threads) as pool:
                # Add parameters to function process_chunk
                chunk_processor = partial(process_chunk, info_idxs=info_idxs, annotations=annotations)
                # Run multiprocessing
                results = pool.imap(chunk_processor, chunk_groups)
                # List of lines for VCF
                concatenated_list = list(itertools.chain.from_iterable(results))
                # Write into VCF
                for variant in concatenated_list:
                    vcf_file.write((variant+"\n").encode("utf-8"))

            chunk_idx += 1

    # BCFTools to reheader sort normalize
    #command = f"""{bcftools} reheader --threads={threads} -f {genome}.fai {output_file}.tmp.translation.vcf.gz | {bcftools} sort --max-mem={maxmem} | {bcftools} norm --threads={threads} --check-ref s -f {genome} -Oz -o {output_file}"""
    log.info(f"VCF Sorting and Normalization...")
    command = f"""{bcftools} sort --max-mem={maxmem} {output_file}.tmp.translation.vcf.gz 2>{output_file}.tmp.err  | {bcftools} norm --threads={threads} --check-ref s -f {genome} -Oz -o {output_file} 2>{output_file}.tmp.err """
    log.debug("bcftools command: " + command)
    subprocess.run(command, shell=True)

    clean_command = f""" rm {output_file}.tmp.* """
    subprocess.run(clean_command, shell=True)


def merge_vcf(input_file, output_file):

    # Open input and output files
    with bgzf.open(input_file, "r") as fin, bgzf.open(output_file, "w") as fout:
        # Iterate over the lines in the input file
        last_info = {}
        last_chrom = None
        last_pos = None
        last_ref = None
        last_alt = None
        write_string = None

        line_i = 0

        for line in fin:
            # If line starts with #, write it to output
            if line.startswith("#"):
                fout.write(line)
            else:
                # Split the line into its columns
                cols = line.strip().split("\t")
                chrom, pos, _, ref, alt, _, _, info_str = cols[:8]

                #log.debug(cols)

                line_i += 1

                # Split the info string into key=value pairs
                info_pairs = info_str.split(";")

                

                # Check if the current variant is the same as the last one
                if last_chrom == chrom and last_pos == pos and last_ref == ref and last_alt == alt:
                    #log.debug("continue")
                    continue  # skip if it's the same variant as before

                else:
                    #log.debug("write")
                    
                    # Write the merged info string to output
                    write_string = f"{last_chrom}\t{last_pos}\t.\t{last_ref}\t{last_alt}\t.\tPASS\t" + ";".join([f"{key}={','.join(values)}" for key, values in last_info.items()]) + "\n"
                    if write_string:
                        fout.write(write_string)
                    
                    # Clear last info dictionary
                    last_info = {}

                # Update the last info dictionary with new values
                for info_pair in info_pairs:
                    key, value = info_pair.split("=")
                    if key in last_info:
                        #last_info[key].append(value)
                        if value not in last_info[key]:
                            last_info[key].append(value)
                    else:
                        last_info[key] = [value]

                # Save the current variant as the last variant
                last_chrom, last_pos, last_ref, last_alt = chrom, pos, ref, alt

            # if line_i > 50:
            #     break

        # Write the merged info string for the last variant to output
        fout.write(f"{last_chrom}\t{last_pos}\t.\t{last_ref}\t{last_alt}\t.\tPASS\t")
        fout.write(";".join([f"{key}={','.join(values)}" for key, values in last_info.items()]))
        fout.write("\n")


def parquet_info_explode(input_file:str, output_file:str, threads:int = 1, tmp:str = "/tmp") -> None:
    """
    > It takes a parquet file, splits it by chromosome, explodes the INFO column, and then merges the
    exploded files back together
    
    :param input_file: The input file to be exploded
    :param output_file: The name of the output file
    """
    
    # Config
    config_ro = {
        "access": "RO",
        "threads": threads
    }
    config_threads = {
        "threads": threads
    }
    param_explode = {
        "explode_infos": True,
        "export_extra_infos": True
    }

    # List of exploded parquet files
    list_of_exploded_files = []

    # Create 
    parquet_input = Variants(input=input_file, config=config_ro)
    parquet_chromosomes = parquet_input.get_query_to_df(query=f""" SELECT "#CHROM" FROM '{input_file}' WHERE "#CHROM" NOT NULL GROUP BY "#CHROM" """)
    for chrom in parquet_chromosomes["#CHROM"]:
        
        # Split Chrom
        query_chrom = f""" SELECT * FROM '{input_file}' WHERE "#CHROM"='{chrom}'  """
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

    df = parquet_input.get_query_to_df(f"SELECT * FROM '{query_explode_output}' LIMIT 10 ")
    log.debug(df)
