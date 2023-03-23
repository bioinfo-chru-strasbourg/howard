#!/usr/bin/env python

import argparse
from functools import partial
import itertools
import multiprocessing
import os
import subprocess
import pyarrow.parquet as pq
import duckdb
import pandas as pd
import Bio.bgzf as bgzf
import numpy as np
import concurrent.futures
from multiprocessing import Pool, cpu_count
import dask.dataframe as dd



# Usage
# python3.9 howard/tools/annovar_to_vcf.py --input=/mnt/STARK/databases/annovar/current/hg19_avsnp150.txt --output=/mnt/STARK/tmp/avsnp.vcf.gz --database_name=avsnp --format=vcf --genome=/mnt/STARK/databases/genomes/current/hg19.fa --threads=8 --maxmem=40G

# for annovardb in /mnt/STARK/databases/annovar/current/hg19_*.txt; do \
#  db_name=$(basename $annovardb | sed 's/.txt$//gi' | sed 's/hg19_//gi'); \
#  if echo "$db_name" | grep -v -q "sites"; then \
#   python3.9 howard/tools/annovar_to_vcf.py --input=$annovardb --output=/mnt/STARK/databases/annotations_annovar/$db_name.vcf.gz --database_name=$db_name --format=vcf --genome=/mnt/STARK/databases/genomes/current/hg19.fa --threads=8 --maxmem=40G; \
#  fi; \
# done;
# time python3.9 howard/tools/annovar_to_vcf.py --format=vcf --input=/mnt/STARK/databases/annovar/current/hg19_dbnsfp42a.txt --output=/mnt/STARK/tmp/dbnsfp42a.1000000chunk.vcf.gz --genome=/mnt/STARK/databases/genomes/current/hg19.fa --threads=8 --database_name=dbnsfp
# time python3.9 howard/tools/annovar_to_vcf.py --format=vcf --input=/mnt/STARK/databases/annovar/current/hg19_nci60.txt --output=/mnt/STARK/tmp/test_dask_nci60.vcf.gz --genome=/mnt/STARK/databases/genomes/current/hg19.fa --threads=8 --database_name=nci60
# time python3.9 howard/tools/annovar_to_vcf.py --format=vcf --input=/mnt/STARK/databases/annovar/current/hg19_avsnp150.txt --output=/mnt/STARK/tmp/test_dask_avsnp150.vcf.gz --genome=/mnt/STARK/databases/genomes/current/hg19.fa --threads=8 --database_name=avsnp150
# time python3.9 howard/tools/annovar_to_vcf.py --format=vcf --input=/mnt/STARK/databases/annovar/current/hg19_clinvar_20210123.txt --output=/mnt/STARK/tmp/clinvar_20210123.vcf.gz --genome=/mnt/STARK/databases/genomes/current/hg19.fa --threads=8 --database_name=clinvar_20210123



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



def tsv_to_vcf(input_file, output_file, annotations=None, header_file=None, database_name=None, bcftools="bcftools", genome="hg19.fa", threads=1, maxmem="40G"): #, chrom_col='#CHROM', pos_col='POS', ref_col='REF', alt_col='ALT'):
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
        

    print("")
    print("")
    print("#########################")
    print("input_file: " + input_file)


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
                    print("Found header")
                    header = line.strip().split('\t')
                    if "#CHROM" not in header:
                        chrom_list_for_header = ["#Chr"]
                        for chr in chrom_list_for_header:
                            if chr in header:
                                header[header.index(chr)] = "#CHROM"
                                print(f"found '{chr}' as chromosome column")
                        if "#CHROM" not in header:
                            print(f"default chromosome column is 0")
                            header[0] = "#CHROM"
                    if "POS" not in header:
                        pos_list_for_header = ["Start"]
                        for pos in pos_list_for_header:
                            if pos in header:
                                header[header.index(pos)] = "POS"
                                print(f"found '{pos}' as position column")
                        if "POS" not in header:
                            print(f"default position column is 1")
                            header[1] = "POS"
                    if "REF" not in header:
                        ref_list_for_header = ["Ref"]
                        for ref in ref_list_for_header:
                            if ref in header:
                                header[header.index(ref)] = "REF"
                                print(f"found '{ref}' as reference column")
                        if "REF" not in header:
                            print(f"default reference column is 3")
                            header[3] = "REF"
                    if "ALT" not in header:
                        alt_list_for_header = ["Alt"]
                        for alt in alt_list_for_header:
                            if alt in header:
                                header[header.index(alt)] = "ALT"
                                print(f"found '{alt}' as alternative column")
                        if "ALT" not in header:
                            print(f"default alternative column is 4")
                            header[4] = "ALT"
                    break
                else:
                    print("NO Header")
                    header = line.strip().split('\t')
                    if len(header) >= 5:
                        print(header)
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
        print("format not supported")
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
        df = pd.read_csv(genome_index, delimiter='\t', header=0,  names = ["#CHROM", "Start", "End", "col3", "col4"], dtype = {"#CHROM": str, "Start": int, "End": int})
        genome_chrom = df['#CHROM'].unique()

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
        print("List of chromosomes: " + str(chrom_list_fixed))
    #return

    # read file
    print("Start reading file")

    # open file input in read only, open (temporary) vcf output in write mode
    with open(input_file, 'r') as f, bgzf.open(output_file+".translation.vcf.gz", 'wt') as vcf_file:

        # header fileformat
        vcf_file.write('##fileformat=VCFv4.3\n')

        # header INFO TAGS
        for col, type_ in vcf_types.items():
            if col not in ["#CHROM", "POS", "REF", "ALT"]:
                info_tag = col
                #if col.startswith():
                info_tag.replace("-", "_").replace("+", "_")
                #print(f"""##INFO=<ID={info_tag},Number=.,Type={type_},Description="{col} annotation">""")
                vcf_file.write(f"""##INFO=<ID={info_tag},Number=.,Type={type_},Description="{info_tag} annotation">\n""")

        # header contig
        for chrom in chrom_list_fixed:
            #print(f"""##contig=<ID={chrom}>""")
            vcf_file.write(f"""##contig=<ID={chrom}>\n""")

        # header #CHROM line
        vcf_file.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n')

        # Determine annotations if not annotations in input
        annotations_auto = []
        if not annotations:
            for h in header:
                if h not in ["#CHROM", "POS", "REF", "ALT", "-", ""] and not str(h).startswith("-"):
                    #print(h)
                    annotations_auto.append(h)
            annotations = annotations_auto

        # create index of annotations
        info_idxs = [header.index(a) for a in annotations]


        # Read Process

        # chunk parameters
        chunk_idx = 0
        chunksize = 1000000

        # chunk file into blocks (prevent memory usage)
        for chunk in pd.read_csv(f, skiprows=0, sep='\t', chunksize=chunksize, engine="c", header=nb_header_line, names=header, dtype=dtype_final, na_values=['.'], low_memory=True):

            # log
            print(chunk)

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
                    vcf_file.write(variant+"\n")

            chunk_idx += 1

    # BCFTools to reheader sort normalize
    command = f"""{bcftools} reheader --threads={threads} -f {genome}.fai {output_file}.translation.vcf.gz | {bcftools} sort --max-mem={maxmem} | {bcftools} norm --threads={threads} --check-ref s -f {genome} -Oz -o {output_file}"""
    print("bcftools command: " + command)
    subprocess.run(command, shell=True)


def main():
    parser = argparse.ArgumentParser(description='Transform a TSV file into a VCF or BED file')
    parser.add_argument('--input_file', type=str, help='Path to the input TSV file')
    parser.add_argument('--output_file', type=str, help='Path to the output file')
    parser.add_argument('--header_file', type=str, help='Path to the header file')
    parser.add_argument('--database_name', type=str, help='Path to the header file', default=None)
    parser.add_argument('--format', type=str, default='vcf',
                        help='Output file format (default: vcf)')
    parser.add_argument('--bcftools', type=str, default="bcftools",
                        help='bcftools bin')
    parser.add_argument('--genome', type=str, default="hg19.fa",
                        help='genome path')
    parser.add_argument('--threads', type=int, default=1,
                        help='number of threads')
    parser.add_argument('--maxmem', type=str, default="40G",
                        help='max memory')

    args = parser.parse_args()

    
    if args.format in ["vcf"]:
        tsv_to_vcf(args.input_file, args.output_file, annotations=None, header_file=args.header_file, database_name=args.database_name, bcftools=args.bcftools, genome=args.genome, threads=args.threads, maxmem=args.maxmem)
    else:
        print("format not yet supported")


if __name__ == '__main__':
    main()


