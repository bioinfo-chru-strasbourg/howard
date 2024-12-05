"""
This should work, tell me if it works see you  I sended to you an invitation to access my git repository

#Download CADD from https://cadd.gs.washington.edu/download
#use your own script to transform TSV to VCF (TSVtoVCF.py)

# you will get one vcf then split by chromosomes (on my repository script on my github)
python split_vcf_chr.py inputvcf outputfolder threads

#convert this folder of vcf per chromosome to parquet partition
#From my github script repo vcftoparquet.py (change only input and output folder)
python vcftoparquet.py folderwithvcf outputpathofpartitioncadd.partition.parquet 1000000 10 log &

#it will create a process.log in the outputfolder

#finally you can annotate
howard annotation --input vcffile --output annotated.vcffile --annotations hg19_whole_genome_SNVschr.parquet
"""

from Bio.bgzf import BgzfReader
import pysam
from howard.tools.tools import command
import os
import logging as log
import pandas as pd
import subprocess
import multiprocessing as mp
from os.path import join as osj

CHUNK_SIZE = 1000000
JOBS = 10


def cadd_tsv_to_vcf(file):
    log.info(f"TSV to VCF {file.replace('.tsv.gz', '.vcf')}, {CHUNK_SIZE} chunksize")
    # Define the VCF columns
    vcf_columns = [
        "#CHROM",
        "POS",
        "ID",
        "REF",
        "ALT",
        "QUAL",
        "FILTER",
        "INFO",
        "FORMAT",
        "hg19_whole_genome_InDels",
    ]

    # Open the output file
    with open(file.replace(".tsv.gz", ".vcf"), "wt") as f:
        # Write the header
        with open(
            osj(os.path.dirname(os.path.abspath(__file__)), "config", "gnomad.hdr"),
            "rt",
        ) as header_file:
            f.write(header_file.read())

        # Process the TSV file in chunks
        first_chunk = True
        for chunk in pd.read_csv(
            file, sep="\t", comment="#", header=None, chunksize=CHUNK_SIZE
        ):
            # Rename the columns
            chunk.columns = ["Chrom", "Pos", "Ref", "Alt", "RawScore", "PHRED"]
            chunk.rename(
                columns={
                    "Chrom": "#CHROM",
                    "Pos": "POS",
                    "Ref": "REF",
                    "Alt": "ALT",
                    "RawScore": "CADDR",
                    "PHRED": "CADDP",
                },
                inplace=True,
            )

            # Add necessary columns for VCF format
            chunk["ID"] = "."
            chunk["QUAL"] = "."
            chunk["FILTER"] = "."
            chunk["INFO"] = chunk.apply(
                lambda row: f"CADDR={row.CADDR};CADDP={row.CADDP}", axis=1
            )
            chunk["FORMAT"] = "GT"
            chunk["hg19_whole_genome_InDels"] = "0/1"

            # Reorder columns to match VCF format
            chunk = chunk[vcf_columns]

            # Write the chunk to the VCF file
            if first_chunk:
                chunk.to_csv(f, sep="\t", index=False, header=True, mode="a")
                first_chunk = False
            else:
                chunk.to_csv(f, sep="\t", index=False, header=False, mode="a")

    # Compress with bgzip
    command("bgzip", file.replace(".tsv.gz", ".vcf"))
    return file.replace(".tsv.gz", ".vcf.gz")


def split_vcf_by_chromosome(input_file, output_folder):
    # Create output folder if it doesn't exist
    os.makedirs(output_folder, exist_ok=True)

    # Get unique chromosome values using bcftools query
    chroms = []
    with BgzfReader(input_file, "rt") as inp:
        for lines in inp:
            if lines.startswith("##contig"):
                chroms.append(lines.split(",", 1)[0].split("=")[-1])
            if lines.startswith("#CHROM"):
                break
    log.debug("Splitting vcf by chrom ...")
    for chr in chroms:
        command(
            f"bcftools view -r {chr} -O z -o {osj(output_folder, chr+'.vcf.gz')} {input_file}"
        )
    log.debug(chroms)
    for chromosome in chroms:
        output_file = osj(output_folder, f"{chromosome}.vcf.gz")
        log.info(f"Creating VCF chr{chromosome}")
        # Use bcftools view to extract variants for each chromosome
        command(
            f"bcftools view -o {output_file} -O z -r {chromosome} --threads 4 {input_file}"
        )
        command(f"tabix -p vcf {output_file}")
    return output_folder


def extract_and_write_vcf_header(input_vcf_file, output_header_file):
    with pysam.VariantFile(input_vcf_file, "r") as vcf_reader:
        header = str(vcf_reader.header)
    log.info("Writing parquet database header")
    with open(output_header_file, "w") as output_file:
        output_file.write(header)


def create_vcf_chunks(input_gz_file, chunk_size, output_prefix, output_folder):
    chrom = os.path.basename(output_folder).split("=")[-1]
    pid = os.getpid()
    log.info(f"Worker PID: {pid} CHROM: {chrom}")

    def write_chunk(chunk_count, current_chunk):
        output_file = osj(output_folder, f"{output_prefix}_chunk_{chunk_count}.vcf")
        output_parquet = output_file.replace(".vcf", ".parquet")
        with open(output_file, "w") as output:
            output.write(header)
            output.write("".join(current_chunk))
        log.info(f"Worker PID: {pid} CHROM: {chrom} Launch HOWARD chunk {chunk_count}")
        command(
            f"howard convert --threads 2 --input {output_file} --output {output_parquet} --explode_infos --verbosity 'ERROR'"
        )
        return output_file

    with pysam.VariantFile(input_gz_file, "r") as vcf_reader:
        header = str(vcf_reader.header)
        current_chunk = []
        chunk_count = 1

        chunk_size = int(chunk_size)

        # Iterate over VCF records
        for record in vcf_reader:
            current_chunk.append(str(record))

            # Check if the current chunk size exceeds the specified limit
            if len(current_chunk) >= chunk_size:
                write_chunk(chunk_count, current_chunk)
                current_chunk = []
                chunk_count += 1

        # Write any remaining records to the last chunk
        if current_chunk:
            output_file = write_chunk(chunk_count, current_chunk)
        # Cleaning
        for files in os.listdir(output_folder):
            if files.endswith(".vcf") or files.endswith(".hdr"):
                os.remove(osj(output_folder, files))

        log.info(f"Worker PID: {pid} Done")
        return len([files for files in os.listdir(output_folder)])


def update_cadd(cadd_input: str, input_folder: str, output_folder: str) -> list:
    """
    Main processing step, iterate over vcf dolder splitted by chrom and create parquet folder partitioned
    """
    output_prefix = "VCF"

    # print(cadd_input, input_folder, output_folder)
    # exit()
    if not os.path.exists(cadd_input.replace(".tsv.gz", ".chr.vcf.gz")):
        vcf = cadd_tsv_to_vcf(cadd_input)
    else:
        vcf = cadd_input.replace(".tsv.gz", ".chr.vcf.gz")
        log.debug("CADD input already converted to vcf")
    if not os.path.exists(input_folder):
        split_vcf_by_chromosome(vcf, input_folder)
    else:
        log.debug("VCF cadd already splitted")
    input_args = []
    log.info(f"Generation of CADD partiton database: {output_folder}")
    if not os.path.exists(output_folder):
        os.mkdir(output_folder)
    for vcf in os.listdir(input_folder):
        if vcf.endswith(".vcf.gz"):
            chrom = vcf.split(".")[0].split("_")[-1]
            chrom_folder = osj(output_folder, f"#CHROM={chrom}")
            if not os.path.exists(chrom_folder):
                log.info(f"Create {os.path.basename(chrom_folder)}")
                os.mkdir(chrom_folder)
            else:
                log.warning(f"{os.path.basename(chrom_folder)} already exists PASS")
                continue
            input_args.append(
                [osj(input_folder, vcf), CHUNK_SIZE, output_prefix, chrom_folder, log]
            )

            # Create header file for database
            if not os.path.exists(output_folder + ".hdr"):
                extract_and_write_vcf_header(
                    osj(input_folder, vcf), output_folder + ".hdr", log
                )
    return input_args
