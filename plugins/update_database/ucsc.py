from howard.functions.commons import download_file, compress_file
from howard.tools.convert import convert
from howard.tools.tools import (
    arguments,
    commands_arguments,
    shared_arguments,
    help_generation,
    set_log_level,
)
from utils import get_compiled_pattern

import argparse
import logging as log
import requests
from bs4 import BeautifulSoup
from os.path import join as osj
import pyBigWig
import tqdm
from utils import extract_gz_file, get_md5
import os
import sys
import gzip
import json
import time
import pathlib
from Bio.bgzf import BgzfWriter

# https://docs.bedbase.org/bedboss/tutorials/bedmaker_tutorial/
# transform bigwig
# URL of the public FTP directory (over HTTP)
# ftp_directory_url = 'https://hgdownload.cse.ucsc.edu/goldenPath/hg19/'  #Example URL
# def download_file(
#     url: str,
#     dest_file_path: str,
#     chunk_size: int = 1024 * 1024,
#     try_aria: bool = True,
#     aria_async_dns: bool = False,
#     threads: int = 1,
#     quiet: bool = True,
# ):
# Test https://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.gc5Base.wig.gz


class Databaseucsc:
    def __init__(
        self,
        link=None,
        database=None,
        exclude_link=None,
        databases_folder=None,
        verbosity=None,
        input=None,
        config_json=None,
    ):
        self.link = link
        self.database = database
        self.exclude_link = ["parentDirectory"]
        self.databases_folder = databases_folder
        self.input = input
        self.config_json = self.read_json(config_json)
        if input is not None and (
            os.path.basename(os.path.splitext(self.input)[0])
            in self.config_json["header"][self.database]
        ):
            self.subdatabase = os.path.basename(os.path.splitext(self.input)[0])
        if verbosity is None:
            verbosity = "debug"
        set_log_level(verbosity)

    @staticmethod
    def read_json(configfile: str) -> dict:
        with open(configfile) as js:
            return json.load(js)

    def list_databases(self):
        try:
            # Send an HTTP GET request to the directory URL
            response = requests.get(
                self.link, timeout=10
            )  # Adjust timeout as         necessary
            response.raise_for_status()  # Raise an error for bad responses

            # Parse the HTML content using BeautifulSoup
            soup = BeautifulSoup(response.content, "html.parser")

            # Print the parsed HTML content for debugging
            log.debug("Parsed HTML content:")
            # print(soup.prettify())

            # Find all anchor tags
            file_links = []
            for link in soup.find_all("a"):
                href = link.get("href")
                if (
                    any(href for pattern in self.exclude_link if pattern not in href)
                    and not href.startswith("?")
                    and not href.startswith("/")
                ):
                    file_links.append(href)
            log.debug("Files available for download:")
            for file in file_links:
                log.info(file)
        except requests.exceptions.RequestException as e:
            log.error(f"Error accessing the FTP server: {e}")
            exit()
        return file_links

    def download_folder(self, ends=None, starts=None):
        available = self.list_databases()
        for files in available:
            if (
                ends is not None
                and files.endswith(ends)
                or starts is not None
                and files.startswith(starts)
                and not self.check_exists(files)
            ):
                self.download(osj(self.link, files))

    def check_exists(self, file):
        if os.path.exists(osj(self.databases_folder, file)):
            log.debug(f"{file} already exists in {self.databases_folder}, skipping")
            return True
        else:
            return False

    def download(self, file):
        name = os.path.basename(file)
        download_file(
            url=file,
            dest_file_path=osj(self.databases_folder, name),
            threads=4,
            try_aria=True,
            quiet=False,
        )
        os.chmod(osj(self.databases_folder, name), 0o755)
        return osj(self.databases_folder, name)

    def create_header(self, header_dict: dict, type=str, subdatabase=None):
        header = []
        header.append("##fileformat=VCFv4.3")
        header.append("##fileDate=%s" % time.strftime("%d/%m/%Y"))
        header.append("##InputFile=" + self.input)
        if subdatabase and header_dict.get(self.subdatabase) is not None:
            header.append(
                self.metaheader_rows(
                    "INFO",
                    self.subdatabase,
                    header_dict[self.subdatabase]["Number"],
                    header_dict[self.subdatabase]["Type"],
                    header_dict[self.subdatabase]["Description"],
                )
            )
        else:
            for annotations in header_dict:
                header.append(
                    self.metaheader_rows(
                        "INFO",
                        annotations,
                        header_dict[annotations]["Number"],
                        header_dict[annotations]["Type"],
                        header_dict[annotations]["Description"],
                    )
                )
        if type == "tsv":
            header.append("#CHROM\tSTART\tEND\t" + "\t".join(list(header_dict.keys())))
        elif type == "vcf":
            header.append(
                "\t".join(
                    ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"]
                )
            )
        else:
            raise ValueError(f"Type {type} not allowed, pick between tsv or vcf EXIT")
        return header

    @staticmethod
    def metaheader_rows(fields, id, number, type, description):
        """
        ##INFO=<ID=STRAND,Number=1,Type=String,Description="Gene strand">
        """

        keys = ["ID", "Number", "Type", "Description"]
        values = list(map(str, [id, number, type, '"' + description + '"']))
        return (
            "##"
            + fields
            + "=<"
            + ",".join(["=".join(val) for val in list(zip(keys, values))])
            + ">"
        )

    def bigwig_to_tsv_batch(self):
        """ """
        for files in self.databases_folder():
            if files.endswith(".bb"):
                self.bigwig_to_tsv(files.replace(".bb", ".tsv"))

    def bigbed_to_vcf_parse_entries(self, entry: str):
        # Split position and frequencies
        position, frequency = entry.split(",", 1)

        # Extract chromosome and position
        chrom, pos = position.split(":")
        pos = int(pos)

        freq = {}
        for pattern in ["REF_AF", "ALT_AF"]:
            # Compile regex patterns using lru_cache
            if pattern == "REF_AF":
                allele_pattern = get_compiled_pattern(rf"{pattern}\((\w+)\)=(\d+\.\d+)")
            else:
                allele_pattern = get_compiled_pattern(
                    rf"{pattern}\(([\w,]+)\)=([\d\.,]+)"
                )
            # Find matches for reference and alternate alleles
            allele_match = allele_pattern.search(frequency)
            if allele_match:
                freq[pattern] = {}
                # Multiallelic
                alleles, frequencies = allele_match.groups()
                for j, mono in enumerate(alleles.split(",")):
                    freq[pattern][mono] = frequencies.split(",")[j]
        try:
            for alt, value in freq["ALT_AF"].items():
                yield [
                    chrom,
                    pos,
                    list(freq["REF_AF"].keys())[0],
                    alt,
                    f"{self.subdatabase}={value}",
                ]
        except KeyError:
            log.debug(f"{entry} does not have frequency skip")

    def bigbed_to_vcf_batch(self, type, output=None, subdatabase=False):
        for files in os.listdir(self.databases_folder):
            if files.endswith(".bb"):
                files = osj(self.databases_folder, files)
                self.input = files
                self.subdatabase = os.path.basename(os.path.splitext(files)[0])
                self.bigbed_to_vcf(type, subdatabase=subdatabase)

    def bigbed_to_vcf(self, type, output=None, subdatabase=False) -> str:
        """
        Convert a BigBed file to a bgzipped VCF file.
        howard does not accept bed just rename to tsv after bb to vcf in case of .bed extension

        Parameters:
        output (str): The path to the output gzipped VCF file.

        Only for 1 based vcf bigbed file for now
        """
        log.info(f"BigBed file: {self.input}")
        bw = pyBigWig.open(self.input)
        log.debug(f"BigBed header: {bw.header()}")
        if output is None:
            output = self.input.replace(".bb", ".vcf")
        if not output.endswith(".gz"):
            output = f"{output}.gz"
        with BgzfWriter(output, "wt") as of:
            dict_json = self.config_json.get("header").get(self.database)
            if dict_json is not None:
                header = [
                    of.write(f"{annot}\n")
                    for annot in self.create_header(dict_json, type, subdatabase)
                ]
                log.debug(f"Header contains {len(header)} rows")
            else:
                log.debug(
                    f"Header won't be generated, available annotations configuration {list(self.config_json.get('header').keys())}"
                )
            for chrom, _ in bw.chroms().items():
                intervals = bw.entries(chrom, 0, bw.chroms(chrom))
                for i, interval in enumerate(intervals):
                    window = interval[2].split("\t")
                    try:
                        variants = list(self.bigbed_to_vcf_parse_entries(window[-1]))
                        for allele in variants:
                            chrom, pos, ref, alt, info = allele
                            vcf_row = [
                                chrom,
                                str(pos),
                                window[0],
                                ref,
                                alt,
                                ".",
                                ".",
                                info,
                            ]
                            of.write("\t".join(vcf_row) + "\n")
                    except TypeError:
                        continue
        return output

    def bigwig_to_tsv(self, output=None) -> str:
        """
        Convert a BigWig file to a bgzipped TSV file.
        howard does not accept bed just rename to tsv after bw to bsv in case of .bed extension
         howard convert --input hg19.100way.phastCons.tsv.gz --output hg19.100way.phastCons.parquet --parquet_partitions '#CHROM' &

        Parameters:
        output (str): The path to the output gzipped TSV file.
        """
        log.info(f"BigWig file: {self.input}")
        bw = pyBigWig.open(self.input)
        log.debug(f"BigWig header: {bw.header()}")
        if output is None:
            output = self.input.replace(".bb", ".tsv")
        if not output.endswith(".gz"):
            output = f"{output}.gz"
        with BgzfWriter(output, "wt") as of:
            dict_json = self.config_json.get("header").get(self.database)
            if dict_json is not None:
                header = [
                    of.write(f"{annot}\n") for annot in self.create_header(dict_json)
                ]
                log.debug(f"Header contains {len(header)} rows")
            else:
                log.debug(
                    f"Header won't be generated, available annotations configuration {list(self.config_json.get('header').keys())}"
                )

            for chrom, _ in bw.chroms().items():
                intervals = bw.intervals(chrom)
                for interval in intervals:
                    of.write(
                        "{}\t{}\t{}\t{}\n".format(
                            chrom, interval[0], interval[1], interval[2]
                        )
                    )
        return output

    def vcf_to_parquet(self, file):
        output = os.path.basename(file).replace(".vcf.gz", ".parquet")
        convert(
            argparse.Namespace(
                command="convert",
                input=file,
                output=osj(self.databases_folder, output),
                arguments_dict={
                    "arguments": arguments,
                    "commands_arguments": commands_arguments,
                    "shared_arguments": shared_arguments,
                },
            )
        )

    # def check_databases(self, databases):
    #     for

    def update_clinvar(self):
        variants_file = self.download(self.link)
        md5_file = self.download(self.link + ".md5")
        # compare md5
        md5_download = open(md5_file).readline().strip().split()[0]
        if md5_download != get_md5(variants_file):
            raise ValueError("MD5 are different EXIT")
        else:
            log.info(f"Create parquet file from {variants_file}")
            self.vcf_to_parquet(variants_file)


def main():
    assembly = "hg19"
    database = "phyloP46way"
    # ucsc = Databaseucsc(
    #     f"https://hgdownload.cse.ucsc.edu/goldenPath/{assembly}/{database}",
    #     database,
    #     assembly,
    #     ["parentDirectory", "goldenPath"],
    #     ".",
    # )
    # ucsc.download("https://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.gc5Base.wig.gz")

    # ucsc = Databaseucsc(
    #     "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/clinvar_20240611.vcf.gz",
    #     "clinvar_20240611.vcf.gz",
    #     ["parentDirectory", "goldenPath"],
    #     "/home1/DB/HOWARD",
    #     "debug",
    # )
    # BIG BED TO BED WORK WELL
    # bw = sys.argv[1]
    # databases = sys.argv[2]
    # ucsc = Databaseucsc(input=bw, databases=databases)
    # ucsc.bigwig_to_bed(
    #     bw.replace(".bw", ".bed"),
    #     "/home1/data/WORK_DIR_JB/howard/plugins/update_database/config/update_databases.config_json",
    # )
    Databaseucsc(
        link="https://ftp.ncbi.nih.gov/snp/population_frequency/TrackHub/latest/hg19/",
        database="ALFA",
        input="/home1/DB/HOWARD/ALFA/hg19/ALFA_AFA.bb",
        databases_folder="/home1/DB/HOWARD/ALFA/hg19",
        config_json="/home1/data/WORK_DIR_JB/howard/plugins/update_database/config/update_databases.json",
    ).bigbed_to_vcf_batch(type="vcf", subdatabase=True)


if __name__ == "__main__":
    main()
