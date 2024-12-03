from howard.functions.commons import download_file, compress_file, find, command
from howard.tools.convert import convert
from howard.tools.tools import (
    arguments,
    commands_arguments,
    shared_arguments,
    help_generation,
    set_log_level,
)
from utils import (
    get_compiled_pattern,
    get_md5,
    metaheader_rows,
    find_files,
    now,
    read_md5_file,
    recursive_chmod,
)
import tempfile
import argparse
import logging as log
import requests
from bs4 import BeautifulSoup
from os.path import join as osj
import pyBigWig
import tqdm
import os
import json
import time
import pathlib
from Bio.bgzf import BgzfWriter
import re
import yaml

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


class Database:
    """
    Parent class to update howard database

    :param link: html file link to deal with
    :param database: database from howard update_database utility (check help)
    :param exclude_link:
    """

    def __init__(
        self,
        link=None,
        database=None,
        exclude_link=None,
        databases_folder=None,
        input=None,
        config_json=None,
        current_folder=None,
    ):
        self.link = link
        self.database = database
        self.exclude_link = ["parentDirectory"]
        self.databases_folder = databases_folder
        self.input = input
        self.config_json = self.read_json(config_json)
        self.current_folder = current_folder
        if input is not None and (
            os.path.basename(os.path.splitext(self.input)[0])
            in self.config_json["header"][self.database]
        ):
            self.subdatabase = os.path.basename(os.path.splitext(self.input)[0])

    @staticmethod
    def check_md5(original, downloaded) -> bool:
        """
        Check md5 from file containing md5 on server and the one calculated locally

        :param original: file containing md5 data of a remote file
        :param downloaded: path of file locally
        :return bool:
        """
        if read_md5_file(original) == get_md5(downloaded):
            log.info("MD5 check ok")
            return True
        else:
            raise ValueError("MD5 server and local are different EXIT")

    @staticmethod
    def read_json(configfile: str) -> dict:
        with open(configfile) as js:
            return yaml.safe_load(js)

    def list_databases(self, linkpage=None) -> list:
        """
        List available file to download on a webpage

        :param linkpage: if link is not set in class arg, could list database on this link
        :return: list of available data
        """
        if linkpage is not None:
            link = linkpage
        else:
            link = self.link
        try:
            response = requests.get(link, timeout=10)
            response.raise_for_status()
            soup = BeautifulSoup(response.content, "html.parser")

            log.debug("Parsed HTML content:")
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
                log.debug(file)
        except requests.exceptions.RequestException as e:
            log.error(f"Error accessing the FTP server: {e}")
            exit()
        return file_links

    def download_folder(
        self, link_list: list, output_folder, link, starts=None, skip=None
    ):
        log.info(
            f"Download files from '{link}', prefix: '{starts}', skip link containing: '{skip}'"
        )
        for files in link_list:
            if (starts is not None and files.startswith(starts)) and skip not in files:
                self.download(osj(link, files), output_folder)

    def download(self, link, output_folder=None):
        if output_folder is None:
            output_folder = self.databases_folder
        name = os.path.basename(link)
        download_file(
            url=link,
            dest_file_path=osj(output_folder, name),
            threads=4,
            try_aria=True,
            quiet=True,
        )
        os.chmod(osj(output_folder, name), 0o755)
        return osj(output_folder, name)

    def create_header(self, header_dict: dict, type: str, subdatabase=None) -> list:
        """
         Create VCF or TSV header
        :param header_dict: from configuration json dict containing metaheader informations
        :param type: either vcf or tsv
        :param subdatabase: optionnal process subdatabase, example ALFA, ALFA_EUR
        :return list of header informations
        """
        header = []
        header.append("##fileformat=VCFv4.3")
        header.append("##fileDate=%s" % time.strftime("%d/%m/%Y"))
        header.append("##InputFile=" + self.input)
        if subdatabase and header_dict.get(self.subdatabase) is not None:
            header.append(
                metaheader_rows(
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
                    metaheader_rows(
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

    def bigwig_to_tsv_batch(self):
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
                self.bigbed_to_vcf(type, output=output, subdatabase=subdatabase)

    def bigbed_to_vcf(self, type, output=None, subdatabase=False) -> str:
        """
        Convert a BigBed file to a bgzipped VCF file.
        howard does not accept bed just rename to tsv after bb to vcf in case of .bed extension.         Only for 1 based vcf bigbed file for now

        :param output (str): The path to the output gzipped VCF file.
        :return: path of output vcf
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

        :param output (str): The path to the output gzipped TSV file.
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
                    of.write(f"{annot}\n")
                    for annot in self.create_header(dict_json, "tsv")
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
        output = file.replace(".vcf.gz", ".parquet")
        param = {"input": file, "output": output, "explode": {"explode_infos": True}}
        convert(
            argparse.Namespace(
                command="convert",
                input=file,
                output=output,
                param=json.dumps(param),
                arguments_dict={
                    "arguments": arguments,
                    "commands_arguments": commands_arguments,
                    "shared_arguments": shared_arguments,
                },
            )
        )


# def main():
#     # ucsc = Ucsc(
#     #     f"https://hgdownload.cse.ucsc.edu/goldenPath/{assembly}/{database}",
#     #     database,
#     #     assembly,
#     #     ["parentDirectory", "goldenPath"],
#     #     ".",
#     # )
#     # ucsc.download("https://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.gc5Base.wig.gz")

#     # ucsc = Ucsc(
#     #     "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/clinvar_20240611.vcf.gz",
#     #     "clinvar_20240611.vcf.gz",
#     #     ["parentDirectory", "goldenPath"],
#     #     "/home1/DB/HOWARD",
#     #     "debug",
#     # )
#     # BIG BED TO BED WORK WELL
#     # bw = sys.argv[1]
#     # databases = sys.argv[2]
#     # ucsc = Ucsc(input=bw, databases=databases)
#     # ucsc.bigwig_to_bed(
#     #     bw.replace(".bw", ".bed"),
#     #     "/home1/data/WORK_DIR_JB/howard/plugins/update_database/config/update_databases.config_json",
#     # )
#     # Ucsc(
#     #     link="https://ftp.ncbi.nih.gov/snp/population_frequency/TrackHub/latest/hg19/",
#     #     database="ALFA",
#     #     input="/home1/DB/HOWARD/ALFA/hg19/ALFA_AFA.bb",
#     #     databases_folder="/home1/DB/HOWARD/ALFA/hg19",
#     #     config_json="/home1/data/WORK_DIR_JB/howard/plugins/update_database/config/update_databases.json",
#     # ).bigbed_to_vcf_batch(type="vcf", subdatabase=True)
#     Ucsc(
#         database="clinvar",
#         databases_folder="/home1/DB/HOWARD",
#         config_json="/home1/data/WORK_DIR_JB/howard/plugins/update_database/config/update_databases.json", verbosity="info").update_clinvar()


# if __name__ == "__main__":
#     main()
