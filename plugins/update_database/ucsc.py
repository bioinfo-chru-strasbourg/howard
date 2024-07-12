from howard.functions.commons import download_file
from howard.tools.convert import convert
from howard.tools.tools import (
    arguments,
    commands_arguments,
    shared_arguments,
    help_generation,
    set_log_level,
)

import argparse
import logging as log
import requests
from bs4 import BeautifulSoup
from os.path import join as osj
import pyBigWig
import tqdm
from utils import extract_gz_file, get_md5
import os

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
        link: str,
        databases: list,
        exclude_link: list,
        databases_folder: str,
        verbosity: str,
    ):
        self.link = link
        self.databases = databases
        # self.assembly = assembly
        self.exclude_link = exclude_link
        self.databases_folder = databases_folder
        set_log_level(verbosity)

    def list_databases(self, link):
        try:
            # Send an HTTP GET request to the directory URL
            response = requests.get(
                link, timeout=10
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
                print(file)

        except requests.exceptions.RequestException as e:
            log.error(f"Error accessing the FTP server: {e}")
            exit()

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

    def bigwig_to_bed(self, bigwigfile):
        bw = pyBigWig.open(bigwigfile)
        print(bw.header())
        # check number of scores/annotations columns
        # bigwig_columns_number = bw.values("chr1", 0, 1)[0]
        # print(bigwig_columns_number)
        of = open("regions.bed", "w")  # Change me

        for chrom, len in bw.chroms().items():
            intervals = bw.intervals(chrom)
            for interval in intervals:
                of.write(
                    "{}\t{}\t{}\t{}\n".format(
                        chrom, interval[0], interval[1], interval[2]
                    )
                )
        bw.close()
        of.close()

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
    # ucsc = Databaseucsc(f"https://hgdownload.cse.ucsc.edu/goldenPath/{assembly}/{database}", database, assembly, ["parentDirectory", "goldenPath"], ".")
    # ucsc.download("https://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.gc5Base.wig.gz")

    ucsc = Databaseucsc(
        "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/clinvar_20240611.vcf.gz",
        "clinvar_20240611.vcf.gz",
        ["parentDirectory", "goldenPath"],
        "/home1/DB/HOWARD",
        "debug",
    )
    ucsc.update_clinvar()


if __name__ == "__main__":
    main()
