from howard.functions.commons import compress_file, command
import logging as log
from os.path import join as osj
from plugins.update_database.factory import Database
import pandas as pd
import gzip
import re
from Bio.bgzf import BgzfWriter, BgzfReader
import subprocess
import os


class Omim(Database):
    def __init__(
        self,
        link=None,
        database=None,
        exclude_link=None,
        databases_folder=None,
        input=None,
        config_json=None,
        current_folder=None,
        refseq=None,
    ):
        super().__init__(
            link,
            database,
            exclude_link,
            databases_folder,
            input,
            config_json,
            current_folder,
        )
        self.header = self.create_header(
            self.config_json.get("header").get(self.database), "tsv"
        )
        self.refseq = refseq
        self.refgene = self.get_refgene_list()

    def write_omim(self, df_gb_sorted, output):
        # output = "/home1/BAS/lamouchj/OMIMtest022025_V2.bed.gz"
        with BgzfWriter(output, "wt") as o:
            for lines in self.header:
                o.write(lines)
            df_gb_sorted.to_csv(o, header=True, index=False, sep="\t", mode="a")

    def get_refgene_list(self):
        # there are non coding transcript inside so more than 20K gene name is  normal
        # refseq = "/home1/HUB/DB/refseq/current/hg19/ncbiRefSeq.all-except-ncc.bed"

        if not os.path.exists(self.refseq):
            raise ValueError(f"{self.refseq} does not exist")
        cmd = f"awk '{{print $4}}' {self.refseq} | sort -u"
        refgene = (
            subprocess.run(cmd, shell=True, capture_output=True, text=True)
            .stdout.strip()
            .split("\n")
        )
        log.debug(f"Refseq: {self.refseq}, {len(refgene)} genes (non-coding included)")
        return refgene

    @staticmethod
    def chrom_sort_key(chrom):
        chrom = chrom.replace("chr", "")  # Remove "chr" prefix
        if chrom.isdigit():  # If it's a number, convert to int
            return (int(chrom), "")
        else:  # Otherwise, sort special chromosomes last
            return (float("inf"), chrom)

    @staticmethod
    def agg_specific(x):
        try:
            return list(set("_".join(x)))
        except TypeError:
            return "_".join(map(str, list(set(x))))

    def update_omim(self):
        """
        Take OMIMannotations.bed.gz and remove genes that are not in refseq for alias
        """

        log.debug(self.header)
        df = pd.read_csv(
            self.input,
            skiprows=len(self.header) - 1,
            header=0,
            compression="gzip",
            sep="\t",
        )
        log.debug(f"Original df: {len(df.index)} rows")
        df_gb = df.groupby("transcript", as_index=False)
        df_list = []
        for i, df_tmp in df_gb:
            gene_to_remove = []
            if len(df_tmp["OMIM_ID"].unique().tolist()) > 1:
                gene_list = df_tmp["genes"].unique().tolist()
                # print(gene_list)
                for gene in gene_list:
                    if gene not in self.refgene:
                        gene_to_remove.append(gene)
            # print(gene_to_remove)
            df_final = df_tmp.loc[~(df_tmp["genes"].isin(gene_to_remove))]
            df_final = df_final.groupby("transcript", as_index=False).agg(
                {
                    "genes": lambda x: ",".join(list(set(x))),
                    "#CHROM": "first",
                    "START": "first",
                    "END": "first",
                    "OMIM_phenotype": lambda x: "_".join(list(set(x))),
                    "OMIM_inheritance": lambda x: "_".join(list(set(x))),
                    "OMIM_ID": lambda x: self.agg_specific(x),
                }
            )
            df_list.append(df_final)

        df_gb = pd.concat(df_list, ignore_index=True)

        df_gb = df_gb[
            [
                "#CHROM",
                "START",
                "END",
                "transcript",
                "genes",
                "OMIM_phenotype",
                "OMIM_inheritance",
                "OMIM_ID",
            ]
        ]
        df_gb_sorted = df_gb.sort_values(
            by=["#CHROM", "START", "END"],
            key=lambda x: x.map(self.chrom_sort_key) if x.name == "#CHROM" else x,
        )
        log.debug(f"OMIM sorted: {len(df_gb_sorted)} rows")
        output = osj(self.databases_folder, "OMIMannotations.final.bed.gz")
        with BgzfWriter(output, "wt") as o:
            for lines in self.header[:-1]:
                o.write(lines + "\n")
            df_gb_sorted.to_csv(o, header=True, index=False, sep="\t", mode="a")
        with open(output + ".hdr", "w+") as h:
            for lines in self.header:
                h.write(lines + "\n")
