from plugins.update_database.utils import (
    read_json,
    count_row_file,
    timeit,
    metaheader_rows,
    sort_vcf,
)
from plugins.update_database.factory import Database
from howard.tools.query import query
from howard.tools.tools import arguments, commands_arguments, shared_arguments
from howard.functions.commons import compress_file
import argparse
import sys
import json
import os
import io
import gzip
import re
from collections import OrderedDict
import pprint
import logging as log
from Bio.bgzf import BgzfWriter, BgzfReader
from tqdm import tqdm
import subprocess
from os.path import join as osj

# Globals
KEEP = [
    "nhomalt",
    "AC",
    "AN",
    "AC_male",
    "AC_asj",
    "AC_amr",
    "AC_afr",
    "AC_eas",
    "AC_fin",
    "AC_nfe",
    "AC_oth",
    "AC_sas",
    "AN_asj",
    "AN_amr",
    "AN_afr",
    "AN_eas",
    "AN_fin",
    "AN_nfe",
    "AN_oth",
    "AN_sas",
    "AC_popmax",
    "AN_popmax",
]

KEEP_GRCH38 = [
    "nhomalt_joint",
    "AC_joint",
    "AN_joint",
    "AC_joint_XY",
    "AC_joint_afr",
    "AC_joint_ami",
    "AC_joint_amr",
    "AC_joint_asj",
    "AC_joint_eas",
    "AC_joint_fin",
    "AC_joint_mid",
    "AC_joint_nfe",
    "AC_joint_sas",
    "AC_joint_remaining",
    "AN_joint_afr",
    "AN_joint_ami",
    "AN_joint_amr",
    "AN_joint_asj",
    "AN_joint_eas",
    "AN_joint_fin",
    "AN_joint_mid",
    "AN_joint_nfe",
    "AN_joint_sas",
    "AN_joint_remaining",
    "AC_grpmax_joint",
    "AN_grpmax_joint"
]


class Gnomad(Database):
    """
    class to deal with gnomad data. It's easier to download gnomad chr through gs utils inside a docker container

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
        data_folder=None,
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
        self.data_folder = data_folder

    def write_header(
        self, io_out: io.TextIOWrapper, ordered_dict: OrderedDict, keep: list
    ):
        """
        Write processed metaheader row in VCF-like

        """
        mandatory = ["FILTER", "ALT", "INFO", "FORMAT", "contig"]
        for key, val in ordered_dict.items():
            if key in mandatory:
                for info_field, meta_desc in val.items():
                    if key in ["FILTER", "ALT"]:
                        io_out.write(
                            metaheader_rows(
                                key, meta_desc["Description"], id=info_field
                            )
                            + "\n"
                        )
                    elif key in ["INFO", "FORMAT"]:
                        if info_field in keep:
                            io_out.write(
                                metaheader_rows(
                                    key,
                                    meta_desc["Description"],
                                    info_field,
                                    meta_desc["Number"],
                                    meta_desc["Type"],
                                )
                                + "\n"
                            )
                    elif key == "contig":
                        io_out.write(
                            metaheader_rows(
                                key,
                                meta_desc["assembly"],
                                number=meta_desc["length"],
                                id=info_field,
                            )
                            + "\n"
                        )
                    else:
                        io_out.write(metaheader_rows(key, val) + "\n")
            else:
                io_out.write(metaheader_rows(key, val) + "\n")

    def process_header_rows(
        self, o: io.TextIOWrapper, lines: str, header: list, keep: list
    ) -> OrderedDict:
        """
        Writing header in VCF-like format, keep only info field provided

        :param o: textiowrapper which is the output stream
        :param lines: header row from the input (vcf)
        :param keep: list of info field to split and extract
        :return: an ordered dict with keys to extract and None as value
        """
        header_explode = self.explode_header(header)
        ordered_dict = OrderedDict(
            (key, None) for key in list(header_explode["INFO"].keys()) if key in keep
        )
        # ADD chr before contig name
        new_contig = {}
        for contig, values in header_explode["contig"].items():
            new_contig[f"chr{contig}"] = values
        header_explode["contig"] = new_contig
        self.write_header(o, header_explode, list(ordered_dict.keys()))
        mandatory_cols = lines.strip().split("\t")
        o.write("\t".join(mandatory_cols) + "\n")
        return ordered_dict

    def process_variants_rows(
        self, o: io.TextIOWrapper, lines: str, ordered_dict: OrderedDict
    ):
        """
        Writing variants in VCF-like format with some info field split in separate column

        :param o: textiowrapper which is the output stream
        :param lines: variant row from the input (vcf)
        :param ordered_dict: an ordered dict with keys to extract and None as value
        """
        dico = ordered_dict.copy()
        lines = lines.strip().split("\t")
        if not lines[0].startswith("chr"):
            lines[0] = f"chr{lines[0]}"
        for item in lines[7].split(";"):
            try:
                key, val = item.split("=")
            except ValueError:
                key = item
                val = True
            if key in dico:
                dico[key] = val
        lines[7] = ";".join(
            [
                f"{key}={val}"
                for key, val in list(
                    zip(dico.keys(), list(self.replace_none_with_dot(dico).values()))
                )
            ]
        )
        o.write("\t".join(lines) + "\n")

    @timeit
    def parse_info_field(self, file: str, output: str, keep: list):
        """
        Process VCF info field

        :param file: input gzipped vcf
        :param output: gzipped output VCF-like
        :param keep: list of info field to split and extract
        """
        header = []
        number_rows = count_row_file(file)
        with BgzfReader(file, "rt") as f, BgzfWriter(output, "wt") as o:
            for i, lines in tqdm(
                enumerate(f),
                total=number_rows,
                desc="Filtering info fields",
                leave=False,
            ):
                if lines.startswith("##"):
                    header.append(lines)
                elif lines.startswith("#"):
                    ordered_dict = self.process_header_rows(o, lines, header, keep)
                else:
                    self.process_variants_rows(o, lines, ordered_dict)

    @staticmethod
    def replace_none_with_dot(ordered_dict):
        for key, value in ordered_dict.items():
            if value is None:
                ordered_dict[key] = "."
        return ordered_dict

    def explode_header(self, header: list, notparse=[]):
        """
                Parses the header of a VCF file and returns a dictionary with metaheader information.

                :param header: List of header lines.
                :param notparse: List of prefixes to avoid splitting in misformatted fields.
                :return: Dictionary containing metaheader information.

        ##FILTER=<ID=RF,Description="Failed random forest filtering thresholds of 0.2634762834546574, 0.22213813189901457 (probabilities of being a true positive variant) for SNPs, indels">
        """
        dico = {}
        error = []

        for line in header:
            line = line.strip()
            if not line.startswith("##"):
                continue

            row = line.split("##", 1)[1]

            if row.endswith(">") and not any(row.startswith(item) for item in notparse):
                key, content = row.split("=", 1)

                # Initialize key in dictionary if not already present
                if key not in dico:
                    dico[key] = {}

                match = re.search(r"<(.*)>", content)
                if match:
                    if key == "INFO":
                        fields = match.group(1).split(",", 3)
                    elif key == "contig":
                        fields = match.group(1).split(",", 2)
                    else:
                        fields = match.group(1).split(",", 1)
                    tmp = {}

                    for field in fields:
                        if "=" in field:
                            k, v = field.split("=", 1)
                            tmp[k] = v
                        else:
                            print("Probably wrong header row:", field)

                    # Check and fix Description field
                    if "Description" in tmp:
                        desc = tmp["Description"]
                        if not (desc.startswith('"') and desc.endswith('"')):
                            print(f"Error in Description field: {desc}")
                            error.append(line)

                    # Use ID as the primary key if it exists
                    if "ID" in tmp:
                        id_value = tmp.pop("ID")
                        dico[key][id_value] = tmp
                    else:
                        dico[key] = tmp

            elif not row.startswith("contig"):
                key, value = row.split("=", 1)
                dico[key] = value
        if error:
            log.warning(f"fields not parsed {' '.join(error)}")
        return dico

    def concat_info_field(self, input, output, header=None):
        mandatory = []
        with BgzfReader(input, "rt") as file, BgzfWriter(output, "wt") as out:
            with tqdm(
                total=count_row_file(input),
                desc="Processing lines",
                unit="line",
                leave=False,
            ) as pbar:  # count_row_file(input)
                if header is None:
                    for line in file:
                        line = line.strip()
                        if line.startswith("##"):
                            out.write(line + "\n")
                        elif line.startswith("#CHROM"):
                            mandatory.extend(line.strip().split())
                            out.write(
                                "\t".join(
                                    [
                                        "#CHROM",
                                        "POS",
                                        "ID",
                                        "REF",
                                        "ALT",
                                        "QUAL",
                                        "FILTER",
                                        "INFO",
                                    ]
                                )
                                + "\n"
                            )
                        else:
                            self.merge_info_columns(line, mandatory, out)
                        pbar.update(1)
                else:
                    with open(header, "r") as hdr:
                        for header_row in hdr:
                            if header_row.startswith("#CHROM"):
                                mandatory.extend(header_row.strip().split())
                                out.write(
                                    "\t".join(
                                        [
                                            "#CHROM",
                                            "POS",
                                            "ID",
                                            "REF",
                                            "ALT",
                                            "QUAL",
                                            "FILTER",
                                            "INFO",
                                        ]
                                    )
                                    + "\n"
                                )
                            else:
                                out.write(header_row)
                    for line in file:
                        line = line.strip()
                        if line.startswith("#CHROM"):
                            continue
                        else:
                            self.merge_info_columns(line, mandatory, out)
                        pbar.update(1)
        return output

    def merge_info_columns(self, line: str, mandatory: list, out: BgzfWriter):
        """
        Merge info columns from SQL to TSV request on gnomAD database
        #CHROM  POS     REF     ALT     annot1  annot2
        """
        line = line.strip().split("\t")
        values = line[4:]
        columns = mandatory[4:]
        res_string = ";".join(
            [f"{key}={val}" for key, val in list(zip(columns, values))]
        )

        variant = [line[0], line[1], ".", line[2], line[3], ".", "."]
        variant.append(res_string)
        out.write("\t".join(variant) + "\n")

    def vcf_query(self):
        #gnomad2.1
        # freq = ["amr", "afr", "asj", "eas", "fin", "nfe", "oth", "sas"]
        # coalesce = [
        #     f"COALESCE(SUM(CAST(AC_{val} AS BIGINT))/SUM(AN_{val}), 0) AS gnomadAltFreq_{val}"
        #     for val in freq
        # ]
        # gnomad_query = (
        #     'COPY (SELECT "#CHROM", POS, REF, ALT, CAST(SUM(CAST(AC AS BIGINT)) AS BIGINT) AS AC_all, CAST(SUM(AN) AS BIGINT) AS AN_all, '
        #     + ", ".join(coalesce)
        #     + ", COALESCE(SUM(CAST(REPLACE(AC_popmax, '.', '0') AS BIGINT))/SUM(CAST(REPLACE(AN_popmax, '.', '0') AS BIGINT)), 0) AS gnomadAltFreq_popmax, \
        #         CAST(SUM(CAST(nhomalt AS BIGINT)) AS BIGINT) AS gnomadHomCount_all, \
        #             COALESCE(SUM(CAST(AC AS BIGINT))/SUM(AN), 0) AS gnomadAltFreq_all, \
        #                 CAST(SUM(CAST(AC AS BIGINT)) - (2 * SUM(CAST(nhomalt AS BIGINT))) AS BIGINT) AS gnomadHetCount_all, \
        #                    CAST(SUM(CASE WHEN \"#CHROM\" = 'chrX' THEN CAST(AC_male AS BIGINT) ELSE 0 END) AS BIGINT) AS gnomadHemCount_all \
        #                         FROM parquet_scan('"
        #     + self.data_folder
        #     + "/*.parquet', union_by_name = true) GROUP BY \"#CHROM\", POS, REF, ALT) TO '"
        #     + self.data_folder
        #     + "/exomes.genomes.processed.csv' DELIMITER '\t' CSV HEADER"
        # )
        freq = ["afr", "ami", "amr", "asj", "eas", "fin", "mid", "nfe", "sas", "remaining"]
        coalesce = [
            f"COALESCE(SUM(CAST(CAST(AC_joint_{val} AS BIGINT) AS BIGINT))/SUM(CAST(AN_joint_{val} AS BIGINT)), 0) AS gnomadAltFreq_{val}"
            for val in freq
        ]
        gnomad_query = (
            'COPY (SELECT "#CHROM", POS, REF, ALT, CAST(SUM(CAST(AC_joint AS BIGINT)) AS BIGINT) AS AC_all, CAST(SUM(CAST(AN_joint AS BIGINT)) AS BIGINT) AS AN_all, '
            + ", ".join(coalesce)
            + ", ROUND(COALESCE(SUM(CAST(REPLACE(AC_grpmax_joint, '.', '0') AS BIGINT))/SUM(CAST(REPLACE(AN_grpmax_joint, '.', '0') AS BIGINT)), 0), 7) AS gnomadAltFreq_popmax, \
                ROUND(CAST(SUM(CAST(nhomalt_joint AS BIGINT)) AS BIGINT), 7) AS gnomadHomCount_all, \
                    ROUND(COALESCE(SUM(CAST(AC_joint AS BIGINT))/SUM(CAST(AN_joint AS BIGINT)), 0), 7) AS gnomadAltFreq_all, \
                        ROUND(CAST(SUM(CAST(AC_joint AS BIGINT)) - (2 * SUM(CAST(nhomalt_joint AS BIGINT))) AS BIGINT), 7) AS gnomadHetCount_all, \
                           ROUND(CAST(SUM(CASE WHEN \"#CHROM\" = 'chrX' THEN CAST(AC_joint_XY AS BIGINT) ELSE 0 END) AS BIGINT), 7) AS gnomadHemCount_all \
                                FROM parquet_scan('"
            + self.data_folder
            + "/*.parquet', union_by_name = true) GROUP BY \"#CHROM\", POS, REF, ALT) TO '"
            + self.data_folder
            + "/exomes.genomes.processed.csv' DELIMITER '\t' CSV HEADER"
        )
        commands_arguments["query"]["groups"]["Query"]["query_print_mode"] = True
        commands_arguments["query"]["groups"]["main"]["param"] = {
            "query": {"query": gnomad_query}
        }
        log.debug(gnomad_query)
        query(
            argparse.Namespace(
                command="query",
                output="",
                input=[file for file in os.listdir(self.data_folder) if file.endswith(".parquet")][0],
                arguments_dict={
                    "arguments": arguments,
                    "commands_arguments": commands_arguments,
                    "shared_arguments": shared_arguments,
                },
                param=json.dumps({"query": {"query": gnomad_query}}),
            )
        )
        return osj(self.data_folder, "exomes.genomes.processed.csv")
    
    def __call__(self, vcf_chrom):
        """
        1) Extract required annotations from vcf -> vcf
        2) Convert each vcf to parquet
        3) SQL query: merge exomes and genomes for all contig and calculate HOWARD annotations -> csv
        4) Merge back csv with header and convert it to parquet
        """
        
        log.info(f"Processing {vcf_chrom} (PID: {os.getpid()})")
        cleaned_vcf = osj(self.data_folder, vcf_chrom).replace(
                ".vcf.bgz", ".parsed.vcf.gz"
            )
        self.parse_info_field(osj(self.data_folder, vcf_chrom), cleaned_vcf, self.config_json["gnomad_fields"])
        # if not os.path.exists(
        #     cleaned_vcf.replace(".vcf.gz", ".parquet")
        # ):
        #     log.debug(
        #         f"VCF to parquet {cleaned_vcf}"
        #     )
        #     self.vcf_to_parquet(
        #         cleaned_vcf
        #     )
        return f"Done {cleaned_vcf.replace('.vcf.gz', '.parquet')}"
    
    def update_gnomad(self):
        parquet_list = []
        for parsed in os.listdir(self.data_folder):
            if parsed.endswith(".parsed.vcf.gz"):
                parquet_file = osj(self.data_folder, parsed).replace(".vcf.gz", ".parquet")
                if not os.path.exists(parquet_file):
                    log.debug(
                        f"VCF to parquet {osj(self.data_folder, parsed)}"
                    )
                    self.vcf_to_parquet(
                        osj(self.data_folder, parsed)
                    )
                parquet_list.append(parquet_file)
        # exit()
        log.info(f"Gnomad DB from {' '.join(parquet_list)}")
        if not os.path.exists(osj(self.data_folder, "exomes.genomes.processed.csv.gz")):
            query_output_file = self.vcf_query()
            compress_file(query_output_file, query_output_file + ".gz")
        header = osj(os.path.dirname(os.path.abspath(__file__)), "config", "gnomad.hdr")

        unsorted_vcf = self.concat_info_field(
            osj(self.data_folder, "exomes.genomes.processed.csv.gz"),
            osj(self.data_folder, "gnomad.vcf.gz"),
            header,
        )
        sorted_vcf = sort_vcf(
            unsorted_vcf, unsorted_vcf.replace(".vcf.gz", ".sorted.vcf.gz")
        )
        self.vcf_to_parquet(sorted_vcf)


"""
Matching annotations gnomad 4.1.0
gnomadaltfreq_amr:          AF_joint_amr
gnomadaltfreq_afr:          AF_joint_afr
gnomadaltfreq_ami:          AF_joint_ami
gnomadaltfreq_asj:          AF_joint_asj
gnomadaltfreq_eas:          AF_joint_eas
gnomadaltfreq_fin:          AF_joint_fin
gnomadaltfreq_mid:          AF_joint_mid
gnomadaltfreq_nfe:          AF_joint_nfe
gnomadaltfreq_remaining:    AF_joint_remaining
gnomadaltfreq_sas:          AF_joint_sas
gnomadaltfreq_all:          AF_joint
gnomadhomcount_all:         nhomalt_joint
gnomadhetcount_all          /
gnomadhemcount_all:         AC_joint_XY
gnomadAltFreq_popmax:       AF_grpmax_joint

Je rajoute:
gnomadAlleleCount:          AC_joint
gnomadAlleleNumber:         AN_joint

"""