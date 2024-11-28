from plugins.update_database.utils import read_json, count_row_file, timeit
from plugins.update_database.factory import Database
import argparse
import sys
import json
import io
import gzip
import re
from collections import OrderedDict
import pprint
import logging as log
from Bio.bgzf import BgzfWriter, BgzfReader
from tqdm import tqdm
import subprocess

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


class Gnomad(Database):
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
                            self.metaheader_rows(
                                key, meta_desc["Description"], id=info_field
                            )
                            + "\n"
                        )
                    elif key in ["INFO", "FORMAT"]:
                        if info_field in keep:
                            io_out.write(
                                self.metaheader_rows(
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
                            self.metaheader_rows(
                                key,
                                meta_desc["assembly"],
                                number=meta_desc["length"],
                                id=info_field,
                            )
                            + "\n"
                        )
                    else:
                        io_out.write(self.metaheader_rows(key, val) + "\n")
            else:
                io_out.write(self.metaheader_rows(key, val) + "\n")

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
        with BgzfReader(file, "rt") as f, BgzfWriter(output, "wt") as o:
            for i, lines in enumerate(f):
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
                    fields = match.group(1).split(",")
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
                total=count_row_file(input), desc="Processing lines", unit="line"
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

    def update_gnomad(self):
        print(self.data_folder)

    # def parseargs():
    #     parser = argparse.ArgumentParser(description="VCF preprocessor")
    #     parser.add_argument("-i", "--input", type=str, help="Compressed vcf", required=True)
    #     parser.add_argument(
    #         "-o",
    #         "--output",
    #         type=str,
    #         help="Compressed vcf-like ready to use by HOWARD",
    #         required=True,
    #     )
    #     parser.add_argument(
    #         "--keep", type=str, help="File containing field to split and extract"
    #     )
    #     parser.add_argument("--process", action="store_true")
    #     parser.add_argument(
    #         "--process_header",
    #         help="Use with process, it need to contains vcf header with mandatory",
    #     )
    #     return parser.parse_args()

    # def main():
    #     args = parseargs()
    #     print(args)
    #     print("VCF preprocessor")
    #     input_file = args.input
    #     output_file = args.output
    #     if args.process:
    #         self.concat_info_field(args.input, args.output, args.process_header)
    #         exit()
    #     if args.keep is not None:
    #         keep = open(args.keep, "r").readline().strip().split(" ")
    #         print(f"keep values from {args.keep}, count: {len(keep)}")
    #     else:
    #         keep = KEEP
    #         print(f"keep values default, count: {len(keep)}")
    #     self.parse_info_field(input_file, output_file, keep)
    """
    for F in ../genomes/*.vcf.bgz; do samp=$(echo "$(basename $F)" | cut -f 1-7 -d '.').PARSED.vcf.gz && echo $samp && python ../parse_info.py --input $F --output $samp --keep ../gnomad_fields.txt;done &

    howard vcf to parquet:
    howard convert --input gnomad.exomes.r2.1.1.sites.22.PARSED.vc.gz --output gnomad.exomes.r2.1.1.sites.22.PARSED.parquet
    SQL pour générer les gnomAD
    howard query --debug --query "COPY (SELECT \"#CHROM\", POS, REF, ALT, CAST(SUM(AC) AS BIGINT) AS AC_all, CAST(SUM(AN) AS BIGINT) AS AN_all, COALESCE(SUM(AC_amr)/SUM(AN_amr), 0) AS gnomadAltFreq_amr, COALESCE(SUM(AC_afr)/SUM(AN_afr), 0) AS gnomadAltFreq_afr, COALESCE(SUM(AC_asj)/SUM(AN_asj), 0) AS gnomadAltFreq_asj, COALESCE(SUM(AC_eas)/SUM(AN_eas), 0) AS gnomadAltFreq_eas, COALESCE(SUM(AC_fin)/SUM(AN_fin), 0) AS gnomadAltFreq_fin, COALESCE(SUM(AC_nfe)/SUM(AN_nfe), 0) AS gnomadAltFreq_nfe, COALESCE(SUM(AC_oth)/SUM(AN_oth), 0) AS gnomadAltFreq_oth, COALESCE(SUM(AC_sas)/SUM(AN_sas), 0) AS gnomadAltFreq_sas, COALESCE(SUM(CAST(REPLACE(AC_popmax, '.', '0') AS BIGINT))/SUM(CAST(REPLACE(AN_popmax, '.', '0') AS BIGINT)), 0) AS gnomadAltFreq_popmax, CAST(SUM(nhomalt) AS BIGINT) AS gnomadHomCount_all, COALESCE(SUM(AC)/SUM(AN), 0) AS gnomadAltFreq_all, CAST(SUM(AC) - (2 * SUM(nhomalt)) AS BIGINT) AS gnomadHetCount_all, CAST(SUM(AC_male) AS BIGINT) AS gnomadHemCount_all FROM parquet_scan('/enadisk/gstock/LGM/lamouche/gnomAD/*.parquet', union_by_name = true) GROUP BY \"#CHROM\", POS, REF, ALT) TO 'chr22.duckdb.csv' DELIMITER '\t' CSV HEADER"

    howard query --input gnomad.exomes.r2.1.1.sites.22.PARSED.parquet --query "SELECT TABLE_CATALOG, DATA_TYPE, COLUMN_NAME FROM INFORMATION_SCHEMA.COLUMNS WHERE TABLE_NAME = 'variants'"
    """


# if __name__ == "__main__":
#     main()
