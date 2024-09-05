import os
import time
from utils import metaheader_rows


def create_vcf_header(input_path) -> list:
    header = []
    header.append("##fileformat=VCFv4.3")
    header.append("##fileDate=%s" % time.strftime("%d/%m/%Y"))
    header.append("##InputFile=%s" % os.path.abspath(input_path))
    for fields in ["MISTIC_score", "MISTIC_pred"]:
        if "score" in fields:
            header.append(
                metaheader_rows(
                    "INFO",
                    fields,
                    1,
                    "Float",
                    "MISTIC score (0:1), score close than 1 show a pathogenic missense",
                )
            )
        else:
            header.append(
                metaheader_rows(
                    "INFO",
                    fields,
                    1,
                    "String",
                    "MISTIC prediction either B for Benign or D for damaging",
                )
            )
    header.append(
        "\t".join(["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"])
    )
    return header


def mistic(file: str, output: str) -> str:
    """
    Process mistic raw file to create vcf
    """
    header = create_vcf_header(file)
    write_header(header, output)
    with open(output, "a+") as out:
        with open(file, "r") as file:
            for lines in file:
                if not lines.startswith("#"):
                    lines = lines.strip().split()
                    if not lines[0].startswith("chr"):
                        lines[0] = f"chr{lines[0]}"
                    info = f"MISTIC_score={lines[4]};MISTIC_pred={lines[5]}"
                    lines = lines[0:4]
                    lines.insert(2, ".")
                    lines.extend([".", ".", info])
                    out.write("\t".join(lines) + "\n")


def write_header(header: list, output: str):
    with open(output, "w+") as out:
        for head in header:
            out.write(head + "\n")
    return output


if __name__ == "__main__":
    file = "/home1/DB/HOWARD/MISTIC/current/hg19/MISTIC_GRCh37.tsv"
    output = "/home1/DB/HOWARD/MISTIC/current/hg19/MISTIC_GRCh37.vcf"
    mistic(file, output)
