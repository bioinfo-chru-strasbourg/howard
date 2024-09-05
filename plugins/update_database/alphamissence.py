import os
import sys

mandatory = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"]

def process(header_row: str, row:str) -> str:
    row = row.strip().split()
    index_info = {}
    info_values = {}
    for j, val in enumerate(header_row):
        if val not in mandatory:
            index_info[j] = val
            info_values[val] = None
    for i, field in enumerate(row):
        if i in index_info.keys():
            info_values[index_info[i]] = field
    info_string = ";".join([f"{k}={v}" for k, v in zip(list(info_values.keys()), list(info_values.values()))])
    variants = []
    for values in mandatory:
        if values in header_row:
            variants.append(row[header_row.index(values)])
        elif values == "INFO":
            variants.append(info_string)
        else:
            variants.append(".")
    return "\t".join(variants)

def main():
    vcf = sys.argv[1]
    header_file = sys.argv[2]
    output = sys.argv[3]
    header = []
    
    with open(output, "w+") as out, open(header_file, "r") as head, open(vcf, "r") as inp:
        for lines in head:
            if lines.startswith("#CHROM"):
                header.append(lines)
            else:
                out.write(lines)
        out.write("\t".join(mandatory)+"\n")
        header_row = header[0].strip().split()
        for lines in inp:
            if not lines.startswith("#CHROM"):
                out.write(process(header_row, lines)+"\n")


if __name__ == "__main__":
    main()