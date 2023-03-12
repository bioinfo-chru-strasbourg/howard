import gzip
import io
import argparse
from Bio import bgzf
import sys
import os

def merge_vcf(input_file, output_file):
    # Open input and output files
    with bgzf.open(input_file, "r") as fin, bgzf.open(output_file, "w") as fout:
        # Iterate over the lines in the input file
        last_info = {}
        last_chrom = None
        last_pos = None
        last_ref = None
        last_alt = None

        for line in fin:
            # If line starts with #, write it to output
            if line.startswith("#"):
                fout.write(line)
            else:
                # Split the line into its columns
                cols = line.strip().split("\t")
                chrom, pos, _, ref, alt, _, _, info_str = cols[:8]

                # Split the info string into key=value pairs
                info_pairs = info_str.split(";")

                # Update the last info dictionary with new values
                for info_pair in info_pairs:
                    key, value = info_pair.split("=")
                    if key in last_info:
                        #last_info[key].append(value)
                        if value not in last_info[key]:
                            last_info[key].append(value)
                    else:
                        last_info[key] = [value]

                # Check if the current variant is the same as the last one
                if last_chrom == chrom and last_pos == pos and last_ref == ref and last_alt == alt:
                    continue  # skip if it's the same variant as before

                if last_chrom:
                    # Write the merged info string to output
                    fout.write(f"{last_chrom}\t{last_pos}\t.\t{last_ref}\t{last_alt}\t.\tPASS\t")
                    fout.write(";".join([f"{key}={','.join(values)}" for key, values in last_info.items()]))
                    fout.write("\n")

                # Clear last info dictionary
                last_info = {}

                # Save the current variant as the last variant
                last_chrom, last_pos, last_ref, last_alt = chrom, pos, ref, alt

        # Write the merged info string for the last variant to output
        fout.write(f"{last_chrom}\t{last_pos}\t.\t{last_ref}\t{last_alt}\t.\tPASS\t")
        fout.write(";".join([f"{key}={','.join(values)}" for key, values in last_info.items()]))
        fout.write("\n")

if __name__ == "__main__":
    # Get input and output file names from command-line arguments
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    # Call the merge_vcf function with the input and output file names
    merge_vcf(input_file, output_file)