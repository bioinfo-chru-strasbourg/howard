#!/usr/bin/env python3

import argparse
import logging as log
from os.path import join as osj
import pandas as pd
import json
import os
import time
from howard.functions import commons

mandatory = ["#CHROM", "START", "END"]


# Aim: This script was created to transform VARANK extann file to vcf, in order to be use in HOWARD 1.0.
def create_metaheader(
    df_extann: pd.DataFrame, input: str, config: dict, extra_cols=None
) -> str:
    """
    From extann file in dataframe, create metaheader of pseudo bed file
    input: path of input extann
    config: dict
    extra_cols: list of column from refgene to keep
    """
    header_processed = []
    header_config_dict = config
    header_config_name = (
        os.path.basename(input)
        .replace(".txt.gz", "")
        .replace(".txt", "")
        .replace(".tsv.gz", "")
        .replace(".tsv", "")
    )
    log.debug(f"Available header annotations {header_config_dict.keys()}")
    columns = list(df_extann.columns)
    if extra_cols:
        columns.extend(extra_cols)
    for col in columns:
        if header_config_dict.get(header_config_name) is not None:
            # log.debug(f"config {header_config_name} key {col} found in main config")
            key_dict = header_config_dict[header_config_name]
            try:
                header_processed.append(
                    metaheader_rows(
                        "INFO",
                        col,
                        key_dict[col].get("Number"),
                        key_dict[col].get("Type"),
                        key_dict[col].get("Description"),
                    )
                )
            except TypeError:
                log.warning(
                    "Mandatory metaheader fields bad formatted, adding default field Number 1 Type String Description . "
                )
                header_processed.append(add_default_metaheader("INFO", col))
            except KeyError:
                log.warning(
                    f"Mandatory metaheader fields {col} not in header config, adding default field Number 1 Type String Description . "
                )
                header_processed.append(add_default_metaheader("INFO", col))

        else:
            log.warning(f"{header_config_name} not found in config, add default values")
            header_processed.append(add_default_metaheader("INFO", col))
    header = ["##fileformat=VCFv4.4", "##fileDate=%s" % time.strftime("%d/%m/%Y")]
    try:
        header_processed.append(
            f"##reference={header_config_dict[header_config_name].get('reference')}"
        )
    except KeyError:
        log.warning("Add default reference hg19")
        header_processed.append("##reference=hg19")
    header.extend(header_processed)

    header = "\n".join(header)
    return header


def add_default_metaheader(fields, id):
    return metaheader_rows(fields, id, "1", "String", ".")


def read_json(file: str) -> dict:
    """
    From json file to python dict
    """
    with open(file, "r") as jsonFile:
        return json.load(jsonFile)


def read_refgene(refgene: str) -> pd.DataFrame:
    if refgene.endswith(".gz"):
        df_refgene = pd.read_csv(refgene, header=0, sep="\t", compression="gzip")
    else:
        df_refgene = pd.read_csv(refgene, header=0, sep="\t")
    assert (
        len(df_refgene.columns) == 7
    ), "This refgene not comming from howard database please reformat it before"
    log.debug(df_refgene.head())
    return df_refgene


def metaheader_rows(
    fields: str, id: str, number: str, type: str, description: str
) -> str:
    """
    ##INFO=<ID=STRAND,Number=1,Type=String,Description="Gene strand">
    fields: INFO, FORMAT....
    number: 0, 1, ., ...
    type: String, Float, ....
    description: descprition of the field
    conf https://samtools.github.io/hts-specs/VCFv4.4.pdf
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


def replace_values(input_string: str, config: dict) -> str:
    for key, value in config.items():
        input_string = input_string.replace(str(key), str(value))
    return input_string


def write_extann(
    param,
    header,
    output,
    df_extann,
    df_refgene,
    extra_cols=None,
    mode=None,
    df_transcript=None,
):
    """
    Write ExtAnn into a bed like file and his hdr mate
    """

    # Variants and write to file
    headerfile = output + ".hdr"
    with open(output, "w+") as o, open(headerfile, "w+") as h:
        if extra_cols:
            mandatory.extend(extra_cols)
        mandatory.extend([col for col in df_extann.columns])

        # Header
        def write_header(header, h):
            for lines in header:
                h.write(lines)
            h.write("\n" + "\t".join(mandatory) + "\n")

        write_header(header, h)

        # Annotations + header
        write_header(header, o)
        df_extann.fillna(".", inplace=True)
        for i, rows in df_extann.iterrows():
            # for each transcript
            pos_list = get_gene_coordinate(
                df_refgene, rows, extra_cols, mode, df_transcript
            )
            if pos_list:
                for data in pos_list:
                    for key, val in rows.items():
                        if param.get("replace"):
                            data.append(replace_values(str(val), param.get("replace")))
                        else:
                            data.append(val)
                    o.write("\t".join(list(map(str, data))) + "\n")


def extann_to_info(record: pd.Series) -> str:
    """
    from pandas series (row of dataframe) create the info field of the vcf from extann data per gene
    """
    info_field_list = []
    for key, val in record.items():
        info_field_list.append(f"{key}={val}")
    return ";".join(info_field_list)


def get_longest_transcript(df: pd.DataFrame, extra_col=None) -> any:
    """
    From pandas dataframe containing one gene and many transcript and coordinate return the longest \n
    if there are many same size transcript keep the MANE
    """
    longest = {}
    match = df.groupby("transcript")
    for transcript, data in match:
        start = data.iloc[0]["START"]
        end = data.iloc[-1]["END"]
        length = end - start
        longest[transcript] = length
    longest_transcript = max(longest, key=lambda k: longest[k])
    df_longest_transcript = df.loc[df["transcript"] == longest_transcript]
    tmp = [
        df_longest_transcript.iloc[0]["#CHROM"],
        df_longest_transcript.iloc[0]["START"],
        df_longest_transcript.iloc[-1]["END"],
    ]
    if extra_col:
        tmp.extend([df_longest_transcript.iloc[0][val] for val in extra_col])
    return [tmp]


def get_all_transcript(match: pd.DataFrame, extra_col=None) -> any:
    """
    Get all from trasncript from refgene matching gene name
    """
    # Keep all transcript
    match_gb = match.groupby("transcript")
    match_list = []
    # Iterate over each transcript for this gene
    for grp, grp_data in match_gb:
        tmp = [
            grp_data.iloc[0]["#CHROM"],
            grp_data.iloc[0]["START"],
            grp_data.iloc[-1]["END"],
        ]
        if extra_col:
            tmp.extend([grp_data.iloc[0][val] for val in extra_col])
        match_list.append(tmp)
    return match_list


def get_chosen_transcript(
    match: pd.DataFrame, df_transcript: pd.DataFrame, extra_col=None
):
    """
    From a txt / tsv file with gene and transcript, it will keep only provided transcript for this gene, if gene does not match it will take the longest
    """
    gene = match["name"].unique()
    if df_transcript.empty:
        return get_longest_transcript(match, extra_col)
    if len(df_transcript.loc[df_transcript["genes"] == gene].index) != 0:
        try:
            chosen = df_transcript.loc[df_transcript["genes"] == gene].iloc[0][
                "transcript"
            ]
            df_chosen_transcript = match.loc[match["transcript"] == chosen]
            tmp = [
                df_chosen_transcript.iloc[0]["#CHROM"],
                df_chosen_transcript.iloc[0]["START"],
                df_chosen_transcript.iloc[-1]["END"],
            ]
            if extra_col:
                tmp.extend([df_chosen_transcript.iloc[0][val] for val in extra_col])
            return [tmp]
        except IndexError:
            log.debug(f"Chosen transcript for {' '.join(gene)} not found in refgene")
            return get_longest_transcript(match, extra_col)
    else:
        log.debug(f"{' '.join(gene)} not provided")
        return get_longest_transcript(match, extra_col)


def get_gene_coordinate(
    df_refgene: pd.DataFrame,
    gene_row: pd.Series,
    extra_col=None,
    mode=None,
    df_transcript=None,
) -> any:
    """
    From pandas dataframe containing refgene file, get chr start stop from each gene present in extann
    do the same process for each gene/transcript it will lead to duplicate
    df_refgene: refgene dataframe
    gene_row: pandas series of extann row
    log
    """
    match = df_refgene.loc[df_refgene["name"] == gene_row["genes"]]
    if len(match.index) > 1:
        if mode == "all":
            # Return list of list for each transcript
            return get_all_transcript(match, extra_col)
        elif mode == "longest":
            return get_longest_transcript(match, extra_col)
        elif mode == "chosen":
            return get_chosen_transcript(match, df_transcript, extra_col)
        else:
            return get_all_transcript(match, extra_col)
    else:
        log.debug(f"Can't find {gene_row['genes']} symbol in refgene database")


def from_extann(args: argparse) -> None:
    """
    This function converts an txt or tsv files containing genes-bases information\n
    From a "genes" columns which contains genes symbol it will match gene coordinates in refgene database and create a bed-like output with vcf header

    :param args: `args` is an object with several attributes representing the input parameters for the
    function. These attributes include:
    :type args: argparse
    """
    param = read_json(args.param_extann)

    log.debug(f"Replace {str(param.get('replace'))}")
    log.debug(f"Config mode_extann {param['mode_extann']}")
    log.debug(f"Config extra_col_list {' '.join(param['extra_col_list'])}")
    log.debug(f"Config transcript_extann {param['transcript_extann']}")
    # Inputs
    if isinstance(args.input_extann, str):
        input_file = args.input_extann
    else:
        input_file = args.input_extann.name
    input_file = commons.full_path(input_file)

    # Output
    if isinstance(args.output_extann, str):
        output_file = args.output_extann
    else:
        output_file = args.output_extann.name
    output_file = commons.full_path(output_file)

    # Refgene
    if args.refgene_extann:
        df_refgene = read_refgene(args.refgene_extann)
    else:
        df_refgene = read_refgene(param.get("refgene_extann"))

    if args.input_extann.endswith(".gz"):
        df_extann = pd.read_csv(
            args.input_extann, header=0, sep="\t", compression="gzip"
        )
    else:
        df_extann = pd.read_csv(args.input_extann, header=0, sep="\t")

    # Transcript
    if args.transcript_extann:
        df_transcript = pd.read_csv(args.transcript_extann, header=0, sep="\t")
    else:
        df_transcript = pd.read_csv(param["transcript_extann"], header=0, sep="\t")

    # Header
    log.info("Create metaheader")
    header = create_metaheader(
        df_extann, args.input_extann, param, param["extra_col_list"]
    )
    log.info("Writting extann ...")

    write_extann(
        param,
        header,
        args.output_extann,
        df_extann,
        df_refgene,
        param["extra_col_list"],
        param["mode_extann"],
        df_transcript,
    )
