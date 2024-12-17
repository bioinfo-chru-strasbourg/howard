#!/usr/bin/env python3

import argparse
import logging as log
import pandas as pd
import os
import time
import yaml
import subprocess

from howard.functions import commons

def run_shell(command):
    try:
        # Run the command without check=True to handle return codes manually
        result = subprocess.run(
            command,
            shell=True,
            text=True,
            stdout=subprocess.PIPE,  # Capture stdout
            stderr=subprocess.PIPE,  # Capture stderr
        )

        if result.returncode == 0:
            # Command succeeded and returned output
            if result.stdout.strip():
                return result.stdout
            else:
                return ""
        elif result.returncode == 1 and "grep" in command:
            # Handle grep's no match scenario gracefully
            return ""
        else:
            # Handle other non-zero exit codes as errors
            log.error("Error occurred!")
            log.error(f"Return Code:, {result.returncode}")
            log.error(f"Error Output: {result.stderr}")
    except Exception as e:
        log.error(f"Unexpected error: {str(e)}")


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
        return yaml.safe_load(jsonFile)


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
    df_transcript=None
):
    """
    Write ExtAnn into a bed like file and his hdr mate
    """
    if os.path.exists(param.get("hgnc_file")):
        alias = pd.read_csv(param.get("hgnc_file"), sep="\t", header=0, low_memory=False)
        alias = alias[["symbol", "alias_symbol", "prev_symbol"]]
        alias.fillna("NO", inplace=True)
    else:
        alias = None
        log.warning("From extann will not check alias name to match extann")

    mandatory = ["#CHROM", "START", "END"]
    # Variants and write to file
    headerfile = output + ".hdr"
    with open(f"{output}.tmp", "w+") as o, open(headerfile, "w+") as h:
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
        log.info("Get gene coordinates from refseq")
        for i, rows in df_extann.iterrows():
            pos_list = get_gene_coordinate(
                df_refgene, rows, extra_cols, mode,df_transcript, alias
            )
            if pos_list:
                # for each transcript from same gene
                for data in pos_list:
                    # it means that gene name comming from alias is already in data list
                    if len(data) == 5:
                        data.pop()
                    for key, val in rows.items():
                        if param.get("replace"):
                            data.append(replace_values(str(val),param.get("replace")))
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
    gene = match["name"].unique()[0]
    if df_transcript.empty:
        return get_longest_transcript(match, extra_col)
    assert "gene" in df_transcript.columns, "Gene column is missing (genes) exit"
    assert (
        "transcript" in df_transcript.columns
    ), "Transcript column is missing (transcript) exit"
    if len(df_transcript.loc[df_transcript["gene"] == gene].index) != 0:
        try:
            chosen = df_transcript.loc[df_transcript["gene"] == gene].iloc[0][
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
    alias=None
) -> any:
    """
    From pandas dataframe containing refgene file, get chr start stop from each gene present in extann
    do the same process for each gene/transcript it will lead to duplicate
    df_refgene: refgene dataframe
    gene_row: pandas series of extann row
    log
    """
    gene = gene_row["genes"]
    match = df_refgene.loc[df_refgene["name"] == gene_row["genes"]]

    if len(match.index) >= 1:
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
        if alias is not None:
        #For all transcript mode search in alias file like OMIM or search in HUGO database file
            return get_refseq_match(gene, df_refgene, alias)

def get_transcript_coordinates(data, transcript):
    """
    """
    data['Extracted'] = data['exon'].str.extract(r'(\d+)', expand=False).astype(int)
    max_value = data['Extracted'].max()
    min_value = data['Extracted'].min()
    min_value_df = data.loc[data["Extracted"] == min_value]
    max_value_df = data.loc[data["Extracted"] == max_value]
    # try:
    if data["strand"].unique() == "+":
        yield [max_value_df["#CHROM"].iloc[0], min_value_df["START"].iloc[0], max_value_df["END"].iloc[0], transcript, max_value_df["name"].iloc[0]]
    else:
        yield [max_value_df["#CHROM"].iloc[0], max_value_df["START"].iloc[0], min_value_df["END"].iloc[0], transcript, max_value_df["name"].iloc[0]]

def get_refseq_match(gene: str, refseq: pd.DataFrame, alias: pd.DataFrame) -> list:
    """
    Get from ncbi all transcript coordinate for a gene. Loop over hugo alias gene name

    :param gene: gene name from gnomad gene file
    :param refseq: path of refseq file from ncbi
    :return: cordinate in string of each transcript with transcript name and gene name
    """

    res = refseq.loc[refseq["name"] == gene]
    if len(res.index) == 0:
        alias_list = get_aliases(gene, alias)
        if not alias_list:
            return []
        log.info(f"Check alias gene name {gene}, list {' '.join(alias_list)}")
        found = []
        for snd_name in alias_list:
            res = refseq.loc[refseq["name"] == snd_name]
            if len(res.index) > 1:
                for transcript, data in res.groupby("transcript"):
                   found.extend(list(get_transcript_coordinates(data, transcript)))
            if found:
                log.debug(f"Found match {snd_name}")
                return found
        log.warning(f"{gene} not in refseq")
        return []
    else:
        raise ValueError("Dataframe is not empty EXIT")

def find_rows_with_substring(df: pd.DataFrame, substring:str):
    """
    Catch value in entire dataframe
    """
    mask = df.map(lambda cell: substring in str(cell))
    if mask.any().any():
        matching_rows = df[mask.any(axis=1)]
        return matching_rows
    else:
        return ""


def get_aliases(gene: str, alias: pd.DataFrame) -> list:
    """
    Get alias gene name from HGNC file
    
    :param gene: gene name from extann raw
    :param alias: HNGC dataframe from raw txt file
    """
    try:
        alias_gene = find_rows_with_substring(alias, gene).values.tolist()[0]

    except AttributeError:
        return []
    alias_gene_splitted = []
    for symbol in alias_gene:
        if "|" in symbol:
            symbol = symbol.split("|")
            alias_gene_splitted.extend(symbol)
        else:
            alias_gene_splitted.append(symbol)
    try:
        alias_gene_splitted.remove(gene)
    except ValueError:
        log.warning(
            f"Partial match in aliases for gene: {gene}, list {alias_gene_splitted}, skip record"
        )
        return []
    
    log.debug(alias_gene_splitted)
    return alias_gene_splitted

def from_extann(args: argparse) -> None:
    """
    This function converts an txt or tsv files containing genes-bases information\n
    From a "genes" columns which contains genes symbol it will match gene coordinates in refgene database and create a bed-like output with vcf header

     :param args: `args` is an object with several attributes representing the input parameters for the
     function. These attributes include:
     :type args: argparse
    """
    param = read_json(args.param_extann)
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
    if args.refgene:
        df_refgene = read_refgene(args.refgene)
    else:
        df_refgene = read_refgene(param.get("refgene"))

    if args.input_extann.endswith(".gz"):
        df_extann = pd.read_csv(
            args.input_extann, header=0, sep="\t", compression="gzip"
        )
    else:
        df_extann = pd.read_csv(args.input_extann, header=0, sep="\t")

    # Mode
    if not args.mode_extann:
        assert param.get(
            "mode_extann"
        ), "No mode_extann was provided either in cli or as parameter EXIT"

    # Transcript
    if args.transcripts:
        df_transcript = commons.transcripts_file_to_df(args.transcripts)
    elif param.get("transcripts"):
        df_transcript = commons.transcripts_file_to_df(param.get("transcripts"))
    else:
        if param.get("mode_extann") == "chosen":
            log.error(
                "Extann mode is set to user-specific but no transcript file provided EXIT"
            )
            exit()
        else:
            log.debug("No transcript provided for extann")
            df_transcript = None
    # Param
    log.debug(f"Replace {str(param.get('replace'))}")
    log.debug(f"Config mode_extann {args.mode_extann}")
    log.debug(f"Config extra_col_list {' '.join(param['extra_col_list'])}")
    log.debug(f"Config transcripts {args.transcripts}")
    log.debug(f"Config refgene {args.refgene}")

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
        args.mode_extann,
        df_transcript
    )
    commons.command(
        f"grep '^#' {args.output_extann}.tmp > {args.output_extann} && grep -v '^#' {args.output_extann}.tmp | sort -k1,1V -k2,2n >> {args.output_extann}"
    )
    log.debug("Removing tmp output")
    os.remove(f"{args.output_extann}.tmp")

    if args.output_extann.endswith(".gz"):
        log.debug("Compressing output")
        commons.command(
            f"mv {args.output_extann} {args.output_extann.replace('.gz', '')} && bgzip {args.output_extann.replace('.gz', '')}"
        )
