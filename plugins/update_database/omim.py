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
    """
    :param refseq: Path to refseq file
    """

    def __init__(
        self,
        link=None,
        database=None,
        exclude_link=None,
        databases_folder=None,
        config_json=None,
        current_folder=None,
        refseq=None,
        data_folder=None
    ):
        super().__init__(
            link,
            database,
            exclude_link,
            databases_folder,
            "OMIMannotations",
            config_json,
            current_folder,
        )
        self.header = self.create_header(
            self.config_json.get("header").get(self.database), "tsv"
        )
        self.data_folder = data_folder
        self.refseq = refseq
        self.refgene = self.get_refgene_list()
        if not os.path.exists(osj(self.data_folder, "genemap_filtered.txt")):
            log.error(f"genemap filtered is not in {self.data_folder} folder (grep -v \"^#\" genemap2.txt | cut -f9,13 > genemap_filtered.txt)")
            raise FileExistsError()
        else:
            self.genemap = self.create_genemap()

    def write_omim(self, df_gb_sorted: pd.DataFrame, output: str):
        """
        Write the final OMIM file bgzip compressed

        :param df_gb_sorted: pandas dataframe with OMIM data
        :param output: path to write the file
        """
        with BgzfWriter(output, "wt") as o:
            for lines in self.header:
                o.write(lines)
            df_gb_sorted.to_csv(o, header=True, index=False, sep="\t", mode="a")

    def get_refgene_list(self) -> list:
        """
        there are non coding transcript inside so more than 20K gene name is  normal
        """
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

    @staticmethod
    def agg_morbid_yes(x):
        if len(list(set(x))) > 1:
            raise ValueError(f"{' '.join(x)} are not in allowed values (Yes or empty)")
        return list(set(x))[0]

    @staticmethod
    def refactor_gene_column(row):
        val1, val2 = row['genes_x'], row['genes_y']
        if val1 == val2:  
            return val1
        elif isinstance(val1, float):  
            return val2
        elif isinstance(val2, float):  
            return val1
        elif isinstance(val1, str) and isinstance(val2, str):
            return val1 if len(val1.split(",")) >= len(val2.split(",")) else val2
        else:
            raise ValueError(f"Conflict in row: '{val1}' != '{val2}'")

    @staticmethod
    def get_header(file):
        header = []
        with gzip.open(file, "rt") as  f:
            for i, lines in enumerate(f):
                if lines.startswith("#CHROM"):
                    break
                else:
                    header.append(lines)
        return header, i

    def load_data(self, file_path):
        _, index = self.get_header(file_path)
        df = pd.read_csv(file_path, skiprows=index, header=0, compression="gzip", sep="\t")
        return df
    
    def update_omim_morbid(self, df: pd.DataFrame, name:str) -> pd.DataFrame:
        """
        From morbid or morbid candidate dataframe return dataframe groupby transcript and sorted
        """
        log.debug(f"{name}: {len(df.index)} rows")
        df_gb = df.groupby("transcript", as_index=False)
        df_list = []
        for i, df_tmp in df_gb:
            df_final = df_tmp.groupby("transcript", as_index=False).agg(
             {'genes': lambda x: ','.join(list(set(x))),
              '#CHROM': 'first',
              'START': 'first',
              'END': 'first',
              name: lambda x: self.agg_morbid_yes(x)
              })
            df_list.append(df_final)

        df_gb = pd.concat(df_list, ignore_index=True)
        df_gb = df_gb[["#CHROM", "START", "END", "transcript", "genes", name]]
        df_gb_sorted = df_gb.sort_values(
            by=["#CHROM", "START", "END"],
            key=lambda x: x.map(self.chrom_sort_key) if x.name == "#CHROM" else x
        )
        log.debug(f"{name} processed sorted: {len(df_gb_sorted.index)} rows")
        return df_gb_sorted
    
    def update_omim_raw(self, df:pd.DataFrame):
        log.debug(f"OMIM raw: {len(df.index)} rows")
        log.debug(f"Groupby transcript ...")
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
        log.debug(f"OMIM raw processed sorted: {len(df_gb_sorted)} rows")
        return df_gb_sorted
    
    def is_in_genemap(self, df: pd.DataFrame, row_to_check:dict) -> bool:
        """
        Check whether a gene / phenotype couple is in genemap2 raw file.
        True meaning values need to be keep otherwise discarded from OMIM bundle file
        Be carefull, genes could be a list
        """
        for genes in row_to_check["genes"].split(","):
            pres = ((df["genes"] == genes) & (~df['OMIM_phenotype'].isna())).any()
            if pres:
                return True
        return False
    
    def get_problematic_id(self, df:pd.DataFrame) -> list:
        """
        Get from OMIM bundle database OMIM_ID where ID is associated to different genes
        """
        mim_id = []
        for id, df_g in df.groupby("OMIM_ID"):
            if len(df_g["START"].unique().tolist()) > 1:
                mim_id.append(id)
        return mim_id
    
    def get_problematic_index_to_delete(self, df: pd.DataFrame, genemap: pd.DataFrame, mim_id: list) -> list:
        """
        Get list of index to delete for OMIM id matching two different genes due to alias association and which do not have an associated phenotype
        """
        id_to_delete = []
        df_double = df.loc[df["OMIM_ID"].isin(mim_id)]
        for omim, df_omim in df_double.groupby("OMIM_ID"):
            for start, df_start in df_omim.groupby("START"):
                if self.is_in_genemap(genemap, {"genes": df_start["genes"].iloc[0], "OMIM_phenotype": df_start["OMIM_phenotype"].iloc[0]}):
                    continue
                else:
                    id_to_delete.extend(list(df_omim.loc[df_omim["START"] == start].index.values))
        return id_to_delete
    
    def create_genemap(self):
        """
        
        """
        genemap = pd.read_csv(osj(self.data_folder, "genemap_filtered.txt"), sep="\t", comment="#", header=None)
        genemap.columns = ["genes", "OMIM_phenotype"]
        log.debug(f"genemap filtered {len(genemap.index)} rows")
        return genemap

    def update_omim(self):
        """
        Take OMIMannotations*.bed.gz and remove genes that are not in refseq for alias and merge all files
        """
        files = {}
        for file in ["OMIMannotations.bed.gz", "OMIMannotations.morbid.bed.gz", "OMIMannotations.morbidCandidates.bed.gz"]:
            if not os.path.exists(osj(self.data_folder, file)):
                log.error(f"Missing {file} in {self.data_folder}")
                raise FileNotFoundError()
            else:
                files[file.replace(".bed.gz", "")] = osj(self.data_folder, file)
        log.debug(f"Inputs: {files}")

        #Create dataframes
        log.info("Load dataframe and merge alias gene")
        df_omim_sorted = self.update_omim_raw(self.load_data(files["OMIMannotations"]))
        df_morbid_sorted = self.update_omim_morbid(self.load_data(files["OMIMannotations.morbid"]), "OMIM_morbid")
        df_morbidcandidate_sorted = self.update_omim_morbid(self.load_data(files["OMIMannotations.morbidCandidates"]), "OMIM_morbid_candidate")
        
        #Merge OMIM databases
        log.info("Merging OMIM databases")
        omim_bundle_tmp = pd.merge(df_omim_sorted, df_morbid_sorted, on=["#CHROM", "START", "END", "transcript"], how="outer")
        omim_bundle_tmp["genes"] = omim_bundle_tmp.apply(lambda x: self.refactor_gene_column(x), axis=1)
        omim_bundle_tmp.drop(columns=['genes_x', 'genes_y'], inplace=True)
        omim_bundle = pd.merge(omim_bundle_tmp, df_morbidcandidate_sorted, on=["#CHROM", "START", "END", "transcript"], how="outer")
        omim_bundle["genes"] = omim_bundle.apply(lambda x: self.refactor_gene_column(x), axis=1)
        omim_bundle.drop(columns=['genes_x', 'genes_y'], inplace=True)

        log.debug(f"OMIM bundle {len(omim_bundle.index)}")

        #Drop ambigous row
        problematic_id = self.get_problematic_id(omim_bundle)
        problematic_id_to_delete = self.get_problematic_index_to_delete(omim_bundle, self.genemap, problematic_id)
        log.debug(f"OMIM bundle drop {len(omim_bundle.loc[omim_bundle.index.isin(problematic_id_to_delete)].index)} ambiguous rows")
        omim_bundle_filtered = omim_bundle.drop(problematic_id_to_delete)


        # drop empty phenotype
        log.debug(f"OMIM bundle dropna phenotype: {len(omim_bundle_filtered.loc[omim_bundle_filtered['OMIM_phenotype'].isnull()].index)} rows")
        omim_bundle_filtered.dropna(subset=['OMIM_phenotype'], inplace=True)


        omim_bundle_filtered_ambiguous = omim_bundle_filtered.loc[~omim_bundle_filtered["genes"].apply(lambda x: any(val.strip() in self.refgene for val in x.split(',')))]
        log.debug(f"OMIM bundle {len(omim_bundle_filtered_ambiguous.index)} row with gene not in refseq")


        #Check ambiguous gene from provided list
        omim_ambiguous_check = omim_bundle_filtered_ambiguous.loc[~omim_bundle_filtered_ambiguous["genes"].str.contains('|'.join(self.config_json["omim_ambiguous"]["full"]), case=False, na=False)]
        if len(omim_ambiguous_check.index) >= 1:
            log.warning(f"{' '.join(omim_ambiguous_check['genes'].unique().tolist())} genes not listed in omim_ambiguous, probably due to version changing, check those genes !!!")

        #Keep genes checked manually, Save listed amibuous gene name from omim_ambiguous in update_databases.json
        index_genes_not_refseq = omim_bundle_filtered.index[(~omim_bundle_filtered["genes"].apply(lambda x: any(val.strip() in self.refgene for val in x.split(',')))) & (~omim_bundle_filtered['genes'].str.contains('|'.join(self.config_json["omim_ambiguous"]["keep"]), case=False, na=False))].tolist()
        
        log.debug(f"Saved {' '.join(self.config_json['omim_ambiguous']['keep'])} genes / phenotype")
        log.debug(f"OMIM bundle remove {len(index_genes_not_refseq)} row with gene not in refseq")
        omim_bundle_filtered.drop(index_genes_not_refseq, inplace=True)


        #Writting
        log.debug(f"OMIM bundle write {len(omim_bundle_filtered.index)} rows in final db")

        output = osj(self.databases_folder, "OMIMannotations.final.bed.gz")
        with BgzfWriter(output, "wt") as o:
            for lines in self.header[:-1]:
                o.write(lines + "\n")
            omim_bundle_filtered.to_csv(o, header=True, index=False, sep="\t", mode="a")
        with open(output + ".hdr", "w+") as h:
            for lines in self.header:
                h.write(lines + "\n")
