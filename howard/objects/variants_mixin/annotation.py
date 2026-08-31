import copy
import datetime
import gzip
import io
import os
from pathlib import Path
import random
import re
import subprocess
from tempfile import NamedTemporaryFile, TemporaryDirectory
import tempfile
import json
import yaml  # type: ignore
import Bio.bgzf as bgzf  # type: ignore
import pandas as pd  # type: ignore
import polars as pl  # type: ignore
from pyfaidx import Fasta  # type: ignore
import numpy as np  # type: ignore
import vcf  # type: ignore
import logging as log
import cyvcf2  # type: ignore
import pyBigWig  # type: ignore


from howard.functions.commons import (
    CODE_TYPE_MAP,
    DEFAULT_ANNOTATIONS_FOLDER,
    DEFAULT_CHUNK_SIZE,
    DEFAULT_PARQUET_FOLDER,
    DEFAULT_BCFTOOLS_FOLDER,
    DEFAULT_BIGWIG_FOLDER,
    DEFAULT_ANNOVAR_FOLDER,
    DEFAULT_ASSEMBLY,
    DEFAULT_DATABASE_FOLDER,
    DEFAULT_GENOME_FOLDER,
    DEFAULT_REFSEQ_FOLDER,
    DEFAULT_SNPEFF_FOLDER,
    DEFAULT_SPLICE_FOLDER,
    DEFAULT_TOOLS_BIN,
    DEFAULT_TOOLS_FOLDER,
    annotation_file_find,
    check_docker_image_exists,
    clean_annotation_field,
    command,
    duckdb_has_spilled,
    find,
    find_file_prefix,
    find_genome,
    full_path,
    get_bin_command,
    get_random,
    merge_regions,
    params_string_to_dict,
    remove_if_exists,
    run_parallel_commands,
    folder_config,
    code_type_map,
    get_random,
    ChromMapping
)

from howard.functions.databases import (
    databases_download_annovar,
    databases_download_exomiser,
    databases_download_snpeff,
    databases_infos,
)

from howard.functions.utils import (
    format_hgvs_name,
    get_refseq_table,
    get_transcript,
    read_transcripts,
)

from howard.objects.database import Database
from howard.objects.variants_mixin.annotation_docker import variants_annotation_docker
# from howard.objects.variants import Variants


class variants_annotation(variants_annotation_docker):

    ##############
    # Annotation #
    ##############

    def scan_databases(
        self,
        database_formats: list = ["parquet"],
        database_releases: list = ["current"],
    ) -> dict:
        """
        The function `scan_databases` scans for available databases based on specified formats and
        releases.

        :param database_formats: The `database_formats` parameter is a list that specifies the formats
        of the databases to be scanned. In this case, the accepted format is "parquet"
        :type database_formats: list ["parquet"]
        :param database_releases: The `database_releases` parameter is a list that specifies the
        releases of the databases to be scanned. In the provided function, the default value for
        `database_releases` is set to `["current"]`, meaning that by default, the function will scan
        databases that are in the "current"
        :type database_releases: list
        :return: The function `scan_databases` returns a dictionary containing information about
        databases that match the specified formats and releases.
        """

        # Config
        config = self.get_config()

        # Param
        param = self.get_param()

        # Param - Assembly
        assembly = self.get_assembly()

        # Scan for availabled databases
        log.info(
            f"Annotations - Check annotation parameters - Scan existing databases - Assembly {[assembly]} - Formats {database_formats} - Releases {database_releases}..."
        )
        databases_infos_dict = databases_infos(
            database_folder_releases=database_releases,
            database_formats=database_formats,
            assembly=assembly,
            config=config,
        )
        log.info(
            f"Annotations - Check annotation parameters - Scan existing databases - {len(databases_infos_dict)} databases found"
        )

        return databases_infos_dict

    def annotation(self, section:str = "annotation", param: dict = None) -> None:
        """
        It annotates the VCF file with the annotations specified in the config file.
        """

        # Config
        config = self.get_config()

        # Param
        if param is None:
            param = self.get_param()

        # Param - Assembly
        assembly = self.get_assembly()

        # annotations databases folders
        annotations_databases = set(
            config.get("folders", {})
            .get("databases", {})
            .get("annotations", [DEFAULT_ANNOTATIONS_FOLDER])
            + config.get("folders", {})
            .get("databases", {})
            .get("parquet", [DEFAULT_PARQUET_FOLDER])
            + config.get("folders", {})
            .get("databases", {})
            .get("bcftools", [DEFAULT_BCFTOOLS_FOLDER])
        )

        # Get param annotations
        if param.get("annotations", None) and isinstance(
            param.get("annotations", None), str
        ):
            param_annotation_list = param.get("annotations").split(",")
        else:
            param_annotation_list = []

        # Each tools param
        if param.get("annotation_parquet", None) != None:
            if isinstance(param.get("annotation_parquet", None), list):
                param_annotation_list.append(",".join(param.get("annotation_parquet")))
            elif isinstance(param.get("annotation_parquet", None), str):
                param_annotation_list.append(param.get("annotation_parquet"))
        if param.get("annotation_snpsift", None) != None:
            if isinstance(param.get("annotation_snpsift", None), list):
                param_annotation_list.append(
                    "snpsift:"
                    + "+".join(param.get("annotation_snpsift")).replace(",", "+")
                )
            elif isinstance(param.get("annotation_snpsift", None), str):
                param_annotation_list.append(
                    "snpsift:" + param.get("annotation_snpsift").replace(",", "+")
                )
        if param.get("annotation_snpeff", None) != None and isinstance(param.get("annotation_snpeff", None), str):
            param_annotation_list.append("snpeff:" + param.get("annotation_snpeff"))
        if param.get("annotation_bcftools", None) != None:
            if isinstance(param.get("annotation_bcftools", None), list):
                param_annotation_list.append(
                    "bcftools:"
                    + "+".join(param.get("annotation_bcftools")).replace(",", "+")
                )
            elif isinstance(param.get("annotation_bcftools", None), str):
                param_annotation_list.append(
                    "bcftools:" + param.get("annotation_bcftools").replace(",", "+")
                )
        if param.get("annotation_annovar", None) != None and isinstance(param.get("annotation_annovar", None), str):
            param_annotation_list.append("annovar:" + param.get("annotation_annovar"))
        if param.get("annotation_exomiser", None) != None and isinstance(param.get("annotation_exomiser", None), str):
            param_annotation_list.append("exomiser:" + param.get("annotation_exomiser"))
        if param.get("annotation_splice", None) != None and isinstance(param.get("annotation_splice", None), str):
            param_annotation_list.append("splice:" + param.get("annotation_splice"))

        # Merge param annotations list
        param["annotations"] = ",".join(param_annotation_list)

        if "fast" in param and param["fast"]:
            fast = True
        else:
            fast = False
        log.debug(f"Fast mode: {fast}")

        # # debug
        # log.debug(f"param_annotations={param['annotations']}")

        if param.get("annotations"):

            # Log
            # log.info("Annotations - Check annotation parameters")

            if not section in param:
                param[section] = {}

            # List of annotations parameters
            annotations_list_input = {}
            if isinstance(param.get("annotations", None), str):
                annotation_file_list = list(param.get("annotations", "").split(","))
                for annotation_file in annotation_file_list:
                    annotations_list_input[annotation_file.strip()] = {"INFO": None}
            else:
                annotations_list_input = param.get("annotations", {})

            log.info(f"Quick Annotations:")
            for annotation_key in list(annotations_list_input.keys()):
                log.info(f"   {annotation_key}")

            # List of annotations and associated fields
            annotations_list = {}

            for annotation_file in annotations_list_input:

                # Explode annotations if ALL
                if (
                    annotation_file.upper() == "ALL"
                    or annotation_file.upper().startswith("ALL:")
                ):

                    # check ALL parameters (formats, releases)
                    annotation_file_split = annotation_file.split(":")
                    database_formats = "parquet"
                    database_releases = "current"
                    for annotation_file_option in annotation_file_split[1:]:
                        database_all_options_split = annotation_file_option.split("=")
                        if database_all_options_split[0] == "format":
                            database_formats = database_all_options_split[1].split("+")
                        if database_all_options_split[0] == "release":
                            database_releases = database_all_options_split[1].split("+")

                    # Scan for availabled databases
                    databases_infos_dict = self.scan_databases(
                        database_formats=database_formats,
                        database_releases=database_releases,
                    )

                    # Add found databases in annotation parameters
                    for database_infos in databases_infos_dict.keys():
                        annotations_list[database_infos] = {"INFO": None}

                else:
                    annotations_list[annotation_file] = annotations_list_input[
                        annotation_file
                    ]

            # Check each databases
            if len(annotations_list):

                log.info(
                    f"Annotations - Check annotation parameters - Check {len(annotations_list)} databases..."
                )

                for annotation_file in annotations_list:

                    # Init
                    annotations = annotations_list.get(annotation_file, None)

                    # Annotation snpEff
                    if annotation_file.startswith("snpeff"):

                        log.debug(f"Quick Annotation snpEff")

                        if "snpeff" not in param[section]:
                            param[section]["snpeff"] = {}

                        if "options" not in param[section]["snpeff"]:
                            param[section]["snpeff"]["options"] = ""

                        # snpEff options in annotations
                        param[section]["snpeff"]["options"] = "".join(
                            annotation_file.split(":")[1:]
                        )

                    # Annotation VEP
                    elif annotation_file.startswith("vep"):

                        log.debug(f"Quick Annotation VEP")

                        if "vep" not in param[section]:
                            param[section]["vep"] = {}

                        if "parameters" not in param[section]["vep"]:
                            param[section]["vep"]["parameters"] = ""

                        # VEP options in annotations
                        param[section]["vep"]["parameters"] = [f"--{x}" if not x.startswith("--") else x for x in annotation_file.replace(" ", ":").split(":")[1:] if x]

                    # Annotation Annovar
                    elif annotation_file.startswith("annovar"):

                        log.debug(f"Quick Annotation Annovar")

                        if "annovar" not in param[section]:
                            param[section]["annovar"] = {}

                        if "annotations" not in param[section]["annovar"]:
                            param[section]["annovar"]["annotations"] = {}

                        # Options
                        annotation_file_split = annotation_file.split(":")
                        for annotation_file_annotation in annotation_file_split[1:]:
                            if annotation_file_annotation:
                                param[section]["annovar"]["annotations"][
                                    annotation_file_annotation
                                ] = annotations

                    # Annotation Exomiser
                    elif annotation_file.startswith("exomiser"):

                        log.debug(f"Quick Annotation Exomiser")

                        param[section]["exomiser"] = params_string_to_dict(
                            annotation_file
                        )

                    # Annotation Splice
                    elif annotation_file.startswith("splice"):

                        log.debug(f"Quick Annotation Splice")

                        param[section]["splice"] = params_string_to_dict(
                            annotation_file
                        )

                    # Annotation Parquet or BCFTOOLS
                    else:

                        # Tools detection
                        if annotation_file.startswith("parquet:"):
                            annotation_tool_initial = "parquet"
                            annotation_file = ":".join(annotation_file.split(":")[1:])
                        elif annotation_file.startswith("bcftools:"):
                            annotation_tool_initial = "bcftools"
                            annotation_file = ":".join(annotation_file.split(":")[1:])
                        elif annotation_file.startswith("snpsift:"):
                            annotation_tool_initial = "snpsift"
                            annotation_file = ":".join(annotation_file.split(":")[1:])
                        elif annotation_file.startswith("bigwig:"):
                            annotation_tool_initial = "bigwig"
                            annotation_file = ":".join(annotation_file.split(":")[1:])
                        else:
                            annotation_tool_initial = None

                        # list of files
                        annotation_file_list = annotation_file.replace("+", ":").split(
                            ":"
                        )

                        for annotation_file in annotation_file_list:

                            if annotation_file:

                                # Annotation tool initial
                                annotation_tool = annotation_tool_initial

                                # Find file
                                annotation_file_found = annotation_file_find(
                                    annotation_file=annotation_file,
                                    databases_folders=list(annotations_databases),
                                    assembly=assembly,
                                )

                                # Full path
                                annotation_file_found = full_path(annotation_file_found)

                                if annotation_file_found:

                                    database = Database(database=annotation_file_found)
                                    quick_annotation_format = database.get_format()
                                    quick_annotation_is_compressed = (
                                        database.is_compressed()
                                    )
                                    quick_annotation_is_indexed = os.path.exists(
                                        f"{annotation_file_found}.tbi"
                                    )
                                    bcftools_preference = False

                                    # Check Annotation Tool
                                    if not annotation_tool:
                                        if (
                                            bcftools_preference
                                            and quick_annotation_format
                                            in ["vcf", "bed"]
                                            and quick_annotation_is_compressed
                                            and quick_annotation_is_indexed
                                        ):
                                            annotation_tool = "bcftools"
                                        elif quick_annotation_format in [
                                            "vcf",
                                            "bed",
                                            "tsv",
                                            "tsv",
                                            "csv",
                                            "json",
                                            "tbl",
                                            "parquet",
                                            "duckdb",
                                        ]:
                                            annotation_tool = "parquet"
                                        elif quick_annotation_format.lower() in [
                                            "bw",
                                            "bb",
                                            "bigwig",
                                            "bigbed",
                                        ]:
                                            annotation_tool = "bigwig"
                                        else:
                                            log.error(
                                                f"Quick Annotation File {annotation_file_found} - Format {quick_annotation_format} not supported yet"
                                            )
                                            raise ValueError(
                                                f"Quick Annotation File {annotation_file_found} - Format {quick_annotation_format} not supported yet"
                                            )

                                    log.debug(
                                        f"Quick Annotation File {annotation_file} - Annotation tool: {annotation_tool}"
                                    )

                                    # Annotation Tool dispatch
                                    if annotation_tool:
                                        if annotation_tool not in param[section]:
                                            param[section][annotation_tool] = {}
                                        if (
                                            "annotations"
                                            not in param[section][annotation_tool]
                                        ):
                                            param[section][annotation_tool][
                                                "annotations"
                                            ] = {}
                                        param[section][annotation_tool][
                                            "annotations"
                                        ][annotation_file_found] = annotations

                                else:
                                    log.warning(
                                        f"Quick Annotation File {annotation_file} does NOT exist"
                                    )

                self.set_param(param)

        if param.get(section, None):
            log.info("Annotations")
            if param.get(section, {}).get("parquet", None):
                log.info("Annotations 'parquet'...")
                self.annotation_parquet(section=section)
            if not fast:
                if param.get(section, {}).get("bcftools", None):
                    log.info("Annotations 'bcftools'...")
                    self.annotation_bcftools(section=section)
                if param.get(section, {}).get("snpsift", None):
                    log.info("Annotations 'snpsift'...")
                    self.annotation_snpsift(section=section)
                if param.get(section, {}).get("bigwig", None):
                    log.info("Annotations 'bigwig'...")
                    self.annotation_bigwig(section=section)
                if param.get(section, {}).get("annovar", None):
                    log.info("Annotations 'annovar'...")
                    self.annotation_annovar(section=section)
                if param.get(section, {}).get("snpeff", None):
                    log.info("Annotations 'snpeff'...")
                    self.annotation_snpeff(section=section)
                if param.get(section, {}).get("vep", None):
                    log.info("Annotations 'vep'...")
                    self.annotation_vep(section=section)
                if param.get(section, {}).get("exomiser", None) is not None:
                    log.info("Annotations 'exomiser'...")
                    self.annotation_exomiser(section=section)
                if param.get(section, {}).get("splice", None) is not None:
                    log.info("Annotations 'splice' ...")
                    self.annotation_splice(section=section)
                if param.get(section, {}).get("docker", None) is not None:
                    log.info("Annotations 'docker' ...")
                    self.annotation_docker(section=section)
                if param.get(section, {}).get("hgvs", None) is not None:
                    log.info("Annotations 'hgvs' ...")
                    self.annotation_hgvs(section=section)
            else:
                log.warning("Fast mode - skip some annotations tools")

    def annotation_bigwig(self, section:str = "annotation", threads: int = None) -> None:
        """
        The function `annotation_bigwig` annotates variants in a VCF file using bigwig databases.

        :param threads: The `threads` parameter in the `annotation_bigwig` method is used to specify the
        number of threads to be used for parallel processing during the annotation process. If the
        `threads` parameter is not provided, the method will attempt to determine the optimal number of
        threads to use based on the system configuration
        :type threads: int
        :return: True
        """

        from howard.objects.variants import Variants

        # DEBUG
        log.debug("Start annotation with bigwig databases")

        # Threads
        if not threads:
            threads = self.get_threads()
        log.debug("Threads: " + str(threads))

        # Config
        config = self.get_config()
        log.debug("Config: " + str(config))

        # Config - BCFTools databases folders
        databases_folders = set(
            self.get_config()
            .get("folders", {})
            .get("databases", {})
            .get("annotations", [DEFAULT_ANNOTATIONS_FOLDER])
            + self.get_config()
            .get("folders", {})
            .get("databases", {})
            .get("bigwig", [DEFAULT_BIGWIG_FOLDER])
        )
        log.debug("Databases annotations: " + str(databases_folders))

        # Param
        param = self.get_param()

        # Param - Annotation
        annotations = (
            param
            .get(section, {})
            .get("bigwig", {})
            .get("annotations", None)
        )
        log.debug("Annotations: " + str(annotations))

        # Header fields override (Number/Type/Description), collected from each
        # annotation's own options block
        # (annotation.annovar.annotations.<name>.options.header_fields)
        annotation_header_fields_override = self.get_annotation_header_fields_override(
            annotations=annotations
        )
        log.debug(
            f"annotation_header_fields_override={annotation_header_fields_override}"
        )

        # Param - Annotation param
        annotations_param = (
            param.get(section, {}).get("bigwig", {}).get("param", {})
        )
        log.debug("Annotations param: " + str(annotations_param))

        # Param - Options
        annotations_options = (
            param.get(section, {}).get("bigwig", {}).get("options", {})
        )
        log.debug("Annotations options: " + str(annotations_options))

        # Param - Option chrom_mapping
        chrom_mapping_options = param.get(section, {}).get("bigwig", {}).get("chrom_mapping", None)
        chrom_mapping = ChromMapping(
            mapping=chrom_mapping_options
        )

        # Chunk size
        chunk_size = self.get_config().get("chunk_size", DEFAULT_CHUNK_SIZE)

        # Assembly
        assembly = param.get(
            "assembly", self.get_config().get("assembly", DEFAULT_ASSEMBLY)
        )

        # Data
        table_variants = self.get_table_variants()

        # Check if not empty
        log.debug("Check if not empty")
        sql_query_chromosomes = (
            f"""SELECT count(*) as count FROM {table_variants} as table_variants"""
        )
        sql_query_chromosomes_df = self.get_query_to_df(sql_query_chromosomes)
        if not sql_query_chromosomes_df["count"][0]:
            log.info(f"VCF empty")
            return

        # VCF header
        vcf_reader = self.get_header()
        log.debug("Initial header: " + str(vcf_reader.infos))

        # Existing annotations
        for vcf_annotation in self.get_header().infos:

            vcf_annotation_line = self.get_header().infos.get(vcf_annotation)
            log.debug(
                f"Existing annotations in VCF: {vcf_annotation} [{vcf_annotation_line}]"
            )

        if annotations:

            with TemporaryDirectory(dir=self.get_tmp_dir()) as tmp_dir:

                # Export VCF file
                tmp_vcf_name = os.path.join(tmp_dir, "input.vcf")

                # annotation_bigwig_config
                annotation_bigwig_config_list = []

                for annotation in annotations:

                    # Annotation Name
                    annotation_name = os.path.basename(annotation)

                    # Annotation fields
                    # Add backward compatibility if "annotation_fields" is in annotation dict (instead of directly in annotation value)
                    # Allow to add options for annotations
                    if "annotation_fields" in annotations[annotation]:
                        annotation_fields = annotations[annotation].get(
                            "annotation_fields", {}
                        )
                    else:
                        annotation_fields = annotations[annotation]

                    # Options for annotations
                    annotation_options = annotations[annotation].get("options", {})

                    # Annotation Name
                    annotation_name = os.path.basename(annotation)

                    if not annotation_fields:
                        annotation_fields = {"INFO": None}

                    log.debug(f"Annotation '{annotation_name}'")
                    log.debug(
                        f"Annotation '{annotation_name}' - fields: {annotation_fields}"
                    )

                    # Create Database
                    database = Database(
                        database=annotation,
                        databases_folders=databases_folders,
                        assembly=assembly,
                    )

                    # Find files
                    db_file = database.get_database()
                    db_file = full_path(db_file)
                    db_hdr_file = database.get_header_file()
                    db_hdr_file = full_path(db_hdr_file)
                    db_file_type = database.get_format()

                    # If db_file is http ?
                    if database.get_database().startswith("http"):

                        # Datbase is HTTP URL
                        db_file_is_http = True

                        # DB file keep as URL
                        db_file = database.get_database()
                        log.warning(
                            f"Annotations 'bigwig' database '{db_file}' - is an HTTP URL (experimental)"
                        )

                        # Retrieve automatic annotation field name
                        annotation_field = clean_annotation_field(
                            re.sub(
                                r"\.(bw|bb|bigwig|bigbed)$",
                                "",
                                os.path.basename(db_file),
                                flags=re.IGNORECASE,
                            )
                        )
                        log.debug(
                            f"Create header file with annotation field '{annotation_field}' is an HTTP URL"
                        )

                        # Create automatic header file
                        db_hdr_file = os.path.join(tmp_dir, "header.hdr")
                        with open(db_hdr_file, "w") as f:
                            f.write("##fileformat=VCFv4.2\n")
                            f.write(
                                f"""##INFO=<ID={annotation_field},Number=.,Type=Float,Description="{annotation_field} annotation from {db_file}">\n"""
                            )
                            f.write(f"#CHROM	START	END	{annotation_field}\n")

                    else:

                        # Datbase is NOT HTTP URL
                        db_file_is_http = False

                    # Check index - try to create if not exists
                    if (
                        db_file is None
                        or db_hdr_file is None
                        or (not os.path.exists(db_file) and not db_file_is_http)
                        or not os.path.exists(db_hdr_file)
                        or not db_file_type.lower() in ["bw", "bb", "bigwig", "bigbed"]
                    ):
                        # if False:
                        log.error("Annotation failed: database not valid")
                        log.error(f"Annotation annotation file: {db_file}")
                        log.error(f"Annotation annotation file type: {db_file_type}")
                        log.error(f"Annotation annotation header: {db_hdr_file}")
                        raise ValueError(
                            f"Annotation failed: database not valid - annotation file {db_file} / annotation file type {db_file_type} / annotation header {db_hdr_file}"
                        )
                    else:

                        # Log
                        log.debug(
                            f"Annotation '{annotation}' - file: "
                            + str(db_file)
                            + " and "
                            + str(db_hdr_file)
                        )

                        # Load header as VCF object
                        db_hdr_vcf = Variants(input=db_hdr_file)
                        db_hdr_vcf_header = db_hdr_vcf.get_header()
                        db_hdr_vcf_header_infos = db_hdr_vcf_header.infos
                        log.debug(
                            "Annotation database header: "
                            + str(db_hdr_vcf_header_infos)
                        )

                        # For all fields in database
                        annotation_fields_full = False
                        if "ALL" in annotation_fields or "INFO" in annotation_fields:
                            annotation_fields = {
                                key: key for key in db_hdr_vcf_header_infos
                            }
                            log.debug(
                                "Annotation database header - All annotations added: "
                                + str(annotation_fields)
                            )
                            annotation_fields_full = True

                        # Init
                        cyvcf2_header_rename_dict = {}
                        cyvcf2_header_list = []
                        cyvcf2_header_indexes = {}

                        # process annotation fields
                        for annotation_field in annotation_fields:

                            # New annotation name
                            annotation_field_new = annotation_fields[annotation_field]

                            # Check annotation field and index in header
                            if (
                                annotation_field
                                in db_hdr_vcf.get_header_columns_as_list()
                            ):
                                annotation_field_index = (
                                    db_hdr_vcf.get_header_columns_as_list().index(
                                        annotation_field
                                    )
                                    - 3
                                )
                                cyvcf2_header_indexes[annotation_field_new] = (
                                    annotation_field_index
                                )
                            else:
                                msg_err = f"Database '{db_file}' does NOT contain annotation field '{annotation_field}'"
                                log.error(msg_err)
                                raise ValueError(msg_err)

                            # Append annotation field in cyvcf2 header list
                            cyvcf2_header_rename_dict[annotation_field_new] = (
                                db_hdr_vcf_header_infos[annotation_field].id
                            )

                            # Because values are normalized (mean or join) over the region
                            # Force number to '.' if type is String
                            # Force number to '1' if type is Integer of Float
                            annotation_field_num = "."
                            if (
                                db_hdr_vcf_header_infos[annotation_field].type
                                == "String"
                            ):
                                annotation_field_num = "."
                            elif db_hdr_vcf_header_infos[annotation_field].type in [
                                "Integer",
                                "Float",
                            ]:
                                annotation_field_num = 1
                            else:
                                annotation_field_num = db_hdr_vcf_header_infos[
                                    annotation_field
                                ].num
                            cyvcf2_header_list.append(
                                {
                                    "ID": annotation_field_new,
                                    "Number": annotation_field_num,
                                    "Type": db_hdr_vcf_header_infos[
                                        annotation_field
                                    ].type,
                                    "Description": db_hdr_vcf_header_infos[
                                        annotation_field
                                    ].desc,
                                }
                            )

                            # Add header on VCF
                            db_hdr_vcf_header.infos[annotation_field_new] = db_hdr_vcf_header.infos[annotation_field]

                        # Load bigwig database
                        bw_db = pyBigWig.open(db_file, "r")

                        # Check bigwig format
                        if bw_db.isBigWig():
                            log.debug(f"Database '{db_file}' is in 'BigWig' format")
                        elif bw_db.isBigBed():
                            log.debug(f"Database '{db_file}' is in 'BigBed' format")
                        else:
                            msg_err = f"Database '{db_file}' is NOT in 'BigWig' format"
                            log.error(msg_err)
                            raise ValueError(msg_err)

                        # Annotation param method
                        # Get default method and update with annotation specific method
                        annotations_param_annotation_method = (
                            annotations_options.get("method", {})
                            .copy()
                        )
                        annotations_param_annotation_method.update(
                            annotation_options.get("method", {})
                        )

                        # For each field, determine method by type if not defined
                        for cyvcf2_header_index in cyvcf2_header_indexes.keys():
                            if (
                                cyvcf2_header_index
                                not in annotations_param_annotation_method
                            ):
                                annotations_param_annotation_method[
                                    cyvcf2_header_index
                                ] = annotations_param_annotation_method.get(
                                    db_hdr_vcf_header_infos[
                                        cyvcf2_header_rename_dict[cyvcf2_header_index]
                                    ].type,
                                    (
                                        "mean"
                                        if db_hdr_vcf_header_infos[
                                            cyvcf2_header_rename_dict[
                                                cyvcf2_header_index
                                            ]
                                        ].type
                                        in ["Integer", "Float"]
                                        else "uniq"
                                    ),
                                )

                        # log.debug(
                        #     f"Annotation method param: {annotations_param_annotation_method}"
                        # )

                        # Append annotation config
                        annotation_bigwig_config_list.append(
                            {
                                "db_file": db_file,
                                "bw_db": bw_db,
                                "bw_type": "bw" if bw_db.isBigWig() else "bb",
                                "bw_format": "BigWig" if bw_db.isBigWig() else "BigBed",
                                "vcf_reader": db_hdr_vcf_header, # VCF header with new annotation fields
                                "cyvcf2_header_rename_dict": cyvcf2_header_rename_dict,
                                "cyvcf2_header_list": cyvcf2_header_list,
                                "cyvcf2_header_indexes": cyvcf2_header_indexes,
                                "bigwig_annotation": list(cyvcf2_header_indexes.keys())[
                                    0
                                ],
                                "annotation_method": annotations_param_annotation_method,
                            }
                        )
                        # log.debug(
                        #     f"Annotation config: {annotation_bigwig_config_list[-1]}"
                        # )

                # Annotate
                if annotation_bigwig_config_list:

                    # Export VCF file
                    # No sort and index needed
                    self.export_variant_vcf(
                        vcf_file=tmp_vcf_name,
                        remove_info=True,
                        add_samples=False,
                        sort=False,
                        index=False,
                        threads=threads,
                        chrom_mapping_sql=chrom_mapping.to_tool_sql()
                    )

                    # Load input tmp file
                    log.debug("Load input VCF file...")
                    input_vcf = cyvcf2.VCF(tmp_vcf_name, threads=threads, lazy=True)
                    log.debug("Load input VCF file done.")

                    # Remove tmp VCF file
                    # Not needed anymore, because loaded in memory
                    remove_if_exists(tmp_vcf_name)

                    # Number of variants
                    count_variants = self.get_query_to_df(
                        f"""SELECT count(*) as count FROM {table_variants} as table_variants"""
                    )
                    num_variants = count_variants["count"][0]
                    log.debug(f"Number of variants to annotate: {num_variants}")

                    # Add header in input file
                    for annotation_bigwig_config in annotation_bigwig_config_list:
                        for cyvcf2_header_field in annotation_bigwig_config.get(
                            "cyvcf2_header_list", []
                        ):
                            log.info(
                                f"Annotations '{annotation_bigwig_config.get('bw_format','BigWig')}' database '{os.path.basename(annotation_bigwig_config.get('db_file'))}' - annotation field '{annotation_bigwig_config.get('cyvcf2_header_rename_dict',{}).get(cyvcf2_header_field.get('ID','Unknown'))}' -> '{cyvcf2_header_field.get('ID')}'"
                            )
                            input_vcf.add_info_to_header(cyvcf2_header_field)

                    # Create output VCF file
                    output_vcf_file = os.path.join(tmp_dir, "output.vcf.gz")
                    output_vcf = cyvcf2.Writer(output_vcf_file, input_vcf)

                    # Fetch variants
                    log.info("Annotations start...")
                    variant_i = 0
                    for variant in input_vcf:

                        # Log progress
                        if variant_i and variant_i % chunk_size == 0:
                            log.debug(
                                f"Annotations - Processed {variant_i} variants [{variant_i/num_variants:.2%}]..."
                            )
                        variant_i += 1

                        # For each bigwig annotation config, process annotation
                        for annotation_bigwig_config in annotation_bigwig_config_list:

                            # DB and indexes
                            bw_db = annotation_bigwig_config.get("bw_db", None)
                            bw_db_file = annotation_bigwig_config.get("db_file", None)
                            vcf_header = annotation_bigwig_config.get(
                                "vcf_reader", None
                            )
                            cyvcf2_header_indexes = annotation_bigwig_config.get(
                                "cyvcf2_header_indexes", None
                            )
                            bigwig_annotation = annotation_bigwig_config.get(
                                "bigwig_annotation", None
                            )
                            bw_type = annotation_bigwig_config.get("bw_type", "bw")

                            # Annotation with bigwig
                            if bw_type == "bw":

                                # Annotation method
                                annotation_method = annotation_bigwig_config.get(
                                    "annotation_method", {}
                                ).get(bigwig_annotation, "uniq")

                                # Retrieve value from chrom pos start end
                                # Try to get value over the variant region
                                # Bydefault "mean" is used to retrevie value of a region
                                # TODO: add parameter to choose method (mean, max, min, sum, coverage)
                                # Due to out of bound errors, use TRY/EXCEPT
                                # Location is 0-based in pyBigWig
                                # Length base on REF allele length (TODO: check for SV)
                                try:
                                    variant_value = bw_db.stats(
                                        variant.CHROM,
                                        variant.POS - 1,
                                        variant.POS + len(variant.REF) - 1,
                                        type=annotation_method,
                                    )[0]
                                except RuntimeError as e:
                                    # Output is empty list if error
                                    variant_value = None

                                # If value found and not NaN, add to variant INFO
                                if variant_value is not None and not np.isnan(
                                    variant_value
                                ):
                                    variant.INFO[bigwig_annotation] = variant_value

                            # Annotation with bigbed
                            elif bw_type == "bb":

                                # Retrieve value from chrom pos start end
                                # Try to get value over the variant region
                                # Due to out of bound errors, use TRY/EXCEPT
                                # Location is 0-based in pyBigBed
                                # Length base on REF allele length (TODO: check for SV)
                                # Result "res" is a list of tuples/entry (see pybigwig doc)
                                # example:
                                # [(3112004, 3112046, '59150\t190\t-\t0.057\t0.61\t119'), (3143183, 3143307, '8\t145\t+\t0.033\t0.48\t200'), (3143324, 3143420, '9\t240\t+\t0.101\t0.30\t480')]
                                try:
                                    res = bw_db.entries(
                                        variant.CHROM,
                                        variant.POS - 1,
                                        variant.POS + len(variant.REF) - 1,
                                    )
                                except RuntimeError as e:
                                    res = None

                                # If result found
                                if res is not None:

                                    # Create values dict by header indexes (annotations)
                                    cyvcf2_header_values = {
                                        cyvcf2_header_index: []
                                        for cyvcf2_header_index in cyvcf2_header_indexes
                                    }

                                    # For each entry
                                    for entry in res:

                                        # Extract and split values
                                        # example: (3112004, 3112046, '59150\t190\t-\t0.057\t0.61\t119')
                                        values = entry[2].split("\t")

                                        # For each header index, get value and construct dict
                                        for (
                                            cyvcf2_header_index
                                        ) in cyvcf2_header_indexes:

                                            # Get value index, depend on annotation fields required
                                            value_index = cyvcf2_header_indexes[
                                                cyvcf2_header_index
                                            ]

                                            # Get value type from header infos
                                            cyvcf2_header_index_infos_type = (
                                                vcf_header.infos.get(
                                                    cyvcf2_header_index
                                                ).type
                                            )

                                            # Value is an integer, cast it
                                            if (
                                                cyvcf2_header_index_infos_type
                                                == "Integer"
                                            ):
                                                value = int(values[value_index])

                                            # Value is a float, cast it
                                            elif (
                                                cyvcf2_header_index_infos_type
                                                == "Float"
                                            ):
                                                value = float(values[value_index])

                                            # Value is a string, keep it
                                            else:
                                                value = values[value_index]

                                            # Append value to dict
                                            cyvcf2_header_values[
                                                cyvcf2_header_index
                                            ].append(value)

                                    # Normalize values (mean for Integer/Float, join for String)
                                    # and add to variant INFO
                                    for cyvcf2_header_index in cyvcf2_header_values:

                                        # Get value type from header infos
                                        cyvcf2_header_index_infos_type = (
                                            vcf_header.infos.get(
                                                cyvcf2_header_index
                                            ).type
                                        )

                                        # Get anotation method
                                        annotation_method = (
                                            annotation_bigwig_config.get(
                                                "annotation_method", {}
                                            ).get(cyvcf2_header_index, "uniq")
                                        )

                                        # Normalize value depending on method
                                        if annotation_method == "mean":
                                            # Calculate mean (by default)
                                            value_norm = np.mean(
                                                cyvcf2_header_values[
                                                    cyvcf2_header_index
                                                ]
                                            )
                                        elif annotation_method == "max":
                                            # Calculate max
                                            value_norm = np.max(
                                                cyvcf2_header_values[
                                                    cyvcf2_header_index
                                                ]
                                            )
                                        elif annotation_method == "min":
                                            # Calculate min
                                            value_norm = np.min(
                                                cyvcf2_header_values[
                                                    cyvcf2_header_index
                                                ]
                                            )
                                        elif annotation_method == "sum":
                                            # Calculate sum
                                            value_norm = np.sum(
                                                cyvcf2_header_values[
                                                    cyvcf2_header_index
                                                ]
                                            )
                                        elif annotation_method == "coverage":
                                            # Calculate coverage (number of values)
                                            value_norm = len(
                                                cyvcf2_header_values[
                                                    cyvcf2_header_index
                                                ]
                                            )
                                        elif annotation_method == "join":
                                            # Join values
                                            value_norm = ",".join(
                                                cyvcf2_header_values[
                                                    cyvcf2_header_index
                                                ]
                                            )
                                        else:
                                            # Default to "join"
                                            value_norm = ",".join(
                                                set(
                                                    cyvcf2_header_values[
                                                        cyvcf2_header_index
                                                    ]
                                                )
                                            )

                                        # If Integer, cast to int (otherwise keep float)
                                        # Error occure when trying to add a float value to int with pycvcf2
                                        if cyvcf2_header_index_infos_type in [
                                            "Integer"
                                        ]:
                                            value_norm = int(value_norm)

                                        # Add to variant INFO
                                        variant.INFO[cyvcf2_header_index] = value_norm

                            else:
                                msg_err = f"Annotations '{annotation_bigwig_config.get('bw_format','BigWig')}' database '{bw_db_file}' - type '{bw_type}' not supported yet"
                                log.error(msg_err)
                                # continue
                                raise ValueError(msg_err)

                        # Add record in output file
                        output_vcf.write_record(variant)

                    # Close bw db
                    for annotation_bigwig_config in annotation_bigwig_config_list:

                        # DB and indexes
                        bw_db = annotation_bigwig_config.get("bw_db", None)
                        bw_db_file = annotation_bigwig_config.get("db_file", None)

                        # Try Close bw db
                        try:
                            if bw_db is not None:
                                log.debug(
                                    f"Annotations '{annotation_bigwig_config.get('bw_format','BigWig')}' file '{os.path.basename(bw_db_file)}' closing..."
                                )
                                bw_db.close()
                                log.debug(
                                    f"Annotations '{annotation_bigwig_config.get('bw_format','BigWig')}' file '{os.path.basename(bw_db_file)}' closed"
                                )
                            else:
                                log.debug(
                                    f"Annotations '{annotation_bigwig_config.get('bw_format','BigWig')}' file '{os.path.basename(bw_db_file)}' is already closed or not open"
                                )
                        except RuntimeError as e:
                            log.error(
                                f"RuntimeError while closing '{annotation_bigwig_config.get('bw_format','BigWig')}' file '{os.path.basename(bw_db_file)}': {e}"
                            )
                        except Exception as e:
                            log.error(
                                f"Unexpected error while closing '{annotation_bigwig_config.get('bw_format','BigWig')}' file '{os.path.basename(bw_db_file)}': {e}"
                            )

                    # Log
                    log.debug("Annotations done.")

                    # Close and write file
                    log.info("Annotations write...")
                    output_vcf.close()
                    log.debug("Annotations write done.")

                    # Update variants
                    log.info("Annotations update...")
                    self.update_from_vcf(
                        output_vcf_file,
                        update_header=True,
                        annotation_header_fields_override=annotation_header_fields_override,
                        chrom_mapping_sql=chrom_mapping.from_tool_sql()
                    )
                    remove_if_exists([output_vcf_file])
                    log.debug("Annotations update done.")

        return True

    def annotation_snpsift(self, section:str = "annotation", threads: int = None) -> None:
        """
        This function annotate with bcftools

        :param threads: Number of threads to use
        :return: the value of the variable "return_value".
        """

        from howard.objects.variants import Variants

        # DEBUG
        log.debug("Start annotation with bcftools databases")

        # Threads
        if not threads:
            threads = self.get_threads()
        log.debug("Threads: " + str(threads))

        # Config
        config = self.get_config()
        log.debug("Config: " + str(config))

        # Config - snpSift
        snpsift_bin_command = get_bin_command(
            bin="SnpSift.jar",
            tool="snpsift",
            bin_type="jar",
            config=config,
            default_folder=f"{DEFAULT_TOOLS_FOLDER}/snpeff",
        )
        if not snpsift_bin_command:
            msg_err = f"Annotation failed: no snpsift bin '{snpsift_bin_command}'"
            log.error(msg_err)
            raise ValueError(msg_err)

        # Config - bcftools
        bcftools_bin_command = get_bin_command(
            bin="bcftools",
            tool="bcftools",
            bin_type="bin",
            config=config,
            default_folder=f"{DEFAULT_TOOLS_FOLDER}/bcftools",
        )
        if not bcftools_bin_command:
            msg_err = f"Annotation failed: no bcftools bin '{bcftools_bin_command}'"
            log.error(msg_err)
            raise ValueError(msg_err)

        # Config - BCFTools databases folders
        databases_folders = set(
            self.get_config()
            .get("folders", {})
            .get("databases", {})
            .get("annotations", [DEFAULT_ANNOTATIONS_FOLDER])
            + self.get_config()
            .get("folders", {})
            .get("databases", {})
            .get("bcftools", [DEFAULT_BCFTOOLS_FOLDER])
        )
        log.debug("Databases annotations: " + str(databases_folders))

        # Param
        param = self.get_param()

        # Annotations section from the parameters
        annotations = (
            param
            .get(section, {})
            .get("snpsift", {})
            .get("annotations", None)
        )
        log.debug("Annotations: " + str(annotations))

        # Header fields override (Number/Type/Description), collected from each
        # annotation's own options block
        # (annotation.annovar.annotations.<name>.options.header_fields)
        annotation_header_fields_override = self.get_annotation_header_fields_override(
            annotations=annotations
        )
        log.debug(
            f"annotation_header_fields_override={annotation_header_fields_override}"
        )

        # Param - Option chrom_mapping
        chrom_mapping_options = param.get(section, {}).get("snpsift", {}).get("chrom_mapping", None)
        chrom_mapping = ChromMapping(
            mapping=chrom_mapping_options
        )

        # Assembly
        assembly = param.get(
            "assembly", self.get_config().get("assembly", DEFAULT_ASSEMBLY)
        )

        # Data
        table_variants = self.get_table_variants()

        # Check if not empty
        log.debug("Check if not empty")
        sql_query_chromosomes = (
            f"""SELECT count(*) as count FROM {table_variants} as table_variants"""
        )
        sql_query_chromosomes_df = self.get_query_to_df(sql_query_chromosomes)
        if not sql_query_chromosomes_df["count"][0]:
            log.info(f"VCF empty")
            return

        # VCF header
        vcf_reader = self.get_header()
        log.debug("Initial header: " + str(vcf_reader.infos))

        # Existing annotations
        for vcf_annotation in self.get_header().infos:

            vcf_annotation_line = self.get_header().infos.get(vcf_annotation)
            log.debug(
                f"Existing annotations in VCF: {vcf_annotation} [{vcf_annotation_line}]"
            )

        if annotations:

            with TemporaryDirectory(dir=self.get_tmp_dir()) as tmp_dir:

                # Export VCF file
                tmp_vcf_name = os.path.join(tmp_dir, "input.vcf.gz")

                # Init
                commands = {}

                for annotation in annotations:

                    # Allow to add options for annotations
                    if "annotation_fields" in annotations[annotation]:
                        annotation_fields = annotations[annotation].get(
                            "annotation_fields", {}
                        )
                    else:
                        annotation_fields = annotations[annotation]

                    # Options for annotations (disabled for now, not used)
                    # annotation_options = annotations[annotation].get("options", {})

                    # Annotation Name
                    annotation_name = os.path.basename(annotation)

                    if not annotation_fields:
                        annotation_fields = {"INFO": None}

                    log.debug(f"Annotation '{annotation_name}'")
                    log.debug(
                        f"Annotation '{annotation_name}' - fields: {annotation_fields}"
                    )

                    # Create Database
                    database = Database(
                        database=annotation,
                        databases_folders=databases_folders,
                        assembly=assembly,
                    )

                    # Find files
                    db_file = database.get_database()
                    db_file = full_path(db_file)
                    db_hdr_file = database.get_header_file()
                    db_hdr_file = full_path(db_hdr_file)
                    db_file_type = database.get_format()
                    db_tbi_file = f"{db_file}.tbi"
                    db_file_compressed = database.is_compressed()

                    # Check if compressed
                    if not db_file_compressed:
                        log.error(
                            f"Annotation '{annotation}' - {db_file} NOT compressed file"
                        )
                        raise ValueError(
                            f"Annotation '{annotation}' - {db_file} NOT compressed file"
                        )

                    # Check if indexed
                    if not os.path.exists(db_tbi_file):
                        log.error(
                            f"Annotation '{annotation}' - {db_file} NOT indexed file"
                        )
                        raise ValueError(
                            f"Annotation '{annotation}' - {db_file} NOT indexed file"
                        )

                    # Check index - try to create if not exists
                    if not os.path.exists(db_file) or not os.path.exists(db_hdr_file):
                        log.error("Annotation failed: database not valid")
                        log.error(f"Annotation annotation file: {db_file}")
                        log.error(f"Annotation annotation header: {db_hdr_file}")
                        log.error(f"Annotation annotation index: {db_tbi_file}")
                        raise ValueError(
                            f"Annotation failed: database not valid - annotation file {db_file} / annotation header {db_hdr_file} / annotation index {db_tbi_file} / annotation compression {db_file_compressed}"
                        )
                    else:

                        log.debug(
                            f"Annotation '{annotation}' - file: "
                            + str(db_file)
                            + " and "
                            + str(db_hdr_file)
                        )

                        # Load header as VCF object
                        db_hdr_vcf = Variants(input=db_hdr_file)
                        db_hdr_vcf_header_infos = db_hdr_vcf.get_header().infos
                        log.debug(
                            "Annotation database header: "
                            + str(db_hdr_vcf_header_infos)
                        )

                        # For all fields in database
                        annotation_fields_full = False
                        if "ALL" in annotation_fields or "INFO" in annotation_fields:
                            annotation_fields = {
                                key: key for key in db_hdr_vcf_header_infos
                            }
                            log.debug(
                                "Annotation database header - All annotations added: "
                                + str(annotation_fields)
                            )
                            annotation_fields_full = True

                        # # Create file for field rename
                        # log.debug("Create file for field rename")
                        # tmp_rename = NamedTemporaryFile(
                        #     prefix=self.get_prefix(),
                        #     dir=self.get_tmp_dir(),
                        #     suffix=".rename",
                        #     delete=False,
                        # )
                        # tmp_rename_name = tmp_rename.name
                        # tmp_files.append(tmp_rename_name)

                        # Number of fields
                        nb_annotation_field = 0
                        annotation_list = []
                        annotation_infos_rename_list = []

                        for annotation_field in annotation_fields:

                            # field new name, if parametered SKIPPED !!!!!! not managed actually TODO
                            annotation_fields_new_name = annotation_fields.get(
                                annotation_field, annotation_field
                            )
                            if not annotation_fields_new_name:
                                annotation_fields_new_name = annotation_field

                            # Check if field is in DB and if field is not elready in input data
                            if (
                                annotation_field in db_hdr_vcf.get_header().infos
                                and annotation_fields_new_name
                                not in self.get_header().infos
                            ):

                                log.info(
                                    f"Annotation '{annotation_name}' - '{annotation_field}' -> '{annotation_fields_new_name}'"
                                )

                                # BCFTools annotate param to rename fields
                                if annotation_field != annotation_fields_new_name:
                                    annotation_infos_rename_list.append(
                                        f"{annotation_fields_new_name}:=INFO/{annotation_field}"
                                    )

                                annotation_list.append(annotation_field)

                                nb_annotation_field += 1

                            else:

                                if (
                                    annotation_field
                                    not in db_hdr_vcf.get_header().infos
                                ):
                                    log.warning(
                                        f"Annotation '{annotation_name}' - '{annotation_field}' - not available in vcf/bed file"
                                    )
                                if (
                                    annotation_fields_new_name
                                    in self.get_header().infos
                                ):
                                    log.warning(
                                        f"Annotation '{annotation_name}' - '{annotation_fields_new_name}' - already exists (skipped)"
                                    )

                        log.info(
                            f"Annotation '{annotation_name}' - {nb_annotation_field} annotations available in vcf/bed file"
                        )

                        annotation_infos = ",".join(annotation_list)

                        if annotation_infos != "":

                            # Annotated VCF (and error file)
                            tmp_annotation_vcf_name = os.path.join(
                                tmp_dir, os.path.basename(annotation) + ".vcf.gz"
                            )
                            tmp_annotation_vcf_name_err = (
                                tmp_annotation_vcf_name + ".err"
                            )

                            # Add fields to annotate
                            if not annotation_fields_full:
                                annotation_infos_option = f"-info {annotation_infos}"
                            else:
                                annotation_infos_option = ""

                            # Info fields rename
                            if annotation_infos_rename_list:
                                annotation_infos_rename = " -c " + ",".join(
                                    annotation_infos_rename_list
                                )
                            else:
                                annotation_infos_rename = ""

                            # Annotate command
                            command_annotate = f"{snpsift_bin_command} annotate {annotation_infos_option} {db_file} {tmp_vcf_name} 2>>{tmp_annotation_vcf_name_err} | {bcftools_bin_command} annotate --threads={threads} {annotation_infos_rename} -Oz1 -o {tmp_annotation_vcf_name} 2>>{tmp_annotation_vcf_name_err} "

                            # Add command
                            commands[command_annotate] = tmp_annotation_vcf_name

                if commands:

                    # Export VCF file
                    self.export_variant_vcf(
                        vcf_file=tmp_vcf_name,
                        remove_info=True,
                        add_samples=False,
                        index=True,
                        chrom_mapping_sql=chrom_mapping.to_tool_sql()
                    )

                    # Num command
                    nb_command = 0

                    # Annotate
                    for command_annotate in commands:
                        nb_command += 1
                        log.info(
                            f"Annotation - Annotate [{nb_command}/{len(commands)}]..."
                        )
                        log.debug(f"command_annotate={command_annotate}")
                        run_parallel_commands([command_annotate], threads)

                        # Update variants
                        log.info(
                            f"Annotation - Updating [{nb_command}/{len(commands)}]..."
                        )
                        self.update_from_vcf(
                            commands[command_annotate],
                            update_header=True,
                            annotation_header_fields_override=annotation_header_fields_override,
                            chrom_mapping_sql=chrom_mapping.from_tool_sql()
                        )
                        remove_if_exists(
                            [
                                commands[command_annotate],
                                commands[command_annotate] + ".tbi",
                            ]
                        )

    def annotation_bcftools(self, section:str = "annotation", threads: int = None) -> None:
        """
        This function annotate with bcftools

        :param threads: Number of threads to use
        :return: the value of the variable "return_value".
        """

        from howard.objects.variants import Variants

        # DEBUG
        log.debug("Start annotation with bcftools databases")

        # Threads
        if not threads:
            threads = self.get_threads()
        log.debug("Threads: " + str(threads))

        # Config
        config = self.get_config()
        log.debug("Config: " + str(config))

        # DEBUG
        delete_tmp = True
        if self.get_config().get("verbosity", "warning") in ["debug"]:
            delete_tmp = False
            log.debug("Delete tmp files/folders: " + str(delete_tmp))

        # Config - BCFTools bin command
        bcftools_bin_command = get_bin_command(
            bin="bcftools",
            tool="bcftools",
            bin_type="bin",
            config=config,
            default_folder=f"{DEFAULT_TOOLS_FOLDER}/bcftools",
        )
        if not bcftools_bin_command:
            msg_err = f"Annotation failed: no bcftools bin '{bcftools_bin_command}'"
            log.error(msg_err)
            raise ValueError(msg_err)

        # Config - BCFTools databases folders
        databases_folders = set(
            self.get_config()
            .get("folders", {})
            .get("databases", {})
            .get("annotations", [DEFAULT_ANNOTATIONS_FOLDER])
            + self.get_config()
            .get("folders", {})
            .get("databases", {})
            .get("bcftools", [DEFAULT_BCFTOOLS_FOLDER])
        )
        log.debug("Databases annotations: " + str(databases_folders))

        # Param
        param = self.get_param()

        # Annotations for bcftools
        annotations = (
            self.get_param()
            .get(section, {})
            .get("bcftools", {})
            .get("annotations", None)
        )
        log.debug("Annotations: " + str(annotations))

        # Header fields override (Number/Type/Description), collected from each
        # annotation's own options block
        # (annotation.annovar.annotations.<name>.options.header_fields)
        annotation_header_fields_override = self.get_annotation_header_fields_override(
            annotations=annotations
        )
        log.debug(
            f"annotation_header_fields_override={annotation_header_fields_override}"
        )

        # Param - Option chrom_mapping
        chrom_mapping_options = param.get(section, {}).get("bcftools", {}).get("chrom_mapping", None)
        chrom_mapping = ChromMapping(
            mapping=chrom_mapping_options
        )

        # Assembly
        assembly = self.get_param().get(
            "assembly", self.get_config().get("assembly", DEFAULT_ASSEMBLY)
        )

        # Data
        table_variants = self.get_table_variants()

        # Remove files
        remove_files = []

        # Check if not empty
        log.debug("Check if not empty")
        sql_query_chromosomes = (
            f"""SELECT count(*) as count FROM {table_variants} as table_variants"""
        )
        sql_query_chromosomes_df = self.get_query_to_df(sql_query_chromosomes)
        if not sql_query_chromosomes_df["count"][0]:
            log.info(f"VCF empty")
            return

        # Export in VCF
        log.debug("Create initial file to annotate")
        tmp_vcf = NamedTemporaryFile(
            prefix=self.get_prefix(),
            dir=self.get_tmp_dir(),
            suffix=".vcf.gz",
            delete=False,
        )
        tmp_vcf_name = tmp_vcf.name

        # Remove files
        remove_files.append(tmp_vcf_name)

        # VCF header
        vcf_reader = self.get_header()
        log.debug("Initial header: " + str(vcf_reader.infos))

        # Existing annotations
        for vcf_annotation in self.get_header().infos:

            vcf_annotation_line = self.get_header().infos.get(vcf_annotation)
            log.debug(
                f"Existing annotations in VCF: {vcf_annotation} [{vcf_annotation_line}]"
            )

        if annotations:

            tmp_ann_vcf_list = []
            commands = []
            tmp_files = []
            err_files = []

            for annotation in annotations:
                
                # Allow to add options for annotations
                if "annotation_fields" in annotations[annotation]:
                    annotation_fields = annotations[annotation].get(
                        "annotation_fields", {}
                    )
                else:
                    annotation_fields = annotations[annotation]

                # Options for annotations (disabled for now, not used)
                # annotation_options = annotations[annotation].get("options", {})

                # Annotation Name
                annotation_name = os.path.basename(annotation)

                if not annotation_fields:
                    annotation_fields = {"INFO": None}

                log.debug(f"Annotation '{annotation_name}'")
                log.debug(
                    f"Annotation '{annotation_name}' - fields: {annotation_fields}"
                )

                # Create Database
                database = Database(
                    database=annotation,
                    databases_folders=databases_folders,
                    assembly=assembly,
                )

                # Find files
                db_file = database.get_database()
                db_file = full_path(db_file)
                db_hdr_file = database.get_header_file()
                db_hdr_file = full_path(db_hdr_file)
                db_file_type = database.get_format()
                db_tbi_file = f"{db_file}.tbi"
                db_file_compressed = database.is_compressed()

                # Check if compressed
                if not db_file_compressed:
                    log.error(
                        f"Annotation '{annotation}' - {db_file} NOT compressed file"
                    )
                    raise ValueError(
                        f"Annotation '{annotation}' - {db_file} NOT compressed file"
                    )

                # Check if indexed
                if not os.path.exists(db_tbi_file):
                    log.error(f"Annotation '{annotation}' - {db_file} NOT indexed file")
                    raise ValueError(
                        f"Annotation '{annotation}' - {db_file} NOT indexed file"
                    )

                # Check index - try to create if not exists
                if not os.path.exists(db_file) or not os.path.exists(db_hdr_file):
                    log.error("Annotation failed: database not valid")
                    log.error(f"Annotation annotation file: {db_file}")
                    log.error(f"Annotation annotation header: {db_hdr_file}")
                    log.error(f"Annotation annotation index: {db_tbi_file}")
                    raise ValueError(
                        f"Annotation failed: database not valid - annotation file {db_file} / annotation header {db_hdr_file} / annotation index {db_tbi_file} / annotation compression {db_file_compressed}"
                    )
                else:

                    log.debug(
                        f"Annotation '{annotation}' - file: "
                        + str(db_file)
                        + " and "
                        + str(db_hdr_file)
                    )

                    # Load header as VCF object
                    db_hdr_vcf = Variants(input=db_hdr_file)
                    db_hdr_vcf_header_infos = db_hdr_vcf.get_header().infos
                    log.debug(
                        "Annotation database header: " + str(db_hdr_vcf_header_infos)
                    )

                    # For all fields in database
                    if "ALL" in annotation_fields or "INFO" in annotation_fields:
                        annotation_fields = {
                            key: key for key in db_hdr_vcf_header_infos
                        }
                        log.debug(
                            "Annotation database header - All annotations added: "
                            + str(annotation_fields)
                        )

                    # Number of fields
                    nb_annotation_field = 0
                    annotation_list = []
                    annotation_list_to_remove = []

                    annotation_name_map = {}

                    for annotation_field in annotation_fields:

                        # field new name, if parametered SKIPPED !!!!!! not managed actually TODO
                        annotation_fields_new_name = annotation_fields.get(
                            annotation_field, annotation_field
                        )
                        if not annotation_fields_new_name:
                            annotation_fields_new_name = annotation_field

                        # Check if field is in DB and if field is not elready in input data
                        if (
                            annotation_field in db_hdr_vcf.get_header().infos
                            and annotation_fields_new_name
                            not in self.get_header().infos
                        ):

                            log.info(
                                f"Annotation '{annotation_name}' - '{annotation_field}' -> '{annotation_fields_new_name}'"
                            )

                            # Add INFO field to header
                            db_hdr_vcf_header_infos_number = (
                                db_hdr_vcf_header_infos[annotation_field].num or "."
                            )
                            db_hdr_vcf_header_infos_type = (
                                db_hdr_vcf_header_infos[annotation_field].type
                                or "String"
                            )
                            db_hdr_vcf_header_infos_description = (
                                db_hdr_vcf_header_infos[annotation_field].desc
                                or f"{annotation_field} description"
                            )
                            db_hdr_vcf_header_infos_source = (
                                db_hdr_vcf_header_infos[annotation_field].source
                                or "unknown"
                            )
                            db_hdr_vcf_header_infos_version = (
                                db_hdr_vcf_header_infos[annotation_field].version
                                or "unknown"
                            )

                            vcf_reader.infos[annotation_fields_new_name] = (
                                vcf.parser._Info(
                                    annotation_fields_new_name,
                                    db_hdr_vcf_header_infos_number,
                                    db_hdr_vcf_header_infos_type,
                                    db_hdr_vcf_header_infos_description,
                                    db_hdr_vcf_header_infos_source,
                                    db_hdr_vcf_header_infos_version,
                                    self.code_type_map[db_hdr_vcf_header_infos_type],
                                )
                            )

                            # annotation_list.append(annotation_field)
                            if annotation_field != annotation_fields_new_name:
                                annotation_list.append(
                                    f"{annotation_fields_new_name}:=INFO/{annotation_field}"
                                )
                                annotation_name_map[annotation_field] = (
                                    f"{annotation_fields_new_name}:=INFO/{annotation_field}"
                                )
                                annotation_list_to_remove.append(f"INFO/{annotation_field}")
                            else:
                                annotation_list.append(annotation_field)
                                annotation_name_map[annotation_field] = annotation_field

                            nb_annotation_field += 1

                        else:

                            if annotation_field not in db_hdr_vcf.get_header().infos:
                                log.warning(
                                    f"Annotation '{annotation}' - '{annotation_field}' - not available in vcf/bed file"
                                )
                            if annotation_fields_new_name in self.get_header().infos:
                                log.warning(
                                    f"Annotation '{annotation}' - '{annotation_fields_new_name}' - already exists (skipped)"
                                )

                    log.info(
                        f"Annotation '{annotation_name}' - {nb_annotation_field} annotations available in vcf/bed file"
                    )

                    if db_file_type in ["bed"]:
                        # For bed file, we need to add CHROM, POS, POS in annotation list for bcftools annotate with bed file, and position all in order of columns in BED file (CHROM, POS, POS, then annotations), and with unwanted annotation as "-"
                        annotation_list_bed = []
                        for col in db_hdr_vcf.get_header_columns_as_list()[3:]:
                            if col in annotation_name_map:
                                annotation_list_bed.append(
                                    annotation_name_map.get(col, col)
                                )
                            else:
                                annotation_list_bed.append("-")
                        annotation_list = annotation_list_bed

                    annotation_infos = ",".join(annotation_list)
                    annotation_infos_to_remove = ",".join(annotation_list_to_remove)

                    if annotation_infos != "":

                        # Protect header for bcftools (remove "#CHROM" and variants line)
                        log.debug("Protect Header file - remove #CHROM line if exists")
                        tmp_header_vcf = NamedTemporaryFile(
                            prefix=self.get_prefix(),
                            dir=self.get_tmp_dir(),
                            suffix=".hdr",
                            delete=False,
                        )
                        tmp_header_vcf_name = tmp_header_vcf.name
                        tmp_files.append(tmp_header_vcf_name)
                        # Command
                        if db_hdr_file.endswith(".gz"):
                            command_extract_header = f"zcat < {db_hdr_file} | grep '^##' > {tmp_header_vcf_name}"
                        else:
                            command_extract_header = f"cat < {db_hdr_file} | grep '^##' > {tmp_header_vcf_name}"
                        # Run
                        run_parallel_commands([command_extract_header], 1)

                        # Find chomosomes
                        log.debug("Find chromosomes ")
                        sql_query_chromosomes = f"""SELECT table_variants.\"#CHROM\" as CHROM FROM {table_variants} as table_variants GROUP BY table_variants.\"#CHROM\""""
                        sql_query_chromosomes_df = self.get_query_to_df(
                            sql_query_chromosomes
                        )
                        chomosomes_list = list(sql_query_chromosomes_df["CHROM"])

                        log.debug("Chromosomes found: " + str(list(chomosomes_list)))

                        # BED columns in the annotation file
                        if db_file_type in ["bed"]:
                            annotation_infos = "CHROM,POS,POS," + annotation_infos

                        for chrom in chomosomes_list:

                            # Create BED on initial VCF
                            log.debug("Create BED on initial VCF: " + str(tmp_vcf_name))
                            tmp_bed = NamedTemporaryFile(
                                prefix=self.get_prefix(),
                                dir=self.get_tmp_dir(),
                                suffix=".bed",
                                delete=False,
                            )
                            tmp_bed_name = tmp_bed.name
                            tmp_files.append(tmp_bed_name)

                            # Remove files
                            remove_files.append(tmp_bed_name)
                            remove_files.append(tmp_bed_name + ".tbi")

                            # Detecte regions
                            log.debug(
                                f"Annotation '{annotation}' - Chromosome '{chrom}' - Start detecting regions..."
                            )
                            window = 1000000
                            sql_query_intervals_for_bed = f"""
                                SELECT  \"#CHROM\",
                                        CASE WHEN \"POS\"-{window}-1 < 0 THEN 0 ELSE \"POS\"-{window}-1 END,
                                        \"POS\"+{window}
                                FROM {table_variants} as table_variants
                                WHERE table_variants.\"#CHROM\" = '{chrom}'
                            """
                            regions = self.conn.execute(
                                sql_query_intervals_for_bed
                            ).fetchall()
                            merged_regions = merge_regions(regions)
                            log.debug(
                                f"Annotation '{annotation}' - Chromosome '{chrom}' - Stop detecting regions..."
                            )

                            header = ["#CHROM", "START", "END"]
                            with open(tmp_bed_name, "w") as f:
                                # Write the header with tab delimiter
                                f.write("\t".join(header) + "\n")
                                for d in merged_regions:
                                    # Write each data row with tab delimiter
                                    f.write("\t".join(map(str, d)) + "\n")

                            # Tmp files
                            tmp_annotation_vcf = NamedTemporaryFile(
                                prefix=self.get_prefix(),
                                dir=self.get_tmp_dir(),
                                suffix=".vcf.gz",
                                delete=False,
                            )
                            tmp_annotation_vcf_name = tmp_annotation_vcf.name
                            tmp_files.append(tmp_annotation_vcf_name)
                            tmp_ann_vcf_list.append(tmp_annotation_vcf_name)
                            tmp_annotation_vcf_name_err = (
                                tmp_annotation_vcf_name + ".err"
                            )
                            err_files.append(tmp_annotation_vcf_name_err)

                            # Remove files
                            remove_files.append(tmp_annotation_vcf_name)
                            remove_files.append(tmp_annotation_vcf_name + ".tbi")
                            remove_files.append(tmp_annotation_vcf_name_err)

                            # Annotate Command
                            log.debug(
                                f"Annotation '{annotation}' - add bcftools command"
                            )

                            # Remove INFO fields command
                            if annotation_infos_to_remove:
                                command_remove_infos_fields = f" | {bcftools_bin_command} annotate -x {annotation_infos_to_remove} "
                            else:
                                command_remove_infos_fields = ""

                            # Command
                            command_annotate = f"{bcftools_bin_command} annotate --pair-logic exact --regions-file={tmp_bed_name} -a {db_file} -h {tmp_header_vcf_name} -c {annotation_infos} {tmp_vcf_name} {command_remove_infos_fields} -o {tmp_annotation_vcf_name} -Oz1 2>>{tmp_annotation_vcf_name_err} && tabix {tmp_annotation_vcf_name} 2>>{tmp_annotation_vcf_name_err} "

                            # Add command
                            commands.append(command_annotate)

            # if some commands
            if commands:

                # Export VCF file
                self.export_variant_vcf(
                    vcf_file=tmp_vcf_name,
                    remove_info=True,
                    add_samples=False,
                    index=True,
                    chrom_mapping_sql=chrom_mapping.to_tool_sql()
                )

                # Remove files
                remove_files.append(tmp_vcf_name)
                remove_files.append(tmp_vcf_name + ".tbi")

                # Threads
                # calculate threads for annotated commands
                if commands:
                    threads_bcftools_annotate = round(threads / len(commands))
                else:
                    threads_bcftools_annotate = 1

                if not threads_bcftools_annotate:
                    threads_bcftools_annotate = 1

                # Add threads option to bcftools commands
                if threads_bcftools_annotate > 1:
                    commands_threaded = []
                    for command in commands:
                        commands_threaded.append(
                            command.replace(
                                f"{bcftools_bin_command} annotate ",
                                f"{bcftools_bin_command} annotate --threads={threads_bcftools_annotate} ",
                            )
                        )
                    commands = commands_threaded

                # Command annotation multithreading
                log.debug(f"Annotation - Annotation commands: " + str(commands))
                log.info(
                    f"Annotation - Annotation multithreaded in "
                    + str(len(commands))
                    + " commands"
                )

                run_parallel_commands(commands, threads)

                # Merge
                tmp_ann_vcf_list_cmd = " ".join(tmp_ann_vcf_list)

                if tmp_ann_vcf_list_cmd:

                    # Tmp file
                    tmp_annotate_vcf = NamedTemporaryFile(
                        prefix=self.get_prefix(),
                        dir=self.get_tmp_dir(),
                        suffix=".vcf.gz",
                        delete=True,
                    )
                    tmp_annotate_vcf_name = tmp_annotate_vcf.name
                    tmp_annotate_vcf_name_err = tmp_annotate_vcf_name + ".err"
                    err_files.append(tmp_annotate_vcf_name_err)

                    # Remove files
                    remove_files.append(tmp_annotate_vcf_name + ".tbi")
                    remove_files.append(tmp_annotate_vcf_name_err)

                    # Tmp file remove command
                    tmp_files_remove_command = ""
                    if tmp_files:
                        tmp_files_remove_command = " && rm -f " + " ".join(tmp_files)

                    # Command merge
                    merge_command = f"{bcftools_bin_command} merge --force-samples --threads={threads} {tmp_vcf_name} {tmp_ann_vcf_list_cmd} -o {tmp_annotate_vcf_name} -Oz 2>>{tmp_annotate_vcf_name_err} {tmp_files_remove_command}"
                    log.info(
                        "Annotation - Annotation merging "
                        + str(len(commands))
                        + " annotated files"
                    )
                    log.debug(f"Annotation - merge command: {merge_command}")
                    run_parallel_commands([merge_command], 1)

                    # Error messages
                    error_message_command_all = []
                    error_message_command_warning = []
                    error_message_command_err = []
                    for err_file in err_files:
                        with open(err_file, "r") as f:
                            for line in f:
                                message = line.strip()
                                error_message_command_all.append(message)
                                if line.startswith("[W::"):
                                    error_message_command_warning.append(message)
                                if line.startswith("[E::"):
                                    error_message_command_err.append(
                                        f"{err_file}: " + message
                                    )

                    if len(error_message_command_err):
                        log.error(f"Error messages:")
                        for message in list(set(error_message_command_err)):
                            log.error(f"   {message}")
                    elif len(error_message_command_warning):
                        log.warning(f"Warning messages:")
                        for message in list(set(error_message_command_warning)):
                            log.warning(f"   {message}")
                    # debug info
                    log.debug(f"Warning/Error messages:")
                    for message in list(set(error_message_command_all)):
                        log.debug(f"   {message}")
                    # failed
                    if len(error_message_command_err):
                        log.error("Annotation failed: Error in commands")
                        raise ValueError("Annotation failed: Error in commands")

                    # Update variants
                    log.info("Annotation - Updating...")
                    self.update_from_vcf(
                        tmp_annotate_vcf_name,
                        remove_vcf_file=False,
                        update_header=True,
                        annotation_header_fields_override=annotation_header_fields_override,
                        chrom_mapping_sql=chrom_mapping.from_tool_sql()
                    )

        # Remove files
        remove_if_exists(remove_files)

    def annotation_exomiser(self, section:str = "annotation", threads: int = None) -> None:
        """
        This function annotate with Exomiser

        This function uses args as parameters, in section "annotation" -> "exomiser", with sections:
        - "analysis" (dict/file):
            Full analysis dictionnary parameters (see Exomiser docs).
            Either a dict, or a file in JSON or YAML format.
            These parameters may change depending on other parameters (e.g. phenotipicFeatures/HPO)
            Default : None
        - "preset" (string):
            Analysis preset (available in config folder).
            Used if no full "analysis" is provided.
            Default: "exome"
        - "phenopacket" (dict/file):
            Samples and phenotipic features parameters (see Exomiser docs).
            Either a dict, or a file in JSON or YAML format.
            Default: None
        - "subject" (dict):
            Sample parameters (see Exomiser docs).
            Example:
                "subject":
                    {
                        "id": "ISDBM322017",
                        "sex": "FEMALE"
                    }
            Default: None
        - "sample" (string):
            Sample name to construct "subject" section:
                "subject":
                    {
                        "id": "<sample>",
                        "sex": "UNKNOWN_SEX"
                    }
            Default: None
        - "phenotypicFeatures" (dict)
            Phenotypic features to construct "subject" section.
            Example:
                "phenotypicFeatures":
                    [
                        { "type": { "id": "HP:0001159", "label": "Syndactyly" } },
                        { "type": { "id": "HP:0000486", "label": "Strabismus" } }
                    ]
        - "hpo" (list)
            List of HPO ids as phenotypic features.
            Example:
                "hpo": ['0001156', '0001363', '0011304', '0010055']
            Default: []
        - "outputOptions" (dict):
            Output options (see Exomiser docs).
            Default:
                "output_options" =
                    {
                        "outputContributingVariantsOnly": False,
                        "numGenes": 0,
                        "outputFormats": ["TSV_VARIANT", "VCF"]
                    }
        - "transcript_source" (string):
            Transcript source (either "refseq", "ucsc", "ensembl")
            Default: "refseq"
        - "exomiser_to_info" (boolean):
            Add exomiser TSV file columns as INFO fields in VCF.
            Default: False
        - "release" (string):
            Exomise database release.
            If not exists, database release will be downloaded (take a while).
            Default: None (provided by application.properties configuration file)
        - "exomiser_application_properties" (file):
            Exomiser configuration file (see Exomiser docs).
            Useful to automatically download databases (especially for specific genome databases).

        Notes:
        - If no sample in parameters, first sample in VCF will be chosen
        - If no HPO found, "hiPhivePrioritiser" analysis step will be switch off

        :param threads: The number of threads to use
        :return: None.
        """

        # DEBUG
        log.debug("Start annotation with Exomiser databases")

        # Threads
        if not threads:
            threads = self.get_threads()
        log.debug("Threads: " + str(threads))

        # Config
        config = self.get_config()
        log.debug("Config: " + str(config))

        # Config - Folders - Databases
        databases_folders = (
            config.get("folders", {})
            .get("databases", {})
            .get("exomiser", f"{DEFAULT_DATABASE_FOLDER}/exomiser/current")
        )
        databases_folders = full_path(databases_folders)
        if not os.path.exists(databases_folders):
            log.error(f"Databases annotations: {databases_folders} NOT found")
        log.debug("Databases annotations: " + str(databases_folders))

        # Config - Exomiser
        exomiser_bin_command = get_bin_command(
            bin="exomiser-cli*.jar",
            tool="exomiser",
            bin_type="jar",
            config=config,
            default_folder=f"{DEFAULT_TOOLS_FOLDER}/exomiser",
        )
        log.debug("Exomiser bin command: " + str(exomiser_bin_command))
        if not exomiser_bin_command:
            msg_err = f"Annotation failed: no exomiser bin '{exomiser_bin_command}'"
            log.error(msg_err)
            raise ValueError(msg_err)

        # Param
        param = self.get_param()
        log.debug("Param: " + str(param))

        # Param - Exomiser
        param_exomiser = param.get(section, {}).get("exomiser", {})
        log.debug(f"Param Exomiser: {param_exomiser}")

        # Param - options chrom_mapping
        chrom_mapping_options = param.get(section, {}).get("exomiser", {}).get("chrom_mapping", None)
        chrom_mapping = ChromMapping(
            mapping=chrom_mapping_options
        )

        # Param - Assembly
        assembly = self.get_assembly()
        log.debug("Assembly: " + str(assembly))

        # Data
        table_variants = self.get_table_variants()

        # Check if not empty
        log.debug("Check if not empty")
        sql_query_chromosomes = (
            f"""SELECT count(*) as count FROM {table_variants} as table_variants"""
        )
        if not self.get_query_to_df(f"{sql_query_chromosomes}")["count"][0]:
            log.info(f"VCF empty")
            return False

        # VCF header
        vcf_reader = self.get_header()
        log.debug("Initial header: " + str(vcf_reader.infos))

        # Samples
        samples = self.get_header_sample_list()
        if not samples:
            log.error("No Samples in VCF")
            return False
        log.debug(f"Samples: {samples}")

        # Memory limit
        memory_limit = self.get_memory("8G")
        log.debug(f"memory_limit: {memory_limit}")

        # Exomiser java options
        exomiser_java_options = (
            f" -Xmx{memory_limit} -XX:+UseParallelGC -XX:ParallelGCThreads={threads} "
        )
        log.debug(f"Exomiser java options: {exomiser_java_options}")

        # Download Exomiser (if not exists)
        exomiser_release = param_exomiser.get("release", None)
        exomiser_application_properties = param_exomiser.get(
            "exomiser_application_properties", None
        )
        databases_download_exomiser(
            assemblies=[assembly],
            exomiser_folder=databases_folders,
            exomiser_release=exomiser_release,
            exomiser_phenotype_release=exomiser_release,
            exomiser_application_properties=exomiser_application_properties,
        )

        # Force annotation
        force_update_annotation = True

        if "Exomiser" not in self.get_header().infos or force_update_annotation:
            log.debug("Start annotation Exomiser")

            with TemporaryDirectory(dir=self.get_tmp_dir()) as tmp_dir:

                ### ANALYSIS ###
                ################

                # Create analysis.json through analysis dict
                # either analysis in param or by default
                # depending on preset exome/genome)

                # Init analysis dict
                param_exomiser_analysis_dict = {}

                # analysis from param
                param_exomiser_analysis = param_exomiser.get("analysis", {})
                param_exomiser_analysis = full_path(param_exomiser_analysis)

                # If analysis in param -> load anlaysis json
                if param_exomiser_analysis:

                    # If param analysis is a file and exists
                    if isinstance(param_exomiser_analysis, str) and os.path.exists(
                        param_exomiser_analysis
                    ):
                        # Load analysis file into analysis dict (either yaml or json)
                        with open(param_exomiser_analysis) as json_file:
                            param_exomiser_analysis_dict = yaml.safe_load(json_file)

                    # If param analysis is a dict
                    elif isinstance(param_exomiser_analysis, dict):
                        # Load analysis dict into analysis dict (either yaml or json)
                        param_exomiser_analysis_dict = param_exomiser_analysis

                    # Error analysis type
                    else:
                        log.error(f"Analysis type unknown. Check param file.")
                        raise ValueError(f"Analysis type unknown. Check param file.")

                # Case no input analysis config file/dict
                # Use preset (exome/genome) to open default config file
                if not param_exomiser_analysis_dict:

                    # default preset
                    default_preset = "exome"

                    # Get param preset or default preset
                    param_exomiser_preset = param_exomiser.get("preset", default_preset)

                    # Try to find if preset is a file
                    if os.path.exists(param_exomiser_preset):
                        # Preset file is provided in full path
                        param_exomiser_analysis_default_config_file = (
                            param_exomiser_preset
                        )
                    # elif os.path.exists(full_path(param_exomiser_preset)):
                    #     # Preset file is provided in full path
                    #     param_exomiser_analysis_default_config_file = full_path(param_exomiser_preset)
                    elif os.path.exists(
                        os.path.join(folder_config, param_exomiser_preset)
                    ):
                        # Preset file is provided a basename in config folder (can be a path with subfolders)
                        param_exomiser_analysis_default_config_file = os.path.join(
                            folder_config, param_exomiser_preset
                        )
                    else:
                        # Construct preset file
                        param_exomiser_analysis_default_config_file = os.path.join(
                            folder_config,
                            f"preset-{param_exomiser_preset}-analysis.json",
                        )

                    # If preset file exists
                    param_exomiser_analysis_default_config_file = full_path(
                        param_exomiser_analysis_default_config_file
                    )
                    if os.path.exists(param_exomiser_analysis_default_config_file):
                        # Load prest file into analysis dict (either yaml or json)
                        with open(
                            param_exomiser_analysis_default_config_file
                        ) as json_file:
                            param_exomiser_analysis_dict["analysis"] = yaml.safe_load(
                                json_file
                            )

                    # Error preset file
                    else:
                        log.error(
                            f"No analysis preset config file ({param_exomiser_analysis_default_config_file})"
                        )
                        raise ValueError(
                            f"No analysis preset config file ({param_exomiser_analysis_default_config_file})"
                        )

                # If no analysis dict created
                if not param_exomiser_analysis_dict:
                    log.error(f"No analysis config")
                    raise ValueError(f"No analysis config")

                # Log
                log.debug(f"Pre analysis dict: {param_exomiser_analysis_dict}")

                ### PHENOPACKET ###
                ###################

                # If no PhenoPacket in analysis dict -> check in param
                if "phenopacket" not in param_exomiser_analysis_dict:

                    # If PhenoPacket in param -> load anlaysis json
                    if param_exomiser.get("phenopacket", None):

                        param_exomiser_phenopacket = param_exomiser.get("phenopacket")
                        param_exomiser_phenopacket = full_path(
                            param_exomiser_phenopacket
                        )

                        # If param phenopacket is a file and exists
                        if isinstance(
                            param_exomiser_phenopacket, str
                        ) and os.path.exists(param_exomiser_phenopacket):
                            # Load phenopacket file into analysis dict (either yaml or json)
                            with open(param_exomiser_phenopacket) as json_file:
                                param_exomiser_analysis_dict["phenopacket"] = (
                                    yaml.safe_load(json_file)
                                )

                        # If param phenopacket is a dict
                        elif isinstance(param_exomiser_phenopacket, dict):
                            # Load phenopacket dict into analysis dict (either yaml or json)
                            param_exomiser_analysis_dict["phenopacket"] = (
                                param_exomiser_phenopacket
                            )

                        # Error phenopacket type
                        else:
                            log.error(f"Phenopacket type unknown. Check param file.")
                            raise ValueError(
                                f"Phenopacket type unknown. Check param file."
                            )

                # If no PhenoPacket in analysis dict -> construct from sample and HPO in param
                if "phenopacket" not in param_exomiser_analysis_dict:

                    # Init PhenoPacket
                    param_exomiser_analysis_dict["phenopacket"] = {
                        "id": "analysis",
                        "proband": {},
                    }

                    ### Add subject ###

                    # If subject exists
                    param_exomiser_subject = param_exomiser.get("subject", {})

                    # If subject not exists -> found sample ID
                    if not param_exomiser_subject:

                        # Found sample ID in param
                        sample = param_exomiser.get("sample", None)

                        # Find sample ID (first sample)
                        if not sample:
                            sample_list = self.get_header_sample_list()
                            if len(sample_list) > 0:
                                sample = sample_list[0]
                            else:
                                log.error(f"No sample found")
                                raise ValueError(f"No sample found")

                        # Create subject
                        param_exomiser_subject = {"id": sample, "sex": "UNKNOWN_SEX"}

                    # Add to dict
                    param_exomiser_analysis_dict["phenopacket"][
                        "subject"
                    ] = param_exomiser_subject

                    ### Add "phenotypicFeatures" ###

                    # If phenotypicFeatures exists
                    param_exomiser_phenotypicfeatures = param_exomiser.get(
                        "phenotypicFeatures", []
                    )

                    # If phenotypicFeatures not exists -> Try to infer from hpo list
                    if not param_exomiser_phenotypicfeatures:

                        # Found HPO in param
                        param_exomiser_hpo = param_exomiser.get("hpo", [])

                        # Split HPO if list in string format separated by comma
                        if isinstance(param_exomiser_hpo, str):
                            param_exomiser_hpo = param_exomiser_hpo.split(",")

                        # Create HPO list
                        for hpo in param_exomiser_hpo:
                            hpo_clean = re.sub("[^0-9]", "", hpo)
                            param_exomiser_phenotypicfeatures.append(
                                {
                                    "type": {
                                        "id": f"HP:{hpo_clean}",
                                        "label": f"HP:{hpo_clean}",
                                    }
                                }
                            )

                    # Add to dict
                    param_exomiser_analysis_dict["phenopacket"][
                        "phenotypicFeatures"
                    ] = param_exomiser_phenotypicfeatures

                    # If phenotypicFeatures not exists -> Remove hiPhivePrioritiser step
                    if not param_exomiser_phenotypicfeatures:
                        for step in param_exomiser_analysis_dict.get(
                            "analysis", {}
                        ).get("steps", []):
                            if "hiPhivePrioritiser" in step:
                                param_exomiser_analysis_dict.get("analysis", {}).get(
                                    "steps", []
                                ).remove(step)

                ### Add Input File ###

                # Initial file name and htsFiles
                tmp_vcf_name = os.path.join(tmp_dir, "initial.vcf.gz")
                param_exomiser_analysis_dict["phenopacket"]["htsFiles"] = [
                    {
                        "uri": tmp_vcf_name,
                        "htsFormat": "VCF",
                        "genomeAssembly": assembly,
                    }
                ]

                ### Add metaData ###

                # If metaData not in analysis dict
                if "metaData" not in param_exomiser_analysis_dict:
                    param_exomiser_analysis_dict["phenopacket"]["metaData"] = {
                        "created": f"{datetime.datetime.now()}".replace(" ", "T") + "Z",
                        "createdBy": "howard",
                        "phenopacketSchemaVersion": 1,
                    }

                ### OutputOptions ###

                # Init output result folder
                output_results = os.path.join(tmp_dir, "results")

                # If no outputOptions in analysis dict
                if "outputOptions" not in param_exomiser_analysis_dict:

                    # default output formats
                    defaut_output_formats = ["TSV_VARIANT", "VCF"]

                    # Get outputOptions in param
                    output_options = param_exomiser.get("outputOptions", None)

                    # If no output_options in param -> check
                    if not output_options:
                        output_options = {
                            "outputContributingVariantsOnly": False,
                            "numGenes": 0,
                            "outputFormats": defaut_output_formats,
                        }

                    # Replace outputDirectory in output options
                    output_options["outputDirectory"] = output_results
                    output_options["outputFileName"] = "howard"

                    # Add outputOptions in analysis dict
                    param_exomiser_analysis_dict["outputOptions"] = output_options

                else:

                    # Replace output_results and output format (if exists in param)
                    param_exomiser_analysis_dict["outputOptions"][
                        "outputDirectory"
                    ] = output_results
                    param_exomiser_analysis_dict["outputOptions"]["outputFormats"] = (
                        list(
                            set(
                                param_exomiser_analysis_dict.get(
                                    "outputOptions", {}
                                ).get("outputFormats", [])
                                + ["TSV_VARIANT", "VCF"]
                            )
                        )
                    )

                # log
                log.debug(f"Pre analysis dict: {param_exomiser_analysis_dict}")

                ### ANALYSIS FILE ###
                #####################

                ### Full JSON analysis config file ###

                exomiser_analysis = os.path.join(tmp_dir, "analysis.json")
                with open(exomiser_analysis, "w") as fp:
                    json.dump(param_exomiser_analysis_dict, fp, indent=4)

                ### SPLIT analysis and sample config files

                # Splitted analysis dict
                param_exomiser_analysis_dict_for_split = (
                    param_exomiser_analysis_dict.copy()
                )

                # Phenopacket JSON file
                exomiser_analysis_phenopacket = os.path.join(
                    tmp_dir, "analysis_phenopacket.json"
                )
                with open(exomiser_analysis_phenopacket, "w") as fp:
                    json.dump(
                        param_exomiser_analysis_dict_for_split.get("phenopacket"),
                        fp,
                        indent=4,
                    )

                # Analysis JSON file without Phenopacket parameters
                param_exomiser_analysis_dict_for_split.pop("phenopacket")
                exomiser_analysis_analysis = os.path.join(
                    tmp_dir, "analysis_analysis.json"
                )
                with open(exomiser_analysis_analysis, "w") as fp:
                    json.dump(param_exomiser_analysis_dict_for_split, fp, indent=4)

                ### INITAL VCF file ###
                #######################

                ### Create list of samples to use and include inti initial VCF file ####

                # Subject (main sample)
                # Get sample ID in analysis dict
                sample_subject = (
                    param_exomiser_analysis_dict.get("phenopacket", {})
                    .get("subject", {})
                    .get("id", None)
                )
                sample_proband = (
                    param_exomiser_analysis_dict.get("phenopacket", {})
                    .get("proband", {})
                    .get("subject", {})
                    .get("id", None)
                )
                sample = []
                if sample_subject:
                    sample.append(sample_subject)
                if sample_proband:
                    sample.append(sample_proband)

                # Get sample ID within Pedigree
                pedigree_persons_list = (
                    param_exomiser_analysis_dict.get("phenopacket", {})
                    .get("pedigree", {})
                    .get("persons", {})
                )

                # Create list with all sample ID in pedigree (if exists)
                pedigree_persons = []
                for person in pedigree_persons_list:
                    pedigree_persons.append(person.get("individualId"))

                # Concat subject sample ID and samples ID in pedigreesamples
                samples = list(set(sample + pedigree_persons))

                # Check if sample list is not empty
                if not samples:
                    log.error(f"No samples found")
                    raise ValueError(f"No samples found")

                # Create VCF with sample (either sample in param or first one by default)
                # Export VCF file
                self.export_variant_vcf(
                    vcf_file=tmp_vcf_name,
                    remove_info=True,
                    add_samples=True,
                    list_samples=samples,
                    index=False,
                    chrom_mapping_sql=chrom_mapping.to_tool_sql()
                )

                ### Execute Exomiser ###
                ########################

                # Init command
                exomiser_command = ""

                # Command exomiser options
                exomiser_options = f" --spring.config.location={databases_folders}/{assembly}/application.properties --exomiser.data-directory={databases_folders}/{assembly} "

                # Release
                exomiser_release = param_exomiser.get("release", None)
                if exomiser_release:
                    # phenotype data version
                    exomiser_options += (
                        f" --exomiser.phenotype.data-version={exomiser_release} "
                    )
                    # data version
                    exomiser_options += (
                        f" --exomiser.{assembly}.data-version={exomiser_release} "
                    )
                    # variant white list
                    variant_white_list_file = (
                        f"{exomiser_release}_{assembly}_clinvar_whitelist.tsv.gz"
                    )
                    if os.path.exists(
                        os.path.join(
                            databases_folders, assembly, variant_white_list_file
                        )
                    ):
                        exomiser_options += f" --exomiser.{assembly}.variant-white-list-path={variant_white_list_file} "

                # transcript_source
                transcript_source = param_exomiser.get(
                    "transcript_source", None
                )  # ucsc, refseq, ensembl
                if transcript_source:
                    exomiser_options += (
                        f" --exomiser.{assembly}.transcript-source={transcript_source} "
                    )

                # If analysis contain proband param
                if param_exomiser_analysis_dict.get("phenopacket", {}).get(
                    "proband", {}
                ):
                    exomiser_command_analysis = f" {exomiser_bin_command} --analysis={exomiser_analysis_analysis} --sample={exomiser_analysis_phenopacket} {exomiser_options} "

                # If no proband (usually uniq sample)
                else:
                    exomiser_command_analysis = f" {exomiser_bin_command} --analysis={exomiser_analysis} {exomiser_options}"

                # Log
                log.debug(f"exomiser_command_analysis={exomiser_command_analysis}")

                # Log
                log.info("Starting Exomiser annotation...")

                # Run command
                result = subprocess.call(
                    exomiser_command_analysis.split(), stdout=subprocess.PIPE
                )
                if result:
                    log.error("Exomiser command failed")
                    raise ValueError("Exomiser command failed")

                ### RESULTS ###
                ###############

                ### Annotate with TSV fields ###

                # Init result tsv file
                exomiser_to_info = param_exomiser.get("exomiser_to_info", False)

                # Init result tsv file
                output_results_tsv = os.path.join(output_results, "howard.variants.tsv")

                # Parse TSV file and explode columns in INFO field
                if exomiser_to_info and os.path.exists(output_results_tsv):

                    # Log
                    log.debug("Exomiser columns to VCF INFO field")

                    # Retrieve columns and types
                    query = f""" SELECT * FROM read_csv('{output_results_tsv}', auto_detect=True, delim='\t', sample_size=-1) LIMIT 0 """
                    output_results_tsv_df = self.get_query_to_df(query)
                    output_results_tsv_columns = output_results_tsv_df.columns.tolist()

                    # Init concat fields for update
                    sql_query_update_concat_fields = []

                    # Fields to avoid
                    fields_to_avoid = [
                        "CONTIG",
                        "START",
                        "END",
                        "REF",
                        "ALT",
                        "QUAL",
                        "FILTER",
                        "GENOTYPE",
                    ]

                    # List all columns to add into header
                    for header_column in output_results_tsv_columns:

                        # If header column is enable
                        if header_column not in fields_to_avoid:

                            # Header info type
                            header_info_type = "String"
                            header_column_df = output_results_tsv_df[header_column]
                            header_column_df_dtype = header_column_df.dtype
                            if header_column_df_dtype == object:
                                if (
                                    pd.to_numeric(header_column_df, errors="coerce")
                                    .notnull()
                                    .all()
                                ):
                                    header_info_type = "Float"
                            else:
                                header_info_type = "Integer"

                            # Header info
                            characters_to_validate = ["-"]
                            pattern = "[" + "".join(characters_to_validate) + "]"
                            header_info_name = re.sub(
                                pattern,
                                "_",
                                f"Exomiser_{header_column}".replace("#", ""),
                            )
                            header_info_number = "."
                            header_info_description = (
                                f"Exomiser {header_column} annotation"
                            )
                            header_info_source = "Exomiser"
                            header_info_version = "unknown"
                            header_info_code = CODE_TYPE_MAP[header_info_type]
                            vcf_reader.infos[header_info_name] = vcf.parser._Info(
                                header_info_name,
                                header_info_number,
                                header_info_type,
                                header_info_description,
                                header_info_source,
                                header_info_version,
                                header_info_code,
                            )

                            # Add field to add for update to concat fields
                            sql_query_update_concat_fields.append(
                                f"""
                                CASE
                                    WHEN table_parquet."{header_column}" NOT IN ('','.')
                                    THEN concat(
                                        '{header_info_name}=',
                                        table_parquet."{header_column}",
                                        ';'
                                        )

                                    ELSE ''
                                END
                            """
                            )

                    # Update query
                    sql_query_update = f"""
                        UPDATE {table_variants} as table_variants
                            SET INFO = concat(
                                            CASE
                                                WHEN INFO NOT IN ('', '.')
                                                THEN INFO
                                                ELSE ''
                                            END,
                                            CASE
                                                WHEN table_variants.INFO NOT IN ('','.')
                                                THEN ';'
                                                ELSE ''
                                            END,
                                            (
                                            SELECT 
                                                concat(
                                                    {",".join(sql_query_update_concat_fields)}
                                                )
                                            FROM read_csv('{output_results_tsv}', auto_detect=True, delim='\t', sample_size=-1) as table_parquet
                                                    WHERE concat('chr', CAST(table_parquet.\"CONTIG\" AS STRING)) = table_variants.\"#CHROM\"
                                                    AND table_parquet.\"START\" = table_variants.\"POS\"
                                                    AND table_parquet.\"ALT\" = table_variants.\"ALT\"
                                                    AND table_parquet.\"REF\" = table_variants.\"REF\"
                                            )
                                        )
                            ;
                        """

                    # Update
                    self.conn.execute(sql_query_update)

                ### Annotate with VCF INFO field ###

                # Init result VCF file
                output_results_vcf = os.path.join(output_results, "howard.vcf.gz")

                # If VCF exists
                if os.path.exists(output_results_vcf):

                    # Log
                    log.debug("Exomiser result VCF update variants")

                    # Update variants with VCF
                    self.update_from_vcf(
                        output_results_vcf,
                        update_header=True,
                        chrom_mapping_sql=chrom_mapping.from_tool_sql()
                    )

        return True

    def annotation_snpeff(self, section:str = "annotation", threads: int = None) -> None:
        """
        This function annotate with snpEff

        :param threads: The number of threads to use
        :return: the value of the variable "return_value".
        """

        # DEBUG
        log.debug("Start annotation with snpeff databases")

        # Threads
        if not threads:
            threads = self.get_threads()
        log.debug("Threads: " + str(threads))

        # Delete tmp
        delete_tmp = True
        if self.get_config().get("verbosity", "warning") in ["debug"]:
            delete_tmp = False
            log.debug("Delete tmp files/folders: " + str(delete_tmp))

        # Config
        config = self.get_config()
        log.debug("Config: " + str(config))

        # Config - Folders - Databases
        databases_folders = (
            config.get("folders", {}).get("databases", {}).get("snpeff", ["."])
        )
        log.debug("Databases annotations: " + str(databases_folders))

        # Config - snpEff bin command
        snpeff_bin_command = get_bin_command(
            bin="snpEff.jar",
            tool="snpeff",
            bin_type="jar",
            config=config,
            default_folder=f"{DEFAULT_TOOLS_FOLDER}/snpeff",
        )
        if not snpeff_bin_command:
            msg_err = f"Annotation failed: no snpeff bin '{snpeff_bin_command}'"
            log.error(msg_err)
            raise ValueError(msg_err)

        # Config - snpEff databases
        snpeff_databases = (
            config.get("folders", {})
            .get("databases", {})
            .get("snpeff", DEFAULT_SNPEFF_FOLDER)
        )
        snpeff_databases = full_path(snpeff_databases)
        if snpeff_databases is not None and snpeff_databases != "":
            log.debug(f"Create snpEff databases folder")
            if not os.path.exists(snpeff_databases):
                os.makedirs(snpeff_databases)

        # Param
        param = self.get_param()
        log.debug("Param: " + str(param))

        # Param - options chrom_mapping
        chrom_mapping_options = param.get(section, {}).get("snpeff", {}).get("chrom_mapping", None)
        chrom_mapping = ChromMapping(
            mapping=chrom_mapping_options
        )

        # Param - Assembly
        assembly = self.get_assembly()

        # Param - Options
        snpeff_options = (
            param.get(section, {}).get("snpeff", {}).get("options", "")
        )
        snpeff_stats = param.get(section, {}).get("snpeff", {}).get("stats", None)
        snpeff_csvstats = (
            param.get(section, {}).get("snpeff", {}).get("csvStats", None)
        )
        if snpeff_stats:
            snpeff_stats = snpeff_stats.replace("OUTPUT", self.get_output())
            snpeff_stats = full_path(snpeff_stats)
            snpeff_options += f" -stats {snpeff_stats}"
        if snpeff_csvstats:
            snpeff_csvstats = snpeff_csvstats.replace("OUTPUT", self.get_output())
            snpeff_csvstats = full_path(snpeff_csvstats)
            snpeff_options += f" -csvStats {snpeff_csvstats}"

        # Data
        table_variants = self.get_table_variants()

        # Check if not empty
        log.debug("Check if not empty")
        sql_query_chromosomes = (
            f"""SELECT count(*) as count FROM {table_variants} as table_variants"""
        )
        # if not self.conn.execute(f"{sql_query_chromosomes}").df()["count"][0]:
        if not self.get_query_to_df(f"{sql_query_chromosomes}")["count"][0]:
            log.info(f"VCF empty")
            return

        # Export in VCF
        log.debug("Create initial file to annotate")
        tmp_vcf = NamedTemporaryFile(
            prefix=self.get_prefix(),
            dir=self.get_tmp_dir(),
            suffix=".vcf.gz",
            delete=True,
        )
        tmp_vcf_name = tmp_vcf.name

        # VCF header
        vcf_reader = self.get_header()
        log.debug("Initial header: " + str(vcf_reader.infos))

        # Existing annotations
        for vcf_annotation in self.get_header().infos:

            vcf_annotation_line = self.get_header().infos.get(vcf_annotation)
            log.debug(
                f"Existing annotations in VCF: {vcf_annotation} [{vcf_annotation_line}]"
            )

        # Memory limit
        # if config.get("memory", None):
        #     memory_limit = config.get("memory", "8G")
        # else:
        #     memory_limit = "8G"
        memory_limit = self.get_memory("8G")
        log.debug(f"memory_limit: {memory_limit}")

        # snpEff java options
        snpeff_java_options = (
            f" -Xmx{memory_limit} -XX:+UseParallelGC -XX:ParallelGCThreads={threads} "
        )
        log.debug(f"snpEff java options: {snpeff_java_options}")

        force_update_annotation = True

        if "ANN" not in self.get_header().infos or force_update_annotation:

            # Check snpEff database
            log.debug(f"Check snpEff databases {[assembly]}")
            databases_download_snpeff(
                folder=snpeff_databases, assemblies=[assembly], config=config
            )

            # Export VCF file
            #log.debug(f"chrom_mapping.to_tool_sql(): {chrom_mapping.to_tool_sql()}")
            self.export_variant_vcf(
                vcf_file=tmp_vcf_name,
                remove_info=True,
                add_samples=False,
                index=True,
                chrom_mapping_sql=chrom_mapping.to_tool_sql()
            )

            # Tmp file
            err_files = []
            tmp_annotate_vcf = NamedTemporaryFile(
                prefix=self.get_prefix(),
                dir=self.get_tmp_dir(),
                suffix=".vcf",
                delete=False,
            )
            tmp_annotate_vcf_name = tmp_annotate_vcf.name
            tmp_annotate_vcf_name_err = tmp_annotate_vcf_name + ".err"
            err_files.append(tmp_annotate_vcf_name_err)

            # Command
            snpeff_command = f"{snpeff_bin_command} {assembly} -dataDir {snpeff_databases} {snpeff_options} {tmp_vcf_name} 1>{tmp_annotate_vcf_name} 2>>{tmp_annotate_vcf_name_err}"
            log.debug(f"Annotation - snpEff command: {snpeff_command}")
            run_parallel_commands([snpeff_command], 1)

            # Error messages
            error_message_command_all = []
            error_message_command_warning = []
            error_message_command_err = []
            for err_file in err_files:
                with open(err_file, "r") as f:
                    for line in f:
                        message = line.strip()
                        error_message_command_all.append(message)
                        if line.startswith("[W::"):
                            error_message_command_warning.append(message)
                        if line.startswith("[E::"):
                            error_message_command_err.append(f"{err_file}: " + message)

            # Warning/Error messages
            if len(error_message_command_err):
                log.error("Error messages:")
                for message in set(error_message_command_err):
                    log.error(f"   {message}")
            elif len(error_message_command_warning):
                log.warning("Warning messages:")
                for message in set(error_message_command_warning):
                    log.warning(f"   {message}")
            # debug info
            log.debug("Warning/Error messages:")
            for message in set(error_message_command_all):
                log.debug(f"   {message}")
            # failed
            if len(error_message_command_err):
                log.error("Annotation failed: Error in commands")
                raise ValueError("Annotation failed: Error in commands")

            # Update variants
            log.info(f"Annotation - Updating...")
            #log.debug(f"chrom_mapping.from_tool_sql(): {chrom_mapping.from_tool_sql()}")
            self.update_from_vcf(
                tmp_annotate_vcf_name,
                update_header=True,
                chrom_mapping_sql=chrom_mapping.from_tool_sql()
            )
            list_to_remove = [
                tmp_annotate_vcf_name,
                tmp_annotate_vcf_name_err,
                f"{tmp_annotate_vcf_name}.tbi",
                f"{tmp_vcf_name}.tbi",
            ]
            log.debug(f"tmp_annotate_vcf_name: {list_to_remove}")
            remove_if_exists(list_to_remove)

        else:
            if "ANN" in self.get_header().infos:
                log.debug(f"Existing snpEff annotations in VCF")
            if force_update_annotation:
                log.debug(f"Existing snpEff annotations in VCF - annotation forced")


    def annotation_vep(self, section:str = "annotation", threads: int = None) -> None:
            """
            This function annotate with VEP
    
            :param threads: The number of threads to use
            :return: the value of the variable "return_value".
            """
    
            # DEBUG
            log.debug("Start annotation with VEP databases")
    
            # Threads
            if not threads:
                threads = self.get_threads()
            log.debug("Threads: " + str(threads))
    
            # Config
            config = self.get_config()
            log.debug("Config: " + str(config))

            # Param
            param = self.get_param()

            # Assembly
            assembly = self.get_assembly()
    
            # vep_parameters
            vep_parameters = param.get(section, {}).get("vep", {}).get("parameters", None) or None

            # vep_parameters
            vep_tool = param.get(section, {}).get("vep", {}).get("tool", None) or "vep"

            # Construct vep configuraton (on line) if not exist in config
            if config.get("tools", {}).get(vep_tool) is None:
                if "tools" not in config:
                    config["tools"] = {}
                config["tools"][vep_tool] = {
                    "docker": {
                        "image": "ensemblorg/ensembl-vep:latest",
                        "entrypoint": "vep",
                        "parameters": {
                            "primary": {
                                "input": "--input_file",
                                "output": "--output_file",
                                "threads": "--fork",
                                "assembly": {
                                    "flag": "--assembly",
                                    "source": "GRCH"
                                }
                            },
                            "defaults": {
                                "parameters": [
                                    "--vcf",
                                    "--database"
                                ]
                            }
                        }
                    }
                }

            #log.debug(f"Config: {json.dumps(self.get_config(), indent=4)}")

    
            # Construct param dict for docker
            param_original = copy.deepcopy(param)
            param_vep = copy.deepcopy(param)
            param_vep_section = f"{section}_{get_random()}"
            param_vep["assembly"] = assembly
            param_vep[param_vep_section] = {
                "docker": {
                    "entries": {
                        "annotation_vep": {
                            "tool": vep_tool,
                            "parameters": vep_parameters
                        }
                    },
                }
            }

            # Log
            log.debug(f"param_vep: {param_vep}")

            # Set the modified parameters for the VEP annotation
            self.set_param(param_vep)

            # Run the annotation using the docker container for VEP
            log.info("Running VEP annotation with docker...")
            self.annotation_docker(section=param_vep_section)

            # Restore the original parameters after running the VEP annotation
            self.set_param(param_original)


    def annotation_annovar(self, section:str = "annotation", threads: int = None) -> None:
        """
        It takes a VCF file, annotates it with Annovar, and then updates the database with the new
        annotations

        :param threads: number of threads to use
        :return: the value of the variable "return_value".
        """

        # DEBUG
        log.debug("Start annotation with Annovar databases")

        def create_annovar_command(annotations_parameters, annovar_bin_command, tmp_vcf_name, annovar_databases_assembly, assembly, tmp_annotate_vcf_prefix, threads, tmp_annotate_vcf_name_err, annovar_fields_to_keep, tmp_rename_name):
            """
            Create the annovar command based on the provided parameters.
            """

            char_code_map = {
                "\\x3b": ",",
                "\\x3d": ":"
            }

            # Merged annovar command
            log.debug(f"annotations_parameters: {annotations_parameters}")

            # Create temporary folder
            command_annovar = f"""mkdir -p {tmp_annotate_vcf_prefix} """

            # Command - Annovar
            command_annovar += f""" && {annovar_bin_command} {tmp_vcf_name} {annovar_databases_assembly} --buildver {assembly} --outfile {tmp_annotate_vcf_prefix} --remove --protocol {",".join(p.strip() for p in annotations_parameters.get("protocol",[]))} --operation {",".join(o.strip() for o in annotations_parameters.get("operation",[]))} --argument {",".join(f"' {a.strip()} '" for a in annotations_parameters.get("argument",[]))} {" ".join(list(set(o.strip() for o in annotations_parameters.get("options", []))))} --nastring . --vcfinput --polish --dot2underline --thread {threads}  2>>{tmp_annotate_vcf_name_err} && mv {tmp_annotate_vcf_name_annovar} {tmp_annotate_vcf_name}.tmp.vcf """

            # Command - start pipe
            command_annovar += f""" && {bcftools_bin_command} view --threads={threads} {tmp_annotate_vcf_name}.tmp.vcf 2>>{tmp_annotate_vcf_name_err} """

            # Command - Clean INFO/ANNOVAR_DATE and INFO/ALLELE_END (due to Annovar issue with multiple TAGS!)
            command_annovar += """ | sed "s/ANNOVAR_DATE=[^;\t]*;//gi" """
            command_annovar += """ | sed "s/ALLELE_END=[^;\t]*;//gi" """

            # Command - Special characters (refGene annotation)
            for char, replacement in char_code_map.items():
                command_annovar += f""" | sed "s/\\\\{char}/{replacement}/gi" """

            # Command - Clean empty fields (with value ".")
            command_annovar += """ | awk -F'\\t' -v OFS='\\t' '{if ($0 ~ /^#/) print; else {split($8,a,";");for(i=1;i<=length(a);i++) {split(a[i],b,"=");if(b[2]!=".") {c[b[1]]=b[2]}}; split($8,d,";");for(i=1;i<=length(d);i++) {split(d[i],e,"=");if(c[e[1]]!="") {if(info!="") {info=info";"}; info=info""e[1]"="c[e[1]]}}; if(info!="") {$8=info} else {$8=""}; delete c; info=""; print}}' """

            # Command - Keep only needed fields, and remove ANNOVAR fields, and compress and index final file
            annovar_fields_to_keep = set(annovar_fields_to_keep + annotations_parameters.get("annovar_fields_to_keep", []))
            if len(annovar_fields_to_keep):
                annovar_fields_to_keep_option = "-x " + ",".join(annovar_fields_to_keep)
            else:
                annovar_fields_to_keep_option = ""

            # Command - Annotate with bcftools annotate
            command_annovar += f""" | {bcftools_bin_command} annotate --pair-logic exact --threads={threads} {annovar_fields_to_keep_option} --rename-annots={tmp_rename_name} -o {tmp_annotate_vcf_name} -Oz 2>>{tmp_annotate_vcf_name_err} """

            # Command - indexing
            command_annovar += f"""  && tabix {tmp_annotate_vcf_name} """

            return command_annovar
                
        def err_files_process(err_files, merge_msg:bool = False):
            """
            Process error files and log messages.
            """
            # Error messages
            if merge_msg:
                error_message_command_all = []
                error_message_command_warning = []
                error_message_command_err = []
            for err_file in err_files:

                log.debug(f"Process error file: {err_file}")

                if not merge_msg:
                    error_message_command_all = []
                    error_message_command_warning = []
                    error_message_command_err = []

                with open(err_file, "r") as f:
                    for line in f:
                        message = line.strip()
                        error_message_command_all.append(message)
                        if line.startswith("[W::"):
                            error_message_command_warning.append(message)
                        if line.startswith("[E::"):
                            error_message_command_err.append(f"{err_file}: " + message)

                if not merge_msg:
                    # Warning/Error messages
                    if len(error_message_command_err):
                        log.error("Error messages:")
                        for message in error_message_command_err:
                            log.error(f"   {message}")
                    elif len(error_message_command_warning):
                        log.warning("Warning messages:")
                        for message in error_message_command_warning:
                            log.warning(f"   {message}")
                    # debug info
                    log.debug("Warning/Error messages:")
                    for message in error_message_command_all:
                        log.debug(f"   {message}")
                    # failed
                    if len(error_message_command_err):
                        log.error("Annotation failed: Error in commands")
                        raise ValueError("Annotation failed: Error in commands")

            if merge_msg:
                if not merge_msg:
                    # Warning/Error messages
                    if len(error_message_command_err):
                        log.error("Error messages:")
                        for message in set(error_message_command_err):
                            log.error(f"   {message}")
                    elif len(error_message_command_warning):
                        log.warning("Warning messages:")
                        for message in set(error_message_command_warning):
                            log.warning(f"   {message}")
                    # debug info
                    log.debug("Warning/Error messages:")
                    for message in set(error_message_command_all):
                        log.debug(f"   {message}")
                    # failed
                    if len(error_message_command_err):
                        log.error("Annotation failed: Error in commands")
                        raise ValueError("Annotation failed: Error in commands")

        # Threads
        if not threads:
            threads = self.get_threads()
        log.debug("Threads: " + str(threads))

        # Tmp en Err files
        tmp_files = []
        err_files = []

        # DEBUG
        delete_tmp = True
        if self.get_config().get("verbosity", "warning") in ["debug"]:
            delete_tmp = False
            log.debug("Delete tmp files/folders: " + str(delete_tmp))

        # Config
        config = self.get_config()
        log.debug("Config: " + str(config))

        # Config - Folders - Databases
        databases_folders = (
            config.get("folders", {}).get("databases", {}).get("annovar", ["."])
        )
        log.debug("Databases annotations: " + str(databases_folders))

        # Config - annovar bin command
        annovar_bin_command = get_bin_command(
            bin="table_annovar.pl",
            tool="annovar",
            bin_type="perl",
            config=config,
            default_folder=f"{DEFAULT_TOOLS_FOLDER}/annovar",
        )
        if not annovar_bin_command:
            msg_err = f"Annotation failed: no annovar bin '{annovar_bin_command}'"
            log.error(msg_err)
            raise ValueError(msg_err)

        # Config - BCFTools bin command
        bcftools_bin_command = get_bin_command(
            bin="bcftools",
            tool="bcftools",
            bin_type="bin",
            config=config,
            default_folder=f"{DEFAULT_TOOLS_FOLDER}/bcftools",
        )
        if not bcftools_bin_command:
            msg_err = f"Annotation failed: no bcftools bin '{bcftools_bin_command}'"
            log.error(msg_err)
            raise ValueError(msg_err)

        # Config - annovar databases
        annovar_databases = (
            config.get("folders", {})
            .get("databases", {})
            .get("annovar", DEFAULT_ANNOVAR_FOLDER)
        )
        if annovar_databases is not None:
            if isinstance(annovar_databases, list):
                annovar_databases = full_path(annovar_databases[0])
                log.warning(f"Annovar databases folder '{annovar_databases}' selected")
            annovar_databases = full_path(annovar_databases)
            if not os.path.exists(annovar_databases):
                log.info(f"Annovar databases folder '{annovar_databases}' created")
                Path(annovar_databases).mkdir(parents=True, exist_ok=True)
        else:
            msg_err = f"Annovar databases configuration failed"
            log.error(msg_err)
            raise ValueError(msg_err)

        # Param
        param = self.get_param()
        log.debug("Param: " + str(param))

        # Param - options
        options = param.get(section, {}).get("annovar", {}).get("options", {})
        log.debug("Options: " + str(options))

        # Param - annotations
        annotations = (
            param.get(section, {}).get("annovar", {}).get("annotations", {})
        )
        log.debug("Annotations: " + str(annotations))

        # Header fields override (Number/Type/Description), collected from each
        # annotation's own options block
        # (annotation.annovar.annotations.<name>.options.header_fields)
        annotation_header_fields_override = self.get_annotation_header_fields_override(
            annotations=annotations
        )
        log.debug(
            f"annotation_header_fields_override={annotation_header_fields_override}"
        )

        # Param - options chrom_mapping
        chrom_mapping_options = param.get(section, {}).get("annovar", {}).get("chrom_mapping", None)
        chrom_mapping = ChromMapping(
            mapping=chrom_mapping_options
        )

        # Param - Assembly
        assembly = self.get_assembly()

        # Annovar database assembly
        annovar_databases_assembly = f"{annovar_databases}/{assembly}"
        if annovar_databases_assembly != "" and not os.path.exists(
            annovar_databases_assembly
        ):
            os.makedirs(annovar_databases_assembly)

        # Data
        table_variants = self.get_table_variants()

        # Check if not empty
        log.debug("Check if not empty")
        sql_query_chromosomes = (
            f"""SELECT count(*) as count FROM {table_variants} as table_variants"""
        )
        sql_query_chromosomes_df = self.get_query_to_df(sql_query_chromosomes)
        if not sql_query_chromosomes_df["count"][0]:
            log.info(f"VCF empty")
            return

        # VCF header
        vcf_reader = self.get_header()
        log.debug("Initial header: " + str(vcf_reader.infos))

        # Existing annotations
        for vcf_annotation in self.get_header().infos:

            vcf_annotation_line = self.get_header().infos.get(vcf_annotation)
            log.debug(
                f"Existing annotations in VCF: {vcf_annotation} [{vcf_annotation_line}]"
            )

        force_update_annotation = True

        if annotations:

            with TemporaryDirectory(
                prefix=self.get_prefix(), dir=self.get_tmp_dir(), suffix=".annovar"
            ) as tmp_annotate_vcf_directory:

                #commands = []
                tmp_annotates_vcf_name_list = []

                # Export in VCF
                log.debug("Create initial file to annotate")
                tmp_vcf_name = os.path.join(tmp_annotate_vcf_directory, get_random(10) + "annovar.vcf.gz")

                # Export VCF file
                self.export_variant_vcf(
                    vcf_file=tmp_vcf_name,
                    remove_info=".",
                    add_samples=False,
                    index=True,
                    chrom_mapping_sql=chrom_mapping.to_tool_sql()
                )

                # Create file for field rename
                log.debug("Create file for field rename")
                tmp_rename_name = os.path.join(tmp_annotate_vcf_directory, get_random(10) + "annovar.rename")

                # Check Annovar database
                log.debug(
                    f"Check Annovar databases {[assembly]}: {list(annotations.keys())}"
                )
                databases_download_annovar(
                    folder=annovar_databases,
                    files=list(annotations.keys()),
                    assemblies=[assembly],
                    force_check_dblist=False,
                    force_update=False,
                )

                # Annotations parameters
                annotations_parameters = {
                    "protocol": [],
                    "operation": [],
                    "argument": [],
                    "options": [],
                    "annovar_fields_to_keep": []
                }

                # parallelize annovar command
                parallelize_annovar_command = options.get("parallelize", None)
                if (
                    parallelize_annovar_command is not None
                    and parallelize_annovar_command not in ["multianno", "parallel"]
                ):
                    msg_warning = f"Annotation - annovar parallelize option '{parallelize_annovar_command}' not valid (should be 'multianno' or 'parallel')"
                    #log.warning(msg_warning)
                    if parallelize_annovar_command.startswith("multi"):
                        log.warning(
                            f"{msg_warning} - '{parallelize_annovar_command}' is identified as 'multianno' and parallelization is enabled"
                        )
                        parallelize_annovar_command = "multianno"
                    elif parallelize_annovar_command.startswith("paral"):
                        log.warning(
                            f"{msg_warning} - '{parallelize_annovar_command}' is identified as 'parallel' and parallelization is enabled"
                        )
                        parallelize_annovar_command = "parallel"
                    else:
                        log.warning(
                            f"{msg_warning} - Parallelization is disabled"
                        )
                        parallelize_annovar_command = None

                # command_annovar_list
                command_annovar_list = []

                # err files list
                err_files_list = []

                for annotation in annotations:

                    # Allow to add options for annotations
                    if "annotation_fields" in annotations[annotation]:
                        annotation_fields = annotations[annotation].get(
                            "annotation_fields", {}
                        )
                    else:
                        annotation_fields = annotations[annotation]

                    # Options for annotations
                    annotation_options = annotations[annotation].get("options", {})

                    # Check options for annotation

                    # Operation option
                    operation = annotation_options.get("operation", None)
                    if operation and operation not in ["gx", "g", "f", "r"]:
                        log.error(
                            f"Annotation '{annotation}' - operation option '{operation}' not valid (should be 'gx', 'g', 'f' or 'r')"
                        )
                        raise ValueError(
                            f"Annotation '{annotation}' - operation option '{operation}' not valid (should be 'gx', 'g', 'f' or 'r')"
                        )

                    if not annotation_fields:
                        annotation_fields = {"INFO": None}

                    log.info(f"Annotations Annovar - database '{annotation}'")
                    log.debug(f"Annotation '{annotation}' - fields: {annotation_fields}")

                    # Generate tmp files for annovar

                    # Tmp file for annovar
                    err_files = []
                    tmp_annotate_vcf_random = get_random(10)
                    tmp_annotate_vcf_prefix = tmp_annotate_vcf_directory + "/" + tmp_annotate_vcf_random + "/annovar"
                    tmp_annotate_vcf_name_annovar = (
                        tmp_annotate_vcf_prefix + "." + assembly + "_multianno.vcf"
                    )
                    tmp_annotate_vcf_name_err = tmp_annotate_vcf_directory + "/" + tmp_annotate_vcf_random + "/.err"
                    err_files.append(tmp_annotate_vcf_name_err)
                    tmp_files.append(tmp_annotate_vcf_name_err)

                    # Tmp file final vcf annotated by annovar
                    tmp_annotate_vcf_name = os.path.join(tmp_annotate_vcf_directory, get_random(10) + "annovar.vcf.gz")
                    tmp_annotates_vcf_name_list.append(tmp_annotate_vcf_name)
                    tmp_files.append(tmp_annotate_vcf_name)
                    tmp_files.append(tmp_annotate_vcf_name + ".tbi")

                    # Number of fields
                    annotation_list = []
                    annotation_renamed_list = []

                    # Create file for field rename
                    log.debug("Create file for field rename")
                    annotation_tmp_rename_name = os.path.join(tmp_annotate_vcf_directory, get_random(10) + "annovar.rename")

                    for annotation_field in annotation_fields:

                        # Check if field contains space
                        if " " in annotation_field:
                            msg_err = f"Annotation '{annotation}' - field '{annotation_field}' contains space (not valid)"
                            log.error(msg_err)
                            raise ValueError(msg_err)

                        # field new name, if parametered SKIPPED !!!!!! not managed actually TODO
                        annotation_fields_new_name = annotation_fields.get(
                            annotation_field, annotation_field
                        )
                        if not annotation_fields_new_name:
                            annotation_fields_new_name = annotation_field

                        if (
                            force_update_annotation
                            or annotation_fields_new_name not in self.get_header().infos
                        ):
                            annotation_list.append(annotation_field)
                            annotation_renamed_list.append(annotation_fields_new_name)
                        else:  # annotation_fields_new_name in self.get_header().infos and not force_update_annotation:
                            log.warning(
                                f"Annotation '{annotation}' - '{annotation_fields_new_name}' - already exists (skipped)"
                            )

                        # Add rename info to tmp_rename_name
                        with open(tmp_rename_name, "a") as f:
                            f.write(
                                f"INFO/{annotation_field} {annotation_fields_new_name}\n"
                            )

                        # Add rename info to annotation_tmp_rename_name
                        with open(annotation_tmp_rename_name, "a") as f:
                            f.write(
                                f"INFO/{annotation_field} {annotation_fields_new_name}\n"
                            )

                    # log.debug("fields_to_removed: " + str(fields_to_removed))
                    log.debug("annotation_list: " + str(annotation_list))

                    # protocol
                    if annotation_options.get("protocol", None):
                        protocol = annotation_options.get("protocol", None)
                    else:
                        protocol = annotation

                    # operation
                    if operation is None:
                        log.debug(
                            f"Operation NOT found in options for annotation '{annotation}'"
                        )
                        # Operation for gene-based annotation, region-based for cytoBand, and filter-based for others (by default)
                        if annotation.startswith(("refGene", "ensGene", "knownGene")):
                            operation = "g"
                        elif annotation in ["cytoBand"]:
                            operation = "r"
                        else:
                            operation = "f"
                        log.debug(
                            f"Operation detected in for annotation '{annotation}': {operation}"
                        )
                    else:
                        log.debug(
                            f"Operation found in options for annotation '{annotation}': {operation}"
                        )

                    ### Argument
                    annotation_argument = ""

                    # GeneBase argument for gene-based annotation (refGene, ensGene, knownGene)
                    if operation in ["g", "gx"] and (
                        annotation_options.get("genebase", None)
                        or options.get("genebase", None)
                    ):
                        if annotation_options.get("genebase", None):
                            annotation_argument = (
                                f""" {annotation_options.get("genebase","")} """
                            )
                        elif options.get("genebase", None):
                            annotation_argument = f""" {options.get("genebase","")} """

                    # Additional argument
                    if annotation_options.get("arguments", None):
                        if isinstance(annotation_options.get("arguments"), str):
                            annotation_argument += (
                                f""" {annotation_options.get("arguments","")}"""
                            )
                        elif isinstance(annotation_options.get("arguments"), list):
                            annotation_argument += (
                                f""" {" ".join(annotation_options.get("arguments",""))}"""
                            )

                    ### Options
                    annotation_parameters = ""
                    command_options = []

                    # Xref options for cross information for "gx" operation
                    # Example: -xref ~/hpc/tools/annovar/humandb/gene_fullxref.txt
                    #if operation in ["g", "gx"] and (
                    if operation in ["gx"] and (
                        annotation_options.get("xref", None) or options.get("xref", None)
                    ):
                        #operation = "gx"  # Change operation to "gx" if xref is provided
                        if annotation_options.get("xref", None):
                            annotation_parameters += (
                                f""" --xref {annotation_options.get("xref","")}"""
                            )
                            command_options.append(
                                f""" --xref {annotation_options.get("xref","")}"""
                            )
                        elif options.get("xref", None):
                            annotation_parameters += f""" --xref {options.get("xref","")}"""
                            command_options.append(f""" --xref {options.get("xref","")}""")

                    # command options
                    command_options.append(f""" {annotation_parameters} """)
                    command_options.append(f""" {annotation_options.get("options","")} """)

                    # merge dict of options and annotation_options
                    merged_options = {**options, **annotation_options}

                    for option in merged_options:
                        if option in ["xref"]:
                            if merged_options.get(option, None):
                                command_options.append(
                                    f""" --xref {merged_options[option]}"""
                                )
                        elif option in ["options"]:
                            command_options.append(f""" {merged_options[option]}""")
                        elif option not in [
                            "protocol",
                            "genebase",
                            "xref",
                            "arguments",
                            "operation",
                            "options",
                            "parallelize",
                            "header_fields"
                        ]:
                            command_options.append(
                                f""" --{option}={merged_options[option]}"""
                            )

                    # Command parameters - merge and uniquify
                    command_parameters = " ".join(
                        list(set(o.strip() for o in command_options))
                    )
                    log.debug(
                        f"Annotation '{annotation}' - command_parameters: {command_parameters}"
                    )

                    # Extract only needed fields, and remove ANNOVAR fields, and compress and index final file
                    #annovar_fields_to_keep = ["INFO/ANNOVAR_DATE", "INFO/ALLELE_END"]
                    annovar_fields_to_keep = []
                    if "ALL" not in annotation_list and "INFO" not in annotation_list:
                        # for ann in annotation_renamed_list:
                        for ann in annotation_list:
                            annovar_fields_to_keep.append(f"^INFO/{ann}")

                    # Parallelization of annovar command - multianno
                    if parallelize_annovar_command == "multianno":

                        # Add parameters for parallelization
                        annotations_parameters["protocol"].append(protocol)
                        annotations_parameters["operation"].append(operation)
                        annotations_parameters["argument"].append(
                            annotation_argument.strip()
                        )
                        # annotations_parameters["options"].append(command_parameters)
                        annotations_parameters["options"] += command_options
                        annotations_parameters["annovar_fields_to_keep"] += annovar_fields_to_keep
                        annotations_parameters["tmp_rename_name"] = tmp_rename_name

                    # Either parallelize_annovar_command is "parallel" or None
                    else:

                        # Threads by annotation
                        if parallelize_annovar_command == "parallel":
                            threads_by_annotation = int((threads / len(annotations)) + 1)
                        else:
                            threads_by_annotation = threads

                        # Fixed parameters for Annovar
                        command_parameters += f""" --nastring . --vcfinput --polish --dot2underline --thread {threads_by_annotation} """  # --intronhgvs 10

                        # Command

                        annotations_parameters = {
                            "protocol": [protocol],
                            "operation": [operation],
                            "argument": [annotation_argument.strip()],
                            "options": command_options,
                        }

                        # Create annovar command
                        command_annovar = create_annovar_command(
                            annotations_parameters=annotations_parameters,
                            annovar_bin_command=annovar_bin_command,
                            tmp_vcf_name=tmp_vcf_name,
                            annovar_databases_assembly=annovar_databases_assembly,
                            assembly=assembly,
                            tmp_annotate_vcf_prefix=tmp_annotate_vcf_prefix,
                            threads=threads,
                            tmp_annotate_vcf_name_err=tmp_annotate_vcf_name_err,
                            annovar_fields_to_keep=annovar_fields_to_keep,
                            tmp_rename_name=annotation_tmp_rename_name
                        )

                        # Launch command
                        log.debug(f"Annotation - Annovar command: {command_annovar}")

                        # If parallelization of annovar command is "parallel", add command to list, else run command
                        if parallelize_annovar_command == "parallel":

                            # Append command to list
                            command_annovar_list.append(command_annovar)

                            err_files_list += err_files

                        else:

                            # Run command
                            log.info(
                                f"Annotations Annovar - database '{annotation}' - annotating..."
                            )
                            run_parallel_commands([command_annovar], 1)

                            # Error messages
                            error_message_command_all = []
                            error_message_command_warning = []
                            error_message_command_err = []
                            for err_file in err_files:
                                with open(err_file, "r") as f:
                                    for line in f:
                                        message = line.strip()
                                        error_message_command_all.append(message)
                                        if line.startswith("[W::"):
                                            error_message_command_warning.append(message)
                                        if line.startswith("[E::"):
                                            error_message_command_err.append(
                                                f"{err_file}: " + message
                                            )

                            # Error/Warning messages
                            if len(error_message_command_err):
                                log.error(f"Error messages:")
                                for message in list(set(error_message_command_err)):
                                    log.error(f"   {message}")
                            elif len(error_message_command_warning):
                                log.warning(f"Warning messages:")
                                for message in list(set(error_message_command_warning)):
                                    log.warning(f"   {message}")
                            # debug info
                            log.debug(f"Warning/Error messages:")
                            for message in list(set(error_message_command_all)):
                                log.debug(f"   {message}")
                            # failed
                            if len(error_message_command_err):
                                log.error("Annotation failed: Error in commands")
                                raise ValueError("Annotation failed: Error in commands")

                if parallelize_annovar_command == "multianno":

                    # Generate tmp files for annovar
                    log.debug("Generate temporary files for annovar")

                    # Tmp file for annovar
                    err_files = []
                    tmp_annotate_vcf_prefix = tmp_annotate_vcf_directory + "/annovar"
                    tmp_annotate_vcf_name_annovar = (
                        tmp_annotate_vcf_prefix + "." + assembly + "_multianno.vcf"
                    )
                    tmp_annotate_vcf_name_err = tmp_annotate_vcf_directory + "/.err"
                    err_files.append(tmp_annotate_vcf_name_err)
                    tmp_files.append(tmp_annotate_vcf_name_err)

                    # Tmp file final vcf annotated by annovar
                    tmp_annotate_vcf_name = os.path.join(tmp_annotate_vcf_directory, get_random(10) + "annovar.vcf.gz")
                    tmp_annotates_vcf_name_list.append(tmp_annotate_vcf_name)
                    tmp_files.append(tmp_annotate_vcf_name)
                    tmp_files.append(tmp_annotate_vcf_name + ".tbi")

                    # Merged annovar command
                    log.debug(f"annotations_parameters: {annotations_parameters}")

                    # Create annovar command
                    command_annovar = create_annovar_command(
                        annotations_parameters=annotations_parameters,
                        annovar_bin_command=annovar_bin_command,
                        tmp_vcf_name=tmp_vcf_name,
                        annovar_databases_assembly=annovar_databases_assembly,
                        assembly=assembly,
                        tmp_annotate_vcf_prefix=tmp_annotate_vcf_prefix,
                        threads=threads,
                        tmp_annotate_vcf_name_err=tmp_annotate_vcf_name_err,
                        #annovar_fields_to_keep=annovar_fields_to_keep,
                        annovar_fields_to_keep=[], # annovar_fields_to_keep in annotations_parameters
                        tmp_rename_name=tmp_rename_name
                    )

                    # Launch command
                    log.info(f"Annotations Annovar - annotating...")
                    log.debug(f"Annotation - Annovar command: {command_annovar}")
                    run_parallel_commands([command_annovar], 1)
                    tmp_annotates_vcf_name_list = [tmp_annotate_vcf_name]

                    err_files_process(err_files=err_files, merge_msg=False)

                elif parallelize_annovar_command == "parallel" and command_annovar_list:
                    log.info(
                        f"Annotations Annovar - annotating..."
                    )
                    run_parallel_commands(command_annovar_list, threads)

                    err_files_process(err_files=err_files_list)

                if tmp_annotates_vcf_name_list:

                    # Tmp file
                    tmp_annotate_vcf_name = os.path.join(tmp_annotate_vcf_directory, get_random(10) + "annovar.vcf.gz")
                    tmp_files.append(tmp_annotate_vcf_name)
                    tmp_annotate_vcf_name_err = tmp_annotate_vcf_name + ".err"
                    err_files.append(tmp_annotate_vcf_name_err)
                    tmp_files.append(tmp_annotate_vcf_name_err)

                    # Check if field in vcf header contains space, if yes, raise error (not valid)
                    log.debug("Check if field in vcf header contains space...")
                    for tmp_annotates_vcf_name in tmp_annotates_vcf_name_list:
                        log.debug(f"Check vcf header in file: {tmp_annotates_vcf_name}")
                        with bgzf.open(tmp_annotates_vcf_name, "rt") as f:
                            header_list = self.read_vcf_header(f)
                        annovar_vcf_header = vcf.Reader(io.StringIO("\n".join(header_list)))
                        for ann in annovar_vcf_header.infos:
                            if " " in ann:
                                msg_err = f"Annotations Annovar - field '{ann}' contains space (not valid)"
                                log.error(msg_err)
                                raise ValueError(msg_err)
                    log.debug("Check if field in vcf header contains space - done")

                    if len(tmp_annotates_vcf_name_list) > 1:

                        # List of annotated files
                        tmp_annotates_vcf_name_to_merge = " ".join(
                            tmp_annotates_vcf_name_list
                        )

                        # Command merge
                        merge_command = f"{bcftools_bin_command} merge --force-samples --threads={threads} {tmp_vcf_name} {tmp_annotates_vcf_name_to_merge} -o {tmp_annotate_vcf_name} -Oz  "
                        log.info(
                            f"Annotations Annovar - Annotation merging "
                            + str(len(tmp_annotates_vcf_name_list))
                            + " annotated files"
                        )
                        log.debug(f"Annotations Annovar - merge command: {merge_command}")
                        run_parallel_commands([merge_command], 1)

                    else:

                        # copy tmp_vcf_name to tmp_annotate_vcf_name
                        import shutil

                        shutil.copy(tmp_annotates_vcf_name_list[0], tmp_annotate_vcf_name)

                    # Update variants
                    log.info(f"Annotations Annovar - Updating...")
                    self.update_from_vcf(
                        tmp_annotate_vcf_name,
                        update_header=True,
                        annotation_header_fields_override=annotation_header_fields_override,
                        chrom_mapping_sql=chrom_mapping.from_tool_sql(),
                    )
                    remove_if_exists(
                        [
                            tmp_annotate_vcf_name,
                            tmp_annotate_vcf_name_err,
                            f"{tmp_annotate_vcf_name}.tbi",
                        ]
                    )

                # Clean files
                # Tmp file remove command
                remove_if_exists(tmp_files)

    # Parquet
    def annotation_parquet(self, section:str = "annotation", threads: int = None) -> None:
        """
        It takes a VCF file, and annotates it with a parquet file

        :param threads: number of threads to use for the annotation
        :type threads: int, optional

        :return: the value of the variable "result".
        """

        # DEBUG
        log.debug("Start annotation with parquet databases")

        # Threads
        if not threads:
            threads = self.get_threads()
        log.debug("Threads: " + str(threads))

        # Chunk size
        chunk_size = self.get_config().get("chunk_size", DEFAULT_CHUNK_SIZE)

        # DEBUG
        delete_tmp = True
        if self.get_config().get("verbosity", "warning") in ["debug"]:
            delete_tmp = False
            log.debug("Delete tmp files/folders: " + str(delete_tmp))

        # Config
        databases_folders = set(
            self.get_config()
            .get("folders", {})
            .get("databases", {})
            .get("annotations", [DEFAULT_ANNOTATIONS_FOLDER])
            + self.get_config()
            .get("folders", {})
            .get("databases", {})
            .get("parquet", [DEFAULT_PARQUET_FOLDER])
        )
        log.debug("Databases annotations: " + str(databases_folders))

        # Param
        annotations = (
            self.get_param()
            .get(section, {})
            .get("parquet", {})
            .get("annotations", None)
        )
        log.debug("Annotations: " + str(annotations))

        # Assembly
        assembly = self.get_param().get(
            "assembly", self.get_config().get("assembly", DEFAULT_ASSEMBLY)
        )

        # Force Update Annotation
        force_update_annotation = (
            self.get_param()
            .get(section, {})
            .get("options", {})
            .get("annotations_update", False)
        )
        log.debug(f"force_update_annotation={force_update_annotation}")
        force_append_annotation = (
            self.get_param()
            .get(section, {})
            .get("options", {})
            .get("annotations_append", False)
        )
        log.debug(f"force_append_annotation={force_append_annotation}")

        # Data
        table_variants = self.get_table_variants()

        # Fast strategy
        fast = self.get_param().get("fast", False)

        # Check if not empty
        log.debug("Check if not empty")
        sql_query_chromosomes_df = self.get_query_to_df(
            f"""SELECT count(*) as count FROM {table_variants} as table_variants LIMIT 1"""
        )
        if not sql_query_chromosomes_df["count"][0]:
            log.info("VCF empty")
            return

        # VCF header
        vcf_reader = self.get_header()
        log.debug("Initial header: " + str(vcf_reader.infos))

        # Nb Variants POS
        log.debug("NB Variants Start")
        nb_variants = self.conn.execute(
            "SELECT count(*) AS count FROM variants"
        ).fetchdf()["count"][0]
        log.debug("NB Variants Stop")

        # Existing annotations
        for vcf_annotation in self.get_header().infos:

            vcf_annotation_line = self.get_header().infos.get(vcf_annotation)
            log.debug(
                f"Existing annotations in VCF: {vcf_annotation} [{vcf_annotation_line}]"
            )

        # drop indexes
        log.debug("Drop indexes...")
        self.drop_indexes()

        # Init SQL query chromosomes dict
        sql_query_chromosomes_dict = None

        # Files or folders to remove
        files_or_folders_to_remove = []

        if annotations:

            if "ALL" in annotations:

                all_param = annotations.get("ALL", {})
                all_param_formats = all_param.get("formats", None)
                all_param_releases = all_param.get("releases", None)

                databases_infos_dict = self.scan_databases(
                    database_formats=all_param_formats,
                    database_releases=all_param_releases,
                )
                for database_infos in databases_infos_dict.keys():
                    if database_infos not in annotations:
                        annotations[database_infos] = {"INFO": None}

            # Update sources for all annotations
            update_sources = []

            for annotation in annotations:

                if annotation in ["ALL"]:
                    continue

                # Annotation Name
                annotation_name = os.path.basename(annotation)

                # Annotation fields
                # Add backward compatibility if "annotation_fields" is in annotation dict (instead of directly in annotation value)
                # Allow to add options for annotations
                if "annotation_fields" in annotations[annotation]:
                    annotation_fields = annotations[annotation].get(
                        "annotation_fields", {}
                    )
                else:
                    annotation_fields = annotations[annotation]

                # Options for annotations
                annotation_options = annotations[annotation].get("options", {})

                # Check options
                annotation_option_uniquify = annotation_options.get("uniquify", False)
                annotation_option_force_update_annotation = annotation_options.get(
                    "annotations_update", force_update_annotation
                )
                annotation_option_force_append_annotation = annotation_options.get(
                    "annotations_append", force_append_annotation
                )

                # Check annotation fields
                if not annotation_fields:
                    annotation_fields = {"INFO": None}

                log.info(f"Annotation '{annotation_name}'")
                log.debug(
                    f"Annotation '{annotation_name}' - fields: {annotation_fields}"
                )

                # Find file
                annotation = annotation_file_find(
                    annotation_file=annotation,
                    databases_folders=list(databases_folders),
                    assembly=assembly,
                )

                # Create Database
                database = Database(
                    database=annotation,
                    databases_folders=databases_folders,
                    assembly=assembly,
                )

                # Find files
                parquet_file = database.get_database()
                parquet_hdr_file = database.get_header_file()
                parquet_type = database.get_type()

                # Check if files exists
                if not parquet_file or not parquet_hdr_file:
                    msg_err_list = []
                    if not parquet_file:
                        msg_err_list.append(
                            "Annotation failed: Annotation file not found"
                        )
                    if parquet_file and not parquet_hdr_file:
                        msg_err_list.append(
                            f"Annotation failed: Annotation file '{parquet_file}' header not found. Check for file '{parquet_file}.hdr'"
                        )

                    log.error(". ".join(msg_err_list))
                    raise ValueError(". ".join(msg_err_list))
                else:
                    # Get parquet connexion
                    parquet_sql_attach = database.get_sql_database_attach(
                        output="query"
                    )
                    if parquet_sql_attach:
                        self.conn.execute(parquet_sql_attach)
                    parquet_file_link = database.get_sql_database_link()
                    # Log
                    log.debug(
                        f"Annotation '{annotation_name}' - file: "
                        + str(parquet_file)
                        + " and "
                        + str(parquet_hdr_file)
                    )

                    # Database full header columns
                    parquet_hdr_vcf_header_columns = database.get_header_file_columns(
                        parquet_hdr_file
                    )
                    # Log
                    log.debug(
                        "Annotation database header columns : "
                        + str(parquet_hdr_vcf_header_columns)
                    )

                    # Load header as VCF object
                    parquet_hdr_vcf_header_infos = database.get_header().infos

                    # Get extra infos
                    parquet_columns = database.get_extra_columns()
                    # Log
                    log.debug("Annotation database Columns: " + str(parquet_columns))

                    # Add extra columns if "ALL" in annotation_fields
                    # if "ALL" in annotation_fields:
                    #     allow_add_extra_column = True
                    if "ALL" in annotation_fields and database.get_extra_columns():
                        for extra_column in database.get_extra_columns():
                            if (
                                extra_column not in annotation_fields
                                and extra_column.replace("INFO/", "")
                                not in parquet_hdr_vcf_header_infos
                            ):
                                parquet_hdr_vcf_header_infos[extra_column] = (
                                    vcf.parser._Info(
                                        extra_column,
                                        ".",
                                        "String",
                                        f"{extra_column} description",
                                        "unknown",
                                        "unknown",
                                        self.code_type_map["String"],
                                    )
                                )

                    # For all fields in database
                    annotation_fields_all = False
                    if "ALL" in annotation_fields or "INFO" in annotation_fields:
                        annotation_fields_all = True
                        annotation_fields = {
                            key: key for key in parquet_hdr_vcf_header_infos
                        }

                        log.debug(
                            "Annotation database header - All annotations added: "
                            + str(annotation_fields)
                        )

                    # Init

                    # List of annotation fields to use
                    sql_query_annotation_update_info_sets = []

                    # List of annotation to agregate
                    sql_query_annotation_to_agregate = []

                    # Number of fields
                    nb_annotation_field = 0

                    # Annotation fields processed
                    annotation_fields_processed = []

                    # Columns mapping
                    map_columns = database.map_columns(
                        columns=annotation_fields, prefixes=["INFO/"]
                    )

                    # Query dict for fields to remove (update option)
                    query_dict_remove = {}

                    # Fetch Anotation fields
                    for annotation_field in annotation_fields:

                        # annotation_field_column
                        annotation_field_column = map_columns.get(
                            annotation_field, "INFO"
                        )

                        # field new name, if parametered
                        annotation_fields_new_name = annotation_fields.get(
                            annotation_field, annotation_field
                        )
                        if not annotation_fields_new_name:
                            annotation_fields_new_name = annotation_field

                        # To annotate
                        # force_update_annotation = True
                        # force_append_annotation = True
                        if annotation_field in parquet_hdr_vcf_header_infos and (
                            # force_update_annotation
                            annotation_option_force_update_annotation
                            # or force_append_annotation
                            or annotation_option_force_append_annotation
                            or (
                                annotation_fields_new_name
                                not in self.get_header().infos
                            )
                        ):

                            # Add field to annotation to process list
                            annotation_fields_processed.append(
                                annotation_fields_new_name
                            )

                            # explode infos for the field
                            annotation_fields_new_name_info_msg = ""
                            if (
                                # force_update_annotation
                                annotation_option_force_update_annotation
                                and annotation_fields_new_name
                                in self.get_header().infos
                            ):
                                # Remove field from INFO
                                query = f"""
                                    UPDATE {table_variants} as table_variants
                                    SET INFO = REGEXP_REPLACE(
                                                concat(table_variants.INFO,''),
                                                ';*{annotation_fields_new_name}=[^;]*',
                                                ''
                                                )
                                    WHERE concat(';',table_variants.INFO) LIKE '%;{annotation_fields_new_name}=%'
                                """
                                annotation_fields_new_name_info_msg = " [update]"
                                query_dict_remove[
                                    f"remove 'INFO/{annotation_fields_new_name}'"
                                ] = query

                            # Sep between fields in INFO
                            nb_annotation_field += 1
                            if nb_annotation_field > 1:
                                annotation_field_sep = ";"
                            else:
                                annotation_field_sep = ""

                            log.info(
                                f"Annotation '{annotation_name}' - '{annotation_field}' -> '{annotation_fields_new_name}'{annotation_fields_new_name_info_msg}"
                            )

                            # Add INFO field to header

                            # If regions, force values as list, due to overlap/aggregation
                            if parquet_type in ["regions"]:
                                parquet_hdr_vcf_header_infos_number = "."
                            else:
                                parquet_hdr_vcf_header_infos_number = (
                                    parquet_hdr_vcf_header_infos[annotation_field].num
                                    or "."
                                )
                            parquet_hdr_vcf_header_infos_type = (
                                parquet_hdr_vcf_header_infos[annotation_field].type
                                or "String"
                            )
                            parquet_hdr_vcf_header_infos_description = (
                                parquet_hdr_vcf_header_infos[annotation_field].desc
                                or f"{annotation_field} description"
                            )
                            parquet_hdr_vcf_header_infos_source = (
                                parquet_hdr_vcf_header_infos[annotation_field].source
                                or "unknown"
                            )
                            parquet_hdr_vcf_header_infos_version = (
                                parquet_hdr_vcf_header_infos[annotation_field].version
                                or "unknown"
                            )

                            vcf_reader.infos[annotation_fields_new_name] = (
                                vcf.parser._Info(
                                    annotation_fields_new_name,
                                    parquet_hdr_vcf_header_infos_number,
                                    parquet_hdr_vcf_header_infos_type,
                                    parquet_hdr_vcf_header_infos_description,
                                    parquet_hdr_vcf_header_infos_source,
                                    parquet_hdr_vcf_header_infos_version,
                                    self.code_type_map[
                                        parquet_hdr_vcf_header_infos_type
                                    ],
                                )
                            )

                            # Append
                            # if force_append_annotation:
                            if annotation_option_force_append_annotation:
                                query_case_when_append = f""" AND REGEXP_EXTRACT(concat(';', table_variants.INFO), ';{annotation_fields_new_name}=([^;]*)',1) IN ('','.') """
                            else:
                                query_case_when_append = ""

                            # Annotation/Update query fields
                            # Found in INFO column
                            if (
                                annotation_field_column == "INFO"
                                and "INFO" in parquet_hdr_vcf_header_columns
                            ):
                                sql_query_annotation_update_info_sets.append(
                                    f"""
                                CASE WHEN REGEXP_EXTRACT(concat(';', table_parquet.INFO), ';{annotation_field}=([^;]*)',1) NOT IN ('','.') {query_case_when_append}
                                        THEN concat('{annotation_field_sep}', '{annotation_fields_new_name}=', REGEXP_EXTRACT(concat(';', table_parquet.INFO), ';{annotation_field}=([^;]*)',1))
                                        ELSE ''
                                    END
                                """
                                )
                            # Found in a specific column
                            else:
                                sql_query_annotation_update_info_sets.append(
                                    f"""
                                CASE WHEN CAST(table_parquet."{annotation_field_column}" AS VARCHAR) NOT IN ('','.') {query_case_when_append}
                                        THEN concat('{annotation_field_sep}', '{annotation_fields_new_name}=', replace(CAST(table_parquet."{annotation_field_column}" AS VARCHAR), ';', ','))
                                        ELSE ''
                                    END
                                """
                                )
                                # if annotation_options.get("uniquify", False):
                                if annotation_option_uniquify:
                                    sql_query_annotation_to_agregate.append(
                                        f""" array_to_string(array_sort(array_distinct(string_split(string_agg(DISTINCT table_parquet_from."{annotation_field_column}", ','), ','))), ',') AS "{annotation_field_column}" """
                                    )
                                else:
                                    sql_query_annotation_to_agregate.append(
                                        f""" array_to_string(string_split(string_agg(table_parquet_from."{annotation_field_column}", ','), ','), ',') AS "{annotation_field_column}" """
                                    )

                        # Not to annotate
                        else:

                            # if force_update_annotation:
                            if annotation_option_force_update_annotation:
                                annotation_message = "forced"
                            else:
                                annotation_message = "skipped"

                            if annotation_field not in parquet_hdr_vcf_header_infos:
                                log.warning(
                                    f"Annotation '{annotation_name}' - '{annotation_field}' [{nb_annotation_field}] - not available in parquet file"
                                )
                            if annotation_fields_new_name in self.get_header().infos:
                                log.warning(
                                    f"Annotation '{annotation_name}' - '{annotation_fields_new_name}' [{nb_annotation_field}] - already exists in header ({annotation_message})"
                                )

                    # Check if ALL fields have to be annotated. Thus concat all INFO field
                    # allow_annotation_full_info = True
                    # allow_annotation_full_info = not force_append_annotation
                    allow_annotation_full_info = not annotation_option_force_append_annotation

                    if parquet_type in ["regions"]:
                        allow_annotation_full_info = False

                    if (
                        allow_annotation_full_info
                        and nb_annotation_field == len(annotation_fields)
                        and annotation_fields_all
                        and (
                            "INFO" in parquet_hdr_vcf_header_columns
                            and "INFO" in database.get_extra_columns()
                        )
                    ):
                        log.debug("Column INFO annotation enabled")
                        sql_query_annotation_update_info_sets = []
                        sql_query_annotation_update_info_sets.append(
                            " table_parquet.INFO "
                        )

                    if sql_query_annotation_update_info_sets:

                        # Annotate
                        log.info(
                            f"Annotation '{annotation_name}' - Create annotation table..."
                        )

                        # Join query annotation update info sets for SQL
                        sql_query_annotation_update_info_sets_sql = ",".join(
                            sql_query_annotation_update_info_sets
                        )

                        # Check chromosomes list (and variants infos)
                        if sql_query_chromosomes_dict is None:
                            sql_query_chromosomes_dict
                            sql_query_chromosomes = f"""
                                SELECT table_variants."#CHROM" as CHROM, count(*) AS count_variants, min(POS) AS min_variants, MAX(POS) AS max_variants
                                FROM {table_variants} as table_variants
                                GROUP BY table_variants."#CHROM"
                                ORDER BY table_variants."#CHROM"
                                """
                            sql_query_chromosomes_df = self.conn.execute(
                                sql_query_chromosomes
                            ).df()
                            sql_query_chromosomes_dict = {
                                entry["CHROM"]: {
                                    "count": entry["count_variants"],
                                    "min": entry["min_variants"],
                                    "max": entry["max_variants"],
                                }
                                for index, entry in sql_query_chromosomes_df.iterrows()
                            }

                        # Count total variants to annotate
                        total_variants_to_annotate = sum(
                            [
                                sql_query_chromosomes_dict[chromosome]["count"]
                                for chromosome in sql_query_chromosomes_dict
                            ]
                        )
                        log.debug(
                            f"Annotation '{annotation_name}' - Total variants to annotate: {total_variants_to_annotate}"
                        )

                        # Init
                        query_dict = query_dict_remove

                        # Spilling
                        duckdb_temp_directory = self.get_connexion_config().get(
                            "temp_directory", ".tmp"
                        )
                        duckdb_spilled = duckdb_has_spilled(duckdb_temp_directory)
                        log.debug(
                            f"Annotation '{annotation_name}' - DuckDB spilled: {duckdb_spilled}"
                        )

                        # Choose strategy
                        if duckdb_spilled:
                            annotation_table_strategy = "PARQUET"
                        else:
                            annotation_table_strategy = "TABLE"

                        # DEVEL
                        # TABLE, PARQUET, VIEW_BATCH, TABLE_BATCH, VIEW_SEQ, TABLE_SEQ
                        # annotation_table_strategy = "VIEW_SEQ"

                        # Log
                        log.debug(
                            f"Annotation '{annotation_name}' - Annotation table strategy: {annotation_table_strategy}"
                        )

                        # Create temporary table for batch update
                        sql_query_annotation_chrom_interval_pos_union = (
                            "annotation_chrom_interval_pos_union_" + get_random(10)
                        )

                        if annotation_table_strategy in ["TABLE"]:

                            # Create empty table for batch update
                            sql_create_empty_table = f"""
                                CREATE TABLE {sql_query_annotation_chrom_interval_pos_union} AS
                                    SELECT
                                        "#CHROM",
                                        "POS",
                                        "REF",
                                        "ALT",
                                        "INFO"
                                    FROM {table_variants} AS table_variants
                                    WHERE 1=0
                                ;
                            """
                            # log.debug(f"sql_create_empty_table={sql_create_empty_table}")
                            self.conn.execute(sql_create_empty_table)

                            # Add update source
                            source = {
                                "table": sql_query_annotation_chrom_interval_pos_union,
                                "join_keys": [
                                    "#CHROM",
                                    "POS",
                                    "REF",
                                    "ALT",
                                ],
                                "columns": {
                                    "INFO": {
                                        "columns": ["INFO"],
                                        "mode": "append",
                                        "separator": ";",
                                    }
                                },
                            }
                            update_sources.append(source)

                        elif annotation_table_strategy in ["PARQUET"]:

                            log.debug(
                                "DuckDB has spilled to disk during previous operations. "
                                "CTAS operation will use Parquet intermediate files to avoid out of memory errors. "
                                "This may slow down the operation. "
                                "Consider increasing DuckDB memory limit or adding more RAM to the system."
                            )

                            # Create empty parquet folder for batch update
                            # Folder will be removed after annotation
                            annotation_parquet_folder = os.path.join(
                                self.get_tmp_dir(),
                                # duckdb_temp_directory,
                                f"ctas_parquet_{get_random(10)}.partition.parquet",
                            )
                            os.makedirs(annotation_parquet_folder, exist_ok=True)
                            files_or_folders_to_remove.append(annotation_parquet_folder)

                            # Create temporary view for parquet batch update
                            # Query will be executed after the first batch is created (otherwise view fails)
                            sql_query_annotation_parquet_folder_view = f"""
                                CREATE VIEW {sql_query_annotation_chrom_interval_pos_union} AS SELECT * FROM read_parquet('{annotation_parquet_folder}/*/*.parquet', hive_partitioning = true)
                            """

                            # Add update source
                            source = {
                                "table": sql_query_annotation_chrom_interval_pos_union,
                                "join_keys": [
                                    "#CHROM",
                                    "POS",
                                    "REF",
                                    "ALT",
                                ],
                                "columns": {
                                    "INFO": {
                                        "columns": ["INFO"],
                                        "mode": "append",
                                        "separator": ";",
                                    }
                                },
                            }
                            update_sources.append(source)

                        # For each chromosome, first bacth by chromosome
                        for chrom in sql_query_chromosomes_dict:

                            # --- Number of variant by chromosome
                            nb_of_variant_by_chrom = sql_query_chromosomes_dict.get(
                                chrom, {}
                            ).get("count", 0)

                            # --- Prepare SQL clauses for annotation (regions or variants)

                            # Annotation with regions database
                            if parquet_type in ["regions"]:

                                sql_query_annotation_table = f"""
                                    WITH uniq_variants AS (
                                        SELECT DISTINCT
                                            "#CHROM",
                                            POS
                                        FROM variants
                                        WHERE "#CHROM" = '{chrom}'
                                    )
                                    SELECT
                                        '{chrom}' AS \"#CHROM\",
                                        table_variants_from.\"POS\" AS \"POS\",
                                        {",".join(sql_query_annotation_to_agregate)}
                                    FROM uniq_variants as table_variants_from
                                    LEFT JOIN {parquet_file_link} as table_parquet_from ON (
                                        table_parquet_from."#CHROM" = '{chrom}'
                                        AND table_variants_from.POS BETWEEN table_parquet_from.START AND table_parquet_from.END
                                    )
                                    WHERE table_variants_from.\"#CHROM\" in ('{chrom}')
                                    GROUP BY table_variants_from.\"POS\"
                                """

                                sql_query_annotation_from_clause = f"""
                                    (
                                        {sql_query_annotation_table}
                                    )
                                    as table_parquet
                                """

                                sql_query_annotation_join_using = """
                                    ("#CHROM", "POS")
                                """

                                sql_query_annotation_where_clause = f"""
                                    table_variants."#CHROM" = '{chrom}'
                                    AND table_variants."#CHROM" = table_parquet."#CHROM"
                                    AND table_variants.POS = table_parquet.POS
                                """

                            # Annotation with variants database
                            else:
                                sql_query_annotation_from_clause = f"""
                                    ({parquet_file_link}) as table_parquet
                                """
                                sql_query_annotation_where_clause = f"""
                                    table_variants."#CHROM" = '{chrom}'
                                    AND table_variants."#CHROM" = table_parquet."#CHROM"
                                    AND table_variants.POS = table_parquet.POS
                                    AND table_variants.REF = table_parquet.REF
                                    AND table_variants.ALT = table_parquet.ALT
                                """

                                sql_query_annotation_join_using = """
                                    ("#CHROM", "POS", "REF", "ALT")
                                """

                            # --- Batch by position intervals

                            # Get min/max POS and number of variants by chrom
                            nb_of_variant_by_chrom = sql_query_chromosomes_dict.get(
                                chrom, {}
                            ).get("count", 0)
                            min_of_variant_by_chrom = (
                                sql_query_chromosomes_dict.get(chrom, {}).get("min", 0)
                            ) - 1
                            max_of_variant_by_chrom = sql_query_chromosomes_dict.get(
                                chrom, {}
                            ).get("max", 0)

                            # Create batch queries by position intervals
                            batch_index = 0
                            nb_windows = (nb_of_variant_by_chrom // chunk_size) + 1
                            chunk_size_batch_update = (
                                int(
                                    (max_of_variant_by_chrom - min_of_variant_by_chrom)
                                    / nb_windows
                                )
                                + 1
                            )

                            # --- Create batch queries

                            # Create queries by position intervals
                            for start in range(
                                min_of_variant_by_chrom,
                                max_of_variant_by_chrom,
                                chunk_size_batch_update,
                            ):

                                # chunking positions
                                end = start + chunk_size_batch_update
                                batch_index += 1

                                # where clause in any cases, because it speedup the query
                                where_clause_batch_split = f" AND table_variants.POS >= {start} AND table_variants.POS < {end} "

                                # Init
                                sql_query_annotation_chrom_interval_pos = None

                                # Log
                                log.debug(
                                    f"Annotation '{annotation_name}' - Chromosome '{chrom}' [{nb_of_variant_by_chrom} variants] - batch [{batch_index}/{nb_windows}][{batch_index/nb_windows:.2%}%] - positions [{start}-{end}]..."
                                )

                                # Table strategy
                                # Standard TABLE
                                if annotation_table_strategy in ["TABLE"]:

                                    # Insert into annotation_chrom_interval_pos_union
                                    sql_query_annotation_chrom_interval_pos = f"""
                                        INSERT INTO {sql_query_annotation_chrom_interval_pos_union}
                                            SELECT
                                                table_variants."#CHROM",
                                                table_variants."POS",
                                                table_variants."REF",
                                                table_variants."ALT",
                                                concat(
                                                        {sql_query_annotation_update_info_sets_sql}
                                                        ) AS INFO
                                            FROM {table_variants} AS table_variants
                                            LEFT JOIN {sql_query_annotation_from_clause}
                                            USING {sql_query_annotation_join_using}
                                            WHERE {sql_query_annotation_where_clause}
                                                {where_clause_batch_split}
                                                ;
                                    """
                                    # log.debug(
                                    #     f"Annotation '{annotation_name}' - Chromosome '{chrom}' - batch [{batch_index}/{nb_windows}] - SQL query: {sql_query_annotation_chrom_interval_pos}"
                                    # )

                                # Parquet strategy
                                # TABLE with Parquet intermediate files
                                elif annotation_table_strategy in ["PARQUET"]:

                                    # Chromosome folder for parquet partitioning
                                    annotation_parquet_folder_chromosome = os.path.join(
                                        annotation_parquet_folder, f"#CHROM={chrom}"
                                    )

                                    # Create chromosome folder
                                    os.makedirs(
                                        annotation_parquet_folder_chromosome,
                                        exist_ok=True,
                                    )

                                    # Copy to parquet file
                                    sql_query_annotation_chrom_interval_pos = f"""
                                        COPY (
                                            SELECT
                                                table_variants."#CHROM",
                                                table_variants."POS",
                                                table_variants."REF",
                                                table_variants."ALT",
                                                concat(
                                                        {sql_query_annotation_update_info_sets_sql}
                                                        ) AS INFO
                                            FROM {table_variants} AS table_variants
                                            LEFT JOIN {sql_query_annotation_from_clause}
                                            USING {sql_query_annotation_join_using}
                                            WHERE {sql_query_annotation_where_clause}
                                                {where_clause_batch_split}
                                        ) TO '{annotation_parquet_folder_chromosome}/part_{batch_index}.parquet' (FORMAT 'parquet')
                                    """

                                # Alternative annotation strategies
                                # Not availble for now, use for devel and benchmarking
                                # Disabled for now because not efficient compared to batch parquet or table
                                elif annotation_table_strategy in [
                                    "VIEW_BATCH",
                                    "TABLE_BATCH",
                                    "VIEW_SEQ",
                                    "TABLE_SEQ",
                                ]:

                                    # Create type
                                    if annotation_table_strategy in [
                                        "VIEW_BATCH",
                                        "VIEW_SEQ",
                                    ]:
                                        create_type = "VIEW"
                                    else:
                                        create_type = "TABLE"

                                    # Create batch name
                                    sql_query_annotation_chrom_interval_pos_union_name = f"annotation_chrom_interval_pos_{chrom}_{batch_index}_{get_random(10)}"

                                    # Create
                                    sql_query_annotation_chrom_interval_pos = f"""
                                        CREATE OR REPLACE {create_type} {sql_query_annotation_chrom_interval_pos_union_name} AS
                                            SELECT
                                                table_variants."#CHROM",
                                                table_variants."POS",
                                                table_variants."REF",
                                                table_variants."ALT",
                                                concat(
                                                        {sql_query_annotation_update_info_sets_sql}
                                                        ) AS INFO
                                            FROM {table_variants} AS table_variants
                                            LEFT JOIN {sql_query_annotation_from_clause}
                                            USING {sql_query_annotation_join_using}
                                            WHERE {sql_query_annotation_where_clause}
                                                {where_clause_batch_split}
                                    """

                                    # Add update source
                                    source = {
                                        "table": sql_query_annotation_chrom_interval_pos_union_name,
                                        "join_keys": [
                                            "#CHROM",
                                            "POS",
                                            "REF",
                                            "ALT",
                                        ],
                                        "columns": {
                                            "INFO": {
                                                "columns": ["INFO"],
                                                "mode": "append",
                                                "separator": ";",
                                            }
                                        },
                                    }

                                    if annotation_table_strategy in [
                                        "VIEW_SEQ",
                                        "TABLE_SEQ",
                                    ]:
                                        log.debug(
                                            f"Annotations - Perform with {source} ..."
                                        )
                                        self.update_table(
                                            dest_table=table_variants,
                                            sources=[source],
                                            samples=10000,
                                            force_strategy="update",
                                            chromosomes=None,
                                            only_strategy=False,
                                        )
                                    else:
                                        update_sources.append(source)

                                # Execute query for the batch, if exists
                                if sql_query_annotation_chrom_interval_pos:
                                    self.get_connexion().execute(
                                        sql_query_annotation_chrom_interval_pos
                                    )

                                # Create view for parquet only once
                                if (
                                    annotation_table_strategy in ["PARQUET"]
                                    and sql_query_annotation_parquet_folder_view
                                    is not None
                                ):
                                    # Execute view creation
                                    self.get_connexion().execute(
                                        sql_query_annotation_parquet_folder_view
                                    )
                                    # Reset to None to avoid re-creation
                                    sql_query_annotation_parquet_folder_view = None

                                    # Check content of view

                        # Execute queries one by one to have logging
                        # For queries in query_dict, mainly to remove annotations to update

                        # SET max_expression_depth TO x
                        self.conn.execute("SET max_expression_depth TO 10000")

                        # Execute queries one by one
                        # Queries mainly for update options, to prepare table by removing existing annotations
                        num_query = 0
                        for query_name in query_dict:
                            query = query_dict[query_name]
                            num_query += 1
                            log.info(
                                f"Annotation '{annotation_name}' - Annotation - Query {query_name}..."
                            )
                            self.get_connexion().execute(query)

                    else:

                        # log
                        log.info(
                            f"Annotation '{annotation_name}' - No Annotations available"
                        )

            # Perform update
            if len(update_sources) > 0:
                log.info(f"Annotations - Perform with {len(update_sources)} sources...")

                # Fast strategy
                if fast:
                    strategy = "fast"
                else:
                    strategy = None

                self.update_table(
                    dest_table=table_variants,
                    sources=update_sources,
                    samples=10000,
                    force_strategy=strategy,
                    chromosomes=None,
                    only_strategy=False,
                )
            else:
                log.info("No annotation sources")

        # extend temporary files list
        self.add_tmp_files(files_or_folders_to_remove)

        # Remove temporary files or folders
        # remove_if_exists(files_or_folders_to_remove)

    def annotation_splice(self, section:str = "annotation", threads: int = None) -> None:
        """
        This function annotate with snpEff

        :param threads: The number of threads to use
        :return: the value of the variable "return_value".
        """

        from howard.objects.variants import Variants

        # DEBUG
        log.debug("Start annotation with splice tools")

        # Threads
        if not threads:
            threads = self.get_threads()
        log.debug("Threads: " + str(threads))

        # DEBUG
        delete_tmp = True
        if self.get_config().get("verbosity", "warning") in ["debug"]:
            delete_tmp = False
            log.debug("Delete tmp files/folders: " + str(delete_tmp))

        # Config
        config = self.get_config()
        log.debug("Config: " + str(config))
        splice_config = config.get("tools", {}).get("splice", {})
        if not splice_config:
            splice_config = DEFAULT_TOOLS_BIN.get("splice", {})
            msg_err = "No Splice tool config"
            raise ValueError(msg_err)
        log.debug(f"splice_config: {splice_config}")

        # Config - Folders - Databases
        databases_folders = (
            config.get("folders", {}).get("databases", {}).get("splice", ["."])
        )
        log.debug("Databases annotations: " + str(databases_folders))

        # Splice docker image
        splice_docker_image = splice_config.get("docker").get("image")

        # Pull splice image if it's not already there
        if not check_docker_image_exists(splice_docker_image):
            log.warning(
                f"Annotation: splice docker image {splice_docker_image} not found locally, trying to pull from dockerhub"
            )
            try:
                command(f"docker pull {splice_config.get('docker').get('image')}")
            except subprocess.CalledProcessError:
                msg_err = f"Unable to find docker {splice_docker_image} on dockerhub"
                log.error(msg_err)
                raise ValueError(msg_err)

        # Config - splice databases
        splice_databases = (
            config.get("folders", {})
            .get("databases", {})
            .get("splice", DEFAULT_SPLICE_FOLDER)
        )
        splice_databases = full_path(splice_databases)

        # Param
        param = self.get_param()
        log.debug("Param: " + str(param))

        # Param
        options = param.get(section, {}).get("splice", {}).get("options", {})
        log.debug("Options: " + str(options))

        # param - chrom mapping
        chrom_mapping_options = param.get(section, {}).get("splice", {}).get("chrom_mapping", None)
        chrom_mapping = ChromMapping(
            mapping=chrom_mapping_options
        )

        # Data
        table_variants = self.get_table_variants()

        # Check if not empty
        log.debug("Check if not empty")
        sql_query_chromosomes = (
            f"""SELECT count(*) as count FROM {table_variants} as table_variants"""
        )
        if not self.get_query_to_df(f"{sql_query_chromosomes}")["count"][0]:
            log.info("VCF empty")
            return None

        # Export in VCF
        log.debug("Create initial file to annotate")

        # Create output folder / work folder
        if options.get("output_folder", ""):
            output_folder = options.get("output_folder", "")
            if not os.path.exists(output_folder):
                Path(output_folder).mkdir(parents=True, exist_ok=True)
        else:
            output_folder = os.path.join(self.get_tmp_dir(), f"splice-{get_random()}")
            if not os.path.exists(output_folder):
                Path(output_folder).mkdir(parents=True, exist_ok=True)

        if options.get("workdir", ""):
            workdir = options.get("workdir", "")
        else:
            workdir = "/work"

        # Create tmp VCF file
        tmp_vcf = NamedTemporaryFile(
            prefix=self.get_prefix(),
            dir=output_folder,
            suffix=".vcf",
            delete=False,
        )
        tmp_vcf_name = tmp_vcf.name

        # VCF header
        header = self.get_header()

        # Existing annotations
        for vcf_annotation in self.get_header().infos:

            vcf_annotation_line = self.get_header().infos.get(vcf_annotation)
            log.debug(
                f"Existing annotations in VCF: {vcf_annotation} [{vcf_annotation_line}]"
            )

        # Memory limit
        if config.get("memory", None):
            memory_limit = config.get("memory", "8G").upper()
            # upper()
        else:
            memory_limit = "8G"
        log.debug(f"memory_limit: {memory_limit}")

        # Check number of variants to annotate
        where_clause_regex_spliceai = r"SpliceAI_\w+"
        where_clause_regex_spip = r"SPiP_\w+"
        where_clause = f""" WHERE NOT regexp_matches("INFO", '{where_clause_regex_spliceai}') AND NOT regexp_matches("INFO", '{where_clause_regex_spip}')"""
        df_list_of_variants_to_annotate = self.get_query_to_df(
            query=f""" SELECT * FROM variants {where_clause} """
        )
        if len(df_list_of_variants_to_annotate) == 0:
            log.warning(
                f"No variants to annotate with splice. Variants probably already annotated with splice"
            )
            return None
        else:
            log.info(f"Annotation: {len(df_list_of_variants_to_annotate)} variants")

        # Export VCF file
        self.export_variant_vcf(
            vcf_file=tmp_vcf_name,
            remove_info=True,
            add_samples=True,
            index=False,
            where_clause=where_clause,
            chrom_mapping_sql=chrom_mapping.to_tool_sql()
        )
        mount = [f" -v {path}:{path}:rw" for path in [output_folder]]
        if any(value for value in splice_config.values() if value is None):
            log.warning("At least one splice config parameter is empty")
            # exit annotation_splice
            return None

        # Params in splice nf
        def check_values(dico: dict):
            """
            Ensure parameters for NF splice pipeline
            """
            for key, val in dico.items():
                if key == "genome":
                    if any(
                        assemb in options.get("genome", {})
                        for assemb in ["hg19", "GRCh37", "grch37", "GRCH37"]
                    ):
                        yield f"--{key} hg19"
                    elif any(
                        assemb in options.get("genome", {})
                        for assemb in ["hg38", "GRCh38", "grch38", "GRCH38"]
                    ):
                        yield f"--{key} hg38"
                elif (
                    (isinstance(val, str) and val)
                    or isinstance(val, int)
                    or isinstance(val, bool)
                ):
                    yield f"--{key} {val}"

        # Genome
        genome = options.get("genome", config.get("assembly", DEFAULT_ASSEMBLY))
        options["genome"] = genome
        # NF params
        nf_params = []
        # Add options
        if options:
            log.debug(options)
            nf_params = list(check_values(options))
            log.debug(f"Splice NF params: {' '.join(nf_params)}")
        else:
            log.debug("No NF params provided")
        # Add threads
        if "threads" not in options.keys():
            nf_params.append(f"--threads {threads}")
        # Genome path
        genome_path = find_genome(
            config.get("folders", {})
            .get("databases", {})
            .get("genomes", DEFAULT_GENOME_FOLDER),
            file=f"{genome}.fa",
        )
        # Add genome path
        if not genome_path:
            raise ValueError(
                f"Can't find genome assembly {genome}.fa in {config.get('folders', {}).get('databases', {}).get('genomes', DEFAULT_GENOME_FOLDER)}"
            )
        else:
            log.debug(f"Genome: {genome_path}")
            nf_params.append(f"--genome_path {genome_path}")

        def splice_annotations(options: dict = {}, config: dict = {}) -> list:
            """
            Setting up updated databases for SPiP and SpliceAI
            """

            try:

                # SpliceAI assembly transcriptome
                spliceai_assembly = os.path.join(
                    config.get("folders", {}).get("databases", {}).get("spliceai", {}),
                    options.get("genome"),
                    "transcriptome",
                )
                spip_assembly = options.get("genome")

                spip = find(
                    f"transcriptome_{spip_assembly}.RData",
                    config.get("folders", {}).get("databases", {}).get("spip", {}),
                )
                spliceai = find("spliceai.refseq.txt", spliceai_assembly)
                log.debug(f"SPiP annotations: {spip}")
                log.debug(f"SpliceAI annotations: {spliceai}")
                if spip and spliceai:
                    return [
                        f"--spip_transcriptome {spip}",
                        f"--spliceai_transcriptome {spliceai}",
                    ]
                else:
                    log.warning(
                        "Can't find splice databases in configuration, use annotations file from image"
                    )
            except TypeError:
                log.warning(
                    "Can't find splice databases in configuration, use annotations file from image"
                )
                return []

        # Add options, check if transcriptome option have already beend provided
        if (
            "spip_transcriptome" not in nf_params
            and "spliceai_transcriptome" not in nf_params
        ):
            splice_reference = splice_annotations(options, config)
            if splice_reference:
                nf_params.extend(splice_reference)
        # nf_params.append(f"--output_folder {output_folder}")
        random_uuid = f"HOWARD-SPLICE-{get_random()}"
        cmd = f"nextflow -log {os.path.join(output_folder, f'{random_uuid}.log')} -c /app/SpliceToolBox/src/splicetoolbox/nextflow/nextflow.docker.config run /app/SpliceToolBox/src/splicetoolbox/nextflow/main.nf -entry SPLICE --vcf {tmp_vcf_name} {' '.join(nf_params)} -profile standard,conda,singularity,report,timeline"
        log.debug(cmd)
        splice_config["docker"]["command"] = cmd

        # Ensure proxy is set
        proxy = [
            f"-e {var}={os.getenv(var)}"
            for var in ["https_proxy", "http_proxy", "ftp_proxy"]
            if os.getenv(var) is not None
        ]
        docker_cmd = get_bin_command(
            tool="splice",
            bin_type="docker",
            config=config,
            default_folder=f"{DEFAULT_TOOLS_FOLDER}/docker",
            add_options=f"--name {random_uuid} {' '.join(mount)} -e NXF_DISABLE_CHECK_LATEST=true {' '.join(proxy)}",
        )
        # print(docker_cmd)
        # exit()
        # Docker debug
        # if splice_config.get("rm_container"):
        #     rm_container = "--rm"
        # else:
        #     rm_container = ""
        # docker_cmd = f"docker run {rm_container} --entrypoint '/bin/bash' --name {random_uuid} {' '.join(mount)} {':'.join(splice_config.get('image'))} {cmd}"

        log.info("Starting Splice annotation...")

        log.debug(docker_cmd)
        res = subprocess.run(docker_cmd, shell=True, capture_output=True, text=True)
        log.debug(res.stdout)
        if res.stderr:
            log.warning(f"Splice annotation warning/error:")
            log.warning(res.stderr)
            #log.error(res.stderr)
        res.check_returncode()
        # Update variants
        log.info("Annotation - Updating...")
        # Test find output vcf
        log.debug(
            f"TMP splice output: {os.path.basename(tmp_vcf_name).replace('.vcf', '')}.spip.spliceai.sorted.vcf.gz"
        )
        output_vcf = []
        # Wrong folder to look in
        for files in os.listdir(os.path.dirname(tmp_vcf_name)):
            if (
                files
                == f"{os.path.basename(tmp_vcf_name).replace('.vcf', '')}.spip.spliceai.sorted.vcf.gz"
            ):
                output_vcf.append(os.path.join(os.path.dirname(tmp_vcf_name), files))
        # log.debug(os.listdir(options.get("output_folder")))
        log.debug(f"Splice annotated vcf: {output_vcf[0]}")
        if not output_vcf:
            log.debug(
                f"Splice output was not generated {os.path.basename(tmp_vcf_name)}*.spip.spliceai.sorted.vcf.gz"
            )
        else:

            self.update_from_vcf(
                output_vcf[0],
                update_header=True,
                chrom_mapping_sql=chrom_mapping.from_tool_sql()
            )

        # Remove file
        remove_if_exists(output_vcf)

    ###
    # HGVS
    ###

    def annotation_hgvs(self, section:str = "annotation", threads: int = None) -> None:
        """
        The `annotation_hgvs` function performs HGVS annotation on a set of variants using genomic
        coordinates and alleles.

        :param threads: The `threads` parameter is an optional integer that specifies the number of
        threads to use for parallel processing. If no value is provided, it will default to the number
        of threads obtained from the `get_threads()` method
        :type threads: int
        """

        import dask.dataframe as dd  # type: ignore

        # Function for each partition of the Dask Dataframe
        def partition_function(partition):
            """
            The function `partition_function` applies the `annotation_hgvs_partition` function to
            each row of a DataFrame called `partition`.

            :param partition: The parameter "partition" is a pandas DataFrame that contains the data
            to be processed
            :return: the result of applying the "annotation_hgvs_partition" function to each row of
            the "partition" dataframe along the axis 1.
            """
            return partition.apply(annotation_hgvs_partition, axis=1)

        def annotation_hgvs_partition(row) -> str:
            """
            The function `annotation_hgvs_partition` takes in a row of data and returns a string
            containing a list of HGVS names associated with the given genomic coordinates and alleles.

            :param row: A dictionary-like object that contains the values for the following keys:
            :return: a string that contains the HGVS names associated with the given row of data.
            """

            chr = row["CHROM"]
            pos = row["POS"]
            ref = row["REF"]
            alt = row["ALT"]

            # Find list of associated transcripts
            transcripts_list = list(polars_conn.execute(f"""
                SELECT transcript
                FROM refseq_df
                WHERE CHROM='{chr}'
                AND POS={pos}
            """)["transcript"])

            # Full HGVS annotation in list
            hgvs_full_list = []

            for transcript_name in transcripts_list:

                # Transcript
                transcript = get_transcript(
                    transcripts=transcripts, transcript_name=transcript_name
                )
                # Exon
                if use_exon:
                    exon = transcript.find_exon_number(pos)
                else:
                    exon = None
                # Protein
                transcript_protein = None
                if use_protein or add_protein or full_format:
                    transcripts_protein = list(polars_conn.execute(f"""
                        SELECT protein
                        FROM refseqlink_df
                        WHERE transcript='{transcript_name}'
                        LIMIT 1
                    """)["protein"])
                    if len(transcripts_protein):
                        transcript_protein = transcripts_protein[0]

                # HGVS name
                hgvs_name = format_hgvs_name(
                    chr,
                    pos,
                    ref,
                    alt,
                    genome=genome,
                    transcript=transcript,
                    transcript_protein=transcript_protein,
                    exon=exon,
                    use_gene=use_gene,
                    use_protein=use_protein,
                    full_format=full_format,
                    use_version=use_version,
                    codon_type=codon_type,
                )
                hgvs_full_list.append(hgvs_name)
                if add_protein and not use_protein and not full_format:
                    hgvs_name = format_hgvs_name(
                        chr,
                        pos,
                        ref,
                        alt,
                        genome=genome,
                        transcript=transcript,
                        transcript_protein=transcript_protein,
                        exon=exon,
                        use_gene=use_gene,
                        use_protein=True,
                        full_format=False,
                        use_version=use_version,
                        codon_type=codon_type,
                    )
                    hgvs_full_list.append(hgvs_name)

            # Create liste of HGVS annotations
            hgvs_full = ",".join(hgvs_full_list)

            return hgvs_full

        # Polars connexion
        polars_conn = pl.SQLContext(register_globals=True, eager=True)

        # Config
        config = self.get_config()

        # Databases
        # Genome
        databases_genomes_folders = (
            config.get("folders", {})
            .get("databases", {})
            .get("genomes", DEFAULT_GENOME_FOLDER)
        )
        databases_genome = (
            config.get("folders", {}).get("databases", {}).get("genomes", "")
        )
        # refseq database folder
        databases_refseq_folders = (
            config.get("folders", {})
            .get("databases", {})
            .get("refseq", DEFAULT_REFSEQ_FOLDER)
        )
        # refseq
        databases_refseq = config.get("databases", {}).get("refSeq", None)
        # refSeqLink
        databases_refseqlink = config.get("databases", {}).get("refSeqLink", None)

        # Param
        param = self.get_param()

        # Quick HGVS
        if "hgvs_options" in param and param.get("hgvs_options", ""):
            log.info(f"Quick HGVS Annotation:")
            if not param.get(section, {}).get("hgvs", ""):
                if section not in param:
                    param[section] = {}
                param[section]["hgvs"] = {}
            for option in param.get("hgvs_options", "").split(","):
                option_var_val = option.split("=")
                option_var = option_var_val[0]
                if len(option_var_val) > 1:
                    option_val = option_var_val[1]
                else:
                    option_val = "True"
                if option_val.upper() in ["TRUE"]:
                    option_val = True
                elif option_val.upper() in ["FALSE"]:
                    option_val = False
                log.info(f"   {option_var}={option_val}")
                param[section][option_var] = option_val

        # Check if HGVS annotation enabled
        log.debug(f"HGVS Annotation param: {param}")
        log.debug(f"HGVS Annotation section: {section}")
        if "hgvs" in param.get(section, {}):
            log.info("HGVS Annotation... ")
            for hgvs_option in param.get(section, {}).get("hgvs", {}):
                log.info(f""" {hgvs_option}: {param.get(section, {}).get("hgvs", {}).get(hgvs_option)} """)
        else:
            return

        # HGVS Param
        param_hgvs = param.get(section, {}).get("hgvs", {})
        use_exon = param_hgvs.get("use_exon", False)
        use_gene = param_hgvs.get("use_gene", False)
        use_protein = param_hgvs.get("use_protein", False)
        add_protein = param_hgvs.get("add_protein", False)
        full_format = param_hgvs.get("full_format", False)
        use_version = param_hgvs.get("use_version", False)
        codon_type = param_hgvs.get("codon_type", "3")

        # refSseq refSeqLink
        databases_refseq = param_hgvs.get("refseq", databases_refseq)
        databases_refseqlink = param_hgvs.get("refseqlink", databases_refseqlink)

        # Assembly
        assembly = self.get_assembly()

        # Genome
        genome_file = None
        if find_genome(databases_genome):
            genome_file = find_genome(databases_genome)
        else:
            genome_file = find_genome(
                genome_path=databases_genomes_folders, assembly=assembly
            )
        log.debug("Genome: " + str(genome_file))

        # refSseq
        refseq_file = find_file_prefix(
            input_file=databases_refseq,
            prefix="ncbiRefSeq",
            folder=databases_refseq_folders,
            assembly=assembly,
        )
        log.debug("refSeq: " + str(refseq_file))

        # refSeqLink
        refseqlink_file = find_file_prefix(
            input_file=databases_refseqlink,
            prefix="ncbiRefSeqLink",
            folder=databases_refseq_folders,
            assembly=assembly,
        )
        log.debug("refSeqLink: " + str(refseqlink_file))

        # Threads
        if not threads:
            threads = self.get_threads()
        log.debug("Threads: " + str(threads))

        # Variables
        table_variants = self.get_table_variants(clause="update")

        # Get variants SNV and InDel only
        query_variants = f"""
            SELECT "#CHROM" AS CHROM, POS, REF, ALT
            FROM {table_variants}
            WHERE REF ~ '^[A-Za-z]+$' AND ALT ~ '^[A-Za-z]+$'
            """
        df_variants = self.get_query_to_df(query_variants)

        if len(df_variants) == 0:
            log.debug("No variants found for HGVS annotation")
            return

        # Added columns
        added_columns = []

        # Add hgvs column in variants table
        hgvs_column_name = "hgvs_" + str(get_random())
        added_column = self.add_column(
            table_variants, hgvs_column_name, "STRING", default_value=None
        )
        added_columns.append(added_column)

        log.debug(f"refSeq loading...")
        # refSeq in duckDB
        refseq_table = get_refseq_table(
            conn=self.conn, refseq_table="refseq", refseq_file=refseq_file
        )
        # Loading all refSeq in Dataframe
        refseq_query = f"""
            SELECT df_variants.CHROM, df_variants.POS, {refseq_table}.name AS transcript
            FROM {refseq_table}
            JOIN df_variants ON (
                {refseq_table}.chrom = df_variants.CHROM
                AND {refseq_table}.txStart<=df_variants.POS
                AND {refseq_table}.txEnd>=df_variants.POS
            )
        """
        refseq_df = self.conn.query(refseq_query).pl()

        if refseqlink_file:
            log.debug(f"refSeqLink loading...")
            # refSeqLink in duckDB
            refseqlink_table = get_refseq_table(
                conn=self.conn, refseq_table="refseqlink", refseq_file=refseqlink_file
            )
            # Loading all refSeqLink in Dataframe
            protacc_column = "protAcc_with_ver"
            mrnaacc_column = "mrnaAcc_with_ver"
            refseqlink_query = f"""
                SELECT {refseq_table}.chrom, {protacc_column} AS protein, {mrnaacc_column} AS transcript
                FROM {refseqlink_table} 
                JOIN {refseq_table} ON ({refseq_table}.name = {refseqlink_table}.mrnaAcc_with_ver)
                WHERE protAcc_without_ver IS NOT NULL
            """
            # Polars Dataframe
            refseqlink_df = self.conn.query(f"{refseqlink_query}").pl()

        # Read RefSeq transcripts into a python dict/model.
        log.debug(f"Transcripts loading...")
        with tempfile.TemporaryDirectory() as tmpdir:
            transcripts_query = f"""
                COPY (
                    SELECT {refseq_table}.*
                    FROM {refseq_table}
                    JOIN df_variants ON (
                        {refseq_table}.chrom=df_variants.CHROM
                        AND {refseq_table}.txStart<=df_variants.POS
                        AND {refseq_table}.txEnd>=df_variants.POS
                    )
                )
                TO '{tmpdir}/transcript.tsv' (DELIMITER '\t');
            """
            self.conn.query(transcripts_query)
            with open(f"{tmpdir}/transcript.tsv") as infile:
                transcripts = read_transcripts(infile)

        # Polars connexion
        polars_conn = pl.SQLContext(register_globals=True, eager=True)

        log.debug("Genome loading...")
        # Read genome sequence using pyfaidx.
        genome = Fasta(genome_file)

        log.debug("Start annotation HGVS...")

        # Create
        # a Dask Dataframe from Pandas dataframe with partition as number of threads
        ddf = dd.from_pandas(df_variants, npartitions=threads)

        # Use dask.dataframe.apply() to apply function on each partition
        ddf[hgvs_column_name] = ddf.map_partitions(partition_function)

        # Convert Dask DataFrame to Pandas Dataframe
        df = ddf.compute()

        # Convert Pandas dataframe to parquet (due to error in cast VARCHAR -> NULL ???)
        with tempfile.TemporaryDirectory() as tmpdir:
            df_parquet = os.path.join(tmpdir, "df.parquet")
            df.to_parquet(df_parquet)

            # Update hgvs column
            update_variant_query = f"""
                UPDATE {table_variants}
                SET "{hgvs_column_name}"=df."{hgvs_column_name}"
                FROM read_parquet('{df_parquet}') as df
                WHERE variants."#CHROM" = df.CHROM
                AND variants.POS = df.POS
                AND variants.REF = df.REF
                AND variants.ALT = df.ALT
                AND df."{hgvs_column_name}" NOT IN ('') AND df."{hgvs_column_name}" NOT NULL
                """
            self.execute_query(update_variant_query)

        # Update INFO column
        sql_query_update = f"""
            UPDATE {table_variants}
            SET INFO = 
                concat(
                    CASE 
                        WHEN INFO NOT IN ('','.')
                        THEN concat(INFO, ';')
                        ELSE ''
                    END,
                    'hgvs=',
                    {hgvs_column_name}
                )
            WHERE "{hgvs_column_name}" NOT IN ('') AND "{hgvs_column_name}" NOT NULL
            """
        # log.debug(f"Update INFO column with HGVS: {sql_query_update}")
        self.execute_query(sql_query_update)

        # Add header
        HGVS_INFOS = {
            "hgvs": {
                "ID": "hgvs",
                "Number": ".",
                "Type": "String",
                "Description": f"HGVS annotatation with HOWARD",
            }
        }

        for field in HGVS_INFOS:
            field_ID = HGVS_INFOS[field]["ID"]
            field_description = HGVS_INFOS[field]["Description"]
            self.get_header().infos[field_ID] = vcf.parser._Info(
                field_ID,
                HGVS_INFOS[field]["Number"],
                HGVS_INFOS[field]["Type"],
                field_description,
                "unknown",
                "unknown",
                code_type_map[HGVS_INFOS[field]["Type"]],
            )

        # Remove added columns
        for added_column in added_columns:
            self.drop_column(column=added_column)
