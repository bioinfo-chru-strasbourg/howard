import gc
import os
import random
import re
import json
import yaml  # type: ignore
import pandas as pd  # type: ignore
import vcf  # type: ignore
import logging as log

from howard.functions.commons import (
    add_value_into_dict,
    barcode,
    clean_annotation_field,
    findbypipeline,
    full_path,
    genotype_stats,
    genotypeconcordance,
    transcripts_file_to_df,
    trio,
    vaf_normalization,
    code_type_map_to_vcf,
    get_random,
)


class variants_calculation:

    ###############
    # Calculation #
    ###############

    def get_config_calculations_default(self) -> dict:
        """
        The function `get_config_calculations_default` returns a dictionary containing default configurations for
        various calculations.

        :return: The function `get_config_calculations_default` returns a dictionary containing default configuration
        settings for different calculations.
        """

        config_default = {
            "variant_chr_pos_alt_ref": {
                "type": "sql",
                "name": "variant_chr_pos_alt_ref",
                "description": "Create a variant ID with chromosome, position, alt and ref",
                "comment": "This calculation generates a variant ID with chromosome, position, alt and ref, with format '#CHROM_POS_REF_ALT'. Useful for variant identification and comparison across datasets",
                "available": True,
                "output_column_name": "variant_chr_pos_alt_ref",
                "output_column_number": 1,
                "output_column_type": "String",
                "output_column_description": "variant ID with chromosome, position, alt and ref",
                "operation_query": """ concat("#CHROM", '_', "POS", '_', "REF", '_', "ALT") """,
                "operation_info": True,
            },
            "VARTYPE": {
                "type": "sql",
                "name": "VARTYPE",
                "description": "Variant type (e.g. SNV, INDEL, MNV, BND...)",
                "comment": "Determine the type of variant based on its characteristics, such as the length of the reference and alternate alleles, and the presence of structural variant information. This calculation classifies variants into types like SNV, INDEL, MNV, BND, etc., which is essential for downstream analysis and interpretation.",
                "available": True,
                "table": "variants",
                "output_column_name": "VARTYPE",
                "output_column_number": 1,
                "output_column_type": "String",
                "output_column_description": "Variant type: SNV if X>Y, MOSAIC if X>Y,Z or X,Y>Z, INDEL if XY>Z or X>YZ",
                "operation_query": """
                        CASE
                            WHEN "SVTYPE" NOT NULL THEN "SVTYPE"
                            WHEN LENGTH(REF) = 1 AND LENGTH(ALT) = 1 THEN 'SNV'
                            WHEN REF LIKE '%,%' OR ALT LIKE '%,%' THEN 'MOSAIC'
                            WHEN LENGTH(REF) == LENGTH(ALT) AND LENGTH(REF) > 1 THEN 'MNV'
                            WHEN LENGTH(REF) <> LENGTH(ALT) THEN 'INDEL'
                            ELSE 'UNDEFINED'
                        END
                        """,
                "info_fields": ["SVTYPE"],
                "operation_info": True,
            },
            "MERGED_HGVS": {
                "type": "sql",
                "name": "merged_hgvs",
                "description": "Merge HGVS nomenclatures from snpEff (snpeff_hgvs) and ANNOVAR (AAChange_refGene) into merged_hgvs field",
                "comment": "This calculation merges HGVS nomenclatures from snpEff (snpeff_hgvs) and ANNOVAR (AAChange_refGene) into a single field called merged_hgvs. It aggregates distinct HGVS annotations from both sources, concatenating them with a comma separator. This unified field allows for easier access to HGVS information regardless of the annotation source, facilitating downstream analysis and interpretation.",
                "available": True,
                "table": "variants",
                "output_column_name": "merged_hgvs",
                "output_column_number": ".",
                "output_column_type": "String",
                "output_column_description": "Merge HGVS nomenclatures from snpEff (snpeff_hgvs) and ANNOVAR (AAChange_refGene) into merged_hgvs field",
                "operation_query": """
                        list_aggregate(
                            list_distinct(
                                list_reduce(
                                    [
                                        "snpeff_hgvs",
                                        "AAChange_refGene",
                                    ],
                                    (a, b) -> list_concat(a, b)
                                )
                            ),
                            'string_agg',
                            ','
                        )
                        """,
                "info_fields": ["snpeff_hgvs", "AAChange_refGene"],
                "operation_info": True,
            },
            "snpeff_extract": {
                "type": "python",
                "name": "snpeff_extract",
                "description": "HGVS nomenclatures from snpEff annotation",
                "comment": "Extract HGVS nomenclatures from snpEff annotation (field ANN) and create new INFO fields with prefix 'snpeff_' (e.g. snpeff_hgvs, snpeff_impact, snpeff_gene_name...). This calculation parses the ANN field from snpEff annotations, extracting relevant information such as HGVS nomenclatures, impact, gene name, etc., and creates new INFO fields with a 'snpeff_' prefix for easier access and downstream analysis.",
                "available": True,
                "function_name": "calculation_snpeff_extract",
                "function_params": {
                    "snpeff_field": "ANN",
                    "snpeff_hgvs": "snpeff_hgvs",
                    "snpeff_explode": "snpeff_",
                    "snpeff_json": "snpeff_json",
                    "uniquify": False,
                },
            },
            # "vep_extract": {
            #     "type": "python",
            #     "name": "vep_extract",
            #     "description": "HGVS nomenclatures from VEP annotation",
            #     "comment": "Extract HGVS nomenclatures from VEP annotation (field CSQ) and create new INFO fields with prefix 'vep_' (e.g. vep_hgvs, vep_impact, vep_gene_name...). This calculation parses the CSQ field from VEP annotations, extracting relevant information such as HGVS nomenclatures, impact, gene name, etc., and creates new INFO fields with a 'vep_' prefix for easier access and downstream analysis.",
            #     "available": True,
            #     "function_name": "calculation_snpeff_extract",
            #     "function_params": {
            #         "snpeff_field": "CSQ",
            #         "snpeff_hgvs": "vep_hgvs",
            #         "snpeff_explode": "vep_",
            #         "snpeff_json": "vep_json",
            #         "uniquify": False,
            #     },
            # },
            "snpeff_hgvs": {
                "type": "python",
                "name": "snpeff_hgvs",
                "description": "HGVS nomenclatures from snpEff annotation",
                "comment": "Extract HGVS nomenclatures from snpEff annotation (field ANN) and create new INFO field 'snpeff_hgvs'. This calculation specifically targets the extraction of HGVS nomenclatures from the ANN field of snpEff annotations, creating a new INFO field named 'snpeff_hgvs' that contains the extracted HGVS information for easier access and downstream analysis.",
                "available": True,
                "function_name": "calculation_snpeff_extract",
                "function_params": {
                    "snpeff_field": "ANN",
                    "snpeff_hgvs": "snpeff_hgvs",
                    "snpeff_explode": None,
                    "snpeff_json": None,
                    "uniquify": False,
                }
            },
            "snpeff_ann_explode": {
                "type": "python",
                "name": "snpeff_ann_explode",
                "description": "Explode snpEff annotations",
                "comment": "Explode snpEff annotations from the ANN field, creating new INFO fields with prefix 'snpeff_' for each annotation. This calculation takes the ANN field from snpEff annotations, which may contain multiple annotations for a single variant, and explodes it so that each annotation is represented in separate INFO fields with a 'snpeff_' prefix. This allows for easier access to individual annotations and facilitates downstream analysis.",
                "available": True,
                "function_name": "calculation_snpeff_extract",
                "function_params": {
                    "snpeff_field": "ANN",
                    "snpeff_hgvs": None,
                    "snpeff_explode": "snpeff_",
                    "snpeff_json": None,
                    "uniquify": False,
                }
            },
            "snpeff_ann_explode_uniquify": {
                "type": "python",
                "name": "snpeff_ann_explode_uniquify",
                "description": "Explode snpEff annotations with uniquify values",
                "comment": "Explode snpEff annotations from the ANN field, creating new INFO fields with prefix 'snpeff_uniquify_' for each annotation and ensuring unique values. This calculation takes the ANN field from snpEff annotations, which may contain multiple annotations for a single variant, and explodes it so that each annotation is represented in separate INFO fields with a 'snpeff_uniquify_' prefix, ensuring that only unique values are retained. This allows for easier access to individual annotations and facilitates downstream analysis.",
                "available": True,
                "function_name": "calculation_snpeff_extract",
                "function_params": {
                    "snpeff_field": "ANN",
                    "snpeff_hgvs": None,
                    "snpeff_explode": "snpeff_uniquify_",
                    "snpeff_json": None,
                    "uniquify": True,
                }
            },
            "snpeff_ann_explode_json": {
                "type": "python",
                "name": "snpeff_ann_explode_json",
                "description": "Explode snpEff annotations in JSON format",
                "comment": "Explode snpEff annotations from the ANN field, creating a new INFO field 'snpeff_json' in JSON format that contains all annotations. This calculation takes the ANN field from snpEff annotations, which may contain multiple annotations for a single variant, and explodes it into a structured JSON format stored in a new INFO field named 'snpeff_json'. This allows for easier access to all annotations in a single field and facilitates downstream analysis.",
                "available": True,
                "function_name": "calculation_snpeff_extract",
                "function_params": {
                    "snpeff_field": "ANN",
                    "snpeff_hgvs": None,
                    "snpeff_explode": None,
                    "snpeff_json": "snpeff_json",
                    "uniquify": True,
                },
            },
            "NOMEN": {
                "type": "python",
                "name": "NOMEN",
                "description": "NOMEN information (e.g. NOMEN, CNOMEN, PNOMEN...) from HGVS nomenclature field (see parameters help)",
                "comment": "Extract NOMEN information (e.g. NOMEN, CNOMEN, PNOMEN...) from HGVS nomenclature field (see parameters help). This calculation parses the HGVS nomenclature field to extract specific NOMEN information such as NOMEN, CNOMEN, PNOMEN, etc., based on the parameters provided. It creates new INFO fields with the extracted NOMEN information for easier access and downstream analysis.",
                "available": True,
                "function_name": "calculation_extract_nomen",
                "function_params": {}
            },
            "NOMEN_SNPEFF": {
                "type": "python",
                "name": "NOMEN_SNPEFF",
                "description": "NOMEN information (e.g. NOMEN, CNOMEN, PNOMEN...) from HGVS nomenclature field (see parameters help)",
                "comment": "Extract NOMEN information (e.g. NOMEN, CNOMEN, PNOMEN...) from HGVS nomenclature field (see parameters help). This calculation parses the HGVS nomenclature field to extract specific NOMEN information such as NOMEN, CNOMEN, PNOMEN, etc., specifically based on snpEff HGVS annotations field 'snpeff_hgvs'. It creates new INFO fields with the extracted NOMEN information for easier access and downstream analysis.",
                "available": True,
                "function_name": "calculation_extract_nomen",
                "function_params": {"hgvs_field": "snpeff_hgvs"},
            },
            "RECREATE_INFO_FIELDS": {
                "type": "python",
                "name": "RECREATE_INFO_FIELDS",
                "description": "Recreate INFO_tags, rename or remove tags",
                "comment": "Recreate INFO_tags, rename or remove tags. This calculation allows for the recreation of INFO fields by renaming existing tags or removing them based on specified parameters. It provides flexibility in managing INFO fields, enabling users to customize their VCF annotations according to their analysis needs.",
                "available": True,
                "function_name": "calculation_recreate_info_fields",
                "function_params": [],
            },
            "RENAME_INFO_FIELDS": {
                "type": "python",
                "name": "RENAME_INFO_FIELDS",
                "description": "Rename or remove INFO/tags",
                "comment": [
                    "Rename or remove INFO/tags. This calculation allows for the renaming of existing INFO fields or the removal of specific tags based on provided parameters. It helps in standardizing or customizing INFO fields in the VCF according to user requirements and facilitates downstream analysis. Use paramter section 'calculation_RENAME_INFO_FIELDS' in param.json to specify the INFO fields to rename or remove.",
                    "",
                    "This example will rename INFO field 'ENOMEN' to 'exon' and remove INFO fields 'SPiP', 'SPiP_Alt' and 'SPiP_distSS' from the 'variants' table:",
                    "",
                    "```json",
                    """{"fields_to_rename": {"ENOMEN": "exon", "SPiP": null, "SPiP_Alt": null, "SPiP_distSS": null }, "table": "variants"}""",
                    "```",
                ],
                "available": True,
                "function_name": "calculation_rename_info_fields",
                "function_params": [],
            },
            "FIND_SAMPLES": {
                "type": "python",
                "name": "FIND_SAMPLES",
                "description": "Number of sample that have a genotype for the variant (for multi sample VCF)",
                "comment": [
                    "This calculation counts the number of samples that have a genotype for a given variant in a multi-sample VCF file. It helps in assessing the presence of the variant across different samples.\n",
                    "Options for this calculation can be specified in the JSON parameter file, directly in the calculation parameters (see help.parameters.md). The parameter 'tags' allows for the specification of the INFO tags to be created, with the key being the tag name and the value being either 'count' or 'list' to indicate whether to count the samples or list them. For example, to create an INFO tag 'count_samples' that counts the number of samples and an INFO tag 'list_samples' that lists the samples, you can specify:\n",
                    "```json",
                    """{ "tags": {"count_samples": "count", "list_samples": "list"} }""",
                    "```\n",
                    "Otherwise, change the default tags in the parameter file, like in this example:\n", 
                    "```json",
                    """{ "tags": {"count_samples_for_variants": "count", "list_samples_for_variants": "list", "another_list_samples_tag": "list"} }""",
                    "```",
                ],
                "available": True,
                "function_name": "calculation_find_samples",
                "function_params": {},
            },
            
            "COUNT_SAMPLES": {
                "type": "python",
                "name": "COUNT_SAMPLES",
                "description": "Number of sample that have a genotype for the variant (for multi sample VCF)",
                "comment": "Number of sample that have a genotype for the variant (for multi sample VCF). This calculation counts the number of samples (fixed 'count_samples' field) that have a genotype for a given variant in a multi-sample VCF file. It helps in assessing the presence of the variant across different samples.",
                "available": True,
                "function_name": "calculation_find_samples",
                "function_params": {"tags": {"count_samples": "count"}},
            },
            "LIST_SAMPLES": {
                "type": "python",
                "name": "LIST_SAMPLES",
                "description": "List of samples that have a genotype for the variant (for multi sample VCF)",
                "comment": "List of samples that have a genotype for the variant (for multi sample VCF). This calculation provides a list of samples (fixed 'list_samples' field) that have a genotype for a given variant in a multi-sample VCF file. It helps in understanding which samples support the presence of the variant.",
                "available": True,
                "function_name": "calculation_find_samples",
                "function_params": {"tags": {"list_samples": "list"}},
            },
            "GENOTYPECONCORDANCE": {
                "type": "python",
                "name": "GENOTYPECONCORDANCE",
                "description": "Concordance of genotype for multi caller VCF",
                "comment": "Concordance of genotype for multi caller VCF. This calculation assesses the concordance of genotypes for a given variant across multiple callers in a multi-caller VCF file. It helps in evaluating the consistency of genotype calls and can be used to identify variants with high confidence based on agreement among different callers.",
                "available": True,
                "function_name": "calculation_genotype_concordance",
                "function_params": [],
            },
            "BARCODE": {
                "type": "python",
                "name": "BARCODE",
                "description": "BARCODE as VaRank tool",
                "comment": "BARCODE as VaRank tool. This calculation generates a BARCODE for each sample based on sample genotype (0 for reference, 1 for heterozygous, 2 for homozygous alternate).",
                "available": True,
                "function_name": "calculation_barcode",
                "function_params": [],
                "output_column_number": "1",
                "output_column_type": "Integer",
            },
            "BARCODEFAMILY": {
                "type": "python",
                "name": "BARCODEFAMILY",
                "description": "BARCODEFAMILY as VaRank tool",
                "comment": "BARCODEFAMILY as VaRank tool. This calculation generates a BARCODE for each family based on sample genotype (0 for reference, 1 for heterozygous, 2 for homozygous alternate) and family relationships.",
                "available": True,
                "function_name": "calculation_barcode_family",
                "function_params": ["BCF"],
            },
            "TRIO": {
                "type": "python",
                "name": "TRIO",
                "description": "Inheritance for a trio family",
                "comment": [
                    "Inheritance for a trio family. This calculation assesses the inheritance pattern of a variant within a trio family pedigree (father, mother, and child). It helps in understanding the transmission of genetic variants and identifying potential de novo mutations.",
                    "Options for this calculation can be specified in the JSON parameter file, directly in the calculation parameters (see help.parameters.md). The parameter 'trio_samples' allows for the specification of the sample names for father, mother, and child:\n",
                    "```json",
                    """{"father": "sample_father", "mother": "sample_mother", "child": "sample_child"}""",
                    "```",
                    "This parameter can also be specified in a JSON file (containing the trio sample names):\n",
                    "```bash",
                    """ "/path/to/trio_samples.json" """,
                    "```",
                    "If no pedigree information is provided, the calculation will use the first 3 samples in the VCF file."
                ],
                "available": True,
                "function_name": "calculation_trio",
                "function_params": [],
            },
            "VAF": {
                "type": "python",
                "name": "VAF",
                "description": "Variant Allele Frequency (VAF) harmonization",
                "comment": "Variant Allele Frequency (VAF) harmonization. This calculation normalizes the Variant Allele Frequency (VAF) across different samples and callers to ensure consistency in VAF representation. It helps in comparing VAF values across samples and callers by applying a harmonization process.",
                "available": True,
                "function_name": "calculation_vaf_normalization",
                "function_params": [],
            },
            "VAF_stats": {
                "type": "python",
                "name": "VAF_stats",
                "description": "Variant Allele Frequency (VAF) statistics",
                "comment": "Variant Allele Frequency (VAF) statistics. This calculation provides statistical measures for the Variant Allele Frequency (VAF) across different samples. It helps in understanding the distribution and variability of VAF values.",
                "available": True,
                "function_name": "calculation_genotype_stats",
                "function_params": {"info": "VAF"},
            },
            "DP_stats": {
                "type": "python",
                "name": "DP_stats",
                "description": "Depth (DP) statistics",
                "comment": "Depth (DP) statistics. This calculation provides statistical measures for the sequencing depth (DP) across different samples. It helps in understanding the coverage and reliability of variant calls.",
                "available": True,
                "function_name": "calculation_genotype_stats",
                "function_params": {"info": "DP"},
            },
            "INFO_TO_FORMAT": {
                "type": "python",
                "name": "INFO_TO_FORMAT",
                "description": "INFO to FORMAT conversion",
                "comment": [
                    "INFO to FORMAT conversion. This calculation converts INFO fields to FORMAT fields in the VCF file.\n",
                    "Options for this calculation can be specified in the JSON parameter file, directly in the calculation parameters (see help.parameters.md):\n",
                    "- 'annotation_fields': allows for the specification of the INFO fields to be converted to FORMAT fields, with the key being the INFO field name and the value being the FORMAT field name. For example, to convert INFO field 'count_samples' to FORMAT field 'CS', INFO field 'list_samples' to FORMAT field 'LS', and INFO field 'calling_quality' to FORMAT field 'CQ', you can specify:\n",
                    "- 'samples': allows for the specification of the samples to be included in the conversion. If not specified, all samples will be included.\n",
                    "- 'remove_info_fields': allows for the removal of the original INFO fields after conversion to FORMAT fields. If set to true, the original INFO fields will be removed from the VCF file.\n",
                    "Options example: \n",
                    "```json",
                    "{",
                    "    \"annotation_fields\": {",
                    "        \"count_samples\": null,",
                    "        \"list_samples\": \"LS\",",
                    "        \"calling_quality\": \"CQ\"",
                    "    },",
                    "    \"samples\": [\"sample1\", \"sample2\"],",
                    "    \"remove_info_fields\": true",
                    "}",
                    "```"
                    ],
                "available": True,
                "function_name": "calculation_info_to_format",
                "function_params": {},
            },
            "variant_id": {
                "type": "python",
                "name": "variant_id",
                "description": "Variant ID generated from variant position and type",
                "comment": [
                    "Variant ID generated from variant position and variation (ref and alt) and type (SVTYPE). This calculation creates a unique identifier for each variant, facilitating variant tracking and comparison across datasets.",
                    "Options for this calculation can be specified in the JSON parameter file, directly in the calculation parameters (see help.parameters.md). The parameter 'variant_id_tag' allows for the specification of the INFO tag to be used for the variant ID. By default, it uses 'variant_id', but it can be changed to any other INFO tag name as needed. For example, to use 'varid' as the INFO tag for the variant ID, you can specify:\n",
                    "```json",
                    """{ "variant_id_tag": "varid" }""",
                    "```",
                    "Other options as 'variant_id_tag_info' allows for the specification of the description in the VCF header, and 'keep_variant_id_tag_column' allows for keeping the variant ID tag column in the 'variants' table for further analysis.\n",
                    "As an axample of the JSON parameter file, you can specify:\n",
                    "```json",
                    """{ "variant_id_tag": "varid", "variant_id_tag_info": "Variant ID generated from variant position and type", "keep_variant_id_tag_column": true }""",
                    "```",
                ],
                "available": True,
                "function_name": "calculation_variant_id",
                "function_params": {},
            },
            "variant_id_varid": {
                "type": "python",
                "name": "variant_id_varid",
                "description": "Variant ID generated from variant position and type, using 'varid' as INFO tag",
                "comment": "Variant ID generated from variant position and variation (ref and alt) and type (SVTYPE). This calculation creates a unique identifier for each variant, facilitating variant tracking and comparison across datasets.",
                "available": True,
                "function_name": "calculation_variant_id",
                "function_params": {"variant_id_tag": "varid"},
            },
            # "transcripts_annotations": {
            #     "type": "python",
            #     "name": "transcripts_annotations",
            #     "description": "Perform transcripts annotations and generate a transcripts table/view (using JSON parameters file)",
            #     "comment": "Perform transcripts annotations and generate a transcripts table/view. This calculation allows for the inclusion of transcript annotations in both JSON and structured formats, providing flexibility in how transcript information is stored and accessed within the variant analysis pipeline. The specific format(s) used can be determined based on parameters specified in the JSON parameter file in 'transcripts' section (see help.parameters.md), or directly in the calculation parameters.",
            #     "available": True,
            #     "function_name": "calculation_transcripts_annotation",
            #     "function_params": [None, None],
            # },
            # "transcripts_annotations_json_format": {
            #     "type": "python",
            #     "name": "transcripts_json",
            #     "description": "Perform transcripts annotations and export into INFO field in JSON format (field 'transcripts_json')",
            #     "comment": "Perform transcripts annotations and generate a transcripts table/view. This calculation allows for the inclusion of transcript annotations in both JSON and structured formats, providing flexibility in how transcript information is stored and accessed within the variant analysis pipeline. The specific format(s) used can be determined based on parameters specified in the JSON parameter file in 'transcripts' section (see help.parameters.md). Then and add transcripts fields into INFO fields in structured format (field 'transcripts_ann'). This calculation allows for the inclusion of transcript annotations in a JSON format, facilitating the integration of transcript information into the variant analysis pipeline.",
            #     "available": True,
            #     "function_name": "calculation_transcripts_annotation",
            #     "function_params": ["transcripts_json", None],
            # },
            # "transcripts_annotations_struct_format": {
            #     "type": "python",
            #     "name": "transcripts_ann",
            #     "description": "Perform transcripts annotations and export into INFO field in structured format (field 'transcripts_ann')",
            #     "comment": "Perform transcripts annotations and generate a transcripts table/view. This calculation allows for the inclusion of transcript annotations in both JSON and structured formats, providing flexibility in how transcript information is stored and accessed within the variant analysis pipeline. The specific format(s) used can be determined based on parameters specified in the JSON parameter file in 'transcripts' section (see help.parameters.md). Then and add transcripts fields into INFO fields in structured format (field 'transcripts_ann'). This calculation allows for the inclusion of transcript annotations in a structured format, facilitating the integration of transcript information into the variant analysis pipeline.",
            #     "available": True,
            #     "function_name": "calculation_transcripts_annotation",
            #     "function_params": [None, "transcripts_ann"],
            # },
            "transcripts_annotations": {
                "type": "python",
                "name": "transcripts_annotations",
                "description": "Perform transcripts annotations and generate a transcripts table/view (using JSON parameters file)",
                "comment": "Perform transcripts annotations and generate a transcripts table/view. This calculation allows for the inclusion of transcript annotations in both JSON and structured formats, providing flexibility in how transcript information is stored and accessed within the variant analysis pipeline. The specific format(s) used can be determined based on parameters specified in the JSON parameter file in 'transcripts' section (see help.parameters.md), or directly in the calculation parameters.",
                "available": True,
                "function_name": "calculation_transcripts_annotation",
                "function_params": {},
            },
            "transcripts_annotations_json_format": {
                "type": "python",
                "name": "transcripts_json",
                "description": "Perform transcripts annotations and export into INFO field in JSON format (field 'transcripts_json')",
                "comment": "Perform transcripts annotations and generate a transcripts table/view. This calculation allows for the inclusion of transcript annotations in both JSON and structured formats, providing flexibility in how transcript information is stored and accessed within the variant analysis pipeline. The specific format(s) used can be determined based on parameters specified in the JSON parameter file in 'transcripts' section (see help.parameters.md). Then and add transcripts fields into INFO fields in structured format (field 'transcripts_ann'). This calculation allows for the inclusion of transcript annotations in a JSON format, facilitating the integration of transcript information into the variant analysis pipeline.",
                "available": True,
                "function_name": "calculation_transcripts_annotation",
                "function_params": {"info_json": "transcripts_json", "info_format": None},
            },
            "transcripts_annotations_struct_format": {
                "type": "python",
                "name": "transcripts_ann",
                "description": "Perform transcripts annotations and export into INFO field in structured format (field 'transcripts_ann')",
                "comment": "Perform transcripts annotations and generate a transcripts table/view. This calculation allows for the inclusion of transcript annotations in both JSON and structured formats, providing flexibility in how transcript information is stored and accessed within the variant analysis pipeline. The specific format(s) used can be determined based on parameters specified in the JSON parameter file in 'transcripts' section (see help.parameters.md). Then and add transcripts fields into INFO fields in structured format (field 'transcripts_ann'). This calculation allows for the inclusion of transcript annotations in a structured format, facilitating the integration of transcript information into the variant analysis pipeline.",
                "available": True,
                "function_name": "calculation_transcripts_annotation",
                "function_params": {"info_json": None, "info_format": "transcripts_ann"},
            },
            "transcripts_prioritization": {
                "type": "python",
                "name": "transcripts_prioritization",
                "description": "Prioritize transcripts with a prioritization profile (using JSON parameters file)",
                "comment": "Prioritize transcripts with a prioritization profile. This calculation allows for the prioritization of transcripts based on a predefined profile specified in the JSON parameter file, either in 'transcripts' section or directly in the calculation parameters (see help.parameters.md), helping to identify the most relevant transcripts for further analysis.",
                "available": True,
                "function_name": "calculation_transcripts_prioritization",
                "function_params": {"strict": False},
            },
            "transcripts_prioritization_strict": {
                "type": "python",
                "name": "transcripts_prioritization_strict",
                "description": "Prioritize transcripts with a prioritization profile (using JSON parameters file)",
                "comment": "Prioritize transcripts with a prioritization profile, ensuring that all needed annotations are present and considered. This calculation allows for the prioritization of transcripts based on a predefined profile specified in the JSON parameter file, either in 'transcripts' section or directly in the calculation parameters (see help.parameters.md), helping to identify the most relevant transcripts for further analysis.",
                "available": True,
                "function_name": "calculation_transcripts_prioritization",
                "function_params": {"strict": True},
            },
            "transcripts_export": {
                "type": "python",
                "name": "transcripts_export",
                "description": "Export transcripts table/view as a file (using JSON parameters file)",
                "comment": [
                    "Export transcripts table/view as a file. This calculation allows for the export of the transcripts table or view generated from transcript annotations to an external file format such as TSV or CSV. The parameters for this calculation can be specified in the JSON parameter file, either in 'transcripts' section or directly in the calculation parameters (see help.parameters.md), including options for the output file path, format, and any filters to apply during export."
                ],
                "available": True,
                "function_name": "calculation_transcripts_export",
                "function_params": {},
            },
            "variant_filter": {
                "type": "python",
                "name": "variant_filter",
                "description": "Filter variants based on specified criteria (using SQL parameters)",
                "comment": [
                    "Filter variants based on specified criteria. This calculation allows for the filtering of variants using SQL-like parameters specified in the JSON parameter file, either in the 'variants' section or directly in the calculation parameters (see help.parameters.md), including options for the filtering conditions and any additional criteria to apply during the filtering process."
                ],
                "available": True,
                "function_name": "calculation_variant_filter",
                "function_params": {},
            },
        }

        return config_default


    def get_operations_help(
        self,
        section: str = "calculation",
        operations_config_dict: dict = {},
        operations_config_file: str = None,
        show_calculations_md: str = None,
        format: str = "text",
    ) -> list:

        # Param
        param = self.get_param()

        # Param calculation config file
        if show_calculations_md is None:
            show_calculations_md = param.get(section, {}).get(
                "show_calculations_md", None
            )

        # Init
        operations_help = []
        operations_help_md = [
            "# Calculations",
            "List of available calculations",
        ]

        # operations
        operations = self.get_config_json(
            name="calculations",
            config_dict=operations_config_dict,
            config_file=operations_config_file,
        )

        # Sort operations by key
        operations = dict(sorted(operations.items(), key=lambda x: x[0].upper()))

        for op in operations:
            op_name = operations[op].get("name", op).upper()
            op_description = operations[op].get("description", op_name)
            op_comment = operations[op].get("comment", op_name)
            op_available = operations[op].get("available", False)
            if op_available:
                if isinstance(op_comment, list):
                    op_comment = "\n".join(op_comment).strip()
                operations_help.append(f"   {op_name}: {op_description}")
                operations_help_md.append(f"## {op_name}")
                operations_help_md.append(f"{op_description}")
                operations_help_md.append(f"{op_comment}")

        # insert header
        operations_help.insert(0, "Available calculation operations:")

        # Write operations help in markdown file
        if show_calculations_md is not None:
            with open(show_calculations_md, "w") as f:
                f.write("\n\n".join(operations_help_md))

        # Return
        if format == "text":
            return operations_help
        elif format == "markdown":
            return [
                f"Markdown file generated with operations help: {show_calculations_md}"
            ]
        else:
            return None

    def calculation(
        self,
        section: str = "calculation",
        operations: dict = {},
        operations_config_dict: dict = {},
        operations_config_file: str = None,
    ) -> None:
        """
        It takes a list of operations, and for each operation, it checks if it's a python or sql
        operation, and then calls the appropriate function

        param json example:
            "calculation": {
                "NOMEN": {
                    "options": {
                        "hgvs_field": "hgvs"
                    },
                "middle" : null
            }
        """

        # Param
        param = self.get_param()

        # CHeck operations config file
        if operations_config_file is None:
            operations_config_file = param.get(section, {}).get(
                "calculation_config", None
            )

        # operations config
        operations_config = self.get_config_json(
            name="calculations",
            config_dict=operations_config_dict,
            config_file=operations_config_file,
        )

        # Upper keys
        operations_config = {k.upper(): v for k, v in operations_config.items()}

        # Calculations

        # Operations from param
        operations = param.get(section, {}).get("calculations", operations)

        # Quick calculation - add
        if param.get("calculations", None):

            # List of operations
            calculations_list = [
                value.strip() for value in param.get("calculations", "").split(",")
            ]

            # Log
            log.info(f"Quick Calculations:")
            for calculation_key in calculations_list:
                log.info(f"   {calculation_key}")

            # Create tmp operations (to keep operation order)
            operations_tmp = {}
            for calculation_operation in calculations_list:
                if calculation_operation.upper() not in operations_tmp:
                    # log.debug(
                    #     f"{calculation_operation}.upper() not in {operations_tmp}"
                    # )
                    operations_tmp[calculation_operation.upper()] = {}
                    add_value_into_dict(
                        dict_tree=operations_tmp,
                        sections=[
                            calculation_operation.upper(),
                        ],
                        value=operations.get(calculation_operation.upper(), {}),
                    )
            # Add operations already in param
            for calculation_operation in operations:
                if calculation_operation not in operations_tmp:
                    operations_tmp[calculation_operation] = operations.get(
                        calculation_operation, {}
                    )

            # Update operations in param
            operations = operations_tmp

        # Operations for calculation
        if not operations:
            operations = param.get(section, {}).get("calculations", {})

        if operations:
            log.info(f"Calculations...")

        # Count number of variants
        try:
            nb_variants = self.get_query_to_df(
                f"SELECT count(1) AS count FROM (SELECT 1 FROM {self.get_table_variants()} LIMIT 1)"
            )["count"].tolist()[0]
        except Exception as e:
            log.debug(f"Error counting variants: {e}")
            return None

        # Init
        # To store operation params for each dest table for update merge
        operation_params = {}
        perform_update_all_calculations = False  # Disabled due to operation dependances

        # For each operations
        for operation_name in operations:

            operation_name = operation_name.upper()

            if operation_name not in [""]:

                if operation_name in operations_config:

                    # Log
                    log.info(f"Calculation '{operation_name}'...")

                    # Get operation config
                    operation = operations_config[operation_name]
                    operation_type = operation.get("type", "sql")
                    operation_allow_empty = operation.get("allow_empty", False)

                    if operation_allow_empty or nb_variants > 0:

                        # Python process
                        if operation_type == "python":
                            self.calculation_process_function(section=section, operation=operation)

                        # SQL process
                        elif operation_type == "sql":

                            # Retrive parrams for operation
                            operation_dest_table, operation_param = (
                                self.calculation_process_sql(
                                    section=section, operation=operation, operation_name=operation_name
                                )
                            )

                            if perform_update_all_calculations:

                                # Create list of params for each dest table
                                if operation_dest_table not in operation_params:
                                    operation_params[operation_dest_table] = []
                                # Append param
                                operation_params[operation_dest_table].append(
                                    operation_param.copy()
                                )

                            else:

                                # Perform update for all calculations
                                log.debug(
                                    f"Process calculations for '{operation_name}'..."
                                )
                                self.update_table(
                                    dest_table=operation_dest_table,
                                    sources=[operation_param],
                                    physical_order=True,
                                    # chunk_size=10000000,
                                    force_strategy=None,
                                    chromosomes=None,
                                )

                                # Clean temp operation view
                                self.remove_tables_or_views(
                                    tables=[operation_param.get("table")]
                                )

                        # Fail process
                        else:
                            log.error(
                                f"Operations config: Type '{operation_type}' NOT available"
                            )
                            raise ValueError(
                                f"Operations config: Type '{operation_type}' NOT available"
                            )

                    else:
                        log.info(
                            f"Calculation '{operation_name}': aborded - no variants"
                        )

                else:
                    log.error(
                        f"Operations config: Calculation '{operation_name}' NOT available"
                    )
                    raise ValueError(
                        f"Operations config: Calculation '{operation_name}' NOT available"
                    )

        # Perform update for all calculations
        if perform_update_all_calculations:
            log.info("Process calculations...")
            for operation_table_dest in operation_params:
                self.update_table(
                    dest_table=operation_table_dest,
                    sources=operation_params[operation_table_dest],
                    physical_order=True,
                    # chunk_size=10000000,
                    force_strategy=None,
                )

    def calculation_process_sql(
        self, section:str = "calculation", operation: dict = {}, operation_name: str = "unknown"
    ) -> tuple:
        """
        The `calculation_process_sql` function takes in a mathematical operation as a string and
        performs the operation, updating the specified table with the result.

        :param operation: The `operation` parameter is a dictionary that contains information about the
        mathematical operation to be performed. It includes the following keys:
        :type operation: dict
        :param operation_name: The `operation_name` parameter is a string that represents the name of
        the mathematical operation being performed. It is used for logging and error handling purposes,
        defaults to unknown
        :type operation_name: str (optional)
        """

        # Operation infos
        operation_name = operation.get("name", "unknown")
        log.debug(f"process SQL {operation_name}")
        output_column_name = operation.get("output_column_name", operation_name)
        output_column_number = operation.get("output_column_number", ".")
        output_column_type = operation.get("output_column_type", "String")
        output_column_description = operation.get(
            "output_column_description", f"{operation_name} operation"
        )
        operation_query = operation.get("operation_query", None)
        if isinstance(operation_query, list):
            operation_query = " ".join(operation_query)
        operation_info_fields = operation.get("info_fields", [])
        operation_info_fields_check = operation.get("info_fields_check", False)
        operation_table = operation.get(
            "table", self.get_table_variants(clause="alter")
        )
        operation_table_source = operation.get("table_source", operation_table)
        operation_table_dest = operation.get("table_dest", operation_table)
        operation_table_key = operation.get(
            "table_key", ["#CHROM", "POS", "REF", "ALT"]
        )

        if operation_query:

            # Info fields check
            operation_info_fields_check_result = True
            if operation_info_fields_check:
                header_infos = self.get_header().infos
                for info_field in operation_info_fields:
                    operation_info_fields_check_result = (
                        operation_info_fields_check_result
                        and info_field in header_infos
                    )

            # If info fields available
            if operation_info_fields_check_result:

                # Create VCF header field
                vcf_reader = self.get_header()
                vcf_reader.infos[output_column_name] = vcf.parser._Info(
                    output_column_name,
                    output_column_number,
                    output_column_type,
                    output_column_description,
                    "howard calculation",
                    "0",
                    self.code_type_map.get(output_column_type),
                )

                # Create view
                table_view_name = "calculation_view_" + str(get_random())
                table_view_name = self.create_annotations_view(
                    table=operation_table_source,
                    view=table_view_name,
                    view_type="view",
                    view_mode="explore",
                    fields=operation_info_fields + ["INFO"],
                    fields_needed=operation_table_key,
                    info_prefix_column="",
                    info_struct_column="INFOS",
                    drop_view=True,
                )

                # Table key construct
                clause_key = []
                for key in operation_table_key:
                    clause_key.append(
                        f""" {operation_table_dest}."{key}" = table_view."{key}" """
                    )

                # Create table with calculation
                calculation_view_name = "calculation_view_" + str(
                    get_random()
                )

                try:
                    # # OLD Query
                    # query_create_view = f"""
                    #     CREATE TABLE {calculation_view_name} AS
                    #     SELECT * FROM (
                    #         SELECT {", ".join([f'"{k}"' for k in operation_table_key])},
                    #             CASE
                    #                 WHEN TRY_CAST(({operation_query}) AS VARCHAR) IS NOT NULL
                    #                 THEN
                    #                     concat(
                    #                             '{output_column_name}=',
                    #                             TRY_CAST(({operation_query}) AS VARCHAR)
                    #                         )
                    #                 ELSE NULL 
                    #                 END AS INFO
                    #         FROM {table_view_name} AS v
                    #         )
                    #     WHERE INFO IS NOT NULL
                    # """
                    # NEW Query
                    query_create_view = f"""
                        CREATE TABLE {calculation_view_name} AS
                        WITH base AS (
                            SELECT
                                {", ".join([f'"{k}"' for k in operation_table_key])},
                                TRY_CAST(({operation_query}) AS VARCHAR) AS operation_value
                            FROM {table_view_name} v
                        )
                        SELECT
                            {", ".join([f'"{k}"' for k in operation_table_key])},
                            concat('{output_column_name}=', operation_value) AS INFO
                        FROM base
                        WHERE operation_value IS NOT NULL;
                    """
                    # log.debug(f"query_create_view={query_create_view}")
                    log.debug("Create calculation view...")
                    self.get_connexion().execute(query_create_view)

                    # Clean temp annotation view
                    self.remove_tables_or_views(tables=[table_view_name])

                except Exception as e:
                    log.error(f"Error creating calculation view: {e}")
                    msg_err = f"Operations config: Calculation '{operation_name}' query failed"
                    log.error(msg_err)
                    raise ValueError(msg_err)

                # update table
                return operation_table_dest, {
                    "table": calculation_view_name,
                    "join_keys": operation_table_key,
                    "columns": {
                        "INFO": {
                            "columns": ["INFO"],
                            "mode": "append",
                            "separator": ";",
                        }
                    },
                }

            else:
                msg_err = f"Operations config: Calculation '{operation_name}' DOES NOT contain all mandatory fields {operation_info_fields}"
                log.error(msg_err)
                raise ValueError(msg_err)

        else:
            msg_err = (
                f"Operations config: Calculation '{operation_name}' query NOT defined"
            )
            log.error(msg_err)
            raise ValueError(msg_err)

    def calculation_process_function(
        self, section: str = "calculation", operation: dict = None, operation_name: str = None
    ) -> None:
        """
        The `calculation_process_function` takes in an operation dictionary and performs the specified
        function with the given parameters.

        :param operation: The `operation` parameter is a dictionary that contains information about the
        operation to be performed. It has the following keys:
        :type operation: dict
        :param operation_name: The `operation_name` parameter is a string that represents the name of
        the operation being performed. It is used for logging purposes, defaults to unknown
        :type operation_name: str (optional)
        """

        if operation_name is None:
            operation_name = operation["name"]
        log.debug(f"process Python {operation_name}")
        function_name = operation["function_name"]
        function_params = operation["function_params"]

        # getattr(self, function_name)(*function_params)
        log.debug(f"function_params={function_params}")
        if function_params is not None:
            if isinstance(function_params, dict):
                # Add operation_name in function_params
                if "operation_name" in function_params:
                    log.warning(
                        f"operation_name already in function_params for {operation_name}, it will be overwritten"
                    )
                function_params["operation_name"] = operation_name
                function_params["section"] = section
                # Dict -> **
                getattr(self, function_name)(**function_params)
            elif isinstance(function_params, list):
                # List -> *
                getattr(self, function_name)(*function_params)
            else:
                # Not dict neither list
                log.error(
                    f"Param function type not supported for {operation_name}: {type(function_params)}"
                )
                raise TypeError(f"Function parameters must be a dict or a list.")
        else:
            # No parameters
            getattr(self, function_name)()

    def get_operation_params(
        self,
        section: str = "calculation",
        operation_params: dict = None,
        operation_name: str = None,
    ):

        # Operation name
        operation_name = operation_params.get("operation_name", operation_name)

        # Operation param for test
        operation_params_test = operation_params.copy()
        operation_params_test.pop("operation_name", None)
        operation_params_test.pop("section", None)

        # Param from JSON parameters file
        param = self.get_param()

        # Param calculations with lower keys for case-insensitive access
        param_calculations = {
            key.lower(): value
            for key, value in param.get(section, {}).get("calculations", {}).items()
        }

        # Param from operation name
        try:
            param_operation = (
                param_calculations
                .get(operation_name.lower(), {})
            )
        except Exception as e:
            # log.error(f"Error retrieving operation parameters for {operation_name}: {e}")
            param_operation = None

        # If no operation params, try to get from JSON parameters file
        if operation_params_test is None or len(operation_params_test) == 0:
            operation_params = param_operation

        if operation_params is None:
            operation_params = {}

        # Return
        return operation_params, operation_name

    def calculation_variant_id(
        self,
        section: str = "calculation",
        variant_id_tag: str = None,
        variant_id_tag_info: str = None,
        keep_variant_id_tag_column: bool = None,
        **kwargs,
    ) -> None:
        """
        The function `calculation_variant_id` adds a variant ID annotation to a VCF file header and
        updates the INFO field of a variants table with the variant ID.
        """

        log.debug(f"Calculation variant ID...")

        operation_params, _ = self.get_operation_params(
            section=section, operation_params=kwargs, operation_name="variant_id"
        )

        ### Parameters for variant ID calculation

        # variant_id annotation field
        variant_id_tag = (
            operation_params.get("variant_id_tag")
            or variant_id_tag
            or self.get_variant_id_column()
        )

        # variant_id_tag_info
        if variant_id_tag_info is None:
            variant_id_tag_info = (
                operation_params.get("variant_id_tag_info")
                or variant_id_tag_info
                or "howard variant ID annotation"
            )

        # keep_variant_id_tag_column
        keep_variant_id_tag_column = (
            operation_params.get("keep_variant_id_tag_column")
            or keep_variant_id_tag_column
            or False
        )

        # variant_id_tag_number
        variant_id_tag_number = 1

        # variant_id_tag_type
        variant_id_tag_type = "String"

        # variant_id_tag_source
        variant_id_tag_source = "howard calculation"

        # variant_id_tag_version
        variant_id_tag_version = "0"

        # variant_id_tag_type_code
        variant_id_tag_type_code = self.code_type_map.get(variant_id_tag_type)

        ### Get variant ID column name and generate ids if not exist

        # get variant_id_column
        variant_id_column = self.get_variant_id_column(variant_id_column=variant_id_tag)

        # Add columns for calculation
        added_columns = [variant_id_column]

        ### Get header and add variant_id to header

        # Header
        vcf_reader = self.get_header()

        # Add variant_id to header
        vcf_reader.infos[variant_id_tag] = vcf.parser._Info(
            variant_id_tag,
            variant_id_tag_number,
            variant_id_tag_type,
            variant_id_tag_info,
            variant_id_tag_source,
            variant_id_tag_version,
            variant_id_tag_type_code,
        )

        ### Prepare update query to add variant_id into INFO field of variants table

        # Variants table
        table_variants = self.get_table_variants()

        # Update
        sql_update = f"""
            UPDATE {table_variants}
            SET "INFO" = 
                concat(
                    CASE
                        WHEN "INFO" IS NULL OR "INFO" IN ('','.')
                        THEN ''
                        ELSE concat("INFO", ';')
                    END,
                    '{variant_id_tag}=',
                    "{variant_id_column}"
                )
        """
        self.get_connexion().execute(sql_update)

        # Remove added columns
        if not keep_variant_id_tag_column:
            for added_column in added_columns:
                self.drop_column(column=added_column)

    def calculation_snpeff_extract(
        self,
        section: str = "calculation",
        snpeff_field: str = "ANN",
        snpeff_hgvs: str = "snpeff_hgvs",
        snpeff_explode: str = "snpeff_",
        snpeff_json: str = "snpeff_json",
        uniquify: bool = True,
        **kwargs
    ) -> None:
        """
        This function extracts SnpEff annotations from the specified field in the VCF file and processes them according to the provided parameters. The annotations can be exploded into separate rows, converted into JSON format, and/or ensured to be unique. The processed annotations are then added to the VCF file with the specified prefixes.

        Args:
            snpeff_field (str): The annotation field in the VCF file to extract SnpEff annotations from. Default is "ANN".
            snpeff_hgvs (str): The prefix for the HGVS annotations extracted from SnpEff. Default is "snpeff_hgvs".
            snpeff_explode (str): The prefix for the exploded annotations. Default is "snpeff_".
            snpeff_json (str): The prefix for the JSON annotations. Default is "snpeff_json".
            uniquify (bool): Whether to ensure unique annotations. Default is True.

        Returns:
            None

        """

        log.debug(f"Calculation snpEff extract...")

        operation_params, _ = self.get_operation_params(
            section=section, operation_params=kwargs, operation_name="snpeff_extract"
        )

        ### Parameters for snpEff extraction

        # snpeff field
        snpeff_field = (
            operation_params.get("snpeff_field")
            or snpeff_field
            or "ANN"
        )

        # snpeff_hgvs
        snpeff_hgvs = (
            operation_params.get("snpeff_hgvs")
            or snpeff_hgvs
            or "snpeff_hgvs"
        )

        # snpeff_explode
        snpeff_explode = (
            operation_params.get("snpeff_explode")
            or snpeff_explode
            or "snpeff_"
        )

        # snpeff_json
        snpeff_json = (
            operation_params.get("snpeff_json")
            or snpeff_json
            or "snpeff_json"
        )

        # uniquify
        uniquify = (
            operation_params.get("uniquify")
            or uniquify
            or True
        )

        # Variants table
        table_variants = self.get_table_variants()

        # Header
        vcf_reader = self.get_header()

        # Log
        log.info(f"Extract snpEff annotations")

        # If snpeff_field exists
        if snpeff_field in vcf_reader.infos:

            # Log
            log.info(f"Extract snpEff annotations - from INFO/Tag '{snpeff_field}'")

            # Create view
            view_name = "snpeff_hgvs_" + str(random.randint(1000, 9999))
            view_infos = self.annotation_format_to_table(
                annotation_field=snpeff_field,
                annotation_id="Feature_ID",
                view_name=view_name,
                column_rename={},
                column_clean=False,
                column_case=None,
            )
            view_name = view_infos[0]

            # Describe
            sql_describe = f"""
                SELECT *
                FROM (
                    DESCRIBE {view_name}
                )
                WHERE column_name NOT IN ('#CHROM', 'POS', 'REF', 'ALT', 'INFO')
            """
            sql_describe_result = self.get_query_to_df(sql_describe)

            # Create dict of snpEff annotations
            annotation_dict = {}
            for _, annotation in sql_describe_result.iterrows():

                # Process values for dict
                annotation_name = annotation.column_name
                annotation_clean = clean_annotation_field(name=annotation_name)
                annotation_type = annotation.column_type
                annotation_type_vcf = code_type_map_to_vcf.get(
                    annotation_type, "String"
                )
                annotation_column = f'"{annotation_name}"'
                annotation_number = 1
                if annotation_type_vcf in ["Flag"]:
                    annotation_number = 0
                elif annotation_name in ["Annotation"]:
                    annotation_number = "."
                    annotation_column = (
                        f"""replace(CAST("{annotation_name}" AS VARCHAR), '&', ',')"""
                    )
                elif annotation_name in ["Distance"]:
                    annotation_column = f"""string_split(CAST("{annotation_name}" AS VARCHAR), '.')[1]"""

                    annotation_number = 1
                annotation_desc = f"snpEff annotation '{annotation_name}'"

                # Create dict
                annotation_dict[annotation_name] = {
                    "name": annotation_name,
                    "id": annotation_clean,
                    "number": annotation_number,
                    "type": annotation_type_vcf,
                    "desc": annotation_desc,
                    "column": annotation_column,
                }

            # update clauses
            sql_clauses = []

            # Prepare sql update
            if snpeff_json is not None:

                # Log
                log.info(
                    f"Extract snpEff annotations - into INFO/tag '{snpeff_json}' in JSON format"
                )

                # Add snpeff_hgvs to header
                vcf_reader.infos[snpeff_json] = vcf.parser._Info(
                    snpeff_json,
                    1,
                    "String",
                    "snpEff annotation in JSON format",
                    "howard calculation",
                    "0",
                    self.code_type_map.get("String"),
                )

                # Prepare annotations
                sql_from_select_annotation_list = []
                for annotation in annotation_dict.values():
                    sql_from_select_annotation_list.append(
                        f""" '{annotation.get("id")}', {annotation.get("column")} """
                    )

                # Add snpeff JSON to header
                vcf_reader.infos[snpeff_hgvs] = vcf.parser._Info(
                    snpeff_hgvs,
                    ".",
                    "String",
                    "HGVS nomenclatures from snpEff annotation",
                    "howard calculation",
                    "0",
                    self.code_type_map.get("String"),
                )

                # Clause for INFO concat
                sql_info_concat = f"""
                    CASE
                        WHEN (INFO IS NULL OR INFO IN ('', '.')) OR (
                            SNPEFF_HGVS.json_data IS NULL
                        )
                        THEN INFO
                        ELSE concat(INFO, ';')
                    END,
                    CASE
                        WHEN SNPEFF_HGVS.json_data IS NOT NULL
                        THEN concat(
                            '{snpeff_json}=',
                            SNPEFF_HGVS.json_data
                        )
                    END
                    
                """

                # Clause for subquery
                sql_from_select = f"""
                    CASE
                        WHEN string_agg("Allele") IS NOT NULL
                        THEN
                            concat(
                                '[',
                                string_agg(
                                    json_object(
                                        {",".join(sql_from_select_annotation_list)}
                                    )::JSON
                                ),
                                ']'
                            )
                        ELSE NULL
                    END AS json_data
                """

                # Append clauses
                sql_clauses.append(
                    {
                        "sql_info_concat": sql_info_concat,
                        "sql_from_select": sql_from_select,
                    }
                )

            if snpeff_explode is not None:

                # Log
                log.info(
                    f"Extract snpEff annotations - into INFO/Tags separately with '{snpeff_explode}' prefix"
                )

                # Prepare annotations
                sql_info_concat_annotation_list = []
                sql_from_select_annotation_list = []
                for annotation in annotation_dict.values():

                    # Add snpeff_hgvs to header
                    annotation_id = f'{snpeff_explode}{annotation.get("id")}'
                    vcf_reader.infos[annotation_id] = vcf.parser._Info(
                        annotation_id,
                        annotation.get("number"),
                        annotation.get("type"),
                        annotation.get("desc"),
                        "howard calculation",
                        "0",
                        self.code_type_map.get(annotation.get("type")),
                    )

                    # Log
                    log.info(
                        f"Extract snpEff annotations - into INFO/Tags separately with '{snpeff_explode}' prefix - '{annotation_id}'"
                    )

                    # Clause for INFO concat for each annotation
                    sql_info_concat_annotation_list.append(
                        f""" 
                            CASE
                                WHEN SNPEFF_HGVS.{annotation.get("id")} IS NOT NULL AND CAST(SNPEFF_HGVS.{annotation.get("id")} AS STRING) NOT IN ('','.')
                                THEN concat('{snpeff_explode}{annotation.get("id")}=', CAST(SNPEFF_HGVS.{annotation.get("id")} AS STRING))
                            END
                        """
                    )

                    # Clause for subquery for each annotation
                    if uniquify:
                        sql_from_select_annotation_list.append(
                            f""" string_agg(DISTINCT CAST({annotation.get("column")} AS STRING), ',') AS '{annotation.get("id")}' """
                        )
                    else:
                        sql_from_select_annotation_list.append(
                            f""" string_agg(COALESCE(CAST({annotation.get("column")} AS STRING), '.'), ',') AS {annotation.get("id")} """
                        )

                # Clause for INFO concat
                sql_info_concat = f"""
                    CASE
                        WHEN (INFO IS NULL OR INFO IN ('', '.')) OR ("Allele" IS NULL)
                        THEN INFO
                        ELSE concat(INFO, ';')
                    END,
                    concat_ws(
                        ';',
                        {" , ".join(sql_info_concat_annotation_list)}
                    )
                """

                # Clause for subquery
                sql_from_select = " , ".join(
                    {" , ".join(sql_from_select_annotation_list)}
                )

                # Append clauses
                sql_clauses.append(
                    {
                        "sql_info_concat": sql_info_concat,
                        "sql_from_select": sql_from_select,
                    }
                )

            if snpeff_hgvs is not None:

                log.info(
                    f"Extract snpEff annotations - into INFO/Tags '{snpeff_hgvs}' with list of HGVS nomenclature"
                )

                # Add snpeff_hgvs to header
                vcf_reader.infos[snpeff_hgvs] = vcf.parser._Info(
                    snpeff_hgvs,
                    ".",
                    "String",
                    "HGVS nomenclatures from snpEff annotation",
                    "howard calculation",
                    "0",
                    self.code_type_map.get("String"),
                )

                # Clause for INFO concat
                sql_info_concat = f"""
                    CASE
                        WHEN (INFO IS NULL OR INFO IN ('', '.')) OR (SNPEFF_HGVS.hgvs IS NULL OR SNPEFF_HGVS.hgvs IN (''))
                        THEN INFO
                        ELSE concat(INFO, ';')
                    END,
                    CASE
                        WHEN SNPEFF_HGVS.hgvs IS NOT NULL AND SNPEFF_HGVS.hgvs NOT IN ('')
                        THEN concat('{snpeff_hgvs}=', SNPEFF_HGVS.hgvs)
                    END
                """

                # Clause for subquery
                sql_from_select = f"""
                    string_agg(
                        concat_ws(
                            ':',
                            "Gene_ID",
                            "Feature_ID",
                            CASE 
                                WHEN "Rank" IS NOT NULL
                                THEN concat('exon', split(CAST("Rank" AS VARCHAR), '/')[1])
                                ELSE NULL
                            END,
                            "HGVS.c",
                            "HGVS.p"
                        ),
                    ',') AS hgvs
                """

                # Append clauses
                sql_clauses.append(
                    {
                        "sql_info_concat": sql_info_concat,
                        "sql_from_select": sql_from_select,
                    }
                )

            # Update
            nb_update = 0
            for sql_clause_item in sql_clauses:

                # Nb update
                nb_update += 1

                # Query
                sql_update = f"""
                    UPDATE variants
                    SET INFO = concat(
                                    {sql_clause_item.get("sql_info_concat")}
                                    )
                    FROM (
                        SELECT "#CHROM", "POS", "REF", "ALT",
                            {sql_clause_item.get("sql_from_select")}
                        FROM {view_name}
                        GROUP BY "#CHROM", "POS", "REF", "ALT"
                        ) AS SNPEFF_HGVS
                    WHERE {table_variants}."#CHROM" = SNPEFF_HGVS."#CHROM"
                    AND {table_variants}."POS" = SNPEFF_HGVS."POS"
                    AND {table_variants}."REF" = SNPEFF_HGVS."REF"
                    AND {table_variants}."ALT" = SNPEFF_HGVS."ALT"
                """

                # Log
                log.info(
                    f"Extract snpEff annotations - Process [{nb_update}/{len(sql_clauses)}]"
                )

                # Process query
                self.conn.execute(sql_update)

            # Delete view
            sql_drop_view = f"""
                DROP VIEW {view_name}
            """
            self.conn.execute(sql_drop_view)

        else:

            log.warning(
                f"Extract snpEff annotations - No snpEff annotation '{snpeff_field}'. Please Anotate with snpEff before use this calculation option"
            )

    def calculation_extract_nomen(
        self,
        section: str = "calculation",
        hgvs_field: str = None,
        hgvs_fields: list = None,
        uniquify_hgvs: bool = None,
        **kwargs
    ) -> None:
        """
        Extracts the HGVS nomenclature from the provided field and calculates the NOMEN patterns.

        This function performs the following steps:
        1. Retrieves extra information fields and constructs the NOMEN pattern.
        2. Splits the NOMEN pattern based on predefined separators.
        3. Constructs SQL queries to parse and extract various components of the NOMEN pattern.
        4. Calculates scores for each variant based on the extracted NOMEN components.
        5. Creates a temporary table to store the results of the NOMEN extraction and scoring.
        6. Updates the main variants table with the extracted NOMEN information.

        Args:
            hgvs_field (str, optional): The field containing the HGVS nomenclature to be extracted. Defaults to None.
            hgvs_fields (list, optional): A list of fields containing HGVS nomenclature to be extracted. Defaults to None.
            uniquify_hgvs (bool, optional): Whether to order the HGVS nomenclature in the output. Defaults to False.

        Returns:
            None

        Raises:
            Any exceptions raised during the execution of the SQL queries or file operations.

        Example:
            self.calculation_extract_nomen(hgvs_field="hgvs_column")
            self.calculation_extract_nomen(hgvs_fields=["hgvs_column1", "hgvs_column2"])
        """

        # NOMEN structure
        nomen_dict = {
            "NOMEN": "NOMEN hgvs nomenclature considered as reference hgvs",
            "CNOMEN": "CNOMEN hgvs nomenclature at DNA level related to a transcript (TNOMEN)",
            "RNOMEN": "RNOMEN hgvs nomenclature at RNA level related to a transcript (TNOMEN)",
            "NNOMEN": "NNOMEN hgvs nomenclature for non-coding variant",
            "PNOMEN": "PNOMEN hgvs nomenclature at Protein level related to a transcript (TNOMEN)",
            "UPNOMEN": "UPNOMEN hgvs nomenclature at Protein level as uncertain related to a transcript (TNOMEN)",
            "TVNOMEN": "TVNOMEN hgvs transcript with version (if any) used (e.g. for CNOMEN and PNOMEN)",
            "TNOMEN": "TNOMEN hgvs transcript used (e.g. for CNOMEN and PNOMEN)",
            "VNOMEN": "VNOMEN hgvs transcript version used (e.g. for CNOMEN and PNOMEN)",
            "TPVNOMEN": "TPVNOMEN hgvs protein transcript with version (if any) used (e.g. for CNOMEN and PNOMEN)",
            "TPNOMEN": "TNOMEN hgvs protein transcript used (e.g. for CNOMEN and PNOMEN)",
            "TPVVNOMEN": "VNOMEN hgvs protein transcript version used (e.g. for CNOMEN and PNOMEN)",
            "ENOMEN": "ENOMEN hgvs exon nomenclature related to a transcript (TNOMEN)",
            "GNOMEN": "GNOMEN hgvs gene nomenclature related to a transcript (TNOMEN)",
        }

        # Param
        param = self.get_param()

        # Header
        vcf_reader = self.get_header()

        # Get HGVS fields
        # Keep compatibility with previous versions of Howard where only one hgvs_field could be specified
        if hgvs_field is None:
            hgvs_field = (
                param.get(section, {})
                .get("calculations", {})
                .get("NOMEN", {})
                .get("options", {})
                .get("hgvs_field", "hgvs")
            )
        if hgvs_fields is None:
            hgvs_fields = (
                param.get(section, {})
                .get("calculations", {})
                .get("NOMEN", {})
                .get("options", {})
                .get("hgvs_fields", hgvs_field)
            )
        if isinstance(hgvs_fields, str):
            hgvs_fields = hgvs_fields.split(",")

        # Get uniquify_hgvs option
        if uniquify_hgvs is None:
            uniquify_hgvs = (
                param.get(section, {})
                .get("calculations", {})
                .get("NOMEN", {})
                .get("options", {})
                .get("uniquify_hgvs", False)
            )

        # Get NOMEN pattern
        nomen_pattern = (
            param.get(section, {})
            .get("calculations", {})
            .get("NOMEN", {})
            .get("options", {})
            .get("pattern", None)
        )
        # default NOMEN pattern
        if nomen_pattern is None:
            nomen_pattern = "GNOMEN:TNOMEN:ENOMEN:CNOMEN:RNOMEN:NNOMEN:PNOMEN"

        if isinstance(nomen_pattern, str):
            nomen_patterns = {"NOMEN": nomen_pattern}
        elif isinstance(nomen_pattern, dict):
            nomen_patterns = nomen_pattern
        else:
            msg_err = f"NOMEN pattern '{nomen_pattern}' is not well formed"
            log.error(msg_err)
            raise ValueError(msg_err)

        # Get NOMEN pattern
        nomen_fields = (
            param.get(section, {})
            .get("calculations", {})
            .get("NOMEN", {})
            .get("options", {})
            .get("fields", None)
        )

        # default NOMEN pattern
        if nomen_fields is None:  # or nomen_fields == []:
            nomen_fields = list(nomen_dict.keys())

        # Remove "NOMEN" as patterns as separetly processed
        for nomen_pattern in nomen_patterns.keys():
            if nomen_pattern in nomen_fields:
                nomen_fields.remove(nomen_pattern)

        # transcripts list of preference sources
        transcripts_sources = {}

        # Get transcripts
        transcripts_file = (
            param.get(section, {})
            .get("calculations", {})
            .get("NOMEN", {})
            .get("options", {})
            .get("transcripts", None)
        )
        transcripts_file = full_path(transcripts_file)
        if transcripts_file:
            if os.path.exists(transcripts_file):
                transcripts_dataframe = transcripts_file_to_df(transcripts_file)
                transcripts_from_file = transcripts_dataframe.iloc[:, 0].tolist()
                transcripts_sources["file"] = transcripts_from_file
            else:
                msg_err = f"Transcript file '{transcripts_file}' does NOT exist"
                log.error(msg_err)
                raise ValueError(msg_err)

        # Get transcripts table
        transcripts_table = (
            param.get(section, {})
            .get("calculations", {})
            .get("NOMEN", {})
            .get("options", {})
            .get("transcripts_table", self.get_table_variants())
        )
        # Get transcripts column
        transcripts_column = (
            param.get(section, {})
            .get("calculations", {})
            .get("NOMEN", {})
            .get("options", {})
            .get("transcripts_column", None)
        )

        # Transcripts of preference source order
        transcripts_order = (
            param.get(section, {})
            .get("calculations", {})
            .get("NOMEN", {})
            .get("options", {})
            .get("transcripts_order", ["column", "file"])
        )

        # Transcripts from file
        transcripts = transcripts_sources.get("file", [])

        # Log
        log.info(f"Start NOMEN calculation configuration...")

        # Create annotation view
        annotations_view = "annotations_view_for_extract_nomen_" + str(
            get_random()
        )
        if transcripts_column is None:
            transcripts_column = "transcript"
        annotations_view = self.create_annotations_view(
            table=transcripts_table,
            view=annotations_view,
            view_type="table",
            view_mode="explore",
            fields=[transcripts_column] + hgvs_fields,
            fields_forced_as_varchar=True,
            info_prefix_column="",
        )

        # Construct NOMEN Pattern
        separators = [",", "(", ")", "|", ":", "[", "]", "{", "}"]
        regex_pattern = "|".join(map(re.escape, separators))

        # Init
        nomen_patterns_sql = {}

        for nomen_pattern_name, nomen_pattern in nomen_patterns.items():

            # Split NOMEN pattern
            split_nomen_pattern = re.split(rf"({regex_pattern})", nomen_pattern)

            # Construct SQL NOMEN Pattern
            nomen_pattern_sql_list = []
            nomen_info_previous = ""
            inside_parentheses = False
            inside_brackets = False
            inside_braces = False

            # Parse NOMEN pattern
            for i, nomen_info in enumerate(split_nomen_pattern):
                if nomen_info == "(":
                    inside_parentheses = True
                    nomen_info_previous += nomen_info
                elif nomen_info == ")":
                    inside_parentheses = False
                    if nomen_info_previous:
                        nomen_pattern_sql_list.append(nomen_info)
                    nomen_info_previous = ""
                elif nomen_info == "[":
                    inside_brackets = True
                    nomen_info_previous += nomen_info
                elif nomen_info == "]":
                    inside_brackets = False
                    if nomen_info_previous:
                        nomen_pattern_sql_list.append(nomen_info)
                    nomen_info_previous = ""
                elif nomen_info == "{":
                    inside_braces = True
                    nomen_info_previous += nomen_info
                elif nomen_info == "}":
                    inside_braces = False
                    if nomen_info_previous:
                        nomen_pattern_sql_list.append(nomen_info)
                    nomen_info_previous = ""
                elif nomen_info in separators:
                    nomen_info_previous += nomen_info
                else:
                    if nomen_info != "":
                        next_info = (
                            split_nomen_pattern[i + 1]
                            if i + 1 < len(split_nomen_pattern)
                            else ""
                        )
                        if next_info in separators:
                            if inside_parentheses or inside_brackets or inside_braces:
                                nomen_pattern_sql_list.append(f"""
                                        CASE
                                            WHEN {nomen_info} IS NOT NULL
                                            THEN concat('{nomen_info_previous}', {nomen_info}, '{next_info}')
                                        END
                                    """)
                            else:
                                nomen_pattern_sql_list.append(f"""
                                        CASE
                                            WHEN {nomen_info} IS NOT NULL
                                            THEN concat('{nomen_info_previous}', {nomen_info})
                                        END
                                    """)
                        else:
                            nomen_pattern_sql_list.append(f"""
                                    CASE
                                        WHEN {nomen_info} IS NOT NULL
                                        THEN concat('{nomen_info_previous}', {nomen_info})
                                    END
                                """)
                        nomen_info_previous = ""

                # Construcut NOMEN pattern for SQL
                nomen_pattern_sql = ", ".join(nomen_pattern_sql_list)

                # Add NOMEN pattern for SQL
                nomen_patterns_sql[nomen_pattern_name] = nomen_pattern_sql

        # Transcript source order and index window
        transcripts_order_length = len(transcripts_order)
        transcripts_order_window = 10000000
        try:
            index_transcript_selected = (
                transcripts_order_length - transcripts_order.index("column")
            ) * transcripts_order_window
        except:
            index_transcript_selected = 0
        try:
            index_transcript_prefered = (
                transcripts_order_length - transcripts_order.index("file")
            ) * transcripts_order_window
        except:
            index_transcript_prefered = 0

        # Transcripts rank
        if len(transcripts) >= 1:
            transcripts_pond = {
                transcript: len(transcripts) - rank
                for rank, transcript in enumerate(transcripts, start=0)
            }

            # Construct transcripts pond table
            transcripts_pond_table = "transcripts_pond_" + str(
                get_random()
            )
            transcripts_pond_df = pd.DataFrame(
                list(transcripts_pond.items()), columns=["transcript", "rank"]
            )
            self.execute_query(
                f"CREATE TABLE {transcripts_pond_table} AS SELECT * FROM transcripts_pond_df"
            )

            transcripts_pond_score_sql = f"""
                + CASE
                        WHEN TVNOMEN in (SELECT transcript FROM {transcripts_pond_table})
                            OR TNOMEN IN (SELECT transcript FROM {transcripts_pond_table})
                        THEN {index_transcript_prefered} + (
                                SELECT {transcripts_pond_table}.rank
                                FROM {transcripts_pond_table}
                                WHERE {transcripts_pond_table}.transcript = TVNOMEN
                                    OR {transcripts_pond_table}.transcript = TNOMEN
                                LIMIT 1
                            )
                        ELSE 0
                    END
            """
        else:
            transcripts_pond_score_sql = ""

        # NOMEN Patterns
        pattern_tvnomen = r".*[:]*([NX][MR]_[^:]*).*"
        pattern_tpvnomen = r".*[:]*([NX]P_[^:]*).*"
        pattern_cnomen = r".*[:]*([cgm]\.[^:]*).*"
        pattern_pnomen = r".*[:]*([p]\.[^:]*).*"
        pattern_nnomen = r".*[:]*([n]\.[^:]*).*"
        pattern_rnomen = r".*[:]*([r]\.[^:]*).*"
        pattern_enomen = r".*[:]*(exon[^:]*).*"

        # Check NOMEN fields length
        nomen_fields_select_sql = ""
        if len(nomen_fields) >= 1:
            nomen_fields_select_sql = ", ".join(nomen_fields) + ","
        else:
            nomen_fields_select_sql = ""

        # NOMEN patterns SQL select
        nomen_patterns_sql_select_list = []
        for nomen_pattern_name, nomen_pattern_sql in nomen_patterns_sql.items():
            nomen_patterns_sql_select_list.append(f"""
                concat({nomen_pattern_sql}) AS "{nomen_pattern_name}",
            """)
        nomen_patterns_sql_select = " ".join(nomen_patterns_sql_select_list)

        # NEW query nomen variants select to handle multiple HGVS fields and avoid duplicates with LIST_DISTINCT and FLATTEN
        hgvs_lists = ", ".join(
            [
                f"STRING_SPLIT(NULLIF(\"{hgvs_field}\"::VARCHAR, ''), ',')"
                for hgvs_field in hgvs_fields
            ]
        )

        # Construct the expression to combine HGVS lists from multiple fields, flatten them, and optionally ensure uniqueness
        hgvs_list_expr = f"""
            FLATTEN(
                LIST_VALUE(
                    {hgvs_lists}
                )
            )
        """

        # If uniquify_hgvs is True, wrap the expression with LIST_DISTINCT to ensure unique HGVS nomenclatures
        if uniquify_hgvs:
            hgvs_list_expr = f"""
                LIST_DISTINCT(
                    {hgvs_list_expr}
                )
            """

        query_nomen_variants_select = f"""
        SELECT
            "#CHROM",
            "POS",
            "REF",
            "ALT",
            "{transcripts_column}"::VARCHAR AS transcript,
            UNNEST(
                {hgvs_list_expr}
            ) AS nomen
        FROM {annotations_view}
        """

        # Query find NOMEN
        query_find_nomen = f"""
            WITH
            nomen_variants AS (
                {query_nomen_variants_select}
            ),
            decomposed_variants AS (
                SELECT
                    "#CHROM", "POS", "REF", "ALT",
                    "transcript",
                    -- TVNOMEN
                    NULLIF(regexp_extract(nomen, '{pattern_tvnomen}', 1), '') AS 'TVNOMEN',
                    CASE
                        WHEN array_length(string_split(regexp_extract(nomen, '{pattern_tvnomen}', 1), '.'), 1) >= 1
                        THEN NULLIF(string_split(regexp_extract(nomen, '{pattern_tvnomen}', 1), '.')[1], '')
                        ELSE NULL
                    END AS 'TNOMEN',
                    CASE
                        WHEN array_length(string_split(regexp_extract(nomen, '{pattern_tvnomen}', 1), '.'), 1) >= 2
                        THEN NULLIF(string_split(regexp_extract(nomen, '{pattern_tvnomen}', 1), '.')[2], '')
                        ELSE NULL
                    END AS 'VNOMEN',
                    -- TPVNOMEN
                    NULLIF(regexp_extract(nomen, '{pattern_tpvnomen}', 1), '') AS 'TPVNOMEN',
                    CASE
                        WHEN array_length(string_split(regexp_extract(nomen, '{pattern_tpvnomen}', 1), '.'), 1) >= 1
                        THEN NULLIF(string_split(regexp_extract(nomen, '{pattern_tpvnomen}', 1), '.')[1], '')
                        ELSE NULL
                    END AS 'TPNOMEN',
                    CASE
                        WHEN array_length(string_split(regexp_extract(nomen, '{pattern_tpvnomen}', 1), '.'), 1) >= 2
                        THEN IFNULL(string_split(regexp_extract(nomen, '{pattern_tpvnomen}', 1), '.')[2], '')
                        ELSE NULL
                    END AS 'TPVVNOMEN',
                    -- CPNR-NOMEN
                    NULLIF(regexp_extract(nomen, '{pattern_cnomen}', 1), '') AS 'CNOMEN',
                    NULLIF(regexp_extract(nomen, '{pattern_pnomen}', 1), '') AS 'PNOMEN',
                    NULLIF(regexp_extract(nomen, '{pattern_nnomen}', 1), '') AS 'NNOMEN',
                    NULLIF(regexp_extract(nomen, '{pattern_rnomen}', 1), '') AS 'RNOMEN',
                    -- Uncertain p.
                    NULLIF(
                        CASE
                            WHEN NULLIF(regexp_extract(nomen, '{pattern_pnomen}', 1), '') IS NOT NULL
                            THEN concat(
                                'p.',
                                '(',
                                string_split(regexp_extract(nomen, '{pattern_pnomen}', 1), '.')[2],
                                ')'
                            )
                        END
                    , '') AS 'UPNOMEN',
                    -- exon
                    NULLIF(regexp_extract(nomen, '{pattern_enomen}', 1), '') AS 'ENOMEN',
                    -- gene
                    CASE
                        WHEN NULLIF(regexp_extract(string_split(nomen, ':')[1], '{pattern_tvnomen}', 1), '') IS NOT NULL
                        OR NULLIF(regexp_extract(string_split(nomen, ':')[1], '{pattern_tpvnomen}', 1), '') IS NOT NULL
                        OR NULLIF(regexp_extract(string_split(nomen, ':')[1], '{pattern_cnomen}', 1), '') IS NOT NULL
                        OR NULLIF(regexp_extract(string_split(nomen, ':')[1], '{pattern_pnomen}', 1), '') IS NOT NULL
                        OR NULLIF(regexp_extract(string_split(nomen, ':')[1], '{pattern_nnomen}', 1), '') IS NOT NULL
                        OR NULLIF(regexp_extract(string_split(nomen, ':')[1], '{pattern_rnomen}', 1), '') IS NOT NULL
                        OR NULLIF(regexp_extract(string_split(nomen, ':')[1], '{pattern_enomen}', 1), '') IS NOT NULL
                        THEN NULL
                        ELSE NULLIF(string_split(nomen, ':')[1], '')
                    END AS 'GNOMEN'
                FROM nomen_variants
            ),
            scored_variants AS (
                SELECT
                    "#CHROM", "POS", "REF", "ALT",
                    "TNOMEN", "TVNOMEN", "VNOMEN",
                    "TPVNOMEN", "TPNOMEN", "TPVVNOMEN",
                    "CNOMEN", "PNOMEN", "NNOMEN", "RNOMEN", "UPNOMEN",
                    "ENOMEN", "GNOMEN",
                    -- Score calculation
                    0
                    + CASE WHEN CNOMEN IS NOT NULL THEN 1 ELSE 0 END
                    + CASE WHEN NNOMEN IS NOT NULL THEN 1 ELSE 0 END
                    + CASE WHEN RNOMEN IS NOT NULL THEN 1 ELSE 0 END
                    + CASE WHEN ENOMEN IS NOT NULL THEN 1 ELSE 0 END
                    + CASE WHEN PNOMEN IS NOT NULL THEN 1 ELSE 0 END
                    + CASE WHEN TPVNOMEN IS NOT NULL THEN 1 ELSE 0 END
                    + CASE WHEN regexp_matches(TVNOMEN, '^NM_.*') THEN 2 ELSE 0 END
                    + CASE WHEN regexp_matches(TVNOMEN, '^NR_.*') THEN 1 ELSE 0 END
                    -- Selected transcript
                    + CASE WHEN transcript IS NOT NULL AND (TVNOMEN == transcript OR TNOMEN == transcript) THEN {index_transcript_selected} ELSE 0 END
                    -- Preferend transcripts
                    {transcripts_pond_score_sql}
                    AS 'SCORE'
                FROM decomposed_variants
            )
                SELECT
                    "#CHROM", "POS", "REF", "ALT",
                    {nomen_fields_select_sql}
                    {nomen_patterns_sql_select}
                FROM (
                    SELECT *,
                        ROW_NUMBER() OVER (PARTITION BY "#CHROM", "POS", "REF", "ALT" ORDER BY SCORE DESC) AS rn
                    FROM scored_variants
                )
                WHERE rn = 1
        """
        nomen_annotations_view = "annotations_vies_for_extract_nomen_" + str(
            get_random()
        )
        query_find_nomen_create = f"""
        CREATE TABLE {nomen_annotations_view} AS (
            {query_find_nomen}
        )
        """
        # log.debug(f"query_devel={query_find_nomen}")
        log.info("Start NOMEN calculation...")
        self.execute_query(query=query_find_nomen_create)
        log.debug("Stop NOMEN calculation")

        # Explode NOMEN Structure and create SQL set for update
        sql_nomen_fields = []
        for nomen_field in list(nomen_patterns.keys()) + nomen_fields:

            # Description
            nomen_field_desc = nomen_dict.get(nomen_field, "howard calculation NOMEN")
            if nomen_field in list(nomen_patterns.keys()):
                nomen_field_desc = (
                    nomen_dict.get("NOMEN", "howard calculation NOMEN")
                    + f""". Format '{nomen_patterns.get(nomen_field)}'"""
                )

            # Create VCF header field
            vcf_reader.infos[nomen_field] = vcf.parser._Info(
                nomen_field,
                1,
                "String",
                nomen_field_desc,
                "howard calculation",
                "0",
                self.code_type_map.get("String"),
            )

            # Add field to SQL query update
            sql_nomen_fields.append(f"""
                    CASE 
                        WHEN {nomen_annotations_view}."{nomen_field}" NOT NULL AND {nomen_annotations_view}."{nomen_field}" NOT IN ('')
                        THEN concat(
                                ';{nomen_field}=',
                                {nomen_annotations_view}."{nomen_field}"
                            )
                        ELSE ''
                    END
                """)

        # SQL set for update
        sql_nomen_fields_set = ", ".join(sql_nomen_fields)

        # Update
        sql_update = f"""
            UPDATE {transcripts_table}
            SET "INFO" = 
                concat(
                    CASE
                        WHEN "INFO" IS NULL
                        THEN ''
                        ELSE concat("INFO", ';')
                    END,
                    regexp_replace(
                        concat(
                            {sql_nomen_fields_set}
                        )
                        ,'^;', ''
                    )
                )
            FROM {nomen_annotations_view}
            WHERE {transcripts_table}."#CHROM" = {nomen_annotations_view}."#CHROM"
                AND {transcripts_table}."POS" = {nomen_annotations_view}."POS" 
                AND {transcripts_table}."REF" = {nomen_annotations_view}."REF"
                AND {transcripts_table}."ALT" = {nomen_annotations_view}."ALT"
        """
        log.debug(f"Start NOMEN update...")
        self.conn.execute(sql_update)
        log.debug(f"Stop NOMEN update...")

        # Remove tables and view
        self.remove_tables_or_views(tables=[annotations_view, nomen_annotations_view])

    def calculation_find_samples(
        self, section: str = "calculation", tags: dict = None, **kwargs
    ) -> None:
        """
        The function `calculation_find_samples` performs a calculation to find samples in a VCF file and update the variant information in the database.
        It calculates the number of samples and the list of samples for each variant, and updates the corresponding fields in the variants table.

        :param section: The `section` parameter is a string that represents the section of the calculation to be performed, defaults to "calculation".
        :type section: str (optional)

        :param tags: The `tags` parameter is a dictionary that represents the annotation fields in the VCF file.
        It is used to create the annotation fields in the VCF header and to update the corresponding fields in the variants table, defaults to
        {"count_samples": "count", "list_samples": "list"}
        :type tags: dict (optional)

        :param kwargs: The `kwargs` parameter is a dictionary that allows you to pass additional

        """

        log.debug("Start find samples calculation...")

        # Get operation parameters
        operation_params, _ = self.get_operation_params(
            section=section, operation_params=kwargs, operation_name="find_samples"
        )

        ### Parameters for find samples calculation

        # find samples tags config
        if tags is None or tags == {}:
            tags = (
                operation_params.get("tags", None)
                or {"count_samples": "count", "list_samples": "list"}
            )

        # Process each tag in the find_samples_tags_config
        for tag in tags:

            # Method
            method = tags.get(tag, "count")

            # if FORMAT and samples
            if (
                "FORMAT" in self.get_header_columns_as_list()
                and self.get_header_sample_list()
            ):

                # findbypipeline annotation field
                findbypipeline_tag = tag

                # VCF infos tags
                vcf_infos_tags = {
                    findbypipeline_tag: f"Samples of a variant ({method})",
                }

                # Prefix
                prefix = self.get_explode_infos_prefix()

                # Field
                findbypipeline_infos = prefix + findbypipeline_tag

                # Variants table
                table_variants = self.get_table_variants()

                # Header
                vcf_reader = self.get_header()

                # Create variant id
                variant_id_column = self.get_variant_id_column()
                added_columns = [variant_id_column]

                # variant_id, FORMAT and samples
                samples_fields = f" {variant_id_column}, FORMAT , " + " , ".join(
                    [f""" "{sample}" """ for sample in self.get_header_sample_list()]
                )

                # Create dataframe
                dataframe_findbypipeline = self.get_query_to_df(
                    f""" SELECT {samples_fields} FROM {table_variants} """
                )

                # Create findbypipeline column
                dataframe_findbypipeline[findbypipeline_infos] = (
                    dataframe_findbypipeline.apply(
                        lambda row: findbypipeline(
                            row, samples=self.get_header_sample_list(), method=method
                        ),
                        axis=1,
                    )
                )

                # Add snpeff_hgvs to header
                vcf_reader.infos[findbypipeline_tag] = vcf.parser._Info(
                    findbypipeline_tag,
                    ".",
                    "String",
                    vcf_infos_tags.get(findbypipeline_tag, "Find in pipeline/sample"),
                    "howard calculation",
                    "0",
                    self.code_type_map.get("String"),
                )

                # Update
                sql_update = f"""
                    UPDATE variants
                    SET "INFO" = 
                        concat(
                            CASE
                                WHEN "INFO" IS NULL OR "INFO" IN ('','.')
                                THEN ''
                                ELSE concat("INFO", ';')
                            END,
                            CASE 
                                WHEN dataframe_findbypipeline."{findbypipeline_infos}" NOT IN ('','.')
                                    AND dataframe_findbypipeline."{findbypipeline_infos}" NOT NULL
                                THEN concat(
                                        '{findbypipeline_tag}=',
                                        dataframe_findbypipeline."{findbypipeline_infos}"
                                    )
                                ELSE ''
                            END
                        )
                    FROM dataframe_findbypipeline
                    WHERE variants."{variant_id_column}" = dataframe_findbypipeline."{variant_id_column}"
                """
                self.conn.execute(sql_update)

                # Remove added columns
                for added_column in added_columns:
                    self.drop_column(column=added_column)

                # Delete dataframe
                del dataframe_findbypipeline
                gc.collect()

    def calculation_genotype_concordance(self) -> None:
        """
        The function `calculation_genotype_concordance` calculates the genotype concordance for
        multi-caller VCF files and updates the variant information in the database.
        """

        # if FORMAT and samples
        if (
            "FORMAT" in self.get_header_columns_as_list()
            and self.get_header_sample_list()
        ):

            # genotypeconcordance annotation field
            genotypeconcordance_tag = "genotypeconcordance"

            # VCF infos tags
            vcf_infos_tags = {
                genotypeconcordance_tag: "Concordance of genotype for multi caller VCF",
            }

            # Prefix
            prefix = self.get_explode_infos_prefix()

            # Field
            genotypeconcordance_infos = prefix + genotypeconcordance_tag

            # Variants table
            table_variants = self.get_table_variants()

            # Header
            vcf_reader = self.get_header()

            # Create variant id
            variant_id_column = self.get_variant_id_column()
            added_columns = [variant_id_column]

            # variant_id, FORMAT and samples
            samples_fields = f" {variant_id_column}, FORMAT , " + " , ".join(
                [f""" "{sample}" """ for sample in self.get_header_sample_list()]
            )

            # Create dataframe
            dataframe_genotypeconcordance = self.get_query_to_df(
                f""" SELECT {samples_fields} FROM {table_variants} """
            )

            # Create genotypeconcordance column
            dataframe_genotypeconcordance[genotypeconcordance_infos] = (
                dataframe_genotypeconcordance.apply(
                    lambda row: genotypeconcordance(
                        row, samples=self.get_header_sample_list()
                    ),
                    axis=1,
                )
            )

            # Add genotypeconcordance to header
            vcf_reader.infos[genotypeconcordance_tag] = vcf.parser._Info(
                genotypeconcordance_tag,
                ".",
                "String",
                vcf_infos_tags.get(genotypeconcordance_tag, "snpEff hgvs annotations"),
                "howard calculation",
                "0",
                self.code_type_map.get("String"),
            )

            # Update
            sql_update = f"""
                UPDATE variants
                SET "INFO" = 
                    concat(
                        CASE
                            WHEN "INFO" IS NULL OR "INFO" IN ('','.')
                            THEN ''
                            ELSE concat("INFO", ';')
                        END,
                        CASE
                            WHEN dataframe_genotypeconcordance."{genotypeconcordance_infos}" NOT IN ('','.')
                                AND dataframe_genotypeconcordance."{genotypeconcordance_infos}" NOT NULL
                            THEN concat(
                                    '{genotypeconcordance_tag}=',
                                    dataframe_genotypeconcordance."{genotypeconcordance_infos}"
                                )
                            ELSE ''
                        END
                    )
                FROM dataframe_genotypeconcordance
                WHERE variants."{variant_id_column}" = dataframe_genotypeconcordance."{variant_id_column}"
            """
            self.conn.execute(sql_update)

            # Remove added columns
            for added_column in added_columns:
                self.drop_column(column=added_column)

            # Delete dataframe
            del dataframe_genotypeconcordance
            gc.collect()

    def calculation_barcode(self, tag: str = "barcode") -> None:
        """
        The `calculation_barcode` function calculates barcode values for variants in a VCF file and
        updates the INFO field in the file with the calculated barcode values.

        :param tag: The `tag` parameter in the `calculation_barcode` function is used to specify the tag
        name that will be used for the barcode calculation in the VCF file. If no tag name is provided,
        the default tag name is set to "barcode", defaults to barcode
        :type tag: str (optional)
        """

        # if FORMAT and samples
        if (
            "FORMAT" in self.get_header_columns_as_list()
            and self.get_header_sample_list()
        ):

            # barcode annotation field
            if not tag:
                tag = "barcode"

            # VCF infos tags
            vcf_infos_tags = {
                tag: "barcode calculation (VaRank)",
            }

            # Prefix
            prefix = self.get_explode_infos_prefix()

            # Field
            barcode_infos = prefix + tag

            # Variants table
            table_variants = self.get_table_variants()

            # Header
            vcf_reader = self.get_header()

            # Create variant id
            variant_id_column = self.get_variant_id_column()
            added_columns = [variant_id_column]

            # variant_id, FORMAT and samples
            samples_fields = f" {variant_id_column}, FORMAT , " + " , ".join(
                [f""" "{sample}" """ for sample in self.get_header_sample_list()]
            )

            # Create dataframe
            dataframe_barcode = self.get_query_to_df(
                f""" SELECT {samples_fields} FROM {table_variants} """
            )

            # Create barcode column
            dataframe_barcode[barcode_infos] = dataframe_barcode.apply(
                lambda row: barcode(row, samples=self.get_header_sample_list()), axis=1
            )

            # Add barcode to header
            vcf_reader.infos[tag] = vcf.parser._Info(
                tag,
                "1",
                "String",
                vcf_infos_tags.get(tag, vcf_infos_tags.get(tag)),
                "howard calculation",
                "0",
                self.code_type_map.get("String"),
            )

            # Update
            sql_update = f"""
                UPDATE {table_variants}
                SET "INFO" = 
                    concat(
                        CASE
                            WHEN "INFO" IS NULL OR "INFO" IN ('','.')
                            THEN ''
                            ELSE concat("INFO", ';')
                        END,
                        CASE
                            WHEN dataframe_barcode."{barcode_infos}" NOT IN ('','.')
                            AND dataframe_barcode."{barcode_infos}" NOT NULL
                            THEN concat(
                                    '{tag}=',
                                    dataframe_barcode."{barcode_infos}"
                                )
                            ELSE ''
                        END
                    )
                FROM dataframe_barcode
                WHERE {table_variants}."{variant_id_column}" = dataframe_barcode."{variant_id_column}"
            """
            self.conn.execute(sql_update)

            # Remove added columns
            for added_column in added_columns:
                self.drop_column(column=added_column)

            # Delete dataframe
            del dataframe_barcode
            gc.collect()

    def calculation_barcode_family(
        self, section: str = "calculation", tag: str = None, tag_samples: str = None
    ) -> None:
        """
        The `calculation_barcode_family` function calculates barcode values for variants in a VCF file
        and updates the INFO field in the file with the calculated barcode values.

        :param tag: The `tag` parameter in the `calculation_barcode_family` function is used to specify
        the barcode tag that will be added to the VCF file during the calculation process. If no value
        is provided for the `tag` parameter, the default value used is "BCF", defaults to BCF
        :type tag: str (optional)
        :param tag_samples: The `tag_samples` parameter in the `calculation_barcode_family` function is
        used to specify the barcode tag that will be added to the VCF file for samples during the
        calculation process. If no value is provided for the `tag_samples` parameter, the default value
        used is "BCFS", defaults to BCFS

        """

        # if FORMAT and samples
        if (
            "FORMAT" in self.get_header_columns_as_list()
            and self.get_header_sample_list()
        ):

            # barcode annotation field
            if not tag:
                tag = "BCF"

            # barcode annotation field for samples
            if not tag_samples:
                tag_samples = f"{tag}S"

            # VCF infos tags
            vcf_infos_tags = {
                "tag": "barcode family calculation",
                "tag_samples": "barcode family samples",
            }

            # Param
            param = self.get_param()
            log.debug(f"param={param}")

            # Prefix
            prefix = self.get_explode_infos_prefix()

            # PED param
            ped = (
                param.get(section, {})
                .get("calculations", {})
                .get("BARCODEFAMILY", {})
                .get("family_pedigree", None)
            )
            log.debug(f"ped={ped}")

            # Load PED
            if ped:

                # Pedigree is a file
                if isinstance(ped, str) and os.path.exists(full_path(ped)):
                    log.debug("Pedigree is file")
                    with open(full_path(ped)) as ped:
                        ped = yaml.safe_load(ped)

                # Pedigree is a string
                elif isinstance(ped, str):
                    log.debug("Pedigree is str")
                    try:
                        ped = json.loads(ped)
                        log.debug("Pedigree is json str")
                    except ValueError as e:
                        ped_samples = ped.split(",")
                        ped = {}
                        for ped_sample in ped_samples:
                            ped[ped_sample] = ped_sample

                # Pedigree is a dict
                elif isinstance(ped, dict):
                    log.debug("Pedigree is dict")

                # Pedigree is not well formatted
                else:
                    msg_error = "Pedigree not well formatted"
                    log.error(msg_error)
                    raise ValueError(msg_error)

                # Construct list
                ped_samples = list(ped.values())

            else:
                log.debug("Pedigree not defined. Take all samples")
                ped_samples = self.get_header_sample_list()
                ped = {}
                for ped_sample in ped_samples:
                    ped[ped_sample] = ped_sample

            # Check pedigree
            if not ped or len(ped) == 0:
                msg_error = f"Error in pedigree: samples {ped_samples}"
                log.error(msg_error)
                raise ValueError(msg_error)

            # Log
            log.info(
                "Calculation 'BARCODEFAMILY' - Samples: "
                + ", ".join([f"{member}='{ped[member]}'" for member in ped])
            )
            log.debug(f"ped_samples={ped_samples}")

            # Header
            vcf_reader = self.get_header()

            # Check for other tag names starting with 'tag'
            log.debug(f"tag={tag}")
            log.debug(f"tag_samples={tag_samples}")
            if tag in vcf_reader.formats:
                # Create a new tag name with a suffix based on the number of tags match with the same 'tag' pattern '{tag}_<integer>'
                tag_new = f"{tag}_" + str(
                    len(
                        [
                            t
                            for t in vcf_reader.formats
                            if (t == tag or re.match(rf"^{tag}_\d+$", t))
                        ]
                    )
                )

                tag = tag_new
            if tag_samples in vcf_reader.formats:
                # Create a new tag name with a suffix based on the number of tags match with the same 'tag' pattern '{tag}_<integer>'
                tag_samples_new = f"{tag_samples}_" + str(
                    len(
                        [
                            t
                            for t in vcf_reader.formats
                            if (
                                t == tag_samples or re.match(rf"^{tag_samples}_\d+$", t)
                            )
                        ]
                    )
                )

                tag_samples = tag_samples_new

            # Create vcf_infos_tags for the tags
            vcf_infos_tags[tag] = vcf_infos_tags.get(
                "tag", "barcode family calculation"
            )
            vcf_infos_tags[tag_samples] = vcf_infos_tags.get(
                "tag_samples", "barcode family samples"
            )
            log.debug(f"vcf_infos_tags={vcf_infos_tags}")

            # Field
            barcode_infos = prefix + tag

            # Variants table
            table_variants = self.get_table_variants()

            # Create variant id
            variant_id_column = self.get_variant_id_column()
            added_columns = [variant_id_column]

            # variant_id, FORMAT and samples
            samples_fields = f" {variant_id_column}, FORMAT , " + " , ".join(
                [f""" "{sample}" """ for sample in ped_samples]
            )

            # Create dataframe
            dataframe_barcode = self.get_query_to_df(
                f""" SELECT {samples_fields} FROM {table_variants} """
            )

            # Create barcode column
            dataframe_barcode[barcode_infos] = dataframe_barcode.apply(
                lambda row: barcode(row, samples=ped_samples), axis=1
            )

            # Add barcode family to header
            # Add vaf_normalization to header
            vcf_reader.formats[tag] = vcf.parser._Format(
                id=tag,
                num=1,
                type="String",
                desc=vcf_infos_tags.get(
                    tag, f"barcode family calculation for {ped_samples}"
                )
                + f" for {ped_samples}",
                type_code=self.code_type_map.get("String"),
            )
            vcf_reader.formats[f"{tag}S"] = vcf.parser._Format(
                id=tag_samples,
                num=str(len(ped_samples)),
                type="String",
                desc=vcf_infos_tags.get(
                    tag_samples, f"barcode family samples for {ped_samples}"
                )
                + f" for {ped_samples}",
                type_code=self.code_type_map.get("String"),
            )

            # Update
            # for sample in ped_samples:
            sql_update_set = []
            for sample in self.get_header_sample_list() + ["FORMAT"]:
                if sample in ped_samples:
                    value = f'dataframe_barcode."{barcode_infos}"'
                    value_samples = (
                        "'"
                        + ",".join([f"""{sample}""" for sample in ped_samples])
                        + "'"
                    )
                    ped_samples
                elif sample == "FORMAT":
                    value = f"'{tag}'"
                    value_samples = f"'{tag_samples}'"
                else:
                    value = "'.'"
                    value_samples = "'.'"

                # Format regex
                format_regex = r"[a-zA-Z0-9\s]"

                # Update query
                sql_update_set.append(
                    f"""
                        "{sample}" = 
                        concat(
                            CASE
                                WHEN {table_variants}."{sample}" = './.'
                                THEN concat('./.',regexp_replace(regexp_replace({table_variants}.FORMAT, '{format_regex}', '', 'g'), ':', ':.', 'g'))
                                ELSE {table_variants}."{sample}"
                            END,
                            ':',
                            {value},
                            ':',
                            {value_samples}
                        )
                    """
                )

            sql_update_set_join = ", ".join(sql_update_set)
            sql_update = f"""
                UPDATE {table_variants}
                SET {sql_update_set_join}
                FROM dataframe_barcode
                WHERE {table_variants}."{variant_id_column}" = dataframe_barcode."{variant_id_column}"
            """
            self.conn.execute(sql_update)

            # Remove added columns
            for added_column in added_columns:
                self.drop_column(column=added_column)

            # Delete dataframe
            del dataframe_barcode
            gc.collect()

    def calculation_trio(self, section: str = "calculation") -> None:
        """
        The `calculation_trio` function performs trio calculations on a VCF file by adding trio
        information to the INFO field of each variant.
        """

        # if FORMAT and samples
        if (
            "FORMAT" in self.get_header_columns_as_list()
            and self.get_header_sample_list()
        ):

            # trio annotation field
            trio_tag = "trio"

            # VCF infos tags
            vcf_infos_tags = {
                "trio": "trio calculation",
            }

            # Param
            param = self.get_param()

            # Prefix
            prefix = self.get_explode_infos_prefix()

            # Trio param
            trio_ped = (
                param.get(section, {})
                .get("calculations", {})
                .get("TRIO", {})
                .get("trio_pedigree", None)
            )

            # Load trio
            if trio_ped:

                # Trio pedigree is a file
                if isinstance(trio_ped, str) and os.path.exists(full_path(trio_ped)):
                    log.debug("TRIO pedigree is file")
                    with open(full_path(trio_ped)) as trio_ped:
                        trio_ped = yaml.safe_load(trio_ped)

                # Trio pedigree is a string
                elif isinstance(trio_ped, str):
                    log.debug("TRIO pedigree is str")
                    try:
                        trio_ped = json.loads(trio_ped)
                        log.debug("TRIO pedigree is json str")
                    except ValueError as e:
                        trio_samples = trio_ped.split(",")
                        if len(trio_samples) == 3:
                            trio_ped = {
                                "father": trio_samples[0],
                                "mother": trio_samples[1],
                                "child": trio_samples[2],
                            }
                            log.debug("TRIO pedigree is list str")
                        else:
                            msg_error = "TRIO pedigree not well formatted"
                            log.error(msg_error)
                            raise ValueError(msg_error)

                # Trio pedigree is a dict
                elif isinstance(trio_ped, dict):
                    log.debug("TRIO pedigree is dict")

                # Trio pedigree is not well formatted
                else:
                    msg_error = "TRIO pedigree not well formatted"
                    log.error(msg_error)
                    raise ValueError(msg_error)

                # Construct trio list
                trio_samples = [
                    trio_ped.get("father", ""),
                    trio_ped.get("mother", ""),
                    trio_ped.get("child", ""),
                ]

            else:
                log.debug("TRIO pedigree not defined. Take the first 3 samples")
                samples_list = self.get_header_sample_list()
                if len(samples_list) >= 3:
                    trio_samples = self.get_header_sample_list()[0:3]
                    trio_ped = {
                        "father": trio_samples[0],
                        "mother": trio_samples[1],
                        "child": trio_samples[2],
                    }
                else:
                    msg_error = f"Error in TRIO pedigree: only {len(samples_list)} samples {samples_list}"
                    log.error(msg_error)
                    raise ValueError(msg_error)

            # Check trio pedigree
            if not trio_ped or len(trio_ped) != 3:
                msg_error = f"Error in TRIO pedigree: {trio_ped}"
                log.error(msg_error)
                raise ValueError(msg_error)

            # Log
            log.info(
                f"Calculation 'TRIO' - Samples: "
                + ", ".join([f"{member}='{trio_ped[member]}'" for member in trio_ped])
            )

            # Field
            trio_infos = prefix + trio_tag

            # Variants table
            table_variants = self.get_table_variants()

            # Header
            vcf_reader = self.get_header()

            # Create variant id
            variant_id_column = self.get_variant_id_column()
            added_columns = [variant_id_column]

            # variant_id, FORMAT and samples
            samples_fields = f" {variant_id_column}, FORMAT , " + " , ".join(
                [f""" "{sample}" """ for sample in self.get_header_sample_list()]
            )

            # Create dataframe
            dataframe_trio = self.get_query_to_df(
                f""" SELECT {samples_fields} FROM {table_variants} """
            )

            # Create trio column
            dataframe_trio[trio_infos] = dataframe_trio.apply(
                lambda row: trio(row, samples=trio_samples), axis=1
            )

            # Add trio to header
            vcf_reader.infos[trio_tag] = vcf.parser._Info(
                trio_tag,
                ".",
                "String",
                vcf_infos_tags.get(trio_tag, "snpEff hgvs annotations"),
                "howard calculation",
                "0",
                self.code_type_map.get("String"),
            )

            # Update
            sql_update = f"""
                UPDATE {table_variants}
                SET "INFO" = 
                    concat(
                        CASE
                            WHEN "INFO" IS NULL OR "INFO" IN ('','.')
                            THEN ''
                            ELSE concat("INFO", ';')
                        END,
                        CASE
                            WHEN dataframe_trio."{trio_infos}" NOT IN ('','.')
                             AND dataframe_trio."{trio_infos}" NOT NULL
                            THEN concat(
                                    '{trio_tag}=',
                                    dataframe_trio."{trio_infos}"
                                )
                            ELSE ''
                        END
                    )
                FROM dataframe_trio
                WHERE {table_variants}."{variant_id_column}" = dataframe_trio."{variant_id_column}"
            """
            self.conn.execute(sql_update)

            # Remove added columns
            for added_column in added_columns:
                self.drop_column(column=added_column)

            # Delete dataframe
            del dataframe_trio
            gc.collect()

    def calculation_vaf_normalization(self) -> None:
        """
        The `calculation_vaf_normalization` function calculates the VAF (Variant Allele Frequency)
        normalization for each sample in a VCF file and updates the FORMAT and INFO fields accordingly.
        :return: The function does not return anything.
        """

        # if FORMAT and samples
        if (
            "FORMAT" in self.get_header_columns_as_list()
            and self.get_header_sample_list()
        ):

            # vaf_normalization annotation field
            vaf_normalization_tag = "VAF"

            # VCF infos tags
            vcf_infos_tags = {
                "VAF": "VAF Variant Frequency",
            }

            # Prefix
            prefix = self.get_explode_infos_prefix()

            # Variants table
            table_variants = self.get_table_variants()

            # Header
            vcf_reader = self.get_header()

            # Do not calculate if VAF already exists
            if "VAF" in vcf_reader.formats:
                log.debug("VAF already on genotypes")
                return

            # Create variant id
            variant_id_column = self.get_variant_id_column()
            added_columns = [variant_id_column]

            # variant_id, FORMAT and samples
            samples_fields = f" {variant_id_column}, FORMAT , " + " , ".join(
                f""" "{sample}" """ for sample in self.get_header_sample_list()
            )

            # Create dataframe
            query = f""" SELECT {variant_id_column}, FORMAT, {samples_fields} FROM {table_variants} """
            log.debug(f"query={query}")
            dataframe_vaf_normalization = self.get_query_to_df(query=query)

            vaf_normalization_set = []

            # for each sample vaf_normalization
            for sample in self.get_header_sample_list():
                dataframe_vaf_normalization[sample] = dataframe_vaf_normalization.apply(
                    lambda row: vaf_normalization(row, sample=sample), axis=1
                )
                vaf_normalization_set.append(
                    f""" "{sample}" = dataframe_vaf_normalization."{sample}" """
                )

            # Add VAF to FORMAT
            dataframe_vaf_normalization["FORMAT"] = dataframe_vaf_normalization[
                "FORMAT"
            ].apply(lambda x: str(x) + ":VAF")
            vaf_normalization_set.append(
                f""" "FORMAT" = dataframe_vaf_normalization."FORMAT" """
            )

            # Add vaf_normalization to header
            vcf_reader.formats[vaf_normalization_tag] = vcf.parser._Format(
                id=vaf_normalization_tag,
                num=1,
                type="Float",
                desc=vcf_infos_tags.get(vaf_normalization_tag, "VAF Variant Frequency"),
                type_code=self.code_type_map.get("Float"),
            )

            # Create fields to add in INFO
            sql_vaf_normalization_set = " , ".join(vaf_normalization_set)

            # Update
            sql_update = f"""
                UPDATE {table_variants}
                SET {sql_vaf_normalization_set}
                FROM dataframe_vaf_normalization
                WHERE variants."{variant_id_column}" = dataframe_vaf_normalization."{variant_id_column}"

            """
            self.conn.execute(sql_update)

            # Remove added columns
            for added_column in added_columns:
                self.drop_column(column=added_column)

            # Delete dataframe
            del dataframe_vaf_normalization
            gc.collect()


    def calculation_info_to_format(self, section: str = "calculation", annotation_fields: dict = None, samples:list = None, remove_info_fields: bool = None, **kwargs) -> None:
        """
        The `calculation_info_to_format` function converts INFO fields to FORMAT fields in a VCF file.

        :param section: The `section` parameter is a string that specifies the section of the configuration file to use for the calculation. It is used to retrieve the relevant parameters for the operation, defaults to "calculation"
        :type section: str (optional)
        :param annotation_fields: The `annotation_fields` parameter is a dictionary that maps INFO field names to FORMAT field names. It specifies which INFO fields should be converted to FORMAT fields and what the corresponding FORMAT field names should be, defaults to None
        :type annotation_fields: dict (optional)
        :param samples: The `samples` parameter is a list of sample names for which the INFO fields should be converted to FORMAT fields. If not provided, all samples in the VCF file will be used, defaults to None
        :type samples: list (optional)
        :param remove_info_fields: The `remove_info_fields` parameter is a boolean that determines whether the original INFO fields should be removed from the VCF file after they have been converted to FORMAT fields. If set to True, the INFO fields will be removed; if set to False, they will be retained, defaults to None
        :type remove_info_fields: bool (optional)

        :return: The function does not return anything.

        Example kwargs:

        {
            "annotation_fields": {
                "count_samples": None,
                "list_samples": "LS",
                "calling_quality": "CQ"
            },
            "samples": ["sample1", "sample2"],
            "remove_info_fields": True
        }
    

        """

        # if FORMAT and samples
        if (
            "FORMAT" in self.get_header_columns_as_list()
            and self.get_header_sample_list()
        ):

            log.debug("Calculation INFO to FORMAT...")

            operation_params, _ = self.get_operation_params(
                section=section, operation_params=kwargs, operation_name="info_to_format"
            )

            ### Parameters for variant ID calculation

            # variant_id annotation field
            if annotation_fields is None:
                annotation_fields = (
                    operation_params.get("annotation_fields")
                    or annotation_fields
                    or {}
                )

            # variant_id_tag_info
            if samples is None:
                samples = (
                    operation_params.get("samples")
                    or samples
                    or self.get_header_sample_list()
                )

            # keep_variant_id_tag_column
            if remove_info_fields is None: 
                remove_info_fields = (
                    operation_params.get("remove_info_fields")
                    or remove_info_fields
                    or False
                )

            # Variants table
            table_variants = self.get_table_variants()

            # Check info and format fiel to determine which are in info and not in format
            header = self.get_header()
            header_infos = header.infos
            header_formats = header.formats
            annotation_fields_filtered = {}
            for info_field, format_field in annotation_fields.items():
                if format_field in header_formats:
                    msg_error = f"FORMAT field '{format_field}' already in VCF header format"
                    log.error(msg_error)
                    raise ValueError(msg_error)
                elif info_field not in header_infos:
                    msg_error = f"INFO field '{info_field}' not found in VCF header"
                    log.error(msg_error)
                    raise ValueError(msg_error)
                else:
                    format_field = format_field or info_field
                    annotation_fields_filtered[info_field] = format_field
                    # Add format field to header
                    header.formats[format_field] = vcf.parser._Format(
                        id=format_field,
                        num=header_infos[info_field].num,
                        type=header_infos[info_field].type,
                        desc=header_infos[info_field].desc,
                        type_code=header_infos[info_field].type_code,
                    )

            # Create view from table_variants
            annotation_view_name = f"variants_view_{str(get_random())}"
            annotation_view_name = self.create_annotations_view(
                table=table_variants,
                view=annotation_view_name,
                view_type="view",
                view_mode="explore",
                info_prefix_column="",
                detect_type_list=False,
                fields=annotation_fields_filtered.keys(),
                fields_not_exists=True,
                fields_forced_as_varchar=True,
                fields_needed_all=False,
                fields_to_rename=None,
                drop_view=True,
            )

            # Prepare set clause for samples
            set_clause = []
            for sample in self.get_header_sample_list():

                # For each sample, prepare the set clause to add the INFO fields values to the FORMAT fields values
                sample_set = []
                for info_field, format_field in annotation_fields_filtered.items():
                    if header_infos.get(info_field).type != "Flag":
                        if sample in samples:
                            sample_set.append(
                                f"""
                                    CASE
                                        WHEN src."{info_field}" IS NOT NULL OR src."{info_field}"::STRING NOT IN ('','.')
                                        THEN replace(src."{info_field}"::STRING, ':', '=')
                                        ELSE '.'
                                    END
                                """
                                )
                        else:
                            sample_set.append(" '.' ")
                    else:
                        raise NotImplementedError(f"Flag type for INFO field '{info_field}' not implemented yet for INFO to FORMAT calculation")
                
                # If there are INFO fields to add to FORMAT for this sample, prepare the set clause to add them to the sample FORMAT field
                if len(sample_set):
                    set_clause.append(f""" "{sample}" = concat("{sample}", ':', {", ':', ".join(sample_set)}) """)

            # Prepare set clause for FORMAT
            if len(annotation_fields_filtered):
                set_clause.append(f""" "FORMAT" = concat("FORMAT", ':', {", ':', ".join([f"'{format_field}'" if format_field is not None else f"'{info_field}'" for info_field, format_field in annotation_fields_filtered.items() ])}) """)

            # Update
            if len(set_clause):
                sql_update = f"""
                    UPDATE {table_variants} v
                    SET {', '.join(set_clause)}
                    FROM {annotation_view_name} src
                    WHERE v."#CHROM" = src."#CHROM"
                    AND v."POS"    = src."POS"
                    AND v."REF"    = src."REF"
                    AND v."ALT"    = src."ALT"

                """
                # log.debug(f"sql_update={sql_update}")
                self.get_connexion().execute(sql_update)

                if remove_info_fields:
                    self.rename_info_fields(
                        fields_to_rename=dict.fromkeys(annotation_fields_filtered.keys(), None), table=table_variants
                    )


    def calculation_genotype_stats(self, info: str = "VAF", **kwargs) -> None:
        """
        The `calculation_genotype_stats` function calculates genotype statistics for a given information
        field in a VCF file and updates the INFO column of the variants table with the calculated
        statistics.

        :param info: The `info` parameter is a string that represents the type of information for which
        genotype statistics are calculated. It is used to generate various VCF info tags for the
        statistics, such as the number of occurrences, the list of values, the minimum value, the
        maximum value, the mean, the median, defaults to VAF
        :type info: str (optional)
        """

        # if FORMAT and samples
        if (
            "FORMAT" in self.get_header_columns_as_list()
            and self.get_header_sample_list()
        ):

            # vaf_stats annotation field
            vaf_stats_tag = info + "_stats"

            # VCF infos tags
            vcf_infos_tags = {
                info + "_stats_nb": f"genotype {info} Statistics - number of {info}",
                info + "_stats_list": f"genotype {info} Statistics - list of {info}",
                info + "_stats_min": f"genotype {info} Statistics - min {info}",
                info + "_stats_max": f"genotype {info} Statistics - max {info}",
                info + "_stats_mean": f"genotype {info} Statistics - mean {info}",
                info + "_stats_mediane": f"genotype {info} Statistics - mediane {info}",
                info
                + "_stats_stdev": f"genotype {info} Statistics - standard deviation {info}",
            }

            # Prefix
            prefix = self.get_explode_infos_prefix()

            # Field
            vaf_stats_infos = prefix + vaf_stats_tag

            # Variants table
            table_variants = self.get_table_variants()

            # Header
            vcf_reader = self.get_header()

            # Create variant id
            variant_id_column = self.get_variant_id_column()
            added_columns = [variant_id_column]

            # variant_id, FORMAT and samples
            samples_fields = f" {variant_id_column}, FORMAT , " + " , ".join(
                [f""" "{sample}" """ for sample in self.get_header_sample_list()]
            )

            # Create dataframe
            dataframe_vaf_stats = self.get_query_to_df(
                f""" SELECT {samples_fields} FROM {table_variants} """
            )

            # Create vaf_stats column
            dataframe_vaf_stats[vaf_stats_infos] = dataframe_vaf_stats.apply(
                lambda row: genotype_stats(
                    row, samples=self.get_header_sample_list(), info=info
                ),
                axis=1,
            )

            # List of vcf tags
            sql_vaf_stats_fields = []

            # Check all VAF stats infos
            for stat in vcf_infos_tags:

                # Extract stats
                dataframe_vaf_stats[stat] = dataframe_vaf_stats[vaf_stats_infos].apply(
                    lambda x: dict(x).get(stat, "")
                )

                # Add snpeff_hgvs to header
                vcf_reader.infos[stat] = vcf.parser._Info(
                    stat,
                    ".",
                    "String",
                    vcf_infos_tags.get(stat, "genotype statistics"),
                    "howard calculation",
                    "0",
                    self.code_type_map.get("String"),
                )

                if len(sql_vaf_stats_fields):
                    sep = ";"
                else:
                    sep = ""

                # Create fields to add in INFO
                sql_vaf_stats_fields.append(
                    f"""
                        CASE
                            WHEN dataframe_vaf_stats."{stat}" NOT NULL
                            THEN concat(
                                    '{sep}{stat}=',
                                    dataframe_vaf_stats."{stat}"
                                )
                            ELSE ''
                        END
                    """
                )

            # SQL set for update
            sql_vaf_stats_fields_set = ",  ".join(sql_vaf_stats_fields)

            # Update
            sql_update = f"""
                UPDATE {table_variants}
                SET "INFO" = 
                    concat(
                        CASE
                            WHEN "INFO" IS NULL OR "INFO" IN ('','.')
                            THEN ''
                            ELSE concat("INFO", ';')
                        END,
                        {sql_vaf_stats_fields_set}
                    )
                FROM dataframe_vaf_stats
                WHERE {table_variants}."{variant_id_column}" = dataframe_vaf_stats."{variant_id_column}"

            """
            self.conn.execute(sql_update)

            # Remove added columns
            for added_column in added_columns:
                self.drop_column(column=added_column)

            # Delete dataframe
            del dataframe_vaf_stats
            gc.collect()

    def calculation_transcripts_annotation(
        self, section: str = "calculation", info_json: str = None, info_format: str = None, **kwargs
    ) -> None:
        """
        The `calculation_transcripts_annotation` function creates a transcripts table and adds an info
        field to it if transcripts are available.

        :param section: The `section` parameter in the `calculation_transcripts_annotation` method
        is a string parameter that represents the section of the calculation to be used. It is used to
        specify the section for the transcripts annotation. If no value is provided when calling the
        method, it defaults to "calculation".
        :type section: str
        :param info_json: The `info_json` parameter in the `calculation_transcripts_annotation` method
        is a string parameter that represents the information field to be used in the transcripts JSON.
        It is used to specify the JSON format for the transcripts information. If no value is provided
        when calling the method, it defaults to "
        :type info_json: str
        :param info_format: The `info_format` parameter in the `calculation_transcripts_annotation`
        method is a string parameter that specifies the format of the information field to be used in
        the transcripts JSON. It is used to define the format of the information field
        :type info_format: str
        """

        # Param from calculation
        param_calculation = self.get_param().get(section, {}).copy()
        param_transcripts_annotation = {
            k.upper(): v for k, v in param_calculation.get("calculations", {}).items()
        }.get("TRANSCRIPTS_ANNOTATIONS", None)

        if param_transcripts_annotation:
            param_calculation_transcripts_annotation = {
                "transcripts": param_transcripts_annotation
            }
        else:
            param_calculation_transcripts_annotation = None

        # Create transcripts table
        transcripts_table = self.create_transcript_view(
            param=param_calculation_transcripts_annotation
        )

        # Add info field
        if transcripts_table:
            self.transcript_view_to_variants(
                transcripts_table=transcripts_table,
                transcripts_info_field_json=info_json,
                transcripts_info_field_format=info_format,
            )
        else:
            log.info("No Transcripts to process. Check param.json file configuration")

    def calculation_transcripts_prioritization(self, section: str = "calculation", strict: bool = False, **kwargs) -> None:
        """
        The function `calculation_transcripts_prioritization` creates a transcripts table and
        prioritizes transcripts based on certain criteria.
        """

        # Param from calculation
        param_calculation = self.get_param().get(section, {}).copy()
        param_transcripts_prioritization = {
            k.upper(): v for k, v in param_calculation.get("calculations", {}).items()
        }.get("TRANSCRIPTS_PRIORITIZATION", None)

        if param_transcripts_prioritization:
            param_calculation_transcripts_prioritization = {
                "transcripts": param_transcripts_prioritization
            }
        else:
            param_calculation_transcripts_prioritization = None

        # Create transcripts table
        transcripts_table = self.create_transcript_view(
            param=param_calculation_transcripts_prioritization
        )

        # Add info field
        if transcripts_table:
            self.transcripts_prioritization(
                transcripts_table=transcripts_table,
                strict=strict,
                param=param_calculation_transcripts_prioritization,
            )
        else:
            log.info("No Transcripts to process. Check param.json file configuration")

    def calculation_transcripts_export(self, section: str = "calculation", **kwargs) -> None:
        """ """

        # Param from calculation
        param_calculation = self.get_param().get(section, {}).copy()
        param_transcripts_export = {
            k.upper(): v for k, v in param_calculation.get("calculations", {}).items()
        }.get("TRANSCRIPTS_EXPORT", None)

        if param_transcripts_export:
            param_calculation_transcripts_export = {
                "transcripts": param_transcripts_export
            }
        else:
            param_calculation_transcripts_export = None

        if (
            param_calculation_transcripts_export
            and param_calculation_transcripts_export.get("transcripts", None)
        ):
            param_calculation_transcripts_export_export = (
                param_calculation_transcripts_export.get("transcripts", {}).get(
                    "export", None
                )
            )
            param_calculation_transcripts_export_explode = (
                param_calculation_transcripts_export.get("transcripts", {}).get(
                    "explode", None
                )
            )
        else:
            param_calculation_transcripts_export_export = None
            param_calculation_transcripts_export_explode = None

        # Create transcripts table
        transcripts_table = self.create_transcript_view(
            param=param_calculation_transcripts_export
        )

        # Add info field
        if transcripts_table:
            self.transcripts_export(
                transcripts_table=transcripts_table,
                param_export=param_calculation_transcripts_export_export,
                param_explode=param_calculation_transcripts_export_explode,
            )
        else:
            log.info("No Transcripts to process. Check param.json file configuration")


    def calculation_variant_filter(
        self,
        section="calculation",
        where_clause: str = None,
        **kwargs
    ) -> None:
        """
        Filter variants based on specified criteria (using SQL parameters)
        """

        log.debug("Calculation variant_filter...")

        operation_params, _ = self.get_operation_params(
            section=section, operation_params=kwargs, operation_name="variant_filter"
        )

        ### Parameters for variant filter calculation

        # variant filter where_clause
        where_clause = (
            operation_params.get("where_clause")
            or where_clause
            or None
        )

        # Check where_clause
        if where_clause is None:
            log.warning(
                f"variant_filter calculation: No where clause specified in parameters. No filtering will be applied."
            )
            return None
        else:
            log.debug(f"variant_filter calculation: Applying where_clause: {where_clause}")

            table_variants = self.get_table_variants()

            # Create annotation view
            annotation_view_name = f"variants_view_{str(get_random())}"
            annotation_view_name = self.create_annotations_view(
                table=table_variants,
                view=annotation_view_name,
                view_type="view",
                view_mode="explore",
                info_prefix_column="",
                fields=None,
                fields_needed_all=True,
                info_struct_column="INFOS",
                sample_struct_column="SAMPLES",
                detect_type_list=True,
            )
            
            # Replace table variants by a SQL query
            query_replace_variants = f"""
                CREATE OR REPLACE TABLE {table_variants} AS
                SELECT variants_old.*
                FROM {table_variants} as variants_old
                JOIN (SELECT "#CHROM", "POS", "REF", "ALT" FROM {annotation_view_name} WHERE {where_clause}) as variants
                ON variants_old."#CHROM" = variants."#CHROM"
                AND variants_old.POS = variants.POS
                AND variants_old.REF = variants.REF
                AND variants_old.ALT = variants.ALT
                
            """
            #log.debug(f"variant_filter calculation: SQL query to replace variants table: {query_replace_variants}")
            self.execute_query(query_replace_variants)

        return None
