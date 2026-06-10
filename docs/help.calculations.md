# Calculations

List of available calculations

## BARCODE

BARCODE as VaRank tool

BARCODE as VaRank tool. This calculation generates a BARCODE for each sample based on sample genotype (0 for reference, 1 for heterozygous, 2 for homozygous alternate).

## BARCODEFAMILY

BARCODEFAMILY as VaRank tool

BARCODEFAMILY as VaRank tool. This calculation generates a BARCODE for each family based on sample genotype (0 for reference, 1 for heterozygous, 2 for homozygous alternate) and family relationships.

## DP_STATS

Depth (DP) statistics

Depth (DP) statistics. This calculation provides statistical measures for the sequencing depth (DP) across different samples. It helps in understanding the coverage and reliability of variant calls.

## FINDBYPIPELINE

Number of pipeline that identify the variant (for multi pipeline VCF)

Number of pipeline that identify the variant (for multi pipeline VCF). This calculation counts the number of pipelines that identify a given variant in a multi-pipeline VCF file. It helps in assessing the confidence of variant calls by determining how many different pipelines support the presence of the variant.

## FINDBYSAMPLE

Number of sample that have a genotype for the variant (for multi sample VCF)

Number of sample that have a genotype for the variant (for multi sample VCF). This calculation counts the number of samples that have a genotype for a given variant in a multi-sample VCF file. It helps in assessing the presence of the variant across different samples.

## GENOTYPECONCORDANCE

Concordance of genotype for multi caller VCF

Concordance of genotype for multi caller VCF. This calculation assesses the concordance of genotypes for a given variant across multiple callers in a multi-caller VCF file. It helps in evaluating the consistency of genotype calls and can be used to identify variants with high confidence based on agreement among different callers.

## LISTBYPIPELINE

List of pipelines that identify the variant (for multi pipeline VCF)

List of pipelines that identify the variant (for multi pipeline VCF). This calculation provides a list of pipelines that identify a given variant in a multi-pipeline VCF file. It helps in understanding which pipelines support the presence of the variant.

## LISTBYSAMPLE

List of samples that have a genotype for the variant (for multi sample VCF)

List of samples that have a genotype for the variant (for multi sample VCF). This calculation provides a list of samples that have a genotype for a given variant in a multi-sample VCF file. It helps in understanding which samples support the presence of the variant.

## MERGED_HGVS

Merge HGVS nomenclatures from snpEff (snpeff_hgvs) and ANNOVAR (AAChange_refGene) into merged_hgvs field

This calculation merges HGVS nomenclatures from snpEff (snpeff_hgvs) and ANNOVAR (AAChange_refGene) into a single field called merged_hgvs. It aggregates distinct HGVS annotations from both sources, concatenating them with a comma separator. This unified field allows for easier access to HGVS information regardless of the annotation source, facilitating downstream analysis and interpretation.

## NOMEN

NOMEN information (e.g. NOMEN, CNOMEN, PNOMEN...) from HGVS nomenclature field (see parameters help)

Extract NOMEN information (e.g. NOMEN, CNOMEN, PNOMEN...) from HGVS nomenclature field (see parameters help). This calculation parses the HGVS nomenclature field to extract specific NOMEN information such as NOMEN, CNOMEN, PNOMEN, etc., based on the parameters provided. It creates new INFO fields with the extracted NOMEN information for easier access and downstream analysis.

## NOMEN_SNPEFF

NOMEN information (e.g. NOMEN, CNOMEN, PNOMEN...) from HGVS nomenclature field (see parameters help)

Extract NOMEN information (e.g. NOMEN, CNOMEN, PNOMEN...) from HGVS nomenclature field (see parameters help). This calculation parses the HGVS nomenclature field to extract specific NOMEN information such as NOMEN, CNOMEN, PNOMEN, etc., specifically based on snpEff HGVS annotations field 'snpeff_hgvs'. It creates new INFO fields with the extracted NOMEN information for easier access and downstream analysis.

## RECREATE_INFO_FIELDS

Recreate INFO_tags, rename or remove tags

Recreate INFO_tags, rename or remove tags. This calculation allows for the recreation of INFO fields by renaming existing tags or removing them based on specified parameters. It provides flexibility in managing INFO fields, enabling users to customize their VCF annotations according to their analysis needs.

## RENAME_INFO_FIELDS

Rename or remove INFO/tags

Rename or remove INFO/tags. This calculation allows for the renaming of existing INFO fields or the removal of specific tags based on provided parameters. It helps in standardizing or customizing INFO fields in the VCF according to user requirements and facilitates downstream analysis. Use paramter section 'calculation_RENAME_INFO_FIELDS' in param.json to specify the INFO fields to rename or remove.

This example will rename INFO field 'ENOMEN' to 'exon' and remove INFO fields 'SPiP', 'SPiP_Alt' and 'SPiP_distSS' from the 'variants' table:

{"fields_to_rename": {"ENOMEN": "exon", "SPiP": null, "SPiP_Alt": null, "SPiP_distSS": null }, "table": "variants"}

## SNPEFF_ANN_EXPLODE

Explode snpEff annotations

Explode snpEff annotations from the ANN field, creating new INFO fields with prefix 'snpeff_' for each annotation. This calculation takes the ANN field from snpEff annotations, which may contain multiple annotations for a single variant, and explodes it so that each annotation is represented in separate INFO fields with a 'snpeff_' prefix. This allows for easier access to individual annotations and facilitates downstream analysis.

## SNPEFF_ANN_EXPLODE_JSON

Explode snpEff annotations in JSON format

Explode snpEff annotations from the ANN field, creating a new INFO field 'snpeff_json' in JSON format that contains all annotations. This calculation takes the ANN field from snpEff annotations, which may contain multiple annotations for a single variant, and explodes it into a structured JSON format stored in a new INFO field named 'snpeff_json'. This allows for easier access to all annotations in a single field and facilitates downstream analysis.

## SNPEFF_ANN_EXPLODE_UNIQUIFY

Explode snpEff annotations with uniquify values

Explode snpEff annotations from the ANN field, creating new INFO fields with prefix 'snpeff_uniquify_' for each annotation and ensuring unique values. This calculation takes the ANN field from snpEff annotations, which may contain multiple annotations for a single variant, and explodes it so that each annotation is represented in separate INFO fields with a 'snpeff_uniquify_' prefix, ensuring that only unique values are retained. This allows for easier access to individual annotations and facilitates downstream analysis.

## SNPEFF_EXTRACT

HGVS nomenclatures from snpEff annotation

Extract HGVS nomenclatures from snpEff annotation (field ANN) and create new INFO fields with prefix 'snpeff_' (e.g. snpeff_hgvs, snpeff_impact, snpeff_gene_name...). This calculation parses the ANN field from snpEff annotations, extracting relevant information such as HGVS nomenclatures, impact, gene name, etc., and creates new INFO fields with a 'snpeff_' prefix for easier access and downstream analysis.

## SNPEFF_HGVS

HGVS nomenclatures from snpEff annotation

Extract HGVS nomenclatures from snpEff annotation (field ANN) and create new INFO field 'snpeff_hgvs'. This calculation specifically targets the extraction of HGVS nomenclatures from the ANN field of snpEff annotations, creating a new INFO field named 'snpeff_hgvs' that contains the extracted HGVS information for easier access and downstream analysis.

## TRANSCRIPTS_ANNOTATIONS

Perform transcripts annotations and generate a transcripts table/view (using JSON parameters file)

Perform transcripts annotations and generate a transcripts table/view. This calculation allows for the inclusion of transcript annotations in both JSON and structured formats, providing flexibility in how transcript information is stored and accessed within the variant analysis pipeline. The specific format(s) used can be determined based on parameters specified in the JSON parameter file in 'transcripts' section (see help.parameters.md), or directly in the calculation parameters.

## TRANSCRIPTS_JSON

Perform transcripts annotations and export into INFO field in JSON format (field 'transcripts_json')

Perform transcripts annotations and generate a transcripts table/view. This calculation allows for the inclusion of transcript annotations in both JSON and structured formats, providing flexibility in how transcript information is stored and accessed within the variant analysis pipeline. The specific format(s) used can be determined based on parameters specified in the JSON parameter file in 'transcripts' section (see help.parameters.md). Then and add transcripts fields into INFO fields in structured format (field 'transcripts_ann'). This calculation allows for the inclusion of transcript annotations in a JSON format, facilitating the integration of transcript information into the variant analysis pipeline.

## TRANSCRIPTS_ANN

Perform transcripts annotations and export into INFO field in structured format (field 'transcripts_ann')

Perform transcripts annotations and generate a transcripts table/view. This calculation allows for the inclusion of transcript annotations in both JSON and structured formats, providing flexibility in how transcript information is stored and accessed within the variant analysis pipeline. The specific format(s) used can be determined based on parameters specified in the JSON parameter file in 'transcripts' section (see help.parameters.md). Then and add transcripts fields into INFO fields in structured format (field 'transcripts_ann'). This calculation allows for the inclusion of transcript annotations in a structured format, facilitating the integration of transcript information into the variant analysis pipeline.

## TRANSCRIPTS_EXPORT

Export transcripts table/view as a file (using JSON parameters file)

Export transcripts table/view as a file. This calculation allows for the export of the transcripts table or view generated from transcript annotations to an external file format such as TSV or CSV. The parameters for this calculation can be specified in the JSON parameter file, either in 'transcripts' section or directly in the calculation parameters (see help.parameters.md), including options for the output file path, format, and any filters to apply during export.

## TRANSCRIPTS_PRIORITIZATION

Prioritize transcripts with a prioritization profile (using JSON parameters file)

Prioritize transcripts with a prioritization profile. This calculation allows for the prioritization of transcripts based on a predefined profile specified in the JSON parameter file, either in 'transcripts' section or directly in the calculation parameters (see help.parameters.md), helping to identify the most relevant transcripts for further analysis.

## TRANSCRIPTS_PRIORITIZATION_STRICT

Prioritize transcripts with a prioritization profile (using JSON parameters file)

Prioritize transcripts with a prioritization profile, ensuring that all needed annotations are present and considered. This calculation allows for the prioritization of transcripts based on a predefined profile specified in the JSON parameter file, either in 'transcripts' section or directly in the calculation parameters (see help.parameters.md), helping to identify the most relevant transcripts for further analysis.

## TRIO

Inheritance for a trio family

Inheritance for a trio family. This calculation assesses the inheritance pattern of a variant within a trio family (father, mother, and child). It helps in understanding the transmission of genetic variants and identifying potential de novo mutations.

## VAF

Variant Allele Frequency (VAF) harmonization

Variant Allele Frequency (VAF) harmonization. This calculation normalizes the Variant Allele Frequency (VAF) across different samples and callers to ensure consistency in VAF representation. It helps in comparing VAF values across samples and callers by applying a harmonization process.

## VAF_STATS

Variant Allele Frequency (VAF) statistics

Variant Allele Frequency (VAF) statistics. This calculation provides statistical measures for the Variant Allele Frequency (VAF) across different samples. It helps in understanding the distribution and variability of VAF values.

## VARIANT_CHR_POS_ALT_REF

Create a variant ID with chromosome, position, alt and ref

This calculation generates a variant ID with chromosome, position, alt and ref, with format '#CHROM_POS_REF_ALT'. Useful for variant identification and comparison across datasets

## VARIANT_ID

Variant ID generated from variant position and type

Variant ID generated from variant position and variation (ref and alt) and type (SVTYPE). This calculation creates a unique identifier for each variant, facilitating variant tracking and comparison across datasets.

## VARTYPE

Variant type (e.g. SNV, INDEL, MNV, BND...)

Determine the type of variant based on its characteristics, such as the length of the reference and alternate alleles, and the presence of structural variant information. This calculation classifies variants into types like SNV, INDEL, MNV, BND, etc., which is essential for downstream analysis and interpretation.