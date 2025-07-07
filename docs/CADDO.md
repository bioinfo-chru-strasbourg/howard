# CADD Optimized pathogenicity thresholds for Non Coding regions of the genome

CONCERN : CADD Optimized threshold for Non Coding regions to Evaluate Risk of pathogeNicity

## Introduction

This documentation provides a step-by-step guide for setting up and running the CADD Optimized (CADDO) pathogenicity prediction workflow. The CADDO score is designed to classify genetic variants as Pathogenic or Benign on CADD PHRED scores and functional annotations from ANNOVAR. The guide covers requirements, database setup, configuration, and execution of the workflow, enabling users to efficiently annotate and interpret variants using optimized CADD thresholds.

The commands in this documentation will help you:

- Set up environment variables for paths and assemblies.
- Download and prepare CADD and ANNOVAR databases.
- Configure the annotation and calculation parameters.
- Run the annotation workflow on your VCF files using Howard.
- Query and inspect the annotated results.

## Requirements

To follow this guide and run the CADDO workflow, you will need:

- **Operating System:** Linux or macOS recommended.
- **Software:**
  - [Howard](https://github.com/lebechea/howard): The main tool for annotation and processing.
  - [curl](https://curl.se/): For downloading files from the internet.
  - [awk](https://www.gnu.org/software/gawk/): For text processing.
  - [gzip](https://www.gnu.org/software/gzip/), [bgzip/tabix](http://www.htslib.org/): For file compression and indexing.

- **Databases:**
  - CADD v1.7 for the appropriate genome assembly (e.g. GRCh38/hg38, GRCh37/hg19).
  - ANNOVAR reference gene database for the appropriate genome assembly (e.g. GRCh38/hg38, GRCh37/hg19).

**Note:** Ensure all tools are installed and available in your system's PATH. Sufficient disk space is required for database downloads and processing.

## Init

The following commands initialize key environment variables required for the workflow. These variables define the genome assembly version, the location of your databases, and the working data directory. Adjust these paths as needed for your system.

- `ASSEMBLY`: Specifies the genome assembly to use (e.g., `hg38` or `hg19`).
- `DATABASES`: Sets the directory where CADD and ANNOVAR databases will be stored.
- `DATA`: Defines the directory for input, output, and intermediate files.

Set these variables before proceeding with database downloads and workflow execution.

```bash
ASSEMBLY="hg38"
ASSEMBLY_GRCH="GRCh38"
DATABASES="${HOME}/howard/databases"
DATA="${HOME}/howard/data"
```

## Databases

This section describes how to download and prepare the required CADD and ANNOVAR databases for use in the CADDO workflow. You will set up the directory structure, retrieve the necessary files, and perform any required preprocessing steps to ensure compatibility with Howard. Follow the instructions carefully to ensure that your annotation and prediction steps have access to the correct and properly formatted data sources.

### CADD

These steps ensure that the CADD database is properly formatted and accessible for use in the CADDO variant annotation workflow.

### Set CADD Variables

Define environment variables to specify the storage location, version, and genome assembly for the CADD database.

```bash
CADD_FOLDER="$DATABASES/cadd/current"
CADD_VERSION="v1.7"
```

#### Download CADD Files

Create the appropriate directory structure, then download the CADD SNV data and its index file for your selected genome assembly. The SNV file is processed to add the `chr` prefix to chromosome names and is recompressed for compatibility with downstream tools.

```bash
mkdir -p $CADD_FOLDER/$ASSEMBLY
wget --retry-connrefused --waitretry=10 --tries=10 --continue "https://kircherlab.bihealth.org/download/CADD/$CADD_VERSION/$ASSEMBLY_GRCH/whole_genome_SNVs.tsv.gz" -O - | gzip -d | cut -f1-4,6 | awk 'BEGIN{OFS="\t"} {if($1 ~ /^#/) print $0; else print "chr"$0}' | bgzip -l1 -c > "$CADD_FOLDER/$ASSEMBLY/whole_genome.tsv.gz"
wget --retry-connrefused --waitretry=10 --tries=10 --continue "https://kircherlab.bihealth.org/download/CADD/$CADD_VERSION/$ASSEMBLY_GRCH/gnomad.genomes.r4.0.indel.tsv.gz" -O - | gzip -d | cut -f1-4,6 | awk 'BEGIN{OFS="\t"} $1 !~ /^#/{print "chr"$0}' | bgzip -l1 -c >> "$CADD_FOLDER/$ASSEMBLY/whole_genome.tsv.gz"
```

#### Convert CADD to Parquet

Generate a header file containing metadata and column definitions to ensure the processed CADD data is correctly interpreted during annotation and analysis. Generate the corresponding Parquet file to enable efficient CADD annotation.

```bash
echo "##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description=\"All filters passed\">
##INFO=<ID=PHRED,Number=1,Type=Float,Description=\"CADD PHRED Score\">
##contig=<ID=chr1>
##contig=<ID=chr10>
##contig=<ID=chr11>
##contig=<ID=chr12>
##contig=<ID=chr13>
##contig=<ID=chr14>
##contig=<ID=chr15>
##contig=<ID=chr16>
##contig=<ID=chr17>
##contig=<ID=chr18>
##contig=<ID=chr19>
##contig=<ID=chr2>
##contig=<ID=chr20>
##contig=<ID=chr21>
##contig=<ID=chr22>
##contig=<ID=chr3>
##contig=<ID=chr4>
##contig=<ID=chr5>
##contig=<ID=chr6>
##contig=<ID=chr7>
##contig=<ID=chr8>
##contig=<ID=chr9>
##contig=<ID=chrX>
##contig=<ID=chrY>
#CHROM	POS	REF	ALT	PHRED
" > $CADD_FOLDER/$ASSEMBLY/whole_genome.tsv.gz.hdr

howard query --input=$CADD_FOLDER/$ASSEMBLY/whole_genome.tsv.gz --output="$CADD_FOLDER/$ASSEMBLY/whole_genome.PHRED.parquet" --query="SELECT \"#Chrom\" AS '#CHROM', Pos AS POS, Ref as 'REF', Alt AS 'ALT', PHRED AS 'PHRED' FROM variants " --config='{"access": "RO"}' --parquet_partitions='#CHROM'
```

### ANNOVAR

Set up the directory for ANNOVAR databases and use the HOWARD tool to download the required reference gene annotation files for your selected genome assembly. This ensures that functional annotations (such as gene regions and variant effects) are available for downstream analysis.

- Create the ANNOVAR database directory.
- Download the `refGene` annotation files for the specified assembly.

The following command will automatically fetch and prepare the necessary ANNOVAR files:

```bash
ANNOVAR_FOLDER="$DATABASES/annovar/current"
howard databases --assembly="$ASSEMBLY" --download-annovar="$ANNOVAR_FOLDER" --download-annovar-files='refGene' --assembly="$ASSEMBLY"
```

## Process

This section outlines the complete workflow for processing a VCF file using the CADDO pipeline, from parameter file creation to annotation and calculation.

Workflow Overview:

1. **Create Calculation Configuration**:  
    Define the logic for classifying variants as Pathogenic, Benign, or Unknown based on CADD PHRED scores and ANNOVAR functional annotations. This is saved in a JSON file (e.g., `calculations_config.json`).

2. **Create Parameter File**:  
    Prepare a parameter file (e.g., `param.json`) that specifies the annotation sources (CADD Parquet and ANNOVAR refGene), the mapping of annotation fields, and the calculation configuration to be applied.

3. **Prepare Input VCF**:  
    Ensure your input VCF file is available and formatted correctly. You can create a test VCF as shown in the example.

4. **Run Annotation and Calculation**:  
    Use the `howard process` command to annotate the VCF with CADD and ANNOVAR data, and apply the CADDO calculation to classify each variant.

5. **Inspect Results**:  
    Query the annotated VCF to review the results and verify that variants have been classified according to the defined logic.

This workflow enables automated, reproducible variant annotation and classification using optimized CADD thresholds and functional context.

### Configuration

```bash
CALCULATIONS="$DATA/calculations_config.json"
echo '{
    "CADDO": {
      "type": "sql",
      "name": "CADDO",
      "description": "CADD Optimized pathogenicity prediction",
      "available": true,
      "output_column_name": "CADDO",
      "output_column_number": "1",
      "output_column_type": "String",
      "output_column_description": "Depending on CADD PHRED score (PHRED) and ANNOVAR variant location (Func_refGene), and CADD optimized thresholds, the variant is classified as Pathogenic 'P' or Benin 'B' or Unknown 'U'",
      "operation_query": [
        "CASE",
        "WHEN PHRED IS NULL THEN '\''U'\''",
        "WHEN list_contains(Func_refGene, '\''splicing'\'') THEN CASE WHEN PHRED >= 25.5 THEN '\''P'\'' ELSE '\''B'\'' END",
        "WHEN list_contains(Func_refGene, '\''ncRNA_exonic'\'') THEN CASE WHEN PHRED >= 11.44 THEN '\''P'\'' ELSE '\''B'\'' END",
        "WHEN list_contains(Func_refGene, '\''intronic'\'')  THEN CASE WHEN PHRED >= 12.76 THEN '\''P'\'' ELSE '\''B'\'' END",
        "WHEN list_contains(Func_refGene, '\''UTR5'\'') THEN CASE WHEN PHRED >= 16.79 THEN '\''P'\'' ELSE '\''B'\'' END",
        "WHEN list_contains(Func_refGene, '\''UTR3'\'') THEN CASE WHEN PHRED >= 11.08 THEN '\''P'\'' ELSE '\''B'\'' END",
        "WHEN PHRED >=25 THEN '\''P'\''",
        "ELSE '\''B'\''",
        "END"
      ],
      "info_fields": ["PHRED", "Func_refGene"],
      "operation_info": true
    }
}' > "$CALCULATIONS"
```

```bash
PARAM="$DATA/param.json"
echo '{
  "annotation": {
    "parquet": {
      "annotations": {
        "'$CADD_FOLDER/$ASSEMBLY/whole_genome.PHRED.parquet'": {
          "PHRED": "PHRED"
        }
      }
    },
    "annovar": {
      "annotations": {
        "refGene": {
          "Func_refGene": "Func_refGene"
        }
      },
      "options": {
        "genebase": " -splicing_threshold 3 "
      }
    }
  },
  "calculation": {
    "calculations": {
      "CADDO": null
    },
    "calculation_config": "'$CALCULATIONS'"
  }
}
' > "$PARAM"
```

### Run

```bash
MY_VCF="$DATA/my.vcf"
MY_ANOTATED_VCF="$DATA/my.annotated.vcf"

echo '##fileformat=VCFv4.1
##FILTER=<ID=PASS,Description="All filters passed">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##contig=<ID=chr1,length=249250621,assembly=hg38>
##contig=<ID=chr12,length=133851895,assembly=hg38>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	sample
chr1	28736	.	C	A	100	PASS		GT	1/1
chr1	35144	.	A	C	100	PASS		GT	0/1
chr12	120291839	.	T	TA	100	PASS		GT	0/1
' > $MY_VCF
```

```bash
howard process --input="$MY_VCF" --output="$MY_ANOTATED_VCF" --assembly="$ASSEMBLY" --param="$PARAM"
```

```bash
howard query --input="$MY_ANOTATED_VCF" --query="SELECT Func_refGene AS Location, CADDO AS Prediction, COUNT(*) AS Variant_Count FROM variants_view GROUP BY Func_refGene, CADDO ORDER BY Location, Prediction"
```
