<style>body {text-align: justify}</style>

# HOWARD User Guide

Highly Open and Valuable tool for Variant Annotation & Ranking for Discovery

HOWARD annotates and prioritizes genetic variations, calculates and normalizes annotations, translates files in multiple formats (e.g. vcf, tsv, parquet) and generates variants statistics.

HOWARD annotation is mainly based on a build-in Parquet annotation method, and external tools such as BCFTOOLS, ANNOVAR, snpEff and Exomiser (see docs, automatically downloaded if needed). Parquet annotation uses annotation database in VCF or BED format, in mutliple file format: Parquet/duckdb, VCF, BED, TSV, CSV, TBL, JSON.

HOWARD calculation processes variants information to calculate new information, such as: harmonizes allele frequency (VAF), extracts Nomen (transcript, cNomen, pNomen...) from HGVS fields with an optional list of personalized transcripts, generates VaRank format barcode.

HOWARD prioritization algorithm uses profiles to flag variants (as passed or filtered), calculate a prioritization score, and automatically generate a comment for each variants (example: 'polymorphism identified in dbSNP. associated to Lung Cancer. Found in ClinVar database').Prioritization profiles are defined in a configuration file. A profile is defined as a list of annotation/value, using wildcards and comparison options (contains, lower than, greater than, equal...). Annotations fields may be quality values (usually from callers, such as 'GQ', 'DP') or other annotations fields provided by annotations tools, such as HOWARD itself (example: COSMIC, Clinvar, 1000genomes, PolyPhen, SIFT). Multiple profiles can be used simultaneously, which is useful to define multiple validation/prioritization levels (example: 'standard', 'stringent', 'rare variants', 'low allele frequency').

HOWARD translates VCF format into multiple formats (e.g. VCF, TSV, Parquet), by sorting variants using specific fields (example : 'prioritization score', 'allele frequency', 'gene symbol'), including/excluding annotations/fields, including/excluding variants, adding fixed columns.

HOWARD generates statistics files with a specific algorithm, snpEff and BCFTOOLS.

HOWARD is multithreaded through the number of variants and by database (data-scaling).

## Table of Contents
- [Installation](#installation)
  - [Python](#python)
  - [Docker](#docker)
- [Databases](#databases)
    - [Databases tool](#databases-tool)
    - [Home-made databases](#home-made-databases)
- [Configuration](#configuration)
- [Parameters](#parameters)
- [Tools](#tools)
  - [Stats](#stats)
  - [Convert](#convert)
  - [Query](#query)
  - [Annotation](#annotation)
    - [Quick Annotation](#quick-annotation)
      - [Parquet annotation method](#parquet-annotation-method)
      - [External tools annotation](#external-tools-annotation)
        - [Annovar annotation](#annovar-annotation)
        - [snpEff annotation](#snpeff-annotation)
        - [Exomiser annotation](#exomiser-annotation)
        - [BCFTools annotation](#bcftools-annotation)
      - [Annotation combination](#annotation-combination)
    - [Annotation parameters](#annotation-parameters)
  - [Calculation](#calculation)
  - [Prioritization](#prioritization)
  - [Process](#process)

<!-- - [Annotation](#annotation)
- [Calculation](#calculation)
- [Prioritization](#prioritization)
- [Process](#process)
- [Query](#query)
- [Stats](#stats)
- [Convert](#convert) -->

# Installation

## Python

### Quick install

Install HOWARD using Python Pip tool:
```
python -m pip install -e .
```

Run HOWARD for help options:
```
howard --help
```

### GUI install

Install HOWARD Graphical User Interface using Python Pip tool with supplementary packages:
```
python -m pip install -r requirements-gui.txt
```

Run HOWARD Graphical User Interface as a tool:
```
howard gui
```

![HOWARD Graphical User Interface](../images/howard-gui.png "HOWARD Graphical User Interface")


## Docker

### Quick Start

In order to build images, launch default setup and create a persitent CLI (Command Line Inferface container), docker-compose command build images and launch services as containers.

```
docker-compose up -d
```

The persitent CLI contains external tools, such as:
| External tool | Description |
| -- | -- |
| [BCFTools](https://samtools.github.io/bcftools/) | Utilities for variant calling and manipulating VCFs and BCFs |
| [snpEff](https://pcingola.github.io/SnpEff/) | Genomic variant annotations, and functional effect prediction toolbox |
| [Annovar](https://annovar.openbioinformatics.org/) | Efficient software tool to utilize update-to-date information to functionally annotate genetic variants |
| [Exomiser](https://www.sanger.ac.uk/tool/exomiser/) | Program that finds potential disease-causing variants from whole-exome or whole-genome sequencing data |




### Setup container

Docker service HOWARD-setup creates HOWARD image and download useful databases to start with HOWARD tools. 

List of databases downloaded in HOWARD setup for hg19 assembly (see [Databases section](#databases) for more information):
- Genome
- Annovar ([refGene](https://genome.ucsc.edu/), [COSMIC](https://cancer.sanger.ac.uk/cosmic))
- snpEff
- refSeq
- dbNSFP
- AlphaMissense
- dnSNP

To avoid databases download (see [Databases section](#databases) to download manually), just start [Command Line Interface](#command-line-interface)

### Command Line Interface

A Command Line Interface container (HOWARD-CLI) is started with host data and databases folders mounted (by default both in ${HOME}/HOWARD folder). To manually start CLI container:

```
docker-compose up -d HOWARD-CLI
```

To use HOWARD tools within HOWARD-CLI container:

```
docker exec -ti HOWARD-CLI bash
howard --help
```

To run a command into HOWARD-CLI container:

```
docker exec HOWARD-CLI <howard command>
```

> Example: Query of an existing VCF:
> 
> ```
> docker exec HOWARD-CLI howard query --input=/tool/tests/data/example.vcf.gz --query='SELECT * FROM variants'
> ```

> Example: VCF annotation using HOWARD-CLI (snpEff and ANNOVAR databases will be automatically downloaded), and query list of genes with HGVS:
> 
> ```
> docker exec --workdir=/tool HOWARD-CLI howard process --config=config/config.json --param=config/param.json --input=tests/data/example.vcf.gz --output=/tmp/example.process.tsv --query='SELECT "NOMEN", "PZFlag", "PZScore", "PZComment" FROM variants ORDER BY "PZScore" DESC' --explode_infos
> ```

See [HOWARD Help](docs/help.md) for more options.

Let's play within Docker HOWARD-CLI service!

### Tests

In order to test HOWARD within Docker, use this command:
```
docker exec -ti HOWARD-CLI bash
cd /tool
# Basic test
coverage run -m pytest .
# Debug test
coverage run -m pytest . -x -v --log-cli-level=DEBUG --capture=tee-sys
```

# Databases

## Databases tool

Multiple databases can be automatically downloaded with databases tool, such as:

| database                                                          | description                                                                                                                                                                                                                                                                    |
| ----------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| [Genome](https://genome.ucsc.edu/cgi-bin/hgGateway)               | Genome Reference Consortium Human                                                                                                                                                                                                                                              |
| [Annovar](https://annovar.openbioinformatics.org/en/latest/)      | ANNOVAR is an efficient software tool to utilize update-to-date information to functionally annotate genetic variants detected from diverse genomes                                                                                                                            |      |
| [snpEff](https://pcingola.github.io/SnpEff/)                      | Genetic variant annotation, and functional effect prediction toolbox                                                                                                                                                                                                           |
| [refSeq](https://www.ncbi.nlm.nih.gov/refseq/)                    | A comprehensive, integrated, non-redundant, well-annotated set of reference sequences including genomic, transcript, and protein                                                                                                                                               |
| [dbSNP](https://www.ncbi.nlm.nih.gov/snp/)                        | dbSNP contains human single nucleotide variations, microsatellites, and small-scale insertions and deletions along with publication, population frequency, molecular consequence, and genomic and RefSeq mapping information for both common variations and clinical mutations |
| [dbNSFP](https://sites.google.com/site/jpopgen/dbNSFP)            | dbNSFP is a database developed for functional prediction and annotation of all potential non-synonymous single-nucleotide variants (nsSNVs) in the human genome                                                                                                                |
| [AlphaMissense](https://github.com/google-deepmind/alphamissense) | AlphaMissense model implementation                                                                                                                                                                                                                                             |
| [Exomiser](https://www.sanger.ac.uk/tool/exomiser/)               | The Exomiser is a Java program that finds potential disease-causing variants from whole-exome or whole-genome sequencing data                                                                                                                                                  |

> Example: Download Multiple databases in the same time for assembly 'hg19' (can take a while):
> 
> ```
> howard databases \
>    --assembly=hg19 \
>    --download-genomes=/databases/genomes/current --download-genomes-provider=UCSC --download-genomes-contig-regex='chr[0-9XYM]+$' \
>    --download-annovar=/databases/annovar/current --download-annovar-files='refGene,cosmic70,nci60' \
>    --download-snpeff=/databases/snpeff/current \
>    --download-refseq=/databases/refseq/current --download-refseq-format-file='ncbiRefSeq.txt' \
>    --download-dbnsfp=/databases/dbnsfp/current --download-dbnsfp-release='4.4a' --download-dbnsfp-subdatabases \
>    --download-alphamissense=/databases/alphamissense/current \
>    --download-exomiser=/databases/exomiser/current \
>    --download-dbsnp=/databases/dbsnp/current --download-dbsnp-vcf --threads=8
> ```

See [HOWARD Help Databases tool](docs/help.md#databases-tool) for more information.

## Home-made Databases

Databases can be generated using an home-made existing annotation file and [HOWARD convert](#convert) tool. The home-made annotation file need to contain specific fields (depending on the annotation type):

- variant annotation: '#CHROM', 'POS', 'ALT', 'REF'
- region annotation: '#CHROM', 'START', 'STOP'

An home-made existing annotation file can be converted into multiple formats (e.g. Parquet, VCF, TSV), but it's strongly suggested to use Parquet format.

After convertion, the database file is associated with a 'header' file ('.hdr'), in VCF header format, to describe annotations within the database. Use the 'header' file to describe annotation fields/columns present in the existing file. An Home-made annotation file in VCF format which is converted in another format will keep all annotation information from the initial VCF header.

Note that a VCF can be directly used as a database (annotation field information within the header of the VCF file). Also, an home-made existing annotation file can be used as a database, but will not be totaly compliant due to the lack of annotation information ('header' will be generated by default).

See [HOWARD Help Convert tool](help.md#convert-tool) for more information.

## VCF and Parquet from Annovar

See [HOWARD Help From ANNOVAR tool](help.md#from_annovar-tool) tool for more information (under development).

# Configuration

HOWARD Configuration JSON file defined default configuration regarding resources (e.g. threads, memory), settings (e.g. verbosity, temporary files), default folders (e.g. for databases) and paths to external tools.

See [HOWARD Configuration JSON](help.config.md) for more information.

Configuration file example:
```
{
  "threads": 8,
  "memory": null,
  "verbosity": "warning",
  "folders": {
    "databases": {
      "annotations": [
        "/databases/annotations/current/",
        "/databases/dbnsfp/current/",
        "/databases/dbsnp/current/"
      ],
      "parquet": ["/databases/annotations/current/"],
      "bcftools": ["/databases/annotations/current/"],
      "annovar": "/databases/annovar/current/",
      "snpeff": "/databases/snpeff/current/",
      "exomiser": "/databases/exomiser/current/"
    }
  },
  "tools": {
    "bcftools": "/tools/bcftools/current/bin/bcftools",
    "bgzip": "/tools/htslib/current/bin/bgzip",
    "snpeff": "/tools/snpeff/current/bin/snpEff.jar",
    "annovar": "/tools/annovar/current/bin/table_annovar.pl"
    "exomiser": "/tools/exomiser/current/bin/exomiser-cli-13.2.0.jar",
    "java": "/usr/bin/java",
  }
}

```

# Parameters

HOWARD Parameters JSON file defined parameters to process annotations, prioritization, calculations, convertions and queries.

See [HOWARD Parameters JSON](help.param.md) for more information.

Parameters file example:
```
{
  "annotation": {
    "snpeff": {
      "options": "-lof -hgvs -oicr -noShiftHgvs -spliceSiteSize 3 "
    },
    "annovar": {
      "annotations": {
        "refGene": {
          "INFO": null
        }
      },
      "options": {
        "genebase": "-hgvs -splicing_threshold 3 ",
        "intronhgvs": 10
      }
    },
    "parquet": {
      "annotations": {
        "tests/databases/annotations/current/hg19/avsnp150.parquet": {
          "INFO": null
        },
        "tests/databases/annotations/current/hg19/dbnsfp42a.parquet": {
          "INFO": null
        },
        "tests/databases/annotations/current/hg19/gnomad211_genome.parquet": {
          "INFO": null
        }
      }
    },
    "bcftools": {
      "annotations": {
        "tests/databases/annotations/current/hg19/cosmic70.vcf.gz": {
          "INFO": null
        }
      }
    }
  },
  "calculation": {
    "vartype": null,
    "snpeff_hgvs": null,
    "NOMEN": {
      "options": {
        "hgvs_field": "snpeff_hgvs",
        "transcripts": "tests/data/transcripts.tsv"
      }
    },
    "VAF": ""
  },
  "prioritization": {
    "config_profiles": "config/prioritization_profiles.json",
    "pzfields": ["PZScore", "PZFlag", "PZComment"],
    "profiles": ["default", "GERMLINE"]
  },
  "explode_infos": True,
  "threads": 8
}
```

# Tools

## Stats

Statistics on genetic variations, such as: number of variants, number of samples, statistics by chromosome, genotypes by samples, annotations...
Theses statsitics can be applied to VCF files and all database annotation files.

> Example: Show example VCF statistics and brief overview
> ```
> howard stats --input=tests/data/example.vcf.gz
> ```

See [HOWARD Help Stats tool](help.md#stats-tool) for more information.

## Convert

Convert genetic variations file to another format. Multiple format are available, such as usual and official VCF format, but also other formats such as TSV, CSV, TBL, JSON and Parquet/duckDB. These formats need a header '.hdr' file to take advantage of the power of howard (especially through INFO/tag definition), and using howard convert tool automatically generate header file fo futher use (otherwise, an default '.hdr' file is generated).

> Example: Translate VCF into TSV, export INFO/tags into columns, and show output file
> ```
> howard convert --input=tests/data/example.vcf.gz --output=/tmp/example.tsv --explode_infos && \
> cat /tmp/example.tsv
> ```

> Example: Translate VCF into parquet
> ```
> howard convert --input=tests/data/example.vcf.gz --output=/tmp/example.parquet
> ```

See [HOWARD Help Convert tool](help.md#convert-tool) for more information.

## Query

Query genetic variations in SQL format. Data can be loaded into 'variants' table from various formats (e.g. VCF, TSV, Parquet...). Using --explode_infos allow query on INFO/tag annotations. SQL query can also use external data within the request, such as a Parquet file(s).

> Example: Select variants in VCF with REF and POS fields
> ```
> howard query --input=tests/data/example.vcf.gz \
>    --query="SELECT * FROM variants WHERE REF = 'A' AND POS < 100000"
> ```

> Example: Select variants in VCF with INFO Tags criterions
> ```
> howard query --input=tests/data/example.vcf.gz --explode_infos \
>    --query='SELECT "#CHROM", POS, REF, ALT, DP, CLNSIG, sample2, sample3 
>             FROM variants 
>             WHERE DP >= 50 OR CLNSIG NOT NULL 
>             ORDER BY CLNSIG DESC, DP DESC'
> ```

> Example: Query a Parquet file with specific columns (e.g. from VCF convertion to Parquet)
> ```
> howard query \
>    --query="SELECT * \
>             FROM 'tests/databases/annotations/current/hg19/dbnsfp42a.parquet' \
>             WHERE \"INFO/Interpro_domain\" NOT NULL \
>             ORDER BY \"INFO/SiPhy_29way_logOdds_rankscore\" DESC"
> ```

> Example: Query multiple Parquet files, merge INFO columns, and extract as TSV (in VCF format)
> ```
> howard query \
>    --query="
>       SELECT \
>          \"#CHROM\" AS \"#CHROM\", \
>          POS AS POS, \
>          '' AS ID, \
>          REF AS REF, \
>          ALT AS ALT, \
>          '' AS QUAL, \
>          '' AS FILTER, \
>          STRING_AGG(INFO, ';') AS INFO \
>       FROM 'tests/databases/annotations/current/hg19/*.parquet' \
>       GROUP BY \"#CHROM\", POS, REF, ALT" \
>    --output=/tmp/full_annotation.tsv \
>    --include_header
> ```


See [HOWARD Help Query tool](help.md#query-tool) for more information.

## Annotation

Annotation is mainly based on a build-in Parquet annotation method, using annotation database file (in multiple format such as Parquet, duckdb, VCF, BED, TSV, JSON).

These annotation databases can be automatically downloaded using [HOWARD Databases tool](docs/help.md#databases-tool) and manually generated using existing annotation files and [HOWARD Convert tool](help.md#convert-tool). Annotation databases need a header file (`.hdr`) to describe each annotation in the database. However, a default header will be generated if no header file is associated to the annotation database file.

Moreover, some external annotation tools are integrated into HOWARD to extend annotation with their own options and databases.

HOWARD annotation tool can use annotation databases files in 2 differents ways: [Quick annotation](#quick-annotation) and [Annotation Parameters](#annotation-parameters) JSON file.

### Quick annotation

#### Parquet annotation method

Quick annotation allows to annotates by simply listing annotation databases (in multiple format). 

##### Parquet annotation with path

These annotation databases are defined with their full path (e.g. `/full/path/to/my.database.parquet`) or relative path (e.g. `databases/my.database.parquet`). A list (separator `:` or `+`) of annotation databases files can be used.

> Example: VCF annotation with full path Parquet databases
> ```
> howard annotation \
>    --input=tests/data/example.vcf.gz \
>    --annotations='/tool/tests/databases/annotations/current/hg19/dbnsfp42a.parquet' \
>    --output=/tmp/example.howard.vcf.gz
> ```

> Example: VCF annotation with relative path VCF databases
> ```
> cd /tool && howard annotation --input=tests/data/example.vcf.gz --annotations='tests/databases/annotations/current/hg19/cosmic70.vcf.gz' --output=/tmp/example.howard.vcf.gz
> ```

> Example: VCF annotation with relative path BED databases
> ```
> cd /tool && \
> howard annotation \
>    --input=tests/data/example.vcf.gz \
>    --annotations='tests/databases/annotations/current/hg19/refGene.bed.gz' \
>    --output=/tmp/example.howard.vcf.gz
> ```


> Example: VCF annotation with 3 annotation databases files
> ```
> howard annotation \
>    --input=tests/data/example.vcf.gz \
>    --annotations='/tool/tests/databases/annotations/current/hg19/dbnsfp42a.parquet+tests/databases/annotations/current/hg19/cosmic70.vcf.gz+tests/databases/annotations/current/hg19/refGene.bed.gz' \
>    --output=/tmp/example.howard.vcf.gz
> ```

##### Parquet annotation with annotation folder

If annotation folder is configured in [HOWARD Configuration JSON](help.config.md), just mention the annotation database basename file. Annotation database file will be found (depending of the assembly).

> Example: VCF annotation with Parquet and VCF databases, with annotation database defined in JSON configuration (as a string)
> ```
> howard annotation \
>    --input=tests/data/example.vcf.gz \
>    --annotations='dbnsfp42a.parquet,cosmic70.vcf.gz' \
>    --config='{"folders": {"databases": {"annotations": ["/tool/tests/databases/annotations/current"]}}}' \
>    --output=/tmp/example.howard.vcf.gz 
> ```

##### Full annotation 

In order to annotate with all available annotation databases, the keyword `ALL` will auto-detect files in the databases annotation folder. The option `format` (defaut `parquet`) can filter annotation databases by listing (separator `+`) desired formats (such as `parquet`, `vcf`). The option `release` (default `current`) is able to scan annotation databases in one or more specific releases in a list (separator `+`). See [HOWARD Configuration JSON - Folders - Databases](help.config.md#foldersdatabases) for more information about databases structure.

> Example: VCF annotation with all available database annotation files in Parquet format (within the database annotation folder in configuration):
> ```
> howard annotation \
>    --input=tests/data/example.vcf.gz \
>    --assembly='hg19' \
>    --annotations='ALL:format=parquet+vcf:release=current' \
>    --config='{"folders": {"databases": {"annotations": ["/tool/tests/databases/annotations/current"]}}}' \
>    --output=/tmp/example.howard.tsv
> ```

#### External tools annotation

External annotation tools are also available, such as BCFTOOLS, Annovar, snpEff and Exomiser. Annovar, snpEff and Exomiser databases are automatically downloaded (see [HOWARD Help Databases tool](docs/help.md#databases-tool)). Quick annotation allows to annotates by simply defining external tools keywords.

##### BCFTools annotation

For BCFTools, use HOWARD keyword `bcftools` and list (separator `:` or `+`) annotation databases with format such as VCF or BED. More options are available using [HOWARD Parameters JSON](docs/help.param.md) file.

> Example: VCF annotation with Cosmic VCF databases and refGene BED database
> ```
> howard annotation \
>    --input=tests/data/example.vcf.gz \
>    --annotations='bcftools:tests/databases/annotations/current/hg19/cosmic70.vcf.gz+tests/databases/annotations/current/hg19/refGene.bed.gz' \
>    --output=/tmp/example.howard.vcf.gz
> ```

##### Annovar annotation

For Annovar tool, use HOWARD keyword `annovar` and mention specific Annovar database keywords. More options are available using [HOWARD Parameters JSON](docs/help.param.md) file.

> Example: VCF annotation with Annovar refGene and cosmic70
> ```
> howard annotation \
>    --input=tests/data/example.vcf.gz \
>    --annotations='annovar:refGene:cosmic70' \
>    --output=/tmp/example.howard.tsv
> ```

##### snpEff annotation

For snpEff tool, use HOWARD keyword `snpeff`. No options are available for quick annotation with snpEff, see [HOWARD Parameters JSON](docs/help.param.md) for more options.

> Example: VCF annotation with snpEff 
> ```
> howard annotation \
>    --input=tests/data/example.vcf.gz \
>    --annotations='snpeff' \
>    --output=/tmp/example.howard.tsv
> ```

##### Exomiser Annotation

For Exomiser tool, use HOWARD keyword `exomiser`. A list of options can be provided as key-value format. More options are available using [HOWARD Parameters JSON](docs/help.param.md) file.

> Example: VCF annotation with Exomiser (exome preset, list of HPO terms, transcript as refseq and release 2109)
> ```
> howard annotation \
>    --input=tests/data/example.vcf.gz \
>    --annotations='exomiser:preset=exome:hpo=0001156+0001363+0011304+0010055:transcript_source=refseq:release=2109' \
>    --output=/tmp/example.howard.tsv
> ```

#### Annotation combination

Quick annotation allows to combine annotations, from build-in Parquet method and external tools. Simply use a list with a comma separator.

> Example: VCF annotation with build-in Parquet method and external tools (Annovar, snpEff and Exomiser)
> ```
> howard annotation \
>    --input=tests/data/example.vcf.gz \
>    --annotations='tests/databases/annotations/current/hg19/dbnsfp42a.parquet,tests/databases/annotations/current/hg19/cosmic70.vcf.gz,annovar:cosmic70,snpeff,exomiser:preset=exome|hpo=0001156-0001363-0011304-0010055' \
>    --output=/tmp/example.howard.tsv
> ```

See [HOWARD Help Annotation tool](help.md#annotation-tool) tool for more information.

### Annotation parameters

All annotation parameters can be defined in [HOWARD Parameters JSON](help.param.md) file. All annotations can be combined (bild-in parquet method and external tools annotation), and options can be detailed in a full JSON file format, including selection of annotation database columns (and rename them) and specific options for external tools.

## Calculation

See [HOWARD Help Calculation tool](help.md#calculation-tool) tool for more information.

## Prioritization

See [HOWARD Help Prioritization tool](help.md#prioritization-tool) tool for more information.

## Process

See [HOWARD Help Process tool](help.md#process-tool) tool for more information (under development).
