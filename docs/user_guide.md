# HOWARD User Guide

![HOWARD Graphical User Interface](../images/icon.ico "HOWARD Graphical User Interface")

Highly Open and Valuable tool for Variant Annotation & Ranking toward genetic Discovery

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
    - [CSV, TSV and JSON convert](#csv-tsv-and-json-convert)
    - [Parquet and duckDB convert](#parquet-and-duckdb-convert)
    - [BED convert](#bed-convert)
    - [Partitioning](#partitioning)
    - [Explode INFO Tags](#explode-info-tags)
  - [Query](#query)
    - [Variants file](#variants-file)
      - [Loading data](#loading-data)
      - [External file](#external-file)
    - [Other files](#other-files)
      - [Parquet format](#parquet-format)
      - [CSV format](#csv-format)
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

See [HOWARD Help](help.md) for more options.

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
>    --download-genomes=~/howard/databases/genomes/current --download-genomes-provider=UCSC --download-genomes-contig-regex='chr[0-9XYM]+$' \
>    --download-annovar=~/howard/databases/annovar/current --download-annovar-files='refGene,cosmic70,nci60' \
>    --download-snpeff=~/howard/databases/snpeff/current \
>    --download-refseq=~/howard/databases/refseq/current --download-refseq-format-file='ncbiRefSeq.txt' \
>    --download-dbnsfp=~/howard/databases/dbnsfp/current --download-dbnsfp-release='4.4a' --download-dbnsfp-subdatabases \
>    --download-alphamissense=~/howard/databases/alphamissense/current \
>    --download-exomiser=~/howard/databases/exomiser/current \
>    --download-dbsnp=~/howard/databases/dbsnp/current --download-dbsnp-vcf --threads=8
> ```

See [HOWARD Help Databases tool](help.md#databases-tool) for more information.

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
        "~/howard/databases/annotations/current/",
        "~/howard/databases/dbnsfp/current/",
        "~/howard/databases/dbsnp/current/"
      ],
      "parquet": ["~/howard/databases/annotations/current/"],
      "bcftools": ["~/howard/databases/annotations/current/"],
      "annovar": "~/howard/databases/annovar/current/",
      "snpeff": "~/howard/databases/snpeff/current/",
      "exomiser": "~/howard/databases/exomiser/current/"
    }
  },
  "tools": {
    "bcftools": "~/howard/tools/bcftools/current/bin/bcftools",
    "bgzip": "~/howard/tools/htslib/current/bin/bgzip",
    "snpeff": "~/howard/tools/snpeff/current/bin/snpEff.jar",
    "annovar": "~/howard/tools/annovar/current/bin/table_annovar.pl"
    "exomiser": "~/howard/tools/exomiser/current/bin/exomiser-cli-13.2.0.jar",
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
    "parquet": {
      "annotations": {
        "tests/databases/annotations/current/hg19/avsnp150.parquet": {
          "INFO": null
        },
        "tests/databases/annotations/current/hg19/dbnsfp42a.parquet": {
          "Polyphen2_HDIV_pred": "PolyPhen",
          "Polyphen2_HDIV_score": "PolyPhen_score",
          "REVEL_score": null
        },
        "tests/databases/annotations/current/hg19/gnomad211_genome.parquet": {
          "INFO": null
        }
      }
    },
    "bcftools": {
      "annotations": {
        "tests/databases/annotations/current/hg19/cosmic70.vcf.gz": {
          "cosmic70": "COSMIC"
        }
      }
    },
    "annovar": {
      "annotations": {
        "refGene": {
          "INFO": null
        },
        "cosmic70": {
          "cosmic70": "COSMIC_annovar"
        }
      },
      "options": {
        "genebase": "-hgvs -splicing_threshold 3 ",
        "intronhgvs": 10
      }
    },
    "snpeff": {
      "options": "-lof -hgvs -oicr -noShiftHgvs -spliceSiteSize 3 "
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

Generates statistics on genetic variations, such as number of variants, number of samples, statistics by chromosome, genotypes by samples, annotations. Theses statsitics can be applied to VCF files from all database annotation file formats. Statistics can be wrote into files in Markdown and JSON format (resp. `--stats_md` and `--stats_json` parameter).

See [HOWARD Help Stats tool](help.md#stats-tool) for more information.

> Example: Show example VCF statistics and brief overview
> ```
> howard stats \
>    --input=tests/data/example.vcf.gz
> ```

> Example: Show example VCF statistics and generate a file in JSON and Markdown formats
> ```
> howard stats \
>    --input=tests/data/example.vcf.gz \
>    --stats_json=/tmp/stats.json \
>    --stats_md=/tmp/stats.md
> 
> cat /tmp/stats.json /tmp/stats.md
> ```
> ```
> {
>     "Infos": {
>         "Input file": "tests/data/example.vcf.gz",
>         "Number of variants": 7,
>         "Number of samples": 4,
>         "Number of INFO fields": 5,
>         "Number of FORMAT fields": 7
>     },
>     "Variants": {
>         "Number of variants by chromosome": {
>             "1": {
>                 "CHROM": "chr1",
>                 "count": 6,
>                 "percent": 0.8571428571428571
>             },
>             "0": {
>                 "CHROM": "chr7",
>                 "count": 1,
>                 "percent": 0.14285714285714285
>             }
>         },
> ...
> ```
> ```
> ...
> ## Variants
> ### Number of variants by chromosome
> | CHROM   |   count |   percent |
> |:--------|--------:|----------:|
> | chr1    |       6 |  0.857143 |
> | chr7    |       1 |  0.142857 |
> ### Counts
> | Type   |   count |
> |:-------|--------:|
> | Total  |       7 |
> | SNV    |       7 |
> | MNV    |       0 |
> | InDel  |       0 |
> ...
> ```

> Example of statistics in Markdown output

- Input file: tests/data/example.vcf.gz
- Number of variants: 7
- Number of samples: 4
- Number of INFO fields: 5
- Number of FORMAT fields: 7

| CHROM   |   count |   percent |
|:--------|--------:|----------:|
| chr1    |       6 |  0.857143 |
| chr7    |       1 |  0.142857 |

| Type   |   count |
|:-------|--------:|
| Total  |       7 |
| SNV    |       7 |
| MNV    |       0 |
| InDel  |       0 |


## Convert

Convert genetic variations file to another format. Multiple format are available, such as usual and official VCF format, but also other formats such as TSV, CSV, TBL, JSON and Parquet/duckDB. These formats need a header '.hdr' file to take advantage of the power of howard (especially through INFO/tag definition), and using howard convert tool automatically generate header file fo futher use (otherwise, an default '.hdr' file is generated).

Multiple options are available, such as explode VCF INFO/tags (parameter `--explode_infos`, see [HOWARD Help query - Explode infos](help.md#explode-infos)), order by columns, include header within file (only TSV and CSV format), or use partitioning into multiple files within a folder. See [HOWARD Help Convert tool](help.md#convert-tool) for more information.

### CSV, TSV and JSON convert

To convert a file (multiple formats) into another flat file, such as CSV (tab-delimiter) and TSV (comma-delimiter), or JSON format, simply name output file with desired extension. Use `.gz` extension to compress file.

> Example: Convert VCF into TSV and show output file
> ```
> howard convert \
>    --input=tests/data/example.vcf.gz \
>    --output=/tmp/example.tsv
> cat /tmp/example.tsv
> ```

> Example: Convert TSV into VCF and show output file
> ```
> howard convert \
>    --input=tests/data/example.tsv \
>    --output=/tmp/example.vcf
> cat /tmp/example.vcf
> ```

> Example: Convert VCF into CSV and compress file
> ```
> howard convert \
>    --input=tests/data/example.vcf.gz \
>    --output=/tmp/example.csv.gz
> ```

> Example: Convert VCF into JSON and compress file
> ```
> howard convert \
>    --input=tests/data/example.vcf.gz \
>    --output=/tmp/example.json.gz
> ```

### Parquet and duckDB convert

Files can be format into Parquet and duckDB format. For duckDB format, a duckDB database will be created with a `variants` table.

> Example: Convert VCF into parquet
> ```
> howard convert \
>    --input=tests/data/example.vcf.gz \
>    --output=/tmp/example.parquet
> ```

> Example: Convert VCF into duckDB
> ```
> howard convert \
>    --input=tests/data/example.vcf.gz \
>    --output=/tmp/example.duckdb
> ```

### BED convert

To convert into BED format, input file needs mandatory columns chrommosome "#CHROM", and positions columns "START"/"END" or "POS" (corresponfing to regions with a uniq nucleotide). Header will be automatically included to decribe all columns. 

> Example: Convert BED in TSV format into BED format
> ```
> howard convert \
>    --input=tests/data/example.bed.tsv \
>    --output=/tmp/example.bed
> ```


> Example: Convert BED in Parquet format into BED format
> ```
> howard convert \
>    --input=tests/data/example.bed.parquet \
>    --output=/tmp/example.bed
> ```

HOWARD input file format does not allow BED format. To read a BED file and export into a CSV or Parquet file format, see [Query BED format](#bed-format) section.

### Partitioning

Partitioning (or Hive partitioning) is a partitioning strategy that is used to split a table into multiple files based on partition keys. The files are organized into folders. Within each folder, the partition key has a value that is determined by the name of the folder (see [duckDB hive partitioning](https://duckdb.org/docs/data/partitioning/hive_partitioning.html)).

Simply list columns as keys to process partitioning. Use 'None' (string) for NO partition but split parquet files into a folder. The partitioning is available for all format (e.g. Parquet, TSV, JSON, except duckDB format), by naming output file with desired extension. This option is faster parallel writing, but memory consuming, and also is faster reading.

> Example: Convert VCF into partitioned Parquet and show tree structure
> ```
> howard convert \
>    --input=tests/data/example.vcf.gz \
>    --output=/tmp/example.partitioned.parquet \
>    --parquet_partitions="#CHROM"
> 
> tree /tmp/example.partitioned.parquet
> ```
> ```
> /tmp/example.partitioned.parquet
> ├── #CHROM=chr1
> │   └── fe61bf182de640f8840270a527aa1582-0.parquet
> └── #CHROM=chr7
>     └── fe61bf182de640f8840270a527aa1582-0.parquet
> ```

> Example: Convert VCF into partitioned TSV (compressed) and show tree structure (files are named with `.csv` extension, but are tab-delimited and compressed)
> ```
> howard convert \
>    --input=tests/data/example.vcf.gz \
>    --output=/tmp/example.partitioned.tsv.gz \
>    --parquet_partitions="#CHROM,REF"
> 
> tree /tmp/example.partitioned.tsv.gz
> ```
> ```
> /tmp/example.partitioned.tsv
> ├── #CHROM=chr1
> │   └── REF=A
> │       └── data_0.csv
> └── #CHROM=chr7
>     └── REF=G
>         └── data_0.csv
> ```

### Explode INFO tags

Use `--explode_infos` parameter to extract all INFO tags (i.e. annotations) into columns (see [HOWARD Help query - Explode infos](help.md#explode-infos)).

> Example: Convert VCF into TSV, export INFO/tags into columns, and show output file
> ```
> howard convert \
>    --input=tests/data/example.vcf.gz \
>    --explode_infos \
>    --output=/tmp/example.tsv
> 
> cut /tmp/example.tsv -f1-4,7,15
> ```
> ```
> #CHROM  POS       REF  ALT  FILTER  CLNSIG
> chr1    28736     A    C    PASS    pathogenic
> chr1    35144     A    C    PASS    non-pathogenic
> chr1    69101     A    G    PASS
> chr1    768251    A    G    PASS
> ...
> ```

## Query

Query tool provides a simple way to query genetic variations in SQL format. Data are loaded into 'variants' table from various formats (e.g. VCF, TSV, Parquet...). Using `--explode_infos` allows querying on INFO/tag annotations. SQL query can also use external file within the request, such as a Parquet file(s), or TSV files.

See [HOWARD Help Query tool](help.md#query-tool) for more information.

### Variants file

#### Loading data

Query tool is able to read variants (i.e. VCF) or regions files (i.e. BED) files, in various format (e.g. VCF, BED, Parquet, TSV, JSON), using `--input` parameter. This allows to load data to perfom actions, such as explode VCF INFO/tags (parameter `--explode_infos`, see [HOWARD Help query - Explode infos](help.md#explode-infos)) in columns to be easier querying. Each columns format (e.g. string, integer) are automatically detected to be used in a SQL query.

> Example: Select variants in VCF with REF and POS fields filter
> ```
> howard query \
>    --input=tests/data/example.vcf.gz \
>    --query="SELECT * FROM variants WHERE REF = 'A' AND POS < 100000"
> ```

> Example: Select variants in VCF with INFO Tags criteria filters
> ```
> howard query \
>    --input=tests/data/example.vcf.gz \
>    --explode_infos \
>    --query='SELECT "#CHROM", POS, REF, ALT, DP, CLNSIG, sample2, sample3 
>             FROM variants 
>             WHERE DP >= 50 OR CLNSIG NOT NULL 
>             ORDER BY CLNSIG DESC, DP DESC'
> ```

> Example: Select variants in VCF and generate VCF output with variants
> ```
> howard query \
>    --input=tests/data/example.vcf.gz \
>    --output=/tmp/example.filtered.vcf \
>    --explode_infos \
>    --query='SELECT "#CHROM", POS, REF, ALT, QUAL, FILTER, INFO 
>             FROM variants 
>             WHERE DP >= 50 OR CLNSIG NOT NULL'
> ```

#### External file

Variants files can be used directly within the query, espacially if they already contain variants information (e.g. "#CHROM", "POS", "REF", "ALT") and annotations as columns.

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

### Other files

Whatever the external file format, if it is compatible with duckDB, query tool is able to query data (see [duckDB Select Statement](https://duckdb.org/docs/sql/statements/select)).

#### Parquet format

Simply use Parquet file path within the query (as descibe above).

> Example: Query a Parquet file with specific columns (e.g. from VCF convertion to Parquet)
> ```
> howard query \
>    --query="SELECT * \
>             FROM 'tests/databases/annotations/current/hg19/dbnsfp42a.parquet' \
>             WHERE \"INFO/Interpro_domain\" NOT NULL \
>             ORDER BY \"INFO/SiPhy_29way_logOdds_rankscore\" DESC"
> ```

#### CSV format

Use duckDB function `read_csv_auto` to read a TSV file format as a table. See [duckDB CSV import](https://duckdb.org/docs/data/csv/overview.html) for more information.

> Example: Query a TSV file
> ```
> howard query \
>    --query="SELECT * FROM read_csv_auto('tests/data/transcripts.tsv')"
> ```

> Example: Query a TSV file with columns struct
> ```
> howard query \
>    --query="SELECT * FROM read_csv_auto('tests/data/transcripts.tsv', columns={'transcript': 'VARCHAR','gene': 'VARCHAR'})"
> ```

#### BED format

In order to read a BED file, create a query (using appropiate columns), and export file into a desired format.

> Example: Read a BED file and export in Parquet format
> ```
> howard query \
>    --query="SELECT * FROM read_csv_auto('tests/data/example.bed', columns={'#CHROM': 'VARCHAR', 'START': 'INTEGER', 'END': 'INTEGER'})" \
>    --output=/tmp/example.bed.parquet
> ```

> Example: Convert a BED file in a Parquet format into a BED file format
> ```
> howard convert \
>    --input=tests/data/example.bed.parquet \
>    --output=/tmp/example.bed
> 
> cat /tmp/example.bed
> ```
> ```
> #CHROM	START	END
> chr1	28735	69101
> chr1	768250	768253
> chr7	55249060	55249069
> ```

A BED file can be filtered using positions or other columns such as gene names, or transcripts.

> Example: Filter a BED using positions
> ```
> howard query \
>    --input=tests/data/example.bed.parquet \
>    --query="SELECT * \
>             FROM variants \
>             WHERE \
>                \"#CHROM\" = 'chr1' and \
>                ( \
>                   (\"START\">28000 and \"END\"<70000) or \
>                   (\"START\">760000 and \"END\"<770000) \
>                )"
> ```
> ```
>   #CHROM   START     END
> 0   chr1   28735   69101
> 1   chr1  768250  768253
> ```

### Extract variants

In order to extract variants from a VCF file, without annotations and samples, use a query to construct the VCF.

> Example: Extract variants onnly
> ```
> howard query \
>    --input=tests/data/example.vcf.gz \
>    --output=/tmp/example.vcf.gz \
>    --query="SELECT \"#CHROM\", POS, ID, REF, ALT, QUAL, FILTER, '.' AS INFO \
>             FROM variants"

## Annotation

Annotation is mainly based on a build-in Parquet annotation method, using annotation database file (in multiple format such as Parquet, duckdb, VCF, BED, TSV, JSON).

See [HOWARD Help Annotation tool](help.md#annotation-tool) for more information.

These annotation databases can be automatically downloaded using [HOWARD Databases tool](help.md#databases-tool) and manually generated using existing annotation files and [HOWARD Convert tool](help.md#convert-tool). Annotation databases need a header file (`.hdr`) to describe each annotation in the database. However, a default header will be generated if no header file is associated to the annotation database file.

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
> cd /tool
> howard annotation \
>    --input=tests/data/example.vcf.gz \
>    --annotations='tests/databases/annotations/current/hg19/cosmic70.vcf.gz' \
>    --output=/tmp/example.howard.vcf.gz
> ```

> Example: VCF annotation with relative path BED databases
> ```
> cd /tool
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

External annotation tools are also available, such as BCFTOOLS, Annovar, snpEff and Exomiser. Annovar, snpEff and Exomiser databases are automatically downloaded (see [HOWARD Help Databases tool](help.md#databases-tool)). Quick annotation allows to annotates by simply defining external tools keywords.

##### BCFTools annotation

For BCFTools, use HOWARD keyword `bcftools` and list (separator `:` or `+`) annotation databases with format such as VCF or BED (compressed). More options are available using [HOWARD Parameters JSON](help.param.md) file.

> Example: VCF annotation with Cosmic VCF databases and refGene BED database
> ```
> howard annotation \
>    --input=tests/data/example.vcf.gz \
>    --annotations='bcftools:tests/databases/annotations/current/hg19/cosmic70.vcf.gz+tests/databases/annotations/current/hg19/refGene.bed.gz' \
>    --output=/tmp/example.howard.vcf.gz
> ```

##### Annovar annotation

For Annovar tool, use HOWARD keyword `annovar` and mention specific Annovar database keywords (separator `:`). More options are available using [HOWARD Parameters JSON](help.param.md) file.

> Example: VCF annotation with Annovar refGene and cosmic70
> ```
> howard annotation \
>    --input=tests/data/example.vcf.gz \
>    --annotations='annovar:refGene:cosmic70' \
>    --output=/tmp/example.howard.tsv
> ```

##### snpEff annotation

For snpEff tool, use HOWARD keyword `snpeff`. No options are available for quick annotation with snpEff, see [HOWARD Parameters JSON - snpEff](help.param.md#snpeff) for more options.

> Example: VCF annotation with snpEff 
> ```
> howard annotation \
>    --input=tests/data/example.vcf.gz \
>    --annotations='snpeff' \
>    --output=/tmp/example.howard.tsv
> ```

##### Exomiser Annotation

For Exomiser tool, use HOWARD keyword `exomiser`. A list of options can be provided as key-value format, such as exomiser release, a preset (pre-configured options), source of transcripts (e.g. 'refseq', 'ucsc'), a llist of HPO terms (do not use ':' separator, e.g. '0001156', 'HP0001156', 'hpo0001156'). More options are available using [HOWARD Parameters JSON](help.param.md) file.

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
>    --annotations='tests/databases/annotations/current/hg19/dbnsfp42a.parquet,bcftools:tests/databases/annotations/current/hg19/cosmic70.vcf.gz,annovar:refGene:cosmic70,snpeff,exomiser:preset=exome:hpo=0001156+0001363+0011304+0010055' \
>    --output=/tmp/example.howard.tsv
> ```

See [HOWARD Help Annotation tool](help.md#annotation-tool) tool for more information.

### Annotation parameters

All annotation parameters can be defined in [HOWARD Parameters JSON](help.param.md) file. All annotations can be combined (bild-in parquet method and external tools annotation), and options can be detailed in a full JSON file format, including selection of annotation database columns (and rename them) and specific options for external tools.

## Calculation

Calculation processes variants annotations to generate new annotation, such as: identify variation type (VarType), harmonizes allele frequency (VAF) and calculate sttistics (VAF_stats), extracts Nomen (transcript, cNomen, pNomen...) from an HGVS field (e.g. snpEff, Annovar) with an optional list of personalized transcripts, generates VaRank format barcode, identify trio inheritance. These calculations are based on existing annotations of variants (and genotypes).

Calculations are either provided by HOWARD within code, or configured into a JSON file. Calculations are either an inner HOWARD Python code, or a SQL query.

See [HOWARD Help Calculation tool](help.md#calculation-tool) tool for more information.

To process a calculation, use is keyword with the `--calculations` parameter.

> Example: calculation of the variant type with `vartype` keyword
> ```
> howard calculation \
>    --input=tests/data/example.full.vcf \
>    --calculations='vartype' \
>    --output=/tmp/example.calculation.tsv
> ```

### Available calculations

To list all available calculations, from HOWARD default configuration or with a homemade Calculation configuration JSON file, use the `--show_calculations` parameter.

> Example: List of build-in calculation
> ```
> howard calculation \
>    --show_calculations
> ```
> ```
> #[2024-03-05 10:53:30] [INFO] Start
> #[2024-03-05 10:53:30] [INFO] Available calculation operations:
> #[2024-03-05 10:53:30] [INFO]    BARCODE: BARCODE as VaRank tool
> #[2024-03-05 10:53:30] [INFO]    DP_STATS: Depth (DP) statistics
> #[2024-03-05 10:53:30] [INFO]    FINDBYPIPELINE: Number of pipeline that identify the variant (for multi pipeline VCF)
> #[2024-03-05 10:53:30] [INFO]    FINDBYSAMPLE: Number of sample that have a genotype for the variant (for multi sample VCF)
> #[2024-03-05 10:53:30] [INFO]    GENOTYPECONCORDANCE: Concordance of genotype for multi caller VCF
> #[2024-03-05 10:53:30] [INFO]    NOMEN: NOMEN information (e.g. NOMEN, CNOMEN, PNOMEN...) from HGVS nomenclature field
> #[2024-03-05 10:53:30] [INFO]    SNPEFF_HGVS: HGVS nomenclatures from snpEff annotation
> #[2024-03-05 10:53:30] [INFO]    TRIO: Inheritance for a trio family
> #[2024-03-05 10:53:30] [INFO]    VAF: Variant Allele Frequency (VAF) harmonization
> #[2024-03-05 10:53:30] [INFO]    VAF_STATS: Variant Allele Frequency (VAF) statistics
> #[2024-03-05 10:53:30] [INFO]    VARIANT_ID: Variant ID generated from variant position and type
> #[2024-03-05 10:53:30] [INFO]    VARTYPE: Variant type (e.g. SNV, INDEL, MNV, BND...)
> #[2024-03-05 10:53:30] [INFO] End
> ```

### Calculation configuration JSON file

All calculations are configured in a JSON file. A default configuration is provided with default calculations.

Basically, a calculation is defined by:
- Type: either 'sql' for a SQL query or 'python' for a Python function
- Name/keyword: a keyword that is used with `--show_calculations` parameter (case unsensitive)
- Description: a description of the calculation
- Output column information: Name, type and decription of the new annotation calculated
- Query and fields: an SQL query (for 'sql' type) with parameters such as mandatory INFO fields
- Function name and parameters: a existing Python function and parameters (for 'python' type)

> Example: Calculation of variant type using an SQL query
> ```
> "VARTYPE": {
>     "type": "sql",
>     "name": "VARTYPE",
>     "description": "Variant type (e.g. SNV, INDEL, MNV, BND...)",
>     "available": true,
>     "output_column_name": "VARTYPE",
>     "output_column_type": "String",
>     "output_column_description": "Variant type: SNV if X>Y, MOSAIC if X>Y,Z or X,Y>Z, INDEL if XY>Z or X>YZ",
>     "operation_query": [
>       "CASE",
>       "WHEN \"SVTYPE\" NOT NULL THEN \"SVTYPE\"",
>       "WHEN LENGTH(REF) = 1 AND LENGTH(ALT) = 1 THEN 'SNV'",
>       "WHEN REF LIKE '%,%' OR ALT LIKE '%,%' THEN 'MOSAIC'",
>       "WHEN LENGTH(REF) == LENGTH(ALT) AND LENGTH(REF) > 1 THEN 'MNV'",
>       "WHEN LENGTH(REF) <> LENGTH(ALT) THEN 'INDEL'",
>       "ELSE 'UNDEFINED'",
>       "END"
>     ],
>     "info_fields": ["SVTYPE"],
>     "operation_info": true
>   }
> ```

> Example: Calculation of variant id using an existing Python function `calculation_variant_id`
> ```
> "variant_id": {
>     "type": "python",
>     "name": "variant_id",
>     "description": "Variant ID generated from variant position and type",
>     "available": true,
>     "function_name": "calculation_variant_id",
>     "function_params": []
>   }
> ```

See [Calculation configuration JSON file example](../config/calculations_config.json).

See  [Calculation configuration JSON file](calculations.config.md) for more information.

### Build-in calculations examples

#### Variant type

Variant type calculation `vartype` detect the type of variant (e.g. SNV, INDEL, MNV). Variant type are claculated with these criteria: SNV if X>Y, MOSAIC if X>Y,Z or X,Y>Z, INDEL if XY>Z or X>YZ.

> Example: Identify variant types and generate a table of variant type count
> ```
> howard calculation \
>    --input=tests/data/example.full.vcf \
>    --calculations='vartype' \
>    --output=/tmp/example.calculation.tsv
> 
> howard query \
>    --input=/tmp/example.calculation.tsv \
>    --explode_infos \
>    --query='SELECT
>                "VARTYPE" AS 'VariantType',
>                count(*) AS 'Count'
>             FROM variants
>             GROUP BY "VARTYPE"
>             ORDER BY count DESC'
> ```
> ```
>   VariantType  Count
> 0         BND      7
> 1         DUP      6
> 2         INS      5
> 3         SNV      4
> 4         CNV      3
> 5         DEL      3
> 6         INV      3
> 7      MOSAIC      2
> 8       INDEL      2
> 9         MNV      1
> ```

#### HGVS and NOMEN from snpEff

NOMEN can be extracted from snpEff annotation (see [HOWARD Parameters JSON - snpEff](help.param.md#snpeff)). The first calculation extract list of HGVS annotations from snpEff annotation (`snpeff_hgvs` keyword), the second calculation choose the NOMEN from snpEff HGVS annotations using a list of reference transcripts (`NOMEN` keyword, `--hgvs_field` and `--transcripts` parameters). More options are available (see [HOWARD Parameters JSON](help.param.md)).

> Example: Calculate NOMEN by extracting hgvs from snpEff annotation and identifying transcripts from a list
>
> ```
> howard calculation \
>    --input=tests/data/example.ann.vcf.gz \
>    --calculations='snpeff_hgvs,NOMEN' \
>    --hgvs_field='snpeff_hgvs' \
>    --transcripts=tests/data/transcripts.tsv \
>    --output=/tmp/example.NOMEN.vcf.gz
> 
> howard query \
>    --input=/tmp/example.NOMEN.vcf.gz \
>    --explode_infos \
>    --query="SELECT GNOMEN, NOMEN, snpeff_hgvs \
>             FROM variants \
>             WHERE GNOMEN='EGFR'"
> ```
> ```
>   GNOMEN  NOMEN                                           snpeff_hgvs
> 0   EGFR  EGFR:NM_001346897:exon19:c.2226G>A:p.Gln742Gln  EGFR:NM_005228.5:exon20:c.2361G>A:p.Gln787Gln,...
> ```

## Prioritization

See [HOWARD Help Prioritization tool](help.md#prioritization-tool) tool for more information.

## Process

See [HOWARD Help Process tool](help.md#process-tool) tool for more information (under development).