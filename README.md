# HOWARD

![HOWARD Graphical User Interface](images/icon.ico "HOWARD Graphical User Interface")

Highly Open and Valuable tool for Variant Annotation & Ranking toward genetic Discovery

HOWARD annotates and prioritizes genetic variations, calculates and normalizes annotations, translates files in multiple formats (e.g. vcf, tsv, parquet) and generates variants statistics.

HOWARD annotation is mainly based on a build-in Parquet annotation method, and external tools such as BCFTOOLS, ANNOVAR, snpEff and Exomiser (see docs, automatically downloaded if needed). Parquet annotation uses annotation database in VCF or BED format, in mutliple file format: Parquet/duckdb, VCF, BED, TSV, CSV, TBL, JSON.

HOWARD calculation processes variants information to calculate new information, such as: harmonizes allele frequency (VAF), extracts Nomen (transcript, cNomen, pNomen...) from HGVS fields with an optional list of personalized transcripts, generates VaRank format barcode.

HOWARD prioritization algorithm uses profiles to flag variants (as passed or filtered), calculate a prioritization score, and automatically generate a comment for each variants (example: 'polymorphism identified in dbSNP. associated to Lung Cancer. Found in ClinVar database').Prioritization profiles are defined in a configuration file. A profile is defined as a list of annotation/value, using wildcards and comparison options (contains, lower than, greater than, equal...). Annotations fields may be quality values (usually from callers, such as 'GQ', 'DP') or other annotations fields provided by annotations tools, such as HOWARD itself (example: COSMIC, Clinvar, 1000genomes, PolyPhen, SIFT). Multiple profiles can be used simultaneously, which is useful to define multiple validation/prioritization levels (example: 'standard', 'stringent', 'rare variants', 'low allele frequency').

HOWARD translates VCF format into multiple formats (e.g. VCF, TSV, Parquet), by sorting variants using specific fields (example : 'prioritization score', 'allele frequency', 'gene symbol'), including/excluding annotations/fields, including/excluding variants, adding fixed columns.

HOWARD generates statistics files with a specific algorithm, snpEff and BCFTOOLS.

HOWARD is multithreaded through the number of variants and by database (data-scaling).

## Table of contents

- [Installation](#installation)
  - [Python](#python)
  - [Docker](#docker)
  - [Databases](#databases)
- [Configuration](#configuration)
- [Parameters](#parameters)
- [Tools](#tools)
  - [Stats](#stats)
  - [Convert](#convert)
  - [Query](#query)
  - [Annotation](#annotation)
  - [Calculation](#calculation)
  - [Prioritization](#prioritization)
  - [Process](#process)
- [Docker HOWARD-CLI](#docker-howard-cli)
- [Documentation](#documentation)
- [Contact](#contact)

# Installation

## Python

Install HOWARD using Python Pip tool, and run HOWARD for help options:
```
python -m pip install -e .
howard --help
```

Install HOWARD Graphical User Interface using Python Pip tool with supplementary packages, and run as a tool:
```
python -m pip install -r requirements-gui.txt
howard gui
```

![HOWARD Graphical User Interface](images/howard-gui.png "HOWARD Graphical User Interface")

## Docker

In order to build, setup and create a persitent CLI (running container with all useful external tools such as [BCFTools](https://samtools.github.io/bcftools/), [snpEff](https://pcingola.github.io/SnpEff/), [Annovar](https://annovar.openbioinformatics.org/), [Exomiser](https://www.sanger.ac.uk/tool/exomiser/)), docker-compose command build images and launch services as containers.

```
docker-compose up -d
```

A setup container (HOWARD-setup) will download useful databases (take a while). To avoid databases download (see [Databases section](#databases) to download manually), just start:

```
docker-compose up -d HOWARD-CLI
```

A Command Line Interface container (HOWARD-CLI) is started with host data and databases folders mounted (by default in ${HOME}/HOWARD folder). Let's play within Docker HOWARD-CLI service!

```
docker exec -ti HOWARD-CLI bash
howard --help
```

## Databases

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

- Download Multiple databases in the same time for assembly 'hg19' (can take a while):

```
howard databases --assembly=hg19 --download-genomes=~/howard/databases/genomes/current --download-genomes-provider=UCSC --download-genomes-contig-regex='chr[0-9XYM]+$' --download-annovar=~/howard/databases/annovar/current --download-annovar-files='refGene,cosmic70,nci60' --download-snpeff=~/howard/databases/snpeff/current --download-refseq=~/howard/databases/refseq/current --download-refseq-format-file='ncbiRefSeq.txt' --download-dbnsfp=~/howard/databases/dbnsfp/current --download-dbnsfp-release='4.4a' --download-dbnsfp-subdatabases --download-alphamissense=~/howard/databases/alphamissense/current --download-exomiser=~/howard/databases/exomiser/current --download-dbsnp=~/howard/databases/dbsnp/current --download-dbsnp-vcf --threads=8
```

See [HOWARD Help Databases tool](docs/help.md#databases-tool) for more information.

Databases can be home-made generated, starting with a existing annotation file, especially using [HOWARD convert](#convert) tool. These files need to contain specific fields (depending on the annotation type):

- variant annotation: '#CHROM', 'POS', 'ALT', 'REF'
- region annotation: '#CHROM', 'START', 'STOP'

Each database annotation file is associated with a 'header' file ('.hdr'), in VCF header format, to describe annotations within the database.

# Configuration

HOWARD Configuration JSON file defined default configuration regarding resources (e.g. threads, memory), settings (e.g. verbosity, temporary files), default folders (e.g. for databases) and paths to external tools.

See [HOWARD Configuration JSON](docs/help.config.md) for more information.

# Parameters

HOWARD Parameters JSON file defined parameters to process annotations, prioritization, calculations, convertions and queries.

See [HOWARD Parameters JSON](docs/help.param.md) for more information.

# Tools

## Stats

Statistics on genetic variations, such as: number of variants, number of samples, statistics by chromosome, genotypes by samples, annotations...
Theses statsitics can be applied to VCF files and all database annotation files.

> Example: Show example VCF statistics and brief overview
> ```
> howard stats \
>    --input=tests/data/example.vcf.gz
> ```

See [HOWARD Help Stats tool](docs/help.md#stats-tool) for more information.

## Convert

Convert genetic variations file to another format. Multiple format are available, such as usual and official VCF format, but also other formats such as TSV, CSV, TBL, JSON and Parquet/duckDB. These formats need a header '.hdr' file to take advantage of the power of howard (especially through INFO/tag definition), and using howard convert tool automatically generate header file fo futher use (otherwise, an default '.hdr' file is generated).

> Example: Translate VCF into TSV, export INFO/tags into columns, and show output file
> ```
> howard convert \
>    --input=tests/data/example.vcf.gz \
>    --explode_infos \
>    --output=/tmp/example.tsv
> cat /tmp/example.tsv
> ```

See [HOWARD Help Convert tool](docs/help.md#convert-tool) for more options.

## Query

Query genetic variations in SQL format. Data can be loaded into 'variants' table from various formats (e.g. VCF, TSV, Parquet...). Using --explode_infos allow query on INFO/tag annotations. SQL query can also use external data within the request, such as a Parquet file(s).

> Example: Select variants in VCF with INFO Tags criterions
> ```
> howard query \
>    --input=tests/data/example.vcf.gz \
>    --explode_infos \
>    --query='SELECT "#CHROM", POS, REF, ALT, DP, CLNSIG, sample2, sample3 
>             FROM variants 
>             WHERE DP >= 50 OR CLNSIG NOT NULL 
>             ORDER BY CLNSIG DESC, DP DESC'
> ```

See [HOWARD Help Query tool](docs/help.md#query-tool) for more options.

## Annotation

Annotation is mainly based on a build-in Parquet annotation method, using database format such as Parquet, duckdb, VCF, BED, TSV, JSON. External annotation tools are also available, such as BCFTOOLS, Annovar, snpEff and Exomiser. It uses available databases and homemade databases. Annovar and snpEff databases are automatically downloaded (see [HOWARD Help Databases tool](docs/help.md#databases-tool)). All annotation parameters are defined in [HOWARD Parameters JSON](docs/help.param.md) file.

Quick annotation allows to annotates by simply listing annotation databases, or defining external tools keywords. These annotations can be combined.

> Example: VCF annotation with Parquet and VCF databases, output as VCF format
> ```
> howard annotation \
>    --input=tests/data/example.vcf.gz \
>    --annotations='tests/databases/annotations/current/hg19/dbnsfp42a.parquet,tests/databases/annotations/current/hg19/cosmic70.vcf.gz' \
>    --output=/tmp/example.howard.vcf.gz
> ```

> Example: VCF annotation with external tools (Annovar refGene and snpEff databases), output as TSV format
> ```
> howard annotation \
>    --input=tests/data/example.vcf.gz \
>    --annotations='annovar:refGene,snpeff' \
>    --output=/tmp/example.howard.tsv
> ```

See [HOWARD Help Annotation tool](docs/help.md#annotation-tool) for more options.

## Calculation

Calculation processes variants information to generate new information, such as: identify variation type (VarType), harmonizes allele frequency (VAF) and calculate sttistics (VAF_stats), extracts Nomen (transcript, cNomen, pNomen...) from an HGVS field (e.g. snpEff, Annovar) with an optional list of personalized transcripts, generates VaRank format barcode, identify trio inheritance.

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


See [HOWARD Help Calculation tool](docs/help.md#calculation-tool) for more options.

## Prioritization

Prioritization algorithm uses profiles to flag variants (as passed or filtered), calculate a prioritization score, and automatically generate a comment for each variants (example: 'polymorphism identified in dbSNP. associated to Lung Cancer. Found in ClinVar database'). Prioritization profiles are defined in a configuration file in JSON format. A profile is defined as a list of annotation/value, using wildcards and comparison options (contains, lower than, greater than, equal...). Annotations fields may be quality values (usually from callers, such as 'DP') or other annotations fields provided by annotations tools, such as HOWARD itself (example: COSMIC, Clinvar, 1000genomes, PolyPhen, SIFT). Multiple profiles can be used simultaneously, which is useful to define multiple validation/prioritization levels (example: 'standard', 'stringent', 'rare variants', 'low allele frequency').

- Prioritize variants from criteria on INFO annotations for profiles 'default' and 'GERMLINE' (see 'prioritization_profiles.json'), and query variants on prioritization tags

```
howard prioritization --input=tests/data/example.vcf.gz --prioritizations=config/prioritization_profiles.json --profiles='default,GERMLINE' --pzfields='PZFlag,PZScore,PZComment' --output=/tmp/example.prioritized.vcf.gz
```

- and query variants passing filters

```
howard query --input=/tmp/example.prioritized.vcf.gz --explode_infos --query="SELECT \"#CHROM\", POS, ALT, REF, \"PZFlag\", \"PZScore\", \"DP\", \"CLNSIG\" FROM variants WHERE \"PZScore\" > 0 AND \"PZFlag\" == 'PASS' ORDER BY \"PZScore\" DESC"
```

- and query variants with different prioritization flag between profiles

```
howard query --input=/tmp/example.prioritized.vcf.gz --explode_infos --query="SELECT \"#CHROM\", POS, ALT, REF, \"PZFlag_default\", \"PZFlag_GERMLINE\" FROM variants WHERE \"PZFlag_default\" != \"PZFlag_GERMLINE\" ORDER BY \"PZScore\" DESC"
```

- and showing propritization comments of variants, with flags and scores

```
howard query --input=/tmp/example.prioritized.vcf.gz --explode_infos --query="SELECT \"#CHROM\", POS, ALT, REF, \"PZComment\", \"PZFlag\", \"PZScore\" FROM variants WHERE \"PZComment\" IS NOT NULL ORDER BY \"PZScore\" DESC"
```

See [HOWARD Help Prioritization tool](docs/help.md#prioritization-tool) for more options.

## Process

howard process tool manage genetic variations to:

- annotates genetic variants with multiple annotation databases/files and tools
- calculates and normalizes annotations
- prioritizes variants with profiles (list of citeria) to calculate scores and flags
- translates into various formats
- query genetic variants and annotations
- generates variants statistics

This process tool combines all other tools to pipe them in a uniq command, through a parameter file in JSON format.

- Full process command

```
howard process --config=config/config.json --param=config/param.json --input=tests/data/example.vcf.gz --output=/tmp/example.process.tsv --query='SELECT "NOMEN", "PZFlag", "PZScore", "PZComment" FROM variants ORDER BY "PZScore" DESC' --explode_infos
```

- to obtain this kind of table

```
                                            NOMEN    PZFlag  PZScore                                        PZComment
0              WASH7P:NR_024540:exon1:n.50+585T>G      PASS       15      Described on CLINVAR database as pathogenic
1      OR4F5:NM_001005484:exon1:c.11A>G:p.Glu4Gly      PASS        5                                DP higher than 50
2  EGFR:NM_001346897:exon19:c.2226G>A:p.Gln742Gln      PASS        5                                DP higher than 50
3         LINC01128:NR_047519:exon2:n.287+3767A>G      PASS        0                                              NaN
4         LINC01128:NR_047519:exon2:n.287+3768A>G      PASS        0                                              NaN
5         LINC01128:NR_047519:exon2:n.287+3769A>G      PASS        0                                              NaN
6                  MIR1302-9:NR_036266:n.*4641A>C  FILTERED     -100  Described on CLINVAR database as non-pathogenic
```

- with parameter JSON file

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

See [HOWARD Help Process tool](docs/help.md#process-tool) for more options.

# Command Line Interface

Docker HOWARD-CLI container (Command Line Interface) can be used to execute commands.

- Query of an existing VCF:

```
docker exec HOWARD-CLI howard query --input=/tool/tests/data/example.vcf.gz --query='SELECT * FROM variants'
```

- VCF annotation using HOWARD-CLI (snpEff and ANNOVAR databases will be automatically downloaded), and query list of genes with HGVS:

```
docker exec --workdir=/tool HOWARD-CLI howard process --config=config/config.json --param=config/param.json --input=tests/data/example.vcf.gz --output=/tmp/example.process.tsv --query='SELECT "NOMEN", "PZFlag", "PZScore", "PZComment" FROM variants ORDER BY "PZScore" DESC' --explode_infos
```

- Let's play within Docker HOWARD-CLI service!

```
docker exec -ti HOWARD-CLI bash
howard --help
```

# Documentation

[HOWARD User Guide](docs/user_guide.md) is available to assist users for particular commands, such as software installation, databases download, annotation command, and so on.

HOWARD provides multiple tools to manage variants and databases, with many options. [HOWARD Help](docs/options.md) is availalbe in markdown format, and compiles all information also available in help command (for each tool):

```
howard --help               # Help to list all tools
howard query --help         # Help for 'query' tool
howard annotation --help    # Help for 'annotation' tool
...
```

# Contact

:hospital: [Medical Bioinformatics applied to Diagnosis Lab](https://www.chru-strasbourg.fr/service/bioinformatique-medicale-appliquee-au-diagnostic-unite-de/) @ Strasbourg Univerty Hospital

:email: [bioinfo@chru-strasbourg.fr](bioinfo@chru-strasbourg.fr)

:bookmark: [GitHub](https://github.com/bioinfo-chru-strasbourg)