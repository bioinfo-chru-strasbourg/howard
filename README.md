HOWARD
===

HOWARD annotates and prioritizes genetic variations, calculates and normalizes annotations, translates files in multiple formats (e.g. vcf, tsv, parquet) and generates variants statistics.

HOWARD annotation is mainly based on a build-in Parquet annotation method, and tools such as BCFTOOLS, ANNOVAR and snpEff. It using available databases (see ANNOVAR and snpEff) and homemade databases. Format of databases are: parquet, duckdb, vcf, bed, ANNOVAR and snpEff (ANNOVAR and snpEff databases are automatically downloaded if needed). 

HOWARD calculation processes variants information to calculate new information, such as: harmonizes allele frequency (VAF), extracts Nomen (transcript, cNomen, pNomen...) from HGVS fields with an optional list of personalized transcripts, generates VaRank format barcode.

HOWARD prioritization algorithm uses profiles to flag variants (as passed or filtered), calculate a prioritization score, and automatically generate a comment for each variants (example: 'polymorphism identified in dbSNP. associated to Lung Cancer. Found in ClinVar database').Prioritization profiles are defined in a configuration file. A profile is defined as a list of annotation/value, using wildcards and comparison options (contains, lower than, greater than, equal...). Annotations fields may be quality values (usually from callers, such as 'GQ', 'DP') or other annotations fields provided by annotations tools, such as HOWARD itself (example: COSMIC, Clinvar, 1000genomes, PolyPhen, SIFT). Multiple profiles can be used simultaneously, which is useful to define multiple validation/prioritization levels (example: 'standard', 'stringent', 'rare variants', 'low allele frequency').

HOWARD translates VCF format into multiple formats (e.g. VCXF, TSV, Parquet), by sorting variants using specific fields (example : 'prioritization score', 'allele frequency', 'gene symbol'), including/excluding annotations/fields, including/excluding variants, adding fixed columns.

HOWARD generates statistics files with a specific algorithm, snpEff and BCFTOOLS.

HOWARD is multithreaded through the number of variants and by database (data-scaling).


# Getting Started

## Installation

### Python

Install using Python Pip tool:
```
python -m pip install -e .
```

### Docker

In order to build, setup and create a persitent CLI (running container), docker-compose command build images and launch services as containers.

```
docker-compose up -d
```

A setup container (HOWARD-setup) automatically downloads required databases according to an HOWARD VCF example annotation using ANNOVAR and snpEff. Configuration of host data and databases folders (default ${HOME}/HOWARD), assembly and databases to download in `.env` file. See HOWARD, ANNOVAR and snpEff documentation for custom databases download.

A Command Line Interface container (HOWARD-CLI) is started with host data and databases folders mounted.

## Quick HOWARD commands

### Show VCF stats and overview

Show example VCF statistics and brief overview
```
howard --input=tests/data/example.vcf.gz --stats --overview
```

### Translate VCF into other format

Translate VCF into CSV and show output file
```
howard --input=tests/data/example.vcf.gz --output=tests/data/example.csv && cat tests/data/example.csv
```

Translate VCF into parquet, and show statistics on output file
```
howard --input=tests/data/example.vcf.gz --output=tests/data/example.parquet && howard --input=tests/data/example.parquet --stats
```

### Query VCF

Select variants in VCF with REF and POS fields
```
howard --input=tests/data/example.vcf.gz --query="SELECT * FROM variants WHERE REF = 'A' AND POS < 100000"
```

Select variants in VCF with INFO Tags criterions
```
howard --input=tests/data/example.vcf.gz --param='{"explode_infos": true}' --query="SELECT \"#CHROM\", POS, REF, ALT, \"INFO/DP\", \"INFO/CLNSIG\", sample2, sample3 FROM variants WHERE \"INFO/DP\" >=50 OR \"INFO/CLNSIG\" = 'pathogenic'"
```

### Annotation

VCF annotation with Parquet databases, output as VCF format, and show INFO
```
howard --verbose --input=tests/data/example.vcf.gz --output=/tmp/example.howard.vcf.gz --annotation=tests/data/annotations/avsnp150.parquet,tests/data/annotations/dbnsfp42a.parquet,tests/data/annotations/gnomad211_genome.parquet --query='SELECT INFO FROM variants'
```

VCF annotation with Clinvar Parquet databases, output as TSV format and INFO tags as column, and query Clinvar results
```
howard --verbose --input=tests/data/example.vcf.gz --output=/tmp/example.howard.tsv --annotation=tests/data/annotations/clinvar_20210123.parquet --param='{"explode_infos": true, "export_extra_infos": true}' --query='SELECT "INFO/CLNDN", count(*) AS count FROM variants WHERE "INFO/CLNDN" NOT NULL GROUP BY "INFO/CLNDN"'
```

VCF annotation (Parquet, BCFTOOLS, ANNOVAR and snpEff) using HOWARD-CLI (snpEff and ANNOVAR databases will be automatically downloaded), and query list of genes with HGVS

```
docker exec HOWARD-CLI howard --verbose --input=/tool/tests/data/example.vcf.gz --output=/data/example.howard.vcf.gz --annotation=snpeff,annovar:refGene,/tool/tests/data/annotations/refGene.bed.gz,/tool/tests/data/annotations/avsnp150.vcf.gz,tests/data/annotations/dbnsfp42a.parquet --param='{"explode_infos": true}' --query='SELECT "INFO/symbol", "INFO/AAChange_refGene" FROM variants WHERE "INFO/symbol" NOT NULL ORDER BY "INFO/symbol"'
```

Let's play within Docker HOWARD-CLI service!
```
$ docker exec -ti HOWARD-CLI bash
[data]# howard --help
```

# Documentation

More documentation in [docs/howard.md](docs/howard.md)



# Docker 

HOWARD image presents a container that runs on AlmaLinux, and includes yum modules and other tools dependencies (see config.cfg and Dockerfile for releases):
- Python
- Java
- bcftools/htslib
- ANNOVAR
- snpEff

## Docker Build - Image


The `Dockerfile` provided with this package provides everything that is needed to build the image. The build system must have Docker installed in
order to build the image.

```
$ cd ${HOME}/HOWARD
$ docker build -t howard:latest .
```

## Running Run - Container

The container host must have Docker installed in order to run the image as a container. Then the image can be pulled and a container can be started directly. Any standard Docker switches may be provided on the command line when running a container.

```
$ docker run howard:latest
```


### Mount Data and Databases volumes

In order to make data and databases persistent, host volumes can be mounted. Content may also be copied directly into the running container using a
`docker cp ...`.

```
-v ${HOME}/HOWARD/data:/data
-v ${HOME}/HOWARD/databases:/databases
```

### Run as a terminal

In order to execute command directly to an container, start HOWARD container with terminal interface:

```
$ docker run --name howard --entrypoint=bash -ti howard:latest
```

### Example

Run HOWARD as a uniq command.

```
$ docker run --rm -v ${HOME}/HOWARD/data:/data -v ${HOME}/HOWARD/databases:/databases howard:latest --input=/tool/docs/example.vcf --output=/data/example.howard.tsv --annotation=snpeff,annovar:refGene,/databases/gnomad211_exome.parquet,/databases/cosmic70.vcf.gz --calculation=VAF,BARCODE,NOMEN --prioritization=GERMLINE
```

Database download (TODO)
-------------------

Databases are downloaded automatically by using annotation configuratin file, or options in command line (--annovar_databases, --snpeff_databases, assembly...).

Use a vcf file, such as HOWARD VCF example, to download ANNOVAR and snpEff databases (WITHOUT multithreading, "ALL" for all databases, "core" for core databases, "snpeff" for snpEff database, or a list of databases, or ANNOVAR code). Use this command multiple times for all needed databases and assembly (such as hg19, hg38, mm9).

```
$ docker run howard:latest --input=/tool/docs/example.vcf --output=/tool/docs/example.annotated.vcf --annotation=ALL,snpeff --thread=1 --verbose
```


> Note: For home made databases, refer to ```config.annotation.ini``` file to construct and configure your own database (deprecated).

> Note: Beware of proxy configuration!


