HOWARD
===

HOWARD annotates and prioritizes genetic variations, calculates and normalizes annotations, translates vcf format and generates variants statistics.

HOWARD annotation is mainly based on ANNOVAR and snpEff tools to annotate, using available databases (see ANNOVAR and snpEff) and home made databases. It also uses BCFTOOLS to annotate variants with a VCF file. ANNOVAR and snpEff databases are automatically downloaded if needed.

HOWARD calculation harmonizes allele frequency (VAF), extracts Nomen (transcript, cNomen, pNomen...) from HGVS fields with an optional list of personalized transcripts, generates VaRank format barcode.

HOWARD prioritization algorithm uses profiles to flag variants (as passed or filtered), calculate a prioritization score, and automatically generate a comment for each variants (example: 'polymorphism identified in dbSNP. associated to Lung Cancer. Found in ClinVar database').Prioritization profiles are defined in a configuration file. A profile is defined as a list of annotation/value, using wildcards and comparison options (contains, lower than, greater than, equal...). Annotations fields may be quality values (usually from callers, such as 'GQ', 'DP') or other annotations fields provided by annotations tools, such as HOWARD itself (example: COSMIC, Clinvar, 1000genomes, PolyPhen, SIFT). Multiple profiles can be used simultaneously, which is useful to define multiple validation/prioritization levels (example: 'standard', 'stringent', 'rare variants', 'low allele frequency').

HOWARD translates VCF format into TSV format, by sorting variants using specific fields (example : 'prioritization score', 'allele frequency', 'gene symbol'), including/excluding annotations/fields, including/excluding variants, adding fixed columns.

HOWARD generates statistics files with a specific algorithm, snpEff and BCFTOOLS.

HOWARD is multithreaded through the number of variants and by database (data-scaling).


# Getting Started


In order to build, setup and create a persitent CLI (running container), docker-compose command build images and launch services as containers.

```
$ docker-compose up
```

A setup container (HOWARD-setup) automatically downloads required databases according to an HOWARD VCF example annotation using ANNOVAR and snpEff. Configuration of host data and databases folders (default ${HOME}/HOWARD), assembly and databases to download in `.env` file. See HOWARD, ANNOVAR and snpEff documentation for custom databases download.

A Command Line Interface container (HOWARD-CLI) is started with host data and databases folders mounted. Execute a command, or connect to the CLI as a terminal, and let's start with HOWARD!

Using an HOWARD VCF example, this command:
- 1- annotates with HGVS (variation identification), outcome and location (fonctionnal annotation), and clinical databases (ClinVar and Cosmic),
- 2- calculates the Variant Allele Frquency (VAF), a genotype barcode (BARCODE), and process HGVS to extract NOMEN information,
- 3- prioritizes variations according to priorization rules specific to somatic focus (quality, functionnal and clinical annotations),
- 4- translates into TSV format, with specific fields order for the first 3 columns (ALL for the rest), and a sorting to focus on intersting variations (Flag as PASS, with best score) 
- 5- generates final file into host data folder (e.g. ${HOME}/HOWARD/data/example.howard.tsv)

```
$ docker exec HOWARD --input=/tool/docs/example.vcf --output=/data/example.howard.tsv --annotation=hgvs,symbol,outcome,location,CLINVAR,CLINVAR_CLNDN,COSMIC --calculation=VAF,BARCODE,NOMEN --prioritization=SOMATIC --translation=TSV --fields=NOMEN,PZFlag,PZScore,ALL --sort=PZFlag::DESC,PZScore:n:DESC
```

```
$ docker exec -ti HOWARD-CLI bash
[data]# HOWARD --help
```

# Docker 

HOWARD image presents a container that runs on CentOS, and includes yum modules and other tools dependencies:
- Java [1.8]
- bcftools/htslib [1.12]
- ANNOVAR [2019Oct24]
- snpEff [5.0e]

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
$ docker run --rm -v ${HOME}/HOWARD/data:/data -v ${HOME}/HOWARD/databases:/databases howard:latest --input=/tool/docs/example.vcf --output=/data/example.howard.tsv --annotation=hgvs,symbol,outcome,location,CLINVAR,CLINVAR_CLNDN,COSMIC --calculation=VAF,BARCODE,NOMEN --prioritization=SOMATIC --translation=TSV --fields=NOMEN,PZFlag,PZScore,ALL --sort=PZFlag::DESC,PZScore:n:DESC
```

Database download
-------------------

Databases are downloaded automatically by using annotation configuratin file, or options in command line (--annovar_databases, --snpeff_databases, assembly...).

Use a vcf file, such as HOWARD VCF example, to download ANNOVAR and snpEff databases (WITHOUT multithreading, "ALL" for all databases, "core" for core databases, "snpeff" for snpEff database, or a list of databases, or ANNOVAR code). Use this command multiple times for all needed databases and assembly (such as hg19, hg38, mm9).

```
$ docker run howard:latest --input=/tool/docs/example.vcf --output=/tool/docs/example.annotated.vcf --annotation=ALL,snpeff --thread=1 --verbose
```


> Note: For home made databases, refer to ```config.annotation.ini``` file to construct and configure your own database.

> Note: Beware of proxy configuration!


