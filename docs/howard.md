# HOWARD User Guide

## Sommaire
- [Installation](#installation)
  - [Python](#python)
  - [Docker](#docker)
  - [Databases](#databases)
    - [snpEff](#snpeff)
    - [Annovar](#annovar)
    - [VCF and BED](#vcf-and-bed)
    - [Parquet and DuckDB](#parquet-and-duckdb)
- [Annotation](#annotation)
- [Calculation](#calculation)
- [Prioritization](#prioritization)



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

Docker service HOWARD-setup will create main HOWARD image, and download useful databases to start with HOWARD analysis.

A Command Line Interface container (HOWARD-CLI) is started with host data and databases folders mounted (by default in ${HOME}/HOWARD folder)

List of databases downloaded in HOWARD setup (for assembly hg19 and hg38):
- Annovar (see [Annovar documentation](https://annovar.openbioinformatics.org/en/latest/user-guide/download/#annovar-main-package) for more information, especially about releases)
  - refGene: FASTA sequences for all annotated transcripts in RefSeq Gene, HGVS annotation, genetic variation impact [[UCSC](https://genome.ucsc.edu/)]
  - gnomad_exome: gnomAD genome collection, aggregating and harmonizing exome sequencing data from a wide variety of large-scale sequencing projects [[gnomAD](https://gnomad.broadinstitute.org/)]
  - dbnsfp42a: functional prediction and annotation of all potential non-synonymous single-nucleotide variants [[dbNSFP](https://sites.google.com/site/jpopgen/dbNSFP)]
  - clinvar_202*: CLINVAR database with Variant Clinical Significance (unknown, untested, non-pathogenic, probable-non-pathogenic, probable-pathogenic, pathogenic, drug-response, histocompatibility, other) and Variant disease name [[ClinVar](https://www.ncbi.nlm.nih.gov/clinvar/)]
  - cosmic70: Catalogue Of Somatic Mutations In Cancer [[COSMIC](https://cancer.sanger.ac.uk/cosmic)]
- snpEff: Genetic variant annotation and functional effect prediction, effects of genetic variants on genes and proteins (such as amino acid changes) [[snpEff](http://pcingola.github.io/SnpEff/)]


### Databases

#### snpEff

Download snpEff databases:
```
howard download --download-snpeff=/databases/snpeff/current --assembly=hg19,hg38
```

#### Annovar

Download Annovar databases:
```
howard download --download-annovar=/databases/annovar/current --download-annovar-files=refGene,gnomad211_exome,cosmic70,clinvar_202*,nci60 --assembly=hg19,hg38
```

#### VCF and BED

#### Parquet and DuckDB

## Query

## Annotation

## Calculation

## Prioritization


