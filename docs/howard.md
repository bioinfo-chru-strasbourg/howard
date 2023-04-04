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

A Command Line Interface container (HOWARD-CLI) is started with host data and databases folders mounted (by default in ${HOME}/HOWARD folder)

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

## Annotation

## Calculation

## Prioritization


