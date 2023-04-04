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

- [Installation](#installation)
  - [Python](#python)
  - [Docker](#docker)
- [Quick HOWARD commands](#quick-howard-commands)
  - [Show VCF stats and overview](#show-vcf-stats-and-overview)
  - [Translate VCF into other format](#translate-vcf-into-other-format)
  - [Query](#query)
  - [Annotation](#annotation)
  - [Calculation](#calculation)
  - [Prioritization](#prioritization)
  - [Docker HOWARD-CLI](#docker-howard-cli)
- [Documentation](#documentation)



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

## Quick HOWARD commands

### Show VCF stats and overview

Show example VCF statistics and brief overview
```
howard analysis --input=tests/data/example.vcf.gz --stats --overview
```

### Translate VCF into other format

Translate VCF into CSV and show output file
```
howard analysis --input=tests/data/example.vcf.gz --output=tests/data/example.csv
cat tests/data/example.csv
```

Translate VCF into parquet, and show statistics on output parquet file (same as VCF)
```
howard analysis --input=tests/data/example.vcf.gz --output=tests/data/example.parquet
howard analysis --input=tests/data/example.parquet --stats
```

### Query

Select variants in VCF with REF and POS fields
```
howard analysis --input=tests/data/example.vcf.gz --query="SELECT * FROM variants WHERE REF = 'A' AND POS < 100000"
```

Select variants in VCF with INFO Tags criterions
```
howard analysis --input=tests/data/example.vcf.gz --param='{"explode_infos": true}' --query='SELECT "#CHROM", POS, REF, ALT, "INFO/DP", "INFO/CLNSIG", sample2, sample3 FROM variants WHERE "INFO/DP" >=50 OR "INFO/CLNSIG" NOT NULL'
```

### Annotation

VCF annotation with Parquet databases, output as VCF format, and show INFO
```
howard analysis --input=tests/data/example.vcf.gz --output=/tmp/example.howard.vcf.gz --annotation=tests/data/annotations/avsnp150.parquet,tests/data/annotations/dbnsfp42a.parquet,tests/data/annotations/gnomad211_genome.parquet --query='SELECT INFO FROM variants'
```

VCF annotation with Clinvar Parquet databases, output as TSV format and INFO tags as column, and query Clinvar results
```
howard analysis --input=tests/data/example.vcf.gz --output=/tmp/example.howard.tsv --annotation=tests/data/annotations/clinvar_20210123.parquet --param='{"explode_infos": true, "export_extra_infos": true}' --query='SELECT "INFO/CLNDN", count(*) AS count FROM variants WHERE "INFO/CLNDN" NOT NULL GROUP BY "INFO/CLNDN"'
```

### Calculation

Count number of variants by type
```
howard analysis --input=tests/data/example.full.vcf --calculations=vartype  --query='SELECT "INFO/VARTYPE", count(*) AS count FROM variants GROUP BY "INFO/VARTYPE" ORDER BY count DESC'
```

Extract hgvs from snpEff annotation and calculate NOMEN with default transcripts list
```
howard analysis --input=tests/data/example.ann.vcf.gz --param=tests/data/param.snpeff_hgvs.json --output=/tmp/example.snpeff_hgvs.vcf.gz --query='SELECT "#CHROM", POS, REF, ALT, "INFO/ANN" AS snpEff, "INFO/NOMEN" AS NOMEN FROM variants'
```
with 'param.snpeff_hgvs.json':
```
{
  "calculation": {
    "snpeff_hgvs": null,
    "NOMEN": {
        "options": {
            "hgvs_field": "snpeff_hgvs",
            "transcripts": "tests/data/tanscripts.tsv"
        }
    }
  },
  "explode_infos": "INFO/"
}
```
and 'transcripts.tsv':
```
NR_024540	WASH7P
NR_036266	MIR1302-9
NM_001346897	EGFR
NM_005228	EGFR
```

### Prioritization

Prioritize variants from criteria on INFO annotations (see 'prioritization_profiles.json')
```
howard analysis --input=tests/data/example.vcf.gz --prioritizations=tests/data/prioritization_profiles.json --output=/tmp/test.vcf --query='SELECT "#CHROM", POS, REF, ALT, "INFO/PZFlag", "INFO/PZScore" FROM variants ORDER BY "INFO/PZFlag" DESC, "INFO/PZScore" DESC' --param='{"explode_infos": "INFO/"}'
```

### Docker HOWARD-CLI

VCF annotation (Parquet, BCFTOOLS, ANNOVAR and snpEff) using HOWARD-CLI (snpEff and ANNOVAR databases will be automatically downloaded), and query list of genes with HGVS

```
docker exec HOWARD-CLI howard analysis --input=/tool/tests/data/example.vcf.gz --output=/data/example.howard.vcf.gz --annotation=snpeff,annovar:refGene,/tool/tests/data/annotations/refGene.bed.gz,/tool/tests/data/annotations/avsnp150.vcf.gz,tests/data/annotations/dbnsfp42a.parquet --param='{"explode_infos": true}' --query='SELECT "INFO/symbol", "INFO/AAChange_refGene" FROM variants WHERE "INFO/symbol" NOT NULL ORDER BY "INFO/symbol"'
```


Let's play within Docker HOWARD-CLI service!
```
docker exec -ti HOWARD-CLI bash
howard --help
```

## Documentation

More documentation in [docs/howard.md](docs/howard.md)


