HOWARD
===

HOWARD annotates and prioritizes genetic variations, calculates and normalizes annotations, translates files in multiple formats (e.g. vcf, tsv, parquet) and generates variants statistics.

HOWARD annotation is mainly based on a build-in Parquet annotation method, and external tools such as BCFTOOLS, ANNOVAR, snpEff and Exomiser (see docs, automatically downloaded if needed). Parquet annotation uses annotation database in VCF or BED format, in mutliple file format: Parquet/duckdb, VCF, BED, TSV, CSV, TBL, JSON. 

HOWARD calculation processes variants information to calculate new information, such as: harmonizes allele frequency (VAF), extracts Nomen (transcript, cNomen, pNomen...) from HGVS fields with an optional list of personalized transcripts, generates VaRank format barcode.

HOWARD prioritization algorithm uses profiles to flag variants (as passed or filtered), calculate a prioritization score, and automatically generate a comment for each variants (example: 'polymorphism identified in dbSNP. associated to Lung Cancer. Found in ClinVar database').Prioritization profiles are defined in a configuration file. A profile is defined as a list of annotation/value, using wildcards and comparison options (contains, lower than, greater than, equal...). Annotations fields may be quality values (usually from callers, such as 'GQ', 'DP') or other annotations fields provided by annotations tools, such as HOWARD itself (example: COSMIC, Clinvar, 1000genomes, PolyPhen, SIFT). Multiple profiles can be used simultaneously, which is useful to define multiple validation/prioritization levels (example: 'standard', 'stringent', 'rare variants', 'low allele frequency').

HOWARD translates VCF format into multiple formats (e.g. VCF, TSV, Parquet), by sorting variants using specific fields (example : 'prioritization score', 'allele frequency', 'gene symbol'), including/excluding annotations/fields, including/excluding variants, adding fixed columns.

HOWARD generates statistics files with a specific algorithm, snpEff and BCFTOOLS.

HOWARD is multithreaded through the number of variants and by database (data-scaling).


# Getting Started

- [Installation](#installation)
  - [Python](#python)
  - [Docker](#docker)
- [Quick HOWARD commands](#quick-howard-commands)
  - [Stats](#stats)
  - [Convert](#convert)
  - [Query](#query)
  - [Annotation](#annotation)
  - [Calculation](#calculation)
  - [Prioritization](#prioritization)
  - [Docker HOWARD-CLI](#docker-howard-cli)
- [Documentation](#documentation)



# Installation

## Python

Install using Python Pip tool:
```
python -m pip install -e .
```

## Docker

In order to build, setup and create a persitent CLI (running container), docker-compose command build images and launch services as containers.

```
docker-compose up -d
```

A Command Line Interface container (HOWARD-CLI) is started with host data and databases folders mounted (by default in ${HOME}/HOWARD folder)

# Databases

Multiple databases can be automatically downloaded with databases tool, such as:

| database | description | Link |
|-----|-----|---|
| [Genome](https://genome.ucsc.edu/cgi-bin/hgGateway) | Genome Reference Consortium Human  |
| [Annovar](https://annovar.openbioinformatics.org/en/latest/) | ANNOVAR is an efficient software tool to utilize update-to-date information to functionally annotate genetic variants detected from diverse genomes |  |
| [snpEff](https://pcingola.github.io/SnpEff/) | Genetic variant annotation, and functional effect prediction toolbox |
| [refSeq](https://www.ncbi.nlm.nih.gov/refseq/) | A comprehensive, integrated, non-redundant, well-annotated set of reference sequences including genomic, transcript, and protein |
| [dbSNP](https://www.ncbi.nlm.nih.gov/snp/) | dbSNP contains human single nucleotide variations, microsatellites, and small-scale insertions and deletions along with publication, population frequency, molecular consequence, and genomic and RefSeq mapping information for both common variations and clinical mutations |
| [dbNSFP](https://sites.google.com/site/jpopgen/dbNSFP) | dbNSFP is a database developed for functional prediction and annotation of all potential non-synonymous single-nucleotide variants (nsSNVs) in the human genome |
| [AlphaMissense](https://github.com/google-deepmind/alphamissense) | AlphaMissense model implementation |
| [Exomiser](https://www.sanger.ac.uk/tool/exomiser/) | The Exomiser is a Java program that finds potential disease-causing variants from whole-exome or whole-genome sequencing data |


- Download Multiple databases in the same time for assembly 'hg19' (can take a while):
```
howard databases --assembly=hg19 --download-genomes=/databases/genomes/current --download-genomes-provider=UCSC --download-genomes-contig-regex='chr[0-9XYM]+$' --download-annovar=/databases/annovar/current --download-annovar-files='refGene,cosmic70,nci60' --download-snpeff=/databases/snpeff/current --download-refseq=/databases/refseq/current --download-refseq-format-file='ncbiRefSeq.txt' --download-dbnsfp=/databases/dbnsfp/current --download-dbnsfp-release='4.4a' --download-dbnsfp-subdatabases --download-alphamissense=/databases/alphamissense/current --download-exomiser=/databases/exomiser/current --download-dbsnp=/databases/dbsnp/current --download-dbsnp-vcf --threads=8
```

Databases can be home-made generated, starting with a existing annotation file, especially using [HOWARD convert](#convert) tool. These files need to contain specific fields (depending on the annotation type):
- variant annotation: '#CHROM', 'POS', 'ALT', 'REF'
- region annotation: '#CHROM', 'START', 'STOP'

Each database annotation file is associated with a 'header' file ('.hdr'), in VCF header format, to describe annotations within the database.


# Quick HOWARD commands

## Stats

Statistics on genetic variations, such as: number of variants, number of samples, statistics by chromosome, genotypes by samples, annotations...
Theses statsitics can be applied to VCF files and all database annotation files.

- Show example VCF statistics and brief overview
```
howard stats --input=tests/data/example.vcf.gz
```

## Convert

Convert genetic variations file to another format. Multiple format are available, such as usual and official VCF format, but also other formats such as TSV, CSV, TBL, JSON and Parquet/duckDB. These formats need a header '.hdr' file to take advantage of the power of howard (especially through INFO/tag definition), and using howard convert tool automatically generate header file fo futher use (otherwise, an default '.hdr' file is generated).

- Translate VCF into TSV, export INFO/tags into columns, and show output file
```
howard convert --input=tests/data/example.vcf.gz --output=/tmp/example.tsv --explode_infos && cat /tmp/example.tsv
```

- Translate VCF into parquet
```
howard convert --input=tests/data/example.vcf.gz --output=/tmp/example.parquet
```

- and show statistics on output parquet file (same as VCF)
```
howard stats --input=/tmp/example.parquet
```

## Query

Query genetic variations in SQL format. Data can be loaded into 'variants' table from various formats (e.g. VCF, TSV, Parquet...). Using --explode_infos allow query on INFO/tag annotations. SQL query can also use external data within the request, such as a Parquet file(s).

- Select variants in VCF with REF and POS fields
```
howard query --input=tests/data/example.vcf.gz --query="SELECT * FROM variants WHERE REF = 'A' AND POS < 100000"
```

- Select variants in VCF with INFO Tags criterions
```
howard query --input=tests/data/example.vcf.gz --explode_infos --query='SELECT "#CHROM", POS, REF, ALT, "DP", "CLNSIG", sample2, sample3 FROM variants WHERE "DP" >= 50 OR "CLNSIG" NOT NULL ORDER BY "CLNSIG" DESC, "DP" DESC'
```

- Query a Parquet file with specific columns (e.g. from VCF convertion to Parquet)
```
howard query --query="SELECT * FROM 'tests/databases/annotations/hg19/dbnsfp42a.parquet' WHERE \"INFO/Interpro_domain\" NOT NULL ORDER BY \"INFO/SiPhy_29way_logOdds_rankscore\" DESC"
```

- Query multiple Parquet files, merge INFO columns, and extract as TSV (in VCF format)
```
howard query --query="SELECT \"#CHROM\" AS \"#CHROM\", POS AS POS, '' AS ID, REF AS REF, ALT AS ALT, '' AS QUAL, '' AS FILTER, STRING_AGG(INFO, ';') AS INFO FROM 'tests/databases/annotations/hg19/*.parquet' GROUP BY \"#CHROM\", POS, REF, ALT" --output=/tmp/full_annotation.tsv
```


## Annotation

Annotation is mainly based on a build-in Parquet annotation method, and tools such as BCFTOOLS, Annovar, snpEff and Exomiser. It uses available databases and homemade databases. Format of databases are: Parquet/duckdb, VCF, BED, Annovar and snpEff (Annovar and snpEff databases are automatically downloaded, see howard databases tool).

- VCF annotation with Parquet and VCF databases, output as VCF format
```
howard annotation --input=tests/data/example.vcf.gz --output=/tmp/example.howard.vcf.gz --annotations='tests/databases/annotations/hg19/dbnsfp42a.parquet,tests/databases/annotations/hg19/gnomad211_genome.parquet,tests/databases/annotations/hg19/cosmic70.vcf.gz'
```

- VCF annotation with Clinvar Parquet, Annovar refGene and snpEff databases, output as TSV format
```
howard annotation --input=tests/data/example.vcf.gz --output=/tmp/example.howard.tsv --annotations='annovar:refGene,snpeff,tests/databases/annotations/hg19/clinvar_20210123.parquet'
```

- VCF annotation with all available database annotation files in Parquet format (within the database annotation folder in configuration):
```
howard annotation --input=tests/data/example.vcf.gz --output=/tmp/example.howard.tsv --assembly='hg19' --annotations='ALL:parquet'
```

## Calculation

Calculation processes variants information to generate new information, such as: identify variation type (VarType), harmonizes allele frequency (VAF) and calculate sttistics (VAF_stats), extracts Nomen (transcript, cNomen, pNomen...) from an HGVS field (e.g. snpEff, Annovar) with an optional list of personalized transcripts, generates VaRank format barcode, identify trio inheritance.

- Identify variant types
```
howard calculation --input=tests/data/example.full.vcf --output=/tmp/example.calculation.tsv --calculations='vartype'
```
- and generate a table of variant type count
```
howard query --input=/tmp/example.calculation.tsv --explode_infos --query='SELECT "VARTYPE" AS 'VariantType', count(*) AS 'Count' FROM variants GROUP BY "VARTYPE" ORDER BY count DESC'
```

- Calculate NOMEN by extracting hgvs from snpEff annotation and identifying default transcripts from a list
```
howard calculation --input=tests/data/example.ann.vcf.gz --output=/tmp/example.NOMEN.vcf.gz --calculations='snpeff_hgvs,NOMEN' --hgvs_field='snpeff_hgvs' --transcripts=tests/data/transcripts.tsv && gzip -dc /tmp/example.NOMEN.vcf.gz | grep "##" -v | head -n2
```
- and query NOMEN for gene 'EGFR'
```
howard query --input=/tmp/example.NOMEN.vcf.gz --explode_infos --query="SELECT \"NOMEN\" AS 'NOMEN' FROM variants WHERE \"GNOMEN\" == 'EGFR'"
```

## Prioritization

Prioritization algorithm uses profiles to flag variants (as passed or filtered), calculate a prioritization score, and automatically generate a comment for each variants (example: 'polymorphism identified in dbSNP. associated to Lung Cancer. Found in ClinVar database'). Prioritization profiles are defined in a configuration file in JSON format. A profile is defined as a list of annotation/value, using wildcards and comparison options (contains, lower than, greater than, equal...). Annotations fields may be quality values (usually from callers, such as 'DP') or other annotations fields provided by annotations tools, such as HOWARD itself (example: COSMIC, Clinvar, 1000genomes, PolyPhen, SIFT). Multiple profiles can be used simultaneously, which is useful to define multiple validation/prioritization levels (example: 'standard', 'stringent', 'rare variants', 'low allele frequency').

- Prioritize variants from criteria on INFO annotations for profiles 'default' and 'GERMLINE' (see 'prioritization_profiles.json'), and query variants on prioritization tags
```
howard prioritization --input=tests/data/example.vcf.gz --output=/tmp/example.prioritized.vcf.gz --prioritizations=config/prioritization_profiles.json --profiles='default,GERMLINE' --pzfields='PZFlag,PZScore,PZComment'
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
        "tests/databases/annotations/hg19/avsnp150.parquet": {
          "INFO": null
        },
        "tests/databases/annotations/hg19/dbnsfp42a.parquet": {
          "INFO": null
        },
        "tests/databases/annotations/hg19/gnomad211_genome.parquet": {
          "INFO": null
        }
      }
    },
    "bcftools": {
      "annotations": {
        "tests/databases/annotations/hg19/cosmic70.vcf.gz": {
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


## Docker HOWARD-CLI

VCF annotation (Parquet, BCFTOOLS, ANNOVAR, snpEff and Exomiser) using HOWARD-CLI (snpEff and ANNOVAR databases will be automatically downloaded), and query list of genes with HGVS

```
docker exec --workdir=/tool HOWARD-CLI howard process --config=config/config.json --param=config/param.json --input=tests/data/example.vcf.gz --output=/tmp/example.process.tsv --query='SELECT "INFO/NOMEN" AS "NOMEN", "INFO/PZFlag" AS "PZFlag", "INFO/PZScore" AS "PZScore", "INFO/PZComment" AS "PZComment" FROM variants ORDER BY "INFO/PZScore" DESC' --explode_infos
```

Let's play within Docker HOWARD-CLI service!
```
docker exec -ti HOWARD-CLI bash
howard --help
```

## Documentation

More documentation in [docs/howard.md](docs/howard.md)


