# HOWARD User Guide

## Sommaire
- [Installation](#installation)
  - [Python](#python)
  - [Docker](#docker)
- [Databases](#databases)
    - [Download snpEff and Annovar](#download-snpeff-and-annovar)
    - [VCF and Parquet from Annovar](#vcf-and-parquet-from-annovar)
    - [VCF, CSV, Parquet and DuckDB](#vcf-csv-parquet-and-duckdb)
- [Process](#process)
- [Annotation](#annotation)
- [Calculation](#calculation)
- [Prioritization](#prioritization)
- [Query](#query)
- [Stats](#stats)
- [Convert](#convert)



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


# Databases

## Download snpEff and Annovar

```
usage: howard databases [-h] [--assembly ASSEMBLY] [--download-snpeff FOLDER] [--download-annovar FOLDER] [--download-annovar-files CODE]
                        [--download-annovar-url URL] [--config JSON] [--threads INTEGER] [--memory FLOAT[kMG]] [--verbosity LEVEL]

Download databases and needed files for howard and associated tools

Options:
  -h, --help            show this help message and exit
  --assembly ASSEMBLY   Genome Assembly
                        Default: 'hg19'

snpEff options:
  --download-snpeff FOLDER
                        Download snpEff databases within snpEff folder

Annovar options:
  --download-annovar FOLDER
                        Download Annovar databases within Annovar folder
  --download-annovar-files CODE
                        Download Annovar databases for a list of Annovar file code (see Annovar Doc)
                        Default: All available files
                        Example: refGene,gnomad211_exome,cosmic70,clinvar_202*,nci60
                        Note: refGene will be at leaset downloaded
                        Note2: Only file that not exists or with a different size will be downloaded
  --download-annovar-url URL
                        Download Annovar databases URL (see Annovar Doc)
                        Default: 'http://www.openbioinformatics.org/annovar/download/'

Shared options:
  --config JSON         Configuration file
                        Format: JSON
                        Default: {}
  --threads INTEGER     Number of threads to use
                        Format: INTEGER
                        Default: 1
  --memory FLOAT[kMG]   Memory to use
                        Format: FLOAT[kMG]
                        Default: 8G
  --verbosity LEVEL     Verbosity level (CRITICAL, ERROR, WARNING, INFO, DEBUG or NOTSET)
                        Default: INFO

Usage examples:
   howard databases --assembly=hg19 --download-annovar=/databases/annovar/current --download-annovar-files='refGene,gnomad_exome,dbnsfp42a,cosmic70,clinvar_202*,nci60' --download-snpeff=/databases/annovar/current
```

## VCF and Parquet from Annovar


```
usage: howard from_annovar [-h] --input FILE --output FILE --genome GENOME [--annovar-code CODE] [--to_parquet FILE] [--reduce_memory BOOL]
                           [--multi_variant BOOL] [--config JSON] [--threads INTEGER] [--memory FLOAT[kMG]] [--verbosity LEVEL]

(beta) Formatting Annovar database file to other format (VCF and Parquet). Exported Parquet file includes INFO/tags columns as VCF INFO columns had been exploded

Options:
  -h, --help            show this help message and exit
  --input FILE          Input file path
                        Format: BCF, VCF, TSV, CSV, PSV, Parquet or duckDB
                        Files can be compressesd (e.g. vcf.gz, tsv.gz)
  --output FILE         Output file path
                        Format: BCF, VCF, TSV, CSV, PSV, Parquet or duckDB
                        Files can be compressesd (e.g. vcf.gz, tsv.gz)
  --genome GENOME       Genome file in fasta format
                        Default: 'hg19.fa'

Annovar options:
  --annovar-code CODE   Annovar code, or database name. Usefull to name databases columns

Parquet options:
  --to_parquet FILE     Parquet file conversion

Modes options:
  --reduce_memory BOOL  Reduce memory option
                        Values: 'auto' (auto-detection), 'enable', 'disable'
                        default: 'auto'
  --multi_variant BOOL  Variant with multiple annotation lines
                        Values: 'auto' (auto-detection), 'enable', 'disable'
                        default: 'auto'

Shared options:
  --config JSON         Configuration file
                        Format: JSON
                        Default: {}
  --threads INTEGER     Number of threads to use
                        Format: INTEGER
                        Default: 1
  --memory FLOAT[kMG]   Memory to use
                        Format: FLOAT[kMG]
                        Default: 8G
  --verbosity LEVEL     Verbosity level (CRITICAL, ERROR, WARNING, INFO, DEBUG or NOTSET)
                        Default: INFO

Usage examples:
   howard from_annovar --input=/databases/annovar/current/hg19_nci60.txt --output=/databases/annotations/current/hg19/nci60.vcf.gz --to_parquet=/databases/annotations/current/hg19/nci60.parquet --annovar-code=nci60 --genome=/databases/genomes/current/hg19.fa --config=/tool/config/config.json --threads=8 
```

## VCF, CSV, Parquet and DuckDB

Use [Convert](#convert) tool

# Process

```
usage: howard process [-h] --input FILE --output FILE [--param JSON] [--annotations ANNOTATIONS] [--calculations OPERATIONS] [--prioritizations JSON]
                      [--query QUERY] [--config JSON] [--threads INTEGER] [--memory FLOAT[kMG]] [--verbosity LEVEL]

howard process tool manage genetic variations to:
- annotates genetic variants with multiple annotation databases/files and tools
- calculates and normalizes annotations
- prioritizes variants with profiles (list of citeria) to calculate scores and flags
- translates into various formats
- query genetic variants and annotations
- generates variants statistics

Options:
  -h, --help            show this help message and exit
  --input FILE          Input file path
                        Format: BCF, VCF, TSV, CSV, PSV, Parquet or duckDB
                        Files can be compressesd (e.g. vcf.gz, tsv.gz)
  --output FILE         Output file path
                        Format: BCF, VCF, TSV, CSV, PSV, Parquet or duckDB
                        Files can be compressesd (e.g. vcf.gz, tsv.gz)
  --param JSON          Parameters file
                        Format: JSON
                        Default: {}

Quick Processes options:
  --annotations ANNOTATIONS
                        Annotation with databases files, or with tools
                        Format: list of files in Parquet, VCF, BED, or keywords
                        For a Parquet/VCF/BED file, use file path (e.g. '/path/to/file.parquet')
                        For snpeff annotation, use keyword 'snpeff'
                        For Annovar annotation, use keyword 'annovar' with annovar code (e.g. 'annovar:refGene', 'annovar:cosmic70')
  --calculations OPERATIONS
                        Calculations on genetic variants information and genotype information
                        List of available calculations (see doc for more information):
                         VARTYPE  snpeff_hgvs  FINDBYPIPELINE  GENOTYPECONCORDANCE  BARCODE  TRIO  VAF  VAF_STATS  DP_STATS 
  --prioritizations JSON
                        Prioritization file in JSON format (defines profiles, see doc)
  --query QUERY         Query in SQL format
                        Format: SQL
                        Example: 'SELECT * FROM variants LIMIT 5'

Shared options:
  --config JSON         Configuration file
                        Format: JSON
                        Default: {}
  --threads INTEGER     Number of threads to use
                        Format: INTEGER
                        Default: 1
  --memory FLOAT[kMG]   Memory to use
                        Format: FLOAT[kMG]
                        Default: 8G
  --verbosity LEVEL     Verbosity level (CRITICAL, ERROR, WARNING, INFO, DEBUG or NOTSET)
                        Default: INFO

Usage examples:
   howard process --input=tests/data/example.vcf.gz --output=/tmp/example.annotated.vcf.gz --param=config/param.json 
   howard process --input=tests/data/example.vcf.gz --annotations='snpeff' --calculations='snpeff_hgvs,NOMEN' --prioritizations=config/prioritization_profiles.json --query='SELECT "INFO/NOMEN" FROM variants'
```


# Annotation

```
usage: howard annotation [-h] --input FILE --output FILE --annotations ANNOTATIONS [--config JSON] [--threads INTEGER] [--memory FLOAT[kMG]]
                         [--verbosity LEVEL]

Annotation is mainly based on a build-in Parquet annotation method, and tools such as BCFTOOLS, Annovar and snpEff. It uses available databases (see Annovar and snpEff) and homemade databases. Format of databases are: parquet, duckdb, vcf, bed, Annovar and snpEff (Annovar and snpEff databases are automatically downloaded, see howard databases tool). 

Options:
  -h, --help            show this help message and exit
  --input FILE          Input file path
                        Format: BCF, VCF, TSV, CSV, PSV, Parquet or duckDB
                        Files can be compressesd (e.g. vcf.gz, tsv.gz)
  --output FILE         Output file path
                        Format: BCF, VCF, TSV, CSV, PSV, Parquet or duckDB
                        Files can be compressesd (e.g. vcf.gz, tsv.gz)
  --annotations ANNOTATIONS
                        Annotation with databases files, or with tools
                        Format: list of files in Parquet, VCF, BED, or keywords
                        For a Parquet/VCF/BED file, use file path (e.g. '/path/to/file.parquet')
                        For snpeff annotation, use keyword 'snpeff'
                        For Annovar annotation, use keyword 'annovar' with annovar code (e.g. 'annovar:refGene', 'annovar:cosmic70')

Shared options:
  --config JSON         Configuration file
                        Format: JSON
                        Default: {}
  --threads INTEGER     Number of threads to use
                        Format: INTEGER
                        Default: 1
  --memory FLOAT[kMG]   Memory to use
                        Format: FLOAT[kMG]
                        Default: 8G
  --verbosity LEVEL     Verbosity level (CRITICAL, ERROR, WARNING, INFO, DEBUG or NOTSET)
                        Default: INFO

Usage examples:
   howard annotation --input=tests/data/example.vcf.gz --output=/tmp/example.howard.vcf.gz --annotations='tests/databases/annotations/hg19/avsnp150.parquet,tests/databases/annotations/hg19/dbnsfp42a.parquet,tests/databases/annotations/hg19/gnomad211_genome.parquet' 
   howard annotation --input=tests/data/example.vcf.gz --output=/tmp/example.howard.tsv --annotations='annovar:refGene,snpeff,tests/databases/annotations/hg19/clinvar_20210123.parquet'
```

# Calculation

```
usage: howard calculation [-h] --input FILE --output FILE --calculations OPERATIONS [--show_calculations] [--hgvs_field TAG] [--transcripts FILE]
                          [--trio_pedigree JSON] [--config JSON] [--threads INTEGER] [--memory FLOAT[kMG]] [--verbosity LEVEL]

Calculation processes variants information to generate new information, such as: identify variation type (VarType), harmonizes allele frequency (VAF) and calculate sttistics (VAF_stats), extracts Nomen (transcript, cNomen, pNomen...) from an HGVS field (e.g. snpEff, Annovar) with an optional list of personalized transcripts, generates VaRank format barcode, identify trio inheritance.

Options:
  -h, --help            show this help message and exit
  --input FILE          Input file path
                        Format: BCF, VCF, TSV, CSV, PSV, Parquet or duckDB
                        Files can be compressesd (e.g. vcf.gz, tsv.gz)
  --output FILE         Output file path
                        Format: BCF, VCF, TSV, CSV, PSV, Parquet or duckDB
                        Files can be compressesd (e.g. vcf.gz, tsv.gz)
  --calculations OPERATIONS
                        Calculations on genetic variants information and genotype information
                        List of available calculations (see doc for more information):
                         VARTYPE  snpeff_hgvs  FINDBYPIPELINE  GENOTYPECONCORDANCE  BARCODE  TRIO  VAF  VAF_STATS  DP_STATS 
  --show_calculations   Show available calculation operations

NOMEN calculation options:
  --hgvs_field TAG      HGVS INFO/tag containing a list o HGVS annotations
                        default: 'hgvs'
  --transcripts FILE    Transcripts file in TSV format
                        Format: Transcript in first column, optional Gene in second column 
                        default: None

TRIO calculation options:
  --trio_pedigree JSON  Pedigree Trio for trio inheritance calculation
                        Format: JSON file or dict (e.g. 'trio.ped.json', '{"father":"sample1", "mother":"sample2", "child":"sample3"}') 
                        default: None

Shared options:
  --config JSON         Configuration file
                        Format: JSON
                        Default: {}
  --threads INTEGER     Number of threads to use
                        Format: INTEGER
                        Default: 1
  --memory FLOAT[kMG]   Memory to use
                        Format: FLOAT[kMG]
                        Default: 8G
  --verbosity LEVEL     Verbosity level (CRITICAL, ERROR, WARNING, INFO, DEBUG or NOTSET)
                        Default: INFO

Usage examples:
   howard calculation --input=tests/data/example.full.vcf --output=/tmp/example.calculation.tsv --calculations='vartype' 
   howard calculation --input=tests/data/example.ann.vcf.gz --output=/tmp/example.calculated.tsv --calculations='snpeff_hgvs,NOMEN' --hgvs_field=snpeff_hgvs --transcripts=tests/data/transcripts.tsv 
```

# Prioritization

```
usage: howard prioritization [-h] --input FILE --output FILE --prioritizations JSON [--profiles PROFILES] [--default_profile PROFILE] [--pzfields PZFIELD]
                             [--prioritization_score_mode MODE] [--config JSON] [--threads INTEGER] [--memory FLOAT[kMG]] [--verbosity LEVEL]

Prioritization algorithm uses profiles to flag variants (as passed or filtered), calculate a prioritization score, and automatically generate a comment for each variants (example: 'polymorphism identified in dbSNP. associated to Lung Cancer. Found in ClinVar database'). Prioritization profiles are defined in a configuration file in JSON format. A profile is defined as a list of annotation/value, using wildcards and comparison options (contains, lower than, greater than, equal...). Annotations fields may be quality values (usually from callers, such as 'DP') or other annotations fields provided by annotations tools, such as HOWARD itself (example: COSMIC, Clinvar, 1000genomes, PolyPhen, SIFT). Multiple profiles can be used simultaneously, which is useful to define multiple validation/prioritization levels (example: 'standard', 'stringent', 'rare variants', 'low allele frequency').

Options:
  -h, --help            show this help message and exit
  --input FILE          Input file path
                        Format: BCF, VCF, TSV, CSV, PSV, Parquet or duckDB
                        Files can be compressesd (e.g. vcf.gz, tsv.gz)
  --output FILE         Output file path
                        Format: BCF, VCF, TSV, CSV, PSV, Parquet or duckDB
                        Files can be compressesd (e.g. vcf.gz, tsv.gz)
  --prioritizations JSON
                        Prioritization file in JSON format (defines profiles, see doc)

Prioritization options:
  --profiles PROFILES   Prioritization profiles to use (based on file in JSON)
                        default: all profiles available
  --default_profile PROFILE
                        Prioritization profile by default (see doc)
                        default: First profile in JSON file
  --pzfields PZFIELD    Prioritization fields to provide (see doc)
                        available: PZScore, PZFlag, PZTags, PZComment, PZInfos
                        default: PZScore,PZFlag
  --prioritization_score_mode MODE
                        Prioritization Score mode (see doc)
                        available: HOWARD (increment score), VaRank (max score)
                        default: HOWARD

Shared options:
  --config JSON         Configuration file
                        Format: JSON
                        Default: {}
  --threads INTEGER     Number of threads to use
                        Format: INTEGER
                        Default: 1
  --memory FLOAT[kMG]   Memory to use
                        Format: FLOAT[kMG]
                        Default: 8G
  --verbosity LEVEL     Verbosity level (CRITICAL, ERROR, WARNING, INFO, DEBUG or NOTSET)
                        Default: INFO

Usage examples:
   howard prioritization --input=tests/data/example.vcf.gz --output=/tmp/example.prioritized.vcf.gz --prioritizations=config/prioritization_profiles.json --profiles='default,GERMLINE'
```


# Query

```
usage: howard query [-h] [--input FILE] [--output FILE] --query QUERY [--explode_infos] [--config JSON] [--threads INTEGER] [--memory FLOAT[kMG]]
                    [--verbosity LEVEL]

Query genetic variations in SQL format. Data can be loaded into 'variants' table from various formats (e.g. VCF, TSV, Parquet...). Using --explode_infos allow query on INFO/tag annotations. SQL query can also use external data within the request, such as a Parquet file(s).  

Options:
  -h, --help           show this help message and exit
  --input FILE         Input file path
                       Format: BCF, VCF, TSV, CSV, PSV, Parquet or duckDB
                       Files can be compressesd (e.g. vcf.gz, tsv.gz)
  --output FILE        Output file path
                       Format: BCF, VCF, TSV, CSV, PSV, Parquet or duckDB
                       Files can be compressesd (e.g. vcf.gz, tsv.gz)
  --query QUERY        Query in SQL format
                       Format: SQL
                       Example: 'SELECT * FROM variants LIMIT 5'
  --explode_infos      Explode VCF INFO/Tag into 'variants' table columns
                       default: False

Shared options:
  --config JSON        Configuration file
                       Format: JSON
                       Default: {}
  --threads INTEGER    Number of threads to use
                       Format: INTEGER
                       Default: 1
  --memory FLOAT[kMG]  Memory to use
                       Format: FLOAT[kMG]
                       Default: 8G
  --verbosity LEVEL    Verbosity level (CRITICAL, ERROR, WARNING, INFO, DEBUG or NOTSET)
                       Default: INFO

Usage examples:
   howard query --input=tests/data/example.vcf.gz --query="SELECT * FROM variants WHERE REF = 'A' AND POS < 100000" 
   howard query --input=tests/data/example.vcf.gz --explode_infos --query='SELECT "#CHROM", POS, REF, ALT, "INFO/DP", "INFO/CLNSIG", sample2, sample3 FROM variants WHERE "INFO/DP" >= 50 OR "INFO/CLNSIG" NOT NULL ORDER BY "INFO/DP" DESC' 
   howard query --query="SELECT * FROM 'tests/databases/annotations/hg19/dbnsfp42a.parquet' WHERE \"INFO/Interpro_domain\" NOT NULL ORDER BY \"INFO/SiPhy_29way_logOdds_rankscore\" DESC" 
   howard query --query="SELECT \"#CHROM\" AS \"#CHROM\", POS AS POS, '' AS ID, REF AS REF, ALT AS ALT, '' AS QUAL, '' AS FILTER, STRING_AGG(INFO, ';') AS INFO FROM 'tests/databases/annotations/hg19/*.parquet' GROUP BY \"#CHROM\", POS, REF, ALT" --output=/tmp/full_annotation.tsv 
```

# Stats

```
usage: howard stats [-h] --input FILE [--config JSON] [--threads INTEGER] [--memory FLOAT[kMG]] [--verbosity LEVEL]

Statistics on genetic variations, such as: number of variants, number of samples, statistics by chromosome, genotypes by samples...

Options:
  -h, --help           show this help message and exit
  --input FILE         Input file path
                       Format: BCF, VCF, TSV, CSV, PSV, Parquet or duckDB
                       Files can be compressesd (e.g. vcf.gz, tsv.gz)

Shared options:
  --config JSON        Configuration file
                       Format: JSON
                       Default: {}
  --threads INTEGER    Number of threads to use
                       Format: INTEGER
                       Default: 1
  --memory FLOAT[kMG]  Memory to use
                       Format: FLOAT[kMG]
                       Default: 8G
  --verbosity LEVEL    Verbosity level (CRITICAL, ERROR, WARNING, INFO, DEBUG or NOTSET)
                       Default: INFO

Usage examples:
   howard stats --input=tests/data/example.vcf.gz 
```

# Convert

```
usage: howard convert [-h] --input FILE --output FILE [--export_infos] [--export_infos_prefix PREFIX] [--config JSON] [--threads INTEGER]
                      [--memory FLOAT[kMG]] [--verbosity LEVEL]

Convert genetic variations file to another format. Multiple format are available, such as usual and official VCF and BCF format, but also other formats such as TSV, CSV, PSV and Parquet/duckDB. These formats need a header '.hdr' file to take advantage of the power of howard (especially through INFO/tag definition), and using howard convert tool automatically generate header file fo futher use. 

Options:
  -h, --help            show this help message and exit
  --input FILE          Input file path
                        Format: BCF, VCF, TSV, CSV, PSV, Parquet or duckDB
                        Files can be compressesd (e.g. vcf.gz, tsv.gz)
  --output FILE         Output file path
                        Format: BCF, VCF, TSV, CSV, PSV, Parquet or duckDB
                        Files can be compressesd (e.g. vcf.gz, tsv.gz)
  --export_infos        Export VCF INFO/Tag into columns file
                        Available only for non VCF format file (e.g. TSV, Parquet...) 
                        default: False
  --export_infos_prefix PREFIX
                        VCF INFO/Tag prefix for exported columns
                        default: 'INFO/'

Shared options:
  --config JSON         Configuration file
                        Format: JSON
                        Default: {}
  --threads INTEGER     Number of threads to use
                        Format: INTEGER
                        Default: 1
  --memory FLOAT[kMG]   Memory to use
                        Format: FLOAT[kMG]
                        Default: 8G
  --verbosity LEVEL     Verbosity level (CRITICAL, ERROR, WARNING, INFO, DEBUG or NOTSET)
                        Default: INFO

Usage examples:
   howard convert --input=tests/data/example.vcf.gz --output=/tmp/example.tsv --export_infos   howard convert --input=tests/data/example.vcf.gz --output=/tmp/example.parquet
```

