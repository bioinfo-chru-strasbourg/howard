# HOWARD Parameters

HOWARD Parameters JSON file defined parameters to process annotations, prioritization, calculations, convertions and queries.

## Table of contents

- [HOWARD Parameters](#howard-parameters)
   - [annotation](#annotation)
      - [parquet](#annotationparquet)
         - [annotations](#annotationparquetannotations)
      - [bcftools](#annotationbcftools)
         - [annotations](#annotationbcftoolsannotations)
      - [annovar](#annotationannovar)
         - [annotations](#annotationannovarannotations)
         - [options](#annotationannovaroptions)
      - [snpeff](#annotationsnpeff)
         - [options](#annotationsnpeffoptions)
         - [stats](#annotationsnpeffstats)
         - [csvStats](#annotationsnpeffcsvstats)
      - [exomiser](#annotationexomiser)
         - [release](#annotationexomiserrelease)
      - [options](#annotationoptions)
         - [annotations_update](#annotationoptionsannotations_update)
         - [annotations_append](#annotationoptionsannotations_append)
   - [calculation](#calculation)
   - [prioritization](#prioritization)
      - [prioritizations](#prioritizationprioritizations)
      - [profiles](#prioritizationprofiles)
      - [default_profile](#prioritizationdefault_profile)
      - [pzfields](#prioritizationpzfields)
      - [prioritization_score_mode](#prioritizationprioritization_score_mode)
   - [hgvs](#hgvs)
      - [use_gene](#hgvsuse_gene)
      - [use_exon](#hgvsuse_exon)
      - [use_protein](#hgvsuse_protein)
      - [add_protein](#hgvsadd_protein)
      - [full_format](#hgvsfull_format)
      - [codon_type](#hgvscodon_type)
      - [refgene](#hgvsrefgene)
      - [refseqlink](#hgvsrefseqlink)
   - [stats](#stats)
      - [stats_md](#statsstats_md)
      - [stats_json](#statsstats_json)
   - [query](#query)
   - [explode_infos](#explode_infos)
   - [explode_infos_prefix](#explode_infos_prefix)
   - [explode_infos_fields](#explode_infos_fields)
   - [parquet_partitions](#parquet_partitions)
   - [order_by](#order_by)
   - [query_limit](#query_limit)
   - [query_print_mode](#query_print_mode)
   - [threads](#threads)


Examples: 
```json
# Parameters for annotation, calculation, prioritization, HGVS annotation...
{
   "annotation": {
      "parquet": {
         "annotations": {
            "/path/to/database1.parquet": {
               ...
            },
            ...
         }
      }
      "annovar": {
         "annotations": {
            "annovar_keyword1": {
               ...
            },
            ...
         }
      }
   },
   "calculation": {
      "operation1": null,
      "operation2": {
        "options": {
          "option1": "value1",
          "option2": "value2"
          ...
        }
        ...
      }
   },
   "prioritization": {
      "prioritizations": "config/prioritization_profiles.json",
      "profiles": ["GENOME", "GERMLINE"],
      "default_profile": "GERMLINE",
      "pzfields": ["PZScore", "PZFlag", "PZComment"],
      "prioritization_score_mode": "VaRank"
   },
   "hgvs": {
      "full_format": true,
      "use_exon": true
   }
}
```

## annotation

Annotation process using HOWARD algorithms or external tools.

For HOWARD Parquet algorithm, specify the list of database files available (formats such as Parquet, VCF, TSV, duckDB, JSON). This parameter enables users to select specific database fields and optionally rename them. Use 'INFO' keyword to select all fields within the database. If a full path is not provided, the system will automatically detect files within database folders (see Configuration doc) and assembly (see Parameter option).

For external tools, such as Annovar, snpEff and Exomiser, specify parameters such as annotation keywords (Annovar) and options (depending on the tool).

Examples: 
```json
# Annotation with multiple tools in multiple formats with multiple options
"annotation": {
   "parquet": {
      "annotations": {
         "/path/to/database1.parquet": {
            ...
         },
         ...,
      }
   }
   "annovar": {
      "annotations": {
         "annovar_keyword1": {
            ...
         },
         ...
      }
   }
}
```

### annotation::parquet

Annotation process using HOWARD parquet algorithm. Provide a list of database files and annotation fields.

Examples: 
```json
# Annotation with multiple databases in multiple formats
"parquet": {
   "annotations": {
      "/path/to/database1.parquet": {
         "field1": null
         "field2": "field2_renamed",
         ...
      },
      "database2.vcf.gz": {
         ...
      },
      ...
   }
}
```

#### annotation::parquet::annotations

Specify the list of database files available in multiple formats such as Parquet, VCF, BED, TSV, duckDB, JSON.

Type: ```str```

Default: ```None```

Examples: 
```json
# Annotation with dbSNP database  with all fields
"annotations": {
   "tests/databases/annotations/current/hg19/avsnp150.parquet": {
      "INFO": null
   }
}

# Annotation with dbNSFP (only PolyPhen, ClinVar and REVEL score), and rename fields
"annotations": {
   "tests/databases/annotations/current/hg19/dbnsfp42a.parquet": {
      "Polyphen2_HDIV_pred": "PolyPhen",
      "ClinPred_pred": "ClinVar",
      "REVEL_score": null
   }
}

# Annotation with refSeq as a BED file
"annotations": {
   "tests/databases/annotations/current/hg19/refGene.bed": {
      "INFO": null
   }
}

# Annotation with dbNSFP REVEL annotation (as a VCF file) within configured annotation databases folders (default: '~/howard/databases/annotations/current') and assembly (default: 'hg19')
"annotations": {
   "dbnsfp42a.REVEL.vcf.gz": {
      "REVEL_score": null,
      "REVEL_rankscore": null
   }
}

```

### annotation::bcftools

Annotation process using BCFTools. Provide a list of database files and annotation fields.

Examples: 
```json
# Annotation with multiple databases in multiple formats
"parquet": {
   "bcftools": {
      "/path/to/database1.vcf.gz": {
         "field1": null
         "field2": "field2_renamed",
         ...
      },
      "database2.bed.gz": {
         ...
      },
      ...
   }
}
```

#### annotation::bcftools::annotations

Specify the list of database files in formats VCF or BED. Files need to be compressed and indexed.

Type: ```str```

Default: ```None```

Examples: 
```json
# Annotation with dbSNP database  with all fields
"annotations": {
   "tests/databases/annotations/current/hg19/avsnp150.vcf.gz": {
      "INFO": null
   }
}

# Annotation with dbNSFP (only PolyPhen, ClinVar and REVEL score), and rename fields
"annotations": {
   "tests/databases/annotations/current/hg19/dbnsfp42a.vcf.gz": {
      "Polyphen2_HDIV_pred": "PolyPhen",
      "ClinPred_pred": "ClinVar",
      "REVEL_score": null
   }
}

# Annotation with refSeq as a BED file
"annotations": {
   "tests/databases/annotations/current/hg19/refGene.bed": {
      "INFO": null
   }
}

# Annotation with dbNSFP REVEL annotation (as a VCF file) within configured annotation databases folders (default: '~/howard/databases/annotations/current') and assembly (default: 'hg19')
"annotations": {
   "dbnsfp42a.REVEL.vcf.gz": {
      "REVEL_score": null,
      "REVEL_rankscore": null
   }
}

```

### annotation::annovar

Annotation process using Annovar tool. Provides a list of keywords to select Annovar databases, and defines Annovar options (see [Annovar documentation](https://annovar.openbioinformatics.org)).

Examples: 
```json
# Annotation with multiple Annovar databases, with fields selection, and Annovar options
"annovar": {
   "annotations": {
      "annotation1": {
         "field1": null
         "field2": "field2_renamed",
         ...
      },
      "annotation2": {
         "INFO": null
      },
      ...
   }
   "options": {
      "option1": "value1",
      ...,
   }
}
```

#### annotation::annovar::annotations

List of keywords refering to Annovar databases (see [Annovar Databases documentation](https://annovar.openbioinformatics.org/en/latest/user-guide/download/)), with a list of selected fields for each of them (rename available)

Type: ```str```

Default: ```None```

Examples: 
```json
# Annotation with ClinVar (fields CLNSIG and CLNDN renamed) and Cosmic (all fields)
"annotations": {
   "clinvar_20221231": {
      "CLNSIG": "ClinVar_class"
      "CLNDN": "ClinVar_desease",
      ...
   },
   "cosmic70": {
      "INFO": null
   },
   ...
}
```

#### annotation::annovar::options

List of options available with Annovar tool (see Annovar documentation). As example, these options allows to define splicing threshold or HGVS annotation with refGene database

Examples: 
```json
# HGVS Annotation with refGene (add 'refGene' to 'annotations') and a splicing threshold as 3
"options": {
   "splicing_threshold": 3
   "argument": "'-hgvs'"
   }
}
```

### annotation::snpeff

Annotation process using snpEff tool and options (see [snpEff documentation](https://pcingola.github.io/SnpEff/snpeff/commandline/)).

Examples: 
```json
# Annotation with snpEff databases, with options for HGVS annotation and additional tags.
"snpeff": {
   "options": {
      " -hgvs -noShiftHgvs -spliceSiteSize 3 -lof -oicr "}
   }
}
```

#### annotation::snpeff::options

String (as command line) of options available such as:

 - filters on variants (regions filter, specific changes as intronic or downstream)

 - annotation (e.g. HGVS, loss of function) 

 - database (e.g. only protein coding transcripts, splice sites size)

Examples: 
```json
# Annotation with snpEff databases, with options to generate HGVS annotation, specify to not shift variants according to HGVS notation, define splice sites size to 3, add loss of function (LOF), Nonsense mediated decay and OICR tags.
"options": " -hgvs -noShiftHgvs -spliceSiteSize 3 -lof -oicr "
```

#### annotation::snpeff::stats

HTML file for snpEff stats. Use keyword 'OUTPUT' to generate file according to output file.

Default: ```False```

Examples: 
```json
# Annotation with snpEff databases, and generate a specific stats in HTML format.
"stats": "/path/to/stats.html"
# Annotation with snpEff databases, and generate stats in HTML format associated with output file.
"stats": "OUTPUT.html"
```

#### annotation::snpeff::csvStats

CSV file for snpEff stats. Use keyword 'OUTPUT' to generate file according to output file.

Examples: 
```json
# Annotation with snpEff databases, and generate a specific stats in CSV format.
"csvStats": "/path/to/stats.csv"
# Annotation with snpEff databases, and generate stats in CSV format associated with output file.
"csvStats": "OUTPUT.csv"
```

### annotation::exomiser

Annotation process using Exomiser tool and options (see [Exomiser website documentation](https://www.sanger.ac.uk/tool/exomiser/)).

Examples: 
```json
# Annotation with Exomiser, using database relse '2109', transcripts source as UCSC and a list of HPO terms.
"exomiser": {
   "release": "2109"
   "transcript_source": "refseq"
   "hpo": ['HP:0001156', 'HP:0001363', 'HP:0011304', 'HP:0010055']
}
```

#### annotation::exomiser::release

Release of Exomiser database. This option replace the release variable in 'application.properties' file (see 'exomiser_application_properties' option). The release will be downloaded if it is not available locally. 

Examples: 
```json
# Annotation with release '2109' of Exomiser database.
"release": "2109"
```

### annotation::options

Options for annotations, such as annotation strategy (skip if exists, update, append)

Examples: 
```json
# Annotation with Parquet databases, with update annotation strategy.
"options": {
   "update": True
}
```

#### annotation::options::annotations_update

Update option for annotation (only for Parquet annotation). If True, annotation fields will be removed and re-annotated. These options will be applied to all annotation databases.

Default: ```False```

Examples: 
```json
# Apply update on all annotation fields for all databases.
"update": True
```

#### annotation::options::annotations_append

Append option for annotation (only for Parquet annotation). If True, annotation fields will be annotated only if not annotation exists for the variant. These options will be applied to all annotation databases.

Default: ```False```

Examples: 
```json
# Apply append on all annotation fields for all databases.
"append": True
```

## calculation

Calculation process operations that are defiend in a Calculation Configuration JSON file. List available calculation operations with possible options (see [Calculation JSON file](help.calculation.md) help).

Examples: 
```json
# Calculation with operations for generate variant_id and variant type, extract HGVS from snpEff annotation, select NOMEN from snpEff HGVS with a list of transcripts of preference
"calculation": {
  "variant_id": null,
  "vartype": null,
  "snpeff_hgvs": null,
  "NOMEN": {
    "options": {
      "hgvs_field": "snpeff_hgvs",
      "transcripts": "tests/data/transcripts.tsv"
    }
  }
}
```

## prioritization

Prioritization process use a JSON configuration file that defines all profiles that can be used. By default, all profiles will be calculated from the JSON configuration file, and the first profile will be considered as default. Proritization annotations (INFO/tags) will be generated depending of a input list (default 'PZScore' and 'PZFlag'), for all profiles (e.g. 'PZScore_GERMLINE' for 'GERMLINE' profile) and for default profile (e.g. 'PZScore' for default). Prioritization score mode is 'HOWARD' by default.

Examples: 
```json
# Prioritization with 'GENOME' and 'GERMLINE' (default) profiles, from a list of configured profiles, only 3 prioritization fields returned, and score calculated in 'VaRank' mode
"prioritization": {
  "prioritizations": "config/prioritization_profiles.json",
  "profiles": ["GENOME", "GERMLINE"],
  "default_profile": "GERMLINE",
  "pzfields": ["PZScore", "PZFlag", "PZComment"],
  "prioritization_score_mode": "VaRank"
}
```

### prioritization::prioritizations

Prioritization configuration profiles JSON file defining profiles to calculate. All configured profiles will be calculated by default (see 'profiles' parameter). First profile will be considered as 'default' if none are provided (see 'default_profile' parameter). Default score calculation mode is 'HOWARD'. This option refers to the quick prioritization command line parameter `--prioritizations`.

Type: ```Path```

Default: ```None```

Examples: 
```json
# Prioritization configuration profile JSON file
"prioritizations": "config/prioritization_profiles.json"
```

### prioritization::profiles

Prioritization profiles to consider, from the list of configured profiles. If empty, all configured profiles will be calculated. First profile will be considered as 'default' if none are provided (see 'default_profile' parameter). Prioritization annotations (INFO/tags) will be generated for all these profiles (e.g. 'PZScore_GERMLINE' for 'GERMLINE' profile).

Type: ```str```

Default: ```None```

Examples: 
```json
# Prioritization with 'GERMLINE' profile only
"profiles": ["GERMLINE"]
# Prioritization with 'GENOME' and 'GERMLINE' profiles
"profiles": ["GENOME", "GERMLINE"]
```

### prioritization::default_profile

Prioritization default profile from the list of processed profiles. Prioritization annotations (INFO/tags) will be generated for this default profile (e.g. 'PZScore', 'PZFlags').

Type: ```str```

Default: ```None```

Examples: 
```json
# Prioritization default profile 'GERMLINE'
"default_profile": "GERMLINE"
```

### prioritization::pzfields

Prioritization annotations (INFO/tags) to generate. By default 'PZScore', 'PZFlags'.

Prioritization fields can be selected from:

- PZScore: calculated score from all passing filters, depending of the mode

- PZFlag: final flag ('PASS' or 'FILTERED'), with strategy that consider a variant is filtered as soon as at least one filter do not pass. By default, the variant is considered as 'PASS' (no filter pass)

- PZComment: concatenation of all passing filter comments

- PZTags: combinason of score, flags and comments in a tags format (e.g. 'PZFlag#PASS|PZScore#15|PZComment#Described on ...')

- PZInfos: information about passing filter criteria

Type: ```str```

Default: ```PZScore,PZFlag```

Examples: 
```json
# Prioritization annotations (INFO/tags) list
"pzfields": ["PZScore", "PZFlag", "PZComment"]
```

### prioritization::prioritization_score_mode

Prioritization score can be calculated following multiple mode. The HOWARD mode will increment scores of all passing filters (default). The VaRank mode will select the maximum score from all passing filters.

Type: ```str```

Choices: ```['HOWARD', 'VaRank']```

Default: ```HOWARD```

Examples: 
```json
# Prioritization score calculation mode 'HOWARD'
"prioritization_score_mode": "HOWARD"
# Prioritization score calculation mode 'VaRank'
"prioritization_score_mode": "VaRank"
```

## hgvs

HOWARD annotates variants with HGVS annotation using HUGO HGVS internation Sequence Variant Nomenclature (http://varnomen.hgvs.org/). Annotation refere to refGene and genome to generate HGVS nomenclature for all available transcripts. This annotation add 'hgvs' field into VCF INFO column of a VCF file. Several options are available, to add gene, exon and protein information, to generate a 'full format' detailed annotation, to choose codon format.

Type: ```str```

Default: ```None```

Examples: 
```json
# HGVS annotation  with operations for generate variant_id and variant type, extract HGVS from snpEff annotation, select NOMEN from snpEff HGVS with a list of transcripts of preference
"hgvs": {
  "full_format": true,
  "use_exon": true
}
```

### hgvs::use_gene

Add Gene information to generate HGVS annotation (e.g. 'NM_152232**(TAS1R2)**:c.231T>C').

Default: ```False```

Examples: 
```json
# Use Gene in HGVS annotation
"use_gene": true
```

### hgvs::use_exon

Add Exon information to generate HGVS annotation (e.g. 'NM_152232(exon2):c.231T>C'). Only if 'use_gene' is not enabled.

Default: ```False```

Examples: 
```json
# Use Exon in HGVS annotation
"use_exon": true
```

### hgvs::use_protein

Use Protein level to generate HGVS annotation (e.g. 'NP_689418:p.Cys77Arg'). Can be used with 'use_exon' or 'use_gene'.

Default: ```False```

Examples: 
```json
# Use Protein in HGVS annotation
"use_protein": true
```

### hgvs::add_protein

Add Protein level to DNA HGVS annotation (e.g. 'NM_152232:c.231T>C,NP_689418:p.Cys77Arg').

Default: ```False```

Examples: 
```json
# Add Protein level to DNA HGVS annotation
"add_protein": true
```

### hgvs::full_format

Generates HGVS annotation in a full format (non-standard, e.g. 'TAS1R2:NM_152232:NP_689418:c.231T>C:p.Cys77Arg', 'TAS1R2:NM_152232:NP_689418:exon2:c.231T>C:p.Cys77Arg'). Full format use all information to generates an exhaustive annotation. Use specifically 'use_exon' to add exon information.

Default: ```False```

Examples: 
```json
# Use full format for HGVS annotation
"full_format": true
```

### hgvs::codon_type

Amino Acide Codon format type to use to generate HGVS annotation (default '3'):

- '1': codon in 1 character (e.g. 'C', 'R')

- '3': codon in 3 character (e.g. 'Cys', 'Arg')

- 'FULL': codon in full name (e.g. 'Cysteine', 'Arginine')

Type: ```str```

Choices: ```['1', '3', 'FULL']```

Default: ```3```

Examples: 
```json
# Amino Acide Codon format with 1 character
"codon_type": '1'
# Amino Acide Codon format with 3 character
"codon_type": '3'
# Amino Acide Codon format with full name
"codon_type": 'FULL'
```

### hgvs::refgene

Path to refGene annotation file (see [HOWARD User Guide](user_guide.md#databases-tool)).

Type: ```Path```

Default: ```None```

Examples: 
```json
# Path to refSeq file
"refgene": '~/howard/databases/refseq/current/hg19/ncbiRefSeq.txt'
```

### hgvs::refseqlink

Path to refGeneLink annotation file (see [HOWARD User Guide](user_guide.md#databases-tool)).

Type: ```Path```

Default: ```None```

Examples: 
```json
# Path to refSeq file
"refseqlink": '~/howard/databases/refseq/current/hg19/ncbiRefSeqLink.txt'
```

## stats

Statistics after loading data. 

Default: ```False```

### stats::stats_md

Stats Output file in MarkDown format. 

Type: ```Path```

Default: ```None```

### stats::stats_json

Stats Output file in JSON format. 

Type: ```Path```

Default: ```None```

## query

Query in SQL format (e.g. 'SELECT * FROM variants LIMIT 50'). 

Type: ```str```

Default: ```None```

Examples: 
```
# Simple query to show all variants file
SELECT "#CHROM", POS, REF, ALT, INFO 
FROM variants
```

## explode_infos

Explode VCF INFO/Tag into 'variants' table columns. 

Default: ```False```

## explode_infos_prefix

Explode VCF INFO/Tag with a specific prefix. 

Type: ```str```

Default: ``` ```

## explode_infos_fields

Explode VCF INFO/Tag specific fields/tags. Keyword `*` specify all available fields, except those already specified. Pattern (regex) can be used, such as `.*_score` for fields named with '_score' at the end. Examples:

- 'HGVS,SIFT,Clinvar' (list of fields)

- 'HGVS,*,Clinvar' (list of fields with all other fields at the end)

- 'HGVS,.*_score,Clinvar' (list of 2 fields with all scores in the middle)

- 'HGVS,.*_score,*' (1 field, scores, all other fields)

- 'HGVS,*,.*_score' (1 field, all other fields, all scores) 

Type: ```str```

Default: ```*```

## parquet_partitions

Parquet partitioning using hive (available for any format). This option is faster parallel writing, but memory consuming. Use 'None' (string) for NO partition but split parquet files into a folder. Examples: '#CHROM', '#CHROM,REF', 'None'. 

Type: ```str```

Default: ```None```

## order_by

List of columns to sort the result-set in ascending or descending order. Use SQL format, and keywords ASC (ascending) and DESC (descending). If a column is not available, order will not be considered. Order is enable only for compatible format (e.g. TSV, CSV, JSON). Examples: 'ACMG_score DESC', 'PZFlag DESC, PZScore DESC'. 

Type: ```str```

Default: ``` ```

## query_limit

Limit of number of row for query (only for print result, not output). 

Type: ```int```

Default: ```10```

## query_print_mode

Print mode of query result (only for print result, not output). Either None (native), 'markdown' or 'tabulate'. 

Type: ```str```

Choices: ```[None, 'markdown', 'tabulate']```

Default: ```None```

## threads

Specify the number of threads to use for processing HOWARD. It determines the level of parallelism, either on python scripts, duckdb engine and external tools. It and can help speed up the process/tool. Use -1 to use all available CPU/cores. Either non valid value is 1 CPU/core. 

Type: ```int```

Default: ```-1```

