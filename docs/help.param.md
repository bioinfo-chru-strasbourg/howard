# HOWARD Parameters

HOWARD Parameters JSON file defines parameters to process annotations, calculations, prioritizations, convertions and queries.

## Table of contents

- [HOWARD Parameters](#howard-parameters)
   - [hgvs](#hgvs)
      - [use_gene](#hgvsuse_gene)
      - [use_exon](#hgvsuse_exon)
      - [use_protein](#hgvsuse_protein)
      - [add_protein](#hgvsadd_protein)
      - [full_format](#hgvsfull_format)
      - [codon_type](#hgvscodon_type)
      - [refgene](#hgvsrefgene)
      - [refseqlink](#hgvsrefseqlink)
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
      - [calculations](#calculationcalculations)
      - [calculation_config](#calculationcalculation_config)
   - [prioritization](#prioritization)
      - [prioritizations](#prioritizationprioritizations)
      - [profiles](#prioritizationprofiles)
      - [default_profile](#prioritizationdefault_profile)
      - [pzfields](#prioritizationpzfields)
      - [prioritization_score_mode](#prioritizationprioritization_score_mode)
   - [stats](#stats)
      - [stats_md](#statsstats_md)
      - [stats_json](#statsstats_json)
   - [query](#query)
      - [query](#queryquery)
      - [query_limit](#queryquery_limit)
      - [query_print_mode](#queryquery_print_mode)
   - [export](#export)
      - [order_by](#exportorder_by)
      - [include_header](#exportinclude_header)
      - [parquet_partitions](#exportparquet_partitions)
   - [explode](#explode)
      - [explode_infos](#explodeexplode_infos)
      - [explode_infos_prefix](#explodeexplode_infos_prefix)
      - [explode_infos_fields](#explodeexplode_infos_fields)
   - [threads](#threads)


Examples: 
> Parameters for annotation, calculation, prioritization, HGVS annotation and export

```json
{
   "hgvs": {
      "full_format": true,
      "use_exon": true
   }
   "annotation": {
      "parquet": {
         "annotations": {
            "/path/to/database3.parquet": {
               "field1": null,
               "field2": "field2_renamed"
            },
            "/path/to/database4.vcf.gz": {
               "INFO": null
            },
            "/path/to/database5.bed.gz": {
               "INFO": null
            }
         }
      },
      "bcftools": {
         "annotations": {
            "/path/to/database6.vcf.gz": {
               "field1": null,
               "field2": "field2_renamed"
            },
            "/path/to/database7.bed": {
               "INFO": null
            }
         }
      },
      "annovar": {
         "annotations": {
            "annovar_keyword2": {
               "field1": null,
               "field2": "field2_renamed"
            },
            "annovar_keyword3": {
               "INFO": null
            }
         }
      },
      "snpeff": {
         "options": " -hgvs -noShiftHgvs -spliceSiteSize 3 -lof -oicr "
      },
      "exomiser": {
         "release": "2109",
         "transcript_source": "refseq",
         "hpo": ["HP:0001156", "HP:0001363", "HP:0011304", "HP:0010055"]
      },
      "options": {
         "append": true
      }
   },
   "calculation": {
      "operation1": null,
      "operation2": {
        "options": {
          "option1": "value1",
          "option2": "value2"
        }
      }
   },
   "prioritization": {
      "prioritizations": "config/prioritization_profiles.json",
      "profiles": ["GENOME", "GERMLINE"],
      "default_profile": "GERMLINE",
      "pzfields": ["PZScore", "PZFlag", "PZComment"],
      "prioritization_score_mode": "VaRank"
   },
   "export": {
      "include_header": true
   }
}
```

## hgvs

HOWARD annotates variants with HGVS annotation using HUGO HGVS internation Sequence Variant Nomenclature (http://varnomen.hgvs.org/). Annotation refere to refGene and genome to generate HGVS nomenclature for all available transcripts. This annotation add 'hgvs' field into VCF INFO column of a VCF file. Several options are available, to add gene, exon and protein information, to generate a 'full format' detailed annotation, to choose codon format.

Examples: 

> HGVS annotation  with operations for generate variant_id and variant type, extract HGVS from snpEff annotation, select NOMEN from snpEff HGVS with a list of transcripts of preference

```json
"hgvs": {
  "full_format": true,
  "use_exon": true
}
```

### hgvs::use_gene

Add Gene information to generate HGVS annotation (e.g. 'NM_152232**(TAS1R2)**:c.231T>C').

Default: ```False```

Examples: 

> Use Gene in HGVS annotation

```json
"use_gene": true
```

### hgvs::use_exon

Add Exon information to generate HGVS annotation (e.g. 'NM_152232(exon2):c.231T>C'). Only if 'use_gene' is not enabled.

Default: ```False```

Examples: 

> Use Exon in HGVS annotation

```json
"use_exon": true
```

### hgvs::use_protein

Use Protein level to generate HGVS annotation (e.g. 'NP_689418:p.Cys77Arg'). Can be used with 'use_exon' or 'use_gene'.

Default: ```False```

Examples: 

> Use Protein in HGVS annotation

```json
"use_protein": true
```

### hgvs::add_protein

Add Protein level to DNA HGVS annotation (e.g. 'NM_152232:c.231T>C,NP_689418:p.Cys77Arg').

Default: ```False```

Examples: 

> Add Protein level to DNA HGVS annotation

```json
"add_protein": true
```

### hgvs::full_format

Generates HGVS annotation in a full format (non-standard, e.g. 'TAS1R2:NM_152232:NP_689418:c.231T>C:p.Cys77Arg', 'TAS1R2:NM_152232:NP_689418:exon2:c.231T>C:p.Cys77Arg'). Full format use all information to generates an exhaustive annotation. Use specifically 'use_exon' to add exon information.

Default: ```False```

Examples: 

> Use full format for HGVS annotation

```json
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

> Amino Acide Codon format with 1 character

```json
"codon_type": "1"
```
> Amino Acide Codon format with 1 character# Amino Acide Codon format with 3 character

```json
"codon_type": "3"
```
> Amino Acide Codon format with 1 character# Amino Acide Codon format with 3 character# Amino Acide Codon format with full name

```json
"codon_type": "FULL"
```

### hgvs::refgene

Path to refGene annotation file (see [HOWARD User Guide](user_guide.md#databases-tool)).

Type: ```Path```

Default: ```None```

Examples: 

> Path to refSeq file

```json
"refgene": "~/howard/databases/refseq/current/hg19/ncbiRefSeq.txt"
```

### hgvs::refseqlink

Path to refGeneLink annotation file (see [HOWARD User Guide](user_guide.md#databases-tool)).

Type: ```Path```

Default: ```None```

Examples: 

> Path to refSeq file

```json
"refseqlink": "~/howard/databases/refseq/current/hg19/ncbiRefSeqLink.txt"
```

## annotation

Annotation process using HOWARD algorithms or external tools.


For HOWARD Parquet algorithm, provide the list of database files available (formats such as Parquet, VCF, TSV, duckDB, JSON) and select fields (rename possible, 'INFO' keyword for all fields), or use 'ALL' keyword to detect available databases.


For external tools, such as Annovar, snpEff and Exomiser, specify parameters such as annotation keywords (Annovar) and options (depending on the tool), and select fields (BCFtools and Annovar, field rename available).

Examples: 

> Annotation with multiple tools in multiple formats with multiple options

```json
"annotation": {
   "parquet": {
      "annotations": {
         "/path/to/database3.parquet": {
            "field1": null,
            "field2": "field2_renamed"
         },
         "/path/to/database4.vcf.gz": {
            "INFO": null
         }
         "/path/to/database5.bed.gz": {
            "INFO": null
         }
      }
   }
   "bcftools": {
      "annotations": {
         "/path/to/database6.vcf.gz": {
            "field1": null,
            "field2": "field2_renamed"
         },
         "/path/to/database7.bed": {
            "INFO": null
         }
      }
   }
   "annovar": {
      "annotations": {
         "annovar_keyword2": {
            "field1": null,
            "field2": "field2_renamed"
         },
         "annovar_keyword3": {
            "INFO": null
         }
      }
   }
   "snpeff": {
      "options": {
         " -hgvs -noShiftHgvs -spliceSiteSize 3 -lof -oicr "}
      }
   }
   "exomiser": {
      "release": "2109",
      "transcript_source": "refseq",
      "hpo": ["HP:0001156", "HP:0001363", "HP:0011304", "HP:0010055"]
   }
   "options": {
      "append": true
   }
}
```

### annotation::parquet

Annotation process using HOWARD Parquet algorithm, for the list of databases available (formats such as Parquet, VCF, TSV, duckDB, JSON).

Examples: 
> Annotation with multiple databases in multiple formats

```json
"parquet": {
   "annotations": {
      "/path/to/database3.parquet": {
         "field1": null,
         "field2": "field2_renamed",
      },
      "/path/to/database4.vcf.gz": {
         "INFO": null
      }
      "/path/to/database5.bed.gz": {
         "INFO": null
      }
   }
}
```

#### annotation::parquet::annotations

Specify a list of databases files available (formats such as Parquet, VCF, TSV, duckDB, JSON). This parameter enables users to select specific database fields and optionally rename them (e.g. '"field": null' to keep field name, '"field": "new_name"' to rename field). Use 'INFO' keyword to select all fields within the database INFO/Tags header (e.g. '"INFO": null'). Use 'ALL' keyword to select all fields within the database regardless INFO/Tags header (e.g. '"ALL": null').


For add all availalbe databases files, use 'ALL' keyword, with filters on type and release (e.g. 'ALL', 'ALL:parquet:current', 'ALL:parquet,vcf:current,devel').


If a full path is not provided, the system will automatically detect files within database folders (see Configuration doc) and assembly (see Parameter option).

Examples: 
> Annotation with dbSNP database with INFO/tags fields, and dbNSFP databases with all fields

```json
"annotations": {
   "tests/databases/annotations/current/hg19/avsnp150.parquet": {
      "INFO": null
   }
   "tests/databases/annotations/current/hg19/dbnsfp42a.parquet": {
      "ALL": null
   }
}
```
> Annotation with dbNSFP (only PolyPhen, ClinVar and REVEL score), and rename fields

```json
"annotations": {
   "tests/databases/annotations/current/hg19/dbnsfp42a.parquet": {
      "Polyphen2_HDIV_pred": "PolyPhen",
      "ClinPred_pred": "ClinVar",
      "REVEL_score": null
   }
}

```
> Annotation with refSeq as a BED file

```json
"annotations": {
   "tests/databases/annotations/current/hg19/refGene.bed": {
      "INFO": null
   }
}
```
> Annotation with dbNSFP REVEL annotation (as a VCF file) within configured annotation databases folders (default: '~/howard/databases/annotations/current') and assembly (default: 'hg19')

```json
"annotations": {
   "dbnsfp42a.REVEL.vcf.gz": {
      "REVEL_score": null,
      "REVEL_rankscore": null
   }
}

```
> Annotation with all available databases in Parquet for ''current release

```json
"parquet": {
   "annotations": {
      "ALL": {
         "formats": ["parquet"],
         "releases": ["current"]
      }
   }
}
```

### annotation::bcftools

Annotation process using BCFTools. Provide a list of database files and annotation fields.

Examples: 

> Annotation with multiple databases in multiple formats

```json
"parquet": {
   "bcftools": {
      "/path/to/database1.vcf.gz": {
         "field1": null
         "field2": "field2_renamed",
      },
      "database2.bed.gz": {
         "INFO": null
      }
   }
}
```

#### annotation::bcftools::annotations

Specify the list of database files in formats VCF or BED. Files need to be compressed and indexed.

This parameter enables users to select specific database fields and optionally rename them (e.g. '"field": null' to keep field name, '"field": "new_name"' to rename field). Use 'INFO' or 'ALL' keyword to select all fields within the database INFO/Tags header (e.g. '"INFO": null', '"ALL": null').


If a full path is not provided, the system will automatically detect files within database folders (see Configuration doc) and assembly (see Parameter option).

Examples: 

> Annotation with dbSNP database  with all fields

```json
"annotations": {
   "tests/databases/annotations/current/hg19/avsnp150.vcf.gz": {
      "INFO": null
   }
}
```
> Annotation with dbSNP database  with all fields# Annotation with dbNSFP (only PolyPhen, ClinVar and REVEL score), and rename fields

```json
"annotations": {
   "tests/databases/annotations/current/hg19/dbnsfp42a.vcf.gz": {
      "Polyphen2_HDIV_pred": "PolyPhen",
      "ClinPred_pred": "ClinVar",
      "REVEL_score": null
   }
}
```
> Annotation with dbSNP database  with all fields# Annotation with dbNSFP (only PolyPhen, ClinVar and REVEL score), and rename fields# Annotation with refSeq as a BED file

```json
"annotations": {
   "tests/databases/annotations/current/hg19/refGene.bed": {
      "INFO": null
   }
}
```
> Annotation with dbSNP database  with all fields# Annotation with dbNSFP (only PolyPhen, ClinVar and REVEL score), and rename fields# Annotation with refSeq as a BED file# Annotation with dbNSFP REVEL annotation (as a VCF file) within configured annotation databases folders (default: '~/howard/databases/annotations/current') and assembly (default: 'hg19')

```json
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

> Annotation with multiple Annovar databases, with fields selection, and Annovar options

```json
"annovar": {
   "annotations": {
      "annovar_keyword2": {
         "field1": null,
         "field2": "field2_renamed",
      },
      "annovar_keyword3": {
         "INFO": null
      }
   }
}
```

#### annotation::annovar::annotations

List of keywords refering to Annovar databases (see [Annovar Databases documentation](https://annovar.openbioinformatics.org/en/latest/user-guide/download/)), with a list of selected fields for each of them (rename available)

Examples: 

> Annotation with ClinVar (fields CLNSIG and CLNDN renamed) and Cosmic (all fields)

```json
"annotations": {
   "clinvar_20221231": {
      "CLNSIG": "ClinVar_class"
      "CLNDN": "ClinVar_desease",
   },
   "cosmic70": {
      "INFO": null
   },
}
```

#### annotation::annovar::options

List of options available with Annovar tool (see Annovar documentation). As example, these options allows to define splicing threshold or HGVS annotation with refGene database

Examples: 

> HGVS Annotation with refGene (add 'refGene' to 'annotations') and a splicing threshold as 3

```json
"options": {
   "splicing_threshold": 3,
   "argument": "'-hgvs'"
   }
}
```

### annotation::snpeff

Annotation process using snpEff tool and options (see [snpEff documentation](https://pcingola.github.io/SnpEff/snpeff/commandline/)).

Examples: 

> Annotation with snpEff databases, with options for HGVS annotation and additional tags.

```json
"snpeff": {
   "options": {
      " -hgvs -noShiftHgvs -spliceSiteSize 3 -lof -oicr "
   }
}
```

#### annotation::snpeff::options

String (as command line) of options available such as:

 - filters on variants (regions filter, specific changes as intronic or downstream)

 - annotation (e.g. HGVS, loss of function) 

 - database (e.g. only protein coding transcripts, splice sites size)

Examples: 

> Annotation with snpEff databases, with options to generate HGVS annotation, specify to not shift variants according to HGVS notation, define splice sites size to 3, add loss of function (LOF), Nonsense mediated decay and OICR tags.

```json
"options": " -hgvs -noShiftHgvs -spliceSiteSize 3 -lof -oicr "
```

#### annotation::snpeff::stats

HTML file for snpEff stats. Use keyword 'OUTPUT' to generate file according to output file.

Examples: 

> Annotation with snpEff databases, and generate a specific stats in HTML format.

```json
"stats": "/path/to/stats.html"
```
> Annotation with snpEff databases, and generate a specific stats in HTML format.# Annotation with snpEff databases, and generate stats in HTML format associated with output file.

```json
"stats": "OUTPUT.html"
```

#### annotation::snpeff::csvStats

CSV file for snpEff stats. Use keyword 'OUTPUT' to generate file according to output file.

Examples: 

> Annotation with snpEff databases, and generate a specific stats in CSV format.

```json
"csvStats": "/path/to/stats.csv"
```
> Annotation with snpEff databases, and generate a specific stats in CSV format.# Annotation with snpEff databases, and generate stats in CSV format associated with output file.

```json
"csvStats": "OUTPUT.csv"
```

### annotation::exomiser

Annotation process using Exomiser tool and options (see [Exomiser website documentation](https://www.sanger.ac.uk/tool/exomiser/)).

Examples: 

> Annotation with Exomiser, using database relse '2109', transcripts source as UCSC and a list of HPO terms.

```json
"exomiser": {
   "release": "2109"
   "transcript_source": "refseq"
   "hpo": ["HP:0001156", "HP:0001363", "HP:0011304", "HP:0010055"]
}
```

#### annotation::exomiser::release

Release of Exomiser database. This option replace the release variable in 'application.properties' file (see 'exomiser_application_properties' option). The release will be downloaded if it is not available locally. 

Examples: 

> Annotation with release '2109' of Exomiser database.

```json
"release": "2109"
```

### annotation::options

Options for annotations, such as annotation strategy (skip if exists, update, append)

Examples: 

> Annotation with Parquet databases, with update annotation strategy.

```json
"options": {
   "update": true
}
```

#### annotation::options::annotations_update

Update option for annotation (only for Parquet annotation). If True, annotation fields will be removed and re-annotated. These options will be applied to all annotation databases.

Default: ```False```

Examples: 

> Apply update on all annotation fields for all databases.

```json
"update": true
```

#### annotation::options::annotations_append

Append option for annotation (only for Parquet annotation). If True, annotation fields will be annotated only if not annotation exists for the variant. These options will be applied to all annotation databases.

Default: ```False```

Examples: 

> Apply append on all annotation fields for all databases.

```json
"append": true
```

## calculation

Calculation process operations that are defiend in a Calculation Configuration JSON file. List available calculation operations with possible options (see [Calculation JSON file](help.calculation.md) help).

Examples: 

> Calculation of operations 'operation1' and 'operation2' (with options) defined in 'calculation_config.json' file

```json
"calculation": {
  "calculations": {
    "operation1": null,
    "operation2": {
      "options": {
        "option1": "value1",
        "option2": "value2"
      }
    }
  },
  "calculation_config": "calculation_config.json"
}
```

### calculation::calculations

List of operations to process with possible options (see [Calculation JSON file](help.calculation.md) help).

Examples: 

> Calculation with operations for generate variant_id and variant type, extract HGVS from snpEff annotation, select NOMEN from snpEff HGVS with a list of transcripts of preference

```json
"calculations": {
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

### calculation::calculation_config

Calculation configuration JSON file. 

Type: ```Path```

Default: ```None```

Examples: 
> Calculation configuration JSON file as an option

```json
"calculation_config": "calculation_config.json" 
```

## prioritization

Prioritization process use a JSON configuration file that defines all profiles that can be used. By default, all profiles will be calculated from the JSON configuration file, and the first profile will be considered as default. Proritization annotations (INFO/tags) will be generated depending of a input list (default 'PZScore' and 'PZFlag'), for all profiles (e.g. 'PZScore_GERMLINE' for 'GERMLINE' profile) and for default profile (e.g. 'PZScore' for default). Prioritization score mode is 'HOWARD' by default.

Examples: 

> Prioritization with 'GENOME' and 'GERMLINE' (default) profiles, from a list of configured profiles, only 3 prioritization fields returned, and score calculated in 'VaRank' mode

```json
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

Type: ```str```

Default: ```None```

Examples: 

> Prioritization configuration profile JSON file

```json
"prioritizations": "config/prioritization_profiles.json"
```

### prioritization::profiles

Prioritization profiles to consider, from the list of configured profiles. If empty, all configured profiles will be calculated. First profile will be considered as 'default' if none are provided (see 'default_profile' parameter). Prioritization annotations (INFO/tags) will be generated for all these profiles (e.g. 'PZScore_GERMLINE' for 'GERMLINE' profile).

Type: ```str```

Default: ```None```

Examples: 

> Prioritization with 'GERMLINE' profile only

```json
"profiles": ["GERMLINE"]
```
> Prioritization with 'GERMLINE' profile only# Prioritization with 'GENOME' and 'GERMLINE' profiles

```json
"profiles": ["GENOME", "GERMLINE"]
```

### prioritization::default_profile

Prioritization default profile from the list of processed profiles. Prioritization annotations (INFO/tags) will be generated for this default profile (e.g. 'PZScore', 'PZFlags').

Type: ```str```

Default: ```None```

Examples: 

> Prioritization default profile 'GERMLINE'

```json
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

> Prioritization annotations (INFO/tags) list

```json
"pzfields": ["PZScore", "PZFlag", "PZComment"]
```

### prioritization::prioritization_score_mode

Prioritization score can be calculated following multiple mode. The HOWARD mode will increment scores of all passing filters (default). The VaRank mode will select the maximum score from all passing filters.

Type: ```str```

Choices: ```['HOWARD', 'VaRank']```

Default: ```HOWARD```

Examples: 

> Prioritization score calculation mode 'HOWARD'

```json
"prioritization_score_mode": "HOWARD"
```
> Prioritization score calculation mode 'HOWARD'# Prioritization score calculation mode 'VaRank'

```json
"prioritization_score_mode": "VaRank"
```

## stats

Statistics on loaded variants.

### stats::stats_md

Stats Output file in MarkDown format. 

Type: ```Path```

Default: ```None```

Examples: 
> Export statistics in Markdown format

```json
"stats_md": "/tmp/stats.md" 
```

### stats::stats_json

Stats Output file in JSON format. 

Type: ```Path```

Default: ```None```

Examples: 
> Export statistics in JSON format

```json
"stats_json": "/tmp/stats.json" 
```

## query

Query options tools. Mainly a SQL query, based on 'variants' table corresponding on input file data, or a independant query. Print options for 'query' tool allow limiting number of lines and choose printing mode.

Type: ```str```

Default: ```None```

### query::query

Query in SQL format (e.g. 'SELECT * FROM variants LIMIT 50'). 

Type: ```str```

Default: ```None```

Examples: 

> Simple query to show all variants file

```
SELECT "#CHROM", POS, REF, ALT, INFO 
FROM variants
```

### query::query_limit

Limit of number of row for query (only for print result, not output). 

Type: ```int```

Default: ```10```

### query::query_print_mode

Print mode of query result (only for print result, not output). Either None (native), 'markdown' or 'tabulate'. 

Type: ```str```

Choices: ```[None, 'markdown', 'tabulate']```

Default: ```None```

## export

Export options for output files, such as data order, include header in output and hive partitioning.

### export::order_by

List of columns to sort the result-set in ascending or descending order. Use SQL format, and keywords ASC (ascending) and DESC (descending). If a column is not available, order will not be considered. Order is enable only for compatible format (e.g. TSV, CSV, JSON). Examples: 'ACMG_score DESC', 'PZFlag DESC, PZScore DESC'. 

Type: ```str```

Default: ``` ```

Examples: 
> Order by ACMG score in descending order

```json
"order_by": "ACMG_score DESC" 
```
> Order by PZFlag and PZScore in descending order

```json
"order_by": PZFlag DESC, PZScore DESC" 
```

### export::include_header

Include header (in VCF format) in output file. Only for compatible formats (tab-delimiter format as TSV or BED). 

Default: ```False```

### export::parquet_partitions

Parquet partitioning using hive (available for any format). This option is faster parallel writing, but memory consuming. Use 'None' (string) for NO partition but split parquet files into a folder. Examples: '#CHROM', '#CHROM,REF', 'None'. 

Type: ```str```

Default: ```None```

## explode

Explode options for INFO/tags annotations within VCF files.

### explode::explode_infos

Explode VCF INFO/Tag into 'variants' table columns. 

Default: ```False```

### explode::explode_infos_prefix

Explode VCF INFO/Tag with a specific prefix. 

Type: ```str```

Default: ``` ```

### explode::explode_infos_fields

Explode VCF INFO/Tag specific fields/tags. Keyword `*` specify all available fields, except those already specified. Pattern (regex) can be used, such as `.*_score` for fields named with '_score' at the end. Examples:

- 'HGVS,SIFT,Clinvar' (list of fields)

- 'HGVS,*,Clinvar' (list of fields with all other fields at the end)

- 'HGVS,.*_score,Clinvar' (list of 2 fields with all scores in the middle)

- 'HGVS,.*_score,*' (1 field, scores, all other fields)

- 'HGVS,*,.*_score' (1 field, all other fields, all scores) 

Type: ```str```

Default: ```*```

## threads

Specify the number of threads to use for processing HOWARD. It determines the level of parallelism, either on python scripts, duckdb engine and external tools. It and can help speed up the process/tool. Use -1 to use all available CPU/cores. Either non valid value is 1 CPU/core. 

Type: ```int```

Default: ```-1```

Examples: 
> Automatically detect all available CPU/cores

```json
"threads": -1
```
> Define 8 CPU/cores

```json
"threads": 8
```

