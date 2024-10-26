---
title: HOWARD Help Parameters
---

- [<span class="toc-section-number">1</span>
  Introduction](#introduction)
- [<span class="toc-section-number">2</span> hgvs](#hgvs)
  - [<span class="toc-section-number">2.1</span> use_gene](#use_gene)
  - [<span class="toc-section-number">2.2</span> use_exon](#use_exon)
  - [<span class="toc-section-number">2.3</span>
    use_protein](#use_protein)
  - [<span class="toc-section-number">2.4</span>
    add_protein](#add_protein)
  - [<span class="toc-section-number">2.5</span>
    full_format](#full_format)
  - [<span class="toc-section-number">2.6</span>
    codon_type](#codon_type)
  - [<span class="toc-section-number">2.7</span> refgene](#refgene)
  - [<span class="toc-section-number">2.8</span>
    refseqlink](#refseqlink)
- [<span class="toc-section-number">3</span> annotation](#annotation)
  - [<span class="toc-section-number">3.1</span> parquet](#parquet)
    - [<span class="toc-section-number">3.1.1</span>
      annotations](#annotations)
  - [<span class="toc-section-number">3.2</span> bcftools](#bcftools)
    - [<span class="toc-section-number">3.2.1</span>
      annotations](#annotations-1)
  - [<span class="toc-section-number">3.3</span> annovar](#annovar)
    - [<span class="toc-section-number">3.3.1</span>
      annotations](#annotations-2)
    - [<span class="toc-section-number">3.3.2</span> options](#options)
  - [<span class="toc-section-number">3.4</span> snpeff](#snpeff)
    - [<span class="toc-section-number">3.4.1</span>
      options](#options-1)
    - [<span class="toc-section-number">3.4.2</span> stats](#stats)
    - [<span class="toc-section-number">3.4.3</span>
      csvStats](#csvstats)
  - [<span class="toc-section-number">3.5</span> snpsift](#snpsift)
    - [<span class="toc-section-number">3.5.1</span>
      annotations](#annotations-3)
  - [<span class="toc-section-number">3.6</span> bigwig](#bigwig)
    - [<span class="toc-section-number">3.6.1</span>
      annotations](#annotations-4)
  - [<span class="toc-section-number">3.7</span> exomiser](#exomiser)
    - [<span class="toc-section-number">3.7.1</span> release](#release)
    - [<span class="toc-section-number">3.7.2</span>
      transcript_source](#transcript_source)
    - [<span class="toc-section-number">3.7.3</span> hpo](#hpo)
  - [<span class="toc-section-number">3.8</span> splice](#splice)
    - [<span class="toc-section-number">3.8.1</span>
      split_mode](#split_mode)
    - [<span class="toc-section-number">3.8.2</span>
      spliceai_distance](#spliceai_distance)
    - [<span class="toc-section-number">3.8.3</span>
      spliceai_mask](#spliceai_mask)
    - [<span class="toc-section-number">3.8.4</span>
      transcript](#transcript)
    - [<span class="toc-section-number">3.8.5</span> rm_snps](#rm_snps)
    - [<span class="toc-section-number">3.8.6</span>
      rm_annot](#rm_annot)
    - [<span class="toc-section-number">3.8.7</span>
      whitespace](#whitespace)
  - [<span class="toc-section-number">3.9</span> options](#options-2)
    - [<span class="toc-section-number">3.9.1</span>
      annotations_update](#annotations_update)
    - [<span class="toc-section-number">3.9.2</span>
      annotations_append](#annotations_append)
- [<span class="toc-section-number">4</span> calculation](#calculation)
  - [<span class="toc-section-number">4.1</span>
    calculations](#calculations)
  - [<span class="toc-section-number">4.2</span>
    calculation_config](#calculation_config)
- [<span class="toc-section-number">5</span>
  prioritization](#prioritization)
  - [<span class="toc-section-number">5.1</span>
    prioritizations](#prioritizations)
  - [<span class="toc-section-number">5.2</span> profiles](#profiles)
  - [<span class="toc-section-number">5.3</span>
    default_profile](#default_profile)
  - [<span class="toc-section-number">5.4</span> pzfields](#pzfields)
  - [<span class="toc-section-number">5.5</span>
    prioritization_score_mode](#prioritization_score_mode)
  - [<span class="toc-section-number">5.6</span> pzprefix](#pzprefix)
- [<span class="toc-section-number">6</span> stats](#stats-1)
  - [<span class="toc-section-number">6.1</span> stats_md](#stats_md)
  - [<span class="toc-section-number">6.2</span>
    stats_json](#stats_json)
- [<span class="toc-section-number">7</span> query](#query)
  - [<span class="toc-section-number">7.1</span> query](#query-1)
  - [<span class="toc-section-number">7.2</span>
    query_limit](#query_limit)
  - [<span class="toc-section-number">7.3</span>
    query_print_mode](#query_print_mode)
- [<span class="toc-section-number">8</span> export](#export)
  - [<span class="toc-section-number">8.1</span> order_by](#order_by)
  - [<span class="toc-section-number">8.2</span>
    include_header](#include_header)
  - [<span class="toc-section-number">8.3</span>
    parquet_partitions](#parquet_partitions)
- [<span class="toc-section-number">9</span> explode](#explode)
  - [<span class="toc-section-number">9.1</span>
    explode_infos](#explode_infos)
  - [<span class="toc-section-number">9.2</span>
    explode_infos_prefix](#explode_infos_prefix)
  - [<span class="toc-section-number">9.3</span>
    explode_infos_fields](#explode_infos_fields)
- [<span class="toc-section-number">10</span> transcripts](#transcripts)
  - [<span class="toc-section-number">10.1</span> table](#table)
  - [<span class="toc-section-number">10.2</span>
    transcripts_info_field_json](#transcripts_info_field_json)
  - [<span class="toc-section-number">10.3</span>
    transcripts_info_field_format](#transcripts_info_field_format)
  - [<span class="toc-section-number">10.4</span>
    transcripts_info_json](#transcripts_info_json)
  - [<span class="toc-section-number">10.5</span>
    transcripts_info_format](#transcripts_info_format)
  - [<span class="toc-section-number">10.6</span>
    transcript_id_remove_version](#transcript_id_remove_version)
  - [<span class="toc-section-number">10.7</span>
    transcript_id_mapping_file](#transcript_id_mapping_file)
  - [<span class="toc-section-number">10.8</span> Example of transcript
    ID mapping file](#example-of-transcript-id-mapping-file)
  - [<span class="toc-section-number">10.9</span>
    transcript_id_mapping_force](#transcript_id_mapping_force)
  - [<span class="toc-section-number">10.10</span> struct](#struct)
    - [<span class="toc-section-number">10.10.1</span>
      from_column_format](#from_column_format)
    - [<span class="toc-section-number">10.10.2</span>
      from_columns_map](#from_columns_map)
    - [<span class="toc-section-number">10.10.3</span> commons
      parameters](#commons-parameters)
  - [<span class="toc-section-number">10.11</span>
    prioritization](#prioritization-1)
    - [<span class="toc-section-number">10.11.1</span>
      profiles](#profiles-1)
    - [<span class="toc-section-number">10.11.2</span>
      default_profile](#default_profile-1)
    - [<span class="toc-section-number">10.11.3</span>
      prioritization_config](#prioritization_config)
    - [<span class="toc-section-number">10.11.4</span>
      prioritization_score_mode](#prioritization_score_mode-1)
    - [<span class="toc-section-number">10.11.5</span>
      pzprefix](#pzprefix-1)
    - [<span class="toc-section-number">10.11.6</span>
      pzfields](#pzfields-1)
    - [<span class="toc-section-number">10.11.7</span>
      prioritization_transcripts_order](#prioritization_transcripts_order)
    - [<span class="toc-section-number">10.11.8</span>
      prioritization_transcripts](#prioritization_transcripts)
    - [<span class="toc-section-number">10.11.9</span>
      prioritization_transcripts_force](#prioritization_transcripts_force)
    - [<span class="toc-section-number">10.11.10</span>
      prioritization_transcripts_version_force](#prioritization_transcripts_version_force)
  - [<span class="toc-section-number">10.12</span> export](#export-1)
- [<span class="toc-section-number">11</span> threads](#threads)
- [<span class="toc-section-number">12</span> samples](#samples)
  - [<span class="toc-section-number">12.1</span> list](#list)
  - [<span class="toc-section-number">12.2</span> check](#check)
- [<span class="toc-section-number">13</span> databases](#databases)

# Introduction

HOWARD Parameters JSON file defines parameters to process annotations,
calculations, prioritizations, convertions and queries.

Examples:

> Parameters for annotation, calculation, prioritization, HGVS
> annotation and export

> ``` json
> {
>    "hgvs": {
>       "full_format": true,
>       "use_exon": true
>    }
>    "annotation": {
>       "parquet": {
>          "annotations": {
>             "/path/to/database3.parquet": {
>                "field1": null,
>                "field2": "field2_renamed"
>             },
>             "/path/to/database4.vcf.gz": {
>                "INFO": null
>             },
>             "/path/to/database5.bed.gz": {
>                "INFO": null
>             }
>          }
>       },
>       "bcftools": {
>          "annotations": {
>             "/path/to/database6.vcf.gz": {
>                "field1": null,
>                "field2": "field2_renamed"
>             },
>             "/path/to/database7.bed": {
>                "INFO": null
>             }
>          }
>       },
>       "annovar": {
>          "annotations": {
>             "annovar_keyword2": {
>                "field1": null,
>                "field2": "field2_renamed"
>             },
>             "annovar_keyword3": {
>                "INFO": null
>             }
>          }
>       },
>       "snpeff": {
>          "options": " -hgvs -noShiftHgvs -spliceSiteSize 3 -lof -oicr "
>       },
>       "exomiser": {
>          "release": "2109",
>          "transcript_source": "refseq",
>          "hpo": ["HP:0001156", "HP:0001363", "HP:0011304", "HP:0010055"]
>       },
>       "options": {
>          "append": true
>       }
>    },
>    "calculation": {
>       "operation1": null,
>       "operation2": {
>         "options": {
>           "option1": "value1",
>           "option2": "value2"
>         }
>       }
>    },
>    "prioritization": {
>       "prioritizations": "config/prioritization_profiles.json",
>       "profiles": ["GENOME", "GERMLINE"],
>       "default_profile": "GERMLINE",
>       "pzfields": ["PZScore", "PZFlag", "PZComment"],
>       "prioritization_score_mode": "VaRank"
>    },
>    "export": {
>       "include_header": true
>    }
> }
> ```

# hgvs

HOWARD annotates variants with HGVS annotation using HUGO HGVS
internation Sequence Variant Nomenclature (http://varnomen.hgvs.org/).
Annotation refere to refGene and genome to generate HGVS nomenclature
for all available transcripts. This annotation add 'hgvs' field into VCF
INFO column of a VCF file. Several options are available, to add gene,
exon and protein information, to generate a 'full format' detailed
annotation, to choose codon format.

Examples:

> HGVS annotation with operations for generate variant_id and variant
> type, extract HGVS from snpEff annotation, select NOMEN from snpEff
> HGVS with a list of transcripts of preference

> ``` json
> {
>    "hgvs": {
>      "full_format": true,
>      "use_exon": true
>    }
> }
> ```

## use_gene

Add Gene information to generate HGVS annotation (e.g.
'NM_152232**(TAS1R2)**:c.231T\>C').

Default: `False`

Examples:

> Use Gene in HGVS annotation

> ``` json
> {
>    "use_gene": true
> }
> ```

## use_exon

Add Exon information to generate HGVS annotation (e.g.
'NM_152232(exon2):c.231T\>C'). Only if 'use_gene' is not enabled.

Default: `False`

Examples:

> Use Exon in HGVS annotation

> ``` json
> {
>    "use_exon": true
> }
> ```

## use_protein

Use Protein level to generate HGVS annotation (e.g.
'NP_689418:p.Cys77Arg'). Can be used with 'use_exon' or 'use_gene'.

Default: `False`

Examples:

> Use Protein in HGVS annotation

> ``` json
> {
>    "use_protein": true
> }
> ```

## add_protein

Add Protein level to DNA HGVS annotation (e.g.
'NM_152232:c.231T\>C,NP_689418:p.Cys77Arg').

Default: `False`

Examples:

> Add Protein level to DNA HGVS annotation

> ``` json
> {
>    "add_protein": true
> }
> ```

## full_format

Generates HGVS annotation in a full format (non-standard, e.g.
'TAS1R2:NM_152232:NP_689418:c.231T\>C:p.Cys77Arg',
'TAS1R2:NM_152232:NP_689418:exon2:c.231T\>C:p.Cys77Arg'). Full format
use all information to generates an exhaustive annotation. Use
specifically 'use_exon' to add exon information.

Default: `False`

Examples:

> Use full format for HGVS annotation

> ``` json
> {
>    "full_format": true
> }
> ```

## codon_type

Amino Acide Codon format type to use to generate HGVS annotation
(default '3'):

- '1': codon in 1 character (e.g. 'C', 'R')

- '3': codon in 3 character (e.g. 'Cys', 'Arg')

- 'FULL': codon in full name (e.g. 'Cysteine', 'Arginine')

Type: `str`

Choices: `['1', '3', 'FULL']`

Default: `3`

Examples:

> Amino Acide Codon format with 1 character

> ``` json
> {
>    "codon_type": "1"
> }
> ```

> Amino Acide Codon format with 3 character

> ``` json
> {
>    "codon_type": "3"
> }
> ```

> Amino Acide Codon format with full name

> ``` json
> {
>    "codon_type": "FULL"
> }
> ```

## refgene

Path to refGene annotation file (see [HOWARD User
Guide](user_guide.md#databases-tool)).

Type: `Path`

Default: `None`

Examples:

> Path to refSeq file

> ``` json
> {
>    "refgene": "~/howard/databases/refseq/current/hg19/ncbiRefSeq.txt"
> }
> ```

## refseqlink

Path to refGeneLink annotation file (see [HOWARD User
Guide](user_guide.md#databases-tool)).

Type: `Path`

Default: `None`

Examples:

> Path to refSeq file

> ``` json
> {
>    "refseqlink": "~/howard/databases/refseq/current/hg19/ncbiRefSeqLink.txt"
> }
> ```

# annotation

Annotation process using HOWARD algorithms or external tools.

For HOWARD Parquet algorithm, provide the list of database files
available (formats such as Parquet, VCF, TSV, duckDB, JSON) and select
fields (rename possible, 'INFO' keyword for all fields), or use 'ALL'
keyword to detect available databases.

For external tools, such as Annovar, snpEff and Exomiser, specify
parameters such as annotation keywords (Annovar) and options (depending
on the tool), and select fields (BCFtools and Annovar, field rename
available).

Examples:

> Annotation with multiple tools in multiple formats with multiple
> options

> ``` json
> {
>    "annotation": {
>       "parquet": {
>          "annotations": {
>             "/path/to/database3.parquet": {
>                "field1": null,
>                "field2": "field2_renamed"
>             },
>             "/path/to/database4.vcf.gz": {
>                "INFO": null
>             }
>             "/path/to/database5.bed.gz": {
>                "INFO": null
>             }
>          }
>       }
>       "bcftools": {
>          "annotations": {
>             "/path/to/database6.vcf.gz": {
>                "field1": null,
>                "field2": "field2_renamed"
>             },
>             "/path/to/database7.bed": {
>                "INFO": null
>             }
>          }
>       }
>       "annovar": {
>          "annotations": {
>             "annovar_keyword2": {
>                "field1": null,
>                "field2": "field2_renamed"
>             },
>             "annovar_keyword3": {
>                "INFO": null
>             }
>          }
>       }
>       "snpeff": {
>          "options": " -hgvs -noShiftHgvs -spliceSiteSize 3 -lof -oicr "
>       }
>       "exomiser": {
>          "release": "2109",
>          "transcript_source": "refseq",
>          "hpo": ["HP:0001156", "HP:0001363", "HP:0011304", "HP:0010055"]
>       }
>       "options": {
>          "append": true
>       }
>    }
> }
> ```

## parquet

Annotation process using HOWARD Parquet algorithm, for the list of
databases available (formats such as Parquet, VCF, TSV, duckDB, JSON).

Examples:

> Annotation with multiple databases in multiple formats

> ``` json
> {
>    "parquet": {
>       "annotations": {
>          "/path/to/database3.parquet": {
>             "field1": null,
>             "field2": "field2_renamed",
>          },
>          "/path/to/database4.vcf.gz": {
>             "INFO": null
>          }
>          "/path/to/database5.bed.gz": {
>             "INFO": null
>          }
>       }
>    }
> }
> ```

### annotations

Specify a list of databases files available (formats such as Parquet,
VCF, TSV, duckDB, JSON). This parameter enables users to select specific
database fields and optionally rename them (e.g. '"field": null' to keep
field name, '"field": "new_name"' to rename field). Use 'INFO' keyword
to select all fields within the database INFO/Tags header (e.g. '"INFO":
null'). Use 'ALL' keyword to select all fields within the database
regardless INFO/Tags header (e.g. '"ALL": null').

For add all availalbe databases files, use 'ALL' keyword, with filters
on type and release (e.g. 'ALL', 'ALL:parquet:current',
'ALL:parquet,vcf:current,devel').

If a full path is not provided, the system will automatically detect
files within database folders (see Configuration doc) and assembly (see
Parameter option).

Examples:

> Annotation with dbSNP database with INFO/tags fields, and dbNSFP
> databases with all fields

> ``` json
> {
>    "annotations": {
>       "tests/databases/annotations/current/hg19/avsnp150.parquet": {
>          "INFO": null
>       }
>       "tests/databases/annotations/current/hg19/dbnsfp42a.parquet": {
>          "ALL": null
>       }
>    }
> }
> ```

> Annotation with dbNSFP (only PolyPhen, ClinVar and REVEL score), and
> rename fields

> ``` json
> {
>    "annotations": {
>       "tests/databases/annotations/current/hg19/dbnsfp42a.parquet": {
>          "Polyphen2_HDIV_pred": "PolyPhen",
>          "ClinPred_pred": "ClinVar",
>          "REVEL_score": null
>       }
>    }
>    
> }
> ```

> Annotation with refSeq as a BED file

> ``` json
> {
>    "annotations": {
>       "tests/databases/annotations/current/hg19/refGene.bed": {
>          "INFO": null
>       }
>    }
> }
> ```

> Annotation with dbNSFP REVEL annotation (as a VCF file) within
> configured annotation databases folders (default:
> '~/howard/databases/annotations/current') and assembly (default:
> 'hg19')

> ``` json
> {
>    "annotations": {
>       "dbnsfp42a.REVEL.vcf.gz": {
>          "REVEL_score": null,
>          "REVEL_rankscore": null
>       }
>    }
>    
> }
> ```

> Annotation with all available databases in Parquet for ''current
> release

> ``` json
> {
>    "parquet": {
>       "annotations": {
>          "ALL": {
>             "formats": ["parquet"],
>             "releases": ["current"]
>          }
>       }
>    }
> }
> ```

## bcftools

Annotation process using BCFTools. Provide a list of database files and
annotation fields.

Examples:

> Annotation with multiple databases in multiple formats

> ``` json
> {
>    "parquet": {
>       "bcftools": {
>          "/path/to/database1.vcf.gz": {
>             "field1": null,
>             "field2": "field2_renamed"
>          },
>          "database2.bed.gz": {
>             "INFO": null
>          }
>       }
>    }
> }
> ```

### annotations

Specify the list of database files in formats VCF or BED. Files need to
be compressed and indexed.

This parameter enables users to select specific database fields and
optionally rename them (e.g. '"field": null' to keep field name,
'"field": "new_name"' to rename field). Use 'INFO' or 'ALL' keyword to
select all fields within the database INFO/Tags header (e.g. '"INFO":
null', '"ALL": null').

If a full path is not provided, the system will automatically detect
files within database folders (see Configuration doc) and assembly (see
Parameter option).

Examples:

> Annotation with dbSNP database with all fields

> ``` json
> {
>    "annotations": {
>       "tests/databases/annotations/current/hg19/avsnp150.vcf.gz": {
>          "INFO": null
>       }
>    }
> }
> ```

> Annotation with dbNSFP (only PolyPhen, ClinVar and REVEL score), and
> rename fields

> ``` json
> {
>    "annotations": {
>       "tests/databases/annotations/current/hg19/dbnsfp42a.vcf.gz": {
>          "Polyphen2_HDIV_pred": "PolyPhen",
>          "ClinPred_pred": "ClinVar",
>          "REVEL_score": null
>       }
>    }
> }
> ```

> Annotation with refSeq as a BED file

> ``` json
> {
>    "annotations": {
>       "tests/databases/annotations/current/hg19/refGene.bed": {
>          "INFO": null
>       }
>    }
> }
> ```

> Annotation with dbNSFP REVEL annotation (as a VCF file) within
> configured annotation databases folders (default:
> '~/howard/databases/annotations/current') and assembly (default:
> 'hg19')

> ``` json
> {
>    "annotations": {
>       "dbnsfp42a.REVEL.vcf.gz": {
>          "REVEL_score": null,
>          "REVEL_rankscore": null
>       }
>    }
> }
> ```

## annovar

Annotation process using Annovar tool. Provides a list of keywords to
select Annovar databases, and defines Annovar options (see [Annovar
documentation](https://annovar.openbioinformatics.org)).

Examples:

> Annotation with multiple Annovar databases, with fields selection, and
> Annovar options

> ``` json
> {
>    "annovar": {
>       "annotations": {
>          "annovar_keyword2": {
>             "field1": null,
>             "field2": "field2_renamed",
>          },
>          "annovar_keyword3": {
>             "INFO": null
>          }
>       }
>    }
> }
> ```

### annotations

List of keywords refering to Annovar databases (see [Annovar Databases
documentation](https://annovar.openbioinformatics.org/en/latest/user-guide/download/)),
with a list of selected fields for each of them (rename available)

Examples:

> Annotation with ClinVar (fields CLNSIG and CLNDN renamed) and Cosmic
> (all fields)

> ``` json
> {
>    "annotations": {
>       "clinvar_20221231": {
>          "CLNSIG": "ClinVar_class"
>          "CLNDN": "ClinVar_desease",
>       },
>       "cosmic70": {
>          "INFO": null
>       },
>    }
> }
> ```

### options

List of options available with Annovar tool (see Annovar documentation).
As example, these options allows to define splicing threshold or HGVS
annotation with refGene database

Examples:

> HGVS Annotation with refGene (add 'refGene' to 'annotations') and a
> splicing threshold as 3

> ``` json
> {
>    "options": {
>       "splicing_threshold": 3,
>       "argument": "'-hgvs'"
>       }
>    }
> }
> ```

## snpeff

Annotation process using snpEff tool and options (see [snpEff
documentation](https://pcingola.github.io/SnpEff/snpeff/commandline/)).

Examples:

> Annotation with snpEff databases, with options for HGVS annotation and
> additional tags.

> ``` json
> {
>    "snpeff": {
>       "options": " -hgvs -noShiftHgvs -spliceSiteSize 3 -lof -oicr "
>    }
> }
> ```

### options

String (as command line) of options available such as:

- filters on variants (regions filter, specific changes as intronic or
  downstream)

- annotation (e.g. HGVS, loss of function)

- database (e.g. only protein coding transcripts, splice sites size)

Examples:

> Annotation with snpEff databases, with options to generate HGVS
> annotation, specify to not shift variants according to HGVS notation,
> define splice sites size to 3, add loss of function (LOF), Nonsense
> mediated decay and OICR tags.

> ``` json
> {
>    "options": " -hgvs -noShiftHgvs -spliceSiteSize 3 -lof -oicr "
> }
> ```

### stats

HTML file for snpEff stats. Use keyword 'OUTPUT' to generate file
according to output file.

Examples:

> Annotation with snpEff databases, and generate a specific stats in
> HTML format.

> ``` json
> {
>    "stats": "/path/to/stats.html"
> }
> ```

> Annotation with snpEff databases, and generate stats in HTML format
> associated with output file.

> ``` json
> {
>    "stats": "OUTPUT.html"
> }
> ```

### csvStats

CSV file for snpEff stats. Use keyword 'OUTPUT' to generate file
according to output file.

Examples:

> Annotation with snpEff databases, and generate a specific stats in CSV
> format.

> ``` json
> {
>    "csvStats": "/path/to/stats.csv"
> }
> ```

> Annotation with snpEff databases, and generate stats in CSV format
> associated with output file.

> ``` json
> {
>    "csvStats": "OUTPUT.csv"
> }
> ```

## snpsift

Annotation process using snpSift. Provide a list of database files and
annotation fields.

Examples:

> Annotation with multiple databases in multiple formats

> ``` json
> {
>    "snpsift": {
>       "annotations": {
>          "/path/to/database1.vcf.gz": {
>             "field1": null,
>             "field2": null
>          },
>          "/path/to/database2.vcf.gz": {
>             "field1": null,
>             "field2": null
>          }
>       }
>    }
> }
> ```

### annotations

Specify the list of database files in formats VCF. Files need to be
compressed and indexed.

This parameter enables users to select specific database fields and
optionally rename them (e.g. '"field": null' to keep field name,
'"field": "new_name"' to rename field). Use 'INFO' or 'ALL' keyword to
select all fields within the database INFO/Tags header (e.g. '"INFO":
null', '"ALL": null').

If a full path is not provided, the system will automatically detect
files within database folders (see Configuration doc) and assembly (see
Parameter option).

Examples:

> Annotation with dbSNP database with all fields

> ``` json
> {
>    "annotations": {
>       "tests/databases/annotations/current/hg19/avsnp150.vcf.gz": {
>          "INFO": null
>       }
>    }
> }
> ```

> Annotation with dbNSFP (only PolyPhen, ClinVar and REVEL score)

> ``` json
> {
>    "annotations": {
>       "tests/databases/annotations/current/hg19/dbnsfp42a.vcf.gz": {
>          "Polyphen2_HDIV_pred": "PolyPhen",
>          "ClinPred_pred": "ClinVar",
>          "REVEL_score": null
>       }
>    }
> }
> ```

> Annotation with dbNSFP REVEL annotation (as a VCF file) within
> configured annotation databases folders (default:
> '~/howard/databases/annotations/current') and assembly (default:
> 'hg19')

> ``` json
> {
>    "annotations": {
>       "dbnsfp42a.REVEL.vcf.gz": {
>          "REVEL_score": null,
>          "REVEL_rankscore": null
>       }
>    }
> }
> ```

## bigwig

Annotation process using BigWig files. Provide a list of database files
in BigWig format ('.bw') and annotation fields.

Examples:

> Annotation with multiple databases in BigWig format

> ``` json
> {
>    "bigwig": {
>       "annotations": {
>          "/path/to/database1.bw": {
>             "field1": null,
>             "field2": null
>          },
>          "/path/to/database2.bw": {
>             "field1": null,
>             "field2": null
>          }
>       }
>    }
> }
> ```

### annotations

Specify the list of database files in BigWig format.

This parameter enables users to select specific database fields and
optionally rename them (e.g. '"field": null' to keep field name,
'"field": "new_name"' to rename field). Use 'INFO' or 'ALL' keyword to
select all fields within the database INFO/Tags header (e.g. '"INFO":
null', '"ALL": null').

If a full path is not provided, the system will automatically detect
files within database folders (see Configuration doc) and assembly (see
Parameter option).

A URL can be provided as a database file (experimental). In this case,
associated header file will be automatically generated with ua uniq
value as the name of the file (cleaned for avoid special characters, and
'.bw' extension).

Examples:

> Annotation with GERP database with all fields

> ``` json
> {
>    "annotations": {
>       "tests/databases/annotations/current/hg19/gerp.bw": {
>          "INFO": null
>       }
>    }
> }
> ```

> Annotation with GERP (only gerp score)

> ``` json
> {
>    "annotations": {
>       "tests/databases/annotations/current/hg19/gerp.bw": {
>          "gerp": "GERP_score"
>       }
>    }
> }
> ```

> Annotation with GERP from a distante database (experimental)

> ``` json
> {
>    "annotations": {
>       "https://hgdownload.soe.ucsc.edu/gbdb/hg19/bbi/All_hg19_RS.bw": {
>          "INFO": null
>       }
>    }
> }
> ```

## exomiser

Annotation process using Exomiser tool and options (see [Exomiser
website documentation](https://www.sanger.ac.uk/tool/exomiser/)).

Examples:

> Annotation with Exomiser, using database release '2109', transcripts
> source as UCSC and a list of HPO terms.

> ``` json
> {
>    "exomiser": {
>       "release": "2109"
>       "transcript_source": "refseq"
>       "hpo": ["HP:0001156", "HP:0001363", "HP:0011304", "HP:0010055"]
>    }
> }
> ```

### release

Release of Exomiser database. This option replace the 'release' variable
in 'application.properties' file (see 'exomiser_application_properties'
option). The release will be downloaded if it is not available locally.

Examples:

> Annotation with release '2109' of Exomiser database.

> ``` json
> {
>    "release": "2109"
> }
> ```

### transcript_source

Transcription source of Exomiser. This option replace the
'transcript_source' variable in 'application.properties' file (see
'exomiser_application_properties' option). The release will be
downloaded if it is not available locally.

Examples:

> Annotation with transcription source 'refseq' of Exomiser.

> ``` json
> {
>    "transcript_source": "refseq"
> }
> ```

### hpo

List of HPO for Exomiser. This option replace the 'hpo' variable in
'application.properties' file (see 'exomiser_application_properties'
option). The release will be downloaded if it is not available locally.

Examples:

> Annotation with a list of 4 HPO for Exomiser.

> ``` json
> {
>    "hpo": ["HP:0001156", "HP:0001363", "HP:0011304", "HP:0010055"]
> }
> ```

## splice

Annotation process using Splice tool and options. This annotation will
be proccessed only for variants that are not already annotated (i.e.
without annotation like 'SpliceAI\_\*' and 'SPiP\_\*')

Examples:

> Annotation with Splice, using database splice mode ('one'), spliceAI
> distance (500) and spliceAI mask (1).

> ``` json
> {
>    "splice": {
>       "split_mode": "one",
>       "spliceai_distance": 500,
>       "spliceai_mask": 1
>    }
> }
> ```

### split_mode

Split mode of Exomiser database (default 'one'):

- all: report all annotated transcript for one gene.

- one: keep only the transcript with the most pathogenic score (in case
  of identical score, take the first).

- list: keep transcript provided in transcript file, if no matching
  transcript in file 'one' mode is activated.

- mixed: 'one' mode, if identical score, list mode is activated.

Examples:

> Split mode to report all annotated transcript for one gene.

> ``` json
> {
>    "split_mode": "all"
> }
> ```

### spliceai_distance

Maximum distance between the variant and gained/lost splice site
(default: 500).

Examples:

> Maximum distance of '500' between variant and splice site.

> ``` json
> {
>    "spliceai_distance": 500
> }
> ```

### spliceai_mask

Mask scores representing annotated acceptor/donor gain and unannotated
acceptor/donor loss (default: 1).

Examples:

> Mask score of '1' for acceptor/donor gain fain and loss.

> ``` json
> {
>    "spliceai_mask": 1
> }
> ```

### transcript

Path to a list of transcripts of preference (default '').

Examples:

> Path to file of transcripts.

> ``` json
> {
>    "transcript": "tests/data/transcripts.tsv"
> }
> ```

### rm_snps

Do not consider SNV for the analysis, only Indels and MNV (default
'false').

Examples:

> Analysing only non SNV.

> ``` json
> {
>    "rm_snps": "true"
> }
> ```

### rm_annot

Remove existing annotation before analysing (default 'true').

Examples:

> Remove annotation before analysing.

> ``` json
> {
>    "rm_annot": "true"
> }
> ```

### whitespace

Remove spaces in INFO field, 'true' to remove (default 'true').

Examples:

> Remove spaces in INFO field.

> ``` json
> {
>    "whitespace": "true"
> }
> ```

## options

Options for annotations, such as annotation strategy (skip if exists,
update, append)

Examples:

> Annotation with Parquet databases, with update annotation strategy.

> ``` json
> {
>    "options": {
>       "update": true
>    }
> }
> ```

### annotations_update

Update option for annotation (only for Parquet annotation). If True,
annotation fields will be removed and re-annotated. These options will
be applied to all annotation databases.

Default: `False`

Examples:

> Apply update on all annotation fields for all databases.

> ``` json
> {
>    "update": true
> }
> ```

### annotations_append

Append option for annotation (only for Parquet annotation). If True,
annotation fields will be annotated only if not annotation exists for
the variant. These options will be applied to all annotation databases.

Default: `False`

Examples:

> Apply append on all annotation fields for all databases.

> ``` json
> {
>    "append": true
> }
> ```

# calculation

Calculation process operations that are defiend in a Calculation
Configuration JSON file. List available calculation operations with
possible options (see [Calculation JSON file](help.calculation.md)
help).

Examples:

> Calculation of operations 'operation1' and 'operation2' (with options)
> defined in 'calculation_config.json' file

> ``` json
> {
>    "calculation": {
>      "calculations": {
>        "operation1": null,
>        "operation2": {
>          "options": {
>            "option1": "value1",
>            "option2": "value2"
>          }
>        }
>      },
>      "calculation_config": "calculation_config.json"
>    }
> }
> ```

## calculations

List of operations to process with possible options (see [Calculation
JSON file](help.calculation.md) help).

Examples:

> Calculation with operations for generate variant_id and variant type,
> extract HGVS from snpEff annotation, select NOMEN from snpEff HGVS
> with a prioritized transcript (from prioritization transcript
> calculation) and list of transcripts of preference, with a specific
> NOMEN pattern

> ``` json
> {
>    "calculations": {
>      "variant_id": null,
>      "vartype": null,
>      "snpeff_hgvs": null,
>      "NOMEN": {
>        "options": {
>          "hgvs_field": "snpeff_hgvs",
>          "transcripts": "tests/data/transcripts.tsv",
>          "transcripts_table": "variants",
>          "transcripts_column": "PZTTranscript",
>          "transcripts_order", ["column", "file"],
>          "pattern": "GNOMEN:TNOMEN:ENOMEN:CNOMEN:RNOMEN:NNOMEN:PNOMEN"
>        }
>      }
>    }
> }
> ```

## calculation_config

Calculation configuration JSON file.

Type: `Path`

Default: `None`

Examples:

> Calculation configuration JSON file as an option

> ``` json
> {
>    "calculation_config": "calculation_config.json" 
> }
> ```

# prioritization

Prioritization process use a JSON configuration file that defines all
profiles that can be used. By default, all profiles will be calculated
from the JSON configuration file, and the first profile will be
considered as default. Proritization annotations (INFO/tags) will be
generated depending of a input list (default 'PZScore' and 'PZFlag'),
for all profiles (e.g. 'PZScore_GERMLINE' for 'GERMLINE' profile) and
for default profile (e.g. 'PZScore' for default). Prioritization score
mode is 'HOWARD' by default.

Examples:

> Prioritization with 'GENOME' and 'GERMLINE' (default) profiles, from a
> list of configured profiles, only 3 prioritization fields returned,
> and score calculated in 'VaRank' mode

> ``` json
> {
>    "prioritization": {
>      "prioritizations": "config/prioritization_profiles.json",
>      "profiles": ["GENOME", "GERMLINE"],
>      "default_profile": "GERMLINE",
>      "pzfields": ["PZScore", "PZFlag", "PZComment"],
>      "prioritization_score_mode": "VaRank"
>    }
> }
> ```

## prioritizations

Prioritization configuration profiles JSON file defining profiles to
calculate. All configured profiles will be calculated by default (see
'profiles' parameter). First profile will be considered as 'default' if
none are provided (see 'default_profile' parameter). Default score
calculation mode is 'HOWARD'. This option refers to the quick
prioritization command line parameter `--prioritizations`.

Type: `str`

Default: `None`

Examples:

> Prioritization configuration profile JSON file

> ``` json
> {
>    "prioritizations": "config/prioritization_profiles.json"
> }
> ```

## profiles

Prioritization profiles to consider, from the list of configured
profiles. If empty, all configured profiles will be calculated. First
profile will be considered as 'default' if none are provided (see
'default_profile' parameter). Prioritization annotations (INFO/tags)
will be generated for all these profiles (e.g. 'PZScore_GERMLINE' for
'GERMLINE' profile).

Type: `str`

Default: `None`

Examples:

> Prioritization with 'GERMLINE' profile only

> ``` json
> {
>    "profiles": ["GERMLINE"]
> }
> ```

> Prioritization with 'GENOME' and 'GERMLINE' profiles

> ``` json
> {
>    "profiles": ["GENOME", "GERMLINE"]
> }
> ```

## default_profile

Prioritization default profile from the list of processed profiles.
Prioritization annotations (INFO/tags) will be generated for this
default profile (e.g. 'PZScore', 'PZFlags').

Type: `str`

Default: `None`

Examples:

> Prioritization default profile 'GERMLINE'

> ``` json
> {
>    "default_profile": "GERMLINE"
> }
> ```

## pzfields

Prioritization annotations (INFO/tags) to generate. By default
'PZScore', 'PZFlags'.

Prioritization fields can be selected from:

- PZScore: calculated score from all passing filters, depending of the
  mode

- PZFlag: final flag ('PASS' or 'FILTERED'), with strategy that consider
  a variant is filtered as soon as at least one filter do not pass. By
  default, the variant is considered as 'PASS' (no filter pass)

- PZComment: concatenation of all passing filter comments

- PZTags: combinason of score, flags and comments in a tags format (e.g.
  'PZFlag#PASS\|PZScore#15\|PZComment#Described on ...')

- PZInfos: information about passing filter criteria

Type: `str`

Default: `PZScore,PZFlag`

Examples:

> Prioritization annotations (INFO/tags) list

> ``` json
> {
>    "pzfields": ["PZScore", "PZFlag", "PZComment"]
> }
> ```

## prioritization_score_mode

Prioritization score can be calculated following multiple mode. The
HOWARD mode will increment scores of all passing filters (default). The
VaRank mode will select the maximum score from all passing filters.

Type: `str`

Choices: `['HOWARD', 'VaRank']`

Default: `HOWARD`

Examples:

> Prioritization score calculation mode 'HOWARD'

> ``` json
> {
>    "prioritization_score_mode": "HOWARD"
> }
> ```

> Prioritization score calculation mode 'VaRank'

> ``` json
> {
>    "prioritization_score_mode": "VaRank"
> }
> ```

## pzprefix

Prioritization prefix for all annotations generated by prioritization.

Type: `str`

Default: `PZ`

Examples:

> Prioritization prefix by default ('PZ'):

> ``` json
> {
>    "pzprefix": "PZ"
> }
> ```

> Specific prioritization prefix:

> ``` json
> {
>    "pzprefix": "PrioritiZation_"
> }
> ```

> Specific prioritization prefix for transcript (see below):

> ``` json
> {
>    "pzprefix": "PZT"
> }
> ```

# stats

Statistics on loaded variants.

## stats_md

Stats Output file in MarkDown format.

Type: `Path`

Default: `None`

Examples:

> Export statistics in Markdown format

> ``` json
> {
>    "stats_md": "/tmp/stats.md" 
> }
> ```

## stats_json

Stats Output file in JSON format.

Type: `Path`

Default: `None`

Examples:

> Export statistics in JSON format

> ``` json
> {
>    "stats_json": "/tmp/stats.json" 
> }
> ```

# query

Query options tools. Mainly a SQL query, based on 'variants' table
corresponding on input file data, or a independant query. Print options
for 'query' tool allow limiting number of lines and choose printing
mode.

Type: `str`

Default: `None`

## query

Query in SQL format (e.g. 'SELECT \* FROM variants LIMIT 50').

Type: `str`

Default: `None`

Examples:

> Simple query to show all variants file

>     SELECT "#CHROM", POS, REF, ALT, INFO 
>     FROM variants

## query_limit

Limit of number of row for query (only for print result, not output).

Type: `int`

Default: `10`

## query_print_mode

Print mode of query result (only for print result, not output). Either
None (native), 'markdown', 'tabulate' or disabled.

Type: `str`

Choices: `[None, 'markdown', 'tabulate', 'disabled']`

Default: `None`

# export

Export options for output files, such as data order, include header in
output and hive partitioning.

## order_by

List of columns to sort the result-set in ascending or descending order.
Use SQL format, and keywords ASC (ascending) and DESC (descending). If a
column is not available, order will not be considered. Order is enable
only for compatible format (e.g. TSV, CSV, JSON). Examples: 'ACMG_score
DESC', 'PZFlag DESC, PZScore DESC'.

Type: `str`

Default: ``

Examples:

> Order by ACMG score in descending order

> ``` json
> {
>    "order_by": "ACMG_score DESC" 
> }
> ```

> Order by PZFlag and PZScore in descending order

> ``` json
> {
>    "order_by": PZFlag DESC, PZScore DESC" 
> }
> ```

## include_header

Include header (in VCF format) in output file. Only for compatible
formats (tab-delimiter format as TSV or BED).

Default: `False`

## parquet_partitions

Parquet partitioning using hive (available for any format). This option
is faster parallel writing, but memory consuming. Use 'None' (string)
for NO partition but split parquet files into a folder. Examples:
'#CHROM', '#CHROM,REF', 'None'.

Type: `str`

Default: `None`

# explode

Explode options for INFO/tags annotations within VCF files.

## explode_infos

Explode VCF INFO/Tag into 'variants' table columns.

Default: `False`

## explode_infos_prefix

Explode VCF INFO/Tag with a specific prefix.

Type: `str`

Default: ``

## explode_infos_fields

Explode VCF INFO/Tag specific fields/tags. Keyword `*` specify all
available fields, except those already specified. Pattern (regex) can be
used, such as `.*_score` for fields named with '\_score' at the end.
Examples:

- 'HGVS,SIFT,Clinvar' (list of fields)

- 'HGVS,\*,Clinvar' (list of fields with all other fields at the end)

- 'HGVS,.\*\_score,Clinvar' (list of 2 fields with all scores in the
  middle)

- 'HGVS,.\*\_score,\*' (1 field, scores, all other fields)

- 'HGVS,*,.*\_score' (1 field, all other fields, all scores)

Type: `str`

Default: `*`

# transcripts

Transcripts information to create transcript view. Useful to add
transcripts annotations in INFO field, to calculate transcripts specific
scores (see [Calculation JSON file](help.configuration.calculation.md)
help), to merge and map transcript IDs (e.g. from Ensembl to refSeq), or
prioritize transcripts (see [Priorization JSON
file](help.configuration.prioritization.md) help).

Type: `dict`

Default: `{}`

Examples:

> Trancripts information from snpEff and dbNSFP annotation

> ``` json
> {
>    "transcripts": {
>      "table": "transcripts",
>      "transcripts_info_field_json": "transcripts_json",
>      "transcripts_info_field_format": "transcripts_ann",
>      "transcripts_info_json": "transcripts_json",
>      "transcripts_info_format": "transcripts_format",
>      "transcript_id_remove_version": true,
>      "transcript_id_mapping_file": "transcripts.for_mapping.tsv",
>      "transcript_id_mapping_force": false,
>      "struct": {
>          "from_column_format": [
>              {
>                  "transcripts_column": "ANN",
>                  "transcripts_infos_column": "Feature_ID",
>                  "column_clean": true,
>                  "column_case": "lower"
>              }
>          ],
>          "from_columns_map": [
>              {
>                  "transcripts_column": "Ensembl_transcriptid",
>                  "transcripts_infos_columns": [
>                      "genename",
>                      "Ensembl_geneid",
>                      "LIST_S2_score",
>                      "LIST_S2_pred"
>                  ],
>                  "column_rename": {
>                      "LIST_S2_score": "LISTScore",
>                      "LIST_S2_pred": "LISTPred"
>                  },
>              },
>              {
>                  "transcripts_column": "Ensembl_transcriptid",
>                  "transcripts_infos_columns": [
>                      "genename",
>                      "VARITY_R_score",
>                      "Aloft_pred"
>                  ]
>              }
>          ]
>      },
>      "prioritization": {
>         "profiles": ["transcripts"],
>         "prioritization_config": "config/prioritization_transcripts_profiles.json",
>         "pzprefix": "PZT",
>         "pzfields": ["Score", "Flag", "LISTScore", "LISTPred"],
>         "prioritization_transcripts_order": {
>              "PZTFlag": "ASC",
>              "PZTScore": "DESC"
>         }
>         "prioritization_score_mode": "HOWARD",
>         "prioritization_transcripts": null,
>         "prioritization_transcripts_force": false,
>         "prioritization_transcripts_version_force": false
>      },
>      "export": {
>         "output": "/tmp/output.tsv.gz"
>      }
>    }
> }
> ```

## table

Transcripts table name to create.

Type: `str`

Default: `transcripts`

Examples:

> Name of transcript table:

> ``` json
> {
>    "table": "transcripts"
> }
> ```

## transcripts_info_field_json

Transcripts INFO field name to add in VCF INFO field in JSON format.

Type: `str`

Default: `None`

Examples:

> Transcripts INFO field name:

> ``` json
> {
>    "transcripts_info_field_json": "transcripts_json"
> }
> ```

## transcripts_info_field_format

Transcripts INFO field name to add in VCF INFO field in strutured
format.

Type: `str`

Default: `None`

Examples:

> Transcripts INFO field name:

> ``` json
> {
>    "transcripts_info_field_format": "transcripts_ann"
> }
> ```

## transcripts_info_json

Transcripts column name to add to transcripts table in JSON format.

Type: `str`

Default: `None`

Examples:

> Transcripts column name:

> ``` json
> {
>    "transcripts_info_json": "transcripts_json"
> }
> ```

## transcripts_info_format

Transcripts column name to add to transcripts table in structured
format.

Type: `str`

Default: `None`

Examples:

> Transcripts column name:

> ``` json
> {
>    "transcripts_info_format": "transcripts_format"
> }
> ```

## transcript_id_remove_version

When merging and mapping transcript IDs, remove possible version of
transcript (e.g 'NM_123456.2' to 'NM_123456').

Type: `Boolean`

Default: `False`

Examples:

> Remove transcript version when merging and mapping:

> ``` json
> {
>    "transcript_id_remove_version": true
> }
> ```

## transcript_id_mapping_file

When merging and mapping transcript IDs, indicate a transcript mapping
file that provides mapping between transcripts IDs (useful to map refSeq
and Ensembl transcript IDs).

Type: `Path`

Default: `None`

Examples:

> Transcript IDs mapping file:

> ``` json
> {
>    "transcript_id_mapping_file": "My_transcripts_mapping_file.tsv.gz"
> }
> ```

## Example of transcript ID mapping file

Transcript IDs mapping file is a tab-delimited file with 2 columns:

- first column corresponds to the reference transcript ID to use

- second column correspond to an alias of the reference transcript

Second column can be empty (no alias is provided). Transcript IDs can
include version or not (see [transcript_id_remove_version
section](#transcript_id_remove_version))

Examples:

> Example of transcripts mapping file:

> ``` ts
> NM_001005484    ENST00000641515.1
> NM_005228.5     ENST00000275493.7
> NM_001346897    ENSG00000146648.21
> NM_001346941    
> NM_005228       
> ```

## transcript_id_mapping_force

When merging and mapping transcript IDs, allows to filter transcript IDs
only if they are present in first column of the transcript mapping file
(see [transcript_id_mapping_file section](#transcript_id_mapping_file)).

Type: `Boolean`

Default: `False`

Examples:

> Filter transcripts IDs only if present in mapping file:

> ``` json
> {
>    "transcript_id_mapping_force": true
> }
> ```

## struct

Structure of transcripts information, corresponding to columns dedicated
to transcripts, such as:

- 'from_column_format': a uniq annotation field with a specific format,
  like snpEff annotation,

- 'from_columns_map': a list of annotation fields corresponding to
  transcripts in another specific field, like dbNSFP annotation.

Some parameters are commons between these structure (e.g.
'column_rename', 'column_clean' and 'column_case').

Type: `dict`

Default: `{}`

### from_column_format

Structure of transcripts information from a uniq annotation field with a
specific format (such as snpEff annotation):

- 'transcripts_column' correspond to INFO field with annotations

- 'transcripts_infos_column' correspond to transcription ID annotations
  field within INFO field

Column can be renamed, cleaned and/or case changed (see below).

Type: `dict`

Default: `{}`

Examples:

> Structure from snpEff annotation (columns names must be clean or
> changed for standard snpEff annotations):

> ``` json
> {
>    "from_column_format": [
>      {
>        "transcripts_column": "ANN",
>        "transcripts_infos_column": "Feature_ID",
>        "column_rename": null,
>        "column_clean": true,
>        "column_case": null
>      }
>    ]
> }
> ```

### from_columns_map

list of annotation fields corresponding to transcripts in another
specific field (such as dbNSFP annotation):

- 'transcripts_column' correspond to INFO field with transcription ID

- 'transcripts_infos_columns' correspond to a list of INFO fields with
  transcript information

Column can be renamed, cleaned and/or case changed (see below).

Type: `dict`

Default: `{}`

Examples:

> Structure from dbNSFP annotations (with 2 columns renamed):

> ``` json
> {
>    "from_columns_map": [
>      {
>        "transcripts_column": "Ensembl_transcriptid",
>        "transcripts_infos_columns": [
>          "genename",
>          "Ensembl_geneid",
>          "LIST_S2_score",
>          "LIST_S2_pred"
>        ],
>        "column_rename": {
>          "LIST_S2_score": "LISTScore",
>          "LIST_S2_pred": "LISTPred"
>        },
>        "column_clean": false,
>        "column_case": null
>      }
>    ]
> }
> ```

### commons parameters

Some parameters are commons between these structure:

- 'column_rename': dict defining mapping of column name changes

- 'column_clean': if true, clean column name to remove not alphanum
  characters not allowed in VCF (e.g. space, dash)

- 'column_case': rename column into lowercase ('lower') or uppercase
  ('upper')

Combining 'column_clean' and 'column_case' ensure well formed VCF field
name and merging identical columns (e.g. same field 'Gene_Name', 'Gene
name' and 'genename' from multiple sources). However, controling column
names through 'column_rename' is much more efficient.

Examples:

> Commons parameters by default:

> ``` json
> {
>   "column_rename": null,
>   "column_clean": false,
>   "column_case": null
> }
> ```

> Parameters to change 2 column names:

> ``` json
> {
>   "column_rename": {
>     "LIST_S2_score": "LISTScore",
>     "LIST_S2_pred": "LISTPred"
>   },
>   "column_clean": false,
>   "column_case": null
> }
> ```

> Parameters to ensure well-named column from format annotations field
> (such as snpEff):

> ``` json
> {
>   "column_rename": null,
>   "column_clean": true,
>   "column_case": null
> }
> ```

## prioritization

Prioritization parameters for transcripts (see [Prioritization
section](#prioritization) and [Priorization JSON
file](help.configuration.prioritization.md) help), defining
prioritization criteria with a configuration file of all available
profiles and which profiles to use, which prioritization method to use
(e.g. 'HOWARD', 'VaRank').

Prioritized transcripts fields can be defined to provide VCF fields
specific to the choosen transcript, using a specific prefix (e.g.
'PZTScore', 'PZTFlag'). The selected 'best' transcript ID is always
provided (e.g. 'PZTTranscript'). Extra annotation fields can also be
defined (e.g. 'LISTScore', 'LISTPred'), as well as prioritizations
informations from multiple profiles (e.g. 'PZTScore_myprofile',
'PZTScore_myotherprofile' in 'pzfields' section with 2 profiles
'myprofile' and 'myotherprofile' in 'profiles' section).

In order to choose the 'best' transcript, parameteres can define order
of transcripts (by annotation columns or a preference transcripts file),
by dealing with transcript versions.

Type: `dict`

Default: `{}`

Examples:

> Prioritization of transcripts in 'HOWARD' mode with 'transcripts'
> profiles available in a configuration JSON file, with 'PZT' as prefix:

> ``` json
> {
>    "prioritization": {
>       "profiles": ["transcripts"],
>       "default_profile": "transcripts",
>       "prioritization_config": "config/prioritization_transcripts_profiles.json",
>       "prioritization_score_mode": "HOWARD",
>       "pzprefix": "PZT",
>       "pzfields": ["Score", "Flag", "LISTScore", "LISTPred"],
>       "prioritization_transcripts_order": {
>          "PZTFlag": "ASC",
>          "PZTScore": "DESC"
>       },
>       "prioritization_transcripts": null,
>       "prioritization_transcripts_force": false,
>       "prioritization_transcripts_version_force": false
>    }
> }
> ```

### profiles

See [Prioritization section](#prioritization)

Type: `str`

Default: `None`

### default_profile

See [Prioritization section](#prioritization)

Type: `str`

Default: `None`

### prioritization_config

See [Prioritization section](#prioritization)

Type: `Path`

Default: `None`

Examples:

> Prioritization configuration JSON file as an option

> ``` json
> {
>    "prioritization_config": "prioritization_config.json" 
> }
> ```

### prioritization_score_mode

See [Prioritization section](#prioritization)

Type: `str`

Choices: `['HOWARD', 'VaRank']`

Default: `HOWARD`

### pzprefix

See [Prioritization section](#prioritization)

### pzfields

See [Prioritization section](#prioritization)

Type: `str`

Default: `PZScore,PZFlag`

### prioritization_transcripts_order

Defines the order of transcripts to determine which one is chosen (by
default PZTFlag and PZTScore). All available annotation can be used
(e.g. scores, length of transcript, predictions...). The first
transcript will be choosen in case of equal order.

Type: `dict`

Default: `{}`

Examples:

> Default order of transcript using Flag and Score:

> ``` json
> {
>    "prioritization_transcripts_order": {
>       "PZTFlag": "ASC",
>       "PZTScore": "DESC"
>    }
> }
> ```

> Order of transcript using Flag, Score and additional specific
> spliceAI_score:

> ``` json
> {
>    "prioritization_transcripts_order": {
>       "PZTFlag": "ASC",
>       "PZTScore": "DESC",
>       "spliceAI_score": "DESC"
>    }
> }
> ```

### prioritization_transcripts

Defines a file with an ordered list of transcripts of preference. The
first transcript in this list will be chosen, if no order is define (see
above), or if this list is not forced (see below).

Type: `Path`

Default: `None`

Examples:

> File with transcripts of preference:

> ``` json
> {
>    "prioritization_transcripts": "transcripts.tsv"
> }
> ```

### prioritization_transcripts_force

Force to use the list of transcripts of preference define in the
provided file (see above).

Type: `Boolean`

Default: `False`

Examples:

> Force using transcripts of preference in file:

> ``` json
> {
>    "prioritization_transcripts_force": true
> }
> ```

### prioritization_transcripts_version_force

By default, versions of transcripts are not considered when comparison
is needed (e.g. transcript ID and list of transcript of preference). If
true, all transcript ID will be considered with their version.

Type: `Boolean`

Default: `False`

Examples:

> Force using transcripts version:

> ``` json
> {
>    "prioritization_transcripts_version_force": true
> }
> ```

## export

Options to export transcripts view/table into a file ('output'
parameter). All HOWARD format are available (e.g. VCF, Parquet, TSV).
For VCF format, all columns will be concatenate into INFO column,
otherwise each column will be exported.

Type: `Dict`

Default: `{}`

Examples:

> Export as compressed TSV:

> ``` json
> {
>    "export": {
>      "output": "/tmp/output.tsv.gz"
>    }
> }
> ```

> Export as VCF:

> ``` json
> {
>    "export": {
>      "output": "/tmp/output.vcf"
>    }
> }
> ```

> Export as Parquet:

> ``` json
> {
>    "export": {
>      "output": "/tmp/output.parquet"
>    }
> }
> ```

# threads

Specify the number of threads to use for processing HOWARD. It
determines the level of parallelism, either on python scripts, duckdb
engine and external tools. It and can help speed up the process/tool.
Use -1 to use all available CPU/cores. Either non valid value is 1
CPU/core.

Type: `int`

Default: `-1`

Examples:

> Automatically detect all available CPU/cores

> ``` json
> {
>    "threads": -1
> }
> ```

> Define 8 CPU/cores

> ``` json
> {
>    "threads": 8
> }
> ```

# samples

Samples parameters to defined a list of samples or use options to check
samples. Only for export in VCF format. By default, if no samples are
listed, all existing samples are checked if they contain well-formed
genotype annotations (based on 'FORMAT' VCF column).

Type: `dict`

Examples:

> Export only a list of samples:

> ``` json
> {
>    "samples": {
>       "list": ["sample1", "sample2"]
>    }
> }
> ```

> Do not check existing samples (all VCF columns after FORMAT column):

> ``` json
> {
>    "samples": {
>       "check": false
>    }
> }
> ```

> Default configuration, with all samples are considered (null) and
> checked (true):

> ``` json
> {
>    "samples": {
>       "list": null,
>       "check": true
>    }
> }
> ```

## list

List of columns that correspond to samples (with well formed genotype,
based on 'FORMAT' VCF column). Only for export in VCF format. Only these
samples are exported in VCF format file.

Type: `dict`

Examples:

> Export only a list of samples:

> ``` json
> {
>    "list": ["sample1", "sample2"]
> }
> ```

## check

Check if samples (either provided in 'list' parameters, or all existing
column after 'FORMAT' column) according to 'FORMAT' VCF column. Only for
export in VCF format. By default, samples are checked (beware of format
if check is disabled) and removed if they are not well-formed.

Type: `dict`

Examples:

> Do not check existing samples:

> ``` json
> {
>    "check": false
> }
> ```

# databases

[HOWARD Parameters Databases JSON](help.parameters.databases.md)
describes configuration JSON file for databases download and convert.
