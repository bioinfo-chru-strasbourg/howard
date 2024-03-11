# HOWARD Parameters

HOWARD Parameters JSON file defined parameters to process annotations, prioritization, calculations, convertions and queries.

## Table of contents

- [HOWARD Parameters](#howard-parameters)
   - [annotations](#annotations)
      - [parquet](#annotationsparquet)
         - [annotations](#annotationsparquetannotations)
      - [bcftools](#annotationsbcftools)
         - [annotations](#annotationsbcftoolsannotations)
      - [annovar](#annotationsannovar)
         - [annotations](#annotationsannovarannotations)
         - [options](#annotationsannovaroptions)
      - [snpeff](#annotationssnpeff)
         - [options](#annotationssnpeffoptions)
         - [stats](#annotationssnpeffstats)
         - [csvStats](#annotationssnpeffcsvstats)
      - [exomiser](#annotationsexomiser)
         - [release](#annotationsexomiserrelease)
      - [options](#annotationsoptions)
         - [update](#annotationsoptionsupdate)
         - [append](#annotationsoptionsappend)


## annotations

Annotation process using HOWARD algorithms or external tools.

For HOWARD Parquet algorithm, specify the list of database files available (formats such as Parquet, VCF, TSV, duckDB, JSON). This parameter enables users to select specific database fields and optionally rename them. Use 'INFO' keyword to select all fields within the database. If a full path is not provided, the system will automatically detect files within database folders (see Configuration doc) and assembly (see Parameter option).

For external tools, such as Annovar, snpEff and Exomiser, specify parameters such as annotation keywords (Annovar) and options (depending on the tool).

Examples: 
```
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
         "/path/to/database1.txt": {
            ...
         },
         ...
      }
   }
}
```

### annotations::parquet

Annotation process using HOWARD parquet algorithm. Provide a list of database files and annotation fields.

Examples: 
```
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

#### annotations::parquet::annotations

Specify the list of database files available in multiple formats such as Parquet, VCF, BED, TSV, duckDB, JSON.

Examples: 
```
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

### annotations::bcftools

Annotation process using BCFTools. Provide a list of database files and annotation fields.

Examples: 
```
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

#### annotations::bcftools::annotations

Specify the list of database files in formats VCF or BED. Files need to be compressed and indexed.

Examples: 
```
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

### annotations::annovar

Annotation process using Annovar tool. Provides a list of keywords to select Annovar databases, and defines Annovar options (see [Annovar documentation](https://annovar.openbioinformatics.org)).

Examples: 
```
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

#### annotations::annovar::annotations

List of keywords refering to Annovar databases (see [Annovar Databases documentation](https://annovar.openbioinformatics.org/en/latest/user-guide/download/)), with a list of selected fields for each of them (rename available)

Examples: 
```
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

#### annotations::annovar::options

List of options available with Annovar tool (see Annovar documentation). As example, these options allows to define splicing threshold or HGVS annotation with refGene database

Examples: 
```
# HGVS Annotation with refGene (add 'refGene' to 'annotations') and a splicing threshold as 3
"options": {
   "splicing_threshold": 3
   "argument": "'-hgvs'"
   }
}
```

### annotations::snpeff

Annotation process using snpEff tool and options (see [snpEff documentation](https://pcingola.github.io/SnpEff/snpeff/commandline/)).

Examples: 
```
# Annotation with snpEff databases, with options for HGVS annotation and additional tags.
"snpeff": {
   "options": {
      " -hgvs -noShiftHgvs -spliceSiteSize 3 -lof -oicr "}
   }
}
```

#### annotations::snpeff::options

String (as command line) of options available such as:

 - filters on variants (regions filter, specific changes as intronic or downstream)

 - annotation (e.g. HGVS, loss of function) 

 - database (e.g. only protein coding transcripts, splice sites size)

Examples: 
```
# Annotation with snpEff databases, with options to generate HGVS annotation, specify to not shift variants according to HGVS notation, define splice sites size to 3, add loss of function (LOF), Nonsense mediated decay and OICR tags.
"options": " -hgvs -noShiftHgvs -spliceSiteSize 3 -lof -oicr "
```

#### annotations::snpeff::stats

HTML file for snpEff stats. Use keyword 'OUTPUT' to generate file according to output file.

Examples: 
```
# Annotation with snpEff databases, and generate a specific stats in HTML format.
"stats": "/path/to/stats.html"
# Annotation with snpEff databases, and generate stats in HTML format associated with output file.
"stats": "OUTPUT.html"
```

#### annotations::snpeff::csvStats

CSV file for snpEff stats. Use keyword 'OUTPUT' to generate file according to output file.

Examples: 
```
# Annotation with snpEff databases, and generate a specific stats in CSV format.
"csvStats": "/path/to/stats.csv"
# Annotation with snpEff databases, and generate stats in CSV format associated with output file.
"csvStats": "OUTPUT.csv"
```

### annotations::exomiser

Annotation process using Exomiser tool and options (see [Exomiser website documentation](https://www.sanger.ac.uk/tool/exomiser/)).

Examples: 
```
# Annotation with Exomiser, using database relse '2109', transcripts source as UCSC and a list of HPO terms.
"exomiser": {
   "release": "2109"
   "transcript_source": "refseq"
   "hpo": ['HP:0001156', 'HP:0001363', 'HP:0011304', 'HP:0010055']
}
```

#### annotations::exomiser::release

Release of Exomiser database. This option replace the release variable in 'application.properties' file (see 'exomiser_application_properties' option). The release will be downloaded if it is not available locally. 

Examples: 
```
# Annotation with release '2109' of Exomiser database.
"release": "2109"
```

### annotations::options

Options for annotations, such as annotation strategy (skip if exists, update, append)

Examples: 
```
# Annotation with Parquet databases, with update annotation strategy.
"options": {
   "update": True
}
```

#### annotations::options::update

Update option for annotation (only for Parquet annotation). If True, annotation fields will be removed and re-annotated. These options will be applied to all annotation databases.

Examples: 
```
# Apply update on all annotation fields for all databases.
"update": True
```

#### annotations::options::append

Append option for annotation (only for Parquet annotation). If True, annotation fields will be annotated only if not annotation exists for the variant. These options will be applied to all annotation databases.

Examples: 
```
# Apply append on all annotation fields for all databases.
"append": True
```

