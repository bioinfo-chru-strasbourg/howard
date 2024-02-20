# HOWARD Parameters

HOWARD Parameters JSON file defined parameters to process annotations, prioritization, calculations, convertions and queries.

## Table of contents

- [HOWARD Parameters](#howard-parameters)
   - [annotations](#annotations)
      - [parquet](#annotationsparquet)
         - [annotations](#annotationsparquetannotations)
      - [annovar](#annotationsannovar)
         - [annotations](#annotationsannovarannotations)
         - [options](#annotationsannovaroptions)


## annotations

Annotation process using HOWARD algorithms or external tools.

For HOWARD Parquet algorithm, specify the list of database files available (formats such as Parquet, VCF, TSV, duckDB, JSON). This parameter enables users to select specific database fields and optionally rename them. Ue 'INFO' keyword to select all fields within the database. If a full path is not provided, the system will automatically detect files within database folders (see Configuration doc) and assembly (see Parameter option).

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

Annotation process using HOWARD parquet algorithm. Provide a list of database files and annotation fields

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
   "tests/databases/annotations/hg19/avsnp150.parquet": {
      "INFO": null
   }
}

# Annotation with dbNSFP (only PolyPhen, ClinVar and REVEL score), and rename fields
"annotations": {
   "tests/databases/annotations/hg19/dbnsfp42a.parquet": {
      "Polyphen2_HDIV_pred": "PolyPhen",
      "ClinPred_pred": "ClinVar",
      "REVEL_score": null
   }
}

# Annotation with refSeq as a BED file
"annotations": {
   "tests/databases/annotations/hg19/refGene.bed": {
      "INFO": null
   }
}

# Annotation with dbNSFP REVEL annotation (as a VCF file) within configured annotation databases folders (default: '/databases/annotations/current') and assembly (default: 'hg19')
"annotations": {
   "dbnsfp42a.REVEL.vcf.gz": {
      "REVEL_score": null,
      "REVEL_rankscore": null
   }
}

# Annotation with dbNSFP REVEL annotation (as a VCF file) within configured annotation databases folders (default: '/databases/annotations/current') and assembly (default: 'hg19')
"annotations": {
   "dbnsfp42a.REVEL.vcf.gz": {
      "REVEL_score": null,
      "REVEL_rankscore": null
   }
}
```

### annotations::annovar

Annotation process using Annovar tool. Provides a list of keywords to select Annovar databases, and defines Annovar options (see Annovar documentation).

Examples: 
```
# Annotation with multiple Annovar databases, with fields selection, and Annovar options
"annovar": {
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
   "options": {
      "splicing_threshold": 3,
      ...,
   }
}
```

#### annotations::annovar::annotations

List of keywords refering to Annovar databases (see Annovar documentation), with a list of selected fields for each of them (rename available)

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

List of options available with Annovar tool (see Annovar documentation). As example, these options allows to define splicing threashold or HGVS annotation with refGene database

Examples: 
```
# HGVS Annotation with refGene (add 'refGene' to 'annotations') and a splicing threshold as 3
"options": {
   "splicing_threshold": 3
   "argument": "'-hgvs'"
   }
}
```

