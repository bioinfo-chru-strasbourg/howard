# HOWARD Help

<!--TOC-->

- [HOWARD Help](#howard-help)
  - [QUERY tool](#query-tool)
    - [Main options](#main-options)
    - [Explode infos](#explode-infos)
    - [Query](#query)
    - [Output](#output)
  - [STATS tool](#stats-tool)
    - [Main options](#main-options-1)
  - [CONVERT tool](#convert-tool)
    - [Main options](#main-options-2)
  - [ANNOTATION tool](#annotation-tool)
    - [Main options](#main-options-3)
  - [CALCULATION tool](#calculation-tool)
    - [Main options](#main-options-4)
    - [NOMEN calculation](#nomen-calculation)
    - [TRIO calculation](#trio-calculation)
  - [PRIORITIZATION tool](#prioritization-tool)
    - [Main options](#main-options-5)
    - [Prioritization](#prioritization)
  - [PROCESS tool](#process-tool)
    - [Main options](#main-options-6)
    - [Quick Processes](#quick-processes)
  - [HGVS tool](#hgvs-tool)
    - [Main options](#main-options-7)
    - [HGVS](#hgvs)
    - [Databases](#databases)
  - [DATABASES tool](#databases-tool)
    - [Main options](#main-options-8)
    - [Genomes](#genomes)
    - [snpEff](#snpeff)
    - [Annovar](#annovar)
    - [refSeq](#refseq)
    - [dbNSFP](#dbnsfp)
    - [AlphaMissense](#alphamissense)
    - [Exomiser](#exomiser)
    - [dbSNP](#dbsnp)
    - [HGMD](#hgmd)
    - [Parameters file](#parameters-file)
  - [FROM_ANNOVAR tool](#from_annovar-tool)
    - [Main options](#main-options-9)
    - [Annovar](#annovar-1)
    - [Parquet](#parquet)
    - [Modes](#modes)
  - [GUI tool](#gui-tool)
  - [HELP tool](#help-tool)
    - [Main options](#main-options-10)
  - [Shared arguments](#shared-arguments)

<!--TOC-->

HOWARD:1.0.0

Highly Open and Valuable tool for Variant Annotation & Ranking toward genetic Discovery

HOWARD annotates and prioritizes genetic variations, calculates and normalizes annotations, convert on multiple formats, query variations and generates statisticsUsage examples:

> howard process --input=tests/data/example.vcf.gz --output=/tmp/example.annotated.vcf.gz --param=config/param.json 

> howard annotation --input=tests/data/example.vcf.gz --output=/tmp/example.howard.vcf.gz --annotations='tests/databases/annotations/current/hg19/dbnsfp42a.parquet,tests/databases/annotations/current/hg19/gnomad211_genome.parquet' 

> howard calculation --input=tests/data/example.full.vcf --output=/tmp/example.calculation.tsv --calculations='vartype' 

> howard prioritization --input=tests/data/example.vcf.gz --output=/tmp/example.prioritized.vcf.gz --prioritizations=config/prioritization_profiles.json --profiles='default,GERMLINE' 

> howard query --input=tests/data/example.vcf.gz --explode_infos --query='SELECT "#CHROM", POS, REF, ALT, "DP", "CLNSIG", sample2, sample3 FROM variants WHERE "DP" >= 50 OR "CLNSIG" NOT NULL ORDER BY "CLNSIG" DESC, "DP" DESC' 

> howard stats --input=tests/data/example.vcf.gz 

> howard convert --input=tests/data/example.vcf.gz --output=/tmp/example.tsv --explode_infos && cat /tmp/example.tsv 



## QUERY tool
Query genetic variations in SQL format. Data can be loaded into 'variants' table from various formats (e.g. VCF, TSV, Parquet...). Using --explode_infos allow query on INFO/tag annotations. SQL query can also use external data within the request, such as a Parquet file(s).  

Usage examples:

> howard query --input=tests/data/example.vcf.gz --query="SELECT * FROM variants WHERE REF = 'A' AND POS < 100000" 

> howard query --input=tests/data/example.vcf.gz --explode_infos --query='SELECT "#CHROM", POS, REF, ALT, DP, CLNSIG, sample2, sample3 FROM variants WHERE DP >= 50 OR CLNSIG NOT NULL ORDER BY DP DESC' 

> howard query --query="SELECT \"#CHROM\", POS, REF, ALT, \"INFO/Interpro_domain\" FROM 'tests/databases/annotations/current/hg19/dbnsfp42a.parquet' WHERE \"INFO/Interpro_domain\" NOT NULL ORDER BY \"INFO/SiPhy_29way_logOdds_rankscore\" DESC LIMIT 10" 

> howard query --explode_infos --explode_infos_prefix='INFO/' --query="SELECT \"#CHROM\", POS, REF, ALT, STRING_AGG(INFO, ';') AS INFO FROM 'tests/databases/annotations/current/hg19/*.parquet' GROUP BY \"#CHROM\", POS, REF, ALT" --output=/tmp/full_annotation.tsv  && head -n2 /tmp/full_annotation.tsv 

> howard query --input=tests/data/example.vcf.gz --param=config/param.json

### Main options
```
--input=<input>

Input file path.
Format file must be either VCF, Parquet, TSV, CSV, PSV or duckDB.
Files can be compressesd (e.g. vcf.gz, tsv.gz).

```

```
--output=<output>

Output file path.
Format file must be either VCF, Parquet, TSV, CSV, PSV or duckDB.
Files can be compressesd (e.g. vcf.gz, tsv.gz).

```

```
--query=<query>

Query in SQL format
(e.g. 'SELECT * FROM variants LIMIT 50').

```

```
--param=<param> ({})

Parameters JSON file or JSON string.

```

### Explode infos
```
--explode_infos

Explode VCF INFO/Tag into 'variants' table columns.

```

```
--explode_infos_prefix=<explode infos prefix>

Explode VCF INFO/Tag with a specific prefix.

```

```
--explode_infos_fields=<explode infos list> (*)

Explode VCF INFO/Tag specific fields/tags.
Keyword `*` specify all available fields, except those already specified.
Pattern (regex) can be used, such as `.*_score` for fields named with '_score' at the end.
Examples:
- 'HGVS,SIFT,Clinvar' (list of fields)
- 'HGVS,*,Clinvar' (list of fields with all other fields at the end)
- 'HGVS,.*_score,Clinvar' (list of 2 fields with all scores in the middle)
- 'HGVS,.*_score,*' (1 field, scores, all other fields)
- 'HGVS,*,.*_score' (1 field, all other fields, all scores)

```

### Query
```
--query_limit=<query limit> (10)

Limit of number of row for query (only for print result, not output).

```

```
--query_print_mode=<print mode> [None, 'markdown', 'tabulate']

Print mode of query result (only for print result, not output).
Either None (native), 'markdown' or 'tabulate'.

```

### Output
```
--include_header

Include header (in VCF format) in output file.
Only for compatible formats (tab-delimiter format as TSV or BED).

```



## STATS tool
Statistics on genetic variations, such as: number of variants, number of samples, statistics by chromosome, genotypes by samples...

Usage examples:

> howard stats --input=tests/data/example.vcf.gz 

> howard stats --input=tests/data/example.vcf.gz --stats_md=/tmp/stats.md 

> howard stats --input=tests/data/example.vcf.gz --param=config/param.json 

### Main options
```
--input=<input> | required

Input file path.
Format file must be either VCF, Parquet, TSV, CSV, PSV or duckDB.
Files can be compressesd (e.g. vcf.gz, tsv.gz).

```

```
--stats_md=<stats markdown>

Stats Output file in MarkDown format.

```

```
--stats_json=<stats json>

Stats Output file in JSON format.

```

```
--param=<param> ({})

Parameters JSON file or JSON string.

```



## CONVERT tool
Convert genetic variations file to another format. Multiple format are available, such as usual and official VCF and BCF format, but also other formats such as TSV, CSV, PSV and Parquet/duckDB. These formats need a header '.hdr' file to take advantage of the power of howard (especially through INFO/tag definition), and using howard convert tool automatically generate header file fo futher use. 

Usage examples:

> howard convert --input=tests/data/example.vcf.gz --output=/tmp/example.tsv 

> howard convert --input=tests/data/example.vcf.gz --output=/tmp/example.parquet 

> howard convert --input=tests/data/example.vcf.gz --output=/tmp/example.tsv --explode_infos --explode_infos_fields='CLNSIG,SIFT,DP' --order_by='CLNSIG DESC, DP DESC' 

> howard convert --input=tests/data/example.vcf.gz --output=/tmp/example.tsv --explode_infos --explode_infos_prefix='INFO/' --explode_infos_fields='CLNSIG,SIFT,DP,*' --order_by='"INFO/CLNSIG" DESC, "INFO/DP" DESC' --include_header 

> howard convert --input=tests/data/example.vcf.gz --output=/tmp/example.tsv --param=config/param.json 

### Main options
```
--input=<input> | required

Input file path.
Format file must be either VCF, Parquet, TSV, CSV, PSV or duckDB.
Files can be compressesd (e.g. vcf.gz, tsv.gz).

```

```
--output=<output> | required

Output file path.
Format file must be either VCF, Parquet, TSV, CSV, PSV or duckDB.
Files can be compressesd (e.g. vcf.gz, tsv.gz).

```

```
--explode_infos

Explode VCF INFO/Tag into 'variants' table columns.

```

```
--explode_infos_prefix=<explode infos prefix>

Explode VCF INFO/Tag with a specific prefix.

```

```
--explode_infos_fields=<explode infos list> (*)

Explode VCF INFO/Tag specific fields/tags.
Keyword `*` specify all available fields, except those already specified.
Pattern (regex) can be used, such as `.*_score` for fields named with '_score' at the end.
Examples:
- 'HGVS,SIFT,Clinvar' (list of fields)
- 'HGVS,*,Clinvar' (list of fields with all other fields at the end)
- 'HGVS,.*_score,Clinvar' (list of 2 fields with all scores in the middle)
- 'HGVS,.*_score,*' (1 field, scores, all other fields)
- 'HGVS,*,.*_score' (1 field, all other fields, all scores)

```

```
--order_by=<order by>

List of columns to sort the result-set in ascending or descending order.
Use SQL format, and keywords ASC (ascending) and DESC (descending).
If a column is not available, order will not be considered.
Order is enable only for compatible format (e.g. TSV, CSV, JSON).
Examples: 'ACMG_score DESC', 'PZFlag DESC, PZScore DESC'.

```

```
--include_header

Include header (in VCF format) in output file.
Only for compatible formats (tab-delimiter format as TSV or BED).

```

```
--parquet_partitions=<parquet partitions>

Parquet partitioning using hive (available for any format).
This option is faster parallel writing, but memory consuming.
Use 'None' (string) for NO partition but split parquet files into a folder.
Examples: '#CHROM', '#CHROM,REF', 'None'.

```

```
--param=<param> ({})

Parameters JSON file or JSON string.

```



## ANNOTATION tool
Annotation is mainly based on a build-in Parquet annotation method, and tools such as BCFTOOLS, Annovar and snpEff. It uses available databases (see Annovar and snpEff) and homemade databases. Format of databases are: parquet, duckdb, vcf, bed, Annovar and snpEff (Annovar and snpEff databases are automatically downloaded, see howard databases tool). 

Usage examples:

> howard annotation --input=tests/data/example.vcf.gz --output=/tmp/example.howard.vcf.gz --annotations='tests/databases/annotations/current/hg19/avsnp150.parquet,tests/databases/annotations/current/hg19/dbnsfp42a.parquet,tests/databases/annotations/current/hg19/gnomad211_genome.parquet' 

> howard annotation --input=tests/data/example.vcf.gz --output=/tmp/example.howard.tsv --assembly=hg19 --annotations='annovar:refGene,annovar:cosmic70,snpeff,tests/databases/annotations/current/hg19/clinvar_20210123.parquet' 

> howard annotation --input=tests/data/example.vcf.gz --output=/tmp/example.howard.tsv --assembly=hg19 --annotations='ALL:parquet' 

> howard annotation --input=tests/data/example.vcf.gz --output=/tmp/example.howard.tsv --param=config/param.json 

### Main options
```
--input=<input> | required

Input file path.
Format file must be either VCF, Parquet, TSV, CSV, PSV or duckDB.
Files can be compressesd (e.g. vcf.gz, tsv.gz).

```

```
--output=<output> | required

Output file path.
Format file must be either VCF, Parquet, TSV, CSV, PSV or duckDB.
Files can be compressesd (e.g. vcf.gz, tsv.gz).

```

```
--annotations=<annotations>

Annotation with databases files, or with tools,
as a list of files in Parquet, VCF, BED, or keywords.
For a Parquet/VCF/BED file, use file path (e.g. '/path/to/file.parquet').
For Annovar annotation, use keyword 'annovar' with annovar code (e.g. 'annovar:refGene', 'annovar:cosmic70').
For add all availalbe databases files, use 'ALL' keyword:
- Use 'ALL:<types>:<releases>'
- e.g. 'ALL', 'ALL:parquet:current', 'ALL:parquet,vcf:devel'
For snpeff annotation, use keyword 'snpeff'

```

```
--assembly=<assembly> (hg19)

Default assembly
(e.g. 'hg19'.

```

```
--annotations_update

Update option for annotation (Only for Parquet annotation).
If True, annotation fields will be removed and re-annotated.
These options will be applied to all annotation databases.

```

```
--annotations_append

Append option for annotation (Only for Parquet annotation).
If True, annotation fields will be annotated only if not annotation exists for the variant.
These options will be applied to all annotation databases.

```

```
--param=<param> ({})

Parameters JSON file or JSON string.

```



## CALCULATION tool
Calculation processes variants information to generate new information, such as: identify variation type (VarType), harmonizes allele frequency (VAF) and calculate sttistics (VAF_stats), extracts Nomen (transcript, cNomen, pNomen...) from an HGVS field (e.g. snpEff, Annovar) with an optional list of personalized transcripts, generates VaRank format barcode, identify trio inheritance.

Usage examples:

> howard calculation --input=tests/data/example.full.vcf --output=/tmp/example.calculation.tsv --calculations='vartype' 

> howard calculation --input=tests/data/example.ann.vcf.gz --output=/tmp/example.calculated.tsv --calculations='snpeff_hgvs,NOMEN' --hgvs_field=snpeff_hgvs --transcripts=tests/data/transcripts.tsv 

> howard calculation --input=tests/data/example.ann.vcf.gz --output=/tmp/example.ann.tsv --param=config/param.json 

> howard calculation --show_calculations 



### Main options
```
--input=<input>

Input file path.
Format file must be either VCF, Parquet, TSV, CSV, PSV or duckDB.
Files can be compressesd (e.g. vcf.gz, tsv.gz).

```

```
--output=<output>

Output file path.
Format file must be either VCF, Parquet, TSV, CSV, PSV or duckDB.
Files can be compressesd (e.g. vcf.gz, tsv.gz).

```

```
--calculations=<operations>

Calculations on genetic variants information and genotype information,
as a list of operations (e.g. 'VARTYPE,variant_id').
List of available calculations (unsensitive case, see doc for more information):
 VARTYPE  snpeff_hgvs  FINDBYPIPELINE  GENOTYPECONCORDANCE  BARCODE  TRIO  VAF  VAF_STATS  DP_STATS 

```

```
--calculation_config=<calculation config>

Calculation configuration JSON file.

```

```
--show_calculations

Show available calculation operations.

```

```
--param=<param> ({})

Parameters JSON file or JSON string.

```

### NOMEN calculation
```
--hgvs_field=<HGVS field> (hgvs)

HGVS INFO/tag containing a list o HGVS annotations.

```

```
--transcripts=<transcripts>

Transcripts TSV file,
with Transcript in first column, optional Gene in second column.

```

### TRIO calculation
```
--trio_pedigree=<trio pedigree>

Pedigree Trio for trio inheritance calculation.
either a JSON file or JSON string(e.g. '{"father": "sample1", "mother": "sample2", "child": "sample3"}').

```



## PRIORITIZATION tool
Prioritization algorithm uses profiles to flag variants (as passed or filtered), calculate a prioritization score, and automatically generate a comment for each variants (example: 'polymorphism identified in dbSNP. associated to Lung Cancer. Found in ClinVar database'). Prioritization profiles are defined in a configuration file in JSON format. A profile is defined as a list of annotation/value, using wildcards and comparison options (contains, lower than, greater than, equal...). Annotations fields may be quality values (usually from callers, such as 'DP') or other annotations fields provided by annotations tools, such as HOWARD itself (example: COSMIC, Clinvar, 1000genomes, PolyPhen, SIFT). Multiple profiles can be used simultaneously, which is useful to define multiple validation/prioritization levels (example: 'standard', 'stringent', 'rare variants', 'low allele frequency').



Usage examples:

> howard prioritization --input=tests/data/example.vcf.gz --output=/tmp/example.prioritized.vcf.gz --prioritizations=config/prioritization_profiles.json --profiles='default,GERMLINE' 

> howard prioritization --input=tests/data/example.vcf.gz --output=/tmp/example.prioritized.tsv --param=config/param.json 

### Main options
```
--input=<input> | required

Input file path.
Format file must be either VCF, Parquet, TSV, CSV, PSV or duckDB.
Files can be compressesd (e.g. vcf.gz, tsv.gz).

```

```
--output=<output> | required

Output file path.
Format file must be either VCF, Parquet, TSV, CSV, PSV or duckDB.
Files can be compressesd (e.g. vcf.gz, tsv.gz).

```

```
--prioritizations=<prioritisations>

Prioritization file in JSON format (defines profiles, see doc).

```

```
--param=<param> ({})

Parameters JSON file or JSON string.

```

### Prioritization
```
--profiles=<profiles>

List of prioritization profiles to process (based on Prioritization JSON file),
such as 'default', 'rare variants', 'low allele frequency', 'GERMLINE'.
By default, all profiles available will be processed.

```

```
--default_profile=<default profile>

Prioritization profile by default (see doc).
Default is the first profile in the list of prioritization profiles.

```

```
--pzfields=<pzfields> (PZScore,PZFlag)

Prioritization fields to provide (see doc).
Available: PZScore, PZFlag, PZTags, PZComment, PZInfos

```

```
--prioritization_score_mode=<prioritization score mode> ['HOWARD', 'VaRank'] (HOWARD)

Prioritization Score mode (see doc).
Available: HOWARD (increment score), VaRank (max score)

```



## PROCESS tool
howard process tool manage genetic variations to:

- annotates genetic variants with multiple annotation databases/files and tools

- calculates and normalizes annotations

- prioritizes variants with profiles (list of citeria) to calculate scores and flags

- translates into various formats

- query genetic variants and annotations

- generates variants statistics

Usage examples:

> howard process --input=tests/data/example.vcf.gz --output=/tmp/example.annotated.vcf.gz --param=config/param.json 

> howard process --input=tests/data/example.vcf.gz --annotations='snpeff' --calculations='snpeff_hgvs' --prioritizations=config/prioritization_profiles.json --explode_infos --output=/tmp/example.annotated.tsv --query='SELECT "#CHROM", POS, ALT, REF, snpeff_hgvs FROM variants' 

> howard process --input=tests/data/example.vcf.gz --hgvs='full_format,use_exon' --explode_infos --output=/tmp/example.annotated.tsv --query='SELECT "#CHROM", POS, ALT, REF, hgvs FROM variants' 

> howard process --input=tests/data/example.vcf.gz --output=/tmp/example.howard.vcf.gz --hgvs='full_format,use_exon' --annotations='tests/databases/annotations/current/hg19/avsnp150.parquet,tests/databases/annotations/current/hg19/dbnsfp42a.parquet,tests/databases/annotations/current/hg19/gnomad211_genome.parquet' --calculation='NOMEN' --explode_infos --query='SELECT NOMEN, REVEL_score, SIFT_score, AF AS 'gnomad_AF', ClinPred_score, ClinPred_pred FROM variants' 



### Main options
```
--input=<input> | required

Input file path.
Format file must be either VCF, Parquet, TSV, CSV, PSV or duckDB.
Files can be compressesd (e.g. vcf.gz, tsv.gz).

```

```
--output=<output> | required

Output file path.
Format file must be either VCF, Parquet, TSV, CSV, PSV or duckDB.
Files can be compressesd (e.g. vcf.gz, tsv.gz).

```

```
--param=<param> ({})

Parameters JSON file or JSON string.

```

### Quick Processes
```
--annotations=<annotations>

Annotation with databases files, or with tools,
as a list of files in Parquet, VCF, BED, or keywords.
For a Parquet/VCF/BED file, use file path (e.g. '/path/to/file.parquet').
For Annovar annotation, use keyword 'annovar' with annovar code (e.g. 'annovar:refGene', 'annovar:cosmic70').
For add all availalbe databases files, use 'ALL' keyword:
- Use 'ALL:<types>:<releases>'
- e.g. 'ALL', 'ALL:parquet:current', 'ALL:parquet,vcf:devel'
For snpeff annotation, use keyword 'snpeff'

```

```
--calculations=<operations>

Calculations on genetic variants information and genotype information,
as a list of operations (e.g. 'VARTYPE,variant_id').
List of available calculations (unsensitive case, see doc for more information):
 VARTYPE  snpeff_hgvs  FINDBYPIPELINE  GENOTYPECONCORDANCE  BARCODE  TRIO  VAF  VAF_STATS  DP_STATS 

```

```
--prioritizations=<prioritisations>

Prioritization file in JSON format (defines profiles, see doc).

```

```
--hgvs=<HGVS options>

Quick HGVS annotation options.
This option will skip all other hgvs options.
Examples:
- 'default' (for default options)
- 'full_format' (for full format HGVS annotation)
- 'use_gene:True,add_protein:true,codon_type:FULL'

```

```
--query=<query>

Query in SQL format
(e.g. 'SELECT * FROM variants LIMIT 50').

```

```
--explode_infos

Explode VCF INFO/Tag into 'variants' table columns.

```

```
--explode_infos_prefix=<explode infos prefix>

Explode VCF INFO/Tag with a specific prefix.

```

```
--explode_infos_fields=<explode infos list> (*)

Explode VCF INFO/Tag specific fields/tags.
Keyword `*` specify all available fields, except those already specified.
Pattern (regex) can be used, such as `.*_score` for fields named with '_score' at the end.
Examples:
- 'HGVS,SIFT,Clinvar' (list of fields)
- 'HGVS,*,Clinvar' (list of fields with all other fields at the end)
- 'HGVS,.*_score,Clinvar' (list of 2 fields with all scores in the middle)
- 'HGVS,.*_score,*' (1 field, scores, all other fields)
- 'HGVS,*,.*_score' (1 field, all other fields, all scores)

```

```
--include_header

Include header (in VCF format) in output file.
Only for compatible formats (tab-delimiter format as TSV or BED).

```



## HGVS tool
HGVS annotation using HUGO HGVS internation Sequence Variant Nomenclature (http://varnomen.hgvs.org/). Annotation refere to refGene and genome to generate HGVS nomenclature for all available transcripts. This annotation add 'hgvs' field into VCF INFO column of a VCF file.

Usage examples:

> howard hgvs --input=tests/data/example.full.vcf --output=/tmp/example.hgvs.vcf 

> howard hgvs --input=tests/data/example.full.vcf --output=/tmp/example.hgvs.tsv --param=config/param.json 

> howard hgvs --input=tests/data/example.full.vcf --output=/tmp/example.hgvs.vcf --full_format --use_exon 

>  



### Main options
```
--input=<input> | required

Input file path.
Format file must be either VCF, Parquet, TSV, CSV, PSV or duckDB.
Files can be compressesd (e.g. vcf.gz, tsv.gz).

```

```
--output=<output>

Output file path.
Format file must be either VCF, Parquet, TSV, CSV, PSV or duckDB.
Files can be compressesd (e.g. vcf.gz, tsv.gz).

```

```
--hgvs=<HGVS options>

Quick HGVS annotation options.
This option will skip all other hgvs options.
Examples:
- 'default' (for default options)
- 'full_format' (for full format HGVS annotation)
- 'use_gene:True,add_protein:true,codon_type:FULL'

```

```
--assembly=<assembly> (hg19)

Default assembly
(e.g. 'hg19'.

```

```
--param=<param> ({})

Parameters JSON file or JSON string.

```

### HGVS
```
--use_gene

Use Gene information to generate HGVS annotation
(e.g. 'NM_152232(TAS1R2):c.231T>C')
```

```
--use_exon

Use Exon information to generate HGVS annotation
(e.g. 'NM_152232(exon2):c.231T>C').
Only if 'use_gene' is not enabled.

```

```
--use_protein

Use Protein level to generate HGVS annotation
(e.g. 'NP_689418:p.Cys77Arg').
Can be used with 'use_exon' or 'use_gene'.

```

```
--add_protein

Add Protein level to DNA HGVS annotation (e.g 'NM_152232:c.231T>C,NP_689418:p.Cys77Arg').

```

```
--full_format

Generates HGVS annotation in a full format
by using all information to generates an exhaustive annotation
(non-standard, e.g. 'TAS1R2:NM_152232:NP_689418:c.231T>C:p.Cys77Arg').
Use 'use_exon' to add exon information
(e.g 'TAS1R2:NM_152232:NP_689418:exon2:c.231T>C:p.Cys77Arg').

```

```
--codon_type=<Codon type> ['1', '3', 'FULL'] (3)

Amino Acide Codon format type to use to generate HGVS annotation.
Available:
- '1': codon in 1 character (e.g. 'C', 'R')
- '3': codon in 3 character (e.g. 'Cys', 'Arg')
-'FULL': codon in full name (e.g. 'Cysteine', 'Arginine')

```

```
--refgene=<refGene>

Path to refGene annotation file.

```

```
--refseqlink=<refSeqLink>

Path to refSeqLink annotation file.

```

### Databases
```
--refseq-folder=<refseq folder> (/Users/lebechea/howard/databases/refseq/current)

Folder containing refSeq files.

```

```
--genomes-folder=<genomes> (/Users/lebechea/howard/databases/genomes/current)

Folder containing genomes.
(e.g. '/Users/lebechea/howard/databases/genomes/current'
```



## DATABASES tool
Download databases and needed files for howard and associated tools

Usage examples:

> howard databases --assembly=hg19 --download-genomes=~/howard/databases/genomes/current --download-genomes-provider=UCSC --download-genomes-contig-regex='chr[0-9XYM]+$' 

> howard databases --assembly=hg19 --download-annovar=~/howard/databases/annovar/current --download-annovar-files='refGene,cosmic70,nci60' 

> howard databases --assembly=hg19 --download-snpeff=~/howard/databases/snpeff/current 

> howard databases --assembly=hg19 --download-refseq=~/howard/databases/refseq/current --download-refseq-format-file='ncbiRefSeq.txt' 

> howard databases --assembly=hg19 --download-dbnsfp=~/howard/databases/dbnsfp/current --download-dbnsfp-release='4.4a' --download-dbnsfp-subdatabases 

> howard databases --assembly=hg19 --download-alphamissense=~/howard/databases/alphamissense/current 

> howard databases --assembly=hg19 --download-exomiser=~/howard/databases/exomiser/current 

> howard databases --assembly=hg19 --download-dbsnp=~/howard/databases/dbsnp/current --download-dbsnp-vcf 

> cd ~/howard/databases && howard databases --assembly=hg19 --download-genomes=genomes/current --download-genomes-provider=UCSC --download-genomes-contig-regex='chr[0-9XYM]+$' --download-annovar=annovar/current --download-annovar-files='refGene,cosmic70,nci60' --download-snpeff=snpeff/current --download-refseq=refseq/current --download-refseq-format-file='ncbiRefSeq.txt' --download-dbnsfp=dbnsfp/current --download-dbnsfp-release='4.4a' --download-dbnsfp-subdatabases --download-alphamissense=alphamissense/current --download-exomiser=exomiser/current --download-dbsnp=dbsnp/current --download-dbsnp-vcf --threads=8 

> howard databases --generate-param=/tmp/param.json --generate-param-description=/tmp/test.description.json --generate-param-formats=parquet 

Notes:

> - Downloading databases can take a while, depending on network, threads and memory

> - Proxy: Beware of network and proxy configuration

> - dbNSFP download: More threads, more memory usage (8 threads ~ 16Gb, 24 threads ~ 32Gb)



### Main options
```
--assembly=<assembly> (hg19)

Default assembly
(e.g. 'hg19'.

```

```
--genomes-folder=<genomes> (/Users/lebechea/howard/databases/genomes/current)

Folder containing genomes.
(e.g. '/Users/lebechea/howard/databases/genomes/current'
```

### Genomes
```
--download-genomes=<genomes>

Path to genomes folder
with Fasta files, indexes,
and all files generated by pygenome module.
(e.g. '/Users/lebechea/howard/databases/genomes/current').

```

```
--download-genomes-provider=<genomes provider> ['GENCODE', 'Ensembl', 'UCSC', 'NCBI'] (UCSC)

Download Genome from an external provider.
Available: GENCODE, Ensembl, UCSC, NCBI.

```

```
--download-genomes-contig-regex=<genomes contig regex>

Regular expression to select specific chromosome
(e.g 'chr[0-9XYM]+$').

```

### snpEff
```
--download-snpeff=<snpEff>

Download snpEff databases within snpEff folder
```

### Annovar
```
--download-annovar=<Annovar>

Path to Annovar databases
(e.g. '/Users/lebechea/howard/databases/annovar/current').

```

```
--download-annovar-files=<Annovar code>

Download Annovar databases for a list of Annovar file code (see Annovar Doc).
Use None to donwload all available files,
or Annovar keyword (e.g. 'refGene', 'cosmic70', 'clinvar_202*').
Note that refGene will at least be downloaded,
and only files that not already exist or changed will be downloaded.

```

```
--download-annovar-url=<Annovar url> (http://www.openbioinformatics.org/annovar/download)

Annovar databases URL (see Annovar Doc).

```

### refSeq
```
--download-refseq=<refSeq>

Path to refSeq databases
(e.g. '/Users/lebechea/howard/databases/refseq/current').

```

```
--download-refseq-url=<refSeq url> (http://hgdownload.soe.ucsc.edu/goldenPath)

refSeq databases URL (see refSeq WebSite)
(e.g. 'http://hgdownload.soe.ucsc.edu/goldenPath')•/n
```

```
--download-refseq-prefix=<refSeq prefix> (ncbiRefSeq)

Check existing refSeq files in refSeq folder.

```

```
--download-refseq-files=<refSeq files> (ncbiRefSeq.txt,ncbiRefSeqLink.txt)

List of refSeq files to download.

```

```
--download-refseq-format-file=<refSeq format file>

Name of refSeq file to convert in BED format
(e.g. 'ncbiRefSeq.txt').
Process only if not None.

```

```
--download-refseq-include-utr5

Formating BED refSeq file including 5'UTR.

```

```
--download-refseq-include-utr3

Formating BED refSeq file including 3'UTR.

```

```
--download-refseq-include-chrM

Formating BED refSeq file including Mitochondiral chromosome 'chrM' or 'chrMT'.

```

```
--download-refseq-include-non-canonical-chr

Formating BED refSeq file including non canonical chromosomes.

```

```
--download-refseq-include-non-coding-transcripts

Formating BED refSeq file including non coding transcripts.

```

```
--download-refseq-include-transcript-version

Formating BED refSeq file including transcript version.

```

### dbNSFP
```
--download-dbnsfp=<dbNSFP>

Download dbNSFP databases within dbNSFP folder(e.g. '/Users/lebechea/howard/databases').

```

```
--download-dbnsfp-url=<dbNSFP url> (https://dbnsfp.s3.amazonaws.com)

Download dbNSFP databases URL (see dbNSFP website)
(e.g. https://dbnsfp.s3.amazonaws.com').

```

```
--download-dbnsfp-release=<dnNSFP release> (4.4a)

Release of dbNSFP to download (see dbNSFP website)
(e.g. '4.4a').

```

```
--download-dbnsfp-parquet-size=<dbNSFP parquet size> (100)

Maximum size (Mb) of data files in Parquet folder.
Parquet folder are partitioned (hive) by chromosome (sub-folder),
which contain N data files.

```

```
--download-dbnsfp-subdatabases

Generate dbNSFP sub-databases.
dbNSFP provides multiple databases which are split onto multiple columns.
This option create a Parquet folder for each sub-database (based on columns names).

```

```
--download-dbnsfp-parquet

Generate a Parquet file for each Parquet folder.

```

```
--download-dbnsfp-vcf

Generate a VCF file for each Parquet folder.
Need genome FASTA file (see --download-genome).

```

```
--download-dbnsfp-no-files-all

Not generate database Parquet/VCF file for the entire database ('ALL').
Only sub-databases files will be generated.
(see '--download-dbnsfp-subdatabases').

```

```
--download-dbnsfp-add-info

Add INFO column (VCF format) in Parquet folder and file.
Useful for speed up full annotation (all available columns).
Increase memory and space during generation of files.

```

```
--download-dbnsfp-row-group-size=<dnNSFP row grooup size> (100000)

Minimum number of rows in a parquet row group (see duckDB doc).
Lower can reduce memory usage and slightly increase space during generation,
speed up highly selective queries, slow down whole file queries (e.g. aggregations).

```

### AlphaMissense
```
--download-alphamissense=<AlphaMissense>

Path to AlphaMissense databases
```

```
--download-alphamissense-url=<AlphaMissense url> (https://storage.googleapis.com/dm_alphamissense)

Download AlphaMissense databases URL (see AlphaMissense website)
(e.g. 'https://storage.googleapis.com/dm_alphamissense').

```

### Exomiser
```
--download-exomiser=<Exomiser>

Path to Exomiser databases
(e.g. /Users/lebechea/howard/databases/exomiser/current).

```

```
--download-exomiser-application-properties=<Exomiser application properties>

Exomiser Application Properties configuration file (see Exomiser website).
This file contains configuration settings for the Exomiser tool.
If this parameter is not provided, the function will attempt to locate
the application properties file automatically based on the Exomiser.
Configuration information will be used to download expected releases (if no other parameters).
CADD and REMM will be downloaded only if 'path' are provided.

```

```
--download-exomiser-url=<Exomiser url> (http://data.monarchinitiative.org/exomiser)

URL where Exomiser database files can be downloaded from
(e.g. 'http://data.monarchinitiative.org/exomiser').

```

```
--download-exomiser-release=<Exomiser release>

Release of Exomiser data to download.
If "default", "auto", or "config", retrieve from Application Properties file.
If not provided (None), from Application Properties file (Exomiser data-version) 
or default '2109'.

```

```
--download-exomiser-phenotype-release=<Exomiser phenoptye release>

Release of Exomiser phenotype to download.
If not provided (None), from Application Properties file (Exomiser Phenotype data-version)
or Exomiser release.

```

```
--download-exomiser-remm-release=<Exomiser remm release>

Release of ReMM (Regulatory Mendelian Mutation) database to download.
If "default", "auto", or "config", retrieve from Application Properties file.

```

```
--download-exomiser-remm-url=<Exomiser remm url> (https://kircherlab.bihealth.org/download/ReMM)

URL where ReMM (Regulatory Mendelian Mutation) database files can be downloaded from
(e.g. 'https://kircherlab.bihealth.org/download/ReMM').

```

```
--download-exomiser-cadd-release=<Exomiser cadd release>

Release of CADD (Combined Annotation Dependent Depletion) database to download.
If "default", "auto", or "config", retrieve from Application Properties file.

```

```
--download-exomiser-cadd-url=<Exomiser cadd url> (https://kircherlab.bihealth.org/download/CADD)

URL where CADD (Combined Annotation Dependent Depletion) database files can be downloaded from
(e.g. 'https://kircherlab.bihealth.org/download/CADD').

```

```
--download-exomiser-cadd-url-snv-file=<Exomiser url snv file> (whole_genome_SNVs.tsv.gz)

Name of the file containing the SNV (Single Nucleotide Variant) data
for the CADD (Combined Annotation Dependent Depletion) database.

```

```
--download-exomiser-cadd-url-indel-file=<Exomiser cadd url indel> (InDels.tsv.gz)

Name of the file containing the INDEL (Insertion-Deletion) data
for the CADD (Combined Annotation Dependent Depletion) database.

```

### dbSNP
```
--download-dbsnp=<dnSNP>

Path to dbSNP databases
(e.g. '/Users/lebechea/howard/databases/exomiser/dbsnp').

```

```
--download-dbsnp-releases=<dnSNP releases> (b156)

Release of dbSNP to download
(e.g. 'b152', 'b152,b156').

```

```
--download-dbsnp-release-default=<dnSNP release default>

Default Release of dbSNP ('default' symlink)
(e.g. 'b156').
If None, first release to download will be assigned as default
only if it does not exists.

```

```
--download-dbsnp-url=<dbSNP url> (https://ftp.ncbi.nih.gov/snp/archive)

URL where dbSNP database files can be downloaded from.
(e.g. 'https://ftp.ncbi.nih.gov/snp/archive').

```

```
--download-dbsnp-url-files=<dbSNP url files>

Dictionary that maps assembly names to specific dbSNP URL files.
It allows you to provide custom dbSNP URL files for specific assemblies
instead of using the default file naming convention.

```

```
--download-dbsnp-url-files-prefix=<dbSNP url files prefix> (GCF_000001405)

String that represents the prefix of the dbSNP file name for a specific assembly.
It is used to construct the full URL of the dbSNP file to be downloaded.

```

```
--download-dbsnp-assemblies-map=<dbSNP assemblies map> ({'hg19': '25', 'hg38': '40'})

dictionary that maps assembly names to their corresponding dbSNP versions.
It is used to construct the dbSNP file name based on the assembly name.

```

```
--download-dbsnp-vcf

Generate well-formatted VCF from downloaded file:
- Add and filter contigs associated to assembly
- Normalize by splitting multiallelics
- Need genome (see --download-genome)

```

```
--download-dbsnp-parquet

Generate Parquet file from VCF.

```

### HGMD
```
--convert-hgmd=<HGMD>

Convert HGMD databases.
Folder where the HGMD databases will be stored.
Fields in VCF, Parquet and TSV will be generated.
If the folder does not exist, it will be created.

```

```
--convert-hgmd-file=<HGMD file>

File from HGMD.
Name format 'HGMD_Pro_<release>_<assembly>.vcf.gz'.

```

```
--convert-hgmd-basename=<HGMD basename>

File output basename.
Generated files will be prefixed by basename
(e.g. 'HGMD_Pro_MY_RELEASE')
By default (None), input file name without '.vcf.gz'.

```

### Parameters file
```
--generate-param=<param>

Parameter file (JSON) with all databases found.
Databases folders scanned are defined in config file.
Structure of databases follow this structure (see doc):
.../<database>/<release>/<assembly>/*.[parquet|vcf.gz|...]

```

```
--generate-param-description=<param description>

Description file (JSON) with all databases found.
Contains all databases with description of format, assembly, fields...

```

```
--generate-param-releases=<param release> (current)

List of database folder releases to check
(e.g. 'current', 'latest').

```

```
--generate-param-formats=<param formats> (parquet)

List of database formats to check
(e.g. 'parquet', 'parquet,vcf,bed,tsv').

```

```
--generate-param-bcftools

Generate parameter JSON file with BCFTools annotation for allowed formats
(i.e. 'vcf', 'bed').

```



## FROM_ANNOVAR tool
(beta) Formatting Annovar database file to other format (VCF and Parquet). Exported Parquet file includes INFO/tags columns as VCF INFO columns had been exploded

Usage examples:

> howard from_annovar --input=tests/databases/others/hg19_nci60.txt --output=/tmp/nci60.from_annovar.vcf.gz --to_parquet=/tmp/nci60.from_annovar.parquet --annovar-code=nci60 --genome=~/howard/databases/genomes/current/hg19.fa --threads=8 

### Main options
```
--input=<input> | required

Input file path.
Format file must be either VCF, Parquet, TSV, CSV, PSV or duckDB.
Files can be compressesd (e.g. vcf.gz, tsv.gz).

```

```
--output=<output> | required

Output file path.
Format file must be either VCF, Parquet, TSV, CSV, PSV or duckDB.
Files can be compressesd (e.g. vcf.gz, tsv.gz).

```

```
--genome=<genome> (~/howard/databases/genomes/current/hg19/hg19.fa) | required

Genome file in fasta format (e.g. 'hg19.fa', 'hg38.fa').

```

### Annovar
```
--annovar-code=<Annovar code>

Annovar code, or database name.
Usefull to name databases columns.

```

### Parquet
```
--to_parquet=<to parquet>

Parquet file conversion.

```

### Modes
```
--reduce_memory=<reduce memory> ['auto', 'enable', 'disable'] (auto)

Reduce memory option,
either 'auto' (auto-detection), 'enable' or 'disable'.

```

```
--multi_variant=<multi variant> ['auto', 'enable', 'disable'] (auto)

Variant with multiple annotation lines.
Either 'auto' (auto-detection), 'enable' or 'disable'.

```



## GUI tool
Graphical User Interface tools

Usage examples:

> howard gui 



## HELP tool
Help tools

Usage examples:

> howard help --help_md=docs/help.md --help_html=docs/help.html

> howard help --help_json_input=docs/help.config.json --help_json_input_title='HOWARD Configuration' --help_md=docs/help.config.md --help_html=docs/help.config.html --code_type='json'

> howard help --help_json_input=docs/help.param.json --help_json_input_title='HOWARD Parameters' --help_md=docs/help.param.md --help_html=docs/help.param.html --code_type='json'

### Main options
```
--help_md=<help markdown>

Help Output file in MarkDown format.

```

```
--help_html=<help html>

Help Output file in HTML format.

```

```
--help_json_input=<help JSON input>

Help input file in JSON format.

```

```
--help_json_input_title=<help JSON input title> (Help)

Help JSON input title.

```

```
--code_type=<example code type>

Help example code type for input JSON format
(e.g. 'json', 'bash').

```



## Shared arguments
```
--config=<config> ({})

Configuration JSON file or JSON string.

```

```
--threads=<threads> (-1)

Specify the number of threads to use for processing HOWARD.
It determines the level of parallelism,
either on python scripts, duckdb engine and external tools.
It and can help speed up the process/tool.
Use -1 to use all available CPU/cores.
Either non valid value is 1 CPU/core.

```

```
--memory=<memory>

Specify the memory to use in format FLOAT[kMG]
(e.g. '8G', '12.42G', '1024M').
It determines the amount of memory for duckDB engine and external tools
(especially for JAR programs).
It can help to prevent 'out of memory' failures.
By default (None) is 80%% of RAM (for duckDB).

```

```
--chunk_size=<chunk size> (1000000)

Number of records in batch to export output file.
The lower the chunk size, the less memory consumption.
For Parquet partitioning, files size will depend on the chunk size.

```

```
--tmp=<Temporary folder>

Temporary folder (e.g. '/tmp').
By default, '.tmp' for duckDB (see doc),external tools and python scripts.

```

```
--duckdb_settings=<duckDB settings>

DuckDB settings (see duckDB doc) as JSON (string or file).
These settings have priority (see options 'threads', 'tmp'...).
Examples: '{"TimeZone": "GMT", "temp_directory": "/tmp/duckdb", "threads": 8}'.

```

```
--verbosity=<verbosity> ['CRITICAL', 'ERROR', 'WARNING', 'INFO', 'DEBUG', 'NOTSET'] (INFO)

Verbosity level
Available: CRITICAL, ERROR, WARNING, INFO, DEBUG or NOTSET
- DEBUG: Detailed information, typically of interest only when diagnosing problems.
- INFO: Confirmation that things are working as expected.
- WARNING: An indication that something unexpected happened.
- ERROR: Due to a more serious problem.
- CRITICAL: A serious error.
- NOTSET: All messages.

```

```
--log=<log>

Logs file
(e.g. 'my.log').

```

```
--quiet

==SUPPRESS==
```

```
--verbose

==SUPPRESS==
```

```
--debug

==SUPPRESS==
```

