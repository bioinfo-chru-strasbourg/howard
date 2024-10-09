# HOWARD Minimalize plugin

Minimalize a VCF file consists in put missing value ('.') on INFO/Tags, ID, QUAL or FILTER fields. Options can also minimalize samples (keep only GT) or remove all samples. INFO/tags can by exploded before minimalize to keep tags into separated columns (useful for Parquet or TSV format to constitute a database).

## Examples

> Minimalize all informations to keep only variants and genotype of samples
```bash
howard minimalize --input="tests/data/example.vcf.gz" --output="/tmp/example.minimal.vcf.gz" --minimalize_info --minimalize_filter --minimalize_qual --minimalize_id --minimalize_samples
```
```ts
#CHROM  POS       ID  REF  ALT  QUAL  FILTER  INFO  FORMAT  sample1  sample2  sample3  sample4
chr1    28736     .   A    C    .     .       GT    0/1     0/1      1/1      1/1
chr1    35144     .   A    C    .     .       GT    ./.     0/1      0/1      0/1
chr1    69101     .   A    G    .     .       GT    0/1     ./.      0/1      0/1
chr1    768251    .   A    G    .     .       GT    0/1     ./.      0/1      0/1
chr1    768252    .   A    G    .     .       GT    0/1     ./.      0/1      0/1
chr1    768253    .   A    G    .     .       GT    0/1     ./.      0/1      0/1
chr7    55249063  .   G    A    .     .       GT    0/1     0/1      ./.      0/1
```
> Remove samples and explode INFOS/Tags into columns (and minimalize INFO field)
```bash
howard minimalize --input="tests/data/example.vcf.gz" --output="/tmp/example.minimal.tsv" --remove_samples --explode_infos --minimalize_info
```
```ts
#CHROM  POS       REF  ALT  ID         QUAL  FILTER  INFO  AA              CLNSIG  DP  NS  SIFT
chr1    28736     A    C    .          100   PASS    .     pathogenic
chr1    35144     A    C    .          100   PASS    .     non-pathogenic
chr1    69101     A    G    .          100   PASS    .     50
chr1    768251    A    G    .          100   PASS    .
chr1    768252    A    G    .          100   PASS    .
chr1    768253    A    G    .          100   PASS    .
chr7    55249063  G    A    rs1050171  5777  PASS    .     125
```

## Main options
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
--param=<param> (default: {})

Parameters JSON file (or string) defines parameters to process 
annotations, calculations, prioritizations, convertions and queries.

```

## Minimalize options

```
--minimalize_info

Minimalize INFO field (e.g. '.' value).
(default: False)

```

```
--minimalize_id

Minimalize ID field (e.g. '.' value).
(default: False)

```

```
--minimalize_qual

Minimalize QUAL field (e.g. '' value).
(default: False)

```

```
--minimalize_filter

Minimalize FILTER field (e.g. '.' value).
(default: False)

```

```
--minimalize_samples

Minimalize samples to keep only genotypes (i.e. 'GT').
(default: False)

```

```
--remove_samples

Remove all samples to keep only variants.
(default: False)

```
