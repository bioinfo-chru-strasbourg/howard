# HOWARD To Excel plugin

Convert input file to Excel '.xlsx'

## Examples

> Convert VCF file into Excel file, with header and variants view as sheets
```bash
howard to_excel --input="tests/data/example.vcf.gz" --output="/tmp/example.xlsx" --add_variants_view --add_header
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

## Add options

```
--add_variants_view   

Create a sheet with all INFO fields exploded.
(default: False)

```

```
--add_header

Create a sheet with all INFO fields header descritions.
(default: False)

```