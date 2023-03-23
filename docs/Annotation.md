# Annotation

HOWARD annotates VCF.

## Tools

Annovar, snpEff, BCFTOOLS and Parquet/DuckDB

## Databases

### Formats

- Parquet / DuckDB

- VCF

- BED

- snpEff bin

- Annovar Databases

### Download

- Annovar: ...

- snpEff: ...

- VCF: ...

- BED: UCSC...

- Parquet: Prepare from VCF

### Preparation

From Annovar to VCF

example:
```
annovar_folder=/databases/annotations/current/hg19/
genome=databases/genomes/current/hg19.fa
python howard/tools/annovar_to_vcf.py --format=vcf --input=$annovar_folder/hg19_avsnp150.txt --output=$annovar_folder/avsnp150.vcf.gz --genome=$genome --database_name=avsnp150 --threads=8;
```

If Annovar (or VCF) with multiple line for variants
example (avsnp150):
```
1	10352	10352	-	A	rs145072688
1	10352	10352	-	A	rs555500075
```
```
python howard/tools/vcf_uniq_variant.py $annovar_folder/avsnp150.vcf.gz $annovar_folder/avsnp150.uniquify.vcf.gz;
```


From VCF to Parquet/DuckDB

```
howard --input=db.vcf.gz --output=db.parquet --param='{"connexion_type": "db.duckdb", "explode_infos": true, "export_extra_infos": true}'
```
example:
```
annovar_folder=/databases/annotations/current/hg19/
python -m howard.main --debug  --input=$annovar_folder/avsnp150.vcf.gz --output=$annovar_folder/avsnp150.parquet --param='{"connexion_type": "$annovar_folder/avsnp150.duckdb", "explode_infos": true, "export_extra_infos": true}' --threads=8;
```


