# HOWARD GeneBe plugin

PyGeneBe: A Python client seamlessly integrating with the GeneBe platform, offering efficient annotation of genetic variants through its API, while supporting pandas, VCF file formats, and HGVS parsing

Using this client, you can easily annotate your DNA variants with the GeneBe API. Annotations include:

Gene, transcript, and effect
ClinVar phenotype
GnomAD frequency
ACMG score

See [GeneBe website](https://genebe.net/) and [GeneBe GitHub](https://github.com/pstawinski/pygenebe) for more information


## Installation

GeneBe needs Python packages. Use Pip to install requirements

```bash
python -m pip install -r plugins/genebe/requirements.txt
```

## Examples

> Annotate with GeneBe
```bash
howard genebe --input="tests/data/example.vcf.gz" --output="/tmp/example.genebe.vcf.gz"
```
```ts
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	sample1	sample2	sample3	sample4
chr1	28736	.	A	C	100	PASS	CLNSIG=pathogenic;transcript=NR_024540.1;gene_symbol=WASH7P;phylop100way_score=-3.6419999599456787;acmg_score=0;acmg_classification=Uncertain_significance;hgvs_c=n.;consequences=intron_variant	GT:AD:DP:GQ	0/1:525,204:729:99	0/1:12659,4994:17664:99	1/1:12658,4995:17663:99	1/1:401,175:576:99
chr1	35144	.	A	C	100	PASS	CLNSIG=non-pathogenic;transcript=NR_026818.1;gene_symbol=FAM138A;phylop100way_score=-1.2669999599456787;acmg_score=-4;acmg_classification=Likely_benign;hgvs_c=n.597T>G;consequences=non_coding_transcript_exon_variant;bayesdelnoaf_score=-0.9700000286102295;acmg_criteria=BP4_Strong	GT:AD:DP:GQ	./.	0/1:12659,4994:17664:99	0/1:12658,4995:17663:99	0/1:401,175:576:99
chr1	69101	.	A	G	100	PASS	DP=50;transcript=NM_001005484.2;gene_symbol=OR4F5;phylop100way_score=0.6230000257492065;acmg_score=-1;acmg_classification=Likely_benign;gene_hgnc_id=14825.0;hgvs_c=c.74A>G;consequences=missense_variant;bayesdelnoaf_score=-0.6000000238418579;acmg_criteria=BP4;revel_score=0.07599999755620956;alphamissense_score=0.210099995136261	GT:AD:DP:GQ	0/1:525,204:729:99	./.	0/1:12658,4995:17663:99	0/1:401,175:576:99
chr1	768251	.	A	G	100	PASS	transcript=NR_047519.1;gene_symbol=LINC01128;phylop100way_score=0.289000004529953;acmg_score=-3;acmg_classification=Likely_benign;hgvs_c=n.287+3767T>G;consequences=intron_variant;bayesdelnoaf_score=-0.9300000071525574;acmg_criteria=PM2_Supporting,BP4_Strong;dbsnp=1557557080	GT:AD:DP:GQ	0/1:525,204:729:99	./.	0/1:12658,4995:17663:99	0/1:401,175:576:99
chr1	768252	.	A	G	100	PASS	transcript=NR_047519.1;gene_symbol=LINC01128;phylop100way_score=0.28700000047683716;acmg_score=-3;acmg_classification=Likely_benign;hgvs_c=n.287+3768C>G;consequences=intron_variant;bayesdelnoaf_score=-0.9300000071525574;acmg_criteria=PM2_Supporting,BP4_Strong	GT:AD:DP:GQ	0/1:525,204:729:99	./.	0/1:12658,4995:17663:99	0/1:401,175:576:99
chr1	768253	.	A	G	100	PASS	transcript=NR_047519.1;gene_symbol=LINC01128;phylop100way_score=-0.6769999861717224;acmg_score=-3;acmg_classification=Likely_benign;hgvs_c=n.287+3769A>G;consequences=intron_variant;bayesdelnoaf_score=-0.9300000071525574;acmg_criteria=PM2_Supporting,BP4_Strong	GT:AD:DP:GQ	0/1:525,204:729:99	./.	0/1:12658,4995:17663:99	0/1:401,175:576:99
chr7	55249063	rs1050171	G	A	5777	PASS	DP=125;transcript=NM_005228.5;gene_symbol=EGFR;phylop100way_score=4.2270002365112305;acmg_score=-18;acmg_classification=Benign;gene_hgnc_id=3236.0;hgvs_c=c.2361G>A;consequences=synonymous_variant;bayesdelnoaf_score=-0.2800000011920929;acmg_criteria=BP4_Moderate,BP6_Very_Strong,BA1;dbsnp=1050171;gnomad_exomes_af=0.562637984752655;gnomad_genomes_af=0.5142030119895935;gnomad_exomes_ac=822502.0;gnomad_genomes_ac=78200.0;gnomad_exomes_homalt=236892.0;gnomad_genomes_homalt=21009.0;clinvar_disease=not specified,Lung carcinoma,Squamous cell lung carcinoma,EGFR-related lung cancer,not provided,Inflammatory skin and bowel disease, neonatal, 2,Lung cancer;clinvar_classification=Benign	GT:AD:DP:GQ	0/1:525,204:729:99	0/1:12659,4994:17664:99	./.	0/1:401,175:576:99
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

```
--assembly=<assembly> (default: hg19)

Genome Assembly (e.g. 'hg19', 'hg38').


```