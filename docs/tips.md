# HOWARD Tips

<details>

<summary>
How to annotate a VCF file with a Parquet database?
</summary>

Annotate a VCF file is simple using `annotation` tool.
Number of threads and memory can be handled.

``` bash
INPUT=tests/data/example.vcf.gz
OUTPUT=/tmp/example.annotated.vcf.gz
PARQUET=tests/databases/annotations/current/hg19/nci60.parquet
howard annotation \
   --input="$INPUT" \
   --output="$OUTPUT" \
   --annotations="$PARQUET" \
   --threads=8 --memory=4G
```

</details>

<details>

<summary>
How to speedup annotation and reduce memory for a huge input VCF and/or a huge Parquet database?
</summary>

For huge input VCF file and/or huge annotation Parquet database, and if it's the uniq process,
two specific parameters can speedup annotation and reduce memory usage:

- `fast`: speedup execution to focus on Parquet annotation without any other process
- `access`: read only strategy will not load input data but read the input file

The strategy consist on converting the input VCF file into Parquet, and then annotate it.

``` bash
INPUT=tests/data/example.vcf.gz
OUTPUT=/tmp/example.annotated.vcf.gz
PARQUET=tests/databases/annotations/current/hg19/nci60.parquet
# 1. Convert input VCF file into Parquet
howard convert \
   --input="$INPUT" \
   --output="$OUTPUT.parquet"
   --access='RO' \
   --threads=2 --memory=4G
# 2. Annotate input Parquet file
howard annotation \
   --input="$OUTPUT.parquet" \
   --output="$OUTPUT" \
   --annotations="$PARQUET" \
   --access='RO' \
   --fast \
   --threads=2 --memory=4G
```

</details>

<details>

<summary>
How to hive partitioning into Parquet a VCF format file?
</summary>

In order to create a database from a VCF file, process partitioning
highly speed up futher annotation process. Simply use HOWARD Convert
tool with `--parquet_partitions` option. Format of input and output
files can be any avaialbe format (e.g. Parquet, VCF, TSV).

This process is not ressource intensive, but can take a while for huge
databases. However, using `--explode_infos` option (annotations in columns)
requires much more memory for huge databases.

``` bash
INPUT=~/howard/databases/dbsnp/current/hg19/b156/dbsnp.vcf.gz
OUTPUT=~/howard/databases/dbsnp/current/hg19/b156/dbsnp.partition.parquet
PARTITIONS="#CHROM" # "#CHROM", "#CHROM,REF,ALT" (for SNV file only) 
howard convert \
   --input=$INPUT \
   --output=$OUTPUT \
   --parquet_partitions="#CHROM" \
   --threads=8
```

</details>

<details>

<summary>
How to process a huge input VCF file?
</summary>
To process a huge input VCF file efficiently, a strategy constits on partitioning it into multiple smaller Parquet files, process each partition individually, and then merge the results into a final output file. This approach helps manage memory usage and speeds up processing for large datasets. The `chunking` strategy is implemented into `process` tool:

```bash
# Initialize variables
input_file=tests/data/example.vcf
output_file=/tmp/example.test.vcf
chunking_size=2 # default 1000000
howard process \
    --input=$input_file \
    --output=$output_file \
    --calculations="VARTYPE" \
    --chunking_enable \
    --chunking_size=$chunk_size
```

This method ensures efficient handling of large files while maintaining flexibility for various processing tasks.

</details>

<details>

<summary>
How to aggregate all INFO annotations from multiple Parquet databases
into one INFO field?
</summary>

In order to merge all annotations in INFO column of multiple databases,
use a SQL query on the list of Parquet databases, and use `GROUP_CONCAT`
duckDB function to aggretate values.
This method can requires huge memory if Parquet files are huge.

``` bash
howard query \
   --explode_infos \
   --explode_infos_prefix='INFO/' \
   --query="SELECT \"#CHROM\", POS, REF, ALT, GROUP_CONCAT(INFO, ';') AS INFO \
            FROM parquet_scan('tests/databases/annotations/current/hg19/*.parquet', union_by_name = true) \
            GROUP BY \"#CHROM\", POS, REF, ALT" \
   --output=/tmp/full_annotation.tsv
   
head -n2 /tmp/full_annotation.tsv
```

> Expected result (depends on Parquet databases)
``` csv
#CHROM  POS     REF     ALT     INFO
chr1    69093   G       T       MCAP13=0.001453;REVEL=0.117;SIFT_score=0.082;SIFT_converted_rankscore=0.333;SIFT_pred=T;SIFT4G_score=0.097;SIFT4G_converted_rankscore=0.392;SIFT4G_pred=T;Polyphen2_HDIV_score=0.0;Polyphen2_HDIV_rankscore=0.029;Polyphen2_HDIV_pred=B;Polyphen2_HVAR_score=0.0;Polyphen2_HVAR_rankscore=0.014;Polyphen2_HVAR_pred=B;LRT_score=0.589;LRT_converted_rankscore=0.056;LRT_pred=N;MutationTaster_score=0.635;MutationTaster_converted_rankscore=0.328;MutationTaster_pred=D;MutationAssessor_score=.;MutationAssessor_rankscore=.;MutationAssessor_pred=.;FATHMM_score=6.74;FATHMM_converted_rankscore=0.005;FATHMM_pred=T;PROVEAN_score=0.27;PROVEAN_converted_rankscore=0.043;PROVEAN_pred=N;VEST4_score=0.12;VEST4_rankscore=0.111;MetaSVM_score=-1.003;MetaSVM_rankscore=0.291;MetaSVM_pred=T;MetaLR_score=0.002;MetaLR_rankscore=0.005;MetaLR_pred=T;MetaRNN_score=0.452;MetaRNN_rankscore=0.666;MetaRNN_pred=T;M_CAP_score=0.001;M_CAP_rankscore=0.022;M_CAP_pred=T;REVEL_score=0.117;REVEL_rankscore=0.332;MutPred_score=0.835;MutPred_rankscore=0.943;MVP_score=0.240;MVP_rankscore=0.236;MPC_score=.;MPC_rankscore=.;PrimateAI_score=.;PrimateAI_rankscore=.;PrimateAI_pred=.;DEOGEN2_score=0.008;DEOGEN2_rankscore=0.072;DEOGEN2_pred=T;BayesDel_addAF_score=-0.359;BayesDel_addAF_rankscore=0.042;BayesDel_addAF_pred=T;BayesDel_noAF_score=-0.754;BayesDel_noAF_rankscore=0.033;BayesDel_noAF_pred=T;ClinPred_score=0.090;ClinPred_rankscore=0.112;ClinPred_pred=T;LIST_S2_score=0.065;LIST_S2_rankscore=0.651;LIST_S2_pred=T;Aloft_pred=.,.,;Aloft_Confidence=.,.,;CADD_raw=0.013;CADD_raw_rankscore=0.049;CADD_phred=2.83;DANN_score=0.602;DANN_rankscore=0.064;fathmm_MKL_coding_score=0.505;fathmm_MKL_coding_rankscore=0.287;fathmm_MKL_coding_pred=D;fathmm_XF_coding_score=0.012;fathmm_XF_coding_rankscore=0.001;fathmm_XF_coding_pred=N;Eigen_raw_coding=-0.958;Eigen_raw_coding_rankscore=0.095;Eigen_PC_raw_coding=-0.968;Eigen_PC_raw_coding_rankscore=0.105;GenoCanyon_score=0;GenoCanyon_rankscore=0.012;integrated_fitCons_score=0.487;integrated_fitCons_rankscore=0.14;integrated_confidence_value=0;LINSIGHT=.;LINSIGHT_rankscore=.;GERP___NR=2.31;GERP___RS=1.36;GERP___RS_rankscore=0.211;phyloP100way_vertebrate=1.139;phyloP100way_vertebrate_rankscore=0.311;phyloP30way_mammalian=0.113;phyloP30way_mammalian_rankscore=0.17;phastCons100way_vertebrate=0.841;phastCons100way_vertebrate_rankscore=0.303;phastCons30way_mammalian=0.552;phastCons30way_mammalian_rankscore=0.281;SiPhy_29way_logOdds=6.575;SiPhy_29way_logOdds_rankscore=0.218;Interpro_domain=.,.;GTEx_V8_gene=.;GTEx_V8_tissue=.;InterVar_automated=Uncertain_significance;PVS1=0;PS1=0;PS2=0;PS3=0;PS4=0;PM1=0;PM2=1;PM3=0;PM4=0;PM5=0;PM6=0;PP1=0;PP2=0;PP3=0;PP4=0;PP5=0;BA1=0;BS1=0;BS2=0;BS3=0;BS4=0;BP1=0;BP2=0;BP3=0;BP4=1;BP5=0;BP6=0;BP7=0
```

</details>

<details>

<summary>
How to explore genetics variations from VCF files?
</summary>

[CuteVariant: A standalone and free application to explore genetics
variations from VCF files](https://cutevariant.labsquare.org/)

Cutevariant is a cross-plateform application dedicated to maniupulate
and filter variation from annotated VCF file. When you create a project,
data are imported into an sqlite database that cutevariant queries
according your needs. Presently, SnpEff and VEP annotations are
supported. Once your project is created, you can query variant using
different gui controller or directly using the VQL language. This Domain
Specific Language is specially designed for cutevariant and try to keep
the same syntax than SQL for an easy use.

Published in [Bioinformatics Advanced - Cutevariant: a standalone
GUI-based desktop application to explore genetic variations from an
annotated VCF
file](https://academic.oup.com/bioinformaticsadvances/article/2/1/vbab028/6440032?login=true)

Documentation available
on [cutevariant.labsquare.org](https://cutevariant.labsquare.org/) and
[GitHub]((https://github.com/labsquare/cutevariant))

<figure>

<img
src="https://raw.githubusercontent.com/labsquare/cutevariant/master/screencast.gif"
title="HOWARD Graphical User Interface" alt="CuteVariant" />
<figcaption aria-hidden="true">

CuteVariant
</figcaption>

</figure>

</details>
