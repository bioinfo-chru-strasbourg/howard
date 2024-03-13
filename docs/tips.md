# :bulb: HOWARD Tips


<details><summary>How to hive partitioning into Parquet a huge VCF of SNV?</summary>

Due to memory usage with duckDB, huge VCF file convertion can fail. This tip describe hive partitioning of huge VCF files into Parquet, to prevent memory resource crash.

The following bash script is spliting VCF by "#CHROM" and chunk VCF files with a defined size ("CHUNK_SIZE"). Depending on input VCF file, the number of chunk VCF files will be determined by the content (usually number of annotation within the VCF file). Moreover, Parquet files can also be chunked ("CHUNK_SIZE_PARQUET"). Default options This will ensure a memory usage around 4-5Gb.

Moreover, an additional partitioning can be applied, on one or more specific VCF column or INFO annotation: "None" for no partitioning; "REF,ALT" for VCF only with SNV (e.g. dbSNP); "CLINSIG" for ClinVar database (values such as "Pathogenic, Moslty-Pathogenic"...). Note that partitioning only works for small values (a value lenght will create a too long partition folder), and for not too many values for a column (partition fragments is maximum of 1024). This additional partition can take a while.

In order to create a high-performance database for HOWARD, INFO annotations can be exploded in columns. This option is slower, but generate a Parquet database that will be used by picking needed columns for annotation. Annotations to explode can be chosen with more options.

```bash
# Files
VCF=/tmp/my.vcf.gz                  # Input VCF file
PARQUET=/tmp/my.partition.parquet   # Output Partitioned Parquet folder

# Tools
BCFTOOLS=bcftools   # BCFTools
BGZIP=bgzip         # BGZip
HOWARD=howard       # HOWARD
TABIX=tabix         # Tabix

# Threads
THREADS=12          # Number of threads

# Param
CHUNK_SIZE=1000000000               # 1000000000 for VCF chunk size around 200Mb
CHUNK_SIZE_PARQUET=10000000         # 10000000 for parquet chunk size around 200Mb
PARQUET_PARTITIONS="None"           # "None" for no more partition, "REF,ALT" (SNV VCF) or "REF" or "CLNSIG"...
CONVERT_OPTIONS=" --explode_infos " # Explode INFO annotations into columns

# Create output folder
rm -r $PARQUET
mkdir -p $PARQUET

# Extract header
$BCFTOOLS view -h $VCF --threads $THREADS > $PARQUET/header.vcf

# VCF indexing (if necessary)
if [ ! -e $VCF.tbi ]; then
    $TABIX $VCF
fi;

# For each chromosome
for chr in $($TABIX -l $VCF | cut -f1); do

    if [ "$chr" != "None" ]; then
        echo "# Chromosome '$chr'"

        # Create chromosome folder
        mkdir -p $PARQUET/#CHROM=$chr;

        echo "# Chromosome '$chr' - BCFTools filter and split file..."
        $BCFTOOLS filter $VCF -r $chr --threads $THREADS | $BCFTOOLS view -H --threads $THREADS | split -a 10 -C $CHUNK_SIZE - $PARQUET/#CHROM=$chr/ --filter="$BGZIP -l1 --threads=$THREADS > \$FILE.gz";
        nb_chunk_files=$(ls $PARQUET/#CHROM=$chr/*.gz | wc -l)

        # Convert chunk VCF to Parquet
        i_chunk_files=0
        for file in $PARQUET/#CHROM=$chr/*.gz; do
            
            # Chunk file to TSV and header
            ((i_chunk_files++))
            chunk=$(basename $file | sed s/.gz$//gi)
            echo "# Chromosome '$chr' - Convert VCF to Parquet '$chunk' ($PARQUET_PARTITIONS) [$i_chunk_files/$nb_chunk_files]..."
            mv $file $file.tsv.gz;
            cp $PARQUET/header.vcf $file.tsv.gz.hdr;

            # Convert with partitioning or not
            $HOWARD convert --input=$file.tsv.gz --output=$file.parquet $CONVERT_OPTIONS --threads=$THREADS --parquet_partitions="$PARQUET_PARTITIONS" --verbosity=ERROR --chunk_size=$CHUNK_SIZE_PARQUET

            # Redirect partitioninf folder if any
            if [ "$PARQUET_PARTITIONS" != "" ]; then
                rsync -a $file.parquet/ $(dirname $file)/
                rm -r $file.parquet*
            fi;

            # Clean
            rm -f $file.tsv*

        done;
    fi;
done;

# Create header
cp $PARQUET/header.vcf $PARQUET.hdr

# Show partitioned Parquet folder
tree -h $PARQUET

```

</details>
