# :bulb: HOWARD Tips


<details><summary>How to hive partitioning into Parquet a huge VCF of SNV?</summary>

This tip describe hive partitioning of huge VCF files into Parquet, to prevent memory resource crash. The following bash script is spliting VCF by "#CHROM" and "chunk" files of a defined size (number of lines "NBLINE"). Depending on input VCF file, the number of lines will determine the size of each "chunk" file. For a "chunk" file round 200-300Mb, a number of line for a VCF with few annotations is around 10.000.000, but for a VCF with a lot of annotations it is around 1.000.000. This will ensure a memory usage around 5-6Gb.

Moreover, an additional partitioning can be applied, on one or more specific VCF column or INFO annotation: "REF,ALT" for VCF only with SNV (e.g. dbSNP); "CLINSIG" for ClinVar database (values such as "Pathogenic, Moslty-Pathogenic"...). Note that partitioning only works only for small values (a value lenght will create a too long partition folder), and for not too many values for a column (partition fragments is maximum of 1024). This additional partition will create an additional "chunk" partition (can be used as a column in Parquet), and can take a while.

In order to create a high-performance database for HOWARD, INFO annotations can be exploded in columns. This option is slower, but generate a Parquet database that will be used by picking needed columns for annotation. Annotations to explode can be chosen with more options.

```bash
# Files
VCF=/tmp/my.vcf.gz                  # Input VCF file
PARQUET=/tmp/my.partition.parquet      # Output Partitioned Parquet folder

# Tools
BCFTOOLS=bcftools   # BCFTools
BGZIP=bgzip         # BGZip
HOWARD=howard       # HOWARD
TABIX=tabix         # Tabix

# Threads
THREADS=12          # Number of threads

# Param
NBLINE="10000000"                   # 10000000 for few annotations, 1000000 for a lot of annotations
CONVERT_OPTIONS=" --explode_infos " # Explode INFO annotations into columns
PARQUET_PARTITIONS=""               # "REF,ALT" (SNV VCF) or "REF" or "CLNSIG"
CHUNK_NAME="chunk"                  # Chunk column name

# Create output folder
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
        $BCFTOOLS filter $VCF -r $chr --threads $THREADS | $BCFTOOLS view -H --threads $THREADS | split -a 10 -l $NBLINE - $PARQUET/#CHROM=$chr/ --filter="$BGZIP -l1 --threads=$THREADS > \$FILE.gz";
        nb_chunk_files=$(ls $PARQUET/#CHROM=$chr/*.gz | wc -l)

        # Convert chunk VCF to Parquet
        i_chunk_files=0
        for file in $PARQUET/#CHROM=$chr/*.gz; do
            
            # Chunk file to TSV and header
            ((i_chunk_files++))
            chunk=$(basename $file | sed s/.gz$//gi)
            echo "# Chromosome '$chr' - Convert VCF to Parquet '$chunk' [$i_chunk_files/$nb_chunk_files]..."
            mv $file $file.tsv.gz;
            cp $PARQUET/header.vcf $file.tsv.gz.hdr;

            # Convert with partitioning or not
            if [ "$PARQUET_PARTITIONS" != "" ]; then
                file_chunck_base=$(dirname $file)/$CHUNK_NAME=$chunk;
                $HOWARD convert --input=$file.tsv.gz --output=$file_chunck_base.parquet $CONVERT_OPTIONS --threads=$THREADS --parquet_partitions="$PARQUET_PARTITIONS" --verbosity="ERROR";
                mv $file_chunck_base.parquet $file_chunck_base;
            else
                $HOWARD convert --input=$file.tsv.gz --output=$file.parquet $CONVERT_OPTIONS --threads=$THREADS --verbosity="ERROR";
            fi;

            # Clean
            rm -f $file.tsv*

        done;
    fi;
done;

# Create header
cp $PARQUET/header.vcf $PARQUET.hdr

# Show partitioned Parquet folder
tree $PARQUET -h

```

</details>
