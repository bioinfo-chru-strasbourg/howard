#################################
## HOWARD environment
#################################


# ANNOVAR and SNPEff configuration
####################################

# see config.ini file


## MANDATORY for multithreading
#################################

# Main tools folder
NGS_TOOLS=/home/TOOLS/tools
NGS_DATABASES=/home/TOOLS/databases

# THREADS
THREADS=$(ls -d /sys/devices/system/cpu/cpu[[:digit:]]* | wc -w)

# TOOLS
########

# HTSlib - A C library for reading/writing high-throughput sequencing data
export HTSLIB_DESCRIPTION="A C library for reading/writing high-throughput sequencing data"
export HTSLIB_REF="http://www.htslib.org/"

# TABIX
export TABIX=$NGS_TOOLS/htslib/1.8/bin/tabix			# BIN
export TABIX_PATH=$(dirname $TABIX)				# BIN
export TABIX_VERSION=1.8					# VER
export TABIX_DESCRIPTION="Indexing VCF files"
export TABIX_REF=$HTSLIB_REF

# BGZIP
export BGZIP=$NGS_TOOLS/htslib/1.8/bin/bgzip			# BIN
export BGZIP_VERSION=1.8					# VER
export BGZIP_DESCRIPTION="Compressing VCF files"
export BGZIP_REF=$HTSLIB_REF

# BCFTOOLS
export BCFTOOLS=$NGS_TOOLS/bcftools/1.8/bin/bcftools		# BIN
export BCFTOOLS_VERSION=1.8					# VER
export BCFTOOLS_DESCRIPTION="Reading/writing BCF2/VCF/gVCF files and calling/filtering/summarising SNP and short indel sequence variants"
export BCFTOOLS_REF="http://www.htslib.org/"

# SNPEFF
export SNPEFF_FOLDER=$NGS_TOOLS/snpeff/4.3t/bin		# FOLDER
export SNPEFF=$SNPEFF_FOLDER/snpEff.jar			# BIN-JAR
export SNPEFF_VERSION=4.3t				# VER
export SNPEFF_CONFIG=$SNPEFF_FOLDER/snpeff.config	# CONFIG # NOT USED !!! # CHANGE CONFIG FILE in SNPEFF TOOL if necessary
export SNPEFF_DATABASES=$NGS_FOLDER/snpeff_sources4.3	# DATA # NOT USED !!! # CHANGE DATABASE location in CONFIG FILE in SNPEFF TOOL if necessary
export SNPEFF_DESCRIPTION="Genetic variant annotation and effect prediction toolbox. It annotates and predicts the effects of variants on genes (such as amino acid changes)"
export SNPEFF_REF="A program for annotating and predicting the effects of single nucleotide polymorphisms, SnpEff: SNPs in the genome of Drosophila melanogaster strain w1118; iso-2; iso-3., Cingolani P, Platts A, Wang le L, Coon M, Nguyen T, Wang L, Land SJ, Lu X, Ruden DM. Fly (Austin). 2012 Apr-Jun;6(2):80-92 "
TOOLS_LIST=$TOOLS_LIST" SNPEFF"

# ANNOVAR
export ANNOVAR=$NGS_TOOLS/annovar/2018Apr16/bin                         # DIR
export ANNOVAR_VERSION=2018Apr16                                                        # VER
export ANNOVAR_DATABASES=$NGS_DATABASES/annovar_sources         # DATA
export ANNOVAR_DESCRIPTION="an efficient software tool to utilize update-to-date information to functionally annotate genetic variants detected from diverse genomes"
export ANNOVAR_REF="Wang K, Li M, Hakonarson H. ANNOVAR: Functional annotation of genetic variants from next-generation sequencing data Nucleic Acids Research, 38:e164, 2010"

# JAVA
export JAVA=$NGS_TOOLS/java/1.8.0/bin/java                                      # BIN
export JAVA_PATH=$NGS_TOOLS/java/1.8.0/bin                                      # BIN
export JAVA_VERSION=1.8.0                                                                       # VER
export JAVA_DESCRIPTION="A high-level programming language developed by Sun Microsystems"
export JAVA_REF="http://java.com"
