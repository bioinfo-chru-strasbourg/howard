#!/bin/bash
#################################
## HOWARD 
#################################

SCRIPT_NAME="HOWARD"
SCRIPT_DESCRIPTION="HOWARD Annotation, Calculation, Prioritization and Translation, based on ANNOVAR and snpEff, allowing multithreading"
SCRIPT_RELEASE="0.9.13.1b"
SCRIPT_DATE="17/12/2018"
SCRIPT_AUTHOR="Antony Le Bechec"
SCRIPT_COPYRIGHT="IRC"
SCRIPT_LICENCE="GNU-GPL"

# Realse note
RELEASE_NOTES=$RELEASE_NOTES"# 0.9b-07/10/2016:\n";
RELEASE_NOTES=$RELEASE_NOTES"#\tScript creation\n";
RELEASE_NOTES=$RELEASE_NOTES"# 0.9.1b-11/10/2016:\n";
RELEASE_NOTES=$RELEASE_NOTES"#\tAdd Prioritization and Translation\n";
RELEASE_NOTES=$RELEASE_NOTES"# 0.9.1b-11/10/2016:\n";
RELEASE_NOTES=$RELEASE_NOTES"#\tAdd snpEff annotation and stats\n";
RELEASE_NOTES=$RELEASE_NOTES"# 0.9.8b-21/03/2017:\n";
RELEASE_NOTES=$RELEASE_NOTES"#\tAdd Multithreading on Prioritization and Translation\n";
RELEASE_NOTES=$RELEASE_NOTES"# 0.9.9b-18/04/2017:\n";
RELEASE_NOTES=$RELEASE_NOTES"#\tAdd Calculation step\n";
RELEASE_NOTES=$RELEASE_NOTES"# 0.9.10b-07/11/2017:\n";
RELEASE_NOTES=$RELEASE_NOTES"#\tAdd generic file annotation through --annotation option\n";
RELEASE_NOTES=$RELEASE_NOTES"#\t\tNo need to be in configuration file\n";
RELEASE_NOTES=$RELEASE_NOTES"#\t\tNeed to be in ANNOVAR database folder (file 'ASSEMBLY_ANN.txt' for annotation 'ANN')\n";
RELEASE_NOTES=$RELEASE_NOTES"#\tAdd options: --force , --split\n";
RELEASE_NOTES=$RELEASE_NOTES"#\tAdd options for VCFanotation.pl: --show_annoataion, --show_annotations_full\n";
RELEASE_NOTES=$RELEASE_NOTES"#\tAdd database download option nowget in VCFanotation.pl\n";
RELEASE_NOTES=$RELEASE_NOTES"#\tFixes: multithreading, VAF calculation, configuration and check dependencies\n";
RELEASE_NOTES=$RELEASE_NOTES"# 0.9.11b-07/05/2018:\n";
RELEASE_NOTES=$RELEASE_NOTES"#\tReplace VCFTOOLS command to BCFTOOLS command\n";
RELEASE_NOTES=$RELEASE_NOTES"#\tRelease added into the output VCF\n";
RELEASE_NOTES=$RELEASE_NOTES"#\tUpdate SNPEff options\n";
RELEASE_NOTES=$RELEASE_NOTES"#\tAdd VARTYPE, CALLING_QUALITY and CALLING_QUALITY_EXPLODE option on calculation\n";
RELEASE_NOTES=$RELEASE_NOTES"#\tAdd description on calculations\n";
RELEASE_NOTES=$RELEASE_NOTES"# 0.9.11.1b-14/05/2018:\n";
RELEASE_NOTES=$RELEASE_NOTES"#\tImprove VCF validation\n";
RELEASE_NOTES=$RELEASE_NOTES"#\tFix snpEff annotation bug\n";
RELEASE_NOTES=$RELEASE_NOTES"# 0.9.11.2b-17/08/2018:\n";
RELEASE_NOTES=$RELEASE_NOTES"#\tAdd --vcf input vcf file option\n";
RELEASE_NOTES=$RELEASE_NOTES"#\tCreate Output file directory automatically\n";
RELEASE_NOTES=$RELEASE_NOTES"#\tImprove Multithreading\n";
RELEASE_NOTES=$RELEASE_NOTES"# 0.9.12b-24/08/2018:\n";
RELEASE_NOTES=$RELEASE_NOTES"#\tImprove Multithreading\n";
RELEASE_NOTES=$RELEASE_NOTES"#\tInput VCF compressed with BGZIP accepted\n";
RELEASE_NOTES=$RELEASE_NOTES"#\tOutput VCF compression level\n";
RELEASE_NOTES=$RELEASE_NOTES"#\tAdd VCF input sorting and multiallele split step (by default)\n";
RELEASE_NOTES=$RELEASE_NOTES"#\tAdd VCF input normalization step with option --norm\n";
RELEASE_NOTES=$RELEASE_NOTES"#\tBug fixes\n";
RELEASE_NOTES=$RELEASE_NOTES"# 0.9.13b-04/10/2018:\n";
RELEASE_NOTES=$RELEASE_NOTES"#\tMultithreading improved\n";
RELEASE_NOTES=$RELEASE_NOTES"#\tChange default output vcf\n";
RELEASE_NOTES=$RELEASE_NOTES"#\tInput vcf without samples allowed\n";
RELEASE_NOTES=$RELEASE_NOTES"#\tVCF Validation with contig check\n";
RELEASE_NOTES=$RELEASE_NOTES"#\tAdd multi VCF in input option\n";
RELEASE_NOTES=$RELEASE_NOTES"#\tAdd --annotate option for BCFTOOLS annotation with a VCF and TAG (beta)\n";
RELEASE_NOTES=$RELEASE_NOTES"#\tRemove no multithreading part code to multithreading with 1 thread\n";
RELEASE_NOTES=$RELEASE_NOTES"#\tRemove --multithreading parameter, only --thread parameter to deal with multithreading\n";
RELEASE_NOTES=$RELEASE_NOTES"#\tReplace --filter and --format parameters by --prioritization and --translation parameters\n";
RELEASE_NOTES=$RELEASE_NOTES"#\tAdd snpeff options to VCFannotation.pl\n";
RELEASE_NOTES=$RELEASE_NOTES"# 0.9.13.1b-13/12/2018:\n";
RELEASE_NOTES=$RELEASE_NOTES"#\tChange Number/Type/Description of new INFO/FORMAT header generated\n";
RELEASE_NOTES=$RELEASE_NOTES"#\tRemove snpEff option --snpeff and --snpeff_hgvs. SnpEff is used through --annotation option\n";
RELEASE_NOTES=$RELEASE_NOTES"#\tAdd '#' to the TAB delimiter format header\n";
RELEASE_NOTES=$RELEASE_NOTES"#\tBug fixed: calculation INFO fields header, snpeff parameters options on multithreading\n";
RELEASE_NOTES=$RELEASE_NOTES"#\tBug fixed: snpeff parameters in command line\n";


# Script folder
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# Header
function header () {
	echo "#######################################";
	echo "# $SCRIPT_NAME [$SCRIPT_RELEASE-$SCRIPT_DATE]";
	echo "# $SCRIPT_DESCRIPTION ";
	echo "# $SCRIPT_AUTHOR @ $SCRIPT_COPYRIGHT © $SCRIPT_LICENCE";
	echo "#######################################";
}

# Release
function release () {
	echo "# RELEASE NOTES:";
	echo -e $RELEASE_NOTES
}

# Usage
function usage {
	echo "# USAGE: $(basename $0) --input=<FILE>  [options...]";
	echo "# Following options are available:";
	echo "# --input|vcf=<FILES>                Input file in VCF format (*vcf or *vcf.gz in BGZIP compression format)";
	echo "# --output=<FILE>                    Output annotated file in defined format (see --format option). Default 'output.vcf'";
	echo "# --annotate=<LIST>                  Annotation with BCFTOOLS: List of VCF files with TAG to annotate.";
	echo "#                                    Format 'VCF:TAG;VCF:TAG...' (e.g. 'annotate1.vcf:ID,QUAL,+TAG;annotate2.vcf:INFO/ANN'). Default TAG '+INFO'";
	echo "# --annotation=<LIST>                Annotation: List of annotation (in the order to add into the input VCF)";
	echo "# --calculation=<LIST>               Calculation: List of calculation";
	echo "# --prioritization=<LIST OF FILE>    Priorization: List of prioritization config file";
	echo "# --translation=<STRING>             Translation: Output format (default VCF)";
	echo "# --config=<FILE>                    Configuration file (ANNOVAR and SNPEff)";
	#echo "# --multithreading                   Multithreading: perfom multithreading";
	echo "# --threads=<INTEGER>                Threads: number of thread to use (default defined in environment variable THREADS, or 1)";
	echo "# --env=<FILE>                       Environment configuration for multithreading (BGZIP, TABIX, BCFTOOLS)";
	echo "# --split=<INTEGER>                  Split by group of variants (default (10000)";
	echo "# --compress=<INTEGER>               Compression level output file *vcf.gz (0 to 9, -1 no compression by default)";
	echo "# --norm=<FILE>                      Genome fasta file to normalize (beware of chromosome identification, either 'x' or 'chrx')";
	echo "# --force                            Force annotation even if already exists in VCF header";
	echo "# --tmp=<FOLDER>                     Temporary folder (default /tmp)";
	echo "# --verbose                          VERBOSE option";
	echo "# --debug                            DEBUG option";
	echo "# --release                          RELEASE option";
	echo "# --help                             HELP option";
	echo "#";
	echo "# More options are also available for each steps:";
	echo "# VCFannotation.pl --help";
	echo "# VCFcalculation.pl --help";
	echo "# VCFprioritization.pl --help";
	echo "# VCFtranslation.pl --help";
	echo "";

	#$SCRIPT_DIR/VCFannotation.pl --release
	#$SCRIPT_DIR/VCFannotation.pl --help
	#$SCRIPT_DIR/VCFcalculation.pl --release
	#$SCRIPT_DIR/VCFcalculation.pl --help
	#$SCRIPT_DIR/VCFprioritization.pl --release
	#$SCRIPT_DIR/VCFprioritization.pl --help
	#$SCRIPT_DIR/VCFtranslation.pl --release
	#$SCRIPT_DIR/VCFtranslation.pl --help

}

# EXAMPLE :
# ./HOWARD.sh --input=test/IND8.200k.ann.vcf --output=test/IND8.200k.ann.test10000.txt --annotation=ALL --calculation=NOMEN,VAF,VAF_STAT  --filter=default --fields=NOMEN,PZFlag,PZScore,VAF_average --column=run:EXOME,sample:IND8 --format=tab --multithreading --transcripts=test/transcripts.txt --force --tmp=/home1/IRC/DATA/DEV/TMP/

# header
header;


ARGS=$(getopt -o "i:o:e:a:f:s:r:xt:m:vdnh" --long "input:,vcf:,output:,env:,annotation:,annotate:,calculation:,filter:,prioritization:,format:,translation:,snpeff_stats:,multithreading,threads:,split:,compress:,norm:,tmp:,verbose,debug,release,help" -- "$@" 2> /dev/null)
#ARGS=$(getopt --long "input:,output:,annotation:,multithreading,threads:,verbose,debug,release,help" -- "$@" 2> /dev/null)
if [ $? -ne 0 ]; then
	:
	#echo $?
	#usage;
	#exit;
fi;

PARAM=$@

eval set -- "$ARGS"
while true
do
	#echo "$1=$2"
	#echo "Eval opts";
	case "$1" in
		--input|--vcf)
			INPUT="$2";
			shift 2
			;;
		--output)
			OUTPUT="$2";
			shift 2
			;;
		--env)
			ENV="$2";
			shift 2
			;;
		--annotation)
			ANNOTATION="$2"
			shift 2
			;;
		--annotate)
			ANNOTATE="$2"
			shift 2
			;;
		--calculation)
			CALCULATION="$2"
			shift 2
			;;
		--config)
			CONFIG="$2"
			shift 2
			;;
		--prioritization|--filter)
			FILTER="$2"
			shift 2
			;;
		--translation|--format)
			FORMAT="$2"
			shift 2
			;;
		--snpeff_stats)
			SNPEFF_STATS="$2"
			shift 2
			;;
		--multithreading)
			MULTITHREADING=1
			shift 1
			;;
		--threads)
			THREADS_INPUT="$2"
			shift 2
			;;
		--split)
			SPLIT_INPUT="$2"
			shift 2
			;;
		--compress)
			COMPRESS="$2"
			shift 2
			;;
		--norm)
			NORM="$2"
			shift 2
			;;
		--tmp)
			TMP_INPUT="$2"
			shift 2
			;;
		--verbose)
			VERBOSE=1
			shift 1
			;;
		--debug)
			VERBOSE=1
			DEBUG=1
			shift 1
			;;
		--release)
			release;
			exit 0
			;;
		--help)
			usage
			exit 0
			;;
		--) shift
			break
			;;
		*) 	echo "# Option $1 is not recognized. " "Use -h or --help to display the help." && \
			exit 1
			;;
	esac
done

# ENV
if [ ! -z $ENV ] && [ -s $ENV ] && [ "$ENV" != "" ] && [ ! -d $ENV ]; then
	ENV=$ENV;
	echo "#[INFO] ENV '$ENV' found."
elif [ -s $SCRIPT_DIR/$ENV ] && [ "$ENV" != "" ] && [ ! -d $ENV ]; then
	ENV=$SCRIPT_DIR/$ENV;
	echo "#[INFO] ENV '$ENV' found."
elif [ "$ENV" == "" ] || [ ! -s $ENV ]; then
	if [ -s $SCRIPT_DIR/"env.sh" ]; then
		ENV=$SCRIPT_DIR/"env.sh";
		echo "#[INFO] Default ENV '$ENV' used."
	else
		ENV="";
		echo "#[WARNING] NO ENV defined. No ENV used."
	fi;
fi;
if [ -e $ENV ] && [ "$ENV" != "" ]; then
	source $ENV
fi;


# Mandatory parameters
if [ -z "$INPUT" ]; then #  || [ -z $OUTPUT ]; then
	echo "#[WARNING] INPUT Missing";
	usage;
	exit;
fi


# Mandatory parameters
if [ -z $OUTPUT ] || [ "$OUTPUT" == "" ]; then #  || [ -z $OUTPUT ]; then
	echo "#[WARNING] OUTPUT Missing";
	#OUTPUT=$(echo $INPUT | sed "s/.vcf$/.annotated.vcf/g"); #/\.vcf$/\.output\.vcf/;
	OUTPUT="output.vcf"; #/\.vcf$/\.output\.vcf/;
	echo "#[INFO] Default OUTPUT '$OUTPUT' used."
fi
# mkdir OUPUT directory 
mkdir -p $(dirname $OUTPUT)


# Multithreading
#if [ -z $MULTITHREADING ]; then
#	MULTITHREADING=1;
#	THREADS=1;
#fi;

#[[ $THREADS == ?(-)+([0-9]) ]] && echo "$1 is an integer"

# Threads
if [ ! -z $THREADS_INPUT ]; then
	THREADS=$THREADS_INPUT;
fi;
if [ -z $THREADS ] || ! [[ $THREADS =~ ^[0-9]+$ ]] || [ $THREADS -lt 1 ]; then
	THREADS=1;
fi;
MULTITHREADING=1;


# TMP
if [ ! -z $TMP_INPUT ] && [ ! -d $TMP_INPUT ]; then
	mkdir -p $TMP_INPUT;
fi;
if [ ! -z $TMP_INPUT ] && [ -d $TMP_INPUT ]; then
	TMP_SYS_FOLDER=$TMP_INPUT;
fi;
if [ -z $TMP_SYS_FOLDER ] || [ ! -d $TMP_SYS_FOLDER ]; then
	TMP_SYS_FOLDER=/tmp
	(($VERBOSE)) && echo "#[WARNING] TMP=$TMP_SYS_FOLDER"
	#usage;
	#exit;
fi


# Split
SPLIT_DEFAULT=10000
if [ ! -z $SPLIT_INPUT ]; then
	SPLIT=$SPLIT_INPUT;
fi;
if [ -z $SPLIT ]; then
	SPLIT=$SPLIT_DEFAULT;
fi;

# Compression
if [ -z $COMPRESS ]; then
	COMPRESS=-1;
fi;

# Normalization
if [ -z $NORM ]; then
	NORM="";
fi;


# ANNOTATIONS_split
if [ -z $ANNOTATIONS_SPLIT ]; then
	ANNOTATIONS_SPLIT=$SPLIT;
	#ANNOTATIONS_SPLIT=2;
fi;

# CALCULATION_split
if [ -z $CALCULATION_SPLIT ]; then
	CALCULATION_SPLIT=$SPLIT;
	#ANNOTATIONS_SPLIT=2;
fi;

# PRIORITIZATION_split
if [ -z $PRIORITIZATION_SPLIT ]; then
	PRIORITIZATION_SPLIT=$SPLIT;
	#ANNOTATIONS_SPLIT=2;
fi;


# TRANSLATION_split
if [ -z $TRANSLATION_SPLIT ]; then
	TRANSLATION_SPLIT=$SPLIT;
	#ANNOTATIONS_SPLIT=2;
fi;

# SPLIT SUFFIX LENGTH
SPLIT_SUFFIT_LENGTH=4




# CHECK dependencies
######################


if (($MULTITHREADING)); then

	echo "#[INFO] Multithreading ($THREADS threads)"

	# BGZIP
	if [ -z $BGZIP ] || ! (( $(command -v $BGZIP | grep ^ -c) )); then
		#echo "file not exists";
		BGZIP="bgzip"
		echo "#[INFO] Default BGZIP 'bgzip' configured";
	fi;
	#echo "$BGZIP"
	#if [ -z $BGZIP ] || [ ! -f $BGZIP ]; then
	if [ -z $BGZIP ] || ! (( $(command -v $BGZIP | grep ^ -c) )); then
		echo "#[ERROR] No BGZIP tool found. Add BGZIP in configuration file";
		exit 1;
	fi;
	echo "#[INFO] BGZIP=$BGZIP";

	# TABIX
	#echo "tabix: $TABIX"
	if [ -z $TABIX ] || ! (( $(command -v $TABIX | grep ^ -c) )); then
		#echo "file not exists";
		TABIX="tabix"
		echo "#[INFO] Default TABIX 'tabix' configured";
	fi;
	#echo "$TABIX"
	#if [ -z $TABIX ] || [ ! -f $TABIX ]; then
	if [ -z $TABIX ] || ! (( $(command -v $TABIX | grep ^ -c) )); then
		echo "#[ERROR] No TABIX tool found. Add TABIX in configuration file";
		exit 1;
	fi;
	echo "#[INFO] TABIX=$TABIX";

	#if [ -z $VCFTOOLS/vcf-sort ] || ! (( $(command -v $VCFTOOLS/vcf-sort | grep ^ -c) )); then
	#	echo "#[ERROR] No VCFTOOLS tool found. Add VCFTOOLS in configuration file";
	#	exit 1;
	#fi;
	#echo "#[INFO] VCFTOOLS=$VCFTOOLS";

	# BCFTOOLS
	#echo "BCFTOOLS: $BCFTOOLS"
	if [ -z $BCFTOOLS ] || ! (( $(command -v $BCFTOOLS | grep ^ -c) )); then
		#echo "file not exists";
		BCFTOOLS="bcftools"
		echo "#[INFO] Default BCFTOOLS 'bcftools' configured";
	fi;
	#echo "$BCFTOOLS"
	#if [ -z $BCFTOOLS ] || [ ! -f $BCFTOOLS ]; then
	if [ -z $BCFTOOLS ] || ! (( $(command -v $BCFTOOLS | grep ^ -c) )); then
		echo "#[ERROR] No BCFTOOLS tool found. Add BCFTOOLS in configuration file";
		exit 1;
	fi;
	echo "#[INFO] BCFTOOLS=$BCFTOOLS";
	
	# SnpEff
	echo "#[INFO] SNPEFF=$SNPEFF";

fi;


# TMP_FILES
TMP_FILES=""
TMP_VAR=$RANDOM$RANDOM
TMP_FOLDER=$TMP_SYS_FOLDER/HOWARD_$TMP_VAR
LOG=$TMP_FOLDER/log
ERR=$TMP_FOLDER/err
mkdir -p $TMP_FOLDER
touch $LOG;
touch $ERR;

PARAM=$PARAM" --tmp=$TMP_FOLDER "

TMP_SORT=$TMP_FOLDER/sort
mkdir -p $TMP_SORT

# CHECK INPUT VCF
######################

# VCF input file VALIDATION
if ((1)); then

	##################
	# VCF VALIDATION #
	##################

	echo -e "\n####################\n# VCF VALIDATION\n####################\n"

	# Mult VCF in input
	if [ $(echo $INPUT | tr "," "\n" | wc -w) -gt 0 ]; then
		#echo "#[INFO] List of VCF as input"
		INPUT_ORIGINAL_LIST=""
		INPUT_LIST=""
		NB_VCF=0;
		for INPUT_VCF in $(echo $INPUT | tr "," "\n"); do

			((NB_VCF++))
			echo -e "\n#[INFO] Validation of VCF #$NB_VCF '$INPUT_VCF'"; 
						
			if [ -e $INPUT_VCF ]; then
			
				# CP file in TMP
				INPUT_VCF_NUM=$TMP_FOLDER/INPUT_VCF_$NB_VCF.vcf
				cp $INPUT_VCF $INPUT_VCF_NUM
			
				# Not empty
				VCF_VALIDATION=$TMP_FOLDER/VCF_VALIDATION.$RANDOM
				>$VCF_VALIDATION.tmp
				>$VCF_VALIDATION
				
				# Simple check
				VCFnotEmpty=$($BCFTOOLS view $INPUT_VCF_NUM -H 2>$VCF_VALIDATION | head -n 1 | wc -l)
				
				if (($(grep ^Failed $VCF_VALIDATION -c))); then
					echo "#[WARNING] Input VCF '$INPUT_VCF' loading failed";
				#elif !(($VCFnotEmpty)); then
				#	echo "#[WARNING] Input VCF '$INPUT_VCF' is empty";
				else
				
					if !(($VCFnotEmpty)); then
						echo "#[WARNING] Input VCF '$INPUT_VCF' is empty";
					fi;
				
					# CONTIG
					if ((1)); then
						existContig=$($BCFTOOLS view -h $INPUT_VCF_NUM 2>>$ERR | grep "##contig=" -c)
						existVariants=$($BCFTOOLS view -H $INPUT_VCF_NUM 2>>$ERR | grep ^ -m1 -c)
						
						if (($existVariants)) && ! (($existContig)) ; then
							echo "#[INFO] Input VCF '$INPUT_VCF': No contig in header. Try to add contig in header"; 
							INPUT_TMP_CONTIG=$INPUT_VCF_NUM.contig.vcf
							INPUT_TMP_CONTIG_HEADER=$INPUT_VCF_NUM.contig.vcf.header
							# Create new header
							$BCFTOOLS view -h $INPUT_VCF_NUM 2>>$ERR | grep ^## >> $INPUT_TMP_CONTIG_HEADER 2>>$ERR
							$BCFTOOLS view -H $INPUT_VCF_NUM 2>>$ERR | cut -f1 | uniq | awk '{print "##contig=<ID="$1">"}' >> $INPUT_TMP_CONTIG_HEADER 2>>$ERR
							$BCFTOOLS view -h $INPUT_VCF_NUM 2>>$ERR | grep ^#CHROM >> $INPUT_TMP_CONTIG_HEADER 2>>$ERR
							$BCFTOOLS view -H $INPUT_VCF_NUM >> $INPUT_TMP_CONTIG 2>>$ERR
							# Reheader
							$BCFTOOLS reheader $INPUT_VCF_NUM -h $INPUT_TMP_CONTIG_HEADER > $INPUT_TMP_CONTIG 2>>$ERR
							# move file
							rm -f $INPUT_VCF_NUM
							mv $INPUT_TMP_CONTIG $INPUT_VCF_NUM
							rm -f $INPUT_TMP_CONTIG $INPUT_TMP_CONTIG_HEADER
						fi;
					fi;
				
					# VCF validation	
			
					VALIDATION_OK=1
				
					# Générate validation file
					TMP_SORT_VALIDATION=$TMP_SORT$RANDOM
					mkdir -p $TMP_SORT_VALIDATION
					$BCFTOOLS view $INPUT_VCF_NUM --threads $THREADS 2>>$VCF_VALIDATION.tmp | $BCFTOOLS sort -T $TMP_SORT_VALIDATION 1>/dev/null 2>>$VCF_VALIDATION.tmp
					rm -rf $TMP_SORT_VALIDATION

					# Extract warning and errors
					grep -v -e "^Writing to " -e "^Merging .* temporary files$" -e "^Cleaning$" -e "^Done$" $VCF_VALIDATION.tmp > $VCF_VALIDATION
					rm $VCF_VALIDATION.tmp

					# Test warning/error
					if (($(grep ^Failed $VCF_VALIDATION -c))) || (($(grep ^Aborded $VCF_VALIDATION -c))); then
						(($VERBOSE)) && echo "#[ERROR] Input VCF $INPUT failed";
						(($DEBUG)) && cat $VCF_VALIDATION;
						VALIDATION_OK=0;
					elif (($(grep "^\[E::" $VCF_VALIDATION -c))); then
						(($VERBOSE)) && echo "#[ERROR] Some Errors on INPUT VCF $INPUT_VCF:";
						(($DEBUG)) && cat $VCF_VALIDATION;
						VALIDATION_OK=0;
					elif (($(cat $VCF_VALIDATION | wc -l))); then
						(($VERBOSE)) && echo "#[WARNING] Some Errors/Warnings on INPUT VCF $INPUT_VCF:";
						(($DEBUG)) && cat $VCF_VALIDATION;
						VALIDATION_OK=1;
					else
						(($VERBOSE)) && echo "#[INFO] INPUT VCF validated";
						VALIDATION_OK=1;
					fi;
			
				
					if (($VALIDATION_OK)); then
				
						# Sorting and normalize
						echo "#[INFO] Input VCF '$INPUT_VCF' sorting and normalization"; 
						TMP_SORT_PREPARATION=$TMP_SORT$RANDOM
						mkdir -p $TMP_SORT_PREPARATION
						$BCFTOOLS norm -m- --threads $THREADS $INPUT_VCF_NUM 2>>$ERR | $BCFTOOLS sort -T $TMP_SORT_PREPARATION 2>>$ERR | $BGZIP -@ $THREADS -l 0 > $INPUT_VCF_NUM.gz 2>>$ERR
						$TABIX $INPUT_VCF_NUM.gz 1>>$LOG 2>>$ERR
						INPUT_VCF_NUM=$INPUT_VCF_NUM.gz
						if [ -e $INPUT_VCF_NUM ] && [ -e $INPUT_VCF_NUM.tbi ]; then
							INPUT_LIST=$INPUT_LIST" "$INPUT_VCF_NUM
							INPUT_ORIGINAL_LIST=$INPUT_ORIGINAL_LIST$INPUT_VCF" "
							#echo "#[INFO] INPUT VCF validated";
						else
							echo "#[WARNING] Input VCF '$INPUT_VCF' loading failed (sort and normalization failed)";
						fi;
					
				
					else
						echo "#[WARNING] Input VCF '$INPUT_VCF' loading failed (VCF not validated)";
				
					fi;
				fi;
			else
				echo "#[WARNING] Input VCF '$INPUT_VCF' loading failed (file does not exist)";
			fi;
		done;
		
		# Check Empty list
		if [ "$INPUT_LIST" == "" ]; then
			echo "#[ERROR] No validated VCF found in the input list: '$INPUT'";
			exit 0;
		fi;
		
		# List of Vaidated VCF
		echo -e "\n#[INFO] List of validated VCF: "$INPUT_ORIGINAL_LIST
		
		# Merging
		INPUT_MULTI_VCF=$TMP_FOLDER/INPUT_MULTI_VCF.vcf
		if [ $(echo $INPUT_LIST | wc -w) -gt 1 ]; then
			echo "#[INFO] Merging list of VCF "
			if ! $BCFTOOLS merge $INPUT_LIST --force-samples -o $INPUT_MULTI_VCF -m none 1>>$LOG 2>>$ERR; then
				echo "#[ERROR] Merging failed ";
				exit 0;	
			fi;
		elif (($(echo $INPUT_LIST | wc -w))); then
			if ! $BCFTOOLS view $INPUT_LIST -o $INPUT_MULTI_VCF 1>>$LOG 2>>$ERR; then
				echo "#[ERROR] Input loading failed ";
				exit 0;	
			fi;
		else
			echo "#[ERROR] Input VCF '$INPUT' does not contain any valide VCF files";
			exit 0;
		fi;
		INPUT=$INPUT_MULTI_VCF
		INPUT_ORIGINAL=$INPUT_ORIGINAL_LIST

	else
		echo "#[ERROR] No VCF file '$INPUT'"
		usage; exit;
		
	fi;
	
	# Check if file is empty
	if !(($(grep ^ $INPUT -c))); then
		echo "#[ERROR] Input VCF '$INPUT' is an empty file";
		exit 0;
	fi;
	
	
	if ((1)); then

		# HEADER
		##########

		#echo -e "\n####################\n# INPUT PARAMETERS\n####################\n"

		echo "#[INFO] INPUT: $INPUT_ORIGINAL"; 
		echo "#[INFO] OUTPUT: $OUTPUT"; 
		echo "#[INFO] MULTITHREADING: $MULTITHREADING ($THREADS threads)"; 
		(($VERBOSE)) && echo "#[INFO] LOG: $LOG";
		(($VERBOSE)) && echo "#[INFO] ERR: $ERR";

		# TMP folder and new input file
		TMP_SORT_PREPARATION=$TMP_SORT$RANDOM
		mkdir -p $TMP_SORT_PREPARATION
		INPUT_TMP=$TMP_FOLDER/TMP_$(basename $INPUT | sed "s/.gz$//gi")
	

		# Normalization
		if [ "$NORM" != "" ] && [ -e "$NORM" ]; then
			INPUT_TMP_NORM=$INPUT.norm.vcf
			$BCFTOOLS norm -m- -f $NORM $INPUT > $INPUT_TMP_NORM 2>>$ERR
			cp $INPUT_TMP_NORM $INPUT
			rm $INPUT_TMP_NORM 
		fi;

		#(($VERBOSE)) && echo "#[INFO] TMP_SORT_PREPARATION=$TMP_SORT_PREPARATION"
		(($VERBOSE)) && echo "#[INFO] INPUT_TMP: $INPUT"

		# NB VARIANTS
		#NB_VARIANT=$(grep ^# -cv $INPUT);
		NB_VARIANT=$($BCFTOOLS view -H $INPUT | grep ^# -cv);
		echo "#[INFO] NB VARIANT: $NB_VARIANT";

	fi;

	
	
	
	
	VCF_VALIDATION=$TMP_FOLDER/VCF_VALIDATION
	>$VCF_VALIDATION.tmp
	>$VCF_VALIDATION
	#echo "$BCFTOOLS view $INPUT"
	TMP_SORT_VALIDATION=$TMP_SORT$RANDOM
	mkdir -p $TMP_SORT_VALIDATION
	$BCFTOOLS view $INPUT 2>>$VCF_VALIDATION.tmp | $BCFTOOLS sort -T $TMP_SORT_VALIDATION 1>/dev/null 2>>$VCF_VALIDATION.tmp
	rm -rf $TMP_SORT_VALIDATION
	#/home1/TOOLS/tools/bcftools/1.3.1/bin/bcftools view $INPUT 1>/dev/null 2>$VCF_VALIDATION
	grep -v -e "^Writing to " -e "^Merging .* temporary files$" -e "^Cleaning$" -e "^Done$" $VCF_VALIDATION.tmp > $VCF_VALIDATION
	rm $VCF_VALIDATION.tmp

	if (($(grep ^Failed $VCF_VALIDATION -c))) || (($(grep ^Aborded $VCF_VALIDATION -c))); then
		echo "#[ERROR] Input VCF $INPUT failed";
		cat $VCF_VALIDATION;
		exit 0;
	elif (($(grep "^\[E::" $VCF_VALIDATION -c))); then
		echo "#[ERROR] Some Errors on INPUT VCF $INPUT:";
		cat $VCF_VALIDATION;
		exit 0;
	elif (($(cat $VCF_VALIDATION | wc -l))); then
		echo "#[WARNING] Some Errors/Warnings on INPUT VCF $INPUT:";
		cat $VCF_VALIDATION;
	else
		echo "#[INFO] INPUT VCF validated";
	fi;

fi;


# CLEANING PARAM
#PARAM=$(echo "$PARAM " | sed "s/--threads=[^ |$]*//gi" | sed "s/--multithreading[^ |$]*//gi")
PARAM=$(echo "$PARAM " | sed "s/--threads=[^ |$]*//gi" | sed "s/--multithreading[^ |$]*//gi" | sed "s/--env[^ |$]*//gi" | sed "s/--compress[^ |$]*//gi" | sed "s/--norm[^ |$]*//gi" | sed "s/--annotate[^ |$]*//gi")

#echo $PARAM; exit 0;


# MK
MK=$TMP_FOLDER/mk_annotation
MK_CALCULATION=$TMP_FOLDER/mk_calculation
MK_PRIORITIZATION=$TMP_FOLDER/mk_prioritization
MK_TRANSLATION=$TMP_FOLDER/mk_translation


# OUTPUT FILES
OUTPUT_ANNOTATION=$TMP_FOLDER/output_annotation
OUTPUT_CALCULATION=$TMP_FOLDER/output_calculation
OUTPUT_PRIORITIZATION=$TMP_FOLDER/output_prioritization
OUTPUT_TRANSLATION=$TMP_FOLDER/output_translation





# Clean INPUT
#INPUT_CLEAN=$TMP_FOLDER/input.vcf
#grep "^##" $INPUT | grep "^##.*=.*" > $INPUT_CLEAN
#grep "^##" $INPUT | grep "^##.*=.*" -v >> $INPUT_CLEAN
#grep "^##" $INPUT -v >> $INPUT_CLEAN
#INPUT=$INPUT_CLEAN

####################
# ANNOTATION
####################

if (($VERBOSE)); then
	tail -f $LOG &
	TAIL_LOG_PID=$!
	#echo "TAIL_LOG_PID=$TAIL_LOG_PID"
fi;


#head $INPUT

if (($NB_VARIANT)); then


	# ANNOTATE with BCFTOOLS
	# example : ./HOWARD.sh --vcf=validation/simple.vcf --annotate="validation/test2.vcf.gz:INFO/Symbol|validation/test2.vcf.gz:INFO/RefSeq"
	# ./HOWARD.sh --vcf=validation/simple.vcf --annotate="validation/test2.vcf.gz:INFO/Symbol;validation/test1.vcf.gz:INFO/RefSeq;validation/test2.vcf:INFO/Ensembl;validation/test2.vcf"  --annotation=Symbol,RefSeq --verbose --force --multithreading --debug
	
	if [ "$ANNOTATE" != "" ] && (($NB_VARIANT)); then
		echo -e "\n####################\n# ANNOTATE with BCFTOOLS\n####################\n"
		$BGZIP $INPUT -c > $INPUT.gz
		$TABIX $INPUT.gz

		# list of annotation "VCF/TAG"		
		for ANNOT in $(echo $ANNOTATE | tr ";" " "); do
		
			#echo $ANNOT
			ANNOT_FILE=$(echo $ANNOT | awk -F":" '{print $1}')
			ANNOT_TAG=$(echo $ANNOT | awk -F":" '{print $2}')
			#echo "-a $ANNOT_FILE -c $ANNOT_TAG"
			
			if [ -e "$ANNOT_FILE" ]; then
			
				# Prepare VCF annot
				TMP_SORT_PREPARATION=$TMP_SORT$RANDOM
				ANNOT_FILE_TMP=$TMP_FOLDER/ANNOT_FILE$RANDOM.vcf
				mkdir -p $TMP_SORT_PREPARATION
				$BCFTOOLS sort -T $TMP_SORT_PREPARATION $ANNOT_FILE 2>>$ERR | $BCFTOOLS norm -m- >$ANNOT_FILE_TMP 2>>$ERR
				$BGZIP $ANNOT_FILE_TMP
				$TABIX $ANNOT_FILE_TMP.gz
				rm -rf $TMP_SORT_PREPARATION
				
				if [ "$ANNOT_TAG" == "" ]; then
					#ANNOT_TAG="+INFO"
					ANNOT_TAG="INFO"
				fi;
				
				if [ -e "$ANNOT_FILE_TMP.gz" ] && [ -e "$ANNOT_FILE_TMP.gz.tbi" ]; then
					
					echo "#[INFO] BCFTOOLS Annotation with VCF file '$ANNOT_FILE' and TAG '$ANNOT_TAG'"
					
					INPUT_TMP=$TMP_FOLDER/INPUT_ANNOT_FILE$RANDOM.vcf.gz
					$BGZIP -c $INPUT > $INPUT_TMP
					$TABIX $INPUT_TMP
			
					INPUT_ANNOTATE_TMP=$TMP_FOLDER/INPUT_ANNOTATE_TMP$RANDOM.vcf
			
					$BCFTOOLS annotate -a $ANNOT_FILE_TMP.gz -c $ANNOT_TAG $INPUT_TMP > $INPUT_ANNOTATE_TMP #2>>$ERR
				
					if [ -e "$INPUT_ANNOTATE_TMP" ]; then
						mv -f $INPUT_ANNOTATE_TMP $INPUT 
					fi;
				
				fi;
			
			else
				echo "#[WARNING] Annotation file '$ANNOT' not found"
			fi;
			
		done;
		
		#exit 0;
		
	fi;



	# ANNOTATION
	#if [ ! -z $ANNOTATION ] && [ "$ANNOTATION$SNPEFF_STATS" != "" ] && (($NB_VARIANT)); then
	if [ "$ANNOTATION$SNPEFF_STATS" != "" ] && (($NB_VARIANT)); then
		echo -e "\n####################\n# ANNOTATION\n####################\n"
		#$SCRIPT_DIR/VCFannotation.pl --header

		#echo $SCRIPT_DIR/VCFannotation.pl  $PARAM --input=$INPUT --show_annotations
		#$SCRIPT_DIR/VCFannotation.pl  $PARAM --input=$INPUT --show_annotations; exit 0;
		CMD_SHOW_ANNOTATION="$SCRIPT_DIR/VCFannotation.pl  $PARAM --input=$INPUT --show_annotations"
		#echo $CMD_SHOW_ANNOTATION
		#eval $CMD_SHOW_ANNOTATION
		#$SCRIPT_DIR/VCFannotation.pl  $PARAM --input=$INPUT --show_annotations
		#ANNOTATIONS=$($SCRIPT_DIR/VCFannotation.pl  $PARAM --input=$INPUT --show_annotations | grep "# ANNOTATIONS: " | cut -d: -f2);
		ANNOTATIONS=$(eval $CMD_SHOW_ANNOTATION | grep "# ANNOTATIONS: " | cut -d: -f2);
		NB_ANNOTATIONS=$(echo $ANNOTATIONS | wc -w);
		#echo "NB_VARIANT=$NB_VARIANT ANNOTATIONS_SPLIT=$ANNOTATIONS_SPLIT ANNOTATIONS=$ANNOTATIONS"; exit 0;

		#if (($MULTITHREADING)) && [ $NB_ANNOTATIONS -gt 1 ]; then
		#if (($MULTITHREADING)) || [ $THREADS -gt 0 ]; then
		if (($THREADS)); then

			if ((1)); then
				# ECHO / VERBOSE / DEBUG
				echo "#[INFO] Multithreading ($THREADS threads)"

			
				# NEEDED
				#if [ ! -z $BGZIP ] && [ -e $BGZIP ]; then
				if [ ! -z $BGZIP ] && (( $(command -v $BGZIP | grep ^ -c) )); then

					echo "#[INFO] BGZIP=$BGZIP";
				else
					echo "#[ERROR] BGZIP '$BGZIP' needed. Configure ENV file";
					exit 0;
				fi;
				#if [ ! -z $VCFTOOLS ] && [ -e $VCFTOOLS ]; then
				#	echo "#[INFO] VCFTOOLS=$VCFTOOLS";
				#else
				#	echo "#[ERROR] VCFTOOLS '$VCFTOOLS' needed. Configure ENV file";
				#	exit 0;
				#fi;
				if [ ! -z $BCFTOOLS ] && [ -e $BCFTOOLS ]; then
					echo "#[INFO] BCFTOOLS=$BCFTOOLS";
				else
					echo "#[ERROR] BCFTOOLS '$BCFTOOLS' needed. Configure ENV file";
					exit 0;
				fi;
				#if [ ! -z $TABIX ] && [ -e $TABIX ]; then
				if [ ! -z $TABIX ] && (( $(command -v $TABIX | grep ^ -c) )); then
					echo "#[INFO] TABIX=$TABIX";
				else
					echo "#[ERROR] TABIX '$TABIX' needed. Configure ENV file";
					exit 0;
				fi;
			fi;

			# Annotation
			echo "#[INFO] Annotations splitted into $NB_ANNOTATIONS analyses"
			(($VERBOSE)) && echo "#[INFO] ANNOTATION=$ANNOTATION" 1>>$LOG 2>$ERR
			(($VERBOSE)) && echo "#[INFO] ANNOTATIONS=$ANNOTATIONS" 1>>$LOG 2>$ERR
			
			#if [ $NB_VARIANT -gt $ANNOTATIONS_SPLIT ]; then
				# SPLIT VCF on ANNOTATIONS_SPLIT variants
				# HEAD
				#grep "^#" $INPUT > $TMP_FOLDER/input.header
				$BCFTOOLS view -h $INPUT > $TMP_FOLDER/input.header
				# VARIANTS
				#grep -v "^#" $INPUT > $TMP_FOLDER/input.variants
				$BCFTOOLS view -H $INPUT > $TMP_FOLDER/input.variants
				# SPLIT
				#echo "SPLIT"
				#grep -v "^#" $INPUT | split - -a $SPLIT_SUFFIT_LENGTH -l $ANNOTATIONS_SPLIT $TMP_FOLDER/input.variants.splitted.
				$BCFTOOLS view -H $INPUT | split - -a $SPLIT_SUFFIT_LENGTH -l $ANNOTATIONS_SPLIT $TMP_FOLDER/input.variants.splitted.
				#echo "SPLIT END"
				VCF_SPLITTED_LIST=""
				#for splited_variants in $TMP_FOLDER/input.variants.splitted.*; do
				for splited_variants in $(ls "$TMP_FOLDER/"input.variants.splitted.* | grep "$TMP_FOLDER/"'input.variants.splitted.[a-z]*$' ); do
					cat $TMP_FOLDER/input.header $splited_variants > $splited_variants.vcf
					VCF_SPLITTED_LIST=$VCF_SPLITTED_LIST" $splited_variants.vcf"
				done
				#echo "SPLITTED: $VCF_SPLITTED_LIST"; exit 0;
				#cat $VCF_SPLITTED_LIST
				#exit 0;
				echo "#[INFO] Input VCF splitted into "$(echo $VCF_SPLITTED_LIST | wc -w)" files";

			#fi;
			#exit 0;


			# PARAM ANNOTATION
			#PARAM_ANNOTATION=$(echo "$PARAM " | sed "s/--snpeff_stats=[^ |$]*//gi");
			#PARAM_ANNOTATION=$(echo "$PARAM " | sed "s/--snpeff_stats=[ |$]*//gi" | sed "s/--snpeff_hgvs[ |$]*//gi"  | sed "s/--snpeff [ |$]*//gi");
			#PARAM_ANNOTATION=$(echo "$PARAM " | sed "s/--snpeff_stats=[^ |$]*//gi" | sed "s/--snpeff_hgvs[^ |$]*//gi"  | sed "s/--snpeff[^ |$]*//gi");
			PARAM_ANNOTATION=$(echo "$PARAM " | sed "s/--snpeff_stats=[^ |$]*//gi");
			#echo $PARAM_ANNOTATION; exit 0;

			# FIND SAMPLE
			if [ $($BCFTOOLS view -h $INPUT | grep "#CHROM" | wc -w) -gt 9 ]; then
			#if [ $($BCFTOOLS view -h $INPUT | grep "#CHROM" | wc -w) -gt 9 ]; then
				SAMPLES=$($BCFTOOLS view -h $INPUT | grep "#CHROM" | cut -f10-$(grep "#CHROM" $INPUT | wc -w) --output-delimiter=,)
				#SAMPLES=$($BCFTOOLS view -h $INPUT | grep "#CHROM" | cut -f10-$($BCFTOOLS view -h $INPUT | grep "#CHROM" | wc -w) --output-delimiter=,)
			else
				SAMPLES=""
			fi;

			SAMPLES_COUNT=$(echo $SAMPLES | tr "," " " | wc -w)
			VCF_HEADER_CHROM_LENGTH=$(($SAMPLES_COUNT+9))
			#(($VERBOSE)) && echo "#[INFO] SAMPLES (#SAMPLES_COUNT)=$SAMPLES" 1>>$LOG 2>$ERR
			(($VERBOSE)) && echo "#[INFO] #SAMPLES=$SAMPLES_COUNT" 1>>$LOG 2>$ERR
			(($VERBOSE)) && echo "#[INFO] SAMPLES=$SAMPLES" 1>>$LOG 2>$ERR
			
			# IF SAMPLE(S)
			if [ "$SAMPLES " != "" ]; then

				# ALL
				echo -e "all: $SNPEFF_STATS $OUTPUT_ANNOTATION
				" >>$MK;

				# SNPEFF STATS
				if [ ! -z $SNPEFF_STATS ]; then
					(($VERBOSE)) && echo "# snpEff Stats '$SNPEFF_STATS'";
					echo "$SNPEFF_STATS: $INPUT
						$SCRIPT_DIR/VCFannotation.pl $PARAM --output=/dev/null  --annotation=snpeff --snpeff_stats=$SNPEFF_STATS  1>>$LOG 2>>$ERR
					" >>$MK;
				fi;


				if ((1)); then
				
					# LOOP on split files
					VCFGZ_SPLITTED_ANNOTATED_LIST=""
					VCFGZ_SPLITTED_ANNOTATED_FILE=$TMP_FOLDER/splited.annotated.list
					for VCF_SPLITTED in $VCF_SPLITTED_LIST; do
						VCFGZ_SPLITTED=$VCF_SPLITTED.vcf.gz
						VCFGZ_SPLITTED_ANNOTATED=$VCF_SPLITTED.ann.vcf.gz
						VCFGZ_SPLITTED_ANN_LIST=""
						VCFGZ_SPLITTED_ANN_FILE=$VCF_SPLITTED.ann.list
						#(($VERBOSE)) && echo "#   VCFGZ_SPLITTED=$VCFGZ_SPLITTED" 1>>$LOG 2>$ERR
						
						for ANN in $(echo $ANNOTATIONS | tr "," " "); do
						
							VCFGZ_SPLITTED_ANN=$VCF_SPLITTED.$ANN.vcf.gz
							
							#(($VERBOSE)) && echo "#      VCFGZ_SPLITTED_ANN=$VCFGZ_SPLITTED_ANN" 1>>$LOG 2>$ERR
							
							echo "$VCFGZ_SPLITTED_ANN: $VCF_SPLITTED
								$SCRIPT_DIR/VCFannotation.pl $PARAM_ANNOTATION --annotation=$ANN --output=\$@.tmp --input=$VCF_SPLITTED  1>>$LOG 2>>$ERR
								# Prevent PL errorin number of fields in BCFTOOLS merge
								cat \$@.tmp | sed s/ID=PL,Number=G/ID=PL,Number=./gi > \$@.tmp2
								$BGZIP -c \$@.tmp2 > \$@ 2>>$ERR
								-rm -f \$@.tmp \$@.tmp2
								$TABIX \$@ 1>>$LOG 2>>$ERR
							" >>$MK;
							
							VCFGZ_SPLITTED_ANN_LIST=$VCFGZ_SPLITTED_ANN_LIST" $VCFGZ_SPLITTED_ANN"
							#echo "$VCFGZ_SPLITTED_AN" >> $VCFGZ_SPLITTED_ANN_FILE
							
						done;
						echo $VCFGZ_SPLITTED_ANN_LIST | tr " " "\n" > $VCFGZ_SPLITTED_ANN_FILE
						
						#(($VERBOSE)) && echo "#      VCFGZ_SPLITTED_ANN_LIST=$VCFGZ_SPLITTED_ANN_LIST" 1>>$LOG 2>$ERR
						
						# MERGE ANNOTATION
						echo "$VCFGZ_SPLITTED_ANNOTATED: $VCFGZ_SPLITTED_ANN_LIST
							$BCFTOOLS merge -l $VCFGZ_SPLITTED_ANN_FILE  --force-samples | cut -f1-$VCF_HEADER_CHROM_LENGTH > \$@.tmp;
							$BGZIP -f \$@.tmp -c > \$@ ;
							$TABIX -f \$@;
							-rm -f $\$@.tmp;
						" >>$MK;
						
						# OLD version of merge using annotate. Error when Number=R !!!
						if ((0)); then
						echo "$VCFGZ_SPLITTED_ANNOTATED: $VCFGZ_SPLITTED_ANN_LIST
							
							$BGZIP $VCF_SPLITTED -c > \$@.merge.vcf.gz ;
							$TABIX -f \$@.merge.vcf.gz;
							for A in $VCFGZ_SPLITTED_ANN_LIST; do \
								echo '# Add annotation '\$\$A ; \
								$BCFTOOLS annotate -a \$\$A -c INFO \$@.merge.vcf.gz -o \$@.merge.vcf --no-version ; \
								$BGZIP -f \$@.merge.vcf ; \
								$TABIX -f \$@.merge.vcf.gz; \
								rm -f $A; \
							done;
							cp \$@.merge.vcf.gz \$@;
							cp \$@.merge.vcf.gz.tbi \$@.tbi;
							-rm -f $\$@.merge.vcf*;
						" >>$MK;
						fi;
						
						
						VCFGZ_SPLITTED_ANNOTATED_LIST=$VCFGZ_SPLITTED_ANNOTATED_LIST" "$VCFGZ_SPLITTED_ANNOTATED
						echo $VCFGZ_SPLITTED_ANNOTATED >> $VCFGZ_SPLITTED_ANNOTATED_FILE
						
					done;
				
					
					# construct file list (when processed)
					echo "$OUTPUT_ANNOTATION.list: $VCFGZ_SPLITTED_ANNOTATED_FILE $VCFGZ_SPLITTED_ANNOTATED_LIST
						cp $VCFGZ_SPLITTED_ANNOTATED_FILE $OUTPUT_ANNOTATION.list
					" >>$MK;
					
					# Construct Header	
					# if [[ $(echo $VCFGZ_SPLITTED_ANNOTATED_LIST | wc -w) -gt 1 ]]; then \
					echo "$OUTPUT_ANNOTATION.header: $OUTPUT_ANNOTATION.list
						if [[ \$\$(grep ^# -cv $OUTPUT_ANNOTATION.list) -gt 1 ]]; then \
							$BCFTOOLS merge -l $OUTPUT_ANNOTATION.list --force-samples --print-header --no-version | cut -f1-$VCF_HEADER_CHROM_LENGTH > \$@ ; \
						else \
							$BCFTOOLS view -h \$\$(cat $OUTPUT_ANNOTATION.list) --no-version > \$@ ; \
						fi ;
					" >>$MK;
					
					# Merge splitted files
					echo "$OUTPUT_ANNOTATION: $OUTPUT_ANNOTATION.list $OUTPUT_ANNOTATION.header
							mkdir -p $TMP_SORT\$@
							$BCFTOOLS concat -f $OUTPUT_ANNOTATION.list -a --no-version | $BCFTOOLS reheader -h $OUTPUT_ANNOTATION.header --threads $THREADS | $BCFTOOLS sort -T $TMP_SORT\$@ > \$@ 2>>$ERR
							
						" >>$MK;
					
				fi;
				
				#cat $MK;
				#exit 0;

				# MAKE
				if ! make -s -B -j $THREADS -f $MK 1>>$LOG 2>>$ERR; then echo "#[FAILED] ANNOTATION Multithreading FAILED!!! see $MK file"; exit 1; fi;

				# ECHO / VERBOSE / DEBUG
				(($DEBUG)) && echo $MK && cat $MK;
				(($DEBUG)) && tail -n 10 $OUTPUT_ANNOTATION.merge.vcf
				(($DEBUG)) && tail -n 10 $OUTPUT_ANNOTATION
				

				# CLEANING
				rm -f $MK

			else
				echo "#[ERROR] no Sample name found in input VCF file"; exit 1;
			fi;
		else
			#echo "$SCRIPT_DIR/VCFannotation.pl $PARAM --input=$INPUT --output=$OUTPUT_ANNOTATION"
			#$SCRIPT_DIR/VCFannotation.pl $PARAM --input=$INPUT --output=$OUTPUT_ANNOTATION 1>>$LOG 2>$ERR
			echo "#[ERROR] threads parameter error"; exit 1;
		fi;

	else

		cp $INPUT $OUTPUT_ANNOTATION

	fi;

	# ECHO / VERBOSE / DEBUG
	(($DEBUG)) && tail -n 10 $OUTPUT_ANNOTATION


	####################
	# CALCULATION
	####################

	if [ ! -z $CALCULATION ] && [ "$CALCULATION" != "" ] && (($NB_VARIANT)) && ((1)); then
		echo -e "\n####################\n# CALCULATION\n####################\n"
		#$SCRIPT_DIR/VCFcalculation.pl --header


		if ((1)); then


			#if (($MULTITHREADING)); then
			if (($THREADS)); then

				if ((1)); then
				
					# ECHO / VERBOSE / DEBUG
					echo "#[INFO] Multithreading ($THREADS threads)"

					# NEEDED
					if [ ! -z $BGZIP ] && [ -e $BZZIP ]; then
						echo "#[INFO] BGZIP=$BGZIP";
					else
						echo "#[ERROR] BGZIP '$BGZIP' needed. Configure ENV file";
						exit 0;
					fi;
					#if [ ! -z $VCFTOOLS ] && [ -e $VCFTOOLS ]; then
					#	echo "#[INFO] VCFTOOLS=$VCFTOOLS";
					#else
					#	echo "#[ERROR] VCFTOOLS '$VCFTOOLS' needed. Configure ENV file";
					#	exit 0;
					#fi;
					if [ ! -z $BCFTOOLS ] && [ -e $BCFTOOLS ]; then
						echo "#[INFO] BCFTOOLS=$BCFTOOLS";
					else
						echo "#[ERROR] BCFTOOLS '$BCFTOOLS' needed. Configure ENV file";
						exit 0;
					fi;
					
				fi;

				#echo $NB_VARIANT
				#echo $CALCULATION_SPLIT
				VCF_SPLITTED_CALCULATION_LIST=""
				if [ $NB_VARIANT -gt $CALCULATION_SPLIT ]; then
					# SPLIT VCF on ANNOTATION_SPLIT variants
					# HEAD
					#grep "^#" $OUTPUT_ANNOTATION > $TMP_FOLDER/input_annotation.header
					$BCFTOOLS view -h $OUTPUT_ANNOTATION > $TMP_FOLDER/input_annotation.header
					# VARIANTS
					#grep -v "^#" $OUTPUT_ANNOTATION > $TMP_FOLDER/input_annotation.variants
					$BCFTOOLS view -H $OUTPUT_ANNOTATION > $TMP_FOLDER/input_annotation.variants
					# SPLIT
					#echo "SPLIT"
					#grep -v "^#" $OUTPUT_ANNOTATION | split - -a $SPLIT_SUFFIT_LENGTH -l $CALCULATION_SPLIT $TMP_FOLDER/input_annotation.variants.splitted.
					$BCFTOOLS view -H $OUTPUT_ANNOTATION | split - -a $SPLIT_SUFFIT_LENGTH -l $CALCULATION_SPLIT $TMP_FOLDER/input_annotation.variants.splitted.
					#echo "SPLIT END"
					for splited_variants in $TMP_FOLDER/input_annotation.variants.splitted.*; do
						cat $TMP_FOLDER/input_annotation.header $splited_variants > $splited_variants.vcf
						VCF_SPLITTED_CALCULATION_LIST=$VCF_SPLITTED_CALCULATION_LIST" $splited_variants.vcf"
					done
					#echo "SPLITTED: $VCF_SPLITTED_LIST"
					#cat $VCF_SPLITTED_LIST
					#exit 0;
					echo "#[INFO] Input VCF splitted into "$(echo $VCF_SPLITTED_CALCULATION_LIST | wc -w)" files";
					#echo $VCF_SPLITTED_CALCULATION_LIST

				else
					cp $OUTPUT_ANNOTATION $TMP_FOLDER/input_annotation.variants.splitted.vcf
					VCF_SPLITTED_CALCULATION_LIST=$TMP_FOLDER/input_annotation.variants.splitted.vcf
				fi;
				#exit 0;


				# PARAM ANNOTATION
				PARAM_CALCULATION=$(echo "$PARAM " | sed "s/--snpeff_stats=[^ |$]*//gi");

				# FIND SAMPLE
				if [ $($BCFTOOLS view -h $INPUT | grep "#CHROM" | wc -w) -gt 9 ]; then
				#if [ $($BCFTOOLS view -h $INPUT | grep "#CHROM" | wc -w) -gt 9 ]; then
					SAMPLES=$($BCFTOOLS view -h $INPUT | grep "#CHROM" | cut -f10-$(grep "#CHROM" $INPUT | wc -w) --output-delimiter=,)
					#SAMPLES=$($BCFTOOLS view -h $INPUT | grep "#CHROM" | cut -f10-$($BCFTOOLS view -h $INPUT | grep "#CHROM" | wc -w) --output-delimiter=,)
				else
					SAMPLES=""
				fi;

				SAMPLES_COUNT=$(echo $SAMPLES | tr "," " " | wc -w)
				VCF_HEADER_CHROM_LENGTH=$(($SAMPLES_COUNT+9))
				#(($VERBOSE)) && echo "#[INFO] SAMPLES (#SAMPLES_COUNT)=$SAMPLES" 1>>$LOG 2>$ERR
				(($VERBOSE)) && echo "#[INFO] #SAMPLES=$SAMPLES_COUNT" 1>>$LOG 2>$ERR
				(($VERBOSE)) && echo "#[INFO] SAMPLES=$SAMPLES" 1>>$LOG 2>$ERR

				# IF SAMPLE(S)
				if [ "$SAMPLES " != "" ]; then

					# ALL
					echo -e "all: $OUTPUT_CALCULATION
					" >>$MK_CALCULATION;

					for VCF_SPLITTED in $VCF_SPLITTED_CALCULATION_LIST; do
						VCFGZ_SPLITTED_CALCULATION=$VCF_SPLITTED.calculation.vcf.gz
						#echo "VCFGZ_SPLITTED_CALCULATION=$VCFGZ_SPLITTED_CALCULATION"
						echo "$VCFGZ_SPLITTED_CALCULATION: $VCF_SPLITTED
						$SCRIPT_DIR/VCFcalculation.pl $PARAM_CALCULATION --output=\$@.tmp --input=$VCF_SPLITTED  1>>$LOG 2>>$ERR
						mkdir -p $TMP_SORT\$@
						$BCFTOOLS sort -T $TMP_SORT\$@ \$@.tmp | $BGZIP -c > \$@ 2>>$ERR
						rm -rf $TMP_SORT\$@
						-rm -f \$@.tmp \$@.tmp2
						$TABIX \$@ 1>>$LOG 2>>$ERR
						" >>$MK_CALCULATION;
						VCFGZ_SPLITTED_CALCULATION_LIST=$VCFGZ_SPLITTED_CALCULATION_LIST" $VCFGZ_SPLITTED_CALCULATION"
					done;
					#echo $VCFGZ_SPLITTED_CALCULATION_LIST

					echo "$OUTPUT_CALCULATION: $VCFGZ_SPLITTED_CALCULATION_LIST
							@$BCFTOOLS concat -a $VCFGZ_SPLITTED_CALCULATION_LIST --threads $THREADS --no-version  >> \$@
						" >>$MK_CALCULATION;


					if ! make -B -j $THREADS -f $MK_CALCULATION 1>>$LOG 2>>$ERR; then echo "#[FAILED] CALCULATION Multithreading FAILED!!! see $MK_CALCULATION file"; exit 1; fi;

					# ECHO / VERBOSE / DEBUG
					(($DEBUG)) && echo $MK_CALCULATION && cat $MK_CALCULATION;
					(($DEBUG)) && tail -n 10 $OUTPUT_CALCULATION.merge.vcf
					(($DEBUG)) && tail -n 10 $OUTPUT_CALCULATION

					# CLEANING
					rm -f $MK_CALCULATION

					#echo " $SCRIPT_DIR/VCFcalculation.pl $PARAM --input=$OUTPUT_ANNOTATION --output=$OUTPUT_CALCULATION"
				fi;

			else

				#echo " $SCRIPT_DIR/VCFcalculation.pl $PARAM --input=$OUTPUT_ANNOTATION --output=$OUTPUT_CALCULATION "
				#echo $LOG
				#echo $ERR
				if ! $SCRIPT_DIR/VCFcalculation.pl $PARAM --input=$OUTPUT_ANNOTATION --output=$OUTPUT_CALCULATION 1>>$LOG 2>$ERR; then echo "#[FAILED] CALCULATION FAILED!!! "; exit 1; fi;


			fi;

			#echo "DEBUG"
			#exit 0;

		fi;


	else

		cp $OUTPUT_ANNOTATION $OUTPUT_CALCULATION
	fi;

	# clear intermediate files
	rm -f $OUTPUT_ANNOTATION

	# ECHO / VERBOSE / DEBUG
	(($DEBUG)) && tail -n 10 $OUTPUT_CALCULATION




	####################
	# PRIORITIZATION
	####################

	if [ ! -z $FILTER ] && [ "$FILTER" != "" ] && (($NB_VARIANT)); then
		echo -e "\n####################\n# PRIORITIZATION\n####################\n"
		#$SCRIPT_DIR/VCFprioritization.pl --header


		if ((1)); then


			#if (($MULTITHREADING)); then
			if (($THREADS)); then


				if ((1)); then
					# ECHO / VERBOSE / DEBUG
					echo "#[INFO] Multithreading ($THREADS threads)"

					# NEEDED
					if [ ! -z $BGZIP ] && [ -e $BZZIP ]; then
						echo "#[INFO] BGZIP=$BGZIP";
					else
						echo "#[ERROR] BGZIP '$BGZIP' needed. Configure ENV file";
						exit 0;
					fi;
					#if [ ! -z $VCFTOOLS ] && [ -e $VCFTOOLS ]; then
					#	echo "#[INFO] VCFTOOLS=$VCFTOOLS";
					#else
					#	echo "#[ERROR] VCFTOOLS '$VCFTOOLS' needed. Configure ENV file";
					#	exit 0;
					#fi;
					if [ ! -z $BCFTOOLS ] && [ -e $BCFTOOLS ]; then
						echo "#[INFO] BCFTOOLS=$BCFTOOLS";
					else
						echo "#[ERROR] BCFTOOLS '$BCFTOOLS' needed. Configure ENV file";
						exit 0;
					fi;
				fi;

				#echo $NB_VARIANT
				#echo $PRIORITIZATION_SPLIT
				VCF_SPLITTED_PRIO_LIST=""
				if [ $NB_VARIANT -gt $PRIORITIZATION_SPLIT ]; then
					# SPLIT VCF on CALCULATION_SPLIT variants
					# HEAD
					#grep "^#" $OUTPUT_CALCULATION > $TMP_FOLDER/input_calculation.header
					$BCFTOOLS view -h $OUTPUT_CALCULATION > $TMP_FOLDER/input_calculation.header
					# VARIANTS
					#grep -v "^#" $OUTPUT_CALCULATION > $TMP_FOLDER/input_calculation.variants
					$BCFTOOLS view -H $OUTPUT_CALCULATION > $TMP_FOLDER/input_calculation.variants
					# SPLIT
					#echo "SPLIT"
					#grep -v "^#" $OUTPUT_CALCULATION | split - -a $SPLIT_SUFFIT_LENGTH -l $PRIORITIZATION_SPLIT $TMP_FOLDER/input_calculation.variants.splitted.
					$BCFTOOLS view -H $OUTPUT_CALCULATION | split - -a $SPLIT_SUFFIT_LENGTH -l $PRIORITIZATION_SPLIT $TMP_FOLDER/input_calculation.variants.splitted.
					#echo "SPLIT END"
					for splited_variants in $TMP_FOLDER/input_calculation.variants.splitted.*; do
						cat $TMP_FOLDER/input_calculation.header $splited_variants > $splited_variants.vcf
						VCF_SPLITTED_PRIO_LIST=$VCF_SPLITTED_PRIO_LIST" $splited_variants.vcf"
					done
					#echo "SPLITTED: $VCF_SPLITTED_LIST"
					#cat $VCF_SPLITTED_LIST
					#exit 0;
					echo "#[INFO] Input VCF splitted into "$(echo $VCF_SPLITTED_PRIO_LIST | wc -w)" files";
					#echo $VCF_SPLITTED_PRIO_LIST

				else
					cp $OUTPUT_CALCULATION $TMP_FOLDER/input_calculation.variants.splitted.vcf
					VCF_SPLITTED_PRIO_LIST=$TMP_FOLDER/input_calculation.variants.splitted.vcf
				fi;
				#exit 0;


				# PARAM CALCULATION
				PARAM_PRIORITIZATION=$(echo "$PARAM " | sed "s/--snpeff_stats=[^ |$]*//gi");

				# FIND SAMPLE
				if [ $($BCFTOOLS view -h $INPUT | grep "#CHROM" | wc -w) -gt 9 ]; then
				#if [ $($BCFTOOLS view -h $INPUT | grep "#CHROM" | wc -w) -gt 9 ]; then
					SAMPLES=$($BCFTOOLS view -h $INPUT | grep "#CHROM" | cut -f10-$(grep "#CHROM" $INPUT | wc -w) --output-delimiter=,)
					#SAMPLES=$($BCFTOOLS view -h $INPUT | grep "#CHROM" | cut -f10-$($BCFTOOLS view -h $INPUT | grep "#CHROM" | wc -w) --output-delimiter=,)
				else
					SAMPLES=""
				fi;

				SAMPLES_COUNT=$(echo $SAMPLES | tr "," " " | wc -w)
				VCF_HEADER_CHROM_LENGTH=$(($SAMPLES_COUNT+9))
				#(($VERBOSE)) && echo "#[INFO] SAMPLES (#SAMPLES_COUNT)=$SAMPLES" 1>>$LOG 2>$ERR
				(($VERBOSE)) && echo "#[INFO] #SAMPLES=$SAMPLES_COUNT" 1>>$LOG 2>$ERR
				(($VERBOSE)) && echo "#[INFO] SAMPLES=$SAMPLES" 1>>$LOG 2>$ERR

				# IF SAMPLE(S)
				if [ "$SAMPLES " != "" ]; then

					#echo "#[INFO] PRIORITIZE..."

					# ALL
					echo -e "all: $OUTPUT_PRIORITIZATION
					" >>$MK_PRIORITIZATION;

					for VCF_SPLITTED in $VCF_SPLITTED_PRIO_LIST; do
						VCFGZ_SPLITTED_PRIO=$VCF_SPLITTED.prio.vcf.gz
						#echo "VCFGZ_SPLITTED_PRIO=$VCFGZ_SPLITTED_PRIO"
						echo "$VCFGZ_SPLITTED_PRIO: $VCF_SPLITTED
						$SCRIPT_DIR/VCFprioritization.pl $PARAM_PRIORITIZATION --output=\$@.tmp --input=$VCF_SPLITTED  1>>$LOG 2>>$ERR
						mkdir -p $TMP_SORT\$@
						$BCFTOOLS sort -T $TMP_SORT\$@ \$@.tmp | $BGZIP -c > \$@ 2>>$ERR
						rm -rf $TMP_SORT\$@
						-rm -f \$@.tmp \$@.tmp2
						$TABIX \$@ 1>>$LOG 2>>$ERR
						" >>$MK_PRIORITIZATION;
						VCFGZ_SPLITTED_PRIO_LIST=$VCFGZ_SPLITTED_PRIO_LIST" $VCFGZ_SPLITTED_PRIO"
					done;
					#echo $VCFGZ_SPLITTED_PRIO_LIST

					echo "$OUTPUT_PRIORITIZATION: $VCFGZ_SPLITTED_PRIO_LIST
							$BCFTOOLS concat $VCFGZ_SPLITTED_PRIO_LIST -a --threads $THREADS --no-version >> \$@
						" >>$MK_PRIORITIZATION;

					#echo $MK_PRIORITIZATION
					#cat $MK_PRIORITIZATION

					if ! make -B -j $THREADS -f $MK_PRIORITIZATION 1>>$LOG 2>>$ERR; then echo "#[FAILED] PRIORITIZATION Multithreading FAILED!!! see $MK_PRIORITIZATION file"; exit 1; fi;

					# ECHO / VERBOSE / DEBUG
					(($DEBUG)) && echo $MK_PRIORITIZATION && cat $MK_PRIORITIZATION;
					(($DEBUG)) && tail -n 10 $OUTPUT_PRIORITIZATION.merge.vcf
					(($DEBUG)) && tail -n 10 $OUTPUT_PRIORITIZATION

					# CLEANING
					rm -f $MK_PRIORITIZATION

					#echo " $SCRIPT_DIR/VCFprioritization.pl $PARAM --input=$OUTPUT_CALCULATION --output=$OUTPUT_PRIORITIZATION"
				fi;

			else

				if ! $SCRIPT_DIR/VCFprioritization.pl $PARAM --input=$OUTPUT_CALCULATION --output=$OUTPUT_PRIORITIZATION 1>>$LOG 2>$ERR; then echo "#[FAILED] PRIORITIZATION FAILED!!! "; exit 1; fi;

			fi;

			#echo "DEBUG"
			#exit 0;

		fi;



		# OLD
		if ((0)); then
			if ! $SCRIPT_DIR/VCFprioritization.pl $PARAM --input=$OUTPUT_CALCULATION --output=$OUTPUT_PRIORITIZATION 1>>$LOG 2>$ERR; then echo "#[FAILED] PRIORITIZATION FAILED!!! "; exit 1; fi;
		fi;
	else

		cp $OUTPUT_CALCULATION $OUTPUT_PRIORITIZATION
	fi;

	# clear intermediate files
	rm -f $OUTPUT_CALCULATION

	# ECHO / VERBOSE / DEBUG
	(($DEBUG)) && tail -n 10 $OUTPUT_PRIORITIZATION


else

	# No variants
	cp $INPUT $OUTPUT_PRIORITIZATION;

fi;


# Cleaning VCF final file
if ((1)); then

	################
	# VCF CLEANING #
	################

	cp $OUTPUT_PRIORITIZATION $OUTPUT_PRIORITIZATION.old

	# Replace PL format as G
	cat $OUTPUT_PRIORITIZATION.old | sed s/ID=PL,Number=\./ID=PL,Number=G/gi > $OUTPUT_PRIORITIZATION

	rm $OUTPUT_PRIORITIZATION.old

fi;


# Add HOWARD HEADER to VCF final file
if ((1)); then

	#################
	# HOWARD HEADER #
	#################

	cp $OUTPUT_PRIORITIZATION $OUTPUT_PRIORITIZATION.old

	grep ^## $OUTPUT_PRIORITIZATION.old > $OUTPUT_PRIORITIZATION

	echo "##"$SCRIPT_NAME"=[$SCRIPT_RELEASE-$SCRIPT_DATE] $SCRIPT_DESCRIPTION | $SCRIPT_AUTHOR @ $SCRIPT_COPYRIGHT © $SCRIPT_LICENCE" >> $OUTPUT_PRIORITIZATION;
	echo "##"$SCRIPT_NAME"_RELEASE=$SCRIPT_RELEASE" >> $OUTPUT_PRIORITIZATION
	echo "##"$SCRIPT_NAME"_PARAM=$PARAM" >> $OUTPUT_PRIORITIZATION

	grep ^## $OUTPUT_PRIORITIZATION.old -v >> $OUTPUT_PRIORITIZATION

	rm $OUTPUT_PRIORITIZATION.old

fi;



# VCF final file COMPRESSION
if ((1)); then

	###################
	# VCF COMPRESSION #
	###################
	
	
	if [ $COMPRESS -gt -1 ]; then
	
		echo -e "\n####################\n# VCF COMPRESSION\n####################\n"
		echo "#[INFO] Compression level='$COMPRESS'";

		$BGZIP $OUTPUT_PRIORITIZATION -c -@ $THREADS -l 9 > $OUTPUT_PRIORITIZATION.gz
		OUTPUT_PRIORITIZATION=$OUTPUT_PRIORITIZATION.gz

	fi;

fi;



# VCF final file VALIDATION
if ((1)); then


	##################
	# VCF VALIDATION #
	##################

	echo -e "\n####################\n# VCF VALIDATION\n####################\n"

	if !(($(grep ^ $INPUT -c))); then
		echo "#[ERROR] Input VCF '$INPUT' is an empty file";
		exit 0;
	fi;

	VCF_VALIDATION=$TMP_FOLDER/VCF_VALIDATION
	>$VCF_VALIDATION.tmp
	>$VCF_VALIDATION
	#echo "$BCFTOOLS view $INPUT"
	TMP_SORT_VALIDATION=$TMP_SORT$RANDOM
	mkdir -p $TMP_SORT_VALIDATION
	$BCFTOOLS view $OUTPUT_PRIORITIZATION 2>>$VCF_VALIDATION.tmp | $BCFTOOLS sort -T $TMP_SORT_VALIDATION 1>/dev/null 2>>$VCF_VALIDATION.tmp
	rm -rf $TMP_SORT_VALIDATION
	#/home1/TOOLS/tools/bcftools/1.3.1/bin/bcftools view $INPUT 1>/dev/null 2>$VCF_VALIDATION
	grep -v -e "^Writing to " -e "^Merging .* temporary files$" -e "^Cleaning$" -e "^Done$" $VCF_VALIDATION.tmp > $VCF_VALIDATION
	rm $VCF_VALIDATION.tmp

	if (($(grep ^Failed $VCF_VALIDATION -c))) || (($(grep ^Aborded $VCF_VALIDATION -c))); then
		echo "#[ERROR] Input VCF $OUTPUT_PRIORITIZATION failed";
		cat $VCF_VALIDATION;
		exit 0;
	elif (($(grep "^\[E::" $VCF_VALIDATION -c))); then
		echo "#[ERROR] Some Errors on Final VCF $OUTPUT_PRIORITIZATION:";
		cat $VCF_VALIDATION;
		exit 0;
	elif (($(cat $VCF_VALIDATION | wc -l))); then
		echo "#[WARNING] Some Errors/Warnings on Final VCF $OUTPUT_PRIORITIZATION:";
		cat $VCF_VALIDATION;
	else
		echo "#[INFO] Final VCF validated";
	fi;


fi;



# Translation in any case (translation if no variant too! )
if ((1)); then

	####################
	# TRANSLATION
	####################

	if [ ! -z $FORMAT ] && [ "$FORMAT" != "" ]; then
		echo -e "\n####################\n# TRANSLATION\n####################\n"
		#$SCRIPT_DIR/VCFtranslation.pl --header
		#echo "$SCRIPT_DIR/VCFtranslation.pl $PARAM --input=$OUTPUT_PRIORITIZATION --output=$OUTPUT"

		$BCFTOOLS view $OUTPUT_PRIORITIZATION > $OUTPUT_PRIORITIZATION.uncompressed.vcf
		rm $OUTPUT_PRIORITIZATION
		OUTPUT_PRIORITIZATION=$OUTPUT_PRIORITIZATION.uncompressed.vcf

		#echo $FIELDS
		if (($THREADS)) && (($MULTITHREADING))  && [ $NB_VARIANT -gt $TRANSLATION_SPLIT ]; then
		#if (($THREADS)) && (($MULTITHREADING)); then


			# ECHO / VERBOSE / DEBUG
				if ((1)); then
				echo "#[INFO] Multithreading ($THREADS threads)"

					# NEEDED
					#if [ ! -z $BGZIP ] && [ -e $BZZIP ]; then
					#	echo "#[INFO] BGZIP=$BGZIP";
					#else
					#	echo "#[ERROR] BGZIP '$BGZIP' needed. Configure ENV file";
					#	exit 0;
					#fi;
					#if [ ! -z $VCFTOOLS ] && [ -e $VCFTOOLS ]; then
					#	echo "#[INFO] VCFTOOLS=$VCFTOOLS";
					#else
					#	echo "#[ERROR] VCFTOOLS '$VCFTOOLS' needed. Configure ENV file";
					#	exit 0;
					#fi;
					if [ ! -z $BCFTOOLS ] && [ -e $BCFTOOLS ]; then
						echo "#[INFO] BCFTOOLS=$BCFTOOLS";
					else
						echo "#[ERROR] BCFTOOLS '$BCFTOOLS' needed. Configure ENV file";
						exit 0;
					fi;
				fi;


				#echo $NB_VARIANT
				#echo $TRANSLATION_SPLIT
				VCF_SPLITTED_TRANS_LIST=""
				if [ $NB_VARIANT -gt $TRANSLATION_SPLIT ]; then
					# SPLIT VCF on ANNOTATIONS_SPLIT variants
					# HEAD
					#grep "^#" $OUTPUT_PRIORITIZATION > $TMP_FOLDER/input_prioritization.header
					$BCFTOOLS view -h $OUTPUT_PRIORITIZATION > $TMP_FOLDER/input_prioritization.header
					# VARIANTS
					#grep -v "^#" $OUTPUT_PRIORITIZATION > $TMP_FOLDER/input_prioritization.variants
					$BCFTOOLS view -H $OUTPUT_PRIORITIZATION > $TMP_FOLDER/input_prioritization.variants
					# SPLIT
					#echo "SPLIT"
					#grep -v "^#" $OUTPUT_PRIORITIZATION | split - -a $SPLIT_SUFFIT_LENGTH -l $TRANSLATION_SPLIT $TMP_FOLDER/input_prioritization.variants.splitted.
					$BCFTOOLS view -H $OUTPUT_PRIORITIZATION | split - -a $SPLIT_SUFFIT_LENGTH -l $TRANSLATION_SPLIT $TMP_FOLDER/input_prioritization.variants.splitted.
					#echo "SPLIT END"
					for splited_variants in $TMP_FOLDER/input_prioritization.variants.splitted.*; do
						cat $TMP_FOLDER/input_prioritization.header $splited_variants > $splited_variants.vcf
						VCF_SPLITTED_TRANS_LIST=$VCF_SPLITTED_TRANS_LIST" $splited_variants.vcf"
					done
					#echo "SPLITTED: $VCF_SPLITTED_TRANS_LIST"
					#cat $VCF_SPLITTED_TRANS_LIST
					#exit 0;
					echo "#[INFO] Input VCF splitted into "$(echo $VCF_SPLITTED_TRANS_LIST | wc -w)" files";
					#echo $VCF_SPLITTED_TRANS_LIST

				else
					# cp $OUTPUT_PRIORITIZATION $TMP_FOLDER/input_prioritization.variants.splitted.vcf
					$BCFTOOLS view $OUTPUT_PRIORITIZATION > $TMP_FOLDER/input_prioritization.variants.splitted.vcf
					VCF_SPLITTED_PRIO_LIST=$TMP_FOLDER/input_prioritization.variants.splitted.vcf
				fi;
				#exit 0;

				# PARAM ANNOTATION
				PARAM_TRANSLATION=$(echo "$PARAM " | sed "s/--snpeff_stats=[^ |$]*//gi");

				# FIND SAMPLE
				#if [ $(grep "#CHROM" $INPUT | wc -w) -gt 9 ]; then
				if [ $($BCFTOOLS view -h $INPUT | grep "#CHROM" | wc -w) -gt 9 ]; then
				#if [ $($BCFTOOLS view -h $INPUT | grep "#CHROM" | wc -w) -gt 9 ]; then
					SAMPLES=$($BCFTOOLS view -h $INPUT | grep "#CHROM" | cut -f10-$(grep "#CHROM" $INPUT | wc -w) --output-delimiter=,)
					#SAMPLES=$($BCFTOOLS view -h $INPUT | grep "#CHROM" | cut -f10-$($BCFTOOLS view -h $INPUT | grep "#CHROM" | wc -w) --output-delimiter=,)
				else
					SAMPLES=""
				fi;

				SAMPLES_COUNT=$(echo $SAMPLES | tr "," " " | wc -w)
				VCF_HEADER_CHROM_LENGTH=$(($SAMPLES_COUNT+9))
				#(($VERBOSE)) && echo "#[INFO] SAMPLES (#SAMPLES_COUNT)=$SAMPLES" 1>>$LOG 2>$ERR
				(($VERBOSE)) && echo "#[INFO] #SAMPLES=$SAMPLES_COUNT" 1>>$LOG 2>$ERR
				(($VERBOSE)) && echo "#[INFO] SAMPLES=$SAMPLES" 1>>$LOG 2>$ERR


				# IF SAMPLE(S)
				if [ "$SAMPLES " != "" ]; then

					#echo "#[INFO] TRANSLATION..."

					# ALL
					echo -e "all: $OUTPUT_TRANSLATION
					" >>$MK_TRANSLATION;

					for VCF_SPLITTED in $VCF_SPLITTED_TRANS_LIST; do
						TAB_SPLITTED_TRANS=$VCF_SPLITTED.trans.tab
						#echo "TAB_SPLITTED_TRANS=$TAB_SPLITTED_TRANS"
						echo "$TAB_SPLITTED_TRANS: $VCF_SPLITTED
						$SCRIPT_DIR/VCFtranslation.pl $PARAM_TRANSLATION --output=\$@ --input=$VCF_SPLITTED 1>>$LOG 2>>$ERR
						" >>$MK_TRANSLATION;
						TAB_SPLITTED_TRANS_LIST=$TAB_SPLITTED_TRANS_LIST" $TAB_SPLITTED_TRANS"
					done;
					#echo $TAB_SPLITTED_TRANS_LIST
					#cat $MK_TRANSLATION;

					echo "$OUTPUT_TRANSLATION: $TAB_SPLITTED_TRANS_LIST
							#echo 'test' > \$@
							head -n 1 \$< > \$@
							awk 'FNR>1' \$^ >>  \$@
							#join \$^ >>  \$@
						" >>$MK_TRANSLATION;

					if ! make -B -j $THREADS -f $MK_TRANSLATION all 1>>$LOG 2>>$ERR; then echo "#[FAILED] TRANSLATION Multithreading FAILED!!! see $MK_TRANSLATION file"; exit 1; fi;

					# ECHO / VERBOSE / DEBUG
					(($DEBUG)) && echo $MK_TRANSLATION && cat $MK_TRANSLATION;
					#(($DEBUG)) && tail -n 10 $OUTPUT_TRANSLATION.merge.vcf
					(($DEBUG)) && tail -n 10 $OUTPUT_TRANSLATION

					#cat $MK_TRANSLATION;
					#exit 0;

					# CLEANING
					rm -f $MK_TRANSLATION


				fi;

				#exit 0;

				#mv $OUTPUT_TRANSLATION $OUTPUT

		else

			#echo "$SCRIPT_DIR/VCFtranslation.pl $PARAM --input=$OUTPUT_PRIORITIZATION --output=$OUTPUT_TRANSLATION"
			if ! $SCRIPT_DIR/VCFtranslation.pl $PARAM --input=$OUTPUT_PRIORITIZATION --output=$OUTPUT_TRANSLATION 1>>$LOG 2>$ERR; then echo "#[FAILED] TRANSLATION FAILED!!! "; exit 1; fi;

		fi;

	else
		cp $OUTPUT_PRIORITIZATION $OUTPUT_TRANSLATION
	fi;

	# clear intermediate files
	rm -f $OUTPUT_PRIORITIZATION

	# Final OUTPUT
	mv $OUTPUT_TRANSLATION $OUTPUT

fi;


# ECHO / VERBOSE / DEBUG
if [ ! -e $OUTPUT ]; then
	echo -e "\n#[ERROR] No output '$OUTPUT' file. Please see LOG '$LOG' and ERR '$ERR' files\n";
fi;
#(($VERBOSE)) && echo -e "\n####################\n# LOG '$LOG'\n####################\n" && cat $LOG;
#(($VERBOSE)) && ! (($MULTITHREADING)) && echo -e "\n####################\n# LOG '$LOG'\n####################\n" && cat $LOG;
(($DEBUG)) && echo -e "\n####################\n# LOG '$LOG'\n####################\n" && cat $LOG;
(($DEBUG)) && echo -e "\n####################\n# ERR '$ERR' \n####################\n" && cat $ERR;
#if [ "$(cat $ERR)" != "" ] && [ ! $DEBUG ]; then cat $ERR; fi;
(($DEBUG)) && echo -e "\n####################\n# OUTOUT '$OUTPUT' \n####################\n" && tail -n 10 $OUTPUT

echo "";


# CLEANING
#rm -rf $TMP_FOLDER

# KILL verbose
if (($VERBOSE)); then
	kill $TAIL_LOG_PID
fi;
exit 0;
