#!/bin/bash
#################################
## HOWARD 
#################################

SCRIPT_NAME="HOWARD Validation"
SCRIPT_DESCRIPTION="Validation of HOWARD tool"
SCRIPT_RELEASE="0.9.3b"
SCRIPT_DATE="07/09/2018"
SCRIPT_AUTHOR="Antony Le Bechec"
SCRIPT_COPYRIGHT="IRC"
SCRIPT_LICENCE="GNU-GPL"

# Realse note
RELEASE_NOTES=$RELEASE_NOTES"# 0.9b-06/04/2018:\n";
RELEASE_NOTES=$RELEASE_NOTES"#\tScript creation\n";
RELEASE_NOTES=$RELEASE_NOTES"# 0.9.1b-17/05/2018:\n";
RELEASE_NOTES=$RELEASE_NOTES"#\tAdd header to report\n";
RELEASE_NOTES=$RELEASE_NOTES"#\tAdd validation CALC_NOMEN_TRANS_ALL and CALC_NOMEN_TRANS_ALL_X\n";
RELEASE_NOTES=$RELEASE_NOTES"# 0.9.2b-24/08/2018:\n";
RELEASE_NOTES=$RELEASE_NOTES"#\tChange input/output default values\n";
RELEASE_NOTES=$RELEASE_NOTES"#\tOption --norm added\n";
RELEASE_NOTES=$RELEASE_NOTES"#\tbug fixes\n";
RELEASE_NOTES=$RELEASE_NOTES"# 0.9.3b-07/09/2018:\n";
RELEASE_NOTES=$RELEASE_NOTES"#\tChange output file name by default for multiple input files\n";
RELEASE_NOTES=$RELEASE_NOTES"#\tAdd log files for each validation command\n";


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
	echo "# --input=<FILE>             File with test definitions (tab-delimiter with columns Test_name, File, Expected, Comment)";
	echo "# --output=<FILE>            Folder to provide results (default <INPUT>_VALIDATION)";
	echo "# --env=<FILE>               Environment configuration for multithreading (BGZIP, TABIX, BCFTOOLS)";
	echo "# --tmp=<FOLDER>             Temporary folder (default /tmp)";
	echo "# --norm=<FILE>              Genome fasta file to normalize (beware of chromosome identification, either 'x' or 'chrx')";
	echo "# --verbose                  VERBOSE option";
	echo "# --debug                    DEBUG option";
	echo "# --release                  RELEASE option";
	echo "# --help                     HELP option";
	echo "#";
	echo "# Reconnized tests:";
	echo "# - Multithreading";
	echo "#";
	echo "";

}

# EXAMPLE :
# ./HOWARD_VALIDATION.sh --vcf=validation --output=validation_results 

# header
header;



ARGS=$(getopt -o "i:o:e:m:vdnh" --long "input:,output:,env:,tmp:,norm:,verbose,debug,release,help" -- "$@" 2> /dev/null)
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
		--input)
			if [ ! -e $2 ] || [ "$2" == "" ]; then
				echo "#[ERROR] No INPUT file '$2'"
				usage; exit;
			else
				INPUT="$2";
			fi;
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
		--tmp)
			TMP_INPUT="$2"
			shift 2
			;;
		--norm)
			NORM="$2"
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
if [ -z $INPUT ] || [ ! -f $INPUT ]; then #  || [ -z $OUTPUT ]; then
	echo "#[WARNING] INPUT Missing";
	if [ -f $SCRIPT_DIR/validation/validation.txt ]; then
		INPUT=$SCRIPT_DIR/validation/validation.txt;
		echo "#[INFO] INPUT '$INPUT' used."
	else
		echo "#[ERROR] No INPUT found (default 'validation/validation.txt' missing)";
		usage;
		exit;
	fi;
	
fi

INPUT_DIR=$(dirname $INPUT)

# Mandatory parameters
if [ -z $OUTPUT ] || [ "$OUTPUT" == "" ]; then #  || [ -z $OUTPUT ]; then
	echo "#[WARNING] OUTPUT Missing";
	#OUTPUT=$INPUT_DIR"_VALIDATION"; #/\.vcf$/\.output\.vcf/;
	DATE="V"`date '+%Y%m%d-%H%M%S'`
	OUTPUT=$INPUT"_VALIDATION_"$DATE; #/\.vcf$/\.output\.vcf/;
	echo "#[INFO] OUTPUT '$OUTPUT' used."
	
fi

# Normalization
if [ -z $NORM ]; then
	NORM="";
fi;



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

TMP_VAR=$RANDOM$RANDOM
TMP_FOLDER=$TMP_SYS_FOLDER/HOWARD_VALIDATION_$TMP_VAR
#PARAM=$PARAM" --tmp=$TMP_FOLDER "


# THREADS
THREADS=$(ls -d /sys/devices/system/cpu/cpu[[:digit:]]* | wc -w)

PARAM_ADD=" --force --tmp=$TMP_FOLDER --compress=9 "
if [ "$NORM" != "" ] && [ -e "$NORM" ]; then
	PARAM_ADD=$PARAM_ADD" --norm=$NORM "
fi;

VALIDATION_RESULTS_DIR=$OUTPUT

mkdir -p $VALIDATION_RESULTS_DIR

VALIDATION_REPORT=$VALIDATION_RESULTS_DIR/REPORT
header > $VALIDATION_REPORT

echo -e "# INPUT=$INPUT" >> $VALIDATION_REPORT
echo -e "# OUTPUT=$OUTPUT" >> $VALIDATION_REPORT
echo -e "# ENV=$ENV" >> $VALIDATION_REPORT
echo -e "# TMP=$TMP_FOLDER" >> $VALIDATION_REPORT
echo -e "#Test Name\tResult\tFile(s)\tError\tTest Description\tComment" >> $VALIDATION_REPORT


LOG=$VALIDATION_RESULTS_DIR/log
ERR=$VALIDATION_RESULTS_DIR/err
LOG_TEST=$VALIDATION_RESULTS_DIR/log_test
ERR_TEST=$VALIDATION_RESULTS_DIR/err_test
#RES_TEST=$VALIDATION_RESULTS_DIR/err_test

(($VERBOSE)) && echo "" 


while IFS=$'\t' read -r -a myArray
do
	TEST_NAME="${myArray[0]}"
	FILES=$INPUT_DIR/"${myArray[1]}"
	#echo $FILES
	FILES=${myArray[1]}
	EXPECTED="${myArray[2]}"
	COMMENT="${myArray[3]}"
	
	# INPUT FILES
	FILE=""
	FILE_OUTPUT=""
	for F in $(echo $FILES | tr "," " "); do
		FILE=$FILE" "$INPUT_DIR/$F
		FILE_OUTPUT=$FILE_OUTPUT"_"$F
	done;
	# Clean files
	FILE=$(echo $FILE | sed "s/^ *//gi" | sed "s/ *$//gi")
	FILE_OUTPUT=$(echo $FILE_OUTPUT | sed "s/^_*//gi" | sed "s/_*$//gi")
	
	CMD=""
	
	#echo "$TEST_NAME|$FILE|$EXPECTED|$COMMENT";
	#continue
	
	#(($VERBOSE)) && echo -e "\n# TEST $TEST_NAME\n" 
	#echo -e "\n# TEST $TEST_NAME\n" >> $LOG
	
	if [ "$TEST_NAME" != "" ]; then
	
		mkdir -p $VALIDATION_RESULTS_DIR/$TEST_NAME
		
		RES_OUTPUT=$(echo $VALIDATION_RESULTS_DIR/$TEST_NAME/$(basename $FILE_OUTPUT) | tr "," "_" | tr " " "_")
		
		case "$TEST_NAME" in
		
			### ANNOTATION test ###
			
			ANN)
				TEST_DESCRIPTION="Check annotation using CORE annotation"
				# COMMAND
				CMD="$SCRIPT_DIR/HOWARD.sh --vcf='$FILE' --output=$RES_OUTPUT --annotation=CORE $PARAM_ADD"
				;;
				
			ANN_CORE)
				TEST_DESCRIPTION="Check annotation using CORE annotation"
				# COMMAND
				CMD="$SCRIPT_DIR/HOWARD.sh --vcf='$FILE' --output=$RES_OUTPUT --annotation=CORE $PARAM_ADD"
				;;
		
			ANN_F)
				TEST_DESCRIPTION="Check annotation using CORE annotation and forcing annotation"
				# COMMAND
				CMD="$SCRIPT_DIR/HOWARD.sh --vcf='$FILE' --output=$RES_OUTPUT --annotation=CORE --force $PARAM_ADD"
				;;
		
		
			ANN_X)
				TEST_DESCRIPTION="Check annotation using CORE annotation with multithreading "
				# COMMAND
				CMD="$SCRIPT_DIR/HOWARD.sh --vcf='$FILE' --output=$RES_OUTPUT --annotation=CORE --multithreading $PARAM_ADD"
				;;
				
			ANN_CORE_X)
				TEST_DESCRIPTION="Check annotation using CORE annotations with multithreading"
				# COMMAND
				CMD="$SCRIPT_DIR/HOWARD.sh --vcf='$FILE' --output=$RES_OUTPUT --annotation=CORE --multithreading $PARAM_ADD"
				;;
			
			ANN_SYMBOL_X)
				TEST_DESCRIPTION="Check annotation using uniq annotation symbol with multithreading"
				# COMMAND
				CMD="$SCRIPT_DIR/HOWARD.sh --vcf='$FILE' --output=$RES_OUTPUT --annotation=symbol --multithreading $PARAM_ADD"
				;;
		
			ANN_ALL_X)
				TEST_DESCRIPTION="Check annotation using ALL annotations with multithreading"
				# COMMAND
				CMD="$SCRIPT_DIR/HOWARD.sh --vcf='$FILE' --output=$RES_OUTPUT --annotation=ALL --multithreading $PARAM_ADD"
				;;
		
			## SNPEff
		
			ANN_SNPEFF)
				TEST_DESCRIPTION="Check annotation SNPEff through --annotation option"
				# COMMAND
				CMD="$SCRIPT_DIR/HOWARD.sh --vcf='$FILE' --output=$RES_OUTPUT --annotation=snpeff $PARAM_ADD"
				;;
				
			ANN_SNPEFF_FULL)
				TEST_DESCRIPTION="Check annotation SNPEff, SNPEff HGVS through --annotation option and SNPEff Stats trhough option"
				# COMMAND
				CMD="$SCRIPT_DIR/HOWARD.sh --vcf='$FILE' --output=$RES_OUTPUT --annotation=snpeff,snpeff_hgvs --snpeff_stats=$RES_OUTPUT.stats $PARAM_ADD"
				;;
				
			ANN_SNPEFF_OPT)
				TEST_DESCRIPTION="Check annotation SNPEff through --snpeff option"
				# COMMAND
				CMD="$SCRIPT_DIR/HOWARD.sh --vcf='$FILE' --output=$RES_OUTPUT --snpeff $PARAM_ADD"
				;;
				
			ANN_SNPEFF_HGVS_OPT)
				TEST_DESCRIPTION="Check annotation SNPEff HGVS through --snpeff option"
				# COMMAND
				CMD="$SCRIPT_DIR/HOWARD.sh --vcf='$FILE' --output=$RES_OUTPUT --snpeff_hgvs $PARAM_ADD"
				;;
				
			ANN_SNPEFF_STATS_OPT)
				TEST_DESCRIPTION="Check annotation SNPEff STATS through --snpeff_stats option"
				# COMMAND
				CMD="$SCRIPT_DIR/HOWARD.sh --vcf='$FILE' --output=$RES_OUTPUT --snpeff_stats=$RES_OUTPUT.stats $PARAM_ADD"
				;;
		
			ANN_SNPEFF_X)
				TEST_DESCRIPTION="Check annotation SNPEff through --annotation option with multithreading"
				# COMMAND
				CMD="$SCRIPT_DIR/HOWARD.sh --vcf='$FILE' --output=$RES_OUTPUT --annotation=snpeff --multithreading $PARAM_ADD"
				;;
				
			ANN_SNPEFF_FULL_X)
				TEST_DESCRIPTION="Check annotation SNPEff, SNPEff HGVS through --annotation option and SNPEff Stats through option with multithreading"
				# COMMAND
				CMD="$SCRIPT_DIR/HOWARD.sh --vcf='$FILE' --output=$RES_OUTPUT --annotation=snpeff,snpeff_hgvs --snpeff_stats=$RES_OUTPUT.stats --multithreading $PARAM_ADD"
				;;
				
			ANN_SNPEFF_OPT_X)
				TEST_DESCRIPTION="Check annotation SNPEff through --snpeff option with multithreading"
				# COMMAND
				CMD="$SCRIPT_DIR/HOWARD.sh --vcf='$FILE' --output=$RES_OUTPUT --snpeff --multithreading $PARAM_ADD"
				;;
				
			ANN_SNPEFF_HGVS_OPT_X)
				TEST_DESCRIPTION="Check annotation SNPEff HGVS through --snpeff option with multithreading"
				# COMMAND
				CMD="$SCRIPT_DIR/HOWARD.sh --vcf='$FILE' --output=$RES_OUTPUT --snpeff_hgvs --multithreading $PARAM_ADD"
				;;
				
			ANN_SNPEFF_STATS_OPT_X)
				TEST_DESCRIPTION="Check annotation SNPEff STATS through --snpeff_stats option with multithreading"
				# COMMAND
				CMD="$SCRIPT_DIR/HOWARD.sh --vcf='$FILE' --output=$RES_OUTPUT --snpeff_stats=$RES_OUTPUT.stats --multithreading $PARAM_ADD"
				;;
		
			### CALCULATION TEST ###
			
			CALC_ALL)
				TEST_DESCRIPTION="Check ALL calculations VAF,VAF_STATS,CALLING_QUALITY,CALLING_QUALITY_EXPLODE,NOMEN,BARCODE,GENOTYPECONCORDANCE,FINDBYPIPELINES"
				# COMMAND
				CMD="$SCRIPT_DIR/HOWARD.sh --vcf='$FILE' --output=$RES_OUTPUT --calculation=VAF,VAF_STATS,CALLING_QUALITY,CALLING_QUALITY_EXPLODE,NOMEN,BARCODE,GENOTYPECONCORDANCE,FINDBYPIPELINES $PARAM_ADD"
				;;
				
			CALC_ALL_X)
				TEST_DESCRIPTION="Check ALL calculations VAF,VAF_STATS,CALLING_QUALITY,CALLING_QUALITY_EXPLODE,NOMEN,BARCODE,GENOTYPECONCORDANCE,FINDBYPIPELINES with multithreading"
				# COMMAND
				CMD="$SCRIPT_DIR/HOWARD.sh --vcf='$FILE' --output=$RES_OUTPUT --calculation=VAF,VAF_STATS,CALLING_QUALITY,CALLING_QUALITY_EXPLODE,NOMEN,BARCODE,GENOTYPECONCORDANCE,FINDBYPIPELINES --multithreading $PARAM_ADD"
				;;

			
			CALC_VAF)
				TEST_DESCRIPTION="Check calculation VAF"
				# COMMAND
				CMD="$SCRIPT_DIR/HOWARD.sh --vcf='$FILE' --output=$RES_OUTPUT --calculation=VAF $PARAM_ADD"
				;;
			
			CALC_VAFSTAT)
				TEST_DESCRIPTION="Check calculation VAF STAT"
				# COMMAND
				CMD="$SCRIPT_DIR/HOWARD.sh --vcf='$FILE' --output=$RES_OUTPUT --calculation=VAF_STATS $PARAM_ADD"
				;;
				
			CALC_VAF_VAFSTAT)
				TEST_DESCRIPTION="Check calculation VAF and VAF STAT"
				# COMMAND
				CMD="$SCRIPT_DIR/HOWARD.sh --vcf='$FILE' --output=$RES_OUTPUT --calculation=VAF,VAF_STATS $PARAM_ADD"
				;;
			
			CALC_CQ)
				TEST_DESCRIPTION="Check calculation CALLING QUALITY"
				# COMMAND
				CMD="$SCRIPT_DIR/HOWARD.sh --vcf='$FILE' --output=$RES_OUTPUT --calculation=CALLING_QUALITY $PARAM_ADD"
				;;
			
			CALC_CQE)
				TEST_DESCRIPTION="Check calculation CALLING QUALITY EXPLODE"
				# COMMAND
				CMD="$SCRIPT_DIR/HOWARD.sh --vcf='$FILE' --output=$RES_OUTPUT --calculation=CALLING_QUALITY_EXPLODE $PARAM_ADD"
				;;
			
			CALC_NOMEN)
				TEST_DESCRIPTION="Check calculation NOMEN"
				# COMMAND
				CMD="$SCRIPT_DIR/HOWARD.sh --vcf='$FILE' --output=$RES_OUTPUT --calculation=NOMEN $PARAM_ADD"
				;;
			
			CALC_BC)
				TEST_DESCRIPTION="Check calculation BARCODE"
				# COMMAND
				CMD="$SCRIPT_DIR/HOWARD.sh --vcf='$FILE' --output=$RES_OUTPUT --calculation=BARCODE $PARAM_ADD"
				;;
			
			CALC_GC)
				TEST_DESCRIPTION="Check calculation GENOTYPECONCORDANCE"
				# COMMAND
				CMD="$SCRIPT_DIR/HOWARD.sh --vcf='$FILE' --output=$RES_OUTPUT --calculation=GENOTYPECONCORDANCE $PARAM_ADD"
				;;
			
			CALC_FBP)
				TEST_DESCRIPTION="Check calculation FINDBYPIPELINES"
				# COMMAND
				CMD="$SCRIPT_DIR/HOWARD.sh --vcf='$FILE' --output=$RES_OUTPUT --calculation=FINDBYPIPELINES $PARAM_ADD"
				;;
			
			CALC_VAF_X)
				TEST_DESCRIPTION="Check calculation VAF with multithreading"
				# COMMAND
				CMD="$SCRIPT_DIR/HOWARD.sh --vcf='$FILE' --output=$RES_OUTPUT --calculation=VAF --multithreading $PARAM_ADD"
				;;
			
			CALC_VAFSTAT_X)
				TEST_DESCRIPTION="Check calculation VAF STAT with multithreading"
				# COMMAND
				CMD="$SCRIPT_DIR/HOWARD.sh --vcf='$FILE' --output=$RES_OUTPUT --calculation=VAF_STATS --multithreading $PARAM_ADD"
				;;
				
			CALC_VAF_VAFSTAT_X)
				TEST_DESCRIPTION="Check calculation VAF and VAF STAT with multithreading"
				# COMMAND
				CMD="$SCRIPT_DIR/HOWARD.sh --vcf='$FILE' --output=$RES_OUTPUT --calculation=VAF,VAF_STATS --multithreading $PARAM_ADD"
				;;
			
			CALC_CQ_X)
				TEST_DESCRIPTION="Check calculation CALLING QUALITY with multithreading"
				# COMMAND
				CMD="$SCRIPT_DIR/HOWARD.sh --vcf='$FILE' --output=$RES_OUTPUT --calculation=CALLING_QUALITY --multithreading $PARAM_ADD"
				;;
			
			CALC_CQE_X)
				TEST_DESCRIPTION="Check calculation CALLING QUALITY EXPLODE with multithreading"
				# COMMAND
				CMD="$SCRIPT_DIR/HOWARD.sh --vcf='$FILE' --output=$RES_OUTPUT --calculation=CALLING_QUALITY_EXPLODE --multithreading $PARAM_ADD"
				;;
			
			CALC_NOMEN_X)
				TEST_DESCRIPTION="Check calculation NOMEN with multithreading"
				# COMMAND
				CMD="$SCRIPT_DIR/HOWARD.sh --vcf='$FILE' --output=$RES_OUTPUT --calculation=NOMEN --multithreading $PARAM_ADD"
				;;
			
			CALC_BC_X)
				TEST_DESCRIPTION="Check calculation BARCODE with multithreading"
				# COMMAND
				CMD="$SCRIPT_DIR/HOWARD.sh --vcf='$FILE' --output=$RES_OUTPUT --calculation=BARCODE --multithreading $PARAM_ADD"
				;;
			
			CALC_GC_X)
				TEST_DESCRIPTION="Check calculation GENOTYPECONCORDANCE with multithreading"
				# COMMAND
				CMD="$SCRIPT_DIR/HOWARD.sh --vcf='$FILE' --output=$RES_OUTPUT --calculation=GENOTYPECONCORDANCE --multithreading $PARAM_ADD"
				;;
			
			CALC_FBP_X)
				TEST_DESCRIPTION="Check calculation FINDBYPIPELINES with multithreading"
				# COMMAND
				CMD="$SCRIPT_DIR/HOWARD.sh --vcf='$FILE' --output=$RES_OUTPUT --calculation=FINDBYPIPELINES --multithreading $PARAM_ADD"
				;;
				
			CALC_NOMEN_TRANS_ALL)
				TEST_DESCRIPTION="Check calculation NOMEN with transcripts file"
				VCF=$(echo $FILE | cut -d" " -f1)
				TRANSCRIPTS=$(echo $FILE | cut -d" " -f2)
				# COMMAND
				CMD="$SCRIPT_DIR/HOWARD.sh --vcf=$VCF --output=$RES_OUTPUT --calculation=NOMEN --transcripts=$TRANSCRIPTS $PARAM_ADD"
				;;
			
			CALC_NOMEN_TRANS_ALL_X)
				TEST_DESCRIPTION="Check calculation NOMEN with transcripts file multithreading"
				VCF=$(echo $FILE | cut -d" " -f1)
				TRANSCRIPTS=$(echo $FILE | cut -d" " -f2)
				# COMMAND
				CMD="$SCRIPT_DIR/HOWARD.sh --vcf=$VCF --output=$RES_OUTPUT --calculation=NOMEN --transcripts=$TRANSCRIPTS --multithreading $PARAM_ADD"
				;;
			
			### INTERGRATION TEST ###
		
			INT_CORE)
				TEST_DESCRIPTION="Check annotation, calculation, prioritization and translation using CORE annotations, ALL calculation, default prioritization, tab format of output, without multithreading"
				# COMMAND
				CMD="$SCRIPT_DIR/HOWARD.sh --vcf='$FILE' --output=$RES_OUTPUT --annotation=CORE --calculation=VAF,VAF_STATS,CALLING_QUALITY,CALLING_QUALITY_EXPLODE,NOMEN,BARCODE,GENOTYPECONCORDANCE,FINDBYPIPELINES --prioritization=default --translation=tab $PARAM_ADD"
				;;
				
			INT_CORE_X)
				TEST_DESCRIPTION="Check annotation, calculation, prioritization and translation using CORE annotations, ALL calculation, default prioritization, tab format of output, with multithreading"
				# COMMAND
				CMD="$SCRIPT_DIR/HOWARD.sh --vcf='$FILE' --output=$RES_OUTPUT --annotation=CORE --calculation=VAF,VAF_STATS,CALLING_QUALITY,CALLING_QUALITY_EXPLODE,NOMEN,BARCODE,GENOTYPECONCORDANCE,FINDBYPIPELINES --prioritization=default --translation=tab --multithreading $PARAM_ADD"
				;;
				
			INT_CORE_SNPEFF)
				TEST_DESCRIPTION="Check annotation, calculation, prioritization and translation using CORE annotations, ALL calculation, default prioritization, tab format of output, snpeff options, without multithreading"
				# COMMAND
				CMD="$SCRIPT_DIR/HOWARD.sh --vcf='$FILE' --output=$RES_OUTPUT --annotation=CORE --calculation=VAF,VAF_STATS,CALLING_QUALITY,CALLING_QUALITY_EXPLODE,NOMEN,BARCODE,GENOTYPECONCORDANCE,FINDBYPIPELINES --snpeff --snpeff_hgvs --snpeff_stats=$RES_OUTPUT.stats --prioritization=default --translation=tab $PARAM_ADD"
				;;
				
			INT_CORE_SNPEFF_X)
				TEST_DESCRIPTION="Check annotation, calculation, prioritization and translation using CORE annotations, ALL calculation, default prioritization, tab format of output, with multithreading"
				# COMMAND
				CMD="$SCRIPT_DIR/HOWARD.sh --vcf='$FILE' --output=$RES_OUTPUT --annotation=CORE --calculation=VAF,VAF_STATS,CALLING_QUALITY,CALLING_QUALITY_EXPLODE,NOMEN,BARCODE,GENOTYPECONCORDANCE,FINDBYPIPELINES --snpeff --snpeff_hgvs --snpeff_stats=$RES_OUTPUT.stats --prioritization=default --translation=tab --multithreading $PARAM_ADD"
				;;
			
			INT_FULL)
				TEST_DESCRIPTION="Check annotation, calculation, prioritization and translation using ALL annotations, ALL calculation, default prioritization, tab format of output, snpeff options, without multithreading"
				# COMMAND
				CMD="$SCRIPT_DIR/HOWARD.sh --vcf='$FILE' --output=$RES_OUTPUT --annotation=ALL --calculation=VAF,VAF_STATS,CALLING_QUALITY,CALLING_QUALITY_EXPLODE,NOMEN,BARCODE,GENOTYPECONCORDANCE,FINDBYPIPELINES --prioritization=default --translation=tab $PARAM_ADD"
				;;
				
			INT_FULL_X)
				TEST_DESCRIPTION="Check annotation, calculation, prioritization and translation using ALL annotations, ALL calculation, default prioritization, tab format of output, with multithreading"
				# COMMAND
				CMD="$SCRIPT_DIR/HOWARD.sh --vcf='$FILE' --output=$RES_OUTPUT --annotation=ALL --calculation=VAF,VAF_STATS,CALLING_QUALITY,CALLING_QUALITY_EXPLODE,NOMEN,BARCODE,GENOTYPECONCORDANCE,FINDBYPIPELINES --prioritization=default --translation=tab --multithreading $PARAM_ADD"
				;;
				
			INT_FULL_SNPEFF)
				TEST_DESCRIPTION="Check annotation, calculation, prioritization and translation using ALL annotations, ALL calculation, default prioritization, tab format of output, snpeff options, without multithreading"
				# COMMAND
				CMD="$SCRIPT_DIR/HOWARD.sh --vcf='$FILE' --output=$RES_OUTPUT --annotation=ALL --calculation=VAF,VAF_STATS,CALLING_QUALITY,CALLING_QUALITY_EXPLODE,NOMEN,BARCODE,GENOTYPECONCORDANCE,FINDBYPIPELINES --snpeff --snpeff_hgvs --snpeff_stats=$RES_OUTPUT.stats --prioritization=default --translation=tab $PARAM_ADD"
				;;
				
			INT_FULL_SNPEFF_X)
				TEST_DESCRIPTION="Check annotation, calculation, prioritization and translation using ALL annotations, ALL calculation, default prioritization, tab format of output, with multithreading"
				# COMMAND
				CMD="$SCRIPT_DIR/HOWARD.sh --vcf='$FILE' --output=$RES_OUTPUT --annotation=ALL --calculation=VAF,VAF_STATS,CALLING_QUALITY,CALLING_QUALITY_EXPLODE,NOMEN,BARCODE,GENOTYPECONCORDANCE,FINDBYPIPELINES --snpeff --snpeff_hgvs --snpeff_stats=$RES_OUTPUT.stats --prioritization=default --translation=tab --multithreading $PARAM_ADD"
				;;
		
			*) 	echo "# Test '$TEST_NAME' is not recognized."
				;;
		esac
	
		(($VERBOSE)) && echo "## $TEST_NAME"
		(($VERBOSE)) && echo "# $TEST_DESCRIPTION"
		#(($VERBOSE)) && echo "# Test name: '$TEST_NAME'"
		(($VERBOSE)) && echo "# File: '$FILE'"
		(($VERBOSE)) && echo "# Expected: '$EXPECTED'"
		(($VERBOSE)) && echo "# Comment: '$COMMENT'"
		echo "## $TEST_NAME"  >> $LOG
		echo "# $TEST_DESCRIPTION"  >> $LOG
		echo -e "# Test name: '$TEST_NAME'"  >> $LOG
		echo -e "# File: '$FILE'"  >> $LOG
		echo -e "# Expected: '$EXPECTED'"  >> $LOG
		echo -e "# Comment: '$COMMENT'"  >> $LOG
	
		if [ "$CMD" != "" ]; then
	
			(($VERBOSE)) && echo "# Command: "$CMD
			echo $CMD >> $LOG
			CMD_TO_EVAL="time  $CMD"
			LOG_TIME=$LOG_TEST"TIME"
			#echo $CMD
			eval "time $CMD  1>$LOG_TEST 2>$ERR_TEST" 2>$LOG_TIME  2>$LOG_TIME #_TO_EVAL 1>$LOG_TEST 2>$ERR_TEST
			
			cat $LOG_TEST >> $LOG; cat $LOG_TIME >> $LOG; cat $ERR_TEST >> $ERR
			cat $LOG_TEST >> $RES_OUTPUT.log; cat $LOG_TIME >> $RES_OUTPUT.time.log; cat $ERR_TEST >> $RES_OUTPUT.err
			(($VERBOSE)) && echo "# Exec Time:: "; grep ^$ -v $LOG_TIME;
			#echo "LOG :"; cat $LOG;
			#echo "ERR :"; cat $ERR;
			#echo "LOG TIME:"; cat $LOG_TIME;
			#echo "";
			# exit 0;
			
			# Validation
			RES_TEST_ERR=$(cat $LOG_TEST $ERR_TEST | grep -c "FAILED")	
			if (($RES_TEST_ERR)); then RES_TEST="FAILED"; else RES_TEST="OK"; fi; 
			if (($(echo $EXPECTED | tr "," "\n" | grep "^[[:space:]]*OUTPUT[[:space:]]*$" - -c))); then
				if [ ! -e $RES_OUTPUT ]; then RES_TEST="FAILED"; fi;
			fi;
			if (($(echo $EXPECTED | tr "," "\n" | grep "^[[:space:]]*FAILED[[:space:]]*$" - -c))); then
				if (($RES_TEST_ERR)); then RES_TEST="OK"; fi;
			fi;
			echo -e "$TEST_NAME\t$FILE\t$RES_TEST\t$RES_TEST_ERR\t$TEST_DESCRIPTION\t$COMMENT" >> $VALIDATION_REPORT
			(($VERBOSE)) && grep FAILED $LOG_TEST $ERR_TEST && echo ""
			
	
		fi;
	
		(($VERBOSE)) && echo ""
		echo ""  >> $LOG
	
	fi;
	
done < <(grep ^# -v $INPUT)


# SORT REPORT
mv $VALIDATION_REPORT $VALIDATION_REPORT.old
sort $VALIDATION_REPORT.old > $VALIDATION_REPORT
rm $VALIDATION_REPORT.old

#(($VERBOSE)) && column -t $VALIDATION_REPORT

(($VERBOSE)) && awk -F'\t' '{printf "%-30s|%-30s|%-6s|%-5s|%-1s - %-50s\n",$1,$2,$3,$4,$5,$6}' $VALIDATION_REPORT # > final_report.txt

(($DEBUG)) && cat $LOG
(($DEBUG)) && cat $ERR


# CLEANING
rm -rf $TMP_FOLDER



exit 0;


