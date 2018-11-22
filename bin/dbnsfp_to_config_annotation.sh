#!/bin/bash
#################################
## HOWARD 
#################################

SCRIPT_NAME="HOWARD_DBNSFP"
SCRIPT_DESCRIPTION="HOWARD DBNSFPx to config annotation ini file"
SCRIPT_RELEASE="0.9"
SCRIPT_DATE="07/11/2017"
SCRIPT_AUTHOR="Antony Le Bechec"
SCRIPT_COPYRIGHT="IRC"
SCRIPT_LICENCE="GNU-GPL"

# Realse note
RELEASE_NOTES=$RELEASE_NOTES"# 0.9b-07/11/2017:\n";
RELEASE_NOTES=$RELEASE_NOTES"#\tScript creation\n";


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
	echo "# --input=<FILE>             Input file in VCF format";
	echo "# --output=<FILE>            Output annotated file in defined format (default VCF)";
	echo "# --input_release=<STRING>   Release of the inpur DBNSFP file";
	echo "# --env=<FILE>               Environment configuration for multithreading (BGZIP, TABIX, BCFTOOLS, VCFTOOLS)";
	#echo "# --force                    Force annotation even if already exists in VCF header";
	echo "# --tmp=<FOLDER>             Temporary folder (default /tmp)";
	echo "# --verbose                  VERBOSE option";
	echo "# --debug                    DEBUG option";
	echo "# --release                  RELEASE option";
	echo "# --help                     HELP option";
	echo "#";
	echo "";

}



# EXAMPLE :
# ./dbnsfp_to_config_annotation.sh --input=dbnsfp33a.txt --output=config.annotation.dbnsfp33a.ini

# header
header;


ARGS=$(getopt -o "i:o:e:a:f:s:r:xt:m:vdnh" --long "input:,output:,input_name:,input_release:,input_date:,env:,annotation:,calculation:,filter:,format:,snpeff_stats:,multithreading,threads:,split:,tmp:,verbose,debug,release,help" -- "$@" 2> /dev/null)
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
				echo "#[ERROR] No VCF file '$2'"
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
		--input_name)
			INPUT_NAME="$2";
			shift 2
			;;
		--input_release)
			INPUT_RELEASE="$2";
			shift 2
			;;
		--input_date)
			INPUT_DATE="$2";
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
		--calculation)
			CALCULATION="$2"
			shift 2
			;;
		--config)
			CONFIG="$2"
			shift 2
			;;
		--filter)
			FILTER="$2"
			shift 2
			;;
		--format)
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
if [ -z $INPUT ]; then #  || [ -z $OUTPUT ]; then
	echo "# INPUT Missing";
	usage;
	exit;
fi

if [ -z $INPUT_NAME ]; then #  || [ -z $OUTPUT ]; then
	INPUT_NAME=$(basename $INPUT)
	echo "#[WARNING] INPUT NAME Missing. default '$INPUT_NAME' used";
fi

if [ -z $INPUT_RELEASE ]; then #  || [ -z $OUTPUT ]; then
	INPUT_RELEASE=$(basename $INPUT)
	echo "#[WARNING] INPUT RELEASE Missing. default '$INPUT_RELEASE' used";
fi

if [ -z $INPUT_DATE ]; then #  || [ -z $OUTPUT ]; then
	INPUT_DATE=$(basename $INPUT)
	echo "#[WARNING] INPUT DATE Missing. default '$INPUT_DATE' used";
fi


# Mandatory parameters
if [ -z $OUTPUT ] || [ "$OUTPUT" == "" ]; then #  || [ -z $OUTPUT ]; then
	echo "# OUTPUT Missing";
	#OUTPUT=$(echo $INPUT | sed "s/.vcf$/.annotated.vcf/g"); #/\.vcf$/\.output\.vcf/;
	OUTPUT=config.annotation.$(basename $INPUT).ini; #/\.vcf$/\.output\.vcf/;
fi


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



# HEADER
echo "
##################
# INPUT=$INPUT
# OUTPUT=$OUTPUT
# MULTITHREADING=$MULTITHREADING ($THREADS threads)"
(($VERBOSE)) && echo "# LOG=$LOG";
(($VERBOSE)) && echo "# ERR=$ERR";
echo "##################"

> $OUTPUT

line=6
for ANNOTATION in $(head -n1 $INPUT | cut -f6- | tr "\t" "\n"); do
	annotation_type="annotation"
	((line++))
	otherinfo=$(($line-6))
	if [[ $ANNOTATION =~ ^.*_score$ ]] || [[ $ANNOTATION =~ ^.*_rankscore$ ]] || [[ $ANNOTATION =~ ^.*_phred$ ]]; then
		annotation_type="score"
	elif [[ $ANNOTATION =~ ^.*_pred$ ]]; then
		annotation_type="prediction"
	fi
	
	echo "ANNOTATION $ANNOTATION [$annotation_type]";
	echo "

[$ANNOTATION]
annovar_code=$INPUT_NAME
annovar_annotation_type=filter
release=$INPUT_RELEASE
available=true
core=false
annotation_type=$annotation_type
date=$INPUT_DATE
otherinfo=$otherinfo
description=$ANNOTATION

" >> $OUTPUT

	
done;


