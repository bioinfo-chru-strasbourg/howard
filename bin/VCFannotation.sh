#!/bin/bash
#################################
##
## NGS environment
##
#################################

SCRIPT_NAME="HOWARDAnnotation"
SCRIPT_DESCRIPTION="HOWARD Annotation based on ANNOVAR and snpEff, allowing multithreading"
SCRIPT_RELEASE="0.9b"
SCRIPT_DATE="07/10/2016"
SCRIPT_AUTHOR="Antony Le Bechec"
SCRIPT_COPYRIGHT="IRC"
SCRIPT_LICENCE="GNU-GPL"

# Realse note
RELEASE_NOTES=$RELEASE_NOTES"# 0.9b-30/05/2016: Script creation\n";


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
	echo "# USAGE: $(basename $0) --input=<VCF>  [options...]";
	echo "# Following options are enable for this script $(basename $0)";
	echo "# -i/--input                              Input VCF";
	echo "# -o/--output                             Output VCF annotated";
	echo "# -a/--annotation                         Annotation";
	echo "# -x/--multithreading                     Multithreading";
	echo "# -t/--threads                            Threads";
	echo "# -v/--verbose                            VERBOSE option";
	echo "# -d/--debug                              DEBUG option";
	echo "# -n/--release                            RELEASE option";
	echo "# -h/--help                               HELP option";
	echo "#";
	echo "# Following options are also anable with the used script VCFannotation.pl";
	$SCRIPT_DIR/VCFannotation.pl --help


: "	echo # -o/--output                             Output VCF annotated;
	echo # -a/--annotation                         Annotation;
	echo # -e/--config                             Configuration file;
	echo # -f/--config_annotation                  Configuratin annotation file;
	echo # -y/--assembly                           Assembly;
	echo # -r/--annovar_folder                     ANNOVAR folder including scripts;
	echo # -b/--annovar_databases                  ANNOVAR databases;
	echo # -s/--snpeff                             snpEff annotation;
	echo # -m/--snpeff_jar                         snpEff JAR;
	echo # -l/--snpeff_database                    snpEff databases;
	echo # -g/--snpeff_hgvs                        snpEff HGVS annotation in INFO:hgvs;
	echo # -z/--snpeff_stats                       snpEff stats;
	echo # -j/--java                               JAVA binary;
	echo # -p/--annovar_annotation_type_default    In case of ANNOVAR code in annotation;
	echo # -q/--AlleleFrequency                    AlleleFrequency (Not available);
	echo # -u/--annotations                        Show annotation;
"

: "		# Main options 
	'help'		=>  	0,	# Help parameter
	'man'		=>  	0,	# Man parameter
	'release'	=>  	0,	# Release parameter
	'debug'		=>  	0,	# Debug parameter
	'verbose'	=>  	0,	# Verbose parameter
	# Configuration
	'config'		=>	'config.ini',			# Configuration file
	'config_annotation'    	=>  	"config.annotation.ini",	# Configuratin annotation file
	'assembly'    		=>  	"",				# Assembly
	'annovar_folder'    	=>  	"",				# ANNOVAR folder including scripts
	'annovar_databases'    	=>  	"",				# ANNOVAR databases
	'snpeff'    		=>  	0,				# snpEff annotation
	'snpeff_jar'    	=>  	"",				# snpEff JAR
	'snpeff_databases'    	=>  	"",				# snpEff databases
	'snpeff_hgvs'    	=>  	0,				# snpEff HGVS annotation in INFO:hgvs
	'snpeff_stats'    	=>  	"",				# snpEff stats
	'java'    		=>  	"java",				# JAVA binary
	# Input
	'input'					=>	undef,		# Input VCF file
	'annotation'				=>	'ALL',		# Annotation sources
	'annovar_annotation_type_default'	=>	'geneanno',	# In case of ANNOVAR code in annotation
	'AlleleFrequency'			=>	0,		# Annotation sources
	# Output
	'output'	=>	undef,		# Output VCF file
	# show
	'annotations'	=>	0,
	
	'help|h|?',		# Help
	'man',			# Man
	'release',		# Release
	'debug',		# Debug
	'verbose',		# Verbose
	# Configuration
	'config|config_file=s',				# Configuration file
	'config_annotation|config_annotation_file=s',	# Configuratin annotation file
	'assembly=s',					# Assembly
	'annovar_folder=s',				# ANNOVAR folder including scripts
	'annovar_databases=s',				# ANNOVAR databases
	'snpeff!',					# snpEff annotation
	'snpeff_jar=s',					# snpEff JAR
	'snpeff_databases=s',				# snpEff databases
	'snpeff_hgvs!',					# snpEff HGVS annotation in INFO:hgvs
	'snpeff_stats=s',				# snpEff stats files
	'java=s',					# JAVA binary
	# Input
	'input|input_file=s',			# Input file
	'annotation=s',				# Annotations
	'annovar_annotation_type_default=s',	# In case of ANNOVAR code in annotation
	'AlleleFrequency!',			# Allele Frequency Calculation
	# output
	'output|output_file=s',	# Output file
	# show
	'annotations!',		# Show annotation
"

}

# header
header;

#ARGS=$(getopt -o "i:o:a:e:f:y:r:b:sm:l:gz:j:p:quxt:vdnh" --long "input:,output:,annotation:,config:,config_annotation:,assembly:,annovar_folder:,annovar_databases:,snpeff,snpeff_jar:,snpeff_database:,snpeff_hgvs,snpeff_stats:,java:,annovar_annotation_type_default:,AlleleFrequency,annotations,multithreading,threads:,verbose,debug,release,help" -- "$@" 2> /dev/null)
ARGS=$(getopt -o "i:o:a:xt:vdnh" --long "input:,output:,annotation:,multithreading,threads:,verbose,debug,release,help" -- "$@" 2> /dev/null)
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
		-i|--input)
			if [ ! -e $2 ] || [ "$2" == "" ]; then
				echo "#[ERROR] No VCF file '$2'"
				usage; exit;	
			else 
				INPUT="$2";
			fi;
			shift 2 
			;;
		-o|--output)
			OUTPUT="$2";
			shift 2 
			;;
		-a|--annotation)
			ANNOTATION="$2"
			shift 2 
			;;
		-e|--config)
			CONFIG="$2"
			shift 2 
			;;
		-f|--config_annotation)
			config_annotation="$2"
			shift 2 
			;;
		-y|--assembly)
			assembly="$2"
			shift 2 
			;;
		-r|--annovar_folder)
			annovar_folder="$2"
			shift 2 
			;;
		-b|--annovar_databases)
			annovar_databases="$2"
			shift 2 
			;;
		-s|--snpeff)
			SNPEFF=1
			shift 1
			;;
		-m|--snpeff_jar)
			snpeff_jar="$2"
			shift 2 
			;;
		-l|--snpeff_database)
			snpeff_database="$2"
			shift 2 
			;;
		-g|--snpeff_hgvs)
			snpeff_hgvs=1
			shift 1
			;;
		-z|--snpeff_stats)
			snpeff_stats=1
			shift 1
			;;
		-j|--java)
			java="$2"
			shift 2
			;;
		-p|--annovar_annotation_type_default)
			annovar_annotation_type_default="$2"
			shift 2
			;;
		-q|--AlleleFrequency)
			AlleleFrequency="$2"
			shift 2
			;;
		-u|--annotations)
			annotations="$2"
			shift 2
			;;
		-x|--multithreading)
			MULTITHREADING=1
			shift 1
			;;
		-t|--threads)
			THREADS_INPUT="$2"
			shift 2
			;;
		-v|--verbose)
			VERBOSE=1
			shift 1
			;;
		-d|--debug)
			VERBOSE=1
			DEBUG=1
			shift 1
			;;
		-n|--release)
			release;
			exit 0
			;;
		-h|--help)
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





# Mandatory parameters
if [ -z $INPUT ]; then #  || [ -z $OUTPUT ]; then
	echo "# INPUT Missing";
	usage;
	exit;
fi

# Mandatory parameters
if [ -z $OUTPUT ] || [ "$OUTPUT" == "" ]; then #  || [ -z $OUTPUT ]; then
	echo "# OUTPUT Missing";
	OUTPUT=$(echo $INPUT | sed "s/.vcf$/.annotated.vcf/g"); #/\.vcf$/\.output\.vcf/;
fi


# Multithreading
if [ -z $MULTITHREADING ]; then
	MULTITHREADING=0;
fi;

# Threads
if [ -z $THREADS_INPUT ]; then
	THREADS=$THREADS_INPUT;
fi;
if [ -z $THREADS ]; then
	THREADS=1;
fi;



#echo "TMP_SYS_FOLDER=$TMP_SYS_FOLDER"; exit 0;
# TMP_FILES
TMP_FILES=""
TMP_VAR=$RANDOM$RANDOM
TMP_FOLDER=/tmp/HOWARD_$TMP_VAR
mkdir $TMP_FOLDER

# input
echo "
##################
# INPUT=$INPUT
# OUTPUT=$OUTPUT
# MULTITHREADING=$MULTITHREADING ($threads threads)"
(($VERBOSE)) && echo "# CONFIG=$CONFIG";
(($VERBOSE)) && echo "";
#(($VERBOSE)) && echo "##################"
#(($VERBOSE)) && echo "";


# Cleaning PARAM
PARAM=$(echo "$PARAM " | sed "s/--threads=.*[ |$]//g" | sed "s/-t=.*[ |$]//g" | sed "s/--threads .*[ |$]//g" | sed "s/-t .*[ |$]//g")
PARAM=$(echo "$PARAM " | sed "s/--multithreading[ |$]//g" | sed "s/-x[ |$]//g")
#echo $PARAM; exit 0;

MK=HOWARD.$TMP_VAR.mk

if (($MULTITHREADING)); then
	#echo "#[WARNING] Multithreading not yet/fully implemented...";
	echo "ANNOTATION=$ANNOTATION"
	# FIND SAMPLE
	#SAMPLE=$(grep "^#CHROM" $INPUT | cut -f10)
	SAMPLES=$(grep "#CHROM" $INPUT | cut -f10-$(grep "#CHROM" $INPUT | wc -w) --output-delimiter=,)
	echo "SAMPLES=$SAMPLES"
	if [ "$SAMPLES " != "" ]; then
		echo "all: $OUTPUT" >>$MK;
		for ANN in $(echo $ANNOTATION | tr "," " "); do
			echo $ANN;
			echo "$SCRIPT_DIR/VCFannotation.pl $PARAM --annotation=$ANN"
			#$SCRIPT_DIR/VCFannotation.pl $PARAM --annotation=$ANN --output=$INPUT.TRUC.vcf
			VCFGZ_ANN=$TMP_FOLDER/$ANN.vcf.gz
			echo "$VCFGZ_ANN: $INPUT
				$SCRIPT_DIR/VCFannotation.pl $PARAM --annotation=$ANN --output=\$@.tmp
				$BGZIP -c \$@.tmp > \$@
				$TABIX \$@
			" >>$MK;
			VCFGZ_LIST=$VCFGZ_LIST" $VCFGZ_ANN"
		done;
		echo "$OUTPUT: $VCFGZ_LIST
			$VCFTOOLS/vcf-merge $VCFGZ_LIST > \$@.merge.vcf
			$VCFTOOLS/vcf-subset -c $SAMPLES \$@.merge.vcf > \$@
			#rm $OUTPUT.merge.vcf
		" >>$MK;
		echo -e "\n\n"
		cat $MK
		make -j $THREADS -f $MK
		(($VERBOSE)) && tail -n 10 $OUTPUT.merge.vcf
		(($VERBOSE)) && tail -n 10 $OUTPUT
		rm -f$MK
		
	else
		echo "#[ERROR] no Smaple name found in input VCF file"
	fi;
else
	$SCRIPT_DIR/VCFannotation.pl $PARAM
fi;

rm -Rf $TMP_FOLDER


exit 0;






