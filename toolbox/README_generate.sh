#!/bin/bash
##############################
## HOWARD README Generation ##
##############################

SCRIPT_NAME="HOWARD Validation"
SCRIPT_DESCRIPTION="README generation of HOWARD tool"
SCRIPT_RELEASE="0.9.3.1b"
SCRIPT_DATE="24/01/2019"
SCRIPT_AUTHOR="Antony Le Bechec"
SCRIPT_COPYRIGHT="HUS"
SCRIPT_LICENCE="GNU AGPL V3"

# Realse note
RELEASE_NOTES=$RELEASE_NOTES"# 0.9b-06/04/2018:\n";
RELEASE_NOTES=$RELEASE_NOTES"#\tScript creation\n";
RELEASE_NOTES=$RELEASE_NOTES"# 0.9.3.1b-24/01/2019:\n";
RELEASE_NOTES=$RELEASE_NOTES"#\tScript generalisation\n";

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
	echo "# USAGE: $(basename $0) [options...]";
	echo "# Following options are available:";
	echo "# --verbose                  VERBOSE option";
	echo "# --debug                    DEBUG option";
	echo "# --release                  RELEASE option";
	echo "# --help                     HELP option";
	echo "#";
	echo "";

}


# header
header;


ARGS=$(getopt -o "i:o:e:m:vdnh" --long "verbose,debug,release,help" -- "$@" 2> /dev/null)
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

README=$SCRIPT_DIR/../docs/README

# README Generation
cat $SCRIPT_DIR/../docs/HEADER > $README
echo "" >> $README
cat $SCRIPT_DIR/../docs/REQUIREMENTS >> $README
echo "" >> $README
cat $SCRIPT_DIR/../docs/INSTALLATION >> $README
echo "" >> $README
$SCRIPT_DIR/../bin/HOWARD.sh --release >> $README
$SCRIPT_DIR/../bin/HOWARD.sh --help >> $README

$SCRIPT_DIR/../bin/VCFannotation.pl --release  >> $README
$SCRIPT_DIR/../bin/VCFannotation.pl --help  >> $README
$SCRIPT_DIR/../bin/VCFcalculation.pl --release  >> $README
$SCRIPT_DIR/../bin/VCFcalculation.pl --help  >> $README
$SCRIPT_DIR/../bin/VCFprioritization.pl --release  >> $README
$SCRIPT_DIR/../bin/VCFprioritization.pl --help  >> $README
$SCRIPT_DIR/../bin/VCFtranslation.pl --release  >> $README
$SCRIPT_DIR/../bin/VCFtranslation.pl --help  >> $README
$SCRIPT_DIR/../bin/VCFtranslation.sh --release  >> $README
$SCRIPT_DIR/../bin/VCFtranslation.sh --help  >> $README
