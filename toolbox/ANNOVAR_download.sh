#!/bin/bash
#################################
##
## ANNOVAR DB Download
##
## version: 0.9b
## date: 02/06/2015
## author: Antony Le Bechec
##
#################################

SCRIPT_NAME="ADaBaDoLa"
SCRIPT_DESCRIPTION="ANNOVAR DB Download"
SCRIPT_RELEASE="0.9"
SCRIPT_DATE="03/06/2015"
SCRIPT_AUTHOR="Antony Le Bechec"
SCRIPT_COPYRIGHT="IRC"

echo "#######################################";
echo "# $SCRIPT_NAME [$SCRIPT_RELEASE-$SCRIPT_DATE]";
echo "# $SCRIPT_DESCRIPTION ";
echo "# $SCRIPT_AUTHOR @ $SCRIPT_COPYRIGHT © ";
echo "#######################################";

# Realse note
RELEASE_NOTES="# Release Notes\n"
RELEASE_NOTES=$RELEASE_NOTES"02/06/2015-0.9b: Create script from dev script.\n"
# 
# Example1: SELF "" "cosmic70" "" "" ""
# Example2: SELF "hg19" "cosmic70" "humandb/" "-webfrom annovar" "annotate_variation.pl"



######################################
# 1. Configuration, Input parameters #
######################################

# 1.1. Configuration
#####################

#USAGE/HELP
USAGE="# USAGE: "`basename $0`" <BUILDVER> <DB_LIST> <DB_FOLDER> <WEBFROM> <ANNOVAR_SCRIPT> <PARALLELIZATION>\n";
EXAMPLE="# EXAMPLE: "`basename $0`" hg19 'Symbol COSMIC CLINVAR' '/media/IRCV2/NGSEnv/annovar_sources/' ' -webfrom annovar' '/media/IRCV2/NGSEnv/tools2/annovar/current/bin/annotate_variation.pl' ''\n";
HELP=$USAGE$EXAMPLE"# BUILDVER: Build/Assembly (default: 'hg19')\n"
HELP=$HELP"# DB_LIST: List of ANNOVAR Table Name (default: 'DEFAULT')\n"
HELP=$HELP"# DB_FOLDER: Database folder (default: 'humandb/')\n"
HELP=$HELP"# WEBFROM: Webfrom option (default: '-webfrom annovar', use webfrom)\n"
HELP=$HELP"# ANNOVAR_SCRIPT: ANNOVAR script location (default: 'annotate_variation.pl')\n"
HELP=$HELP"# PARALLELIZATION: Parallelization of the download (default: "", no parallelisation, or whatever to parallelize)\n"

if [ "${1^^}" == "HELP" ] || [ "${1^^}" == "USAGE" ]; then
	echo -e $HELP;
	exit 0;
fi;



# PARAMETERS
##############

# DEFAULT DB LIST
DB_LIST_DEFAULT=" refGene knownGene ensGene "

# FULL DB LIST
DB_LIST_FULL=$DB_LIST_DEFAULT # DEFAULT
#DB_LIST_FULL=$DB_LIST_FULL" ljb26_all cg46 cg69" # SCORES
DB_LIST_FULL=$DB_LIST_FULL" ljb26_sift ljb26_pp2hdiv ljb26_pp2hvar ljb26_cadd ljb26_phylop46way_placental ljb26_phylop100way_vertebrate ljb26_lrt ljb26_mt ljb26_ma ljb26_fathmm ljb26_siphy ljb26_gerp++ ljb26_metasvm ljb26_metalr ljb26_vest" # SCORES
DB_LIST_FULL=$DB_LIST_FULL" cg46 cg69" # SCORES
DB_LIST_FULL=$DB_LIST_FULL" cosmic68 cosmic68wgs cosmic70 clinvar_20140929 " # CLINICAL KNOWLEDGE
DB_LIST_FULL=$DB_LIST_FULL" esp6500siv2_all esp6500siv2_ea esp6500siv2_aa exac02 1000g2014oct popfreq_max" # POPULATIONS
DB_LIST_FULL=$DB_LIST_FULL" snp130 snp135 snp137 snp138 snp130NonFlagged snp135NonFlagged snp137NonFlagged snp138NonFlagged" # DBSNP


# BUILDVER
BUILDVER=$1;
if [ "$BUILDVER" == "" ]; then
	BUILDVER="hg19";
fi;

# DB_LIST
DB_LIST=$2
if [ "$DB_LIST" == "" ]; then
	DB_LIST="DEFAULT";
fi;
if [ "$DB_LIST" == "DEFAULT" ]; then
	DB_LIST=$DB_LIST_DEFAULT;
	echo "# DEFAULT databases will be downloaded ";
fi;
if [ "$DB_LIST" == "FULL" ] || [ "$DB_LIST" == "ALL" ]; then
	DB_LIST=$DB_LIST_FULL;
	echo "# FULL databases will be downloaded ";
fi;

# DB_FOLDER
DB_FOLDER=$3;
if [ "$DB_FOLDER" == "" ]; then
	DB_FOLDER="humandb/";
fi;
if [ ! -d $DB_FOLDER  ]; then
	mkdir -p $DB_FOLDER;
fi;
if [ ! -d $DB_FOLDER ]; then
	echo "# [ERROR] DB_FOLDER '$DB_FOLDER' doesn't exists!";
	echo "# USAGE: $USAGE"
	exit 1;
fi;


# WEBFROM
WEBFROM=$4;
if [ "$WEBFROM" == "" ]; then
	WEBFROM="-webfrom annovar";
fi;

# ANNOVAR_SCRIPT
ANNOVAR_SCRIPT=$5;
if [ "$ANNOVAR_SCRIPT" == "" ]; then
	ANNOVAR_SCRIPT="annotate_variation.pl";
fi;

# PARALLELIZATION
PARALLELIZATION=$6;
if [ "$PARALLELIZATION" == "" ]; then
	PARALLELIZATION="";
else
	PARALLELIZATION="&";	
fi;


# DOWNLOAD
############

echo "# BUILDVER : $BUILDVER ";
echo "# DATABASES to download:";
echo "# $DB_LIST ";

for BUILD in $BUILDVER;
do
	for DB in $DB_LIST;
	do
		echo "# DOWNLOAD $DB";
		if [ "$PARALLELIZATION" == "" ]; then
			TMP=$BUILD.$DB.log;
			COMMAND="$ANNOVAR_SCRIPT -downdb -buildver $BUILD $DB $DB_FOLDER $WEBFROM 1> $TMP 2>$TMP";
			echo "# command: $COMMAND";
			eval $COMMAND;
			if [ `grep -c "WARNING" $TMP` == "0" ]; then
				echo "# DOWNLOAD $DB OK";
			else
				echo "# DOWNLOAD $DB Failed";				
				echo "# "`grep "WARNING" $TMP`;
			fi;
		else
			COMMAND="$ANNOVAR_SCRIPT -downdb -buildver $BUILD $DB $DB_FOLDER $WEBFROM $PARALLELIZATION ";
			echo "# command: $COMMAND";
			eval $COMMAND;
		fi;
	done;
done;

echo "# DONE."


exit 0;
