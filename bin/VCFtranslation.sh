#!/bin/bash
#################################
## HOWARD
#################################

SCRIPT_NAME="HOWARD"
SCRIPT_DESCRIPTION="HOWARD Translation with awk"
SCRIPT_RELEASE="0.9b"
SCRIPT_DATE="21/01/2019"
SCRIPT_AUTHOR="Antony Le Bechec"
SCRIPT_COPYRIGHT="HUS"
SCRIPT_LICENCE="GNU AGPL V3"

# Realse note
RELEASE_NOTES=$RELEASE_NOTES"# 0.9b-21/01/2019:\n";
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
	echo "# --input|vcf=<FILES>                Input file in VCF format (*vcf or *vcf.gz in BGZIP compression format)";
	echo "# --output=<FILE>                    Output annotated file in defined format (see --format option). Default '<input>.translated.<format>'";
	echo "# --translation=<STRING>             Translation: Output format, either TSV or VCF (default VCF)";
	echo "# --fields=<STRING>                  List of annotations from INFO VCF field to include in the output file.";
	echo "#                                    Use 'INFO' to insert a uniq field INFO in TSV format";
	echo "#                                    Use 'ALL' to insert all other annotations in the list";
	echo "#                                    Annotations considered only if present in the file";
	echo "#                                    Default 'ALL'";
	echo "#                                    Example: 'PZScore,PZFlag,PZComment,Symbol,hgvs,location,outcome,ALL'";
	echo "# --sort=<STRING>                    Sort variants by a field and order (default '')";
	echo "#                                    Format: 'field1:type:order,field2:type:order'";
	echo "#                                       field considered only if present in the file";
	echo "#                                       type 'n' for numeric, '?' for autodetect with VCF header. Default ''.";
  	echo "#                                       order either 'DESC' or 'ASC'. Default 'ASC'";
  	echo "#                                    Example: 'PZFlag::DESC,PZScore:n:DESC' (to sort and order by relevance)";
  	echo "# --sort_by=<STRING>                 Sort variants by a field (if no 'sort' option). Default ''";
	echo "#                                    Example: 'PZFlag,PZScore' (to sort by relevance)";
	echo "# --order_by=<STRING>                Order variants by a field (if no 'sort' option). Default ''.";
	echo "#                                    Example: 'DESC,DESC' (useful to sort by relevance)";
	echo "# --bcftools_expression=<STRING>     bcftools include expression to filter variants (default '').";
	echo "#                                    Example: \"PZFlag='PASS'\" (useful to hard filter), \"QUAL>90\", or \"PZFlag='PASS'&&QUAL>90\" to combine ";
	echo "# --env=<FILE>                       Environment configuration for multithreading (BGZIP, TABIX, BCFTOOLS)";
	#echo "# --compress=<INTEGER>               Compression level output file *vcf.gz (0 to 9, -1 no compression by default)";
	echo "# --force                            Force annotation even if already exists in VCF header";
	echo "# --tmp=<FOLDER>                     Temporary folder (default /tmp)";
	echo "# --verbose                          VERBOSE option";
	echo "# --debug                            DEBUG option";
	echo "# --release                          RELEASE option";
	echo "# --help                             HELP option";
	echo "";

}

# EXAMPLE :
# ./VCFtranslation.sh --input=example.vcf --output=example.tsv --sort="PZFlag:?:DESC,PZScore:n:DESC" --fields="PZFlag,PZScore,NOMEN,Symbol,ALL" --translation=TSV  --bcftools_expression="QUAL>90"
# header
header;


ARGS=$(getopt -o "i:o:e:a:f:s:r:xt:m:vdnh" --long "input:,vcf:,output:,env:,translation:,format:,fields:,sort:,sort_by:,order_by:,bcftools_expression:,tmp:,verbose,debug,release,help" -- "$@" 2> /dev/null)
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
		--translation|--format)
			TRANSLATION="$2"
			shift 2
			;;
		--fields)
			FIELDS="$2"
			shift 2
			;;
		--sort)
			SORT="$2"
			shift 2
			;;
		--sort_by)
			SORT_BY="$2"
			shift 2
			;;
		--order_by)
			ORDER_BY="$2"
			shift 2
			;;
		--bcftools_expression)
			BCFTOOLS_EXPRESSION="$2"
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


# TRANSLATION/FORMAT
# Mandatory parameters
if [ -z "$TRANSLATION" ]; then #  || [ -z $OUTPUT ]; then
	echo "#[INFO] FORMAT Missing. Default FORMAT 'VCF'";
	TRANSLATION="VCF"
fi

FORMAT=$(echo "$TRANSLATION" | tr '[:lower:]' '[:upper:]')

# Mandatory parameters
if [ -z $OUTPUT ] || [ "$OUTPUT" == "" ]; then #  || [ -z $OUTPUT ]; then
	OUTPUT=$INPUT".translated."$(echo "$TRANSLATION" | tr '[:upper:]' '[:lower:]'); #/\.vcf$/\.output\.vcf/;
	echo "#[INFO] OUTPUT Missing. Default OUTPUT '$OUTPUT' used."
fi
# mkdir OUPUT directory
mkdir -p $(dirname $OUTPUT)


if [ -z "$FIELDS" ]; then #  || [ -z $OUTPUT ]; then
	echo "#[INFO] FIELDS Missing. Default FIELDS 'ALL'";
	FIELDS="ALL"
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


# TMPORARY FOLDER and LOG/ERR
################################

TMP_FILES=""
TMP_VAR=$RANDOM$RANDOM
TMP_FOLDER=$TMP_SYS_FOLDER/HOWARD_$TMP_VAR
LOG=$TMP_FOLDER/log
ERR=$TMP_FOLDER/err
mkdir -p $TMP_FOLDER
touch $LOG;
touch $ERR;





if [ "$SORT" != "" ]; then
	ANNS=$SORT
else
	echo $SORT_BY
	echo $ORDER_BY

	echo $SORT_BY | tr ',' ' ' | tr ' ' '\n' | sed '/^$/d' > $TMP_FOLDER/sort_by
	echo $ORDER_BY | tr ',' ' ' | tr ' ' '\n' | sed '/^$/d' > $TMP_FOLDER/order_by

	ANNS=$(paste $TMP_FOLDER/sort_by $TMP_FOLDER/order_by | awk '{print $1":?:"$2}')

fi;

LIMIT_VARIANT=
#BCFTOOLS_INCLUDE_EXPRESSION="PZFlag='PASS'"
BCFTOOLS_INCLUDE_EXPRESSION=$BCFTOOLS_EXPRESSION
#BCFTOOLS_EXCLUDE_EXPRESSION=""
VCFTRANSLATION_AWK=$SCRIPT_DIR/VCFtranslation.awk

# BUG FIX
#[ $BCFTOOLS_EXPRESSION ]

INFINIT=100000000000000000000
if [ -z $LIMIT_VARIANT ]; then LIMIT_VARIANT=$INFINIT; fi;
if [ "$FORMAT" == "" ]; then FORMAT="VCF"; fi;
if [ "$FORMAT" == "TAB" ]; then FORMAT="TSV"; fi;
if [ "$FIELDS" == "" ]; then FIELDS="INFO"; fi;


FORMAT_OUTPUT=$(echo "$FORMAT" | tr '[:upper:]' '[:lower:]')


## TEST BCFTOOLS
#BCFTOOLS=bcftools



	ANN=$ANNS

	# output file
	FA=$TMP_FOLDER/INPUT.vcf
	#echo "FA="$FA

	# Generate sort parameters
	nb_ANN=$(echo $ANN | tr "," " " | wc -w)
	sort_param=""
	nb=0
	sep=""
	ANN_clear_list=""
	for a in $(echo $ANN | tr "," " "); do
		((nb++))
		order=""
		ANN_clear=$(echo $a | awk -F: '{print $1}')
		ANN_clear_list=$ANN_clear_list$sep$ANN_clear
		sep=","

		if [ "$(echo $a | awk -F: '{print $3}' | tr '[:upper:]' '[:lower:]')" == "desc" ]; then
			order="r"
		fi;
		comparaison="V"
		if [ "$(echo $a | awk -F: '{print $2}' | tr '[:upper:]' '[:lower:]')" == "n" ]; then
			comparaison="n"
		elif [ "$(echo $a | awk -F: '{print $2}')" == "?" ]; then
			if [ "$(grep "ID=$ANN_clear," $INPUT -m1 | awk -F"Type=" '{print $2}' | awk -F"," '{print $1}')" == "Integer" ]; then
				comparaison="n"
			fi;
		fi;

		#grep ID=PZScore 1G.1000.vcf | awk -F"Type=" '{print $2}' | awk -F"," '{print $1}'
		sort_param=$sort_param" -k"$nb","$nb$comparaison$order" "
	done
	((nb_ANN++))

	if [ "$sort_param" == "" ]; then
		sort_param=" -k1,1V -k2,2n"
	fi;


	# Reduce VCF file with fields
	###############

	#echo ""
	#echo "# Reduce File"

	echo $FIELDS | tr "," "\n" > $FA.tsv.listINFO.param
	$BCFTOOLS view -h $INPUT | awk -F"##INFO=<ID=" '{print $2}' | awk -F"[,>]" '{print $1}' | sort -u > $FA.tsv.listINFO.ALL
	if (($(grep "^ALL$" $FA.tsv.listINFO.param -c))); then
		echo "ALL" >> $FA.tsv.listINFO.ALL
	fi;
	grep -Fxv -f $FA.tsv.listINFO.ALL $FA.tsv.listINFO.param | grep -Fxv -f - $FA.tsv.listINFO.param > $FA.tsv.listINFO.param.clean #| sort -u
	cat $FA.tsv.listINFO.param.clean > $FA.tsv.listINFO

	# IN all but not in listINFO
	grep -Fxv -f $FA.tsv.listINFO $FA.tsv.listINFO.ALL > $FA.tsv.listINFO.ALL.REST


	if (($(grep "^ALL$" $FA.tsv.listINFO.param -c))); then

		> $FA.tsv.listINFO.tmp
		for FIELD in $(cat $FA.tsv.listINFO); do
			if [ "$FIELD" == "ALL" ]; then
				cat $FA.tsv.listINFO.ALL.REST >> $FA.tsv.listINFO.tmp
			else
				echo $FIELD >> $FA.tsv.listINFO.tmp
			fi;
		done
		cp $FA.tsv.listINFO.tmp $FA.tsv.listINFO
		grep "^ALL$" -v $FA.tsv.listINFO.ALL > $FA.tsv.listINFO.ALL.tmp
		cp $FA.tsv.listINFO.ALL.tmp $FA.tsv.listINFO.ALL
	fi;

	(($DEBUG)) && echo "#[INFO] Selected Fields: "$(cat $FA.tsv.listINFO)

	grep -Fxv -f $FA.tsv.listINFO $FA.tsv.listINFO.ALL | sort -u > $FA.tsv.listINFO.inv #| grep -Fxv -f - $FA.tsv.listINFO.ALL



	# FINAL
	##########

	#echo ""
	#echo "# Generate output"



	FINAL=$OUTPUT #$i.final.$FORMAT_OUTPUT

	#if [ "$FORMAT" == "TSV" ]; then

	# translate to TSV

	# List of ORMAT TAG
	$BCFTOOLS view -h $INPUT | awk -F"##FORMAT=<ID=" '{print $2}' | awk -F"," '{print $1}' | sort -u > $FA.tsv.listFORMAT.ALL
	#cat $FA.tsv.listFORMAT.ALL

	LIST_VARIANTS_IDS="%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL\t%FILTER"

	LIST_FORMAT=""
	sep=""
	for VCFFORMAT in $(cat $FA.tsv.listFORMAT.ALL); do
		if [ "$VCFFORMAT" != "" ]; then
			#echo $INFO
			LIST_FORMAT=$LIST_FORMAT$sep"$VCFFORMAT=%$VCFFORMAT"
		fi;
		sep=":"
	done;
	#echo $LIST_FORMAT
	LIST_INFOS_TOTAL_FORMAT="\t"$(echo -ne $LIST_FORMAT | tr -d "%")
	#echo $LIST_INFOS_TOTAL_FORMAT
	LIST_INFOS_TOTAL_SAMPLES="[\t$LIST_FORMAT]"

	# List of fields to export by order
	LIST_INFOS=""
	sep=""
	for INFO in $(cat $FA.tsv.listINFO); do
		if [ "$INFO" != "" ]; then
			#echo $INFO
			LIST_INFOS=$LIST_INFOS$sep"%INFO/$INFO"
		fi;
		sep="\t"
	done;
	if [ -z $LIST_INFOS ]; then LIST_INFOS="INFO"; fi;

	LIST_INFOS_TOTAL="$LIST_VARIANTS_IDS\t$LIST_INFOS$LIST_INFOS_TOTAL_SAMPLES"

	# List of samples
	#LIST_SAMPLES=$(bcftools query -l $i | tr "\n" "\t"  | sed s/\t$//g)
	LIST_SAMPLES=$($BCFTOOLS view -h $INPUT  | grep ^#CHROM | cut -f10-)
	#echo $LIST_SAMPLES
	#echo $LIST_INFOS


	#FORMAT="VCF"

	> $FINAL
	if [ "$FORMAT" == "TSV" ]; then
		echo -ne "#" > $FINAL
		echo -ne $LIST_INFOS_TOTAL | sed "s/INFO\///g" | tr -d "%"  | sed 's/\[.*\]//gi' >> $FINAL
		#echo -ne $ANN_clear_list | tr "," "\t" | sed "s/INFO\///g" | tr -d "%"  | sed 's/\[.*\]//gi' >> $FINAL
		#cat $FA.tsv.listINFO | tr "\n" "\t" | sed s/\t$//g  >> $FINAL

		echo -ne "\tFORMAT" >> $FINAL
		#echo -ne $LIST_FORMAT | tr -d "%" >> $FINAL

		echo -ne "\t" >> $FINAL
		#echo $LIST_SAMPLES >> $FINAL
		$BCFTOOLS view -h $INPUT  | grep ^#CHROM | cut -f10- >> $FINAL
	else
		$BCFTOOLS view -h $INPUT >> $FINAL
	fi;

	#echo "LIST_INFOS_TOTAL=$LIST_INFOS_TOTAL"
	#echo "ANN_clear_list=$ANN_clear_list"
	#cat $FA.tsv.listINFO
#exit
#-v ANN=$ANN_clear_list
#-v FIELDS=$(cat $FA.tsv.listINFO | tr "\n" "," | sed s/,$//g)
	FIELDS_AWK=$(cat $FA.tsv.listINFO | tr "\n" "," | sed s/,$//g)

	if [ -z $nb_ANN ]; then $nb_ANN=1; fi;
	if [ -z $FIELDS_AWK ]; then FIELDS_AWK="INFO"; fi;


	(($DEBUG)) && echo "nb_ANN=$nb_ANN"
	#grep ^# -v $i | awk -F"\t"  -v FORMAT=$FORMAT -v ANN=$ANN_clear_list -v FIELDS=$FIELDS_AWK -f $VCFTRANSLATION_AWK | sort $sort_param | cut -f$nb_ANN- | head -n $LIMIT_VARIANT >> $FINAL
	#echo "grep ^# -v $i | awk -F"\t"  -v FORMAT=$FORMAT -v ANN=$ANN_clear_list -v FIELDS=$FIELDS_AWK -f $VCFTRANSLATION_AWK | sort $sort_param | cut -f$nb_ANN- | head -n $LIMIT_VARIANT >> $FINAL"
	(($DEBUG)) && echo "ANN_clear_list=$ANN_clear_list"
	#echo "FIELDS=$FIELDS_AWK"
	#grep ^# -v $i | awk -F"\t" -v FORMAT=TSV -v ANN=$ANN_clear_list -v FIELDS= -f $VCFTRANSLATION_AWK | sort $sort_param | cut -f$nb_ANN-  | head -n3
	#exit
	#grep ^# -v $i | awk -F"\t" -v FORMAT=VCF -v ANN=$ANN_clear_list -v FIELDS="INFO" -f $VCFTRANSLATION_AWK
	#exit



	# OUTPUT HEADER
	##################

	echo "#[INFO] INPUT=$INPUT";
	echo "#[INFO] OUTPUT=$OUTPUT";
	echo "#[INFO] FORMAT=$TRANSLATION";
	echo "#[INFO] FIELDS=$FIELDS_AWK";
	echo "#[INFO] SORT=$SORT";
	echo "#[INFO] BCFTOOLS_EXPRESSION=$BCFTOOLS_EXPRESSION";
	echo "#[INFO] TMP=$TMP_FOLDER";

	(($DEBUG)) && echo "
	FILES=$INPUT
	ANNS=$ANNS
	FIELDS=$FIELDS
	FORMAT=$FORMAT
	FORMAT_OUTPUT=$FORMAT_OUTPUT
	LIMIT_VARIANT=$LIMIT_VARIANT
	BCFTOOLS_INCLUDE_EXPRESSION=$BCFTOOLS_INCLUDE_EXPRESSION
	"



		if [ "$BCFTOOLS_INCLUDE_EXPRESSION" != "" ]; then
			BCFTOOLS_INCLUDE_EXPRESSION_param="-i \"$BCFTOOLS_INCLUDE_EXPRESSION\""
		fi;


		if [ $LIMIT_VARIANT == $INFINIT ]; then
			CMD="$BCFTOOLS view -H $INPUT $BCFTOOLS_INCLUDE_EXPRESSION_param | awk -F'\t'  -v FORMAT=$FORMAT -v ANN=$ANN_clear_list -v FIELDS=$FIELDS_AWK -f $VCFTRANSLATION_AWK | sort $sort_param | cut -f$nb_ANN- >> $FINAL;"
			#$BCFTOOLS view -H $INPUT $BCFTOOLS_INCLUDE_EXPRESSION_param | awk -F"\t"  -v FORMAT=$FORMAT -v ANN=$ANN_clear_list -v FIELDS=$FIELDS_AWK -f $VCFTRANSLATION_AWK | sort $sort_param | cut -f$nb_ANN- >> $FINAL;
		else
			CMD="$BCFTOOLS view -H $INPUT $BCFTOOLS_INCLUDE_EXPRESSION_param | awk -F'\t' -v FORMAT=VCF -v ANN=$ANN_clear_list -v FIELDS='INFO' -f $VCFTRANSLATION_AWK | sort $sort_param | cut -f$nb_ANN- | head -n $LIMIT_VARIANT | awk -F'\t'  -v FORMAT=$FORMAT -v FIELDS=$FIELDS_AWK -f $VCFTRANSLATION_AWK >> $FINAL"
			#$BCFTOOLS view -H $INPUT $BCFTOOLS_INCLUDE_EXPRESSION_param  | awk -F"\t" -v FORMAT=VCF -v ANN=$ANN_clear_list -v FIELDS="INFO" -f $VCFTRANSLATION_AWK | sort $sort_param | cut -f$nb_ANN- | head -n $LIMIT_VARIANT | awk -F"\t"  -v FORMAT=$FORMAT -v FIELDS=$FIELDS_AWK -f $VCFTRANSLATION_AWK >> $FINAL
		fi;
		(($DEBUG)) && echo $CMD
		eval $CMD

		#exit

		# copy vcf to final
		#cp $FA.reduced.vcf $FINAL

		#fi;

		# INFOS
		#echo "# OUTPUT: $FINAL"
		#echo "# FORMAT: $FORMAT"
		(($DEBUG)) &&  echo "# OUTPUT"
		(($DEBUG)) &&  grep ^# $FINAL | head -n3 #| column -t
		(($DEBUG)) &&  grep ^# -v $FINAL | head -n3 #| column -t

		# Test output file
		if (($VERBOSE)); then
			nb_var_original=$(grep ^# -v $INPUT -c)
			nb_var_final=$(grep ^# -v $FINAL -c)
			echo "#[INFO] NB variants original: $nb_var_original"
			echo "#[INFO] NB variants final:    $nb_var_final"

			# Destructive test
			if (( $nb_var_original == nb_var_final )); then
				echo "#[INFO] NO destructive options"
			else
				echo "#[INFO] DESTRUCTIVE options"
			fi;
		fi;


exit 0;
