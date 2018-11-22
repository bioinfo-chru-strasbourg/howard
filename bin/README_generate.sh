#!/bin/bash
##############################
## HOWARD README Generation ##
##############################

# Script folder
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# README Generation
cat $SCRIPT_DIR/HEADER > $SCRIPT_DIR/README
cat $SCRIPT_DIR/REQUIREMENTS >> $SCRIPT_DIR/README
cat $SCRIPT_DIR/INSTALLATION >> $SCRIPT_DIR/README
$SCRIPT_DIR/HOWARD.sh --release >> $SCRIPT_DIR/README
$SCRIPT_DIR/HOWARD.sh --help >> $SCRIPT_DIR/README

$SCRIPT_DIR/VCFannotation.pl --release  >> $SCRIPT_DIR/README
$SCRIPT_DIR/VCFannotation.pl --help  >> $SCRIPT_DIR/README
$SCRIPT_DIR/VCFcalculation.pl --release  >> $SCRIPT_DIR/README
$SCRIPT_DIR/VCFcalculation.pl --help  >> $SCRIPT_DIR/README
$SCRIPT_DIR/VCFprioritization.pl --release  >> $SCRIPT_DIR/README
$SCRIPT_DIR/VCFprioritization.pl --help  >> $SCRIPT_DIR/README
$SCRIPT_DIR/VCFtranslation.pl --release  >> $SCRIPT_DIR/README
$SCRIPT_DIR/VCFtranslation.pl --help  >> $SCRIPT_DIR/README


