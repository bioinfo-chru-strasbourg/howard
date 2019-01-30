#!/usr/bin/perl
##############################
# Functions library          #
# Author: Antony Le Béchec   #
# Copyright: IRC             #
##############################

## Modules
############

#use DBI;


## Parameters
###############

if ($parameters{"help"}) {
	pod2usage(1);
};#if

## Man
if ($parameters{"man"}) {
	pod2usage(-exitstatus => 0, -verbose => 2);
};#if

## Release
if ($parameters{"release"}) {
	print "## Script Information\n";
	while ((my $var, $val) = each(%information)){
		print "# $var: $val\n";	
	};#while
	exit 0;
};#if

## Debug
our $DEBUG=1 if $parameters{"debug"};

## Verbose
our $VERBOSE=1 if $parameters{"verbose"};


## Configuration folder

our $config_folder=dirname($basename)."/config";
#warn "No configuration folder $config_folder...\n";

#print $config_folder."\n";
#print $parameters{"config"}."\n";
#print "$config_folder/".$parameters{"config"}."\n";
#print "$config_folder/config.ini"."\n";

## Configuration files

## MAIN CONFIGURATION

# Config file
our $config_file;
if (-e $parameters{"config"}) {
	$config_file=$parameters{"config"};
} elsif (-e "$config_folder/".$parameters{"config"}) {
	$config_file="$config_folder/".$parameters{"config"};
} elsif (-e "$config_folder/config.ini") {
	$config_file="$config_folder/config.ini";
#} elsif (-e "$basename/".$parameters{"config"}) {
#	$config_file="$basename/".$parameters{"config"};
#} elsif (-e "$basename/config.ini") {
#	$config_file="$basename/config.ini";
} else {
	warn "No Configuration file...\n";
	pod2usage(1);
};#if

# Default config

our $assembly="hg19";
our $annovar_databases="annovar";
our $annovar_folder=".";
our $snpeff_databases="snpeff";
our $snpeff_jar="";
our $java="java"; #$parameters{"java"};
our $java_flags=""; #$parameters{"java_flags"};
our $snpeff_threads=1;
our $tmp_folder="/tmp";


# Set parameters from config file
if (-e $config_file) {
	
	# read ini file
	our %config=read_ini($config_file);
	# Folders
	#$data_folder=$config{"folders"}{"data_folder"}."/";
	#$annovar_databases=$config{"folders"}{"annovar_databases"}."/";

	# ANNOVAR databases
	if (defined $config{"folders"}{"annovar_databases"} && -d $config{"folders"}{"annovar_databases"}) {
		$annovar_databases=$config{"folders"}{"annovar_databases"};  
	};#if
	if (trim($annovar_databases) ne "") {$annovar_databases.="/"};

	# ANNOVAR folder
	if (defined $config{"folders"}{"annovar_folder"} && -d $config{"folders"}{"annovar_folder"}) {
		$annovar_folder=$config{"folders"}{"annovar_folder"}; 
	};#if
	if (trim($annovar_folder) ne "") {$annovar_folder.="/"};
	
	# SNPEFF databases
	if (defined $config{"folders"}{"snpeff_databases"} && -d $config{"folders"}{"snpeff_databases"}) {
		$annovar_snpeff=$config{"folders"}{"snpeff_databases"};  
	};#if
	if (trim($annovar_snpeff) ne "") {$annovar_snpeff.="/"};
	
	# SNPEFF JAR
	if (defined $config{"folders"}{"snpeff_jar"} && -e $config{"folders"}{"snpeff_jar"}) {
		$snpeff_jar=$config{"folders"}{"snpeff_jar"}; 
	};#if
	
	# JAVA
	if (defined $config{"folders"}{"java"} && -e $config{"folders"}{"java"}) {
		$java=$config{"folders"}{"java"}; 
	};#if
	# JAVA FLAGS
	if (defined $config{"folders"}{"java_flags"} && -e $config{"folders"}{"java_flags"}) {
		$java_flags=$config{"folders"}{"java_flags"}; 
	};#if
	
	#tmp_folder
	if ((defined $config{"folders"}{"tmp_folder"} && -d $config{"folders"}{"tmp_folder"})) {
		$tmp_folder=$config{"folders"}{"tmp_folder"};    
	};#if

	#Project
	$assembly=$config{"project"}{"assembly"};

};#if

# Assembly in parameter
if ($parameters{"assembly"} ne "") {
	$assembly=$parameters{"assembly"};
};#if

# ANNOVAR databases
if ($parameters{"annovar_databases"} ne "" && -d $parameters{"annovar_databases"}) {
	$annovar_databases=$parameters{"annovar_databases"};
};#if

# ANNOVAR folder
if ($parameters{"annovar_folder"} ne "" && -d $parameters{"annovar_folder"}) {
	$annovar_folder=$parameters{"annovar_folder"};
};#if

# SNPEFF databases
if ($parameters{"snpeff_databases"} ne "" && -d $parameters{"snpeff_databases"}) {
	$snpeff_databases=$parameters{"snpeff_databases"};
};#if

# SNPEFF folder
if ($parameters{"snpeff_jar"} ne "" && -e $parameters{"snpeff_jar"}) {
	$snpeff_jar=$parameters{"snpeff_jar"};
};#if

# JAVA
if ($parameters{"java"} ne "") {
	$java=$parameters{"java"};
};#if

# JAVA FLAGS folder
if ($parameters{"java_flags"} ne "") {
	$java_flags=$parameters{"java_flags"};
};#if

# SNPEFF Threads
if ($parameters{"snpeff_threads"} > 0) {
	$snpeff_threads=$parameters{"snpeff_threads"};
};#if

# SNPEFF Threads
if ($parameters{"tmp_folder"} ne "" && -d $parameters{"tmp_folder"}) {
	$tmp_folder=$parameters{"tmp_folder"};
};#if


# TEST


# SNPEFF database
$snpeff_databases=$config{"folders"}{"snpeff_databases"};
if ($parameters{"snpeff_databases"} ne "" && -d $parameters{"snpeff_databases"}) {
	$snpeff_databases=$parameters{"snpeff_databases"};
} elsif (defined $config{"folders"}{"snpeff_databases"} && -d $config{"folders"}{"snpeff_databases"}) {
	$snpeff_databases=$config{"folders"}{"snpeff_databases"};
} else {
	$snpeff_databases=$basename;
};#if

# SNPEFF jar
$snpeff_jar=$config{"folders"}{"snpeff_jar"};

#print "$snpeff_jar" if $VERBOSE;
if ($parameters{"snpeff_jar"} ne "" && -e $parameters{"snpeff_jar"}) {
	$snpeff_jar=$parameters{"snpeff_jar"};
} elsif (defined $config{"folders"}{"snpeff_jar"} && -e $config{"folders"}{"snpeff_jar"}) {
	$snpeff_jar=$config{"folders"}{"snpeff_jar"};
} else {
	print "# snpEff JAR not found...\n";
	#pod2usage(1);
};#if
#print "snpeff_jar=$snpeff_jar";


## ANNOTATION CONFIGURATION

# Config annotation file
our $config_annotation_file;
if (-e $parameters{"config_annotation"}) {
    $config_annotation_file=$parameters{"config_annotation"};
} elsif (-e "$config_folder/".$parameters{"config_annotation"}) {
    $config_annotation_file="$config_folder/".$parameters{"config_annotation"};
} elsif (-e "$config_folder/config.annotation.ini") {
    $config_annotation_file="$config_folder/config.annotation.ini";
} else {
	print "#[WARNING] No Annotation Configuration file...\n";
	#pod2usage(1);
};#if

## FILTER CONFIGURATION

# Config filter file
our $config_filter_file;
if (-e $parameters{"config_prioritization"}) {
    $config_filter_file=$parameters{"config_prioritization"};
} elsif (-e "$config_folder/".$parameters{"config_prioritization"}) {
    $config_filter_file="$config_folder/".$parameters{"config_prioritization"};
} elsif (-e "$config_folder/config.prioritization.ini") {
    $config_filter_file="$config_folder/config.prioritization.ini";
} else {
	print "#[WARNING] No Prioritization Configuration file...\n";
	#pod2usage(1);
};#if

# Read the config filter file
our %config_filter=read_ini($config_filter_file,1);


## COMMON ##

# FORCE
our $force=$parameters{"force"};



return 1;
