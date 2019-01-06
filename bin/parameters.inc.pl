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

## Configuration files


## MAIN CONFIGURATION

# Config file
my $config_file;
if (-e $parameters{"config"}) {
	$config_file=$parameters{"config"};
} elsif (-e "$basename/".$parameters{"config"}) {
	$config_file="$basename/".$parameters{"config"};
} elsif (-e "$basename/config.ini") {
	$config_file="$basename/config.ini";
} else {
	warn "No Configuration file...\n";
	pod2usage(1);
};#if

# Set parameters from config file
if (-e $config_file) {
	
	# read ini file
	our %config=read_ini($config_file);
	# Folders
	$data_folder=$config{"folders"}{"data_folder"}."/";
	$annovar_databases=$config{"folders"}{"annovar_databases"}."/";

	# ANNOVAR folder (default "$basename")
	if (defined $config{"folders"}{"annovar_folder"} && -d $config{"folders"}{"annovar_folder"}) {
		$annovar_folder=$config{"folders"}{"annovar_folder"};
	} else {
		$annovar_folder=$basename;    
	};#if
	if (trim($annovar_folder) ne "") {$annovar_folder.="/"};

	# VCFTOOLS folder in $PATH by default (default "")
	if ((defined $config{"folders"}{"vcftools_folder"} && -d $config{"folders"}{"vcftools_folder"}) || trim($config{"folders"}{"vcftools_folder"}) eq "") {
		$vcftools_folder=$config{"folders"}{"vcftools_folder"};    
	} else {
		$vcftools_folder=$basename;    
	};#if
	if (trim($vcftools_folder) ne "") {$vcftools_folder.="/"};

	# R folder in $PATH by default (default "")
	if ((defined $config{"folders"}{"R_folder"} && -d $config{"folders"}{"R_folder"}) || trim($config{"folders"}{"R_folder"}) eq "") {
		$R_folder=$config{"folders"}{"R_folder"};    
	} else {
		$R_folder=$basename;    
	};#if
	if (trim($R_folder) ne "") {$R_folder.="/"};
	
	#tmp_folder
	if ((defined $config{"folders"}{"tmp_folder"} && -d $config{"folders"}{"tmp_folder"})) {
		$tmp_folder=$config{"folders"}{"tmp_folder"};    
	};#if
	
	# Database connexion
	$host = $config{"database"}{"host"};
	$driver = $config{"database"}{"driver"};
	$database= $config{"database"}{"database"};
	$user = $config{"database"}{"user"};
	$pw = $config{"database"}{"pw"};
	$port = $config{"database"}{"port"};
	#Project
	$assembly=$config{"project"}{"assembly"};
	$platform=$config{"project"}{"platform"};

};#if


## ANNOTATION CONFIGURATION

# Config annotation file
my $config_annotation_file;
if (-e $parameters{"config_annotation"}) {
    $config_annotation_file=$parameters{"config_annotation"};
} elsif (-e "$basename/".$parameters{"config_annotation"}) {
    $config_annotation_file="$basename/".$parameters{"config_annotation"};
} elsif (-e "$basename/config.annotation.ini") {
    $config_annotation_file="$basename/config.annotation.ini";
} else {
	print "#[ERROR] No Annotation Configuration file...\n";
	pod2usage(1);
};#if

# Read the config annotation file
our %config_annotation=read_ini($config_annotation_file);

# find annotation type
my %annotation_type;
while ((my $annotation_name, my $annotation_infos) = each(%config_annotation)){
	$annotation_type{$$annotation_infos{"annotation_type"}}{$annotation_name}=1;
};#while


## FILTER CONFIGURATION

# Config filter file
my $config_filter_file;
if (-e $parameters{"config_prioritization"}) {
    $config_filter_file=$parameters{"config_prioritization"};
} elsif (-e "$basename/".$parameters{"config_prioritization"}) {
    $config_filter_file="$basename/".$parameters{"config_prioritization"};
} elsif (-e "$basename/config.prioritization.ini") {
    $config_filter_file="$basename/config.prioritization.ini";
} else {
	print "No Filter Configuration file...\n";
	pod2usage(1);
};#if

# Read the config filter file
our %config_filter=read_ini($config_filter_file,1);


## COMMON ##

# FORCE
our $force=$parameters{"force"};




## Database connexion
#$dsn = "dbi:$driver:$database:$host:$port:mysql_local_infile=1:max_allowed_packet=512M";
#$dsn = "dbi:$driver:$database:$host:$port:mysql_local_infile=1;max_allowed_packet=1024M";
#$DBIconnect = DBI->connect($dsn, $user, $pw);

return 1;
