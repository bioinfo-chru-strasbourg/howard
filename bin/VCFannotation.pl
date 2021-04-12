#!/usr/bin/perl
#################################
# VCF Annotation using ANNOVAR  #
# Author: Antony Le Béchec      #
#################################

## Main Information
#####################

our %information = ( #
	'script'	=>  	basename($0),			# Script
	'release'	=>  	"0.9.10b",		# Release
	'date'		=>  	"20210411",		# Release parameter
	'author'	=>  	"Antony Le Béchec",	# Author
	'copyright'	=>  	"HUS",			# Copyright
	'licence'	=>  	"GNU AGPL V3",		# Licence
);

## Release Notes
##################
# 20150602-0.9.3b: deal with empty VCF (i.e. no variant)
# 20160519-0.9.4.2b: bug on VCF header with comma in the descrption of FILTER, INFO FORMAT
# 20160520-0.9.4.3b: bug on VCF without variant
# 20160520-0.9.7.1b: remove transcripts file for snpeff annotation
# 20181004-0.9.8b: add snpeff additional options. add  spliceSiteSize option, default 3. Add snpeff_split options
# 20190110-0.9.9b: add annovar_code_for_downdb parameter on annotation config file, generalise file parameter (no assembly by definition), bug fix
# 20210411-0.9.10b: remove snpEff -t option, fix snpEff -stats option, bug fix


## Modules
############

use Getopt::Long;		# Catch Options
use Pod::Usage;			# Pod
use Time::localtime;		# Time
use Data::Dumper;		# Data
use File::Basename;		# File
use Switch;			# Switch
#use File::Temp qw/ tempfile tempdir tmpnam /;
use File::Temp qw/ tmpnam /;
use File::Copy;
use File::stat;
use lib dirname (__FILE__);	# Add lib in the same folder


require "functions.inc.pl";	# Common functions
require "input.inc.pl";		# Input parameters


## HELP/MAN
#############

=head1 NAME

VCFannotation.pl - Annotate variants in the database or from an input file

=head1 DESCRIPTION

Description

=head1 BUGS

Bugs...

=head1 ACKNOWLEDGEMENTS

Thank U!

=head1 COPYRIGHT

HUS - GNU AGPL V3

=head1 AUTHOR

ALB

=head1 USAGE

$ARGV[0] [options] --input=<VCF> --output=<VCF> [...]

=head1 OPTIONS

=head2 MAIN options

=over 2

=item B<--help|h|?>

Print a brief help message and exits.

=item B<--man>

Prints the manual page and exits.

=back

=head2 CONFIGURATION

=over 2

=item B<--config|config_file=<file>>

Configuration file for main parameters (default 'config.ini')

=item B<--config_annotation|config_annotation_file=<file>>

Configuration file for annotation parameters (default 'config.annotation.ini').

=item B<--annovar_folder=<folder>>

Folder with ANNOVAR scripts.

=item B<--annovar_databases=<folder>>

Folder with ANNOVAR databases.


=item B<--snpeff_jar=<file>>

snpEff JAR file.

=item B<--snpeff_databases=<folder>>

Folder with snpEff databases.

=item B<--java>

java binary (needed if snpeff option on). default "java"

=item B<--java_flags>

java flags  (needed if snpeff option on, especially for configure proxy for databases donwload). default ""

=back

=head2 INPUT/OUPUT

=over 2

=item B<--input|input_file|vcf=<file>>

VCF Input file

=item B<--output=<file>>

VCF Output file

=back

=head2 OPTIONS

=over 2

=item B<--annotation=<string>>

Annotations sources, defined in the file 'config_annotation'. CASE SENSITIVE.
Example: "symbol", "Symbol,HGVS", "snpeff", "snpeff_split"

Specific snpEff options (if snpeff_jar and snpeff_databases defined):

"snpeff_split" to annotate "ANN" field and split snpEff annotation into "snpeff_hgvs", "snpeff_gene_name", "snpeff_annotation" and "snpeff_impact" fieds, and add to specific fields "symbol", "location" and "outcome" if empty

"snpeff" to annotate "ANN" field

"snpeff_hgvs" to annotate "ANN" field and annotate "snpeff_hgvs" field and add into "hgvs" if empty

"snpeff_gene_name" to annotate "ANN" field and annotate "snpeff_gene_name" field and add into "symbol" if empty

"snpeff_annotation" to annotate "ANN" field and annotate "snpeff_annotation" field and add into "location" and "outcome" if empty

"snpeff_impact" to annotate "ANN" field and annotate "snpeff_impact" field

=item B<--assembly=<string>>

Genome assembly to use. Default "hg19", "default" to use the default assembly in the configuration file.

=item B<--snpeff_stats=<file>>

Statistics from snpEff (--snpeff option will be turned on). default false.

=item B<--snpeff_spliceSiteSize=<integer>>

Set size for splice sites (donor and acceptor) in bases. Default: 3

=item B<--snpeff_additional_options=<string>>

additional options for snpEff (format "param1:value1|param2:value2"). default empty.

=back

=item B<--show_annotations>

List of annotations available

=back

=item B<--show_annotations_full>

List of annotations available with details

=back

=cut


## Parameters
###############

## Parameters default values


## Main parameters
$date=timestamp();
$basename = dirname (__FILE__);

## Header
$header="##\n";
$header.="## Script: ".$information{"script"}." (release ".$information{"release"}."/".$information{"date"}." - author ".$information{"author"}." - ".$information{"copyright"}." - ".$information{"licence"}.")\n";
$header.="## Excecution Date: ".$date."\n";
$header.="##\n";

## Parameters
###############

require "parameters.inc.pl";	# Parameters

## PrePorcessing
##################

## Help
if ($parameters{"help"}) {
	pod2usage(1);
};#if

## Man
if ($parameters{"man"}) {
	pod2usage(-exitstatus => 0, -verbose => 2);
};#if

## Header
if ($parameters{"header"}) {
	print $header;
	exit 0;
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
my $DEBUG=1 if $parameters{"debug"};

## Verbose
my $VERBOSE=1 if $parameters{"verbose"};


# Test ANNOVAR folder
if ( ! -d $annovar_folder) {
	print "#[ERROR] ANNOVAR folder '$annovar_folder' not found...\n";
	pod2usage(1);
};#if


# ANNOVAR annotation ype default
our $annovar_annotation_type_default=$parameters{"annovar_annotation_type_default"};


# Read the config annotation file
if (-e $config_annotation_file) {
	our %config_annotation=read_ini($config_annotation_file);
} else {
	print "#[ERROR] No Annotation Configuration file...\n";
	pod2usage(1);
};#if

# find annotation type
my %annotation_type;
while ((my $annotation_name, my $annotation_infos) = each(%config_annotation)){
	$annotation_type{$$annotation_infos{"annotation_type"}}{$annotation_name}=1;
};#while


## Input file
my $input_file;
if (-e $parameters{"input"}) {
	$input_file=$parameters{"input"};
	$header.="## Input VCF file: $input_file\n";
} else {
	print "# ERROR: input file '".$parameters{"input"}."' DOES NOT exist\n";
	pod2usage();
	exit 1;
};#if

## Output file
my $output_file;
my $output_filename_pattern="annotated";
if (-e $parameters{"output"} && 0) {
	print "# ERROR: output file '".$parameters{"output"}."' DOES exist\n";
	pod2usage();
	exit 1;
} else {
	$output_file=$parameters{"output"};
	if (trim($output_file) eq "") {
		$output_file=$input_file;
		$output_file =~ s/\.vcf$/\.$output_filename_pattern\.vcf/g; #/\.vcf$/\.output\.vcf/;
	};#if
	$header.="## Output VCF file: $output_file\n";
};#if


## Annotation
my %annotation_hash = map { $_ => undef } split(/,/, $parameters{"annotation"});
#print Dumper(\%annotation_hash) if $DEBUG;

#foreach annotation
while ((my $annotation_name, my $annotation_infos) = each(%annotation_hash)){
	# snpeff
	if ($annotation_name eq "snpeff") {
		$parameters{"snpeff"}=1;
		delete $annotation_hash{"snpeff"};
	} elsif ($annotation_name eq "snpeff_split") {
		$parameters{"snpeff_split"}=1;
		delete $annotation_hash{"snpeff_split"};
		$parameters{"snpeff_hgvs"}=1;
		delete $annotation_hash{"snpeff_hgvs"};
		$parameters{"snpeff_gene_name"}=1;
		delete $annotation_hash{"snpeff_gene_name"};
		$parameters{"snpeff_annotation"}=1;
		delete $annotation_hash{"snpeff_annotation"};
		$parameters{"snpeff_impact"}=1;
		delete $annotation_hash{"snpeff_impact"};
	} elsif ($annotation_name eq "snpeff_hgvs") {
		$parameters{"snpeff_hgvs"}=1;
		delete $annotation_hash{"snpeff_hgvs"};
	} elsif ($annotation_name eq "snpeff_gene_name") {
		$parameters{"snpeff_gene_name"}=1;
		delete $annotation_hash{"snpeff_gene_name"};
	} elsif ($annotation_name eq "snpeff_annotation") {
		$parameters{"snpeff_annotation"}=1;
		delete $annotation_hash{"snpeff_annotation"};
	} elsif ($annotation_name eq "snpeff_impact") {
		$parameters{"snpeff_impact"}=1;
		delete $annotation_hash{"snpeff_impact"};
	} else {
		# Annotation is configured
		if (exists $config_annotation{$annotation_name}) { #&& $config_annotation{$annotation_name}{"available"} eq "true") {
			$annotation_hash{$annotation_name}=$config_annotation{$annotation_name};
		};#if
		# Annotation is ALL available annotations
		if (uc($annotation_name) eq "ALL") {
			# add all annotation "available"
			while ((my $config_annotation_name, my $config_annotation_infos) = each(%config_annotation)){
				if ($$config_annotation_infos{"available"} eq "true"
					|| $$config_annotation_infos{"core"} eq "true") {
					$annotation_hash{$config_annotation_name}=$config_annotation{$config_annotation_name};
				};#if
			};#while
			# remove "ALL"
			delete $annotation_hash{$annotation_name};
		};#if
		# Annotation is CORE available annotations
		if (uc($annotation_name) eq "CORE") {
			# add all annotation "available"
			while ((my $config_annotation_name, my $config_annotation_infos) = each(%config_annotation)){
				if ($$config_annotation_infos{"core"} eq "true") { # don't need to be available
					$annotation_hash{$config_annotation_name}=$config_annotation{$config_annotation_name};
				};#if
			};#while
			# remove "ALL"
			delete $annotation_hash{$annotation_name};
		};#if
		# Annotation is a Type
		if (exists $annotation_type{$annotation_name}) {
			# add all annotation with this type
			while ((my $config_annotation_name, my $config_annotation_infos) = each(%config_annotation)){
				if ($$config_annotation_infos{"available"} eq "true"
					&& $$config_annotation_infos{"annotation_type"} eq $annotation_name) {
					$annotation_hash{$config_annotation_name}=$config_annotation{$config_annotation_name};
				};#if
			};#while
			# remove the type
			delete $annotation_hash{$annotation_name};
		};#if
	};#if

};#while


if ($parameters{"show_annotations"}) {
	#print Dumper(\%annotation_hash) if $DEBUG;
	#print Dumper(\%annotation_hash) if $DEBUG;
	#my @show_annotation;
	#my @show_annotation; #=keys(%annotation_hash);
	my @show_annotation=keys(%annotation_hash);
	#print @show_annotation if $DEBUG;
	if ($parameters{"snpeff_split"}) {
		#print "snpeff split\n" if $DEBUG;
		push(@show_annotation,"snpeff_split");
	} else {
		
		if ($parameters{"snpeff_hgvs"}) {
			#print "snpeff hgvs\n" if $DEBUG;
			push(@show_annotation,"snpeff_hgvs");
		};#if
		if ($parameters{"snpeff_gene_name"}) {
			push(@show_annotation,"snpeff_gene_name");
		};#if
		if ($parameters{"snpeff_annotation"}) {
			push(@show_annotation,"snpeff_annotation");
		};#if
		if ($parameters{"snpeff_impact"}) {
			push(@show_annotation,"snpeff_impact");
		};#if
		if ($parameters{"snpeff"}
			&& !$parameters{"snpeff_hgvs"}
			&& !$parameters{"snpeff_gene_name"}
			&& !$parameters{"snpeff_annotation"}
			&& !$parameters{"snpeff_impact"}
			) {
			push(@show_annotation,"snpeff");
		};#if
	};#if
	#push @show_annotation, keys(%annotation_hash);
	print "# ANNOTATIONS: @show_annotation";
	exit;
}; #if

if ($parameters{"show_annotations_full"}) {
	#print Dumper(\%annotation_hash) if $DEBUG;
	#print Dumper(\%annotation_hash); # if $DEBUG;
	print "\n";
	print "# List of annotations available\n";
	print "#################################\n";
	print "\n";
	
	my $nb_annotation=0;
	while ((my $annotation_hash_name, my $annotation_hash_infos) = each(%annotation_hash)){
		#print "$annotation_hash_name\n";
		if ($$annotation_hash_infos{"available"} eq "true") {
			my $core="";
			if ($$annotation_hash_infos{"core"} eq "true") { $core=" [core]"; }; 
			print "### $annotation_hash_name\n";
			print "# release: ".$$annotation_hash_infos{"release"}."/".$$annotation_hash_infos{"date"}."\n";
			print "# annotation_type: ".$$annotation_hash_infos{"annotation_type"}.$core."\n";
			#print "# core: ".$$annotation_hash_infos{"core"}."\n";
			print "# description: ".$$annotation_hash_infos{"description"}."\n";
			print "\n";
			$nb_annotation++;
		};#if
	};#while
	
	print "# $nb_annotation annotations available\n";
	print "\n";
	#@show_annotation=keys(%annotation_hash);
	
	#print "# ANNOTATIONS FULL: @show_annotation";
	exit;
}; #if

# SNPEFF
if ($parameters{"snpeff"} || $parameters{"snpeff_split"} || $parameters{"snpeff_hgvs"} || $parameters{"snpeff_gene_name"} || $parameters{"snpeff_annotation"} || $parameters{"snpeff_impact"}) {

	# JAVA
	#print $parameters{"java"}."\n";
	if ($parameters{"java"} ne "" && -e $parameters{"java"}) {
		#print "java1\n";
		$java=$parameters{"java"};
	} elsif (defined $config{"folders"}{"java"} && -e $config{"folders"}{"java"}) {
		#print "java2\n";
		$java=$config{"folders"}{"java"};
	} else {
		#print "java3\n";
		$java=trim(`which java`);
	};#
	if ($java ne "" && -e "$java" && -x _ ) {
		#print "# JAVA: $java\n";
	} else {
		print "# java '$java' not found...\n";
		pod2usage(1);
		#pod2usage(1);
	};#if
	#print "snpeff_jar=$snpeff_jar";



	if ($parameters{"snpeff_split"} || $parameters{"snpeff_hgvs"} || $parameters{"snpeff_gene_name"} || $parameters{"snpeff_annotation"} || $parameters{"snpeff_impact"}) {
		$parameters{"snpeff"}=1;
	};#if
	# TEST JAR (and databases?)
	# SNPEFF JAR
	if ($parameters{"snpeff_jar"} ne "") {
		$snpeff_jar=$parameters{"snpeff_jar"};
	};#if
	# Test SNPEFF JAR
	if ( ! -e $snpeff_jar) {
		print "# snpEff JAR not found...\n";
		pod2usage(1);
	};#if
	# SNPEFF databases
	#$snpeff_databases=$config{"folders"}{"snpeff_databases"};
	#if ($parameters{"snpeff_databases"} ne "") {
	#	$snpeff_databases=$parameters{"snpeff_databases"};
	#};#if
	# JAVA test
	#system("java -version");
	#$java_version = `$java -version`;
	#print "# JAVA version: $java_version" if $VERBOSE;
	#print $snpeff_jar if $VERBOSE;
	if (!-e $snpeff_jar) {
		print "# snpEff JAR not found...\n";
		pod2usage(1);
	};#if


};#if


# Case of ALL annotation
# list all available annotation sources
if ($parameters{"annotations"}) {
	# width column
	$max_annotation_name=0;
	$max_release=0;
	$max_description=0;
	while ((my $annotation_name, my $annotation_infos) = each(%config_annotation)){
		$max_annotation_name=(length($annotation_name)>$max_annotation_name)?length($annotation_name):$max_annotation_name;
		$max_release=(length($$annotation_infos{"release"})>$max_release)?length($$annotation_infos{"release"}):$max_release;
		$max_description=(length($$annotation_infos{"description"})>$max_description)?length($$annotation_infos{"description"}):$max_description;
	};#while
	while ((my $annotation_name, my $annotation_infos) = each(%config_annotation)){
		printf("%-".($max_annotation_name+5)."s%-".($max_release+5)."s%-".($max_description+5)."s","# $annotation_name","[".$$annotation_infos{"release"}."]",$$annotation_infos{"description"});
		print "\n";
		#print "\t".$$annotation_infos{"annovar_code"}." (".$$annotation_infos{"annovar_annotation_type"}.")\n";
		#printf "# $annotation_name\t\t[".$$annotation_infos{"release"}."]\t".$$annotation_infos{"annovar_code"}." (".$$annotation_infos{"annovar_annotation_type"}.")\n", 30;
		#print "#    - ANNOVAR code: ".$config_annotation{"annovar_code"}." (".$config_annotation{"annovar_annotation_type"}.")\n";
	};#while
	exit 0;
};#if


## GRANTHAM Matrix
if (!-e "$annovar_databases/grantham.matrix") {
	if (-e "$annovar_folder/grantham.matrix") {
		my $cmd_cp="cp $annovar_folder/grantham.matrix $annovar_databases/grantham.matrix";
		my $result = `$cmd_cp 2>&1`;
	};#if
};#if

## ANNOVAR scripts
my $annotate_variation=$annovar_folder."/annotate_variation.pl";
my $convert2annovar=$annovar_folder."/convert2annovar.pl";

# Test ANNOVAR scripts
if ( ! -e $annotate_variation) {
	print "# ANNOVAR script '$annotate_variation' not found...\n";
	pod2usage(1);
};#if
if ( ! -e $convert2annovar) {
	print "# ANNOVAR script '$convert2annovar' not found...\n";
	pod2usage(1);
};#if


## temp dir
#my $tempdir = tempdir()."/";
#mkdir($tempdir) if (!-d $tempdir);
#my $tempdir = tempdir()."/";
my $tmpdir = "/tmp/".rand(1000000);
if ($parameters{"tmp"} ne "") {
	$tmpdir = $parameters{"tmp"}."/".rand(1000000);
};#
#if ( ! -d $tmpdir ) {
#	mkdir($tmpdir);
#};#if

## DEBUG
##########

#print Dumper(\%parameters) if $DEBUG;
#print Dumper(\%config_annotation) if $DEBUG;


## MAIN
#########

# Config
if ($VERBOSE) {
	print "#[INFO] Assembly: $assembly\n";
	print "#[INFO] ANNOVAR folder: $annovar_folder\n";
	print "#[INFO] ANNOVAR Databases: $annovar_databases\n";
	print "#[INFO] snpEff JAR: $snpeff_jar\n";
	print "#[INFO] snpEff Databases: $snpeff_databases\n";
	print "#[INFO] java: $java\n";
	print "#[INFO] java flags: $java_flags\n";
	#exit 0;
};#if

# Variables
$output="";

# Verbose
if ($VERBOSE) {
	print "#[INFO] Input file: $input_file\n";
	print "#[INFO] Output file: $output_file\n";
	#print "# Output format: @output_formats\n";
	#print "#[INFO] Filter: ".$parameters{"annotation"}."\n";
	print "#\n";
};#if

## VCF into ANNOVAR

# Test NB of variant
my $nb_variant=`grep -cv ^# $input_file`+0;
if ($VERBOSE) {
	print "#[INFO] Nb of variants: $nb_variant\n#\n";
};#if
if ($nb_variant eq 0) {
	print "#[INFO] NO variant to annotate\n";
	print "#[INFO] Copy Input File into Output File\n#\n";
	copy($input_file,$output_file) or die "#[ERROR] Copy failed: $!";
	exit 0;
};#if

#print "FORCEFFFFFFFFFFFFFFFFFFFFFFF=$force\n" if $DEBUG;
#print "FORCEFFFFFFFFFFFFFFFFFFFFFFF=".$parameters{"force"}."\n" if $DEBUG;


###############
### ANNOVAR ###
###############

print "#\n### ANNOVAR Annotation ###\n#" if $VERBOSE;

# ANNOVAR temporary file
#my $annovar_file=tmpnam();
my $annovar_file=$tmpdir.".".rand(1000000);

#print "#\n### ANNOVAR File TMP = $tmpdir / $annovar_file ###\n#" if $DEBUG;

# ANNOVAR convert command
my $cmd_convert2annovar="perl $convert2annovar --format vcf4old --includeinfo --allallele --outfile $annovar_file $input_file 2>&1";
#my $cmd_convert2annovar="perl $convert2annovar --format vcf4 --includeinfo --outfile $annovar_file $input_file 2>&1";
print "$cmd_convert2annovar\n" if $DEBUG;

# Launch command
my $result = `$cmd_convert2annovar`;


# header
my %vcf_header_hash_test=read_vcf($parameters{"input"},"header");
#my %header_INFO_list;
#while ( my ($ann, $poss) = each(%{$vcf_header_hash_test{"INFO"}}) ) {
#	print "$ann, $poss\n" if $DEBUG;
#	$header_INFO_list{$ann}=$poss;
#};#while
#print Dumper(\%header_INFO_list) if $DEBUG;
#print $header_INFO_list{"hgvs"}."\n" if $DEBUG;
#if (defined $header_INFO_list{"hgvs"}) {
#	print "OK\n" if $DEBUG;
#};#if
#if (defined $vcf_header_hash_test{"INFO"}{"hgvs"}) {
#	print "OK\n" if $DEBUG;
#};#if

# Test NB of variant
#$annovar_file

## Foreach Annotation Sources (defined)
my $annotation_nb=scalar(keys %annotation_hash);
my $annotation_i=0;
my %annotation_output;
my $vcf_header_INFOS="";
my %vcf_header;
my $vcf_header_information;
while ((my $annotation_name, my $annotation_infos) = each(%annotation_hash)){

	#if ($force || !defined $$variant_values{"INFOS"}{"GenotypeConcordance"}) {
	# if no in header or force
	#print "TEST $annotation_name, $annotation_infos\n" if $DEBUG;
	#print "DEFINED ".$header_INFO_list{$annotation_name}."\n" if $DEBUG;
	#if (defined $header_INFO_list{$annotation_name}) {
	#if (!($force || !defined $vcf_header_hash_test{"INFO"}{$annotation_name})) {
	if (!$force && defined $vcf_header_hash_test{"INFO"}{$annotation_name}) {
		#$vcf_header_hash_test{"INFO"}{"hgvs"}
		my $step="# [$annotation_i/$annotation_nb] Annotation '$annotation_name' skipped because already in the VCF header (force to annotate)\n";
		$output.=$step;

		next;
		#print "$annotation_name defined !!!!!!\n" if $DEBUG;
	};#if

	## init
	$annotation_i++;
	my $error=0;

	## Output
	my $step="# [$annotation_i/$annotation_nb] Annotation '$annotation_name'\n";
	$output.=$step;
	print $step if $VERBOSE;
	my $output_verbose;
	my $output_verbose_result;

	## Check if annotation is defined
	my $annovar_code="";
	if (defined($annotation_infos)) {
		$annovar_code=$config_annotation{$annotation_name}{"annovar_code"};
	} else {
		print "# NO configuration for Annotation '$annotation_name'. Try with ANNOVAR code and default parameters\n" if $VERBOSE;
		$annovar_code=$annotation_name;
	};#if

	if ($annovar_code ne "") {

		## init
		my $annovar_annotation_type_cmd;
		my $output_extension_from_input_file;
		my $annovar_file_output=$annovar_file.".".$annotation_name;

		## Command options

		# Mandatory options in the config file
		#my $annovar_code=$config_annotation{$annotation_name}{"annovar_code"};
		my $annovar_annotation_type=$config_annotation{$annotation_name}{"annovar_annotation_type"};
		my $release=$config_annotation{$annotation_name}{"release"};
		#$release.=($release eq "")?"unknown":"";
		my $date=$config_annotation{$annotation_name}{"date"};

		# Possibily calculated option
		my $annotation_type=$config_annotation{$annotation_name}{"annotation_type"};
		my $annotation_description=trim($config_annotation{$annotation_name}{"description"});

		my $file=$config_annotation{$annotation_name}{"file"};
		if ($file eq "") {
			$file="".$assembly."_".$annovar_code.".txt";
		} else {
		    $file=$assembly."_".$file;
		};#if
		my $annovar_code_for_downdb=$config_annotation{$annotation_name}{"annovar_code_for_downdb"};
		if ($annovar_code_for_downdb eq "") {
		    $annovar_code_for_downdb=$annovar_code;
		};#if
		my $online_file=$config_annotation{$annotation_name}{"online_file"};
		my $output_file_extension=$config_annotation{$annotation_name}{"output_file_extension"};
		my $additional_options=$config_annotation{$annotation_name}{"additional_options"};
		my $colsWanted_options=$config_annotation{$annotation_name}{"colsWanted"};
		$additional_options.=" --colsWanted=$colsWanted_options " if (trim($colsWanted_options) ne "");
		my $otherinfo_options=$config_annotation{$annotation_name}{"otherinfo"};
		$additional_options.=" --otherinfo " if (trim($otherinfo_options) ne "");
		# column result
		my $column_result_default=1;
		my $column_result_config=$config_annotation{$annotation_name}{"column_result"};
		my $column_result=((trim($column_result_config) ne "") && (($column_result_config+0) gt -1))?$column_result_config:$column_result_default;
		# column chr
		my $column_chr_default=7;
		my $column_chr_config=$config_annotation{$annotation_name}{"column_chr"};
		my $column_chr=(trim($column_chr_config) ne "" && ($column_chr_config+0) gt -1)?$column_chr_config:$column_chr_default;


		## Command option preprocess

		# Annovar Code
		if (trim($annovar_code) eq "") {
			$output_verbose.="# ERROR: 'annovar_code' ('$annovar_code') not correct\n";
			next;
		};#if
		$output_verbose.="#    - annovar_code='$annovar_code'\n";

		# File
		if ($file eq "" || !-e "$annovar_databases/$file") {
			
			# try to download from ANNOVAR
			if ($file eq "" || !-e "$annovar_databases/$file") {
				#if ( -e "./ANNOVAR_download.sh") {
					#my $cmd="./ANNOVAR_download.sh '$assembly' '$annovar_code_for_downdb' '$annovar_databases' '' '$annotate_variation' ''";
					my $cmd="$annotate_variation -downdb -buildver $assembly $annovar_code_for_downdb $annovar_databases -webfrom annovar";# $ANNOVAR_SCRIPT -downdb -buildver $BUILD $DB $DB_FOLDER $WEBFROM
					$output_verbose.="#    - cmdDL='$cmd'\n";
					print "# Try to download database from ANNOVAR...\n" if $VERBOSE;
					print "# cmd=$cmd...\n" if $VERBOSE;
					#print "# cmd=$cmd\n" if $DEBUG;
					my $result = `$cmd 2>&1`; #
				#};#if
			};#if
			if ($file eq "" || !-e "$annovar_databases/$file") {
				#if ( -e "./ANNOVAR_download.sh") {
					#my $cmd="./ANNOVAR_download.sh '$assembly' '$annovar_code_for_downdb' '$annovar_databases' '' '$annotate_variation' ''";
					my $cmd="$annotate_variation -downdb -buildver $assembly $annovar_code_for_downdb $annovar_databases -webfrom annovar -nowget";# $ANNOVAR_SCRIPT -downdb -buildver $BUILD $DB $DB_FOLDER $WEBFROM
					$output_verbose.="#    - cmdDL='$cmd'\n";
					print "# Try to download database from ANNOVAR without wget...\n" if $VERBOSE;
					print "# cmd=$cmd...\n" if $VERBOSE;
					#print "# cmd=$cmd\n" if $DEBUG;
					my $result = `$cmd 2>&1`; #
				#};#if
			};#if
			# try to download from UCSC
			if ($file eq "" || !-e "$annovar_databases/$file") {
				#if ( -e "./ANNOVAR_download.sh") {
					my $cmd="$annotate_variation -downdb -buildver $assembly $annovar_code_for_downdb $annovar_databases ";# $ANNOVAR_SCRIPT -downdb -buildver $BUILD $DB $DB_FOLDER $WEBFROM
					$output_verbose.="#    - cmdDL='$cmd'\n";
					print "# Try to download database from UCSC...\n" if $VERBOSE;
					print "# cmd=$cmd...\n" if $VERBOSE;
					#print "# cmd=$cmd\n" if $DEBUG;
					my $result = `$cmd 2>&1`; #
				#};#if
			};#if
			# try to download from UCSC
			if ($file eq "" || !-e "$annovar_databases/$file") {
				#if ( -e "./ANNOVAR_download.sh") {
					if ($online_file ne "") {
						my $cmd="wget -O $annovar_databases/$file $online_file ";# $ANNOVAR_SCRIPT -downdb -buildver $BUILD $DB $DB_FOLDER $WEBFROM
						$output_verbose.="#    - cmdDL='$cmd'\n";
						print "# Try to download database from online file '$online_file'...\n" if $VERBOSE;
						print "# cmd=$cmd...\n" if $VERBOSE;
						#print "# cmd=$cmd\n" if $DEBUG;
						my $result = `$cmd 2>&1`; #
					};#if
				#};#if
			};#if
			if ($file eq "" || !-e "$annovar_databases/$file") {
				$output.="# ERROR: 'file' ('$file') not correct\n";
				next;
			};#if
		};#if
		$output_verbose.="#    - file='$file'\n";

		# Creation date
		my $cmd_date="stat -c \%y $annovar_databases/$file | cut -d' ' -f1 | tr -d '-'";
		my $file_creation_date = trim(`$cmd_date 2>&1`);
		$output_verbose.="#    - file_date='$file_creation_date'\n";

		$annovar_annotation_type.=($annovar_annotation_type eq "")?"$annovar_annotation_type_default":"";
		$release.=($release eq "")?$file_creation_date:"";
		$date.=($date eq "")?$file_creation_date:"";
		$annotation_type.=($annotation_type eq "")?"unknown":"";
		$annotation_description.=($annotation_description eq "")?$file:"";


		# output Extension and command type
		my $output_extension_from_input_file=$file;
		$output_extension_from_input_file =~s/\.txt$//;
		switch ($annovar_annotation_type){
			case("geneanno") {
				$output_file_extension=((trim($output_file_extension) eq "")?".variant_function":$output_file_extension);
				$annovar_annotation_type_cmd="-geneanno -dbtype ";
			}
			case("filter") {
				$output_file_extension=((trim($output_file_extension) eq "")?".".$output_extension_from_input_file."_dropped":$output_file_extension);
				$annovar_annotation_type_cmd="-filter -dbtype ";
			}
			case("regionanno") {
				$output_file_extension=((trim($output_file_extension) eq "")?".".$output_extension_from_input_file:$output_file_extension);
				$annovar_annotation_type_cmd="-regionanno -dbtype ";
			}
			case("vcf") {
				$output_file_extension=((trim($output_file_extension) eq "")?".".$assembly."_vcf_dropped":$output_file_extension);
				$annovar_annotation_type_cmd="-filter -dbtype vcf -vcfdbfile ";
			}
			case("generic") {
				$output_file_extension=((trim($output_file_extension) eq "")?".".$assembly."_generic_dropped":$output_file_extension);
				$annovar_annotation_type_cmd="-filter -dbtype generic -genericdbfile ";
			}
			else {
				$output_verbose.="# ERROR: 'annovar_annotation_type' ('$annovar_annotation_type') not correct\n";
				$error=1;
			}# else
		}
		$output_verbose.="#    - output_file_extension='$output_file_extension'\n";

		## Other options
		$output_verbose.="#    - column_result='$column_result'\n";
		$output_verbose.="#    - column_chr='$column_chr'\n";
		$output_verbose.="#    - annovar_annotation_type='$annovar_annotation_type'\n";
		$output_verbose.="#    - additional_options='$additional_options'\n";
		$output_verbose.="#    - annotation_description='$annotation_description'\n";


		## Command
		my $cmd="perl $annotate_variation $annovar_annotation_type_cmd $annovar_code --buildver $assembly $additional_options --outfile $annovar_file_output $annovar_file $annovar_databases";
		$output_verbose.="#    - cmd='$cmd'\n";

		if ($error) {
			$output_verbose.="# [ERROR] next\n";
			next;
		};#if

		print $output_verbose if $VERBOSE;

		## Launch Command
		my $result = `$cmd 2>&1`;

		## Read output file
		my $annovar_file_output_results=$annovar_file_output.$output_file_extension;
		my $line_num=0;
		my $nb_annotation=0;
		my $output_results;
		if ( -e $annovar_file_output_results ) {
			open(ANNOVAR_file_output_results, "$annovar_file_output_results") || die "# ERROR: '$annovar_file_output_results' $!";
			while(<ANNOVAR_file_output_results>)
			{
				#print $_ if $DEBUG;
				# init
				chomp; #delete \n character
				$line_num++;
				my $line=$_;
				my @line_split=split("\t",$line);

				# Results
				my $result=$line_split[$column_result];
				my $chr=$line_split[$column_chr];
				my $pos=$line_split[($column_chr+1)];
				my $id=$line_split[($column_chr+2)];
				my $ref=$line_split[($column_chr+3)];
				my $alt=$line_split[($column_chr+4)];

				# OtherInfo
				if (trim($otherinfo_options) ne "") {
					my @result_split=split(",",$result);
					my $result_concat;
					foreach my $i (split(",",$otherinfo_options)) {
						$result_concat.=(($result_concat ne "")?",":"").@result_split[($i-1)] if (@result_split[($i-1)] ne ".");
					};#foreach
					#$result=$result_concat if (trim($result_concat) ne "");
					$result=$result_concat;
				};#if

				# Results PostProcess
				$result =~ s/;/,/g;
				$result =~ s/=/:/g;
				$result =~ s/ /_/g;
				if ( $annotation_name eq "Symbol" || $annotation_name eq "symbol" ) {
					print "$result\n" if $DEBUG;
					@result_split=split(/[\(,]/,$result);
					print "@result_split\n" if $DEBUG;
					$result=$result_split[0];
					print "$result_split[0]\n" if $DEBUG;
				}; #if

				# Output
				if (trim($result) ne "") {
					$nb_annotation++;
					#my $sep=((trim($annotation_output{$chr}{$pos}{$ref}{$alt}{$annotation_name}) eq "")?"":",");
					#if ($annotation_output{$chr}{$pos}{$ref}{$alt}{$annotation_name} eq "$sep$result") {
					#	#$annotation_output{$chr}{$pos}{$ref}{$alt}{$annotation_name}=$result;
					#};#if
					#$annotation_output{$chr}{$pos}{$ref}{$alt}{$annotation_name}.="$sep$result" if ($annotation_output{$chr}{$pos}{$ref}{$alt}{$annotation_name} ne $result);

					my $sep_default=",";
					#print "#       TEST:\t".trim($annotation_output{$chr}{$pos}{$ref}{$alt}{$annotation_name})." ??? ".$result."\n";
					if (trim($annotation_output{$chr}{$pos}{$ref}{$alt}{$annotation_name}) eq ""
						|| trim($annotation_output{$chr}{$pos}{$ref}{$alt}{$annotation_name}) eq $sep_default.trim($result)) {
						$annotation_output{$chr}{$pos}{$ref}{$alt}{$annotation_name}=$result;
					} elsif (trim($annotation_output{$chr}{$pos}{$ref}{$alt}{$annotation_name}) ne trim($result)) {
						$annotation_output{$chr}{$pos}{$ref}{$alt}{$annotation_name}.=$sep_default.$result;
					};#if
					#print "#       $chr:$pos:$ref:$alt\t".$annotation_output{$chr}{$pos}{$ref}{$alt}{$annotation_name}."\n";
					$output_results.="#       $chr:$pos:$ref:$alt\t".$annotation_output{$chr}{$pos}{$ref}{$alt}{$annotation_name}."\n";
				};#if

			};#while
			$nb_results_to_show=10;
			$index_results_to_show=$nb_results_to_show-1;
			$output_verbose_result.="#    - $nb_annotation Results for '$annotation_name'\n";
			$output_verbose_result.="#    - first $nb_results_to_show results for '$annotation_name':\n";
			#$output_verbose_result.=$output_results;
			#$output_verbose_result.=( split /\n/, $output_results)[0-10]."\n";
			$output_verbose_result.= (join "\n", ( split /\n/, $output_results)[0 .. $index_results_to_show])."\n";

			$output.="#       $nb_annotation annotated variants for '$annotation_name'\n";

			# VCF_header_INFOS
			if ($nb_annotation>0) {
				$annotation_description =~ s/,//g;
				$vcf_header_INFOS.="##INFO=<ID=$annotation_name,Number=.,Type=String,Description=\"Annotation '$annotation_name' (release:$release date:$date type:$annotation_type).".((trim($annotation_description) eq "")?"":" ")."$annotation_description\">\n";
				# Mandatory
				$vcf_header{"INFO"}{$annotation_name}{"Number"}=".";
				$vcf_header{"INFO"}{$annotation_name}{"Type"}="String";
				my $description=((trim($annotation_description) eq "")?"unknown":trim($annotation_description));
				my $additional_description=" [Release=$release;Date=$date;AnnotationType=$annotation_type]";
				$vcf_header{"INFO"}{$annotation_name}{"Description"}="\"$description$additional_description\"";
				# Optional but mandotory for DB	(included in description TAG)
				#my $additional_description="[Release=$release,Date=$date,AnnotationType=$annotation_type]";
				#$vcf_header{"INFO"}{$annotation_name}{"Description"}.=" $additional_description";
				# Optional but mandotory for DB
				#$vcf_header{"INFO"}{$annotation_name}{"Release"}="$release";
				#$vcf_header{"INFO"}{$annotation_name}{"Date"}="$date";
				#$vcf_header{"INFO"}{$annotation_name}{"AnnotationType"}="$annotation_type";
				#$vcf_header{"INFO"}{$annotation_name}{"Description"}="\"Annotation '$annotation_name' (release:$release date:$date type:$annotation_type).".((trim($annotation_description) eq "")?"":" ")."$annotation_description\"";
			};#if

			#UNLINK
			my @tmp_files=<$annovar_file_output*>;
			unlink(@tmp_files);

		} else {

			$output_verbose_result.="#    ERROR: NO ANNOVAR result file '$annovar_file_output_results'\n"; #;

		};#if

	} else {

		## Output
		$output_verbose_result.="#    ERROR: NO configuration for Annotation '$annotation_name'\n"; #;

	};#if

	print $output_verbose_result if $VERBOSE;

};#foreach

# Remove Annovar filename
unlink($annovar_file);

##############
### SNPEFF ###
##############

my %annotation_snpeff_output;
my $annotation_snpeff_output_header;
my @snpeff_fields_definition_split;
my %info_split_index_array;


if ($parameters{"snpeff"} || $parameters{"snpeff_stats"}) {

	print "#\n### snpEff Annotation ###\n#" if $VERBOSE;
	
	my $info_snpeff_name="ANN";
	my $snpeff_spliceSiteSize_default="3";
	
	if (!$force && defined $vcf_header_hash_test{"INFO"}{$info_snpeff_name}) {
	#if (0) {
		#$vcf_header_hash_test{"INFO"}{"hgvs"}
		#my $step="# [$annotation_i/$annotation_nb] Annotation '$annotation_name' skipped because already in the VCF header (force to annotate)\n";
		#$output.=$step;
		print "#\n### snpEff Annotation skipped because already in the header\n#" if $VERBOSE;

		#print "$annotation_name defined !!!!!!\n" if $DEBUG;

	} else {

		# SNPEFF temporary file
		my $input_file_snpeff_stats=tmpnam();
		my $output_file_snpeff=tmpnam();
		my $output_file_snpeff_log=tmpnam();
		#$tmpdir.".".rand(1000000);
		my $input_file_snpeff_stats=$tmpdir.".".rand(1000000);
		my $output_file_snpeff=$tmpdir.".".rand(1000000);
		my $output_file_snpeff_log=$tmpdir.".".rand(1000000);
		my $output_file_snpeff_transcripts=$tmpdir.".".rand(1000000);
		my $output_verbose="";

		@snpeff_additional_options_split=split(/\|/,$parameters{"snpeff_additional_options"});
		my $snpeff_additional_options="";
		foreach my $option (@snpeff_additional_options_split) {
			my @option_split=split(/:/,$option);
			$snpeff_additional_options.=" -".$option_split[0];
			if (defined $option_split[1]) {
				$snpeff_additional_options.=" ".$option_split[1];
			};#if
		};#foreach		
		
		#$parameters{"snpeff_additional_options"}
		
		#my $snpeff_options=" -dataDir $snpeff_databases -spliceSiteSize ".$parameters{"snpeff_spliceSiteSize"}." ".$snpeff_additional_options; #.$parameters{"snpeff_additional_options"};
		my $snpeff_options=" -dataDir $snpeff_databases $snpeff_additional_options"; #.$parameters{"snpeff_additional_options"};
		#-nodownload -noShiftHgvs
		#if ($parameters{"snpeff_stats"} eq "" || ! -e $input_file_snpeff_stats) {

		# Splice Sites Size
		if (1) {
			my $spliceSiteSize_input=$parameters{"snpeff_spliceSiteSize"};
			if ($spliceSiteSize_input ne "") { #|| ! -e $input_file_snpeff_stats) {
				$snpeff_options.=" -spliceSiteSize $spliceSiteSize_input ";
				$output_verbose.="#    - spliceSiteSize='$spliceSiteSize_input'\n";
			} else {
				$snpeff_options.=" ";
			};#if
		}

		# Stats
		if (1) {
			my $stats_input=$parameters{"snpeff_stats"};
			if ($stats_input ne "") { #|| ! -e $input_file_snpeff_stats) {
				$snpeff_options.=" -stats $stats_input ";
				$output_verbose.="#    - stats='$stats_input'\n";
				$output_verbose.="#    - stats='$stats_input'.genes.txt\n";
			} else {
				$snpeff_options.=" -noStats ";
			};#if
		};#if

		# TRANSCRIPTS
		if (0) {
			my $transcripts_input=$parameters{"transcripts"};
			if (-e $transcripts_input) { #if file exist
				#my $cmd="awk -F\"\\t\" '{print \$2\"\t\"\$1}' $transcripts_input > $output_file_snpeff_transcripts";
				my $cmd="awk -F\"\\t\" '(\$1!=\”\"){print \$1}' $transcripts_input > $output_file_snpeff_transcripts";
				#$output_verbose.="#    - cmd='$cmd'\n";
				#print $output_verbose if $VERBOSE;
				#my $result = `$cmd 2>&1`;
				$snpeff_options.=" -onlyTr $output_file_snpeff_transcripts ";
				$output_verbose.="#    - onlyTr '$transcripts_input'\n";
			};#if
		};#if


		## Command
		#my $cmd="$java $java_flags -Xmx4g -jar $snpeff_jar $assembly $input_file -stats $input_file_snpeff_stats $snpeff_options -noLog 1>$output_file_snpeff 2>$output_file_snpeff_log";
		my $cmd="$java $java_flags -Xmx4g -jar $snpeff_jar $assembly $input_file $snpeff_options -noLog 1>$output_file_snpeff 2>$output_file_snpeff_log";
		$output_verbose.="##    - cmd='$cmd'\n";
		print $output_verbose if $VERBOSE;
		
		## Launch Command
		my $result = `$cmd 2>&1`;

		# STATS
		# if ($parameters{"snpeff_stats"} ne "" && -e $input_file_snpeff_stats) {
		# 	#my $cmd="cp $input_file_snpeff_stats ".$parameters{"snpeff_stats"};
		# 	#$output_verbose.="#    - cmd='$cmd'\n";
		# 	#print $output_verbose if $VERBOSE;
		# 	#my $result = `$cmd 2>&1`;
		# 	my $cmd="cp $input_file_snpeff_stats.genes.txt ".$parameters{"snpeff_stats"}.".genes.txt";
		# 	$output_verbose.="#    - cmd='$cmd'\n";
		# 	print $output_verbose if $VERBOSE;
		# 	my $result = `$cmd 2>&1`;
		# };#if


		#print $output_file_snpeff_log if $VERBOSE;

		

		# VARIANT ANNOTATION
		%annotation_snpeff_output=read_vcf($output_file_snpeff);
		#my %vcf_header_hash_test=read_vcf($parameters{"input"},"header");
		my %vcf_header_hash_snpeff=read_vcf($output_file_snpeff,"header");
		#print Dumper(\%annotation_output) if $VERBOSE;
		#print Dumper(\%info_split_index_array) if $VERBOSE;
		#foreach (@%info_split_index_array) { print $_; };
		
		# Split ANN header

		my $header_content_snpeff=$vcf_header_hash_snpeff{"INFO"}{"ANN"}{"Description"};
		#print "header_content_snpeff=$header_content_snpeff\n";
		if ($header_content_snpeff=~ /^(.*)'(.*)'(.*)$/) {
			#print "1:$1 \n 2:$2 \n 3:$3 \n" if $VERBOSE;
			my $snpeff_fields_definition=$2;
			$snpeff_fields_definition =~ s/ //gi;
			#print "snpeff_fields_definition=$snpeff_fields_definition \n" if $VERBOSE;
			@snpeff_fields_definition_split=split(/\|/,$snpeff_fields_definition);
			#print "DEF: @snpeff_fields_definition_split\n" if $VERBOSE;
			my $info_split_index=-1;
			foreach (@snpeff_fields_definition_split) {
				#print "$_\n";
				$info_split_index++;
				#print "$info_split_index $_\n" if $VERBOSE;
				$info_split_index_array{$_}=$info_split_index;
				#if ( $_ eq "Feature_ID" ) { @info_split_index_array };#if
			};#foreach
		};#if
		#$annotation_snpeff_output_header=$line_content;



		while ( my ($chr, $poss) = each(%annotation_snpeff_output) ) {
		while ( my ($pos, $refs) = each(%{$poss}) ) {
		while ( my ($ref, $alts) = each(%{$refs}) ) {
		while ( my ($alt, $variant_values) = each(%{$alts}) ) {
			#print "# {$chr}{$pos}{$ref}{$alt} ".$annotation_output{$chr}{$pos}{$ref}{$alt}{"Symbol"}."\n" if $VERBOSE;
			#print "# {$chr}{$pos}{$ref}{$alt} ".$annotation_snpeff_output{$chr}{$pos}{$ref}{$alt}{"INFOS"}{"ANN"}."\n" if $VERBOSE;
			$annotation_output{$chr}{$pos}{$ref}{$alt}{$info_snpeff_name}=$annotation_snpeff_output{$chr}{$pos}{$ref}{$alt}{"INFOS"}{"ANN"};

			if (1) {
				# HEADER Mandatory
				my $annotation_name=$info_snpeff_name;
				my $annotation_description="snpEff annotation";
				my $additional_description=""; #" [Release=$release;Date=$date;AnnotationType=$annotation_type]";
				$annotation_description =~ s/,//g;
				$additional_description =~ s/,//g;
				#$vcf_header{"INFO"}{$annotation_name}{"Number"}="."; # $vcf_header_hash_snpeff{"INFO"}{"ANN"}{"Number"};
				#$vcf_header{"INFO"}{$annotation_name}{"Type"}="String";
				#$vcf_header{"INFO"}{$annotation_name}{"Description"}="\"$annotation_description$additional_description\"";
				$vcf_header{"INFO"}{$annotation_name}{"Number"}=$vcf_header_hash_snpeff{"INFO"}{"ANN"}{"Number"};
				$vcf_header{"INFO"}{$annotation_name}{"Type"}=$vcf_header_hash_snpeff{"INFO"}{"ANN"}{"Type"};
				$vcf_header{"INFO"}{$annotation_name}{"Description"}=$vcf_header_hash_snpeff{"INFO"}{"ANN"}{"Description"};

			};

			if ($parameters{"snpeff_split"} || $parameters{"snpeff_hgvs"} || $parameters{"snpeff_gene_name"} || $parameters{"snpeff_annotation"} || $parameters{"snpeff_impact"}) {
				#my @snpeff_fields_values_split=split(",",$annotation_snpeff_output{$chr}{$pos}{$ref}{$alt}{"INFOS"}{"ANN"});
				my @snpeff_fields_values_split=split(",",$annotation_output{$chr}{$pos}{$ref}{$alt}{$info_snpeff_name});
				#foreach (@snpeff_fields_values_split) { print "$_\n"; };
				my %snpeff_hgvs_list;
				my %snpeff_gene_id_list;
				my %snpeff_gene_name_list;
				my %snpeff_annotation_list;
				my %snpeff_impact_list;
				foreach my $info (@snpeff_fields_values_split) {
					my @info_split=split(/\|/,$info);
					#foreach (@info_split) { print "$_\n"; };
					#print "INFO: @info_split\n" if $VERBOSE;
					#my $info_split_index=0;
					#foreach (@info_split) { $info_split_index++; print "$info_split_index $_\n"; };
					if ($info_split[$info_split_index_array{"Feature_ID"}] ne ""
						&& $info_split[$info_split_index_array{"HGVS.c"}] ne "") {
						my $nb_infos=0;
						my @snpeff_hgvs_transcript_split=split(/:/,$info_split[$info_split_index_array{"Feature_ID"}]);
						my $snpeff_hgvs_transcript=$snpeff_hgvs_transcript_split[0];
						#my $snpeff_hgvs=$info_split[$info_split_index_array{"Feature_ID"}].":".$info_split[$info_split_index_array{"HGVS.c"}];
						my $snpeff_hgvs=$snpeff_hgvs_transcript.":".$info_split[$info_split_index_array{"HGVS.c"}];
						# Gene ID
						if ($info_split[$info_split_index_array{"Gene_Name"}] ne "") {
							#$snpeff_hgvs=$info_split[$info_split_index_array{"Gene_Name"}].":".$snpeff_hgvs;
							$snpeff_gene_name_list{$info_split[$info_split_index_array{"Gene_Name"}]}=1;
							$nb_infos++;
						};#if
						# Gene ID
						if ($info_split[$info_split_index_array{"Gene_ID"}] ne "") {
							$snpeff_hgvs=$info_split[$info_split_index_array{"Gene_ID"}].":".$snpeff_hgvs;
							$snpeff_gene_id_list{$info_split[$info_split_index_array{"Gene_ID"}]}=1;
							$nb_infos++;
						};#if
						# Protein HGVS annotation
						if ($info_split[$info_split_index_array{"HGVS.p"}] ne "") {
							$snpeff_hgvs.=":".$info_split[$info_split_index_array{"HGVS.p"}];
							$nb_infos++;
						};#
						# HGVS annotation
						$snpeff_hgvs_list{$snpeff_hgvs}=1;
						
						# Annotation
						if ($info_split[$info_split_index_array{"Annotation"}] ne "") {
							$snpeff_annotation_list{$info_split[$info_split_index_array{"Annotation"}]}=1;
						};#
						
						# Impact
						if ($info_split[$info_split_index_array{"Annotation_Impact"}] ne "") {
							$snpeff_impact_list{$info_split[$info_split_index_array{"Annotation_Impact"}]}=1;
						};#
						
						
					};#if
				};#foreach
				my $snpeff_hgvs=join(",",sort { length($b) <=> length($a) } keys %snpeff_hgvs_list);
				my $snpeff_gene_name=join(",",sort { length($b) <=> length($a) } keys %snpeff_gene_name_list);
				my $snpeff_gene_id=join(",",sort { length($b) <=> length($a) } keys %snpeff_gene_id_list);
				my $snpeff_annotation=join(",",sort { length($b) <=> length($a) } keys %snpeff_annotation_list);
				my $snpeff_impact=join(",",sort { length($b) <=> length($a) } keys %snpeff_impact_list);
				
				
				# output annotation
				
				if ($parameters{"snpeff_hgvs"}) {
					my $sep_snpeff_hgvs=",";
					if (!defined $annotation_output{$chr}{$pos}{$ref}{$alt}{"snpeff_hgvs"}
						|| $annotation_output{$chr}{$pos}{$ref}{$alt}{"snpeff_hgvs"} eq "") {
						$sep_snpeff_hgvs="";
					};#if
					$annotation_output{$chr}{$pos}{$ref}{$alt}{"snpeff_hgvs"}.=$sep_snpeff_hgvs.$snpeff_hgvs;
				
					my $sep_hgvs=",";
					if (!defined $annotation_output{$chr}{$pos}{$ref}{$alt}{"hgvs"}
						|| $annotation_output{$chr}{$pos}{$ref}{$alt}{"hgvs"} eq "") {
						$sep_hgvs="";
					} else {
						$annotation_output{$chr}{$pos}{$ref}{$alt}{"hgvs"} =~ s/,$//;
					};#if
					$annotation_output{$chr}{$pos}{$ref}{$alt}{"hgvs"}.=$sep_hgvs.$snpeff_hgvs;
					
					# HEADER Mandatory
					my $annotation_name="snpeff_hgvs";
					my $annotation_description="snpEff HGVS annotation";
					my $additional_description=""; #" [Release=$release;Date=$date;AnnotationType=$annotation_type]";
					$annotation_description =~ s/,//g;
					$additional_description =~ s/,//g;
					$vcf_header{"INFO"}{$annotation_name}{"Number"}=".";
					$vcf_header{"INFO"}{$annotation_name}{"Type"}="String";
					$vcf_header{"INFO"}{$annotation_name}{"Description"}="\"$annotation_description$additional_description\"";
					
					my $annotation_name="hgvs";
					my $annotation_description="snpEff HGVS annotation";
					my $additional_description=""; #" [Release=$release;Date=$date;AnnotationType=$annotation_type]";
					$annotation_description =~ s/,//g;
					$additional_description =~ s/,//g;
					$vcf_header{"INFO"}{$annotation_name}{"Number"}=".";
					$vcf_header{"INFO"}{$annotation_name}{"Type"}="String";
					$vcf_header{"INFO"}{$annotation_name}{"Description"}="\"$annotation_description$additional_description\"";
					
				};#if
				
				if ($parameters{"snpeff_gene_name"}) {
					#print "snpeff_gene_name $snpeff_gene_name OK\n" if $VERBOSE;
					my $sep_snpeff_gene_name=",";
					if (!defined $annotation_output{$chr}{$pos}{$ref}{$alt}{"snpeff_gene_name"}
						|| $annotation_output{$chr}{$pos}{$ref}{$alt}{"snpeff_gene_name"} eq "") {
						$sep_snpeff_gene_name="";
					};#if
					$annotation_output{$chr}{$pos}{$ref}{$alt}{"snpeff_gene_name"}.=$sep_snpeff_gene_name.$snpeff_gene_name;
					#print "snpeff_gene_name ".$annotation_output{$chr}{$pos}{$ref}{$alt}{"snpeff_gene_name"}." ADDED\n" if $VERBOSE;
				
					
					my $sep_symbol=",";
					@snpeff_symbol_split=split(/[\(,]/,$snpeff_gene_name);
					$snpeff_symbol=$snpeff_symbol_split[0];
					if (!defined $annotation_output{$chr}{$pos}{$ref}{$alt}{"symbol"}
						|| $annotation_output{$chr}{$pos}{$ref}{$alt}{"symbol"} eq "") {
						$sep_symbol="";
					} else {
						$annotation_output{$chr}{$pos}{$ref}{$alt}{"symbol"} =~ s/,$//;
					};#if
					$annotation_output{$chr}{$pos}{$ref}{$alt}{"symbol"}.=$sep_symbol.$snpeff_symbol;
					#print "symbol ".$annotation_output{$chr}{$pos}{$ref}{$alt}{"symbol"}." ADDED\n" if $VERBOSE;
					
					# HEADER Mandatory
					my $annotation_name="snpeff_gene_name";
					my $annotation_description="snpEff Gene Name annotation";
					my $additional_description=""; #" [Release=$release;Date=$date;AnnotationType=$annotation_type]";
					$annotation_description =~ s/,//g;
					$additional_description =~ s/,//g;
					$vcf_header{"INFO"}{$annotation_name}{"Number"}=".";
					$vcf_header{"INFO"}{$annotation_name}{"Type"}="String";
					$vcf_header{"INFO"}{$annotation_name}{"Description"}="\"$annotation_description$additional_description\"";
					
					my $annotation_name="symbol";
					my $annotation_description="snpEff Gene Name annotation";
					my $additional_description=""; #" [Release=$release;Date=$date;AnnotationType=$annotation_type]";
					$annotation_description =~ s/,//g;
					$additional_description =~ s/,//g;
					$vcf_header{"INFO"}{$annotation_name}{"Number"}=".";
					$vcf_header{"INFO"}{$annotation_name}{"Type"}="String";
					$vcf_header{"INFO"}{$annotation_name}{"Description"}="\"$annotation_description$additional_description\"";
					
					
				};#if
				
				if ($parameters{"snpeff_annotation"}) {
					#print "snpeff_annotation $snpeff_annotation OK\n" if $VERBOSE;
					my $sep_snpeff_annotation=",";
					if (!defined $annotation_output{$chr}{$pos}{$ref}{$alt}{"snpeff_annotation"}
						|| $annotation_output{$chr}{$pos}{$ref}{$alt}{"snpeff_annotation"} eq "") {
						$sep_snpeff_annotation="";
					};#if
					$annotation_output{$chr}{$pos}{$ref}{$alt}{"snpeff_annotation"}.=$sep_snpeff_annotation.$snpeff_annotation;
					#print "snpeff_annotation ".$annotation_output{$chr}{$pos}{$ref}{$alt}{"snpeff_annotation"}." ADDED\n" if $VERBOSE;
				
					my $sep_location=",";
					if (!defined $annotation_output{$chr}{$pos}{$ref}{$alt}{"location"}
						|| $annotation_output{$chr}{$pos}{$ref}{$alt}{"location"} eq "") {
						$sep_location="";
					} else {
						$annotation_output{$chr}{$pos}{$ref}{$alt}{"location"} =~ s/,$//;
					};#if
					$annotation_output{$chr}{$pos}{$ref}{$alt}{"location"}.=$sep_location.$snpeff_annotation;
					#print "location ".$annotation_output{$chr}{$pos}{$ref}{$alt}{"location"}." ADDED\n" if $VERBOSE;
					
					my $sep_location=",";
					if (!defined $annotation_output{$chr}{$pos}{$ref}{$alt}{"outcome"}
						|| $annotation_output{$chr}{$pos}{$ref}{$alt}{"outcome"} eq "") {
						$sep_location="";
					} else {
						$annotation_output{$chr}{$pos}{$ref}{$alt}{"outcome"} =~ s/,$//;
					};#if
					$annotation_output{$chr}{$pos}{$ref}{$alt}{"outcome"}.=$sep_location.$snpeff_annotation;
					#print "outcome ".$annotation_output{$chr}{$pos}{$ref}{$alt}{"outcome"}." ADDED\n" if $VERBOSE;
					
					# HEADER Mandatory
					my $annotation_name="snpeff_annotation";
					my $annotation_description="snpEff Annotation";
					my $additional_description=""; #" [Release=$release;Date=$date;AnnotationType=$annotation_type]";
					$annotation_description =~ s/,//g;
					$additional_description =~ s/,//g;
					$vcf_header{"INFO"}{$annotation_name}{"Number"}=".";
					$vcf_header{"INFO"}{$annotation_name}{"Type"}="String";
					$vcf_header{"INFO"}{$annotation_name}{"Description"}="\"$annotation_description$additional_description\"";
					
					my $annotation_name="location";
					my $annotation_description="snpEff Annotation";
					my $additional_description=""; #" [Release=$release;Date=$date;AnnotationType=$annotation_type]";
					$annotation_description =~ s/,//g;
					$additional_description =~ s/,//g;
					$vcf_header{"INFO"}{$annotation_name}{"Number"}=".";
					$vcf_header{"INFO"}{$annotation_name}{"Type"}="String";
					$vcf_header{"INFO"}{$annotation_name}{"Description"}="\"$annotation_description$additional_description\"";
					
					my $annotation_name="outcome";
					my $annotation_description="snpEff Annotation";
					my $additional_description=""; #" [Release=$release;Date=$date;AnnotationType=$annotation_type]";
					$annotation_description =~ s/,//g;
					$additional_description =~ s/,//g;
					$vcf_header{"INFO"}{$annotation_name}{"Number"}=".";
					$vcf_header{"INFO"}{$annotation_name}{"Type"}="String";
					$vcf_header{"INFO"}{$annotation_name}{"Description"}="\"$annotation_description$additional_description\"";
					
					
				};#if
				
				if ($parameters{"snpeff_impact"}) {
					#print "snpeff_impact $snpeff_impact OK\n" if $VERBOSE;
					my $sep_snpeff_impact=",";
					if (!defined $annotation_output{$chr}{$pos}{$ref}{$alt}{"snpeff_impact"}
						|| $annotation_output{$chr}{$pos}{$ref}{$alt}{"snpeff_impact"} eq "") {
						$sep_snpeff_impact="";
					};#if
					$annotation_output{$chr}{$pos}{$ref}{$alt}{"snpeff_impact"}.=$sep_snpeff_impact.$snpeff_impact;
					#print "snpeff_impact ".$annotation_output{$chr}{$pos}{$ref}{$alt}{"snpeff_impact"}." ADDED\n" if $VERBOSE;
				
					# HEADER Mandatory
					my $annotation_name="snpeff_impact";
					my $annotation_description="snpEff Impact annotation";
					my $additional_description=""; #" [Release=$release;Date=$date;AnnotationType=$annotation_type]";
					$annotation_description =~ s/,//g;
					$additional_description =~ s/,//g;
					$vcf_header{"INFO"}{$annotation_name}{"Number"}=".";
					$vcf_header{"INFO"}{$annotation_name}{"Type"}="String";
					$vcf_header{"INFO"}{$annotation_name}{"Description"}="\"$annotation_description$additional_description\"";
					
					
				};#if
				
				if ($parameters{"snpeff_split"}) {
					#print "snpeff_split $snpeff_split OK\n" if $VERBOSE;
					my $sep_snpeff_split=",";
					if (!defined $annotation_output{$chr}{$pos}{$ref}{$alt}{"snpeff_split"}
						|| $annotation_output{$chr}{$pos}{$ref}{$alt}{"snpeff_split"} eq "") {
						$sep_snpeff_split="";
					};#if
					$snpeff_split="$snpeff_gene_name|$snpeff_impact|$snpeff_annotation|$snpeff_hgvs";
					$annotation_output{$chr}{$pos}{$ref}{$alt}{"snpeff_split"}.=$sep_snpeff_split.$snpeff_split;
					#print "snpeff_split ".$annotation_output{$chr}{$pos}{$ref}{$alt}{"snpeff_split"}." ADDED\n" if $VERBOSE;
				
					# HEADER Mandatory
					my $annotation_name="snpeff_split";
					my $annotation_description="snpEff splitted annotations: snpEff_gene_name|snpEff_impact|snpEff_annotation|snpEff_hgvs";
					my $additional_description=""; #" [Release=$release;Date=$date;AnnotationType=$annotation_type]";
					$annotation_description =~ s/,//g;
					$additional_description =~ s/,//g;
					$vcf_header{"INFO"}{$annotation_name}{"Number"}=".";
					$vcf_header{"INFO"}{$annotation_name}{"Type"}="String";
					$vcf_header{"INFO"}{$annotation_name}{"Description"}="\"$annotation_description$additional_description\"";
					
					
				};#if
				


			};#if


		}
		}
		}
		}



		#my $annotation_name="ANN";
		#$vcf_header{"INFO"}{$annotation_name}{"Number"}=".";
		#$vcf_header{"INFO"}{$annotation_name}{"Type"}="String";
		#my $description=((trim($annotation_description) eq "")?"unknown":trim($annotation_description));
		#my $additional_description=" [Release=$release;Date=$date;AnnotationType=$annotation_type]";
		#$vcf_header{"INFO"}{$annotation_name}{"Description"}="\"$description$additional_description\"";


		#open(FILE_SNPEFF_LOF_INPUT, $output_file_snpeff_log) || die "Problem to open the file '$output_file_snpeff_log': $!";
		#while(<FILE_INPUT>) {
		#	print $_;
		#};#

		# Remove SNPEff tmp files
		if ($input_file_snpeff_stats ne "") {
			my @tmp_snpeff_stats_files=<$input_file_snpeff_stats*>;
			unlink(@tmp_snpeff_stats_files);
		}; # if
		if ($output_file_snpeff ne "") {
			my @tmp_snpeff_files=<$output_file_snpeff*>;
			unlink(@tmp_snpeff_files);
		}; # if
		if ($output_file_snpeff_log ne "") {
			my @tmp_snpeff_log_files=<$output_file_snpeff_log*>;
			unlink(@tmp_snpeff_log_files);
		}; # if

	};#if

};#if



#################
### Write VCF ###
#################


#open the files
open(FILE_OUTPUT, ">$output_file") || die "Problem to open the file '$output_file': $!";
open(FILE_INPUT, $input_file) || die "Problem to open the file '$input_file': $!";
$line=0;
$line_info_read=0;

#read the file
while(<FILE_INPUT>) {

	# init
	chomp; #delete \n character
	$line++;
	$line_content=$_;
	@line_content_split=split("\t",$line_content);

	if (substr($line_content,0,2) eq "##") { # VCF Header

		if ($line_content=~ /^##(.*)=<(.*)>$/) {
			#print "$_\n";
			#print "type=$1\n";
			my $type=$1;
			if ($type eq "INFO" || $type eq "FORMAT" || $type eq "contig" || $type eq "FILTER") {
				my $split_limit=4;
				if ($type eq "INFO") { $split_limit=4; };
				if ($type eq "FORMAT") { $split_limit=4; };
				if ($type eq "contig") { $split_limit=3; };
				if ($type eq "FILTER") { $split_limit=2; };
				my @infos_split=split(",",$2,$split_limit);
				my %infos;
				foreach my $info (@infos_split) {
					my @info_split=split("=",$info,2);
					$infos{$info_split[0]}=$info_split[1];
				};#foreach
				#print Dumper(\%infos) if $VERBOSE;
				if (!exists($vcf_header{$1}{$infos{"ID"}})) { # Updated!
					while ((my $info_var, my $info_val) = each(%infos)){
						$vcf_header{$1}{$infos{"ID"}}{$info_var}=$info_val if ($info_var ne "ID");
					};#while
				};#if
			};#if

		} else {
			# just info, such as ##fileformat=VCFv4.1
			$vcf_header_information.="$line_content\n";
		};#if



		#%vcf_header = map{split /,/, $_}(split /=</, $line_content);

		#print "########HEADER $line_content\n";
		#check the line in order to insert the description of the new info into the VCF INFO column
		#if (substr($line_content,0,7) eq "##INFO=") {
		#	$line_info_read=1;
		#} elsif ($line_info_read==1) {
		#	print FILE_INPUT_annotated "$vcf_header";
		#	$line_info_read=0;
		#};#if
		#keep the line
		#print FILE_OUTPUT "$line_content\n";
	} elsif (substr($line_content,0,6) eq "#CHROM") { # case of variant header "#CHROM"
		#print "########CHROM... $vcf_header_INFOS\n$line_content\n";
		#print "TEST";
		print FILE_OUTPUT $vcf_header_information;

		foreach my $type (sort { $vcf_header{$a} <=> $vcf_header{$b} or $a cmp $b } keys %vcf_header) {
			foreach my $ID (sort { $vcf_header{$type}{$a} <=> $vcf_header{$type}{$b} or $a cmp $b } keys %{$vcf_header{$type}}) {
				my $line="##$type=<ID=$ID";
				while ((my $info_var, my $info_val) = each(%{$vcf_header{$type}{$ID}})){
					if (trim($info_val) ne "") {
						$line.=",$info_var=$info_val";
					};#if
				};#while
				$line.=">\n";
				print FILE_OUTPUT $line;
			};#foreach
		};#foreach

		# Add snpEff header
		if ($parameters{"snpeff"}) {
			#print "# SNPEFF HEADER: $annotation_snpeff_output_header\n" if $VERBOSE;
			if ($annotation_snpeff_output_header ne "") {
				print FILE_OUTPUT "$annotation_snpeff_output_header\n";
			};#if
		};#if

		#print FILE_OUTPUT $vcf_header_INFOS;

		print FILE_OUTPUT "$line_content\n";
	} else {
		my $chr=$line_content_split[0];
		my $pos=$line_content_split[1];
		my $ref=$line_content_split[3];
		my $alt=$line_content_split[4];
		my $annotation_input_line=$line_content_split[7];
		$annotation_input_line="" if (trim($annotation_input_line) eq ".");

		#if ($parameters{"snpeff"}) {
		#	print "# {$chr}{$pos}{$ref}{$alt}{$var} ".$annotation_snpeff_output{$chr}{$pos}{$ref}{$alt}{"INFOS"}{"ANN"} if $VERBOSE;
		#};#if

		#%{$annotation_input{$chr}{$pos}{$ref}{$alt}} = map{split /=/, $_}(split /;/, $annotation_input_line);
		foreach my $ann (split(/;/,$annotation_input_line)) {
			(my $var, my $val) = split(/=/,$ann);
			#print "$ann\t\t'$var'='$val'\n" if $DEBUG;
			$annotation_input{$chr}{$pos}{$ref}{$alt}{$var}=$val;
		};#foreach
		#print Dumper(\%annotation_input) if $DEBUG;
		print "### {$chr}{$pos}{$ref}{$alt}\n" if $DEBUG;

		while ((my $annotation_name, my $annotation_result) = each(%{$annotation_output{$chr}{$pos}{$ref}{$alt}})){
		#while ((my $annotation_name, my $annotation_result) = each(%{$annotation_input{$chr}{$pos}{$ref}{$alt}})){
			#my $info_field=$annotations{$chrom}{$pos}{$ref}{$alt};


			if ($force) {
				$annotation_input{$chr}{$pos}{$ref}{$alt}{$annotation_name}=$annotation_result;
			};#if

			my $sep_default=",";
			#print "#       TEST:\t".trim($annotation_input{$chr}{$pos}{$ref}{$alt}{$annotation_name})." ??? ".$result."\n";
			if (!exists $annotation_input{$chr}{$pos}{$ref}{$alt}{$annotation_name}
				|| trim($annotation_input{$chr}{$pos}{$ref}{$alt}{$annotation_name}) eq ""
				|| trim($annotation_input{$chr}{$pos}{$ref}{$alt}{$annotation_name}) eq $sep_default.trim($annotation_result)) {
				$annotation_input{$chr}{$pos}{$ref}{$alt}{$annotation_name}=$annotation_result;
			} elsif (trim($annotation_input{$chr}{$pos}{$ref}{$alt}{$annotation_name}) ne trim($annotation_result)) {
				$annotation_input{$chr}{$pos}{$ref}{$alt}{$annotation_name}.=$sep_default.$annotation_result;
			};#if
			#print "#       $chr:$pos:$ref:$alt\t".$annotation_output{$chr}{$pos}{$ref}{$alt}{$annotation_name}."\n";
			#$output_results.="#       $chr:$pos:$ref:$alt\t".$annotation_output{$chr}{$pos}{$ref}{$alt}{$annotation_name}."\n";

			#print "$annotation_name=$annotation_result\n";
			#if (0 && (!exists $annotation_input{$chr}{$pos}{$ref}{$alt}{$annotation_name}
			#	|| $annotation_result ne $annotation_input{$chr}{$pos}{$ref}{$alt}{$annotation_name})) {
			#	#$line_content_split[7].=";$annotation_name=$annotation_result";
			#	my $original_result=trim($annotation_input{$chr}{$pos}{$ref}{$alt}{$annotation_name});
			#	my $sep=((trim($annotation_input{$chr}{$pos}{$ref}{$alt}{$annotation_name}) eq "")?"":",");
			#	$annotation_input{$chr}{$pos}{$ref}{$alt}{$annotation_name}="$original_result$sep$annotation_result";
			#};#if
		};#while
		#
		#my $str = join(", ", map { "$_ X $hash{$_}" } keys %hash);
		#my $variant_annotation=join(";", map { "$_=$annotation_input{$chr}{$pos}{$ref}{$alt}{$_}" } keys %{$annotation_input{$chr}{$pos}{$ref}{$alt}});
		#print Dumper(\%annotation_output) if $DEBUG;

		my $variant_annotation="";
		while ((my $annotation_name, my $annotation_result) = each(%{$annotation_input{$chr}{$pos}{$ref}{$alt}})){
			my $sep=($variant_annotation ne "")?";":"";
			if (trim($annotation_result) eq "" && $vcf_header{"INFO"}{$annotation_name}{"Number"} eq "0") {
				$variant_annotation.="$sep$annotation_name";
			} elsif (!defined $annotation_result || trim($annotation_result) eq "") {
				$variant_annotation.="$sep$annotation_name";
			} else {
				$variant_annotation.="$sep$annotation_name=$annotation_result";
			};#if
		};#if

		print "$variant_annotation\n" if $DEBUG;

		$variant_annotation="." if (trim($variant_annotation) eq "");
		$line_content_split[7]=$variant_annotation;
		$variant_line_join=join("\t",@line_content_split);
		print FILE_OUTPUT "$variant_line_join\n";
	};#if

	#print Dumper(\%vcf_header) if $VERBOSE;


};#while
close(FILE_INPUT);
close(FILE_OUTPUT);


#my %vcf_hash=read_vcf($output_file,"header");
#print Dumper(\%vcf_hash) if $VERBOSE;


## PostProcess
################

$header.="##\n";


## OUTPUT
###########

# Header
print $header;
print $debug;
print $verbose;
print $output;


__END__
