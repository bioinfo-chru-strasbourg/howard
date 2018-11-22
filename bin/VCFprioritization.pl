#!/usr/bin/perl
############################
# VCF prioritization       #
# Author: Antony Le Béchec #
# Copyright: IRC           #
############################

## Main Information
#####################

our %information = ( #
	'script'	=>  	basename($0),		# Script
	'release'	=>  	"0.9.5.1b",		# Release
	'date'		=>  	"20180823",		# Release parameter
	'author'	=>  	"Antony Le Béchec",	# Author
	'copyright'	=>  	"IRC",			# Copyright
	'licence'	=>  	"GNU-GPL",		# Licence
);

## Release Notes
##################
# 20180823-0.9.5.1b: Fix bug with INFO/. for empty INFO field


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
use lib dirname (__FILE__);	# Add lib in the same folder
use Scalar::Util qw(looks_like_number);
use Digest::MD5 qw(md5_hex);

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

IRC - GNU GPL License

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

=item B<--config_filter|config_filter_file=<file>>

Configuration file for filter/prioritization parameters (default 'config.prioritization.ini' or 'config.filter.ini').

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

=item B<--filter|prioritization=string>

Filter profile, defined in the 'config_annotation' file. 'ALL' for all filters defined in the 'config_annotation' file. If the filter doesn't exist, 'ALL' will be used (if exist).

=item B<--hard!>

Remove not PASS variant in default priorisation filter (PZFlag)

=item B<--pzfields=<string>>

List of prioritisation information to show

Default: 'PZScore,PZFlag'

Format: 'field1,field2...'

Example: 'PZScore,PZFlag,PZComment,PZInfos'

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



## Input file
my $input_file;
if (-e $parameters{"input"}) {
	$input_file=$parameters{"input"};
	$header.="#[INFO] Input VCF file: $input_file\n";
} else {
	print "#[ERROR] input file '".$parameters{"input"}."' DOES NOT exist\n";
	pod2usage();
	exit 1;
};#if

## Output file
my $output_file;
my $output_filename_pattern="prioritized";
if (-e $parameters{"output"} && 0) {
	print "#[ERROR] output file '".$parameters{"output"}."' DOES exist\n";
	pod2usage();
	exit 1;
} else {
	$output_file=$parameters{"output"};
	if (trim($output_file) eq "") {
		$output_file=$input_file;
		$output_file =~ s/\.vcf$/\.$output_filename_pattern\.vcf/g; #/\.vcf$/\.output\.vcf/;
	};#if
	$header.="#[INFO] Output VCF file: $output_file\n";
};#if


# Filter

my %annotation_filters;
while ((my $filter_name, my $filter_criteria) = each(%config_filter)){
	#print "$filter_name\n" if $DEBUG;
	while ((my $criterion, my $criterion_values) = each(%{$filter_criteria})){
		my %criterion_list;
		#print "\t$criterion\n" if $DEBUG;
		if (trim($criterion) =~ /.*\[\]$/) { # $criterion_values is a hash
			while ((my $criterion_value, my $criterion_value_infos) = each(%{$criterion_values})){
				$criterion_list{$criterion_value}=$criterion_value;
			};#while
		} else { # $criterion_values is a string
			$criterion_list{$criterion_values}=$criterion_values;
		};#if
		# Clean criterion
		$criterion=~s/\[\]$//g;
		$criterion=lc($criterion);
		# Explode criterions values
		while ((my $one_filter, my $one_filter_infos) = each(%criterion_list)){
			@one_filter_array=split(":",$one_filter);
			#print "\t\t@one_filter_array\n" if $DEBUG;
			my $criterion_source;
			my $criterion_value;
			my $criterion_pond;
			my $criterion_comment;
			# Value of the filter
			my $criterion_value=trim($one_filter_array[0]);
			$criterion_value=~s/^"|"$//g;
			# Ponderation of the filter
			my $criterion_pond=trim($one_filter_array[1]);
			$criterion_pond=~s/^"|"$//g;
			# Ponderation check
			if (!defined $criterion_pond) {
			    $criterion_pond=1;
			} else {
			    if (trim($criterion_pond) eq "+") {
				$criterion_pond=1;
			    };#if
			    if (trim($criterion_pond) eq "-") {
				$criterion_pond=-1;
			    };#if
			};#if
			# Comment of the filter
			my $criterion_comment=trim($one_filter_array[2]);
			$criterion_comment=~s/^"|"$//g;
			# Add filter
			#$annotation_filters{$filter_name}{$criterion}{$criterion_value}=$criterion_pond;
			$annotation_filters{$filter_name}{$criterion}{$criterion_value}{$criterion_pond}=$criterion_comment;
		};#while


	};#while
};#while
#print Dumper(\%parameters);
# filter input
#print "FILTER=".$parameters{"filter"};
my %filter_list;
@filter_list_input=split(",",$parameters{"filter"});
#print @filter_list_input;

foreach my $one_filter (@filter_list_input) {
	if (exists $annotation_filters{$one_filter}) { #  ||
		$filter_list{$one_filter}=1;
		if (!defined $filter) {
			$filter=$one_filter; # First filter is the default filter
		};#if
	};#if
	if (lc($one_filter) eq "all") {
		while ((my $one_filter_defined, my $one_filter_defined_infos) = each(%annotation_filters)){
			$filter_list{$one_filter_defined}=1;
		};#while
	};#if

};#foreach
if (!defined $filter) {
	$filter="default"; # if not defined by the list, default
};#if

$header.="## Filters: ".join( ", " ,(map { $_ } keys %filter_list))."\n";
$header.="## Filter by default: $filter\n";
#print "# Filters: ".join( ", " ,(map { $_ } keys %filter_list))."\n";# ;

# HARD
my $hard_filter=$parameters{"hard"};


# PZFields
#@pzfilter_list_input=split(",",$parameters{"pzfields"});
my %pzfilter_list_input = map { $_ => 1 } split(",",$parameters{"pzfields"});


## DEBUG
##########

#print Dumper(\%pzfilter_list_input) ; #if $DEBUG;
#print Dumper(\%config_annotation) if $DEBUG;
#print Dumper(\%config_filter) if $DEBUG;
#print Dumper(\%annotation_filters) if $DEBUG;


## MAIN
#########

# Variables
$output="";


my %variants;
my %annotations;

if ($VERBOSE) {
	print "#[INFO] Input file: $input_file\n";
	print "#[INFO] Output file: $output_file\n";
	#print "# Output format: @output_formats\n";
	print "#[INFO] Prioritization: $filter\n";
	#print "# Filters: ".(map { $_ } keys %filter_list)."\n";# ;
	print "#\n";
};#if

# READ VCF FILE
my %annotation_output=read_vcf($parameters{"input"});


if ($DEBUG) {
	#print Dumper(\%annotation_output) if $DEBUG;
	#print Dumper(\%annotations) if $DEBUG;filter_list
	#print Dumper(\%filter_list) if $DEBUG;
};#DEBUG

my $PZHeader;


sub prioritize {
	# Update prioritization values
	my $prioritization_flag=$_[0];
	my $prioritization_score=$_[1];
	my $prioritization_comment=$_[2];
	my $prioritization_infos=$_[3];
	my $variant_annotations_name=$_[4];
	my $annotation_name=$_[5];
	my $annotation_value=$_[6];
	my $filter_annotation=$_[7];
	my $filter_value=$_[8];
	my $filter_value_comment=$_[9];

	#$prioritization_flag,$prioritization_score,$prioritization_comment,$prioritization_infos,
	#print "$variant_annotations_name,$annotation_name,$annotation_value,$filter_annotation,$filter_value,$filter_value_comment\n" if $DEBUG;

	#print "START: $prioritization_flag,$prioritization_score,$prioritization_comment,$prioritization_infos\n" if $DEBUG;
	#print "   VALUES: @_,\n" if $DEBUG;

	if (looks_like_number($annotation_value)) { # is a number
		#print "\t\t\tprocessing 'score' filter\n" if $DEBUG;

		my $filter_value_numeric = $filter_annotation;
		$filter_value_numeric =~ s/[<>=]//g;
		#my $filter_value_numeric += 0; # integer cast
		$filter_value_comparison = $filter_annotation;
		$filter_value_comparison =~ s/$filter_value_numeric//g;
		#print "            filter_value_numeric $filter_value_numeric, filter_value_comparison $filter_value_comparison\n" if $DEBUG && $annotation_name eq "FILTER";


		if (	looks_like_number($filter_value_numeric)
			&&
			(
			(($filter_value_comparison eq ">=" || $filter_value_comparison eq ">")
				&& $annotation_value>=$filter_value_numeric)
			|| (($filter_value_comparison eq "<=" || $filter_value_comparison eq "<")
				&& $annotation_value<=$filter_value_numeric)
			|| ($annotation_value==$filter_value_numeric)
			)
		) {
			#print "NUMBER! \n" if $DEBUG;
			# PONDERATION (SCORE and FLAG)
			if (lc(trim($filter_value)) =~ /^p/) {
				$prioritization_flag="PASS";
			} elsif (lc(trim($filter_value)) =~ /^f/) {
				if ($prioritization_flag ne "PASS") { # Filtered only if not 'PASS'
					$prioritization_flag="FILTERED";
				};#if
			} elsif (looks_like_number($filter_value)) {
				$prioritization_score+=$filter_value;
			#} elsif ($filter_value =~ /^COMMENT/) {
			#	print "$comment\n" if $DEBUG;
			#	my $comment=$filter_value; $comment=~ s/^COMMENT//g; $comment=~ s/VALUE/$annotation_value/g;
			#	print "$comment\n" if $DEBUG;
			#	$prioritization_comment.="$comment,";
			} else {

			};#if
			# COMMENT
			if (defined $filter_value_comment && trim($filter_value_comment) ne "") {
				my $comment=$filter_value_comment;
				$comment=~ s/VALUE/$annotation_value/g; #$filter_annotation
				$comment=~ s/ANNOTATION/$annotation_name/g; #
				$comment=~ s/FILTER/$filter_annotation/g; #$filter_annotation
				#print "$prioritization_comment =~ /$comment/\n" if $DEBUG;
				if ($prioritization_comment =~ /$comment/) {
					#print "   => not add\n" if $DEBUG;
				} else {
					$prioritization_comment.=((trim($prioritization_comment) ne "")?". ":"").$comment;
					#print "   => ADD\n" if $DEBUG;
				};#if
			};#if

			# INFOS
			$prioritization_infos.="$variant_annotations_name:$annotation_name:$filter_annotation:$filter_value,";
		} else {

		};#if


	} else { # is not a number
		#print "\t\t\tprocessing 'annotation' filter\n" if $DEBUG;

		my $filter_annotation_copy=lc($filter_annotation);
		$filter_value_string = $filter_annotation_copy;
		$filter_value_string =~ s/[!<>=]//g;
		$filter_value_comparison = $filter_annotation_copy;
		$filter_value_comparison =~ s/$filter_value_string//g;
		$filter_value_comparison = substr $filter_annotation_copy, 0, (length($filter_annotation)-length($filter_value_string));
		#print "\t\t\tfilter_value_comparison=$filter_value_comparison\tfilter_value_string=$filter_value_string\tannotation_value=$annotation_value\t\t\t\t\t\t$filter_value_comparison '$filter_value_string' against '$annotation_value'\n" if $DEBUG && $annotation_name eq "FILTER";
		#print "TEST"; print (grep(/^$filter_value_string$/,lc($annotation_value))!=""); print "\n";
		if ( $filter_value_string ne ""
			&& (
			($filter_value_comparison eq "=" && (grep(/^$filter_value_string$/,lc($annotation_value))!="")) # Exact match
			|| ($filter_value_comparison eq "!=" && (grep(/^$filter_value_string$/,lc($annotation_value))=="")) # Exact match negation
			|| ($filter_value_comparison eq "!" && (grep(/^.*$filter_value_string.*$/,lc($annotation_value))=="")) # negation of something
			|| ($filter_value_comparison eq "" && (grep(/^.*$filter_annotation_copy.*$/,lc($annotation_value))!=""))
			)
			) {
			# PONDERATION (SCORE and FLAG)
			#print "\t\t\tMATCH!!!\n" if $DEBUG;
			#if (lc(trim($filter_value)) eq "p") {
			if (lc(trim($filter_value)) =~ /^p/) {
				#print "\t\t\t\tP\n" if $DEBUG;
				$prioritization_flag="PASS";
			#} elsif (lc(trim($filter_value)) eq "f") {
			} elsif (lc(trim($filter_value)) =~ /^f/) {
				#print "\t\t\t\tF\n" if $DEBUG;
				if ($prioritization_flag ne "PASS") { # Filtered only if not 'PASS'
					$prioritization_flag="FILTERED";
				};#if
			} elsif (looks_like_number($filter_value)) {
				#print "\t\t\t\tS\n" if $DEBUG;
				$prioritization_score+=$filter_value;
			#} elsif ($filter_value =~ /^COMMENT/) {
			#	my $comment=$filter_value; $comment=~ s/^COMMENT//g; $comment=~ s/VALUE/$annotation_value/g;
			#	$prioritization_comment.="$comment,";
			} else {

			};#if
			# COMMENT
			if (defined $filter_value_comment && trim($filter_value_comment) ne "") {
				my $comment=$filter_value_comment;
				$comment=~ s/VALUE/$annotation_value/g;
				$comment=~ s/ANNOTATION/$annotation_name/g; #
				$comment=~ s/FILTER/$filter_annotation/g; #$filter_annotation
				if ($prioritization_comment =~ /$comment/) {

				} else {
					$prioritization_comment.=((trim($prioritization_comment) ne "")?". ":"").$comment;
				};#if
			};#if
			# INFOS
			$prioritization_infos.="$variant_annotations_name:$annotation_name:$filter_annotation:$filter_value,";
		};#if

	};#if

	#print "END: $prioritization_flag,$prioritization_score,$prioritization_comment,$prioritization_infos\n" if $DEBUG;
	return ($prioritization_flag,$prioritization_score,$prioritization_comment,$prioritization_infos);

}


my %annotation_output2=%annotation_output;
# READ HASH Variants
while ( my ($chr, $poss) = each(%annotation_output2) ) {
while ( my ($pos, $refs) = each(%{$poss}) ) {
while ( my ($ref, $alts) = each(%{$refs}) ) {
while ( my ($alt, $variant_values) = each(%{$alts}) ) {

	#print "$chr:$pos|$ref>$alt\n" if $DEBUG || $VERBOSE;

	## List of annotations
	my %variant_annotations_to_check;
	my $variant_annotations_to_check;
	#print Dumper($variant_values) if $DEBUG;
	#print Dumper($variant_annotations_to_check) if $DEBUG;

	my %variant_annotations_list;
	# INFOS
	while ( my ($a, $b) = each(%{$$variant_values{"INFOS"}}) ) {
		$$variant_annotations_to_check{$a}=$b;
		#$$variant_annotations_list{"INFOS"}{$a}=$b;
		#print "     $a}=$b\n" if $DEBUG;
	};#while
	foreach my $variant_value (("CHROM","CHROM_NUM","POS","ID","REF","ALT","QUAL","FILTER",)) {
		#print"$variant_value=".$$variant_values{$variant_value}."\n" if $DEBUG;
		$$variant_annotations_to_check{$variant_value}=$$variant_values{$variant_value};
		#$$variant_annotations_list{"INFOS"}{$variant_value}=$$variant_values{$variant_value};
	};#foreach


	$variant_annotations_list{"INFOS"}=$variant_annotations_to_check;

	# SAMPLES
	#print Dumper(\%variant_annotations_list) if $DEBUG;
	my $samples=$$variant_values{"SAMPLES"};
	while ( my ($sample_name, $sample_quality) = each(%$samples) ) {
		#print "$sample_name\n" if $DEBUG;
		$variant_annotations_list{"SAMPLE.$sample_name"}=$sample_quality;
	};#foreach

	# VARIANT

	#print Dumper(\%variant_annotations_list) if $DEBUG && $pos=="106180796";

	# Filter exists on this annotation
	while ( my ($one_filter, $one_filter_infos) = each(%filter_list) ) {
	#print "$one_filter\n" if $DEBUG;

		# Filter Extension
		my $PZExtension="-$one_filter";
		if ($one_filter eq $filter) { # default filter
			$PZExtension="";
		};#if

		# if annotation exist and not forces
		if (!$force && defined $$variant_values{"INFOS"}{"PZFlag$PZExtension"}) {
			next;
		};#if

		my $prioritization_flag="";
		my $prioritization_score=0;
		my $prioritization_comment="";
		my $prioritization_infos="";


		#foreach my $variant_annotations (@variant_annotations_list) {
		while ( my ($variant_annotations_name, $variant_annotations) = each(%variant_annotations_list) ) {

			# For each variant annotation
			while ( my ($annotation_name, $annotation_value) = each(%{$variant_annotations}) ) {
			#print "\t$annotation_name=$annotation_value\n" if $DEBUG ;#&& $pos=="106180796";

				if (exists $annotation_filters{$one_filter}{lc($annotation_name)}) {
					while (my ($filter_annotation, $filter_values) = each(%{$annotation_filters{$one_filter}{lc($annotation_name)}})) {
					while (my ($filter_value, $filter_value_comment) = each(%{$filter_values})) {
						#print "\t? $filter_annotation $filter_values\n" if $DEBUG;
						($prioritization_flag,
						$prioritization_score,
						$prioritization_comment,
						$prioritization_infos)=prioritize($prioritization_flag,$prioritization_score,$prioritization_comment,$prioritization_infos,
											$variant_annotations_name,$annotation_name,$annotation_value,$filter_annotation,$filter_value,$filter_value_comment);

					};#while
					};#while

				};#fi
			};#while

		};#foreach

		# PASS if not FILTERED or other flag
		if (trim($prioritization_flag) eq "") {
			$prioritization_flag="PASS";
		};#if

		
		#print Dumper(\%pzfilter_list_input) ; #if $DEBUG;

		# PZScore
		if (defined $pzfilter_list_input{"PZScore"}) {

			# If Score not calculated (no filter found)
			$prioritization_score=$prioritization_score+0;
			# Annotations summary
			$annotations{"PZScore$PZExtension"}{$prioritization_score}++;
			# Add Filter infos into variant annotation
			$annotation_output{$chr}{$pos}{$ref}{$alt}{"INFOS"}{"PZScore$PZExtension"}=$prioritization_score;
			# Header
			if (!exists $VCF_header{$header_type}{"PZScore$PZExtension"}) {
				$PZHeader.="##INFO=<ID=PZScore$PZExtension,Number=.,Type=String,Description=\"Prioritization Score for filter '$one_filter', depending on the ponderation of each prioritization criterion\">\n";
				$VCF_header{$header_type}{"PZScore$PZExtension"}=1;
			};#if
			
		};#if

		# PZScore
		if (defined $pzfilter_list_input{"PZFlag"}) {

			# If Score not calculated (no filter found)
			$prioritization_flag=trim($prioritization_flag);
			$prioritization_flag="PASS" if $prioritization_flag eq "";
			# Annotations summary
			$annotations{"PZFlag$PZExtension"}{$prioritization_flag}++;
			# Add Filter infos into variant annotation
			$annotation_output{$chr}{$pos}{$ref}{$alt}{"INFOS"}{"PZFlag$PZExtension"}=$prioritization_flag;
			# Header
			if (!exists $VCF_header{$header_type}{"PZFlag$PZExtension"}) {
				$PZHeader.="##INFO=<ID=PZFlag$PZExtension,Number=.,Type=String,Description=\"Prioritization Flag for filter '$one_filter', either 'PASS' or 'FILTERED'\">\n";
				$VCF_header{$header_type}{"PZFlag$PZExtension"}=1;
			};#if
			
		};#if

		# PZComment
		if (defined $pzfilter_list_input{"PZComment"}) {
			# If Score not calculated (no filter found)
			#if (trim($prioritization_comment."")=="") {
			#	$prioritization_comment="NA";
			#};#if
			$prioritization_comment="NA" if trim($prioritization_comment) eq "";
			#$prioritization_comment=trim($prioritization_comment);
			#$prioritization_comment="";
			#$prioritization_comment="NA" if $prioritization_comment eq "";
			# Annotations summary
			$annotations{"PZComment$PZExtension"}{$prioritization_comment}++;
			# Add Filter infos into variant annotation
			$annotation_output{$chr}{$pos}{$ref}{$alt}{"INFOS"}{"PZComment$PZExtension"}=$prioritization_comment;
			# Header
			if (!exists $VCF_header{$header_type}{"PZComment$PZExtension"}) {
				$PZHeader.="##INFO=<ID=PZComment$PZExtension,Number=.,Type=String,Description=\"Prioritization Comment for filter '$one_filter', auto-comment the annotations\">\n";
				$VCF_header{$header_type}{"PZComment$PZExtension"}=1;
			};#if
			
		};#if

		# PZComment
		if (defined $pzfilter_list_input{"PZInfos"}) {

			# If Score not calculated (no filter found)
			$prioritization_infos="NA" if trim($prioritization_infos) eq "";
			# Annotations summary
			$annotations{"PZInfos$PZExtension"}{$prioritization_infos}++;
			# Add Filter infos into variant annotation
			$annotation_output{$chr}{$pos}{$ref}{$alt}{"INFOS"}{"PZInfos$PZExtension"}=$prioritization_infos;
			# Header
			if (!exists $VCF_header{$header_type}{"PZInfos$PZExtension"}) {
				$PZHeader.="##INFO=<ID=PZInfos$PZExtension,Number=.,Type=String,Description=\"Prioritization Infos for filter '$one_filter', describing the PZScore and PZFlag calculation\">\n";
				$VCF_header{$header_type}{"PZInfos$PZExtension"}=1;
			};#if
			
		};#if



	};#while


}
}
}
}

#print Dumper(\%annotation_output) if $DEBUG;

## Write VCF
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

	if (substr($line_content,0,1) eq "#") { # VCF Header
		#check the line in order to insert the description of the new info into the VCF INFO column
		if (substr($line_content,0,7) eq "##INFO=") {
			$line_info_read=1;
		} elsif ($line_info_read==1) {
			#print FILE_INPUT_annotated "$vcf_header";
			$line_info_read=0;
		};#if
		if (substr($line_content,0,6) eq "#CHROM") {
			# Add Header
			print FILE_OUTPUT $PZHeader;
		};#if
		#keep the line
		print FILE_OUTPUT "$line_content\n";
	} else {
		my $chr=$line_content_split[0];
		my $pos=$line_content_split[1];
		my $ref=$line_content_split[3];
		my $alt=$line_content_split[4];
		my $annotation_input_line=$line_content_split[7];

		#### ERROR !!! ####
		#%{$annotation_input{$chr}{$pos}{$ref}{$alt}{"INFOS"}} = map{split /=/, $_, 2}(split /;/, $annotation_input_line);
		

		#if ($DEBUG) {
		#my @test=split /;/, $annotation_input_line;
		foreach my $varval (split /;/, $annotation_input_line) {
			#print "$varval\n" if $DEBUG;
			(my $var, my $val)=(split /=/, $varval, 2);
			
			#print "   $var = $val\n" if $DEBUG;
			if ($var ne ".") {
				$annotation_input{$chr}{$pos}{$ref}{$alt}{"INFOS"}{$var}=$val;
			};#if
			#foreach my $t2 (split /=/, $varval, 2) {
			#	print "   $t2\n" if $DEBUG;
			#	
			#};#foreach
		};#foreach
		#print "@test\n" if $DEBUG;
		#print Dumper(\%{$annotation_input{$chr}{$pos}{$ref}{$alt}}) if $DEBUG;
			
		#};#if

		while ((my $annotation_name, my $annotation_result) = each(%{$annotation_output{$chr}{$pos}{$ref}{$alt}{"INFOS"}})){

			if ($force) {
				$annotation_input{$chr}{$pos}{$ref}{$alt}{"INFOS"}{$annotation_name}=$annotation_result;
			};#if

			if (!exists $annotation_input{$chr}{$pos}{$ref}{$alt}{"INFOS"}{$annotation_name}
				|| $annotation_result ne $annotation_input{$chr}{$pos}{$ref}{$alt}{"INFOS"}{$annotation_name}) {
				my $original_result=trim($annotation_input{$chr}{$pos}{$ref}{$alt}{"INFOS"}{$annotation_name});
				my $sep=((trim($annotation_input{$chr}{$pos}{$ref}{$alt}{"INFOS"}{$annotation_name}) eq "")?"":",");
				$annotation_input{$chr}{$pos}{$ref}{$alt}{"INFOS"}{$annotation_name}="$original_result$sep$annotation_result";
			};#if
		};#while

		print Dumper(\%{$annotation_input{$chr}{$pos}{$ref}{$alt}}) if $DEBUG;

		# Hard filtering
		my $hard_filter_remove=0;
		if ($hard_filter) {
			print "HARD ".$annotation_input{$chr}{$pos}{$ref}{$alt}{'INFOS'}{"PZFlag"} if $DEBUG;
			#if ($annotation_input{$chr}{$pos}{$ref}{$alt}{'INFOS'}{"PZFlag"} ne "PASS") {
			if ($annotation_input{$chr}{$pos}{$ref}{$alt}{'INFOS'}{"PZFlag"} =~ /.*PASS.*/) {
			} else {
				$hard_filter_remove=1;
			};#if
			print " $hard_filter_remove " if $DEBUG;
			print "\n" if $DEBUG;
		};#if

		# OUTPUT
		if (!$hard_filter_remove) {
			#my $variant_annotation=join(";", map { "$_=$annotation_input{$chr}{$pos}{$ref}{$alt}{'INFOS'}{$_}" } keys %{$annotation_input{$chr}{$pos}{$ref}{$alt}{"INFOS"}	});
			
			my $variant_annotation="";
			my @varval_INFOS=keys(%{$annotation_input{$chr}{$pos}{$ref}{$alt}{"INFOS"}}) if $DEBUG;
			#print "@varval_INFOS\n" if $DEBUG;
			foreach my $var (keys(%{$annotation_input{$chr}{$pos}{$ref}{$alt}{"INFOS"}})) {
				my $val=$annotation_input{$chr}{$pos}{$ref}{$alt}{'INFOS'}{$var};
				#print "   $var = $val\n" if $DEBUG;
				my $sep=""; $sep=";" if $variant_annotation ne "";
				if (defined $val) {
					$variant_annotation.="$sep$var=$val";
				} else {
					$variant_annotation.="$sep$var";
				};#if
				
			};#foreach
			#foreach my $varval (split /;/, @varval_INFOS ) {
			#	print "$varval\n" if $DEBUG;
			#	(my $var, my $val)=(split /=/, $varval, 2);
			#
			#	#print "   $var = $val\n" if $DEBUG;
			#	$annotation_input{$chr}{$pos}{$ref}{$alt}{"INFOS"}{$var}=$val;
			#	#foreach my $t2 (split /=/, $varval, 2) {
			#	#	print "   $t2\n" if $DEBUG;
			#	#	
			#	#};#foreach
			#};#foreach
			print "$variant_annotation\n" if $DEBUG;


			$line_content_split[7]=$variant_annotation;
			$variant_line_join=join("\t",@line_content_split);
			print FILE_OUTPUT "$variant_line_join\n";
		};#if
	};#if

};#while
close(FILE_INPUT);
close(FILE_OUTPUT);


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
