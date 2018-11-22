#!/usr/bin/perl
############################
# VCF translation          #
# Author: Antony Le Béchec #
# Copyright: IRC           #
############################

## Main Information
#####################

our %information = ( #
	'script'	=>  	basename($0),		# Script
	'release'	=>  	"0.9.2.1b",		# Release
	'date'		=>  	"20180918",		# Release parameter
	'author'	=>  	"Antony Le Béchec",	# Author
	'copyright'	=>  	"IRC",			# Copyright
	'licence'	=>  	"GNU-GPL",		# Licence
);

## Release Notes
##################
# 20180918-0.9.2.1b: Fix format input, blank lines in vcf


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

=item B<--config_annotation|config_annotation_file=<file>>

Configuration file for annotation parameters (default 'config.annotation.ini').

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

=item B<--translation|format=<string>>

Output format, such as 'tab' (default 'tab').

=item B<--fields=<string>>

List of annotation from INFOS VCF field to include in the output file

Default: 'ALL'

Format: 'annotation1,annotation2,annotation3,...'

Example: 'PZScore,PZFlag,PZComment,Symbol,hgvs,location,outcome,AlleleFrequency,AD,DP,AF,GQ,Ensembl,ALL' (first fourth annotations are defined and ordered, and all other are included but not ordered)

=item B<--sort_by=<string>>

Sort variants by a field (default 'PZFlag,PZScore'). Only 2 levels

Example: 'PZFlag,PZScore' (to sort by relevance) or 'CHROM_NUM,POS' (to sort by position)

=item B<--order_by=<string>>

Order variants by a field (default '').

Example: 'DESC,DESC' (useful to sort by relevance)

=item B<--columns=<string>>

Additional columns with values. All additional columns will be added with the associated values. Available only for 'tab' format

Example: 'run:runA,sample:sample1'

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

## Format
my %format_allowed=("tab"=>"tab-delimited", "vcf"=>"VCF");
#my $format=lc($parameters{"format"});
#print %format_allowed;
my $format="";
foreach my $tag (keys %format_allowed) {
	if (lc $tag eq lc($parameters{"format"})) {
		$format=$tag;
	    #$format=$format_allowed{$tag};
	}
}
if ($format eq "") {
	$format="tab";
};#if
$header.="#[INFO] Output format: ".$format_allowed{$format}."\n";


## Output file
my $output_file;
my $output_filename_pattern="translated";
if (-e $parameters{"output"} && 0) {
	print "#[ERROR] output file '".$parameters{"output"}."' DOES exist\n";
	pod2usage();
	exit 1;
} else {
	$output_file=$parameters{"output"};
	if (trim($output_file) eq "") {
		$output_file=$input_file;
		$output_file =~ s/\.vcf$/\.$output_filename_pattern\.$format/g; #/\.vcf$/\.output\.vcf/;
	};#if
	$header.="#[INFO] Output VCF file: $output_file\n";
};#if

# columns
my @columns_split=split(/,/, $parameters{"columns"});
my %additional_columns;
foreach my $column_infos (@columns_split) {
	my @columns_split_split=split(/:/, $column_infos);
	print $columns_split_split[0]." = ".$columns_split_split[1]."\n" if $DEBUG;
	$additional_columns{trim($columns_split_split[0])}=trim($columns_split_split[1]);
	#my $sep=((defined $parameters{"fields"} || $parameters{"fields"} eq "")?",":"");
	#$parameters{"fields"}.=$sep.trim($columns_split_split[0]);
	#foreach my $column_infos (@columns_split_split) {
	#	@columns_split=split(/:/, $column_infos);
	#};#foreach
};#foreach
print Dumper(\%additional_columns) if $DEBUG;

## Annotation to Show
my %vcf_header_hash_test=read_vcf($parameters{"input"},"header");
@fields_split=split(/,/, $parameters{"fields"});
#print $parameters{"fields"};
print "@fields_split\n" if $DEBUG;
if ( "ALL" ~~ @fields_split && 1 ) {
	#print "ALLinARRAY\n";
	my $ALLfieldsList;
	while ( my ($ann, $poss) = each(%{$vcf_header_hash_test{"INFO"}}) ) {
		if  ( $ann ~~ @fields_split ) {
			
		} else {
			#if ($ALLfieldsList=="") { print "OK"; print $ALLfieldsList; print "\n"; } else { print "ko"; print $ALLfieldsList; print "\n"; }; #if
			my $sep=((defined $ALLfieldsList)?",":"");
			#print $chr;
			$ALLfieldsList.=$sep.$ann;
		};#if
	};
	#print $ALLfieldsList;
	@fields_split = map {s/ALL/$ALLfieldsList/g; $_; } @fields_split;
	my $fields_split_joined=join(",", @fields_split);
	$parameters{"fields"}=$fields_split_joined;
	@fields_split=split(/,/, $fields_split_joined);
};#if
my @fields_split = do { my %seen; grep { !$seen{$_}++ } @fields_split };
print "@fields_split\n" if $DEBUG;
my %fields_hash = map { $fields_split[$_] => $_ } 0..$#fields_split;
$header.="#[INFO] fields to show (".@fields_split."): @fields_split\n";
#print Dumper(\%fields_hash); #if $DEBUG;

#foreach my $annotation_name (sort { $fields_hash{$a} <=> $fields_hash{$b} } keys %fields_hash) {
#	print "$annotation_name\n";
#};#foreach

#$header.="## fields to show: ".@fields_split."\n";
#print $header;
#exit 0;

# Sort by
my @sort_by=split(/,/, $parameters{"sort_by"});
$header.="#[INFO] Sort by: ".((join("",@sort_by) ne "")?"@sort_by":"-")."\n";

# Order by
my @order_by=split(/,/, $parameters{"order_by"});
$header.="#[INFO] Order: ".((join("",@order_by) ne "")?"@order_by":"-")."\n";


#exit 0;

## DEBUG
##########

#print Dumper(\%config) if $DEBUG;
#print Dumper(\%config_annotation) if $DEBUG;
#print Dumper(\%config_filter) if $DEBUG;
#print Dumper(\%annotation_filters) if $DEBUG;


## MAIN
#########

# Variables
$output="";


my %variants;
my %annotations;

# VERBOSE
if ($VERBOSE) {
	print "#[INFO] Input file: $input_file\n";
	print "#[INFO] Output file: $output_file\n";
	print "#[INFO] Output format: $format\n";
	print "#[INFO] Annotation: $fields\n";
	print "#\n";
};#if

# READ VCF FILE
my %annotation_output=read_vcf($parameters{"input"});
my $vcf_header=read_vcf($parameters{"input"},"header_flat");
my $vcf_header_hash=read_vcf($parameters{"input"},"header");
#print $vcf_header_hash if $DEBUG;
#exit 0;

# DEBUG
if ($DEBUG) {
	#print Dumper($vcf_header_hash) if $DEBUG;
	#print Dumper($vcf_header_hash{"INFO"}) if $DEBUG;
	#while ( my ($chr, $poss) = each(%{$vcf_header_hash{"INFO"}}) ) {
	#	print $chr;
	#};
	#print Dumper(\%annotation_output) if $DEBUG;
	#print Dumper(\%annotations) if $DEBUG;
	#print Dumper(\%filter_list) if $DEBUG;
	#print Dumper(\%fields_hash) if $DEBUG;
	#exit 0;
};#DEBUG

my $PZHeader;
my $nb_variant=0;
my $output_header;
my $output_content;
my %output_content_hash;
my $output_all=0;
my @left_annotations=("CHROM","POS","ID","REF","ALT","QUAL","FILTER");
#my @right_annotations=("FORMAT","SAMPLES","SAMPLE_LIST","SAMPLE_LIST_CONCAT","SAMPLE_LIST_CONCAT2","SAMPLE_LIST_HEADER");
my @right_annotations_initial=("FORMAT","SAMPLE_LIST_CONCAT");
#my @right_annotations=("FORMAT","SAMPLES","SAMPLE_LIST","SAMPLE_LIST_CONCAT","SAMPLE_LIST_CONCAT2","SAMPLE_LIST_HEADER");
#my %excluded=map { $_ => 1 } ("CHROM","POS","ID","REF","ALT","QUAL","FILTER","FORMAT","SAMPLE_LIST");
my %excluded=map { $_ => 1 } @left_annotations, @right_annotations_initial, ("SAMPLE_LIST");
#my %excluded=map { $_ => 1 } @left_annotations, @right_annotations;
#print Dumper(\%excluded) if $DEBUG;

# Excluded
#$excluded{"CHROM","POS","ID","REF","ALT","QUAL","FILTER"}=1;

# HEADER
if ($format eq "vcf") {
	#print "########### VCF\n" if $DEBUG;
	$output_header=$vcf_header
};#if


# READ HASH Variants to construct the list of all annotations
while ( my ($chr, $poss) = each(%annotation_output) ) {
while ( my ($pos, $refs) = each(%{$poss}) ) {
while ( my ($ref, $alts) = each(%{$refs}) ) {
while ( my ($alt, $variant_values) = each(%{$alts}) ) {

	# Annotations fields
	my $variant_annotations=$$variant_values{"INFOS"};
	foreach my $variant_value (("CHROM","POS","ID","REF","ALT","QUAL","FILTER",)) {
		$$variant_annotations{$variant_value}=$$variant_values{$variant_value};
	};#foreach

	# ANNOTATIONS
	while ( my ($annotation_name, $annotation_value) = each(%{$variant_annotations}) ) {
		$annotations{$annotation_name}{$annotation_value}++;
	};#while
}}}}

#print Dumper(\%annotations) if $DEBUG;
#print Dumper(\%annotation_output) if $DEBUG;

# READ HASH Variants
while ( my ($chr, $poss) = each(%annotation_output) ) {
while ( my ($pos, $refs) = each(%{$poss}) ) {
while ( my ($ref, $alts) = each(%{$refs}) ) {
while ( my ($alt, $variant_values) = each(%{$alts}) ) {

	my $variant_annotations=$$variant_values{"INFOS"};
	while ( my ($a, $b) = each(%{$$variant_values{"INFOS"}}) ) {
		if (!defined $b) {
			#print "# UNDEFINED !!!!!!!!!!!!!!!\n" if $DEBUG;
			$b=$a;
		};#if
		$$variant_annotations{$a}=$b;
		#print "#$a=$b\n" if $DEBUG;
		#print "#$a\n" if $DEBUG;
	};#while
	foreach my $variant_value (("CHROM","POS","ID","REF","ALT","QUAL","FILTER","SAMPLE_LIST_CONCAT","SAMPLE_LIST_CONCAT2","SAMPLE_LIST","SAMPLES","FORMAT","SAMPLE_LIST_HEADER")) {
		$$variant_annotations{$variant_value}=$$variant_values{$variant_value};
	};#foreach
	#foreach my $sample ( split(/,/, $$variant_values{"SAMPLE_LIST_CONCAT2"})) {    
	my @samples_names=split(/\t/, $$variant_values{"SAMPLE_LIST_HEADER"});
	my @samples_values=split(/\t/, $$variant_values{"SAMPLE_LIST_CONCAT2"});
	print "###### samples_names @samples_names\n" if $DEBUG; #
	print "###### samples_values @samples_values\n" if $DEBUG; #
	my $i=0;
	foreach my $sample (@samples_names) {    
		#print "$i\n";
		$$variant_annotations{$sample}=$samples_values[$i];
		#@right_annotations=(@right_annotations, ($sample));
		#push(@right_annotations, ($sample));
		$i++;
	};#foreach
	if ($format eq "vcf") {
		@right_annotations=(("FORMAT"), @samples_names);
	} else {
		#push(@right_annotations, @samples_name);
		@right_annotations=(@right_annotations_initial, @samples_names);
	};#if
	# VARIANT
	$nb_variant++;
	#print "$chr:$pos|$ref>$alt (chr_num:".$$variant_annotations{"CHROM_NUM"}.")\n" if $DEBUG;
	
	# Line
	my $line;

	# sort by
	$sort_by_index1=(defined $sort_by[0] && $sort_by[0] ne "" && defined $$variant_annotations{$sort_by[0]})?$$variant_annotations{$sort_by[0]}:$nb_variant;
	$sort_by_index2=(defined $sort_by[1] && $sort_by[1] ne "" && defined $$variant_annotations{$sort_by[1]})?$$variant_annotations{$sort_by[1]}:$nb_variant;
	#print "$sort_by_index1 ($sort_by[0]),$sort_by_index2 ($sort_by[1])\n" if $DEBUG;

	# First annotations LEFT
	if ($format eq "vcf") {
		#print "########### VCF\n" if $DEBUG;
		foreach my $annotation_name ("CHROM","POS","ID","REF","ALT","QUAL","FILTER") {
			my $annotation_value=$$variant_annotations{$annotation_name};
			#print "annotation_value=$annotation_value\n" if $DEBUG;
			$line.=((!defined $line)?"":"\t").$annotation_value;
		};#if
	
	} else {
#print "\n\n\n";
		foreach my $annotation_name (@left_annotations) {
#print "$annotation_name\n";
			my $annotation_value=$$variant_annotations{$annotation_name};

			my $sep_header=((!defined $output_header)?"":"\t");
			my $sep_line=((!defined $line)?"":"\t");

			# HEADER
			if ($nb_variant == 1) {

				# TAB
				if ($format eq "tab") {
					$output_header.=$sep_header.$annotation_name;
				};#if

			};#if

			# LINE
			# TAB
			if ($format eq "tab") {
				$line.=$sep_line.$annotation_value;
			};#if

		};#foreach
	};#if
	
	# list of annotations to show
	
	my $line_vcf_info="";
	
#print "\n\n\n";
	foreach my $annotation_name (sort { $fields_hash{$a} <=> $fields_hash{$b} } keys %fields_hash) {
#print "$annotation_name\n";
		my $annotation_value=$$variant_annotations{$annotation_name};

		my $sep_header=((!defined $output_header)?"":"\t");
		my $sep_line=((!defined $line)?"":"\t");

		if (lc($annotation_name) eq "all") {
			$output_all=1;
		} else {

			# HEADER
			if ($nb_variant == 1) {

				# TAB
				if ($format eq "tab") {
					$output_header.=$sep_header.$annotation_name;
				};#if

			};#if

			# LINE
			# VCF
			if ($format eq "vcf") {
				if (defined $annotation_value && trim($annotation_value) ne "") {
					$line_vcf_info.=((!defined $line_vcf_info || $line_vcf_info eq "")?"":";").$annotation_name."=".$annotation_value;
				};#if
			};#if
			# TAB
			if ($format eq "tab") {
				$line.=$sep_line.$annotation_value;
			};#if

		};#if


	};#if

	#if ($format eq "vcf") {
	#	$line.="\t".$line_vcf_info;
	#};#if
	
	#my $line_vcf_info=0;
	# output all
	if ($output_all && 1) {
		while ( my ($annotation_name, $annotation_values) = each(%annotations) ) { # $annotations{$annotation_name}{$annotation_value}++;
			
			my $annotation_value=$$variant_annotations{$annotation_name};
			#print "    $annotation_name\t\t$annotation_value\n" if $DEBUG;
			
			if (!defined $fields_hash{$annotation_name} && !defined $excluded{$annotation_name}) {
				
				my $sep_header=((!defined $output_header)?"":"\t");
				my $sep_line=((!defined $line)?"":"\t");

				# HEADER
				if ($nb_variant == 1) {

					# TAB
					if ($format eq "tab") {
						$output_header.=$sep_header.$annotation_name;
					};#if

				};#if

				# LINE
				# VCF
				if ($format eq "vcf") {
					if (defined $annotation_value && trim($annotation_value) ne "") {
						$line_vcf_info.=((!defined $line_vcf_info || $line_vcf_info eq "")?"":";").$annotation_name."=".$annotation_value;
					};#if
				};#if
			# TAB
				if ($format eq "tab") {
					$line.=$sep_line.$annotation_value;
				};#if
			};#if
		};#while
	};#if
	
	#if ($format eq "vcf") {
	#	$line.="\t".$line_vcf_info;
	#};#if
	# Last annotations RIGHT

	# Additional columns
	while ( my ($annotation_name, $annotation_values) = each(%additional_columns) ) { # $annotations{$annotation_name}{$annotation_value}++;
		
		print "???    $column_name\t\t$column_values\n" if $DEBUG;
		
		my $annotation_value=$annotation_values;

		#print "$annotation_name (@right_annotations\n";
		#$$variant_annotations{$annotation_name};

		my $sep_header=((!defined $output_header)?"":"\t");
		my $sep_line=((!defined $line)?"":"\t");

		# HEADER
		if ($nb_variant == 1) {

			# TAB
			if ($format eq "tab") {
				$output_header.=$sep_header.$annotation_name;
			};#if

		};#if

		# LINE
		# VCF
		#if ($format eq "vcf") {
		#	$line.=$sep_line.$annotation_value;
		#};#if
		if ($format eq "vcf") {
			if (defined $annotation_value && trim($annotation_value) ne "") {
				$line_vcf_info.=((!defined $line_vcf_info || $line_vcf_info eq "")?"":";").$annotation_name."=".$annotation_value;
			};#if
		};#if
		# TAB
		if ($format eq "tab") {
			$line.=$sep_line.$annotation_value;
		};#if

	};#foreach
	
	if ($line_vcf_info eq "") {
		$line_vcf_info=".";
	};#if
	if ($format eq "vcf" && $line_vcf_info ne "") {
		$line.="\t".$line_vcf_info;
	};#if

	foreach my $annotation_name (@right_annotations) {
		my $annotation_value=$$variant_annotations{$annotation_name};

		#print "$annotation_name (@right_annotations\n";
		#$$variant_annotations{$annotation_name};

		my $sep_header=((!defined $output_header)?"":"\t");
		my $sep_line=((!defined $line)?"":"\t");

		# HEADER
		if ($nb_variant == 1) {

			# TAB
			if ($format eq "tab") {
				$output_header.=$sep_header.$annotation_name;
			};#if

		};#if

		# LINE
		# TAB
		if ($format eq "vcf") {
			$line.=$sep_line.$annotation_value;
		};#if
		# TAB
		if ($format eq "tab") {
			$line.=$sep_line.$annotation_value;
		};#if

	};#foreach
	
	
	# Output Content
	$output_content_hash{$sort_by_index1}{$sort_by_index2}.="$line\n";

	#print "line=$line\n" if $DEBUG;

}
}
}
}


#print Dumper(\%output_content_hash) if $DEBUG;

# create content
foreach my $index1 (($order_by[0] eq "ASC" || $order_by[0] eq "")?sort {$a <=> $b || $a cmp $b} keys %output_content_hash:reverse sort {$a <=> $b || $a cmp $b} keys %output_content_hash) {
	my $indexes=$output_content_hash{$index1};
	#print Dumper(\%indexes) if $DEBUG;
	foreach my $index2 (($order_by[1] eq "ASC" || $order_by[1] eq "")?sort {$a <=> $b || $a cmp $b} keys %{$indexes}:reverse sort {$a <=> $b || $a cmp $b} keys %{$indexes}) {
		#print "$index1 $index2\n" if $DEBUG;
		$output_content.=$output_content_hash{$index1}{$index2};
		#print "$output_content\n" if $DEBUG;
	};#foreach
};#foreach

#print "$output_header\n" if $DEBUG;
#print "$output_content\n" if $DEBUG;

if ($DEBUG) {
	#print Dumper(\%annotation_output) if $DEBUG;
	#print Dumper(\%annotations) if $DEBUG;
	#print Dumper(\%filter_list) if $DEBUG;
	#print Dumper(\%output_content_hash) if $DEBUG;
};#DEBUG


# create full file content
my $output_full="$output_header\n$output_content";
# remove blank lines
$output_full =~ s/\n\s*/\n/g;

#print "OUTPUT '$output_file': $output_header\n$output_content\n" if $DEBUG;
## Write OUTPUT
#open the files

open(FILE_OUTPUT, ">$output_file") || die "Problem to open the file '$output_file': $!";
print FILE_OUTPUT $output_full;
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

