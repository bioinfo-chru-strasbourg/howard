#!/usr/bin/perl
############################
# VCF calculation          #
# Author: Antony Le Béchec #
############################

## Main Information
#####################

our %information = ( #
	'script'	=>  	basename($0),		# Script
	'release'	=>  	"0.9.2.3b",			# Release
	'date'		=>  	"20181217",		# Release parameter
	'author'	=>  	"Antony Le Béchec",	# Author
	'copyright'	=>  	"HUS",			# Copyright
	'licence'	=>  	"GNU AGPL V3",		# Licence
);

## Release Notes
##################
# 20180514-0.9.2.1b: Add help information
# 20180823-0.9.2.2b: Add Other NOMEN to avoid misassignation of GNOMEN. Fix bug with INFO/. for empty INFO field
# 20181217-0.9.2.3b: Change Number/Type/Description of new INFO/FORMAT header generateed. Default VAF '.' if no info.  Bug fixed, Description, Number and Type consolidated within INFO


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
use List::Util qw(sum min max);
#use List::MoreUtils qw(uniq);

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


=item B<--calculation=<string>>

calculation to do.
Example: "VAF", "VAF,CNOMEN,PNOMEN", "VAF,VAF_STATS,NOMEN,VARTYPE"

Currently available :

VAF: add VAF calculation for each sample, depending on information provided by callers (by order: FREQ, DP4, AD)

VAF_STATS: add VAF statistics on INFO field, from VAF calculation for each sample, depending on information provided by callers (by order: VAF, FREQ, DP4, AD)

CALLING_QUALITY: Calling quality (FORMAT/*) of all samples in case of multiSample VCF, or all pipelines in case of multipipeline VCF

CALLING_QUALITY_EXPLODE: Explode all Calling quality (FORMAT/*) in multiple fields in INFOS

NOMEN: Find the NOMEN from HGVS annotation. Depend on transcript of reference. If no transcript of reference, first transcript. This option create annotations on INFO field: NOMEN (full HGVS annotation), CNOMEN (DNA level mutation "c."), PNOMEN (Protein level mutation "p."), TNOMEN (transcript), ENOMEN (exon, if any), GNOMEN (gene, if any)

BARCODE: Calculate VaRank BarCode

GENOTYPECONCORDANCE: If all samples with the same genotype

FINDBYPIPELINES: Number of pipeline calling the variant

VARTYPE: SNV if X>Y, MOSAIC if X>Y,Z or X,Y>Z, INDEL if XY>Z or X>YZ

=item B<--transcripts=<string>>

file containing default transcripts (with or without release) for each gene. If mmultiple transcripts for a gene, priority is assigned by position in the list (e.g. TRANSCRIPT2a has priority over TRANSCRIPT2b).

format :

TRANSCRIPT1	GENE1

TRANSCRIPT2a	GENE2

TRANSCRIPT2b	GENE2

TRANSCRIPT3	GENE3

...

=item B<--nomen_fields=<string>>

List of names of annotation field to use to calculate/extract NOMEN annotation.

Field order determine the priority

Format: field1,field2,...

Examples: "hgvs", "snpeff_hgvs", "snpeff_hgvs,hgvs"

Default: "hgvs"

=item B<--trio=<string>>

List of sample to identify a trio. Will automatically calculate VaRank barcode, and add INFO/trio_variant_type (either "denovo", "dominant", "recessive") as:

case "001": denovo

case "011","101","111","021","201","121","211": dominant

case "112","212","122","222": recessive

format: "father_sample_name,mother_sample_name,child_sample_name"

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
my $output_filename_pattern="calculated";
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

# TRIO
my $trio=$parameters{"trio"};
@trio_array=split(",",$trio);
if (defined $trio && $trio ne "") {
	$parameters{"calculation"}.=",BARCODE";
};#if
print "CALCULATION: '$trio' ".$parameters{"calculation"}."\n" if $DEBUG;


# CALCULATION
my $calculation=$parameters{"calculation"};
if (!defined $calculation) {
	$calculation=""; # if not defined by the list, default
};#if
@calculation_list_input=uniq(split(",",$calculation));

# TRANSCRIPTS
my $transcripts_input=$parameters{"transcripts"};
my @default_transcripts;
if (-e $transcripts_input) { #if file exist
	# Open the file
	open(FILE_TRANSCRIPTS, $transcripts_input) || die "Problem to open the file: $!";
	# Read the file
	while(<FILE_TRANSCRIPTS>) {
		# delete \n character
		chomp;
		# ignore blank lines
		next if /^\s*$/;
		# ignore commented lines
		next if /^\s*\;|#/;
		# read line
		$line=$_;
		# split line
		#my @line_split1=split(".",split("\t",$line));
		my @line_split1=split("\t",$line);
		my @line_split2=split(/\./,$line_split1[0]);
		push @default_transcripts, trim($line_split2[0]);
	};#while
};#if
print "TRANSCRIPTS @default_transcripts" if $DEBUG;
#exit 0;

#print "force=$force\n" if $DEBUG;
#exit 0;




$header.="#[INFO] Calculations: @calculation_list_input\n";

# HARD
my $hard_filter=$parameters{"hard"};


# Annotation type
my $annotation_type="calculation";
my $description_plus=" [Release=".$information{"release"}.";Date=".$information{"date"}.";AnnotationType=$annotation_type]";

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
	print "#[INFO] Input VCF file: $input_file\n";
	print "#[INFO] Output VCF file: $output_file\n";
	#print "# Output format: @output_formats\n";
	print "#[INFO] Calculation: @calculation_list_input\n";
	#print "# Filters: ".(map { $_ } keys %filter_list)."\n";# ;
	print "#\n";
};#if

# READ VCF FILE
my %annotation_output=read_vcf($parameters{"input"});

#my $vcf_header=read_vcf($parameters{"input"},"header_flat");
my %vcf_header_hash=read_vcf($parameters{"input"},"header");


my $PZHeader;

#my %vcf_header_hash_test=read_vcf($parameters{"input"},"header");

if ($DEBUG) {
	#print Dumper(\%annotation_output) if $DEBUG;
	#print Dumper(\%annotations) if $DEBUG;filter_list
	#print $vcf_header if $DEBUG;
	#print Dumper(\%vcf_header_hash) if $DEBUG;
	#print "TRUC\n";
};#DEBUG

my %annotation_output2=%annotation_output;
#print Dumper(\%annotation_output2) if $DEBUG;
# READ HASH Variants
while ( my ($chr, $poss) = each(%annotation_output2) ) {
while ( my ($pos, $refs) = each(%{$poss}) ) {
while ( my ($ref, $alts) = each(%{$refs}) ) {
while ( my ($alt, $variant_values) = each(%{$alts}) ) {

	print "$chr:$pos|$ref>$alt\n" if $DEBUG;

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



	#$$variant_values{"SAMPLE_LIST_CONCAT"}
	# VARIANT

	print "# Calculations: @calculation_list_input\n" if $DEBUG;



	foreach my $one_calculation (@calculation_list_input) {
		print "#\n# Calculation '$one_calculation'\n" if $DEBUG;
		my $is_calculated=0;



		########################
		# GENOTYPE_CONCORDANCE #
		########################

		if (uc(trim($one_calculation)) eq "GENOTYPECONCORDANCE") {
			print "# Calculation '$one_calculation'...\n" if $DEBUG;

			if ($force || !defined $$variant_values{"INFOS"}{"GenotypeConcordance"}) {

				my $GENOTYPE_CONCORDANCE="";
				my $FindByPipelines="";

				my @GT_LIST;
				my $number_pipelines=0;
				my $seen=0;
				while ( my ($sample_name, $sample_quality) = each(%$samples) ) {
					#print "$sample_name\n" if $DEBUG;
					$number_pipelines++;
					if (defined $$sample_quality{"GT"}) {
						my $GT=$$sample_quality{"GT"};
						print "\tGenotype Concordance with GT (".$$sample_quality{"GT"}.") ($GT)\n" if $DEBUG ;
						if ($GT ne "./.") {
							push @GT_LIST, $GT;
						};#if
						#push @GT_LIST, $GT;
					};#if
				};# while
				print "\tGT_LIST @GT_LIST\n" if $DEBUG ;
				my @GT_LIST_UNIQ = do { my %seen; grep { !$seen{$_}++ } @GT_LIST };
				$seen=@GT_LIST;
				#my @GT_LIST_UNIQ= uniq @GT_LIST;
				$GENOTYPE_CONCORDANCE=((@GT_LIST_UNIQ >1)?"FALSE":"TRUE");
				#$FindByPipelines="$seen/$number_pipelines";
				print "\tGT_LIST_UNIQ @GT_LIST_UNIQ\n" if $DEBUG ;
				print "\tnumber_pipelines $number_pipelines\n" if $DEBUG ;
				print "\tseen $seen\n" if $DEBUG ;
				print "\tGENOTYPE_CONCORDANCE $GENOTYPE_CONCORDANCE\n" if $DEBUG ;
				#print "\tFindByPipelines $FindByPipelines\n" if $DEBUG ;


				if ($GENOTYPE_CONCORDANCE ne "") {
					# GenotypeConcordance in INFO
					$annotation_output{$chr}{$pos}{$ref}{$alt}{"INFOS"}{"GenotypeConcordance"}=$GENOTYPE_CONCORDANCE;
					# ADD to HEADER
					$vcf_header{"INFO"}{"GenotypeConcordance"}{"Number"}="1";
					$vcf_header{"INFO"}{"GenotypeConcordance"}{"Type"}="String";
					$vcf_header{"INFO"}{"GenotypeConcordance"}{"Description"}="\"TRUE if genotypes are concordant between all the callers, else FALSE$description_plus\"";
				};#if

			};#if

			#};#while

			$is_calculated=1;
		};#if


		###################
		# FINDBYPIPELINES #
		###################

		if (uc(trim($one_calculation)) eq "FINDBYPIPELINES") {
			print "# Calculation '$one_calculation'...\n" if $DEBUG;

			if ($force || !defined $$variant_values{"INFOS"}{"FindByPipelines"}) {

				my $GENOTYPE_CONCORDANCE="";
				my $FindByPipelines="";

				my @GT_LIST;
				my $number_pipelines=0;
				my $seen=0;
				while ( my ($sample_name, $sample_quality) = each(%$samples) ) {
					#print "$sample_name\n" if $DEBUG;
					$number_pipelines++;
					if (defined $$sample_quality{"GT"}) {
						my $GT=$$sample_quality{"GT"};
						print "\tGenotype Concordance with GT (".$$sample_quality{"GT"}.") ($GT)\n" if $DEBUG ;
						if ($GT ne "./.") {
							push @GT_LIST, $GT;
						};#if
					};#if
				};# while
				print "\tGT_LIST @GT_LIST\n" if $DEBUG ;
				my @GT_LIST_UNIQ = do { my %seen; grep { !$seen{$_}++ } @GT_LIST };
				$seen=@GT_LIST;
				#my @GT_LIST_UNIQ= uniq @GT_LIST;
				#$GENOTYPE_CONCORDANCE=((@GT_LIST_UNIQ >1)?"FALSE":"TRUE");
				$FindByPipelines="$seen/$number_pipelines";
				print "\tGT_LIST_UNIQ @GT_LIST_UNIQ\n" if $DEBUG ;
				print "\tnumber_pipelines $number_pipelines\n" if $DEBUG ;
				print "\tseen $seen\n" if $DEBUG ;
				#print "\tGENOTYPE_CONCORDANCE $GENOTYPE_CONCORDANCE\n" if $DEBUG ;
				print "\tFindByPipelines $FindByPipelines\n" if $DEBUG ;


				if ($FindByPipelines ne "") {
					# FindByPipelines in INFO
					$annotation_output{$chr}{$pos}{$ref}{$alt}{"INFOS"}{"FindByPipelines"}=$FindByPipelines;
					# ADD to HEADER
					$vcf_header{"INFO"}{"FindByPipelines"}{"Number"}="1";
					$vcf_header{"INFO"}{"FindByPipelines"}{"Type"}="String";
					$vcf_header{"INFO"}{"FindByPipelines"}{"Description"}="\"Number of pipelines which find the variant (if value = 0, that to say that the variant was filtered in by all pipelines)$description_plus\"";
				};#if

			};#if


			#};#while

			$is_calculated=1;
		};#if


		#######
		# VAF #
		#######

		if (uc(trim($one_calculation)) eq "VAF" || uc(trim($one_calculation)) eq "VAF_STATS") {
		#if (uc(trim($one_calculation)) eq "VAF") {
			print "# Calculation '$one_calculation'...\n" if $DEBUG;

			my @VAF_LIST; #=(0.99,0.999);
			my $VAF_min="";
			my $VAF_max="";
			my $VAF_average="";
			my $VAF_median="";

			while ( my ($sample_name, $sample_quality) = each(%$samples) ) {
				#print "$sample_name ".$$sample_quality{"VAF"}." \n" if $DEBUG;

				#if ($force || !defined $$sample_quality{"VAF"}) {
				if ($force || !defined $annotation_output{$chr}{$pos}{$ref}{$alt}{"SAMPLES"}{$sample_name}{"VAF"}) {
					my $VAF="";


					# if FREQ
					#$$sample_quality{"FREQ"}="2.36%,5.25%";
					if ($VAF eq "" && defined $$sample_quality{"FREQ"}) {
					#if (defined $$sample_quality{"FREQ"}) {
						#print "\tVAF calculated on FREQ (".$$sample_quality{"FREQ"}.")\n" if $DEBUG ;
						my $FREQ=$$sample_quality{"FREQ"};
						$FREQ =~ s/%//g;
						print "\tVAF calculated on FREQ (".$$sample_quality{"FREQ"}." > ".$FREQ.")\n" if $DEBUG ;
						my @FREQ_split=split(/,/,$FREQ);
						my $FREQ_ALT = eval join '+', @FREQ_split;
						if ($FREQ_ALT) {
							$VAF=sprintf("%.4f", $FREQ_ALT)/100;
							print "\tVAF = $VAF\n" if $DEBUG ;
						};#if
					};#if

					# if DP4
					#$$sample_quality{"DP4"}="36,33,30,5";
					#$$sample_quality{"DP4"}="36,33,30,5,10,5";
					if ($VAF eq "" && defined $$sample_quality{"DP4"}) {
						print "\tVAF calculated on DP4 (".$$sample_quality{"DP4"}.")\n" if $DEBUG ;
						my @DP4_split=split(/,/,$$sample_quality{"DP4"});
						#my $DP4_DP = eval join '+', @DP4_split[0 .. 1];
						my $DP4_DP = eval join '+', @DP4_split[0 .. $#DP4_split];
						my $DP4_ALT = eval join '+', @DP4_split[2 .. $#DP4_split];
						if ($DP4_DP) {
							$VAF=sprintf("%.6f", $DP4_ALT/$DP4_DP);
							print "\tVAF = $VAF\n" if $DEBUG ;
						};#if
					};#if

					# if AD
					if ($VAF eq "" && defined $$sample_quality{"AD"}) {
						print "\tVAF calculated on AD (".$$sample_quality{"AD"}.")\n" if $DEBUG ;
						my @AD_split=split(/,/,$$sample_quality{"AD"});
						my $AD_DP = eval join '+', @AD_split;
						my $AD_ALT = eval join '+', @AD_split[1 .. $#AD_split];
						if ($AD_DP) {
							$VAF=sprintf("%.6f", $AD_ALT/$AD_DP);
							print "\tVAF = $VAF\n" if $DEBUG ;
						};#if
					};#if

					# Check VAF NOT null
					if ($VAF ne "") {
						push @VAF_LIST, $VAF ;
					} else {
					#if ($VAF eq "") {
						$VAF=".";
					};# if

					if (uc(trim($one_calculation)) eq "VAF") {
						if ($VAF ne "") {
							# VAF in FORMAT
							$annotation_output{$chr}{$pos}{$ref}{$alt}{"SAMPLES"}{$sample_name}{"VAF"}=$VAF;
							# ADD to FORMAT
							if (!grep( /^VAF$/, split(/:/,$annotation_output{$chr}{$pos}{$ref}{$alt}{"FORMAT"}))) {
								$annotation_output{$chr}{$pos}{$ref}{$alt}{"FORMAT"}.=":VAF";
							};#if
							# ADD to HEADER
							$vcf_header{"FORMAT"}{"VAF"}{"Number"}="1";
							$vcf_header{"FORMAT"}{"VAF"}{"Type"}="Float";
							$vcf_header{"FORMAT"}{"VAF"}{"Description"}="\"VAF Variant Frequency, calculated from quality$description_plus\"";

						};#if
					};#if

				} else {
					push @VAF_LIST, $$sample_quality{"VAF"} ;
					#@VAF_LIST=split(/,/,$$sample_quality{"VAF"});


				};#if

			};# while

			print "@VAF_LIST\n" if $DEBUG;

			#print "VAF ".$$variant_values{"VAF"}."\n" if $DEBUG;
			if (uc(trim($one_calculation)) eq "VAF_STATS") {

				if ($force || !defined $$variant_values{"INFOS"}{"VAF"}) {

				#if ($force || !defined $$variant_values{"VAF"}) {

				$VAF_min = min @VAF_LIST;
				$VAF_max = max @VAF_LIST;
				$VAF_average = mean(@VAF_LIST);
				$VAF_median = median(@VAF_LIST);
				#print sum(@AD_split)/@AD_split if $DEBUG ;
				print "\tVAF MIN = $VAF_min\n" if $DEBUG ;
				print "\tVAF MAX = $VAF_max\n" if $DEBUG ;
				print "\tVAF AVG = $VAF_average\n" if $DEBUG ;
				print "\tVAF MED = $VAF_median\n" if $DEBUG ;

				#if ($VAF ne "" && $VAF ne ".") {
				if (join("",@VAF_LIST) ne "") {
					my $VAF=join(",",@VAF_LIST);
					# VAF in INFO
					$annotation_output{$chr}{$pos}{$ref}{$alt}{"INFOS"}{"VAF"}=$VAF;
					# ADD to HEADER
					#$vcf_header{"FORMAT"}{"VAF"}{"Number"}=".";
					#$vcf_header{"FORMAT"}{"VAF"}{"Type"}="String";
					#$vcf_header{"FORMAT"}{"VAF"}{"Description"}="\"VAF Variant Frequency, calculated from quality\"";
					$vcf_header{"INFO"}{"VAF"}{"Number"}=".";
					$vcf_header{"INFO"}{"VAF"}{"Type"}="Float";
					$vcf_header{"INFO"}{"VAF"}{"Description"}="\"VAF Variant Frequency, calculated from quality$description_plus\"";
				};#if

				if ($VAF_min ne "") {
					# VAF_min in INFO
					$annotation_output{$chr}{$pos}{$ref}{$alt}{"INFOS"}{"VAF_min"}=$VAF_min;
					# ADD to HEADER
					$vcf_header{"INFO"}{"VAF_min"}{"Number"}="1";
					$vcf_header{"INFO"}{"VAF_min"}{"Type"}="Float";
					$vcf_header{"INFO"}{"VAF_min"}{"Description"}="\"VAF Variant Frequency minimum$description_plus\"";
				};#if

				if ($VAF_max ne "") {
					# VAF_min in INFO
					$annotation_output{$chr}{$pos}{$ref}{$alt}{"INFOS"}{"VAF_max"}=$VAF_max;
					# ADD to HEADER
					$vcf_header{"INFO"}{"VAF_max"}{"Number"}="1";
					$vcf_header{"INFO"}{"VAF_max"}{"Type"}="Float";
					$vcf_header{"INFO"}{"VAF_max"}{"Description"}="\"VAF Variant Frequency max$description_plus\"";
				};#if

				if ($VAF_average ne "") {
					# VAF_min in INFO
					$annotation_output{$chr}{$pos}{$ref}{$alt}{"INFOS"}{"VAF_average"}=$VAF_average;
					# ADD to HEADER
					$vcf_header{"INFO"}{"VAF_average"}{"Number"}="1";
					$vcf_header{"INFO"}{"VAF_average"}{"Type"}="Float";
					$vcf_header{"INFO"}{"VAF_average"}{"Description"}="\"VAF Variant Frequency average$description_plus\"";
				};#if

				if ($VAF_median ne "") {
					# VAF_min in INFO
					$annotation_output{$chr}{$pos}{$ref}{$alt}{"INFOS"}{"VAF_median"}=$VAF_median;
					# ADD to HEADER
					$vcf_header{"INFO"}{"VAF_median"}{"Number"}="1";
					$vcf_header{"INFO"}{"VAF_median"}{"Type"}="Float";
					$vcf_header{"INFO"}{"VAF_median"}{"Description"}="\"VAF Variant Frequency median$description_plus\"";
				};#if

				};#if

			};#if

			$is_calculated=1;
		};#if


		################
		# CALLING_QUAL #
		################

		if (uc(trim($one_calculation)) eq "CALLING_QUALITY") {
			print "# Calculation '$one_calculation'...\n" if $DEBUG;

			# Catch CALLING_QUALITY
			my $CALLING_QUALITY="SAMPLE:".$$variant_values{"FORMAT"}."|".$$variant_values{"SAMPLE_LIST_CONCAT"};
			# cleaning CALLING_QUALITY
			$CALLING_QUALITY =~ tr/;/,/;

			print "CALLING_QUALITY: ".$CALLING_QUALITY."\n" if $DEBUG;

			if ($CALLING_QUALITY ne "") {
				$annotation_output{$chr}{$pos}{$ref}{$alt}{"INFOS"}{"CALLING_QUALITY"}=$CALLING_QUALITY;
				$vcf_header{"INFO"}{"CALLING_QUALITY"}{"Number"}="1";
				$vcf_header{"INFO"}{"CALLING_QUALITY"}{"Type"}="String";
				$vcf_header{"INFO"}{"CALLING_QUALITY"}{"Description"}="\"Calling quality (FORMAT/*) of all samples in case of multiSample VCF, or all pipelines in case of multipipeline VCF$description_plus\"";

			};#if

			$is_calculated=1;
		};#if


		########################
		# CALLING_QUAL_EXPLODE #
		########################

		if (uc(trim($one_calculation)) eq "CALLING_QUALITY_EXPLODE") {
			print "# Calculation '$one_calculation'...\n" if $DEBUG;

			#print Dumper($$variant_values{"SAMPLES"}) if $DEBUG;


			while ( my ($sample_name, $sample_quality) = each(%$samples) ) {
				while ( my ($quality_name, $quality_value) = each(%$sample_quality) ) {
					#print "$sample_name ".$$sample_quality{"VAF"}." \n" if $DEBUG;
					my $CALLING_QUALITY_EXPLODE_name="CALLING_QUALITY_EXPLODE_".$quality_name."_".$sample_name;
					my $CALLING_QUALITY_EXPLODE_value=$quality_value;
					# cleaning CALLING_QUALITY_EXPLODE_value
					$CALLING_QUALITY_EXPLODE_value =~ tr/;/,/;

					if ($CALLING_QUALITY_EXPLODE_value ne "") {
						my $value=$CALLING_QUALITY_EXPLODE_value;
						my $number=($vcf_header_hash{"FORMAT"}{$quality_name}{"Number"} eq "R")?".":$vcf_header_hash{"FORMAT"}{$quality_name}{"Number"};
						my $type=$vcf_header_hash{"FORMAT"}{$quality_name}{"Type"};
						my $description=$vcf_header_hash{"FORMAT"}{$quality_name}{"Description"};

						# FORMAT header Consolidation
						# Number
						if (trim($number) eq "") {
							$number=($vcf_header{"FORMAT"}{$quality_name}{"Number"} eq "R")?".":$vcf_header{"FORMAT"}{$quality_name}{"Number"};
						};#if
						if (trim($number) eq "") { $number="."; };#if
						# Type
						if (trim($type) eq "") {
							$type=($vcf_header{"FORMAT"}{$quality_name}{"Type"} eq "R")?".":$vcf_header{"FORMAT"}{$quality_name}{"Type"};
						};#if
						if (trim($type) eq "") { $type="String"; };#if
						# Description
						if (trim($description) eq "") {
							$description=($vcf_header{"FORMAT"}{$quality_name}{"Description"} eq "R")?".":$vcf_header{"FORMAT"}{$quality_name}{"Description"};
						};#if
						if (trim($description) eq "") { $description="\"unknown\""; };#if

						print "$CALLING_QUALITY_EXPLODE_name number=$number\n" if $DEBUG;
						print "$CALLING_QUALITY_EXPLODE_name type=$type\n" if $DEBUG;
						print "$CALLING_QUALITY_EXPLODE_name description=$description\n" if $DEBUG;

						$annotation_output{$chr}{$pos}{$ref}{$alt}{"INFOS"}{$CALLING_QUALITY_EXPLODE_name}=$value;
						$vcf_header{"INFO"}{$CALLING_QUALITY_EXPLODE_name}{"Number"}=$number;
						$vcf_header{"INFO"}{$CALLING_QUALITY_EXPLODE_name}{"Type"}=$type;
						$vcf_header{"INFO"}{$CALLING_QUALITY_EXPLODE_name}{"Description"}=$description;
					};#if
					#print "$CALLING_QUALITY_EXPLODE_name=$CALLING_QUALITY_EXPLODE_value \n" if $DEBUG;
					#print Dumper(\$vcf_header_hash{"FORMAT"}{$quality_name}) if $DEBUG;
					#print $vcf_header_hash{"FORMAT"}{$quality_name}{"Number"} if $DEBUG;
					#print "TEST".$vcf_header{"FORMAT"}{"VAF"}{"Number"}."\n" if $DEBUG;


				};
			};

			if (0) {
			my $CALLING_QUALITY="SAMPLE:".$$variant_values{"FORMAT"}."|".$$variant_values{"SAMPLE_LIST_CONCAT"};
			print "CALLING_QUALITY: ".$CALLING_QUALITY."\n" if $DEBUG;

			if ($CALLING_QUALITY ne "") {
				$annotation_output{$chr}{$pos}{$ref}{$alt}{"INFOS"}{"CALLING_QUALITY"}=$CALLING_QUALITY;
				$vcf_header{"INFO"}{"CALLING_QUALITY"}{"Number"}=".";
				$vcf_header{"INFO"}{"CALLING_QUALITY"}{"Type"}="String";
				$vcf_header{"INFO"}{"CALLING_QUALITY"}{"Description"}="\"Calling quality (FORMAT/*) of all samples in case of multiSample VCF, or all pipelines in case of multipipeline VCF$description_plus\"";

			};#if
			};#if


			$is_calculated=1;
		};#if


		###########
		# VARTYPE #
		###########

		if (uc(trim($one_calculation)) eq "VARTYPE") {
			print "# Calculation '$one_calculation'...\n" if $DEBUG;

			# Catch VARTYPE
			my $REF=$$variant_values{"REF"};
			my $ALT=$$variant_values{"ALT"};
			my $VARTYPE="UNDEFINED";
			if (length($REF) eq 1 && length($ALT) eq 1) {
				$VARTYPE="SNV";
			} elsif ($REF=~/,/ || $ALT=~/,/) {
				$VARTYPE="MOSAIC";
			} elsif (length($REF) ne length($ALT)) {
				$VARTYPE="INDEL";
			} else {
				$VARTYPE="UNDEFINED";
			};#if
			# cleaning VARTYPE
			#$VARTYPE =~ tr/;/,/;

			#print "VARTYPE: ".$VARTYPE."\n" if $DEBUG;

			if ($VARTYPE ne "") {
				$annotation_output{$chr}{$pos}{$ref}{$alt}{"INFOS"}{"VARTYPE"}=$VARTYPE;
				$vcf_header{"INFO"}{"VARTYPE"}{"Number"}="1";
				$vcf_header{"INFO"}{"VARTYPE"}{"Type"}="String";
				$vcf_header{"INFO"}{"VARTYPE"}{"Description"}="\"Variant type: SNV if X>Y, MOSAIC if X>Y,Z or X,Y>Z, INDEL if XY>Z or X>YZ$description_plus\"";

			};#if

			$is_calculated=1;
		};#if


		#########
		# Nomen #
		#########

		#my @default_transcripts=("NM_001122607");


		#if (uc(trim($one_calculation)) eq "NOMEN" || uc(trim($one_calculation)) eq "CNOMEN" || uc(trim($one_calculation)) eq "PNOMEN") {
		if (uc(trim($one_calculation)) eq "NOMEN") {
			print "# Calculation '$one_calculation'...\n" if $DEBUG;
			#print "N   ".$$variant_values{"INFOS"}{"NOMEN"}."\n" if $DEBUG;
			#while ( my ($sample_name, $sample_quality) = each(%$samples) ) {
			if ($force || !defined $$variant_values{"INFOS"}{"NOMEN"}) {
				#print "$sample_name\n" if $DEBUG;

				#print "".$$variant_values{"INFOS"}{"hgvs"} if $DEBUG;
				#print "".$annotation_output{$chr}{$pos}{$ref}{$alt}{"INFOS"}{"hgvs"} if $DEBUG;
				my $NOMEN=$$variant_values{"INFOS"}{"cnomen"};
				my $CNOMEN=$$variant_values{"INFOS"}{"cnomen"};
				my $PNOMEN=$$variant_values{"INFOS"}{"pnomen"};
				my $TVNOMEN=$$variant_values{"INFOS"}{"tvnomen"};
				my $TNOMEN=$$variant_values{"INFOS"}{"tnomen"};
				my $VNOMEN=$$variant_values{"INFOS"}{"vnomen"};
				my $ENOMEN=$$variant_values{"INFOS"}{"enomen"};
				my $GNOMEN=$$variant_values{"INFOS"}{"gnomen"};

				print "# Calculation '$one_calculation'...\n" if $DEBUG;

				print "# nomen_fields '".$parameters{'nomen_fields'}."'...\n" if $VERBOSE;
				my @nomen_fields_split=split(/,/,$parameters{'nomen_fields'});
				my $value_hgvs;
				foreach my $one_nomen_field (@nomen_fields_split) {
					print "# one_nomen_field '".$one_nomen_field."'...\n" if $VERBOSE;
					print "# one_nomen_field value '".$$variant_values{"INFOS"}{$one_nomen_field}."'...\n" if $VERBOSE;
					if (defined $$variant_values{"INFOS"}{$one_nomen_field} && $$variant_values{"INFOS"}{$one_nomen_field} ne "" && !defined $value_hgvs) {
						$value_hgvs=$$variant_values{"INFOS"}{$one_nomen_field};
					}; #if
				}; # foreach
				#my $value_hgvs=$$variant_values{$parameters{'nomen_fields'}};
				print "# nomen_fields value '".$value_hgvs."'...\n" if $VERBOSE;


				#if (!defined $CNOMEN || !defined $PNOMEN) {
				#if (defined $$variant_values{"INFOS"}{"hgvs"}) {
				if (defined $value_hgvs) {
					print "# nomen_fields value defined ! '".$value_hgvs."'...\n" if $VERBOSE;

					#my @HGVS_split=split(/,/,$$variant_values{"INFOS"}{"hgvs"});
					#my @HGVS_split=$value_hgvs;
					my @HGVS_split=split(/,/,$value_hgvs);
					my $assigned=0;
					foreach my $one_hgvs (@HGVS_split) {
						print "$one_hgvs\n" if $DEBUG;
						@one_hgvs_split=split(/:/,$one_hgvs);
						my $N;
						my $C;
						my $P;
						my $TV;
						my $T;
						my $V;
						my $E;
						my $G;
						my $O;
						my $T_D;
						foreach my $one_hgvs_infos (@one_hgvs_split) {
							#print "   $one_hgvs_infos\n" if $DEBUG;
							my $assigned=0;
							# CNOMEN
							if ($one_hgvs_infos=~ /^c\.(.*)$/) {
								$C=$one_hgvs_infos;
								$assigned=1;
							};#if
							# PNOMEN
							if ($one_hgvs_infos=~ /^p\.(.*)$/) {
								$P=$one_hgvs_infos;
								$assigned=1;
							};#if
							# TVNOMEN
							if ($one_hgvs_infos=~ /^NM_(.*)$/) {
								$TV=$one_hgvs_infos;
								$assigned=1;
							};#if
							# TNOMEN
							if ($one_hgvs_infos=~ /^NM_(.*)$/) {
								#$T=$one_hgvs_infos;
								my @one_hgvs_infos_T_split=split(/\./,$one_hgvs_infos);
								if (defined $one_hgvs_infos_T_split[0]) {
									$T=$one_hgvs_infos_T_split[0];
									$assigned=1;
								}; # if
							};#if
							# VNOMEN
							if ($one_hgvs_infos=~ /^NM_(.*)$/) {
								#$V=$one_hgvs_infos;
								my @one_hgvs_infos_T_split=split(/\./,$one_hgvs_infos);
								if (defined $one_hgvs_infos_T_split[1]) {
									$V=$one_hgvs_infos_T_split[1];
									$assigned=1;
								}; # if
							};#if
							# ENOMEN
							if ($one_hgvs_infos=~ /^exon(.*)$/) {
								$E=$one_hgvs_infos;
								$assigned=1;
							};#if
							# NOMEN other genomic
							if ($one_hgvs_infos=~ /^g\.(.*)$/) {
								$O=$one_hgvs_infos;
								$assigned=1;
							};#if
							# NOMEN other genomic
							if ($one_hgvs_infos=~ /^m\.(.*)$/) {
								$O=$one_hgvs_infos;
								$assigned=1;
							};#if
							# NOMEN other genomic
							if ($one_hgvs_infos=~ /^n\.(.*)$/) {
								$O=$one_hgvs_infos;
								$assigned=1;
							};#if
							# NOMEN other genomic
							if ($one_hgvs_infos=~ /^r\.(.*)$/) {
								$O=$one_hgvs_infos;
								$assigned=1;
							};#if
							# GNOMEN
							if (!$assigned) {
								$G=$one_hgvs_infos;
								$assigned=1;
							};#if
							# DEFAULT TRANSCRIPT
							#if (in_array($one_hgvs_infos,@default_transcripts)) {
							#my @one_hgvs_infos_split=split(/\./,$one_hgvs_infos);
							#my $transcript_NM=$one_hgvs_infos_split[0];
							my $transcript_NM=$T;
							#print "default transcipt to find: $transcript_NM\n" if $DEBUG;
							if (grep( /^$transcript_NM$/, @default_transcripts)) {
								print "default transcipt found: $transcript_NM\n" if $DEBUG;
								$T_D=$one_hgvs_infos;
							};#fi

						};#foreach
						$N="$C$P$TV$E$G";
						if ($N ne "") {
							if ( !$assigned || defined $T_D) {
								print "assi\n" if $DEBUG;
								$NOMEN=$one_hgvs;
								$CNOMEN=$C;
								$PNOMEN=$P;
								$TVNOMEN=$TV;
								$TNOMEN=$T;
								$VNOMEN=$V;
								$ENOMEN=$E;
								$GNOMEN=$G;
								$assigned=1;
							};#if
						};#if
						#if (!defined $PNOMEN || defined $T_D) {
						#	$PNOMEN=$P;
						#};#if
					};#foreach
					#print $HGVS_split[0]."\n" if $DEBUG;

				};#fi
				#};#if

				if (1) {
					print "NOMEN=$NOMEN\n" if $DEBUG;
					print "CNOMEN=$CNOMEN\n" if $DEBUG;
					print "PNOMEN=$PNOMEN\n" if $DEBUG;
					print "TVNOMEN=$TVNOMEN\n" if $DEBUG;
					print "TNOMEN=$TNOMEN\n" if $DEBUG;
					print "VNOMEN=$VNOMEN\n" if $DEBUG;
					print "ENOMEN=$ENOMEN\n" if $DEBUG;
					print "GNOMEN=$GNOMEN\n" if $DEBUG;
				};#if

				# NOMEN
				if ($NOMEN ne "") {
					# VAF in INFO
					$annotation_output{$chr}{$pos}{$ref}{$alt}{"INFOS"}{"NOMEN"}=$NOMEN;
					# ADD to HEADER
					$vcf_header{"INFO"}{"NOMEN"}{"Number"}="1";
					$vcf_header{"INFO"}{"NOMEN"}{"Type"}="String";
					$vcf_header{"INFO"}{"NOMEN"}{"Description"}="\"NOMEN hgvs nomenclature considered as reference hgvs (official transcript, first otherwise)$description_plus\"";
				};#if
				# CNOMEN
				if ($CNOMEN ne "") {
					# VAF in INFO
					$annotation_output{$chr}{$pos}{$ref}{$alt}{"INFOS"}{"CNOMEN"}=$CNOMEN;
					# ADD to HEADER
					$vcf_header{"INFO"}{"CNOMEN"}{"Number"}="1";
					$vcf_header{"INFO"}{"CNOMEN"}{"Type"}="String";
					$vcf_header{"INFO"}{"CNOMEN"}{"Description"}="\"CNOMEN hgvs nomenclature at cDNA level related to a transcript (TNOMEN)$description_plus\"";
				};#if
				# PNOMEN
				if ($PNOMEN ne "") {
					# VAF in INFO
					$annotation_output{$chr}{$pos}{$ref}{$alt}{"INFOS"}{"PNOMEN"}=$PNOMEN;
					# ADD to HEADER
					$vcf_header{"INFO"}{"PNOMEN"}{"Number"}="1";
					$vcf_header{"INFO"}{"PNOMEN"}{"Type"}="String";
					$vcf_header{"INFO"}{"PNOMEN"}{"Description"}="\"PNOMEN hgvs nomenclature at Protein level related to a transcript (TNOMEN)$description_plus\"";
				};#if
				# TVNOMEN
				if ($TVNOMEN ne "") {
					# VAF in INFO
					$annotation_output{$chr}{$pos}{$ref}{$alt}{"INFOS"}{"TVNOMEN"}=$TVNOMEN;
					# ADD to HEADER
					$vcf_header{"INFO"}{"TVNOMEN"}{"Number"}="1";
					$vcf_header{"INFO"}{"TVNOMEN"}{"Type"}="String";
					$vcf_header{"INFO"}{"TVNOMEN"}{"Description"}="\"TVNOMEN hgvs transcript with version (if any) used (e.g. for CNOMEN and PNOMEN)$description_plus\"";
				};#if
				# TNOMEN
				if ($TNOMEN ne "") {
					# VAF in INFO
					$annotation_output{$chr}{$pos}{$ref}{$alt}{"INFOS"}{"TNOMEN"}=$TNOMEN;
					# ADD to HEADER
					$vcf_header{"INFO"}{"TNOMEN"}{"Number"}="1";
					$vcf_header{"INFO"}{"TNOMEN"}{"Type"}="String";
					$vcf_header{"INFO"}{"TNOMEN"}{"Description"}="\"TNOMEN hgvs transcript used (e.g. for CNOMEN and PNOMEN)$description_plus\"";
				};#if
				# VNOMEN
				if ($VNOMEN ne "") {
					# VAF in INFO
					$annotation_output{$chr}{$pos}{$ref}{$alt}{"INFOS"}{"VNOMEN"}=$VNOMEN;
					# ADD to HEADER
					$vcf_header{"INFO"}{"VNOMEN"}{"Number"}="1";
					$vcf_header{"INFO"}{"VNOMEN"}{"Type"}="String";
					$vcf_header{"INFO"}{"VNOMEN"}{"Description"}="\"VNOMEN hgvs transcript version used (e.g. for CNOMEN and PNOMEN)$description_plus\"";
				};#if
				# ENOMEN
				if ($ENOMEN ne "") {
					# VAF in INFO
					$annotation_output{$chr}{$pos}{$ref}{$alt}{"INFOS"}{"ENOMEN"}=$ENOMEN;
					# ADD to HEADER
					$vcf_header{"INFO"}{"ENOMEN"}{"Number"}="1";
					$vcf_header{"INFO"}{"ENOMEN"}{"Type"}="String";
					$vcf_header{"INFO"}{"ENOMEN"}{"Description"}="\"ENOMEN hgvs exon nomenclature related to a transcript (TNOMEN)$description_plus\"";
				};#if
				# GNOMEN
				#if ($GNOMEN ne "" && 0) {
				if ($GNOMEN ne "") {
					# VAF in INFO
					$annotation_output{$chr}{$pos}{$ref}{$alt}{"INFOS"}{"GNOMEN"}=$GNOMEN;
					# ADD to HEADER
					$vcf_header{"INFO"}{"GNOMEN"}{"Number"}="1";
					$vcf_header{"INFO"}{"GNOMEN"}{"Type"}="String";
					$vcf_header{"INFO"}{"GNOMEN"}{"Description"}="\"GNOMEN hgvs gene nomenclature related to a transcript (TNOMEN)$description_plus\"";
				};#if


			};#if

			$is_calculated=1;
		};#if



		###########
		# BARCODE #
		###########

		if (uc(trim($one_calculation)) eq "BARCODE") {
		#if (uc(trim($one_calculation)) eq "VAF") {
			print "# Calculation '$one_calculation'...\n" if $DEBUG;

			my @BARCODE_LIST; #=(0.99,0.999);
			#my $VAF_min="";
			#my $VAF_max="";
			#my $VAF_average="";
			#my $VAF_median="";

			my @sample_list_header=split(/\t/,$$variant_values{"SAMPLE_LIST_HEADER"});
			#print "SAMPLE_LIST_HEADER: ".$$variant_values{"SAMPLE_LIST_HEADER"}."\n" if $DEBUG;
			print "SAMPLE_LIST_HEADER: @sample_list_header\n" if $DEBUG;

			if ($force || !defined $$variant_values{"INFOS"}{"BARCODE"}) {


				print "TRIONAME:@trio_array \n" if $DEBUG;

				my $father_name=$trio_array[0];
				my $mother_name=$trio_array[1];
				my $child_name=$trio_array[2];

				my $father_barcode="";
				my $mother_barcode="";
				my $child_barcode="";
				my $trio_barcode="";

				for my $sample_name (@sample_list_header) {
					print "SAMPLE: $sample_name ".$$samples{$sample_name}{"GT"}."\n" if $DEBUG;

					my $BARCODE="";
					my $GT=$$samples{$sample_name}{"GT"};

					if (defined $GT) {

						switch ($GT) {

							case "1/1"
							{
								$BARCODE="2";
								break;
							} # case
							case ["0/1","1/0"]
							{
								$BARCODE="1";
								break;
							} # case
							case ["0/0","./."]
							{
								$BARCODE="0";
								break;
							} # case
							else
							{
								$BARCODE="?";
							}

						};# switch
					};#if


					if ($BARCODE ne "") {
						push @BARCODE_LIST, $BARCODE ;
					} else {
						$BARCODE="?";
					};# if

					if ("$father_name" eq "$sample_name") {
						$father_barcode=$BARCODE;
					};#if
					if ("$mother_name" eq "$sample_name") {
						$mother_barcode=$BARCODE;
					};#if
					if ("$child_name" eq "$sample_name") {
						$child_barcode=$BARCODE;
					};#if

					print "GT: $GT > ".$BARCODE." \n" if $DEBUG;


				};#for

				# TRIO
				$trio_barcode="$father_barcode$mother_barcode$child_barcode";
				print "TRIO:$father_barcode$mother_barcode$child_barcode  FATHER:$father_barcode  MOTHER:$mother_barcode  CHILD:$child_barcode \n" if $DEBUG;


				if (length($trio_barcode)==3) {

					my $trio_variant_type="unknown";

					switch ($trio_barcode) {

						case "001"
						{
							$trio_variant_type="denovo";
							break;
						} # case
						case ["011","101","111","021","201","121","211"]
						{
							$trio_variant_type="dominant";
							break;
						} # case
						case ["112","212","122","222"]
						{
							$trio_variant_type="recessive";
							break;
						} # case

					};# switch

					print "trio_variant_type=$trio_variant_type \n" if $DEBUG;

					# BARCODE in INFO
					$annotation_output{$chr}{$pos}{$ref}{$alt}{"INFOS"}{"trio_variant_type"}=$trio_variant_type;
					# ADD to HEADER
					$vcf_header{"INFO"}{"trio_variant_type"}{"Number"}="1";
					$vcf_header{"INFO"}{"trio_variant_type"}{"Type"}="String";
					$vcf_header{"INFO"}{"trio_variant_type"}{"Description"}="\"Trio variant type from VaRank BARCODE$description_plus\"";

				};#if

				if (@BARCODE_LIST) {
					# BARCODE in INFO
					$annotation_output{$chr}{$pos}{$ref}{$alt}{"INFOS"}{"BARCODE"}=join("",@BARCODE_LIST);
					# ADD to HEADER
					$vcf_header{"INFO"}{"BARCODE"}{"Number"}="1";
					$vcf_header{"INFO"}{"BARCODE"}{"Type"}="String";
					$vcf_header{"INFO"}{"BARCODE"}{"Description"}="\"VaRank BARCODE$description_plus\"";
				};#if




			};#if

			print "@BARCODE_LIST\n" if $DEBUG;

			$is_calculated=1;
		};#if



		# Is calculated ?
		if ($is_calculated) {
			print "# Calculation '$one_calculation' executed\n\n" if $DEBUG;
		} else {
			print "# Calculation '$one_calculation' NOT executed\n\n" if $DEBUG;
		};#if


	};#foreach


	#print $annotation_output{$chr}{$pos}{$ref}{$alt}{"FORMAT"}."\n" if $DEBUG;


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
				#print "IS @infos_split\n";
				my %infos;
				foreach my $info (@infos_split) {
					#print "$info\n";
					#my @info_split=split("=",$info,2);
					my @info_split=split("=",$info,2);
					$infos{$info_split[0]}=$info_split[1];
				};#foreach
				#print Dumper(\%infos); # if $DEBUG;
				if (!exists($vcf_header{$type}{$infos{"ID"}})) { # Updated!
				#if (!exists($vcf_header{$1}{$infos{"ID"}}) || $force) { # Updated!
					if ($type eq "INFO") {
						# INFO header mandatory fields Consolidation
						if (!exists($vcf_header{$type}{$infos{"ID"}}{"Number"})) { $vcf_header{$type}{$infos{"ID"}}{"Number"}="."; };#if
						if (!exists($vcf_header{$type}{$infos{"ID"}}{"Type"})) { $vcf_header{$type}{$infos{"ID"}}{"Type"}="String"; };#if
						if (!exists($vcf_header{$type}{$infos{"ID"}}{"Description"})) { $vcf_header{$type}{$infos{"ID"}}{"Description"}="\"unknown\""; };#if
					};#if
					print "UPDATE! ".$infos{"ID"}."\n" if $DEBUG;
					while ((my $info_var, my $info_val) = each(%infos)){
						$vcf_header{$type}{$infos{"ID"}}{$info_var}=$info_val if ($info_var ne "ID");
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
		print FILE_OUTPUT $vcf_header_information;

		foreach my $type (sort { $vcf_header{$a} <=> $vcf_header{$b} or $a cmp $b } keys %vcf_header) {
			foreach my $ID (sort { $vcf_header{$type}{$a} <=> $vcf_header{$type}{$b} or $a cmp $b } keys %{$vcf_header{$type}}) {
				my $line="##$type=<ID=$ID";
				if ($type eq "INFO") {
					#print "\n" if $DEBUG;
					#print "$type}{$ID\n" if $DEBUG;
					#print "   ".$vcf_header{$type}{$ID}{"Number"}."\n" if $DEBUG;
					#print "   ".$vcf_header{$type}{$ID}{"Type"}."\n" if $DEBUG;
					#print "   ".$vcf_header{$type}{$ID}{"Description"}."\n" if $DEBUG;
					#print "\n" if $DEBUG;
					# INFO header mandatory fields Consolidation
					if (!exists($vcf_header{$type}{$ID}{"Number"}) || trim($vcf_header{$type}{$ID}{"Number"}) eq "") { $vcf_header{$type}{$ID}{"Number"}="."; };#if
					if (!exists($vcf_header{$type}{$ID}{"Type"}) || trim($vcf_header{$type}{$ID}{"Type"}) eq "") { $vcf_header{$type}{$ID}{"Type"}="String"; };#if
					if (!exists($vcf_header{$type}{$ID}{"Description"}) || trim($vcf_header{$type}{$ID}{"Description"}) eq "") { $vcf_header{$type}{$ID}{"Description"}="\"unknown\""; };#if
					#print " > ".$vcf_header{$type}{$ID}{"Number"}."\n" if $DEBUG;
					#print " > ".$vcf_header{$type}{$ID}{"Type"}."\n" if $DEBUG;
					#print " > ".$vcf_header{$type}{$ID}{"Description"}."\n" if $DEBUG;
					#print "\n" if $DEBUG;
				};#if

				while ((my $info_var, my $info_val) = each(%{$vcf_header{$type}{$ID}})){
					if (trim($info_val) ne "") {
						$line.=",$info_var=$info_val";
						#print ",$info_var=$info_val\n" if $DEBUG;
					};#if
				};#while
				$line.=">\n";
				print FILE_OUTPUT $line;
			};#foreach
		};#foreach

	#if (substr($line_content,0,1) eq "#") { # VCF Header
		#check the line in order to insert the description of the new info into the VCF INFO column
		#if (substr($line_content,0,7) eq "##INFO=") {
		#	$line_info_read=1;
		#} elsif ($line_info_read==1) {
		#	#print FILE_INPUT_annotated "$vcf_header";
		#	$line_info_read=0;
		#};#if
		#if (substr($line_content,0,6) eq "#CHROM") {
		# 	#Add Header
		#	print FILE_OUTPUT $PZHeader;
		#};#if
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
				$annotation_input{$chr}{$pos}{$ref}{$alt}{"INFOS"}{$annotation_name}="$annotation_result";
			};#if

			if (!exists $annotation_input{$chr}{$pos}{$ref}{$alt}{"INFOS"}{$annotation_name}
				|| $annotation_result ne $annotation_input{$chr}{$pos}{$ref}{$alt}{"INFOS"}{$annotation_name}) {
				my $original_result=trim($annotation_input{$chr}{$pos}{$ref}{$alt}{"INFOS"}{$annotation_name});
				my $sep=((trim($annotation_input{$chr}{$pos}{$ref}{$alt}{"INFOS"}{$annotation_name}) eq "")?"":",");
				$annotation_input{$chr}{$pos}{$ref}{$alt}{"INFOS"}{$annotation_name}="$original_result$sep$annotation_result";
			};#if
		};#while

		# FORMAT
		$line_content_split[8]=$annotation_output{$chr}{$pos}{$ref}{$alt}{"FORMAT"};
		#print "10: ".$line_content_split[8]." ".$line_content_split[9]."\n" if $DEBUG;
		#print "11: ".$line_content_split[8]." ".$line_content_split[10]."\n" if $DEBUG;
		#print "12: ".$line_content_split[8]." ".$line_content_split[11]."\n" if $DEBUG;
		my $sample_i=9;
		foreach my $sample_name (split(/\t/,$annotation_output{$chr}{$pos}{$ref}{$alt}{"SAMPLE_LIST_HEADER"})) {
			my $sample_quality="";
			foreach my $quality (split(/:/,$annotation_output{$chr}{$pos}{$ref}{$alt}{"FORMAT"})) {

				#print "$sample_i $sample_name $quality ".$annotation_output{$chr}{$pos}{$ref}{$alt}{"SAMPLES"}{$sample_name}{$quality}."\n" if $DEBUG;
				my $sep=((trim($sample_quality) eq "")?"":":");
				$sample_quality.=$sep.$annotation_output{$chr}{$pos}{$ref}{$alt}{"SAMPLES"}{$sample_name}{$quality};

			}
			#print "$sample_i $sample_name ".$line_content_split[8]." $sample_quality\n" if $DEBUG;
			$line_content_split[$sample_i]=$sample_quality;
			$sample_i++;
		}

		#print Dumper(\%{$annotation_input{$chr}{$pos}{$ref}{$alt}}) if $DEBUG;

		# OUTPUT
		if (!$hard_filter_remove) {
			#my $variant_annotation=join(";", map { "$_=$annotation_input{$chr}{$pos}{$ref}{$alt}{'INFOS'}{$_}" } keys %{$annotation_input{$chr}{$pos}{$ref}{$alt}{"INFOS"}	});

			my $variant_annotation="";
			my @varval_INFOS=keys(%{$annotation_input{$chr}{$pos}{$ref}{$alt}{"INFOS"}}) if $DEBUG;
			#print "@varval_INFOS\n" if $DEBUG;
			foreach my $var (keys(%{$annotation_input{$chr}{$pos}{$ref}{$alt}{"INFOS"}})) {
				my $val=$annotation_input{$chr}{$pos}{$ref}{$alt}{'INFOS'}{$var};
				#print "   $var = $val\n" if $DEBUG;
				#if ($var ne ".") {
					my $sep=""; $sep=";" if $variant_annotation ne "";
					if (defined $val) {
						$variant_annotation.="$sep$var=$val";
					} else {
						$variant_annotation.="$sep$var";
					};#if
				#};#if

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
			#print "$variant_annotation\n" if $DEBUG;


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
