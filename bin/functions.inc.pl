#!/usr/bin/perl
##############################
# Functions library          #
# Author: Antony Le Béchec   #
# Date: 23/08/2018           #
# Copyright: IRC             #
##############################

use POSIX qw(locale_h strtod ceil);
use Time::localtime;

use Scalar::Util qw/reftype/;
use List::Util qw(sum min max);
#use POSIX;

sub read_ini {
# Function read_ini
# Read ini file
#$_[0]: ini file
#$_[1]: header case (either NULL, "lower" or "upper")

	# Parameters
	my $ini_file=$_[0];
	my $case=lc($_[1]);
	my %ini_infos;
	my $header="";

        # Initialisation
	if (-e $ini_file) { #if file exist

		# Open the file
		open(FILE_PROJECT_INFOS, $ini_file) || die "Problem to open the file: $!";

		# Read the file		
		while(<FILE_PROJECT_INFOS>) {

			# delete \n character
			chomp;

			# ignore blank lines
			next if /^\s*$/;

			# ignore commented lines
			next if /^\s*\;|#/;

			# check for section header
			if ( /^\s*\[([^\]]+)\]/ ) {

				# read header
				$header = $1;
				
				# case
				if ($case ne "") {
					$header=lc($header) if ($case eq "lower");
					$header=uc($header) if ($case eq "upper");
					#if ($case eq "lower") {	$header=lc($header); };#if
					#if ($case eq "upper") {	$header=uc($header); };#if
				};#if

				# remove leading and trailing white spaces
				$header =~ s/^\s*|\s*$//g;
				$ini_infos{$header} = {};

			};#if

			# check for section data
			#if (/(\S+)\s*[=:\t]\s*(.*?)\s*$/) {
			#if (/(\S+)(\[\])*\s*[=:\t]\s*(.*?)\s*$/) {
			#if (/((\s)*(?([^\=^\s^\n]+))[\s^\n]*\=(\s)*(?([^\n^\s]+(\n){0,1})))/) {
			#if (/^\W*(\w+)=?(\w+)\W*(#.*)?$/) {
			#if (/(\S+)(\(.|\[\]\)*)[=:\t]\s*(.*?)\s*$/) {
			#if (/(.*)[=:\t]\s*(.*?)\s*$/) {
			#if (/^(.*)[=:\t](.*)$/) {
			#if (/\s*(.+?)\s*=\s*(.+)/) { # GOOD
			#if (/^([^=;\r\n]+)=([^;\r\n]*)/) { # VERY GOOD
			if (/^([^=;\r\n]+)=([^;\r\n]*)/) { # VERY GOOD
				#print "(1)$1\t(2)$2\t(3)$3\t(4)$4\n";

				# read line
				$line=$_;

				# split line
				my @line_split1=split("",$line);
				my @line_split2=split("=",$line);
				my $var=$1;
				my $val=$2;
				# Remove comment
				@val_split=split(/([\t;])/,$val);
				$val=trim($val_split[0]);
				#print "VAR $var\t=\t VAL $val\n";

				# case
				if ($case ne "") {
					$var=lc($var) if ($case eq "lower");
					$var=uc($var) if ($case eq "upper");
				};#if

				# add info
				my $value=$val;
				#my $variable=$var;
				#$variable=~s/\[\]$//g;

				if (trim($var) =~ /.*\[\]$/) { # Value is a hash
					$ini_infos{$header}{$var}{$val}=$val;	
				} else {
					$ini_infos{$header}{$var}=$val;	
				};#if

			};#if

		};#while

		# close file
		close(FILE_PROJECT_INFOS);

		# retunr
		return %ini_infos;

	} else {

		# No ini file
		print "File '$ini_file' DOES NOT exist!\n";
		exit 1;

	};#if

}


# function read_vcf
sub read_vcf {
# $_[0]: VCF file
# $_[]:
# return $annotation_output (hash ref), $annotations (hash ref)

	# input
	my $input_file=$_[0];	# Input file in VCF format
	my $output_type=$_[1];	# Type of output. Default "variants", else "header"
	$output_type=(($output_type eq "")?"variants":$output_type);

	# Output
	my %annotation_output;
	my %annotations;

	# Open file
	open(VCF, $input_file) || die "Problem to open the file '$input_file': $!";
	my $line=0;
	my $data_header=0;
	my $variant_group=0;
	my $variant_index=0;
	my $header_VCF;
	my $header_VCF_samples;
	my $nb_sample=0;
	my %VCF_header;
	my $vcf_header_INFOS="";
	my %vcf_header;
	my $vcf_header_information;
	my %sample_hash;
	my %annotation_output;

	# Read
	while(<VCF>) {

		# Clean line
		chomp; #delete \n character
		my $line++;
		my $line_content=$_;

		#print "$line_content\n" if $DEBUG;
		

		# HEADER
		if(grep(/^##.*$/i,$line_content)) {
			$header_VCF.=trim($line_content)."\n";
			# example: ##contig=<ID=Y,length=59373566,assembly=hg19>
			if ($line_content =~ /##(.*)=<(.*)>/) {
				my $header_type=$1;
				#print "$1\t$2\n" if $DEBUG;	
				my %header_variable_infos = map{split /=/, $_}(split /,/, $2);
				my $header_variable_id=$header_variable_infos{"ID"};
				#print "$header_type\t$header_variable_id\n" if $DEBUG;	
				$VCF_header{$header_type}{$header_variable_id}=\%header_variable_infos;
				#print Dumper(\%header_variable_infos) if $DEBUG;
				#print Dumper($VCF_header{$header_type}{$header_variable_id}) if $DEBUG;
				#print "test=".$$VCF_header{$header_type}{$header_variable_id}{"ID"}."\n" if $DEBUG;
			};#if

			if ($line_content=~ /^##(.*)=<(.*)>$/) {
				my $type=$1;
				if ($type eq "INFO" || $type eq "FORMAT" || $type eq "contig" || $type eq "FILTER") {
					my @infos_split=split(",",$2);
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


		# HEADER Variant line
		} elsif (grep(/^#.*$/i,$line_content)) {
			$header_VCF.=trim($line_content)."\n";
			@line_content_split=split("\t",trim($line_content));
			my $size_array = @line_content_split;
			$nb_sample=$size_array-9;
			for ($i=9;$i<@line_content_split;$i++) {
				my $sep=($i>9)?"\t":"";
				$header_VCF_samples.=$sep.$line_content_split[$i];
				$sample_hash{$i}=$line_content_split[$i];
			};#if

		# VARIANT
		} else {

			# SPLIT line
			my @line_content_split=split("\t",trim($line_content));
			my $chr=$line_content_split[0];
			my $chr_num=$chr; $chr_num =~ s/chr//g;
			my $pos=$line_content_split[1];
			my $id=$line_content_split[2]; $id =~ s/;/,/g;
			my $ref=$line_content_split[3];
			my $alt=$line_content_split[4];
			my $qual=$line_content_split[5];
			my $filter=$line_content_split[6];
			my $infos=$line_content_split[7];
			my $format=$line_content_split[8];
			my $sample=$line_content_split[9];
		
			# SPLIT sample col
			my $sample_list="";
			my $sample_list_concat="";
			my $sample_list_concat2="";
			my $sample_list_concat2_header="";
			for ($i=9;$i<@line_content_split;$i++) {
				my $sep=($i>9)?"\t":"";
				$sample_list.=$sep.$sample_hash{$i}.":".$line_content_split[$i];
				#$sample_list.=$sep.$line_content_split[$i];
				my $sep=($i>9)?"|":"";
				my $sep2=($i>9)?"\t":"";
				$sample_list_concat.=$sep.$sample_hash{$i}.":".$line_content_split[$i];
				$sample_list_concat2.=$sep2.$line_content_split[$i];
				$sample_list_concat2_header.=$sep2.$sample_hash{$i};
				#$sample_list_concat.=$sep.$line_content_split[$i];
				#print "Sample List $i:\t$sample_list\n";
			};#for
			#$format="SAMPLE:".$format;

			# SAMPLES QUALITY
			my @format_split=split(/:/, $format);
			my %sample_quality;
			for ($i=9;$i<@line_content_split;$i++) {
				my $sample_quality=$line_content_split[$i];
				@sample_quality_split=split(/:/, $sample_quality);
				my $nb_quality=0;
				foreach my $quality (@sample_quality_split) {
					$sample_quality{$sample_hash{$i}}{$format_split[$nb_quality]}=$quality;
					$nb_quality++;
				};#foreach
			};#for

			# INDEX
			#$variant_index++;
		
			#foreach $alt (split(",",$alts)) {

				# Ref / Alt transformation

				# HASH VCF
				$variants_VCF{$variant_index}{$chr}{$pos}{$ref}{$alt}{"VCF"}=trim($line_content);
				$variants_VCF{$variant_index}{$chr}{$pos}{$ref}{$alt}{"VCFinfos"}=$infos;
				$variants_VCF{$variant_index}{$chr}{$pos}{$ref}{$alt}{"VCFformat"}=$format;
				$variants_VCF{$variant_index}{$chr}{$pos}{$ref}{$alt}{"VCFsample_list"}=$sample_list;
		
				# HASH Variants
				$annotation_output{$chr}{$pos}{$ref}{$alt}{"CHROM"}=$chr;		# CHR
				$annotation_output{$chr}{$pos}{$ref}{$alt}{"CHROM_NUM"}=$chr_num;	# CHR
				$annotation_output{$chr}{$pos}{$ref}{$alt}{"POS"}=$pos;			# POS
				$annotation_output{$chr}{$pos}{$ref}{$alt}{"ID"}=$id;			# ID
				$annotation_output{$chr}{$pos}{$ref}{$alt}{"REF"}=$ref;			# REF
				$annotation_output{$chr}{$pos}{$ref}{$alt}{"ALT"}=$alt;			# ALT
				$annotation_output{$chr}{$pos}{$ref}{$alt}{"QUAL"}=$qual;		# QUAL
				$annotation_output{$chr}{$pos}{$ref}{$alt}{"FILTER"}=$filter;		# FILTER
				$annotation_output{$chr}{$pos}{$ref}{$alt}{"FORMAT"}=$format;		# FORMAT
				$annotation_output{$chr}{$pos}{$ref}{$alt}{"SAMPLES"}=\%sample_quality;	# SAMPLES
				$annotation_output{$chr}{$pos}{$ref}{$alt}{"SAMPLE_LIST"}=$sample_list;	# SAMPLE_LIST
				$annotation_output{$chr}{$pos}{$ref}{$alt}{"SAMPLE_LIST_CONCAT"}=$sample_list_concat;	# SAMPLE_LIST_CONCAT
				$annotation_output{$chr}{$pos}{$ref}{$alt}{"SAMPLE_LIST_CONCAT2"}=$sample_list_concat2;	# SAMPLE_LIST_CONCAT
				$annotation_output{$chr}{$pos}{$ref}{$alt}{"SAMPLE_LIST_HEADER"}=$sample_list_concat2_header;	# SAMPLE_LIST_CONCAT
				#print "###### SAMPLE_LIST_CONCAT2 $chr}{$pos}{$ref}{$alt: $sample_list_concat2\n" if $DEBUG; # =$info_value
				
			
		
				# INFOS col
				my $infos=$line_content_split[7];
				#print "###### $infos\n" if $DEBUG; # =$info_value
				my @infos_split=split(";",$infos);
				#print "### @infos_split\n" if $DEBUG; # =$info_value
				foreach $info (@infos_split) {
					my @info_split=split("=",$info,2);
					my $info_name=$info_split[0];
					my $info_value=$info_split[1];
					#print "# $info_name\n" if $DEBUG; # =$info_value
					#if (!defined $info_value) {$info_value=$info_name};#if
					if ($info_name ne ".") {
						$annotation_output{$chr}{$pos}{$ref}{$alt}{"INFOS"}{$info_name}=$info_value;
					};#if
				};#foreach
			
				# ANNOTATIONS
				#while ( my ($annotation_name, $annotation_value) = each(%{$annotation_output{$chr}{$pos}{$ref}{$alt}}) ) {
				#	$annotations{$annotation_name}{$annotation_value}++;
				#};#while

			#};#foreach

		};#if
	
	};#while

	#print $header_VCF;

	#my %return;
	#$return{"annotation_output"}=\%annotation_output;
	#$return{"annotations"}=\%annotations;
	#return %return;
	if ($output_type eq "header_flat") {
		return $header_VCF;
	} elsif ($output_type eq "header") {
		#print "truc\n\n\n";
		#print Dumper($VCF_header{"FORMAT"});
		#return %vcf_header;
		return %VCF_header;
	} else {
		return %annotation_output;
	};#if


}


#Function trim
#Perl trim function to remove whitespace from the start and end of the string
sub trim($) {
	my $string = shift;
	$string =~ s/^\s+//;
	$string =~ s/\s+$//;
	return $string;
}


#Function getnum DOESNOTWORK
#tranform an string on integer 
sub getnum {
    #use POSIX;
    my $str = shift;
    $str =~ s/^\s+//;
    $str =~ s/\s+$//;
    
    $! = 0;
    my($num, $unparsed) = strtod($str);
    if (($str eq '') && ($unparsed != 0) && $!) {
        return undef;
    } else {
        return $num;
    }
}

sub extract_number {
	my $str = shift;
	$str =~ s/[^0-9\.]//g;
	if (nifty_number($str)) {
		return $str;
	} else {
		return 0;
	};#if
}

#Function timestamp
#return the date in mutiple formats
sub timestamp {
#parameters:
#*format of the date:
# - "date" as 'yyyy-mm-dd'
# - "datamin" as 'yyyymmdd'
# - default yyy-mm-dd_hh:mm:ss
  my $t = localtime;
	if (lc($_[0]) eq "date") {
		return sprintf( "%04d-%02d-%02d", $t->year + 1900, $t->mon + 1, $t->mday);
	} elsif (lc($_[0]) eq "datemin") {
		return sprintf( "%04d%02d%02d%02d%02d%02d", $t->year + 1900, $t->mon + 1, $t->mday,$t->hour, $t->min, $t->sec);
	} else {
		return sprintf( "%04d-%02d-%02d %02d:%02d:%02d",
			  $t->year + 1900, $t->mon + 1, $t->mday,
			  $t->hour, $t->min, $t->sec );
	};#if
}

#Function in_array
#return true if an element is in an array
sub in_array {
#parameters:
#$_[0]: The array in which you look in
#$_[1]: The element you look in the array
#$_[3]: case sensitivity (default FALSE)
#usage:  if(in_array(\@arr,'Amlan'))...
	my ($arr,$search_for,$case) = @_;
	($case eq "")?$case=1:$case=$case;
	#my %items = map {$_ => 1} @$arr;
	#return (exists($items{$search_for}))?1:0;
	my %items = map {($case)?$_:lc($_) => 1} @$arr;
	#return (exists($items{($case)?$search_for:lc($search_for)}))?1:0;
	return (exists($items{($case)?$search_for:lc($search_for)}))?1:0;
}

#Function in_array_array
#return true if at least one element in an array is inan array
sub in_array_array {
#parameters:
#$_[0]: The array in which you look in
#$_[1]: The array of element you look in the array
#usage:  if(in_array(\@arr,'Amlan'))...
	#my $arr=$_[];
	my ($arr,@arr2) = @_;
	my %items = map {$_ => 1} @$arr;
	#print "TEST"."@arr2"."\n";
	#foreach $search_for (@$arr2) {
	foreach $search_for (@arr2) {
		#print "$search_for == ".$items{$search_for}."\n";
		if (exists($items{$search_for})) {
			return 1;
		};#if
	};#foreach
	return 0;
}

#Function index_array
#return the index of an element in an array
sub index_array {
#parameters
#$_[0]: The array in which you look in
#$_[1]: The element you want the index
#$_[2]: case sentitivity (default false)
	my $array = $_[0];
	my @arr=@$array;
	my $element = $_[1];
	my $case_sensitivity = $_[2];
	if ($case_sensitivity) {
		$element=lc($element);
	};#if
	return first { $arr[$_] eq $element } 0..$#arr;
}

#functions for number checking
sub is_whole_number {
	$_[0] =~ /^\d+$/
}
sub is_number {
	$_[0] =~ /^\d+$/
}	
sub is_integer {
	$_[0] =~ /^[+-]?\d+$/
}
sub is_float {
	$_[0] =~ /^[+-]?\d+\.?\d*$/
}
#from an old paper-tape processing machine from Mark Biggar:
sub nifty_number {
    $_[0] =~ m{   # YANETUT
		^ ( [+-]? )
		  (?= \d 
		    | \.\d 
		  )
		  \d*
		  ( \. \d* ) ?
		  (   [Ee] ( [+-]? \d+ ) ) ?
		$
    }x
}


#function providing all posible combination of element of a list
sub combo {
#Usage:
#my $iter = combo( 5 , 1 .. 50  );
#while ( my @combo = $iter->() ) {
#    print "@combo\n";
#}
    my $by = shift;
    return sub { () } if ! $by || $by =~ /\D/ || @_ < $by;
    my @list = @_;

    my @position = (0 .. $by - 2, $by - 2);
    my @stop     = @list - $by .. $#list;
    my $end_pos  = $#position;
    my $done     = undef;

    return sub {
        return () if $done;
        my $cur = $end_pos;
        {
            if ( ++$position[ $cur ] > $stop[ $cur ] ) {
                $position[ --$cur ]++;
                redo if $position[ $cur ] > $stop[ $cur ];
                my $new_pos = $position[ $cur ];
                @position[ $cur .. $end_pos ] = $new_pos .. $new_pos + $by;
            }
        }
        $done = 1 if $position[0] == $stop[0];
        return @list[ @position ];
    }
}

#Function command_fc
#return execution tim (in second) of a system command
#add the command in a log file, as well as a comment if any
sub command_fc {
#parameters:
#$_[0]: The command line to launch
#$_[1]: The logfile to fill in
#$_[1]: The comment to add to the log file
#usage:  command_fc($command[,$logFile][,$comment])
	my $command=$_[0];	# Command to run
	my $logFile=$_[1];	# Log file
	my $comment=$_[2];	# Comment
	if (trim($command) eq "") {
		return false;
	} else {
		if (trim($logFile) ne "") {
			if (trim($comment) ne "") {
				system "echo \"# $comment\" >> $logFile";
			};#if
			system "echo \"#COMMAND[".timestamp()."]: $command\" >> $logFile";
		};#if
		print "$command\n" if $DEBUG;
		my $start=time;	
		system $command;
		my $end=time;
		if (trim($logFile) ne "") {
			system "echo \"#EXECTIME = $time sec\" >> $logFile";
		};#if
		return ($end-$start);
	};#if
}

sub mean { return @_ ? sum(@_) / @_ : 0 }

sub median {
  sum( ( sort { $a <=> $b } @_ )[ int( $#_/2 ), ceil( $#_/2 ) ] )/2;
}

sub uniq {
    my %seen;
    grep !$seen{$_}++, @_;
}

return 1;
