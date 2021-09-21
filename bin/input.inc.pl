

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

=item B<--config_prioritization|--config_prioritization_file|--config_filter|config_filter_file=<file>>

Configuration file for filter/prioritization parameters (default 'config.prioritization.ini' or 'config.filter.ini').

=item B<--annovar_folder=<folder>>

Folder with ANNOVAR scripts.

=item B<--annovar_databases=<folder>>

Folder with ANNOVAR databases.


=item B<--snpeff_jar=<file>>

snpEff JAR file.

=item B<--snpeff_databases=<folder>>

Folder with snpEff databases.

=item B<--java=<file>>

java binary (needed if snpeff option on). default "java"

=item B<--java_flags=<string>>

java flags  (needed if snpeff option on, especially for configure proxy for databases donwload). default ""

=back

=head2 COMMON PARAMETERS

=over 2

=item B<--force!>

Force calculation, annotation...

=item B<--output=<file>>

VCF Output file

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

=item B<--calculation=<string>>

calculation to do.
Example: "VAF", "VAF,cNomen,pNomen"

Currently available :

VAF: add VAF calculation for each sample, depending on information provided by callers (by order: FREQ, DP4, AD)

VAF_STATS: add VAF statistics on INFO field, from VAF calculation for each sample, depending on information provided by callers (by order: VAF, FREQ, DP4, AD)

CALLING_QUALITY: Calling quality (FORMAT/*) of all samples in case of multiSample VCF, or all pipelines in case of multipipeline VCF

CALLING_QUALITY_EXPLODE: Explode all Calling quality (FORMAT/*) in multiple fields in INFOS

NOMEN: Find the NOMEN from HGVS annotation. Depend on transcript of reference. If no transcript of reference, first transcript. This option create annotations on INFO field: NOMEN (full HGVS annotation), CNOMEN (DNA level mutation "c."), RNOMEN (RNA level mutation "c."), NNOMEN (non-coding mutation "c."), PNOMEN (Protein level mutation "p."), TNOMEN (transcript), TVNOMEN (transcript with version), ENOMEN (exon, if any), GNOMEN (gene, if any)

BARCODE: Calculate VaRank BarCode

GENOTYPECONCORDANCE: If all samples with the same genotype

FINDBYPIPELINES: Number of pipeline calling the variant

VARTYPE: SNV if X>Y, MOSAIC if X>Y,Z or X,Y>Z, INDEL if XY>Z or X>YZ

=item B<--transcripts=<string>>

file containing default transcript for each gene.

format :

TRANSCRIPT1	GENE1

TRANSCRIPT2	GENE2

...

=item B<--nomen_fields=<string>>

List of names of annotation field to use to calculate/extract NOMEN annotation. Field order determine the priority.

Format: field1,field2,...

Examples: "hgvs", "snpeff_hgvs", "snpeff_hgvs,hgvs"

Default: "hgvs"

=item B<--nomen_pattern=<string>>

List of information to generate NOMEN annotation.

Information order determine the NOMEN pattern.

Empty information will not be used.

Format: INFO1:INFO2:...

Available information: see --calculation option

Examples: "GNOMEN:TNOMEN:ENOMEN:CNOMEN:RNOMEN:NNOMEN:PNOMEN", "GNOMEN:TNOMEN:CNOMEN:RNOMEN:PNOMEN", "GNOMEN:CNOMEN:PNOMEN"

Default: "GNOMEN:TNOMEN:ENOMEN:CNOMEN:RNOMEN:NNOMEN:PNOMEN"

=item B<--trio=<string>>

List of sample to identify a trio. Will automatically calculate VaRank barcode, and add INFO/trio_variant_type (either "denovo", "dominant", "recessive") as:

case "001": denovo

case "011","101","111","021","201","121","211": dominant

case "112","212","122","222": recessive

format: "father_sample_name,mother_sample_name,child_sample_name"

=item B<--prioritization|filter=<file>>

Filter profile, defined in the 'config_annotation' file. 'ALL' for all filters defined in the 'config_annotation' file. If the filter doesn't exist, 'ALL' will be used (if exist).

=item B<--hard!>

Remove not PASS variant in default priorisation filter (PZFlag)

=item B<--pzfields=<string>>

List of prioritisation information to show

Default: 'PZScore,PZFlag'

Format: 'field1,field2...'

Example: 'PZScore,PZFlag,PZComment,PZInfos'

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




$annotation_default="ALL";
$calculation_default="";
$sort_by_default="PZFlag,PZScore";
$order_by_default="DESC,DESC";
$PZFields_default="PZScore,PZFlag";

our %parameters = ( #
	# Main options
	'help'			=>  	0,	# Help parameter
	'man'			=>  	0,	# Man parameter
	'release'		=>  	0,	# Release parameter
	'debug'			=>  	0,	# Debug parameter
	'verbose'		=>  	0,	# Verbose parameter
	'header'		=>  	0,	# Verbose parameter

	# Configuration
	'config'		=>	'config.ini',			# Configuration file
	'config_annotation'    	=>  	"config.annotation.ini",	# Configuratin annotation file
	'config_prioritization'    	=>  	"config.prioritization.ini",	# Configuration filter file
	'tmp'			=>	undef,		# Output VCF file

	# COMMON PARAMETERS
	'force'			=>	0,		# Configuration file


	# INPUT/OUTPUT
	'input'			=>	undef,		# Input VCF file
	'output'		=>	undef,		# Output VCF file

	# ANNOTATION
	# Input
	'annotation'				=>	$annotation_default,	# Annotation sources
	#'annovar_annotation_type_default'	=>	'geneanno',		# In case of ANNOVAR code in annotation
	'annovar_annotation_type_default'	=>	'filter',		# In case of ANNOVAR code in annotation
	'AlleleFrequency'			=>	0,			# Annotation sources

	# Configuration
	'assembly'    			=>  	"",		# Assembly
	'annovar_folder'    		=>  	"",		# ANNOVAR folder including scripts
	'annovar_databases'    		=>  	"",		# ANNOVAR databases
	#'snpeff'    			=>  	0,		# snpEff annotation
	'snpeff_threads'		=>  	0,		# snpEff annotation
	'snpeff_jar'    		=>  	"",		# snpEff JAR
	'snpeff_databases'    		=>  	"",		# snpEff databases
	#'snpeff_hgvs'    		=>  	0,		# snpEff HGVS annotation in INFO:hgvs
	'snpeff_spliceSiteSize'    	=>  	3,		# snpEff additional options
	'snpeff_additional_options'    	=>  	"",		# snpEff additional options
	'snpeff_stats'    		=>  	"",		# snpEff stats
	'bcftools_stats'    		=>  	"",		# BCFTOOLS stats
	'stats'    		=>  	"",		# HOWARD stats
	'java'    			=>  	"",		# JAVA binary
	'java_flags'    		=>  	"",		# JAVA binary
	# show
	'show_annotations'		=>	0,		# Show annotations
	'show_annotations_full'		=>	0,		# Show annotations FULL

	# CALCULATION
	# Input
	'calculation'		=>	$calculation_default,	# Calculation to do
	'transcripts'		=>	'',			# File with default transcripts by gene
	'nomen_fields'		=>	'hgvs',			# Fields to calculate NOMEN
	'nomen_pattern'		=>	'GNOMEN:TNOMEN:ENOMEN:CNOMEN:RNOMEN:NNOMEN:PNOMEN',			# Information to generate NOMEN
	'trio'			=>	'',			# Trio identification

	# PRIORITIZATION
	# Input
	'filter'	=>	'ALL',			# Filter
	'pzfields'	=>	$PZFields_default,	# Filter fields
	'hard'		=>	0,			# Hard filtering parameters
	'prioritization_score_mode'		=>	"howard",			# Prioritization score mode parameters

	# TRANSLATION
	# Input
	'format'	=>	'tab',			# Output format
	'fields'	=>	$annotation_default,	# Fields to show
	'sort_by'	=>	$sort_by_default,	# Sort variants by an annotation
	'order_by'	=>	$order_by_default,	# Order variants by an annotation, ASC or DESC
	'columns'	=>	'',			# Colmuns to add in case of tab delimiter format

	# OTHER
	'split'		=>	'10000',		# Split variants, only for HOWARD.sh


);

our @options=(
	# Main options
	'help|h|?',		# Help
	'man',			# Man
	'release',		# Release
	'debug',		# Debug
	'verbose',		# Verbose
	'header',		# Verbose

	# CONFIGURATION
	'config|config_file=s',				# Configuration file
	'config_annotation|config_annotation_file=s',	# Configuratin annotation file
	'config_prioritization|config_prioritization_file|config_filter|config_filter_file=s',		# Configuration filter file
	'tmp=s',		# Configuration filter file

	# COMMON PARAMETERS
	'force!',		# Input file

	# INPUT/OUTPUT
	'input|input_file|vcf=s',	# Input file
	'output|output_file=s',		# Output file

	# ANNOTATION
	# Input
	'annotation=s',				# Annotations
	'annovar_annotation_type_default=s',	# In case of ANNOVAR code in annotation
	'AlleleFrequency!',			# Allele Frequency Calculation
	# Configuration
	'assembly=s',			# Assembly
	'annovar_folder=s',		# ANNOVAR folder including scripts
	'annovar_databases=s',		# ANNOVAR databases
	#'snpeff!',			# snpEff annotation
	'snpeff_threads=i',		# snpEff JAR
	'snpeff_jar=s',			# snpEff JAR
	'snpeff_databases=s',		# snpEff databases
	#'snpeff_hgvs!',			# snpEff HGVS annotation in INFO:hgvs
	'snpeff_spliceSiteSize=i',	# snpEff additional options
	'snpeff_additional_options=s',	# snpEff additional options
	'snpeff_stats=s',		# snpEff stats files
	'bcftools_stats=s',		# BCFTOOLS stats files
	'stats=s',		# HOWARD stats files
	'java=s',			# JAVA binary
	'java_flags=s',			# JAVA binary
	# show
	'show_annotations!',		# Show annotations
	'show_annotations_full!',	# Show annotations FULL

	# CALCULATION
	# Input
	'calculation=s',	# Calculation to do
	'transcripts=s',	# File with default transcripts by gene
	'nomen_fields=s',	# Fields to calculate NOMEN
	'nomen_pattern=s',	# Information to generate NOMEN
	'trio=s',		# Trio identification

	# PRIORITIZATION
	# Input
	'filter|prioritization=s',		# Filter
	'pzfields=s',		# Filter fields
	'hard!',		# Hard Filtering
	'prioritization_score_mode=s',		# # Prioritization score mode
	
	# TRANSLATION
	# Input
	'format|translation=s',		# Output format
	'fields=s',		# Fields to show
	'sort_by=s',		# Sort variants by an annotation
	'order_by=s',		# Order variants by an annotation, ASC or DESC
	'columns=s',		# Colmuns to add in case of tab delimiter format

	# OTHER
	'split=s',		# Output format
);



## Catch Options and put into parameters
GetOptions (\%parameters,@options,@common_options)  or pod2usage();

#print Dumper(\%parameters);
