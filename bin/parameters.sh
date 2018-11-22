



our %parameters = ( #
	# Main options 
	'help'			=>  	0,	# Help parameter
	'man'			=>  	0,	# Man parameter
	'release'		=>  	0,	# Release parameter
	'debug'			=>  	0,	# Debug parameter
	'verbose'		=>  	0,	# Verbose parameter
	
	# Configuration
	'config'		=>	'config.ini',			# Configuration file
	'config_annotation'    	=>  	"config.annotation.ini",	# Configuratin annotation file
	
	# INPUT/OUTPUT
	'input'			=>	undef,		# Input VCF file
	'output'		=>	undef,		# Output VCF file
	
	# ANNOTATION
	# Input
	'annotation'				=>	'ALL',		# Annotation sources
	'annovar_annotation_type_default'	=>	'geneanno',	# In case of ANNOVAR code in annotation
	'AlleleFrequency'			=>	0,		# Annotation sources
	
	# Configuration
	'assembly'    		=>  	"",		# Assembly
	'annovar_folder'    	=>  	"",		# ANNOVAR folder including scripts
	'annovar_databases'    	=>  	"",		# ANNOVAR databases
	'snpeff'    		=>  	0,		# snpEff annotation
	'snpeff_jar'    	=>  	"",		# snpEff JAR
	'snpeff_databases'    	=>  	"",		# snpEff databases
	'snpeff_hgvs'    	=>  	0,		# snpEff HGVS annotation in INFO:hgvs
	'snpeff_stats'    	=>  	"",		# snpEff stats
	'java'    		=>  	"java",		# JAVA binary
	# show
	'show_annotations'	=>	0,
	
	# PRIORITIZATION
	# Input
	'filter'	=>	'ALL',	# Filter sources
	'hard'		=>	0,		# Hard filterind parameters
	
	# TRANSLATION
	# Input
	'format'	=>	'tab',			# Output format
	'annotation'	=>	$annotation_default,	# Annotation to include
	'sort_by'	=>	$sort_by_default,	# Sort variants by an annotation
	'order_by'	=>	$order_by_default,	# Order variants by an annotation, ASC or DESC
	
	
};

our @options=(
	# Main options
	'help|h|?',		# Help
	'man',			# Man
	'release',		# Release
	'debug',		# Debug
	'verbose',		# Verbose
	
	# CONFIGURATION
	'config|config_file=s',				# Configuration file
	'config_annotation|config_annotation_file=s',	# Configuratin annotation file
	
	# INPUT/OUTPUT
	'input|input_file=s',		# Input file
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
	'snpeff!',			# snpEff annotation
	'snpeff_jar=s',			# snpEff JAR
	'snpeff_databases=s',		# snpEff databases
	'snpeff_hgvs!',			# snpEff HGVS annotation in INFO:hgvs
	'snpeff_stats=s',		# snpEff stats files
	'java=s',			# JAVA binary
	# show
	'show_annotations!',		# Show annotation
	
	# PRIORITIZATION
	# Input
	'filter=s',		# Filter
	'hard!',		# Hard Filtering
	
	# TRANSLATION
	# Input
	'format=s',		# Output format
	'annotation=s',		# Annotation to include
	'sort_by=s',		# Sort variants by an annotation
	'order_by=s',		# Order variants by an annotation, ASC or DESC
);    





our %parameters = ( #
	# Main options 
	'help'		=>  	0,	# Help parameter
	'man'		=>  	0,	# Man parameter
	'release'	=>  	0,	# Release parameter
	'debug'		=>  	0,	# Debug parameter
	'verbose'	=>  	0,	# Verbose parameter
	# Configuration
	'config'		=>	'config.ini',			# Configuration file
	#'config_annotation'    	=>  	"config.annotation.ini",	# Configuration annotation file
	#'config_filter'    	=>  	"config.filter.ini",		# Configuration filter file
	# Input
	'input'		=>	undef,			# Input VCF file
	'format'	=>	'tab',			# Output format
	'annotation'	=>	$annotation_default,	# Annotation to include
	'sort_by'	=>	$sort_by_default,	# Sort variants by an annotation
	'order_by'	=>	$order_by_default,	# Order variants by an annotation, ASC or DESC
	# Output
	'output'	=>	undef,		# Output VCF file
);

## Parameters definition

our @options=(
	# Main options
	'help|h|?',		# Help
	'man',			# Man
	'release',		# Release
	'debug',		# Debug
	'verbose',		# Verbose
	# Configuration
	'config|config_file=s',				# Configuration file
	#'config_annotation|config_annotation_file=s',	# Configuration annotation file
	#'config_filter|config_filter_file=s',		# Configuration filter file
	# Input
	'input|input_file=s',	# Input file
	'format=s',		# Output format
	'annotation=s',		# Annotation to include
	'sort_by=s',		# Sort variants by an annotation
	'order_by=s',		# Order variants by an annotation, ASC or DESC
	# output
	'output|output_file=s',	# Output file
);

## Catch Options and put into parameters
GetOptions (\%parameters,@options,@common_options)  or pod2usage();



























our %parameters = ( #
	# Main options 
	'help'		=>  	0,	# Help parameter
	'man'		=>  	0,	# Man parameter
	'release'	=>  	0,	# Release parameter
	'debug'		=>  	0,	# Debug parameter
	'verbose'	=>  	0,	# Verbose parameter
	# Configuration
	'config'		=>	'config.ini',			# Configuration file
	'config_annotation'    	=>  	"config.annotation.ini",	# Configuratin annotation file
	'assembly'    		=>  	"",				# Assembly
	'annovar_folder'    	=>  	"",				# ANNOVAR folder including scripts
	'annovar_databases'    	=>  	"",				# ANNOVAR databases
	'snpeff'    		=>  	0,				# snpEff annotation
	'snpeff_jar'    	=>  	"",				# snpEff JAR
	'snpeff_databases'    	=>  	"",				# snpEff databases
	'snpeff_hgvs'    	=>  	0,				# snpEff HGVS annotation in INFO:hgvs
	'snpeff_stats'    	=>  	"",				# snpEff stats
	'java'    		=>  	"java",				# JAVA binary
	# Input
	'input'					=>	undef,		# Input VCF file
	'annotation'				=>	'ALL',		# Annotation sources
	'annovar_annotation_type_default'	=>	'geneanno',	# In case of ANNOVAR code in annotation
	'AlleleFrequency'			=>	0,		# Annotation sources
	# Output
	'output'	=>	undef,		# Output VCF file
	# show
	'show_annotations'	=>	0,
);

## Parameters definition

our @options=(
	# Main options
	'help|h|?',		# Help
	'man',			# Man
	'release',		# Release
	'debug',		# Debug
	'verbose',		# Verbose
	# Configuration
	'config|config_file=s',				# Configuration file
	'config_annotation|config_annotation_file=s',	# Configuratin annotation file
	'assembly=s',					# Assembly
	'annovar_folder=s',				# ANNOVAR folder including scripts
	'annovar_databases=s',				# ANNOVAR databases
	'snpeff!',					# snpEff annotation
	'snpeff_jar=s',					# snpEff JAR
	'snpeff_databases=s',				# snpEff databases
	'snpeff_hgvs!',					# snpEff HGVS annotation in INFO:hgvs
	'snpeff_stats=s',				# snpEff stats files
	'java=s',					# JAVA binary
	# Input
	'input|input_file=s',			# Input file
	'annotation=s',				# Annotations
	'annovar_annotation_type_default=s',	# In case of ANNOVAR code in annotation
	'AlleleFrequency!',			# Allele Frequency Calculation
	# output
	'output|output_file=s',	# Output file
	# show
	'show_annotations!',		# Show annotation
);    



our %parameters = ( #
	# Main options
	'help'		=>  	0,	# Help parameter
	'man'		=>  	0,	# Man parameter
	'release'	=>  	0,	# Release parameter
	'debug'		=>  	0,	# Debug parameter
	'verbose'	=>  	0,	# Verbose parameter
	# Configuration
	'config'		=>	'config.ini',			# Configuration file
	'config_annotation'    	=>  	"config.annotation.ini",	# Configuration annotation file
	'config_filter'    	=>  	"config.prioritization.ini",		# Configuration filter file
	# Input
	'input'		=>	undef,		# Input VCF file
	'filter'	=>	'ALL',	# Filter sources
	'hard'		=>	0,		# Hard filterind parameters
	# Output
	'output'	=>	undef,		# Output VCF file
);

## Parameters definition

our @options=(
	# Main options
	'help|h|?',		# Help
	'man',			# Man
	'release',		# Release
	'debug',		# Debug
	'verbose',		# Verbose
	# Configuration
	'config|config_file=s',				# Configuration file
	'config_annotation|config_annotation_file=s',	# Configuration annotation file
	'config_filter|config_filter_file|config_prioritization|config_prioritization_file=s',		# Configuration filter file
	# Input
	'input|input_file=s',	# Input file
	'filter=s',		# Filter
	'hard!',		# Hard Filtering
	# output
	'output|output_file=s',	# Output file
);






