# HOWARD Configuration

HOWARD Configuration JSON file defined default configuration regarding resources (e.g. threads, memory), settings (e.g. verbosity, temporary files), default folders (e.g. for databases) and paths to external tools.

## Table of contents

- [HOWARD Configuration](#howard-configuration)
   - [threads](#threads)
   - [memory](#memory)
   - [assembly](#assembly)
   - [verbosity](#verbosity)
   - [tmp](#tmp)
   - [access](#access)
   - [duckdb_settings](#duckdb_settings)
   - [chunk_size](#chunk_size)
   - [log](#log)
   - [folders](#folders)
      - [databases](#foldersdatabases)
         - [annotations](#foldersdatabasesannotations)
         - [parquet](#foldersdatabasesparquet)
         - [bcftools](#foldersdatabasesbcftools)
         - [annovar](#foldersdatabasesannovar)
         - [snpeff](#foldersdatabasessnpeff)
         - [exomiser](#foldersdatabasesexomiser)
   - [tools](#tools)
      - [bcftools](#toolsbcftools)
      - [bgzip](#toolsbgzip)
      - [snpeff](#toolssnpeff)
      - [annovar](#toolsannovar)
      - [exomiser](#toolsexomiser)
      - [java](#toolsjava)


## threads

Number of threads to use for processing HOWARD. It determines the level of parallelism, either on python scripts, duckDB engine and external tools. It and can help speed up the process/tool

Use -1 to use all available CPU/cores.

Default: -1

Examples: 
- -1 (for all available CPU/cores)

- 8 (for 8 CPU/cores)

## memory

Specify the memory to use. It determines the amount of memory for duckDB engine and external tools (especially for JAR prorams). It can help to prevent 'out of memory' failures.

Format: FLOAT[kMG]

Default: None (80% of RAM for duckDB)

Examples: 
- 8G (for 8Go of memory)

- 1024M (for 1Go of memory)

- 24.8G (for 24Go and 800Mo of memory)

## assembly

Default assembly. This parameter will by overwritten in the paramter JSON file.

Default: hg19

Examples: 
- "hg19" (for Homo Sapiens hg19/GRCh37 assembly)

- "hg38" (for Homo Sapiens hg38/GRCh38 assembly)

- "ailMel1.99" (for Giant Panda assembly)

## verbosity

Verbosity level, such as:

- DEBUG: Detailed information, typically of interest only when diagnosing problems.

- INFO: Confirmation that things are working as expected.

- WARNING: An indication that something unexpected happened, or indicative of some problem in the near future (e.g. ‘disk space low’). The software is still working as expected.

- ERROR: Due to a more serious problem, the software has not been able to perform some function.

- CRITICAL: A serious error, indicating that the program itself may be unable to continue running.

- NOTSET: All messages.

Default: INFO

Examples: CRITICAL, ERROR, WARNING, INFO, DEBUG or NOTSET

## tmp

Temporary folder, espacially for duckDB (see doc), external tools and python scripts.

Format: Path

Default: None

Examples: 
- /tmp

- .tmp

## access

Access mode to variants file or database.Either 'RW' for Read and Write, or 'RO' for Read Only.

Format: either 'RO' or 'RW'

Default: RW

Examples: 
- RO

- RW

## duckdb_settings

DuckDB settings (see duckDB doc) as JSON (string or file). These settings have priority (see options 'threads', 'tmp'...).

Default: {}

Examples: 
- {"TimeZone": "GMT", "temp_directory": "/tmp/duckdb", "threads": 8}

- /path/to/duckdb_config.json

## chunk_size

Number of records in batch to export output file. The lower the chunk size, the less memory consumption. For Parquet partitioning, files size will depend on the chunk size.

Default: 1000000

Examples: 
- 1000000

- 100000

## log

Log file.

Default: None

Examples: 
- my.log

- /tmp/my.log

## folders

Configuration for folders such as databases

### folders::databases

Default folders for databases that follows the specific database HOWARD format. These folders will be used in HOWARD tools to autodetect databases by their name and using assembly.

Format: /path/to/databases/db_name/db_release/assembly/database_file

#### folders::databases::annotations

Annotation databases folders that contains databases in various format such as Parquet, VCF, duckDB and TSV.

Format: a list of folder path (without assembly)

Default: ["~/howard/databases/annotations/current"]

Examples: 
- ["~/howard/databases/annotations/current/"]

- ["~/howard/databases/annotations/current/","~/howard/databases/dejavu/current/","~/howard/databases/dbnsfp/current/"]

#### folders::databases::parquet

Annotation databases folders that contains databases in Parquet format.

Format: a list of folder path (without assembly)

Default: ["~/howard/databases/annotations/current"]

Examples: 
- ["~/howard/databases/parquet/current/"]

- ["~/howard/databases/annotations/current/"]

- ["~/howard/databases/parquet/current/","~/howard/databases/dejavu/current/"]

#### folders::databases::bcftools

Annotation databases folders for BCFTools annotation.

Format: a list of folder path (without assembly)

Default: ["~/howard/databases/bcftools/current"]

Examples: 
- ["~/howard/databases/bcftools/current/"]

- ["~/howard/databases/bcftools/current/","~/howard/databases/dejavu/current/"]

#### folders::databases::annovar

Annotation databases folders for Annovar annotation.

Format: a list of folder path (without assembly)

Default: ["~/howard/databases/annovar/current"]

Examples: 
- ["~/howard/databases/annovar/current/"]

- ["~/howard/databases/annovar/current/","~/howard/databases/annovar/homemade/"]

#### folders::databases::snpeff

Annotation databases folders for snpEff annotation.

Format: A folder path (without assembly)

Default: ~/howard/databases/snpeff/current

Examples: 
- "~/howard/databases/snpeff/current/"

#### folders::databases::exomiser

Annotation databases folders for Exomiser annotation.

Format: A folder path (without assembly)

Default: ~/howard/databases/exomiser/current

Examples: 
- "~/howard/databases/exomiser/current/"

## tools

External tools paths that can be defined as path to a binary or a dict including the binary type (such as "bin", "jar", "perl").

Examples: 
- "/path/to/tool/bin/tool.bin"

- {"bin": "/path/to_tool/bin/tool.sh"}

- {"bin": "/path/to_tool/bin/java"}, {"perl": "/path/to_tool/bin/tool.pl"}

### tools::bcftools

BCFTools binary (see https://samtools.github.io/bcftools/).

Default: ~/howard/tools/bcftools/current/bin/bcftools

Examples: 
- "~/howard/tools/bcftools/current/bin/bcftools"

- {"bin": "~/howard/tools/bcftools/current/bin/bcftools"}

### tools::bgzip

BGZip binary (see https://samtools.github.io/bcftools/).

Default: ~/howard/tools/htslib/current/bin/bgzip

Examples: 
- "~/howard/tools/htslib/current/bin/bgzip"

- {"bin": "~/howard/tools/htslib/current/bin/bgzip"}

### tools::snpeff

snpEff binary (see https://pcingola.github.io/SnpEff/).

Default: ~/howard/tools/snpeff/current/bin/snpEff.jar

Examples: 
- "~/howard/tools/snpeff/current/bin/snpEff.jar"

- {"jar": "~/howard/tools/snpeff/current/bin/snpEff.jar"}

### tools::annovar

ANNOVAR perl script (see https://annovar.openbioinformatics.org/).

Default: ~/howard/tools/annovar/current/bin/table_annovar.pl

Examples: 
- "~/howard/tools/annovar/current/bin/table_annovar.pl"

- {"perl": "~/howard/tools/annovar/current/bin/table_annovar.pl"}

### tools::exomiser

Exomiser binary (see https://www.sanger.ac.uk/tool/exomiser/).

Default: ~/howard/tools/exomiser/current/bin/exomiser-cli-13.2.0.jar

Examples: 
- "~/howard/tools/exomiser/current/bin/exomiser-cli-13.2.0.jar"

- {"jar": "~/howard/tools/exomiser/current/bin/exomiser-cli-13.2.0.jar"}

### tools::java

Java binary (see https://www.java.com).

Default: ~/howard/tools/java/current/bin/java

Examples: 
- "~/howard/tools/java/current/bin/java"

- "java"

- {"bin": "/usr/bin/java"}

