---
title: HOWARD Help Configuration
---

- [<span class="toc-section-number">1</span>
  Introduction](#introduction)
- [<span class="toc-section-number">2</span> folders](#folders)
  - [<span class="toc-section-number">2.1</span> databases](#databases)
    - [<span class="toc-section-number">2.1.1</span> genomes](#genomes)
    - [<span class="toc-section-number">2.1.2</span>
      annotations](#annotations)
    - [<span class="toc-section-number">2.1.3</span> parquet](#parquet)
    - [<span class="toc-section-number">2.1.4</span>
      bcftools](#bcftools)
    - [<span class="toc-section-number">2.1.5</span> annovar](#annovar)
    - [<span class="toc-section-number">2.1.6</span> snpeff](#snpeff)
    - [<span class="toc-section-number">2.1.7</span>
      exomiser](#exomiser)
    - [<span class="toc-section-number">2.1.8</span> refseq](#refseq)
- [<span class="toc-section-number">3</span> tools](#tools)
  - [<span class="toc-section-number">3.1</span> bcftools](#bcftools-1)
  - [<span class="toc-section-number">3.2</span> bgzip](#bgzip)
  - [<span class="toc-section-number">3.3</span> java](#java)
  - [<span class="toc-section-number">3.4</span> snpeff](#snpeff-1)
  - [<span class="toc-section-number">3.5</span> annovar](#annovar-1)
  - [<span class="toc-section-number">3.6</span> exomiser](#exomiser-1)
  - [<span class="toc-section-number">3.7</span> splice](#splice)
- [<span class="toc-section-number">4</span> chunking](#chunking)
  - [<span class="toc-section-number">4.1</span>
    chunking_enable](#chunking_enable)
  - [<span class="toc-section-number">4.2</span>
    chunking_size](#chunking_size)
  - [<span class="toc-section-number">4.3</span>
    chunking_partitions](#chunking_partitions)
  - [<span class="toc-section-number">4.4</span>
    chunking_sort](#chunking_sort)
- [<span class="toc-section-number">5</span> threads](#threads)
- [<span class="toc-section-number">6</span> memory](#memory)
- [<span class="toc-section-number">7</span> assembly](#assembly)
- [<span class="toc-section-number">8</span> verbosity](#verbosity)
- [<span class="toc-section-number">9</span> tmp](#tmp)
- [<span class="toc-section-number">10</span> access](#access)
- [<span class="toc-section-number">11</span>
  duckdb_settings](#duckdb_settings)
- [<span class="toc-section-number">12</span> chunk_size](#chunk_size)
- [<span class="toc-section-number">13</span> log](#log)

# Introduction

HOWARD Configuration JSON file defined default configuration regarding
resources (e.g. threads, memory), settings (e.g. verbosity, temporary
files), default folders (e.g. for databases) and paths to external
tools.

Examples:

> Example of a configuration JSON file

> ``` json
> {
>   "threads": 8,
>   "memory": null,
>   "verbosity": "WARNING",
>   "folders": {
>     "databases": {
>       "genomes": "~/howard/databases/genomes/current",
>       "annotations": [
>         "~/howard/databases/annotations/current",
>         "~/howard/databases/dbnsfp/current",
>         "~/howard/databases/dbsnp/current"
>       ],
>       "parquet": [
>         "~/howard/databases/annotations/current"
>       ],
>       "bcftools": [
>         "~/howard/databases/annotations/current"
>       ],
>       "annovar": "~/howard/databases/annovar/current",
>       "snpeff": "~/howard/databases/snpeff/current",
>       "varank": "~/howard/databases/varank/current"
>     }
>   },
>   "tools": {
>     "bcftools": "bcftools",
>     "bgzip": "bgzip",
>     "java": "java",
>     "snpeff": "~/howard/tools/snpeff/current/bin/snpEff.jar",
>     "annovar": "~/howard/tools/annovar/current/bin/table_annovar.pl",
>     "exomiser": "~/howard/tools/exomiser/current/bin/exomiser-cli-13.2.0.jar",
>     "splice": {
>        "docker": {
>          "image": "bioinfochrustrasbourg/splice:0.2.4",
>          "entrypoint": "/bin/bash",
>          "options": null,
>          "command": null
>        }
>     }
>   }
> }
> ```

# folders

Folders configuration such as for databases.

## databases

Default folders for databases that follows the specific database HOWARD
format. These folders will be used in HOWARD tools to autodetect
databases by their name and using assembly. Within database folders,
multiple releases can be provides (e.g. 'b152' and 'b156' for dbSNP and
'hg19' assembly within folder '~/howard/databases/dbsnp/current/hg19',
in 2 subfolders resp.)

Format:
`/path/to/databases/<db_name>/<db_release>/<assembly>/<database_files>`

Examples:

> Example of a configuration for databases folders

> ``` json
> {
>    "databases": {
>       "genomes": "~/howard/databases/genomes/current",
>       "annotations": [
>          "~/howard/databases/annotations/current",
>          "~/howard/databases/dbnsfp/current",
>          "~/howard/databases/dbsnp/current"
>       ],
>       "parquet": [
>          "~/howard/databases/annotations/current"
>       ],
>       "bcftools": [
>          "~/howard/databases/bcftools/current"
>       ],
>       "annovar": "~/howard/databases/annovar/current",
>       "snpeff": "~/howard/databases/snpeff/current",
>       "exomiser": "~/howard/databases/exomiser/current"
>    }
> }
> ```

### genomes

Genome folder with, for each assembly, FASTA files, indexes, and all
files generated by pygenome module.

Type: `Path`

Format: `A folder path (without assembly)`

Default: `~/howard/databases/genomes/current`

Examples:

> Path to genomes folder

> ``` json
> {
>    "genomes": "~/howard/databases/genomes/current"
> }
> ```

### annotations

Annotation databases folders that contains databases in various format
such as Parquet, VCF, duckDB and TSV.

Type: `Path`

Format: `a list of folder path (without assembly)`

Default: `["~/howard/databases/annotations/current"]`

Examples:

> Uniq folder with multiple annotations for Parquet annotation method,
> or other External tools

> ``` json
> {
>    "annotations": [
>       "~/howard/databases/annotations/current"
>    ]
> }
> ```

> Combinason of 2 folders with multiple annotations for Parquet
> annotation method, or other External tools

> ``` json
> {
>    "annotations": [
>       "~/howard/databases/annotations/current",
>       "~/howard/databases/dejavu/current",
>       "~/howard/databases/dbnsfp/current"
>    ]
> }
> ```

### parquet

Annotation databases folders that contains databases in Parquet format.

Format: `a list of folder path (without assembly)`

Default: `["~/howard/databases/annotations/current"]`

Examples:

> Uniq folder with multiple annotations for Parquet annotation method

> ``` json
> {
>    "annotations": [
>       "~/howard/databases/annotations/current"
>    ]
> }
> ```

> Combinason of 2 folders with multiple annotations for Parquet
> annotation method

> ``` json
> {
>    "annotations": [
>       "~/howard/databases/annotations/current",
>       "~/howard/databases/dejavu/current",
>       "~/howard/databases/dbnsfp/current"
>    ]
> }
> ```

### bcftools

Annotation databases folders for BCFTools annotation.

Format: `a list of folder path (without assembly)`

Default: `["~/howard/databases/bcftools/current"]`

Examples:

> Uniq folder with multiple VCF and BED files for BCFTools annotation

> ``` json
> {
>    "bcftools": [
>       "~/howard/databases/bcftools/current"
>    ]
> }
> ```

> Combinason of 2 folders with multiple VCF and BED files for BCFTools
> annotation

> ``` json
> {
>    "bcftools": [
>       "~/howard/databases/bcftools/current",
>       "~/howard/databases/dejavu/current"
>    ]
> }
> ```

### annovar

Annotation databases folder for Annovar annotation.

Format: `a list of folder path (without assembly)`

Default: `["~/howard/databases/annovar/current"]`

Examples:

> Uniq folder with multiple Annovar TXT files for Annovar annotation

> ``` json
> {
>    "annovar": "~/howard/databases/annovar/current/"
> }
> ```

### snpeff

Annotation databases folders for snpEff annotation.

Format: `A folder path (without assembly)`

Default: `~/howard/databases/snpeff/current`

Examples:

> Path to snpEff database folder

> ``` json
> {
>    "snpeff": "~/howard/databases/snpeff/current/"
> }
> ```

### exomiser

Annotation databases folders for Exomiser annotation.

Format: `A folder path (without assembly)`

Default: `~/howard/databases/exomiser/current`

Examples:

> Path to Exomiser database folder

> ``` json
> {
>    "exomiser": "~/howard/databases/exomiser/current/"
> }
> ```

### refseq

Annotation databases folders for refSeq annotation.

Type: `str`

Format: `A folder path (without assembly)`

Default: `~/howard/databases/refseq/current`

Examples:

> Path to refSeq files folder

> ``` json
> {
>    "refseq": "~/howard/databases/refseq/current/"
> }
> ```

# tools

External tools paths that can be defined as path to a binary or a dict
including the binary type (such as "bin", "jar", "perl"). External tools
can be configured with docker, using 'docker' as binary type and options
to define docker 'image' (mandatory), to specify 'entrypoint', 'command'
and docker 'options' (e.g. folder mount '-v
/path/to/folder:/path/to/folder').

Examples:

> Example of a configuration for tools, with env \$PATH, full path and
> path with type

> ``` json
> {
>    "tools": {
>       "bcftools": "bcftools",
>       "bgzip": "bgzip",
>       "java": "/usr/bin/java",
>       "snpeff": "~/howard/tools/snpeff/current/bin/snpEff.jar",
>       "annovar": {"jar": "~/howard/tools/annovar/current/bin/table_annovar.pl"},
>       "exomiser": {"jar": "~/howard/tools/exomiser/current/bin/exomiser-cli-13.2.0.jar"}
>    }
> }
> ```

> Example of a configuration for bcftools with a docker image (example
> with howard docker image)

> ``` json
> {
>    "tools": {
>       "bcftools": {
>          "docker": {
>            "image": "howard:0.13.0",
>            "entrypoint": "bcftools",
>            "options": null,
>            "command": null
>          }
>       }
>    }
> }
> ```

> Example of a configuration for splice with a docker image

> ``` json
> {
>    "tools": {
>       "splice": {
>          "docker": {
>            "image": "bioinfochrustrasbourg/splice:0.2.4",
>            "entrypoint": "/bin/bash",
>            "options": null,
>            "command": null
>          }
>       }
>    }
> }
> ```

## bcftools

BCFTools binary (see https://samtools.github.io/bcftools/).

Default: `bcftools`

Examples:

> Path to binary in \$PATH env variable

> ``` json
> {
>    "bcftools": "bcftools"
> }
> ```

> Path to binary as a dict with binary type 'bin'

> ``` json
> {
>    "bcftools": {"bin": "~/howard/tools/bcftools/current/bin/bcftools"}
> }
> ```

## bgzip

BGZip binary (see https://samtools.github.io/bcftools/).

Default: `bgzip`

Examples:

> Path to binary in \$PATH env variable

> ``` json
> {
>    "bgzip": "bgzip"
> }
> ```

> Path to binary as a dict with binary type 'bin'

> ``` json
> {
>    "bgzip": {"bin": "~/howard/tools/htslib/current/bin/bgzip"}
> }
> ```

## java

Java binary (see https://www.java.com).

Default: `java`

Examples:

> Path to binary in \$PATH env variable

> ``` json
> {
>    "java": "java"
> }
> ```

> Path to binary as a dict with binary type 'bin'

> ``` json
> {
>    "java": {"bin": "/usr/bin/java"}
> }
> ```

## snpeff

snpEff binary (see https://pcingola.github.io/SnpEff/).

Default: `~/howard/tools/snpeff/current/bin/snpEff.jar`

Examples:

> Path to binary as a dict without binary type

> ``` json
> {
>    "snpeff": "~/howard/tools/snpeff/current/bin/snpEff.jar"
> }
> ```

> Path to binary as a dict with binary type 'jar'

> ``` json
> {
>    "snpeff": {"jar": "~/howard/tools/snpeff/current/bin/snpEff.jar"}
> }
> ```

## annovar

ANNOVAR perl script (see https://annovar.openbioinformatics.org/).

Default: `~/howard/tools/annovar/current/bin/table_annovar.pl`

Examples:

> Path to binary as a dict without binary type

> ``` json
> {
>    "annovar": "~/howard/tools/annovar/current/bin/table_annovar.pl"
> }
> ```

> Path to binary as a dict with binary type 'perl'

> ``` json
> {
>    "annovar": {"jar": "~/howard/tools/annovar/current/bin/table_annovar.pl"}
> }
> ```

## exomiser

Exomiser binary (see https://www.sanger.ac.uk/tool/exomiser/).

Default: `~/howard/tools/exomiser/current/bin/exomiser-cli-13.2.0.jar`

Examples:

> Path to binary as a dict without binary type

> ``` json
> {
>    "snpeff": "~/howard/tools/exomiser/current/bin/exomiser-cli-13.2.0.jar"
> }
> ```

> Path to binary as a dict with binary type 'jar'

> ``` json
> {
>    "snpeff": {"jar": "~/howard/tools/exomiser/current/bin/exomiser-cli-13.2.0.jar"}
> }
> ```

## splice

Splice Docker image binary (see
https://hub.docker.com/r/bioinfochrustrasbourg/splice).

Default: `None`

Examples:

> Configuration of Docker image

> ``` json
> {
>    "splice": {
>       "docker": {
>         "image": "bioinfochrustrasbourg/splice:0.2.4",
>         "entrypoint": "/bin/bash",
>         "options": null,
>         "command": null
>       }
>    }
> }
> ```

# chunking

Chunking is a technique used to process large files by dividing them
into smaller, more manageable pieces (chunks), by specifying chunk size
and partitioning scheme. This approach allows each chunk to be processed
independently, reducing memory usage, disk usage (duckDB swap), and can
help prevent crashes or slowdowns due to memory constraints. Once all
chunks are processed, the results are merged to produce the final
output.

This method is particularly useful for handling large datasets or files
that would otherwise be too resource-intensive to process in one go
(such as low memory, or data particularly huge due to extensive
annotation).

However, splitting, processing, and merging chunks can introduce
additional computational overhead, may mislead with some calculations
(due to not covering the full region of interest, such as all variants
of a gene for compound heterozygote calculation), and can not apply full
data features (such as statistics). Moreover, the order of data can
change during the chunking-merging process.

To enable chunking, include the `chunking` section in the parameters
JSON file. This section includes parameters to enable chunking, define
chunk size, partitioning scheme, and sorting of the final output.

Examples:

> Default configuration for chunking

> ``` json
> {
>    "chunking": {
>       "chunking_enable": false,
>       "chunking_size": 1000000
>       "chunking_partitions": "None",
>       "chunking_sort": false
> }
> ```

> Example of a configuration for chunking big files by chromosome with
> 100 million variants per chunk and sorting of final output

> ``` json
> {
>    "chunking": {
>       "chunking_enable": true,
>       "chunking_size": 100000000
>       "chunking_partitions": "#CHROM",
>       "chunking_sort": false
> }
> ```

## chunking_enable

Chunking process by splitting input file into smaller parts. Use
chunking parameters to define chunking size and partitioning fields.

Default: `False`

Examples:

> Enable chunking

> ``` json
> {
>    "chunking_enable": true
> }
> ```

> Disable chunking

> ``` json
> {
>    "chunking_enable": false
> }
> ```

## chunking_size

Chunking size in number of variants to use when chunking is enabled.
Chunking will be applied if number of variants exceeding this chunking
size. For example, setting this to 1,000,000 will chunk input files with
more than 1 million variants into multiple files, each containing up to
1 million variants.

Type: `int`

Default: `1000000`

Examples:

> Set chunking size to 1 million variants

> ``` json
> {
>    "chunking_size": 1000000
> }
> ```

> Set chunking size to 10 million variants

> ``` json
> {
>    "chunking_size": 10000000
> }
> ```

## chunking_partitions

Partitioning scheme to use when chunking is enabled, defining how the
input file is partitioned during chunking. For example, setting this to
'#CHROM' will partition the input file by chromosome.

Type: `str`

Default: `None`

Examples:

> No partitioning

> ``` json
> {
>    "chunking_partitions": "None"
> }
> ```

> Partition by chromosome

> ``` json
> {
>    "chunking_partitions": "#CHROM"
> }
> ```

## chunking_sort

Sort final output after chunking is applied. Prevents issues with
unordered variants after chunked processing.

Default: `False`

Examples:

> Enable sorting

> ``` json
> {
>    "chunking_sort": true
> }
> ```

> Disable sorting

> ``` json
> {
>    "chunking_sort": false
> }
> ```

# threads

Specify the number of threads to use for processing HOWARD. It
determines the level of parallelism, either on python scripts, duckdb
engine and external tools. It and can help speed up the process/tool.
Use -1 to use all available CPU/cores. Either non valid value is 1
CPU/core.

Type: `int`

Default: `-1`

Examples:

> Automatically detect all available CPU/cores

> ``` json
> {
>    "threads": -1
> }
> ```

> Define 8 CPU/cores

> ``` json
> {
>    "threads": 8
> }
> ```

# memory

Specify the memory to use in format FLOAT\[kMG\] (e.g. '8G', '12.42G',
'1024M'). It determines the amount of memory for duckDB engine and
external tools (especially for JAR programs). It can help to prevent
'out of memory' failures. By default (None) is 80%% of RAM (for duckDB).

Type: `str`

Format: `FLOAT[kMG]`

Default: `None`

Examples:

> Automatically detect all available CPU/cores

> ``` json
> {
>    "threads": -1
> }
> ```

> Define 8 CPU/cores

> ``` json
> {
>    "threads": 8
> }
> ```

# assembly

Genome Assembly (e.g. 'hg19', 'hg38').

Type: `str`

Default: `hg19`

Examples:

> Default assembly for all analysis tools

> ``` json
> {
>    "assembly": "hg19" 
> }
> ```

> List of assemblies for databases download tool

> ``` json
> {
>    "assembly": "hg19,hg38" 
> }
> ```

# verbosity

Verbosity level Available: CRITICAL, ERROR, WARNING, INFO, DEBUG or
NOTSET

- DEBUG: Detailed information, typically of interest only when
  diagnosing problems.

- INFO: Confirmation that things are working as expected.

- WARNING: An indication that something unexpected happened.

- ERROR: Due to a more serious problem.

- CRITICAL: A serious error.

- FATAL: A fatal error.

- NOTSET: All messages.

Type: `str`

Choices:
`['CRITICAL', 'ERROR', 'WARNING', 'INFO', 'DEBUG', 'NOTSET', 'WARN', 'FATAL']`

Default: `INFO`

Examples:

> Default verbosity

> ``` json
> {
>    "verbosity": "INFO"
> }
> ```

> ERROR level (quiet mode)

> ``` json
> {
>    "verbosity": "ERROR"
> }
> ```

> For debug

> ``` json
> {
>    "verbosity": "DEBUG"
> }
> ```

# tmp

Temporary folder (e.g. '/tmp'). By default, '.tmp' for duckDB (see
doc),external tools and python scripts.

Type: `Path`

Default: `None`

Examples:

> System temporary folder

> ``` json
> {
>    "tmp": "/tmp"
> }
> ```

> HOWARD work directory

> ``` json
> {
>    "tmp": "~/howard/tmp"
> }
> ```

> Current work directory

> ``` json
> {
>    "tmp": ".tmp"
> }
> ```

# access

Access mode to variants file or database. Either 'RW' for Read and
Write, or 'RO' for Read Only.

Type: `str`

Choices: `['RW', 'RO']`

Default: `RW`

Examples:

> Read and Write mode

> ``` json
> {
>    "access": "RW"
> }
> ```

> Read only mode

> ``` json
> {
>    "access": "RO"
> }
> ```

# duckdb_settings

DuckDB settings (see duckDB doc) as JSON (string or file). These
settings have priority (see options 'threads', 'tmp'...). Examples:
'{"TimeZone": "GMT", "temp_directory": "/tmp/duckdb", "threads": 8}'.

Type: `Path`

Default: `None`

Examples:

> DuckDB settings JSON file

> ``` json
> {
>    "duckdb_settings": "/path/to/duckdb_config.json"
> }
> ```

> JSON string for Time zone, temporary directory and threads for duckDB

> ``` json
> {
>    "duckdb_settings": {
>       "TimeZone": "GMT",
>       "temp_directory": "/tmp/duckdb",
>       "threads": 8
>    }
> }
> ```

# chunk_size

Number of records in batch to export output file. The lower the chunk
size, the less memory consumption. For Parquet partitioning, files size
will depend on the chunk size.

Type: `int`

Default: `1048576`

Examples:

> Chunk size of 1048576 by default

> ``` json
> {
>    "chunk_size": 1048576
> }
> ```

> Smaller chunk size to reduce Parquet file size and memory usage

> ``` json
> {
>    "chunk_size": 100000
> }
> ```

# log

Logs file (e.g. 'my.log').

Type: `Path`

Default: `None`

Examples:

> Relative path to log file

> ``` json
> {
>    "log": "my.log"
> }
> ```

> HOWARD work directory

> ``` json
> {
>    "log": "~/howard/log"
> }
> ```

> Full path to log file

> ``` json
> {
>    "log": "/tmp/my.log"
> }
> ```
