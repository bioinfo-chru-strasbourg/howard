HOWARD
============

HOWARD annotates and prioritizes variants, calculates and normalizes annotations, translates vcf format and generates variants statistics.
HOWARD annotation is mainly based on ANNOVAR and snpEff tools to annotate, using available databases (see ANNOVAR and snpEff) and home made databases. It also uses BCFTOOLS to annotate variants with a VCF file. ANNOVAR and snpEff databases are automatically downloaded if needed.
HOWARD calculation harmonizes allele frequency (VAF), extracts Nomen (transcript, cNomen, pNomen...) from HGVS fields with an optional list of personalized transcripts, generates VaRank format barcode.
HOWARD prioritization algorithm uses profiles to flag variants (as passed or filtered), calculate a prioritization score, and automatically generate a comment for each variants (example: 'polymorphism identified in dbSNP. associated to Lung Cancer. Found in ClinVar database').Prioritization profiles are defined in a configuration file. A profile is defined as a list of annotation/value, using wildcards and comparison options (contains, lower than, greater than, equal...). Annotations fields may be quality values (usually from callers, such as 'GQ', 'DP') or other annotations fields provided by annotations tools, such as HOWARD itself (example: COSMIC, Clinvar, 1000genomes, PolyPhen, SIFT). Multiple profiles can be used simultaneously, which is useful to define multiple validation/prioritization levels (example: 'standard', 'stringent', 'rare variants', 'low allele frequency').
HOWARD translates VCF format into TSV format, by sorting variants using specific fields (example : 'prioritization score', 'allele frequency', 'gene symbol'), including/excluding annotations/fields, including/excluding variants, adding fixed columns.
HOWARD generates statistics files with a specific algorithm, snpEff and BCFTOOLS.
HOWARD is multithreaded through the number of variants and by database (data-scaling).


Getting Started
---------------

The default build of this image presents a container that runs HOWARD App on CentOS.

### Image Layout

Includes yum modules and other tools dependencies.


Building
--------

The `Dockerfile` provided with this package provides everything that is
needed to build the image. The build system must have Docker installed in
order to build the image.

```
$ cd PROJECT_ROOT
$ docker build -t howard:latest .
```
> Note: PROJECT_ROOT should be replaced with the path to where you have
>       cloned this project on the build system


Running a Container
-------------------

The container host must have Docker installed in order to run the image as a
container. Then the image can be pulled and a container can be started directly.

```
$ docker run howard:latest
```

### Swtiches

Any standard Docker switches may be provided on the command line when running
a container. Some specific switches of interest are documented below.

#### Configuration
```
-v HOST_DATA_FOLDER:/data
-v HOST_DATABASES_FOLDER:/databases
```
Content may be copied directly into the running container using a
`docker cp ...` command, alternatively one may choose to simply expose a host
configuration folder to the container.

### Examples

Run HOWARD as a uniq command.

```
$ docker run --rm howard:latest HOWARD --input=sample.vcf --output=sample.tsv --annotation=symbol,hgvs --calculation=VAF,NOMEN --prioritization=default --translation=TSV
```

Start HOWARD container.

```
$ docker run --name howard --entrypoint=bash -ti howard:latest
```


Debugging
---------

You may connect to a running container using the following command
```
$ docker exec -it --user root CONTAINER_NAME /bin/bash
```
> Note: CONTAINER_NAME should be the name provided to Docker when creating the
>       container initially. If not provided explicitly, Docker may have
>       assigned a random name. A container ID may also be used.

You may tail the container logs using the following commands
```
$ docker logs -f CONTAINER_NAME
```
