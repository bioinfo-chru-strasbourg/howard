# Databases
<!-- TOC -->

- [Databases](#databases)
- [UCSC](#ucsc)
  - [Clinvar](#clinvar)
    - [Resume](#resume)
    - [Update](#update)
    - [Format](#format)
    - [Assembly](#assembly)
  - [PhyloP](#phylop)
    - [Resume](#resume-1)
    - [Update](#update-1)
    - [Format](#format-1)
    - [Assembly](#assembly-1)
  - [PhastCons](#phastcons)
    - [Resume](#resume-2)
    - [Update](#update-2)
    - [Format](#format-2)
    - [Assembly](#assembly-2)
  - [GERP++](#gerp)
    - [Resume](#resume-3)
    - [Update](#update-3)
    - [Format](#format-3)
    - [Assembly](#assembly-3)
- [Annovar *TODO conf Antony*](#annovar-todo-conf-antony)
- [Autres](#autres)
  - [Alphamissence](#alphamissence)
    - [Resume](#resume-4)
    - [Update](#update-4)
    - [Format](#format-4)
    - [Assembly](#assembly-4)
  - [1000g](#1000g)
    - [Update](#update-5)
    - [Format](#format-5)
    - [Assembly](#assembly-5)
  - [gnomAD](#gnomad)
    - [Resume](#resume-5)
    - [Update](#update-6)
    - [Format](#format-6)
    - [Assembly](#assembly-6)

<!-- /TOC -->
<!-- /TOC -->
<!-- /TOC -->
<!-- /TOC -->
# UCSC
FTP: https://ftp.ncbi.nlm.nih.gov

## Clinvar

### Resume
Base de données de variants rapportés par différents laboratoire à travers le monde dont le labo de diag de Strasbourg et curé par des experts. Contient des variants majoritairement associé àune patho mais aussi des Classes 1 à 3

### Update
Récupérer de l'UCSC, mise à jour hebdomadaire sur le FTP, https://ftp.ncbi.nlm.nih.gov/pub/clinvar
Devrait rester fixe dans le temps

### Format
Format <b>VCF</b>, à transformer en parquet, base de données de quelques Mo donc parquet unique 

### Assembly
version 37 et 38 disponible et mis à jour en même temps

## PhyloP

### Resume
Score de Conservation nucléotidique, issues d'alignements multiple inter-espèces scores positifs sites conservés au cours de l'évolution et négatifs positions qui évoluent rapidement, peu conservés (threshold iontorrent min -20 et positif +30), le way correspond au nombre d'espèce utilisé pour l'alignement

### Update
Pas d'update fichiers d'il y a 10 ans pour hg19, hg38 plus récent
https://hgdownload.cse.ucsc.edu/goldenPath/

### Format
Format <b>BigWig</b> de l'UCSC, à transformer en bed puis en parquet, script dans le plugin de mise à jour des bases de données.

### Assembly
version dossier hg19 et hg38 disponible

## PhastCons

### Resume
Comme PhyloP score de conservation nucléotidique basé sur un HMM pour établier la proba qu'une base soit conservé ou non

### Update
Pas d'update fichiers d'il y a 10 ans pour hg19, hg38 plus récent

### Format
Format <b>BigWig</b> de l'UCSC, à transformer en bed puis en parquet, script dans le plugin de mise à jour des bases de données.

### Assembly
version dossier hg19 et hg38 disponible
## GERP++
### Resume
GERP score est une mesure de la conservation inter espèce. Difference entre le nombre attendu et observé de substitution a un locus entre plusieurs espèces.

### Update
Compliqué déjà à trouver
peut etre en 38 des updates
hg38: https://ftp.ensembl.org/pub/current_compara/conservation_scores/91_mammals.gerp_conservation_score/gerp_conservation_scores.homo_sapiens.GRCh38.bw
hg19:https://hgdownload.soe.ucsc.edu/gbdb/hg19/bbi/All_hg19_RS.bw
### Format
Format <b>BigWig</b> de l'UCSC, à transformer en bed puis en parquet, script dans le plugin de mise à jour des bases de données.
### Assembly
version 37 et 38 disponible UCSC et ensembl


# Annovar *TODO conf Antony*

# Autres

## Alphamissence

### Resume
Machine Learning model de DeepMind pour prédire la pathogénicité des variants faux sens. Plusieurs annotations disponible, tout les missense possible sur les transcripts canoniques, la moyenne des scores de ceux-ci par gènes. Des prédictions par Acide aminés en fonction des id uniprot et prédiction hg38 uniquement pour des transcripts non canonique input provenant de genecode

### Update
Version de fin 2023
https://zenodo.org/records/8208688

### Format
Format <b>tsv.gz</b>, étape de processing pour les transformer en parquet, mono block ou hive à voir

### Assembly
version dossier hg19 et hg38 disponible

## 1000g
Projet de séquençage de personnes sain, permet de déterminer une fréquence de population de variants

### Update
Ca va pas bouger, j'ai repris de Nirvanna, GRCh"7 de 2013 et 38 de 2019

### Format
VCF donc a transformer en parquet ou laisser tel quel

### Assembly
version dossier hg19 et hg38 disponible
hg19: https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/
hg38: https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20190312_biallelic_SNV_and_INDEL/

## gnomAD
### Resume
<b>gs://gcp-public-data--gnomad/release/</b>
Base de données comprenant le séquençage de milliers de génomes et d'exomes
toutes les données sont disponible sur le cloud google, manière la plus simple et rapide de tout télécharger

### Update
A prévoir peut être une nouvelle version après la 4.1.0, sinon version stable donc pas besoin de prévoir une politique d'update
Certaines valeurs de pLi seront sans doute à calculer d'après les annotations des fichiers constraints

### Format
VCF pour les variants ou TSV pour les constraints
Prévoir un formatage en parquet pour les VCFs et from extann pour les TSV avec descriptif des colonnes
pLi: gnomad_pli
LOEUF: oe_lof_upper
o/e: oe_lof


### Assembly

Version 2.1.1 en GRCh37 avec une version liftover en GRCh38
Version 4.1.0 en GRCh38 only (version Beta d'après gnomAD)


HUSDIAGGEN

dejavu_HUSDIAGGEN_GOMV1_SOMATIC
popfreq
gnomAD
gene_name
transcript_name
c_nomen
PNOMEN
location
COVAR_ANN
COVAR_COM
CADD_phred
DEEP_INTRON
CLINVAR
weighted_variation_class_name
weighted_variation_class_number
chromosome_number
c_dna_start
REF
ALT
variant_name
coding_effect
PZFlag
ENOMEN
ANN
dbscsnv
dbSNP
dbSNPNonFlagged
EXAC
gme
n1000genomesAFR
n1000genomesALL
n1000genomesAMR
n1000genomesEUR
n1000genomesSAS
n6500NHLBIAA
n6500NHLBIALL
n6500NHLBIEA
FindByPipelines
id


TUMSOL
 
    - chr
    - pos
    - ref
    - alt
    - nomen
    - findbypipelines
    - vaf_list
    - pzscore-tumsol
    - count_var
    - dejavu.hustumsol.xths
    - dbsnpnonflagged
    - popfreq
    - gnomad
    - cosmic
    - interpro_domain
    - clinvar
    - mutationtaster_pred
    - sift_pred
    - outcome
    - samples.T22395.gt

DIAG

    - varankvarscore
    - gene
    - cnomen
    - pnomen
    - transcript
    - chr
    - pos
    - ref
    - alt
    - tags
    - genedesc
    - clinvarclinsignifs
    - clinvarphenotypes
    - rsid
    - gnomadaltfreq_all
    - gnomadaltfreq_popmax
    - gnomadhomcount_all
    - gnomadhetcount_all
    - gnomadhemcount_all
    - deltamaxentscorepercent
    - deltassfscorepercent
    - deltannsscorepercent
    - localspliceeffect
    - siftprediction
    - homcount
    - hetcount
    - allelecount
    - samplecount
    - allelefrequency
    - ddd_hi_percent
    - omim_id
    - omim_phenotype
    - omim_inheritance
    - gnomen
    - vartype
    - codingeffect
    - varlocation
    - wtmaxentscore
    - varmaxentscore
    - wtssfscore
    - varssfscore
    - wtnnsscore
    - varnnsscore
    - nearestsschange
    - distnearestss
    - nearestsstype
    - branchpointpos
    - branchpointchange
    - phylop
    - phastcons
    - granthamdist
    - samva
    - loeuf_bin
    - omim_morbid
    - exomiser_gene_pheno_score
    - exon
    - intron
    - gnomad_pli
    - gnomadaltfreq_afr
    - gnomadaltfreq_amr
    - gnomadaltfreq_asj
    - gnomadaltfreq_eas
    - gnomadaltfreq_fin
    - gnomadaltfreq_nfe
    - gnomadaltfreq_oth
    - gnomadaltfreq_sas
    - 1000g_af