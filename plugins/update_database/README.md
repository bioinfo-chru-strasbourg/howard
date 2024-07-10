# Databases
<!-- TOC -->

- [Databases](#databases)
- [dbNSFP](#dbnsfp)
  - [InterproDomain](#interprodomain)
  - [SIFT / SIFT 4G](#sift--sift-4g)
  - [PROVEAN](#provean)
  - [fatHMM](#fathmm)
  - [PrimateAI](#primateai)
- [NCBI](#ncbi)
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
  - [dbSNP](#dbsnp)
    - [Resume](#resume-4)
    - [Update](#update-4)
    - [Format](#format-4)
    - [Assembly](#assembly-4)
- [Annovar *TODO conf Antony*](#annovar-todo-conf-antony)
  - [PopFreq](#popfreq)
    - [Resume](#resume-5)
    - [Update](#update-5)
    - [Assembly](#assembly-5)
    - [Format](#format-5)
  - [gme](#gme)
    - [Resume](#resume-6)
    - [Update](#update-6)
    - [Assembly](#assembly-6)
    - [Format](#format-6)
  - [outcome](#outcome)
    - [Resume](#resume-7)
    - [Update](#update-7)
- [Autres](#autres)
  - [Alphamissence](#alphamissence)
    - [Resume](#resume-8)
    - [Update](#update-8)
    - [Format](#format-7)
    - [Assembly](#assembly-7)
  - [1000g](#1000g)
    - [Update](#update-9)
    - [Format](#format-8)
    - [Assembly](#assembly-8)
  - [gnomAD](#gnomad)
    - [Resume](#resume-9)
    - [Update](#update-10)
    - [Format](#format-9)
    - [Assembly](#assembly-9)
  - [Decipher](#decipher)
    - [Resume](#resume-10)
    - [Update](#update-11)
    - [Format](#format-10)
    - [Assembly](#assembly-10)
  - [CADD](#cadd)
    - [Resume](#resume-11)
    - [Update](#update-12)
    - [Assembly](#assembly-11)
    - [Format](#format-11)

<!-- /TOC -->
<!-- /TOC -->
<!-- /TOC -->
<!-- /TOC -->
# dbNSFP

https://sites.google.com/site/jpopgen/dbNSFP
aggregator de differents outils pour évaluer la pathogeenicité des nsSNV sur le génome. La bdd ne contient en fait que les données d'exomes (hors outil WGSA).

Disponible en hg19 et 38, utiliser les versions académiques 4.6a par exemple.
Mise à jour régulière, plusieurs fois par an.
## InterproDomain
Base de données des familles de protéines, domaines et site importants. EMBL-EBI

## SIFT / SIFT 4G
Prediction de l'effet sur la protéine d'un changement d'AA. SIFT4G plus performant que SIFT.
https://sift.bii.a-star.edu.sg/sift4g/AboutSIFT4G.html#:~:text=SIFT%204G%20is%20a%20faster,predictions%20for%20single%20nucleotide%20variants.
score entre 0 et 1 ou label  Tolerated ou Damaging (T or D)

## PROVEAN
Prediction de l'effet sur la protéine d'un changement d'AA
score entre 0 et 1 ou label  Tolerated ou Neutral (T or N)

## fatHMM
Outil pour évaluer pathogénicité des faux sens avec algo HMM

## PrimateAI
Deep neural network pour évaluer la pathogénicité des mutations faux sesn

# NCBI
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

## dbSNP
### Resume
base de données de public de SNP depuis plus de 20 ans

### Update
Tout les ans / 2 ans a peu près

### Format
VCF mais comprenant les NC comme numéro de chromosome au lieu de chr, donc a process avant d'être utilisable chez nous, eventuellement à transformer en parquet

### Assembly
Release en 37 et 38 avec les noms des assemblages originaux
GCF_000001405.25 => GRCh37
GCF_000001405.40 => GRCh48 

# Annovar *TODO conf Antony*
## PopFreq
<b>DEPRECIATED (bases de données beaucoup plus récente, voir faire le score nous meme avec les nouvelles db)</b>

### Resume
Frequence Allélique max de de plusieures bases de données 1000g, Exac ESP6500

### Update
Aucune c'est 20150413 la version actuelle (modification de 2019 d'apres annovar avdblist)

### Assembly
que pour hg19

### Format
| #Chr | Start | End   | Ref | Alt | PopFreqMax |
|:-----|:------|:------|:----|:----|:-----------|
| 1    | 10177 | 10177 | A   | AC  | 0.49       |

## gme
### Resume
Greater middle east variome, base de donnée de fréquence allélique pour des populations peu représenté dans les grandes db, (afrique ouest et est , peninsule arabique israel...)

### Update
Comme popfreq on a une maj de 2019 on dirait

### Assembly
hg19 et 38

### Format


| #Chr | Start | End   | Ref | Alt | GME_AF   | GME_NWA  | GME_NEA  | GME_AP   | GME_Israel | GME_SD   | GME_TP   | GME_CA   |
|:-----|:------|:------|:----|:----|:---------|:---------|:---------|:---------|:-----------|:---------|:---------|:---------|
| 1    | 69134 | 69134 | A   | G   | 0.049505 | 0.000000 | 0.032787 | 0.000000 | 0          | 0.000000 | 0.181818 | 0.133333 |

## outcome
### Resume
Simplement la position du variant exon, splice UTR downstream etc établis par annovar en fonction du fichier refgene

### Update

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
Projet de séquençage de personnes "saines", permet de déterminer une fréquence de population de variants

### Update
Ca va pas bouger, j'ai repris de Nirvanna, GRCh37 de 2013 et 38 de 2019

### Format
VCF donc a transformer en parquet ou laisser tel quel
Annotation pour DIAG, AF
Pour les autres EUR_AF, EAS_AF etc.

### Assembly
version dossier hg19 et hg38 disponible
hg19: https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/
hg38: https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20190312_biallelic_SNV_and_INDEL/

## gnomAD
### Resume
<b>gs://gcp-public-data--gnomad/release/</b>
Base de données comprenant le séquençage de milliers de génomes et d'exomes
toutes les données sont disponible sur le cloud google, manière la plus simple et rapide de tout télécharger
Genome 450 Go, possibilité de prende gnomAD genome coupé sur l'exome (10Go)

### Update
A prévoir peut être une nouvelle version après la 4.1.0, sinon version stable donc pas besoin de prévoir une politique d'update
Certaines valeurs de pLi seront sans doute à calculer d'après les annotations des fichiers constraints

### Format
VCF pour les variants ou TSV pour les constraints
Prévoir un formatage en parquet pour les VCFs et from extann pour les TSV avec descriptif des colonnes</br>

<b>Constrains</b></br>
pLi: gnomad_pli</br>
LOEUF: oe_lof_upper</br>
o/e: oe_lof

<b>Statistiques</b></br>
gnomAD Exome et Genome
gnomadaltfreq_popmax: AF_popmax</br>
gnomadaltfreq_asj: AF_asj</br>
gnomadaltfreq_amr: AF_amr</br>
gnomadaltfreq_afr: AF_afr</br>
gnomadaltfreq_eas: AF_eas</br>
gnomadaltfreq_fin: AF_fin</br>
gnomadaltfreq_nfe: AF_nfe</br>
gnomadaltfreq_oth: AF_oth</br>
gnomadaltfreq_sas: AF_sas</br>
gnomadhomcount_all: nhomalt</br>

<b>Pour calculer la vraie AF gnomadaltfreq_all</b> représenté sur gnomad, il faut:</br>
requete sql: AC exome + AC genome/ AN exome + AN genome = gnomadaltfreq_all</br>
<b>Pour les heterozygotes:</b></br>
AC exome + AC genome - 2x(nhomalt exome + nhomalt genome) = gnomadhetcount_all</br>
<b>Pour les hemizygotes:</b></br>
AC_male exome + AC_male genome = gnomadhemcount_all

### Assembly

Version 2.1.1 en GRCh37 avec une version liftover en GRCh38
Version 4.1.0 en GRCh38 only (version Beta d'après gnomAD)

## Decipher
### Resume
BDD variants et phenotype / browser /tools etc
on s'interesse nous au score de haploinsuffisance et a la base de données phenotype des troubles neurodéveloppementaux
https://www.deciphergenomics.org/about/downloads/data

### Update
Aucune

### Format
NDD au format csv / extann, à transformer en extann bed
HI format deja bed mais need parsing et ajouter le header

### Assembly
hg19 pour les HI uniquement

## CADD
### Resume
Prediction sur la pathogénicité d'un variant (modèle ML), sur toutes les bases du génomes pour la version 1.7. (exon score > 20 threshold)

### Update
Sans doute plus d'un an pour la prochaine vu la release 1.7

### Assembly
hg19 et hg38

### Format
| #Chrom | Pos   | Ref | Alt | RawScore | PHRED |
|:-------|:------|:----|:----|:---------|:------|
| chr1   | 10001 | T   | A   | 0.767991 | 7.993 |





HUSDIAGGEN

PZFlag
dbscsnv
gme
n6500NHLBIAA
n6500NHLBIALL
n6500NHLBIEA

TUMSOL
 
    - pzscore-tumsol
    - count_var
    - cosmic
    - mutationtaster_pred
    - outcome

DIAG

    - tags
    - genedesc
    - granthamdist
    - samva