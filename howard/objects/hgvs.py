r"""
HGVS language currently implemented.

HGVS = ALLELE
     | PREFIX_NAME : ALLELE

PREFIX_NAME = TRANSCRIPT
            | TRANSCRIPT '(' GENE ')'

TRANSCRIPT = TRANSCRIPT_NAME
           | TRANSCRIPT_NAME '.' TRANSCRIPT_VERSION

TRANSCRIPT_VERSION = NUMBER

ALLELE = 'c.' CDNA_ALLELE    # cDNA
       | 'g.' GENOMIC_ALLELE # genomic
       | 'm.' MIT_ALLELE     # mitochondrial sequence
       | 'n.' NC_ALLELE      # non-coding RNA reference sequence
       | 'r.' RNA_ALLELE     # RNA sequence (like r.76a>u)
       | 'p.' PROTEIN_ALLELE # protein sequence (like  p.Lys76Asn)

NC_ALLELE =
RNA_ALLELE =
CDNA_ALLELE = CDNA_COORD SINGLE_BASE_CHANGE
            | CDNA_COORD_RANGE MULTI_BASE_CHANGE

GENOMIC_ALLELE =
MIT_ALLELE = COORD SINGLE_BASE_CHANGE
           | COORD_RANGE MULTI_BASE_CHANGE

SINGLE_BASE_CHANGE = CDNA_ALLELE = CDNA_COORD BASE '='        # no change
                   | CDNA_COORD BASE '>' BASE                 # substitution
                   | CDNA_COORD 'ins' BASE                    # 1bp insertion
                   | CDNA_COORD 'del' BASE                    # 1bp deletion
                   | CDNA_COORD 'dup' BASE                    # 1bp duplication
                   | CDNA_COORD 'ins'                         # 1bp insertion
                   | CDNA_COORD 'del'                         # 1bp deletion
                   | CDNA_COORD 'dup'                         # 1bp duplication
                   | CDNA_COORD 'del' BASE 'ins' BASE         # 1bp indel
                   | CDNA_COORD 'delins' BASE                 # 1bp indel

MULTI_BASE_CHANGE = COORD_RANGE 'del' BASES             # deletion
                  | COORD_RANGE 'ins' BASES             # insertion
                  | COORD_RANGE 'dup' BASES             # duplication
                  | COORD_RANGE 'del'                   # deletion
                  | COORD_RANGE 'dup'                   # duplication
                  | COORD_RANGE 'del' BASES 'ins' BASES # indel
                  | COORD_RANGE 'delins' BASES          # indel


AMINO1 = [GAVLIMFWPSTCYNQDEKRH]

AMINO3 = 'Gly' | 'Ala' | 'Val' | 'Leu' | 'Ile' | 'Met' | 'Phe' | 'Trp' | 'Pro'
       | 'Ser' | 'Thr' | 'Cys' | 'Tyr' | 'Asn' | 'Gln' | 'Asp' | 'Glu' | 'Lys'
       | 'Arg' | 'His'

PROTEIN_ALLELE = AMINO3 COORD '='               # no peptide change
               | AMINO1 COORD '='               # no peptide change
               | AMINO3 COORD AMINO3 PEP_EXTRA  # peptide change
               | AMINO1 COORD AMINO1 PEP_EXTRA  # peptide change
               | AMINO3 COORD '_' AMINO3 COORD PEP_EXTRA        # indel
               | AMINO1 COORD '_' AMINO1 COORD PEP_EXTRA        # indel
               | AMINO3 COORD '_' AMINO3 COORD PEP_EXTRA AMINO3 # indel
               | AMINO1 COORD '_' AMINO1 COORD PEP_EXTRA AMINO1 # indel

# A genomic range:
COORD_RANGE = COORD '_' COORD

# A cDNA range:
CDNA_COORD_RANGE = CDNA_COORD '_' CDNA_COORD

# A cDNA coordinate:
CDNA_COORD = COORD_PREFIX COORD
           | COORD_PREFIX COORD OFFSET_PREFIX OFFSET
COORD_PREFIX = '' | '-' | '*'
COORD = NUMBER
OFFSET_PREFIX = '-' | '+'
OFFSET = NUMBER

# Primatives:
NUMBER = "\d+"
BASE = [ACGT]
BASES = BASE+

"""

import re

from .cdna import CDNACoord
from .variant import revcomp

CHROM_PREFIX = 'chr'


CODON_1 = {
    "TTT": "F", "TTC": "F", "TCT": "S", "TCC": "S", "TAT": "Y", "TAC": "Y", "TGT": "C", "TGC": "C", "TTA": "L", "TCA": "S", "TAA": "*", "TGA": "*", "TTG": "L", "TCG": "S", "TAG": "*", "TGG": "W", "CTT": "L", "CTC": "L", "CCT": "P", "CCC": "P", "CAT": "H", "CAC": "H", "CGT": "R", "CGC": "R", "CTA": "L", "CTG": "L", "CCA": "P", "CCG": "P", "CAA": "Q", "CAG": "Q", "CGA": "R", "CGG": "R", "ATT": "I", "ATC": "I", "ACT": "T", "ACC": "T", "AAT": "N", "AAC": "N", "AGT": "S", "AGC": "S", "ATA": "I", "ACA": "T", "AAA": "K", "AGA": "R", "ATG": "M", "ACG": "T", "AAG": "K", "AGG": "R", "GTT": "V", "GTC": "V", "GCT": "A", "GCC": "A", "GAT": "D", "GAC": "D", "GGT": "G", "GGC": "G", "GTA": "V", "GTG": "V", "GCA": "A", "GCG": "A", "GAA": "E", "GAG": "E", "GGA": "G", "GGG": "G"
}
CODON_3 = {
    "TTT": "Phe", "TTC": "Phe", "TCT": "Ser", "TCC": "Ser", "TAT": "Tyr", "TAC": "Tyr", "TGT": "Cys", "TGC": "Cys",
    "TTA": "Leu", "TCA": "Ser", "TAA": "*", "TGA": "*", "TTG": "Leu", "TCG": "Ser", "TAG": "*", "TGG": "Trp",
    "CTT": "Leu", "CTC": "Leu", "CCT": "Pro", "CCC": "Pro", "CAT": "His", "CAC": "His", "CGT": "Arg", "CGC": "Arg",
    "CTA": "Leu", "CTG": "Leu", "CCA": "Pro", "CCG": "Pro", "CAA": "Gln", "CAG": "Gln", "CGA": "Arg", "CGG": "Arg",
    "ATT": "Ile", "ATC": "Ile", "ACT": "Thr", "ACC": "Thr", "AAT": "Asn", "AAC": "Asn", "AGT": "Ser", "AGC": "Ser",
    "ATA": "Ile", "ACA": "Thr", "AAA": "Lys", "AGA": "Arg", "ATG": "Met", "ACG": "Thr", "AAG": "Lys", "AGG": "Arg",
    "GTT": "Val", "GTC": "Val", "GCT": "Ala", "GCC": "Ala", "GAT": "Asp", "GAC": "Asp", "GGT": "Gly", "GGC": "Gly",
    "GTA": "Val", "GTG": "Val", "GCA": "Ala", "GCG": "Ala", "GAA": "Glu", "GAG": "Glu", "GGA": "Gly", "GGG": "Gly"
}
CODON_FULL = {
    "TTT": "Phenylalanine", "TTC": "Phenylalanine", "TCT": "Serine", "TCC": "Serine", "TAT": "Tyrosine", "TAC": "Tyrosine",
    "TGT": "Cysteine", "TGC": "Cysteine", "TTA": "Leucine", "TCA": "Serine", "TAA": "Stop", "TGA": "Stop", "TTG": "Leucine",
    "TCG": "Serine", "TAG": "Stop", "TGG": "Tryptophan", "CTT": "Leucine", "CTC": "Leucine", "CCT": "Proline", "CCC": "Proline",
    "CAT": "Histidine", "CAC": "Histidine", "CGT": "Arginine", "CGC": "Arginine", "CTA": "Leucine", "CTG": "Leucine",
    "CCA": "Proline", "CCG": "Proline", "CAA": "Glutamine", "CAG": "Glutamine", "CGA": "Arginine", "CGG": "Arginine",
    "ATT": "Isoleucine", "ATC": "Isoleucine", "ACT": "Threonine", "ACC": "Threonine", "AAT": "Asparagine", "AAC": "Asparagine",
    "AGT": "Serine", "AGC": "Serine", "ATA": "Isoleucine", "ACA": "Threonine", "AAA": "Lysine", "AGA": "Arginine",
    "ATG": "Methionine", "ACG": "Threonine", "AAG": "Lysine", "AGG": "Arginine", "GTT": "Valine", "GTC": "Valine",
    "GCT": "Alanine", "GCC": "Alanine", "GAT": "Aspartic acid", "GAC": "Aspartic acid", "GGT": "Glycine", "GGC": "Glycine",
    "GTA": "Valine", "GTG": "Valine", "GCA": "Alanine", "GCG": "Alanine", "GAA": "Glutamic acid", "GAG": "Glutamic acid",
    "GGA": "Glycine", "GGG": "Glycine"
}

NUCLEOTIDE_TRANSLATE = {
    "T": "A",
    "A": "T",
    "G": "C",
    "C": "G",
}


# The HGVSRegex class is used for working with regular expressions related to the Human Genome
# Variation Society (HGVS) nomenclature.
class HGVSRegex(object):
    """
    All regular expression for HGVS names.
    """

    # DNA syntax
    # http://www.hgvs.org/mutnomen/standards.html#nucleotide
    BASE = r"[acgtbdhkmnrsvwyACGTBDHKMNRSVWY]|\d+"
    BASES = r"[acgtbdhkmnrsvwyACGTBDHKMNRSVWY]+|\d+"
    DNA_REF = "(?P<ref>" + BASES + ")"
    DNA_ALT = "(?P<alt>" + BASES + ")"

    # Mutation types
    EQUAL = "(?P<mutation_type>=)"
    SUB = "(?P<mutation_type>>)"
    INS = "(?P<mutation_type>ins)"
    DEL = "(?P<mutation_type>del)"
    DUP = "(?P<mutation_type>dup)"

    # Simple coordinate syntax
    COORD_START = r"(?P<start>\d+)"
    COORD_END = r"(?P<end>\d+)"
    COORD_RANGE = COORD_START + "_" + COORD_END

    # cDNA coordinate syntax
    CDNA_COORD = (r"(?P<coord_prefix>|-|\*)(?P<coord>\d+)"
                  r"((?P<offset_prefix>-|\+)(?P<offset>\d+))?")
    CDNA_START = (r"(?P<start>(?P<start_coord_prefix>|-|\*)(?P<start_coord>\d+)"
                  r"((?P<start_offset_prefix>-|\+)(?P<start_offset>\d+))?)")
    CDNA_END = (r"(?P<end>(?P<end_coord_prefix>|-|\*)(?P<end_coord>\d+)"
                r"((?P<end_offset_prefix>-|\+)(?P<end_offset>\d+))?)")
    CDNA_RANGE = CDNA_START + "_" + CDNA_END

    # cDNA allele syntax
    CDNA_ALLELE = [
        # No change
        CDNA_START + EQUAL,
        CDNA_START + DNA_REF + EQUAL,

        # Substitution
        CDNA_START + DNA_REF + SUB + DNA_ALT,

        # 1bp insertion, deletion, duplication
        CDNA_START + INS + DNA_ALT,
        CDNA_START + DEL + DNA_REF,
        CDNA_START + DUP + DNA_REF,
        CDNA_START + DEL,
        CDNA_START + DUP,

        # Insertion, deletion, duplication
        COORD_RANGE + EQUAL,
        CDNA_RANGE + INS + DNA_ALT,
        CDNA_RANGE + DEL + DNA_REF,
        CDNA_RANGE + DUP + DNA_REF,
        CDNA_RANGE + DEL,
        CDNA_RANGE + DUP,

        # Indels
        "(?P<delins>" + CDNA_START + 'del' + DNA_REF + 'ins' + DNA_ALT + ")",
        "(?P<delins>" + CDNA_RANGE + 'del' + DNA_REF + 'ins' + DNA_ALT + ")",
        "(?P<delins>" + CDNA_START + 'delins' + DNA_ALT + ")",
        "(?P<delins>" + CDNA_RANGE + 'delins' + DNA_ALT + ")",
    ]

    CDNA_ALLELE_REGEXES = [re.compile("^" + regex + "$")
                           for regex in CDNA_ALLELE]

    # Peptide syntax
    PEP = "([A-Z]([a-z]{2}))+"
    PEP_REF = "(?P<ref>" + PEP + ")"
    PEP_REF2 = "(?P<ref2>" + PEP + ")"
    PEP_ALT = "(?P<alt>" + PEP + ")"

    PEP_EXTRA = r"(?P<extra>(|=|\?)(|fs))"

    # Peptide allele syntax
    PEP_ALLELE = [
        # No peptide change
        # Example: Glu1161=
        PEP_REF + COORD_START + PEP_EXTRA,

        # Peptide change
        # Example: Glu1161Ser
        PEP_REF + COORD_START + PEP_ALT + PEP_EXTRA,

        # Peptide indel
        # Example: Glu1161_Ser1164?fs
        "(?P<delins>" + PEP_REF + COORD_START + "_" + PEP_REF2 + COORD_END +
        PEP_EXTRA + ")",
        "(?P<delins>" + PEP_REF + COORD_START + "_" + PEP_REF2 + COORD_END +
        PEP_ALT + PEP_EXTRA + ")",
    ]

    PEP_ALLELE_REGEXES = [re.compile("^" + regex + "$")
                          for regex in PEP_ALLELE]

    # Genomic allele syntax
    GENOMIC_ALLELE = [
        # No change
        COORD_START + EQUAL,
        COORD_START + DNA_REF + EQUAL,

        # Substitution
        COORD_START + DNA_REF + SUB + DNA_ALT,

        # 1bp insertion, deletion, duplication
        COORD_START + INS + DNA_ALT,
        COORD_START + DEL + DNA_REF,
        COORD_START + DUP + DNA_REF,
        COORD_START + DEL,
        COORD_START + DUP,

        # Insertion, deletion, duplication
        COORD_RANGE + EQUAL,
        COORD_RANGE + INS + DNA_ALT,
        COORD_RANGE + DEL + DNA_REF,
        COORD_RANGE + DUP + DNA_REF,
        COORD_RANGE + DEL,
        COORD_RANGE + DUP,

        # Indels
        "(?P<delins>" + COORD_START + 'del' + DNA_REF + 'ins' + DNA_ALT + ")",
        "(?P<delins>" + COORD_RANGE + 'del' + DNA_REF + 'ins' + DNA_ALT + ")",
        "(?P<delins>" + COORD_START + 'delins' + DNA_ALT + ")",
        "(?P<delins>" + COORD_RANGE + 'delins' + DNA_ALT + ")",
    ]

    GENOMIC_ALLELE_REGEXES = [re.compile("^" + regex + "$")
                              for regex in GENOMIC_ALLELE]


# The RefSeq standard for naming contigs/transcripts/proteins:
# http://www.ncbi.nlm.nih.gov/books/NBK21091/table/ch18.T.refseq_accession_numbers_and_mole/?report=objectonly  # nopep8
REFSEQ_PREFIXES = [
    ('AC_', 'genomic',
     'Complete genomic molecule, usually alternate assembly'),
    ('NC_', 'genomic',
     'Complete genomic molecule, usually reference assembly'),
    ('NG_', 'genomic', 'Incomplete genomic region'),
    ('NT_', 'genomic', 'Contig or scaffold, clone-based or WGS'),
    ('NW_', 'genomic', 'Contig or scaffold, primarily WGS'),
    ('NS_', 'genomic', 'Environmental sequence'),
    ('NZ_', 'genomic', 'Unfinished WGS'),
    ('NM_', 'mRNA', ''),
    ('NR_', 'RNA', ''),
    ('XM_', 'mRNA', 'Predicted model'),
    ('XR_', 'RNA', 'Predicted model'),
    ('AP_', 'Protein', 'Annotated on AC_ alternate assembly'),
    ('NP_', 'Protein', 'Associated with an NM_ or NC_ accession'),
    ('YP_', 'Protein', ''),
    ('XP_', 'Protein', 'Predicted model, associated with an XM_ accession'),
    ('ZP_', 'Protein', 'Predicted model, annotated on NZ_ genomic records'),
]


REFSEQ_PREFIX_LOOKUP = dict(
    (prefix, (kind, description))
    for prefix, kind, description in REFSEQ_PREFIXES
)


def get_refseq_type(name:str) -> str:
    """
    The `get_refseq_type` function returns the RefSeq type for a given RefSeq name.
    
    :param name: The `name` parameter is a string representing a RefSeq name
    :type name: str
    :return: The function `get_refseq_type` returns the RefSeq type for a given RefSeq name.
    """
    
    prefix = name[:3]
    return REFSEQ_PREFIX_LOOKUP.get(prefix, (None, ''))[0]


# The InvalidHGVSName class is a subclass of ValueError and is used to represent an error when an
# invalid HGVS name is encountered.
class InvalidHGVSName(ValueError):

    def __init__(self, name:str = '', part:str = 'name', reason:str = '') -> None:
        """
        The function initializes an InvalidHGVSName object with a message, name, part, and reason.
        
        :param name: The name parameter is a string that represents the invalid HGVS name. It is the
        name that is considered invalid and does not meet the required criteria
        :type name: str
        :param part: The "part" parameter represents the part of the HGVS (Human Genome Variation
        Society) name that is invalid. It is used to provide more specific information about the error
        that occurred, defaults to name
        :type part: str (optional)
        :param reason: The "reason" parameter is an optional argument that provides additional
        information or context for why the HGVS name is considered invalid. It can be used to provide
        specific details about the error or to explain why the name does not meet the required criteria
        :type reason: str
        """
        
        if name:
            message = 'Invalid HGVS %s "%s"' % (part, name)
        else:
            message = 'Invalid HGVS %s' % part
        if reason:
            message += ': ' + reason
        super(InvalidHGVSName, self).__init__(message)

        self.name = name
        self.part = part
        self.reason = reason


# The HGVSName class is a Python class for handling HGVS (Human Genome Variation Society) names.
class HGVSName(object):
    """
    Represents a HGVS variant name.
    """

    def __init__(self, name:str = '', prefix:str = '', chrom:str = '', transcript:str = '', transcript_protein:str = None, gene:str = '', exon:str = None, kind:str = '', mutation_type:str = None, start:int = 0, end:int = 0, ref_allele:str = '', ref2_allele:str = '', alt_allele:str = '', cdna_start:int = None, cdna_end:int = None, pep_extra:str = ''):
        """
        The function is a constructor that initializes various attributes of an object and parses a
        given name to populate those attributes.
        
        :param name: The full HGVS name of the variant
        :type name: str
        :param prefix: The `prefix` parameter is a string that is used as a prefix for the HGVS name. It
        can be used to indicate additional information or context about the variant
        :type prefix: str
        :param chrom: The `chrom` parameter represents the chromosome where the mutation occurs. It is a
        string that specifies the chromosome number or identifier
        :type chrom: str
        :param transcript: The `transcript` parameter represents the transcript ID or name associated
        with the mutation. It is used to specify the specific transcript in which the mutation occurs
        :type transcript: str
        :param transcript_protein: The `transcript_protein` parameter is used to store information about
        the protein associated with the transcript. It can be used to specify the protein variant or
        isoform that is affected by the mutation
        :type transcript_protein: str
        :param gene: The "gene" parameter represents the gene associated with the variant. It is a
        string that specifies the gene name or identifier
        :type gene: str
        :param exon: The `exon` parameter represents the exon number or range in which the mutation
        occurs. It is used to specify the location of the mutation within the transcript
        :type exon: str
        :param kind: The "kind" parameter is used to specify the type of variant or mutation. It can be
        a string that represents the kind of mutation, such as "substitution", "deletion", "insertion",
        etc. This parameter helps to categorize and describe the type of mutation being represented by
        the
        :type kind: str
        :param mutation_type: The `mutation_type` parameter is used to specify the type of mutation. It
        can be a string that represents the type of mutation, such as "SNP" (single nucleotide
        polymorphism), "DEL" (deletion), "INS" (insertion), etc
        :type mutation_type: str
        :param start: The `start` parameter represents the starting position of the mutation or variant
        in the genomic sequence. It is an integer value that indicates the position of the mutation or
        variant on the genomic sequence. If not provided, it defaults to 0, defaults to 0
        :type start: int (optional)
        :param end: The "end" parameter represents the end position of the mutation or variant. It is an
        integer value that indicates the position of the mutation or variant on the genomic sequence,
        defaults to 0
        :type end: int (optional)
        :param ref_allele: The `ref_allele` parameter represents the reference allele in a genetic
        mutation. It is the allele that is present in the reference genome at a specific position
        :type ref_allele: str
        :param ref2_allele: The `ref2_allele` parameter represents the reference allele at the end of a
        peptide indel. In the context of genetic mutations, an indel refers to the insertion or deletion
        of nucleotides in a DNA sequence. The `ref2_allele` specifically represents the reference allele
        that is
        :type ref2_allele: str
        :param alt_allele: The `alt_allele` parameter represents the alternate allele in a genetic
        mutation. In genetics, an allele is one of the possible forms of a gene. In the context of this
        code, `alt_allele` is used to store the alternate allele that is present in a mutation
        :type alt_allele: str
        :param cdna_start: The `cdna_start` parameter is used to specify the start position of the
        mutation in the cDNA sequence. It is an optional parameter and if not provided, it will be set
        to a default value of `CDNACoord()`
        :type cdna_start: int
        :param cdna_end: The `cdna_end` parameter is used to store the end coordinate of the cDNA
        (complementary DNA) sequence. It is an optional parameter and if not provided, it will be
        initialized as a `CDNACoord` object. The `CDNACoord` object is likely a
        :type cdna_end: int
        :param pep_extra: The `pep_extra` parameter is a string that represents any additional
        information related to the protein. It is used in the context of protein-specific fields
        :type pep_extra: str
        """

        # Full HGVS name.
        self.name = name

        # Name parts.
        self.prefix = prefix
        self.chrom = chrom
        self.transcript = transcript
        self.transcript_protein = transcript_protein
        self.gene = gene
        self.exon = exon
        self.kind = kind
        self.mutation_type = mutation_type
        self.start = start
        self.end = end
        self.ref_allele = ref_allele  # reference allele
        self.ref2_allele = ref2_allele  # reference allele at end of pep indel
        self.alt_allele = alt_allele  # alternate allele

        # cDNA-specific fields
        self.cdna_start = cdna_start if cdna_start else CDNACoord()
        self.cdna_end = cdna_end if cdna_end else CDNACoord()

        # Protein-specific fields
        self.pep_extra = pep_extra

        if name:
            self.parse(name)


    def parse(self, name:str) -> None:
        """
        The `parse` function is used to split an HGVS name into a prefix and allele, and then validate
        the parsed components.
        
        :param name: The `name` parameter is a string that represents an HGVS name. It is the input to
        the `parse` function and is used to parse the HGVS name by splitting it into a prefix and allele
        :type name: str
        """
        
        # Does HGVS name have transcript/gene prefix?
        if ':' in name:
            prefix, allele = name.split(':', 1)
        else:
            prefix = ''
            allele = name

        self.name = name

        # Parse prefix and allele.
        self.parse_allele(allele)
        self.parse_prefix(prefix)
        self._validate()


    def parse_prefix(self, prefix:str):
        """
        The `parse_prefix` function is used to parse a HGVS prefix (gene/transcript/chromosome) and
        assign the parsed values to the corresponding attributes of the object.
        
        :param prefix: The `prefix` parameter is a string that represents a HGVS prefix, which can be a
        gene, transcript, or chromosome identifier. It is used to determine the type of prefix and
        assign the parsed values to the corresponding attributes of the object
        :type prefix: str
        :return: The function `parse_prefix` returns the parsed values for the transcript and gene
        attributes, or sets the chrom or gene attributes based on the given prefix.
        """

        self.prefix = prefix

        # No prefix.
        if prefix == '':
            self.chrom = ''
            self.transcript = ''
            self.gene = ''
            return

        # Transcript and gene given with parens.
        # example: NM_007294.3(BRCA1):c.2207A>C
        match = re.match(r"^(?P<transcript>[^(]+)\((?P<gene>[^)]+)\)$", prefix)
        if match:
            self.transcript = match.group('transcript')
            self.gene = match.group('gene')
            return

        # Transcript and gene given with braces.
        # example: BRCA1{NM_007294.3}:c.2207A>C
        match = re.match(r"^(?P<gene>[^{]+){(?P<transcript>[^}]+)}$", prefix)
        if match:
            self.transcript = match.group('transcript')
            self.gene = match.group('gene')
            return

        # Determine using Ensembl type.
        if prefix.startswith('ENST'):
            self.transcript = prefix
            return

        # Determine using LRG type.
        if prefix.startswith('LRG_'):
            self.transcript = prefix
            return

        # Determine using refseq type.
        refseq_type = get_refseq_type(prefix)
        if refseq_type in ('mRNA', 'RNA'):
            self.transcript = prefix
            return

        # Determine using refseq type.
        if prefix.startswith(CHROM_PREFIX) or refseq_type == 'genomic':
            self.chrom = prefix
            return

        # Assume gene name.
        self.gene = prefix

    def parse_allele(self, allele:str) -> None:
        """
        The function `parse_allele` parses a HGVS allele description and determines the kind of HGVS
        name (c., p., g., etc.) and the mutation type.
        
        Some examples include:
          cDNA substitution: c.101A>C,
          cDNA indel: c.3428delCinsTA, c.1000_1003delATG, c.1000_1001insATG
          No protein change: p.Glu1161=
          Protein change: p.Glu1161Ser
          Protein frameshift: p.Glu1161_Ser1164?fs
          Genomic substitution: g.1000100A>T
          Genomic indel: g.1000100_1000102delATG

        :param allele: The `allele` parameter is a string that represents a HGVS allele description. It
        can contain various types of mutations, such as cDNA substitutions, cDNA indels, protein
        changes, protein frameshifts, genomic substitutions, and genomic indels. The purpose of the
        `parse_allele`
        :type allele: str
        """

        if '.' not in allele:
            raise InvalidHGVSName(allele, 'allele',
                                  'expected kind "c.", "p.", "g.", etc')

        # Determine HGVS name kind.
        kind, details = allele.split('.', 1)
        self.kind = kind
        self.mutation_type = None

        if kind in ("c", 'n'):
            self.parse_cdna(details)
            if kind == 'n':  # Ensure no 3'UTR or 5'UTR coords in non-coding
                if self.cdna_start.coord < 0:
                    raise InvalidHGVSName(allele, 'allele',
                                          "Non-coding transcript cannot contain negative (5'UTR) coordinates")
                if self.cdna_start.landmark == 'cdna_stop' or self.cdna_end and self.cdna_end.landmark == 'cdna_stop':
                    raise InvalidHGVSName(allele, 'allele',
                                          "Non-coding transcript cannot contain '*' (3'UTR) coordinates")
        elif kind == "p":
            self.parse_protein(details)
        elif kind in ("g", 'm'):
            self.parse_genome(details)
        else:
            raise NotImplementedError("unknown kind: %s" % allele)


    def parse_cdna(self, details:str) -> None:
        """
        The function `parse_cdna` is used to parse a HGVS cDNA name and extract information such as
        mutation type, coordinates, and alleles.
        
        Some examples include:
          Substitution: 101A>C,
          Indel: 3428delCinsTA, 1000_1003delATG, 1000_1001insATG

        :param details: The `details` parameter is a string that represents a HGVS cDNA name. It
        contains information about a genetic mutation, such as a substitution or an indel, along with
        the specific coordinates and alleles involved in the mutation
        :type details: str
        :return: None.
        """

        for regex in HGVSRegex.CDNA_ALLELE_REGEXES:
            match = re.match(regex, details)
            if match:
                groups = match.groupdict()

                # Parse mutation type.
                if groups.get('delins'):
                    self.mutation_type = 'delins'
                else:
                    self.mutation_type = groups['mutation_type']

                # Parse coordinates.
                self.cdna_start = CDNACoord(string=groups.get('start'))
                if groups.get('end'):
                    self.cdna_end = CDNACoord(string=groups.get('end'))
                else:
                    self.cdna_end = CDNACoord(string=groups.get('start'))

                # Parse alleles.
                self.ref_allele = groups.get('ref', '')
                self.alt_allele = groups.get('alt', '')

                # Convert numerical allelles.
                if self.ref_allele.isdigit():
                    self.ref_allele = "N" * int(self.ref_allele)
                if self.alt_allele.isdigit():
                    self.alt_allele = "N" * int(self.alt_allele)

                # Convert duplication alleles.
                if self.mutation_type == "dup":
                    self.alt_allele = self.ref_allele * 2

                # Convert no match alleles.
                if self.mutation_type == "=":
                    self.alt_allele = self.ref_allele
                return

        raise InvalidHGVSName(details, 'cDNA allele')


    def parse_protein(self, details:str) -> None:
        """
        The function `parse_protein` is used to parse a HGVS protein name and extract information such
        as mutation type, coordinates, alleles, and additional details.
        
        Some examples include:
          No change: Glu1161=
          Change: Glu1161Ser
          Frameshift: Glu1161_Ser1164?fs

        :param details: The `details` parameter is a string that represents a HGVS protein name. It
        contains information about a protein mutation, such as the amino acid change and the position of
        the mutation
        :type details: str
        :return: The method `parse_protein` does not return anything. It updates the instance variables
        of the object it is called on.
        """

        for regex in HGVSRegex.PEP_ALLELE_REGEXES:
            match = re.match(regex, details)
            if match:
                groups = match.groupdict()

                # Parse mutation type.
                if groups.get('delins'):
                    self.mutation_type = 'delins'
                else:
                    self.mutation_type = '>'

                # Parse coordinates.
                self.start = int(groups.get('start'))
                if groups.get('end'):
                    self.end = int(groups.get('end'))
                else:
                    self.end = self.start

                # Parse alleles.
                self.ref_allele = groups.get('ref', '')
                if groups.get('ref2'):
                    self.ref2_allele = groups.get('ref2')
                    self.alt_allele = groups.get('alt', '')
                else:
                    # If alt is not given, assume matching with ref
                    self.ref2_allele = self.ref_allele
                    self.alt_allele = groups.get(
                        'alt', self.ref_allele)

                self.pep_extra = groups.get('extra')
                return

        raise InvalidHGVSName(details, 'protein allele')


    def parse_genome(self, details:str) -> None:
        """
        The function `parse_genome` is used to parse a HGVS genomic name and extract information such as
        mutation type, coordinates, and alleles.
        
        Some examples include:
          Substitution: 1000100A>T
          Indel: 1000100_1000102delATG

        :param details: The `details` parameter is a string that represents a HGVS genomic name. It
        contains information about a genomic mutation, such as a substitution or an indel
        :type details: str
        :return: None.
        """

        for regex in HGVSRegex.GENOMIC_ALLELE_REGEXES:
            match = re.match(regex, details)
            if match:
                groups = match.groupdict()

                # Parse mutation type.
                if groups.get('delins'):
                    self.mutation_type = 'delins'
                else:
                    self.mutation_type = groups['mutation_type']

                # Parse coordinates.
                self.start = int(groups.get('start'))
                if groups.get('end'):
                    self.end = int(groups.get('end'))
                else:
                    self.end = self.start

                # Parse alleles.
                self.ref_allele = groups.get('ref', '')
                self.alt_allele = groups.get('alt', '')

                # Convert numerical alleles.
                if self.ref_allele.isdigit():
                    self.ref_allele = "N" * int(self.ref_allele)
                if self.alt_allele.isdigit():
                    self.alt_allele = "N" * int(self.alt_allele)

                # Convert duplication alleles.
                if self.mutation_type == "dup":
                    self.alt_allele = self.ref_allele * 2

                # Convert no match alleles.
                if self.mutation_type == "=":
                    self.alt_allele = self.ref_allele
                return

        raise InvalidHGVSName(details, 'genomic allele')


    def _validate(self) -> None:
        """
        The function checks for internal inconsistencies in the representation of coordinates.
        """
        
        if self.start > self.end:
            raise InvalidHGVSName(reason="Coordinates are nonincreasing")


    def __repr__(self) -> str:
        """
        The function returns a string representation of an HGVSName object.
        :return: The `__repr__` method is returning a string representation of the object. If the
        `format` method is implemented, it will return a string in the format
        "HGVSName('formatted_string')". If the `format` method is not implemented, it will return a
        string in the format "HGVSName('name')".
        """

        try:
            return "HGVSName('%s')" % self.format()
        except NotImplementedError:
            return "HGVSName('%s')" % self.name


    def __unicode__(self) -> str:
        """
        The function returns a formatted string representation of the object.
        :return: The `__unicode__` method is returning the result of the `format()` method.
        """

        return self.format()


    def format(self, use_prefix:bool = True, use_gene:bool = True, use_exon:bool = False, use_protein:bool = False, full_format=False, use_version:bool = False) -> str:
        """
        The `format` function generates a HGVS name as a string based on various formatting options.
        
        :param use_prefix: A boolean indicating whether to include the prefix in the HGVS name. If set
        to True, the prefix will be included in the HGVS name. If set to False, the prefix will be
        excluded. The default value is True, defaults to True
        :type use_prefix: bool (optional)
        :param use_gene: A boolean indicating whether to include the gene name in the HGVS name. If set
        to True, the gene name will be included in the HGVS name. If set to False, the gene name will
        not be included. The default value is True, defaults to True
        :type use_gene: bool (optional)
        :param use_exon: A boolean indicating whether to include exon information in the HGVS name. If
        set to True, exon information will be included in the HGVS name. If set to False, exon
        information will not be included, defaults to False
        :type use_exon: bool (optional)
        :param use_protein: A boolean indicating whether to include the protein change in the HGVS name.
        If set to True, the protein change will be included in the HGVS name. If set to False, the
        protein change will not be included, defaults to False
        :type use_protein: bool (optional)
        :param full_format: A boolean parameter that determines whether the full format of the allele
        should be included in the output. If set to True, and if the allele is not a protein variant,
        the allele will be appended with ':p.' followed by the formatted protein variant, defaults to
        False (optional)
        :param use_version: A boolean parameter that determines whether to include the version number in
        the formatted HGVS name. If set to True, the version number will be included in the output. If
        set to False, the version number will not be included, defaults to False
        :type use_version: bool (optional)
        :return: a HGVS name as a string.
        """

        if use_protein and self.format_protein():
            allele = 'p.' + self.format_protein()
        else:
            if self.kind in ('c', 'n'):
                allele = self.kind + '.' + self.format_cdna()
            elif self.kind == 'p':
                allele = 'p.' + self.format_protein()
            elif self.kind in ('g', 'm'):
                allele = self.kind + '.' + self.format_genome()
            else:
                raise NotImplementedError("not implemented: '%s'" % self.kind)

        if full_format and not use_protein and self.format_protein():
            allele += ':p.' + self.format_protein()
        
        prefix = self.format_prefix(use_gene=use_gene, use_exon=use_exon, use_protein=use_protein, full_format=full_format, use_version=use_version) if use_prefix else ''
        
        if prefix:
            return prefix + ':' + allele
        else:
            return allele


    def format_prefix(self, use_gene:bool = True, use_exon:bool = False, use_protein:bool = False, full_format:bool = False, use_version:bool = False) -> str:
        """
        The `format_prefix` function generates an HGVS transcript/gene prefix based on various
        parameters.
        
        :param use_gene: A boolean parameter that determines whether to include the gene name in the
        prefix. If set to True, the gene name will be included in the prefix. If set to False, the gene
        name will not be included in the prefix. The default value is True, defaults to True
        :type use_gene: bool (optional)
        :param use_exon: A boolean parameter that determines whether to include the exon information in
        the prefix. If set to True, the exon information will be included in the prefix. If set to
        False, the exon information will not be included, defaults to False
        :type use_exon: bool (optional)
        :param use_protein: A boolean indicating whether to use the protein transcript instead of the
        nucleotide transcript if available. If set to True, the protein transcript will be used. If set
        to False, the nucleotide transcript will be used. The default value is False, defaults to False
        :type use_protein: bool (optional)
        :param full_format: A boolean parameter that determines whether to generate the full HGVS name
        with transcript/gene prefix or not. If set to True, the full format will be generated. If set to
        False, only the transcript/gene prefix will be generated, defaults to False
        :type full_format: bool (optional)
        :param use_version: A boolean parameter that determines whether to include the version number in
        the transcript prefix. If set to True, the version number will be included in the prefix (e.g.,
        NM_007294.3). If set to False, only the transcript ID without the version number will be
        included in the prefix, defaults to False
        :type use_version: bool (optional)
        :return: The function `format_prefix` returns a formatted HGVS transcript/gene prefix as a
        string.
        """

        if full_format:

            prefix = []
            if self.gene:
                prefix.append(self.gene)
            if self.transcript:
                prefix.append(self.transcript)
            if self.transcript_protein:
                prefix.append(self.transcript_protein)
            if self.exon: # TODO
                prefix.append("exon"+str(self.exon))
            return ":".join(prefix)

        else:

            if self.kind in ('g', 'm'):
                if self.chrom:
                    return self.chrom

            if self.transcript:
                if use_protein and self.transcript_protein:
                    transcript = self.transcript_protein
                else:
                    transcript = self.transcript
                if not use_version and transcript:
                    transcript = transcript.split('.')[0]
                if use_gene and self.gene:
                    return '%s(%s)' % (transcript, self.gene)
                elif use_exon and self.exon:
                    return '%s(exon%s)' % (transcript, self.exon)
                else:
                    return transcript
            else:
                if use_gene:
                    return self.gene
                else:
                    return ''


    def format_cdna_coords(self) -> str:
        """
        The function `format_cdna_coords` generates a string representing HGVS cDNA coordinates,
        returning either the start coordinate or a string in the format "start_end" depending on whether
        the start and end coordinates are the same or not.
        :return: a string representing the cDNA coordinates. If the start and end coordinates are the
        same, it returns just the start coordinate. Otherwise, it returns a string in the format
        "start_end".
        """
        
        # Format coordinates.
        if self.cdna_start == self.cdna_end:
            return str(self.cdna_start)
        else:
            return "%s_%s" % (self.cdna_start, self.cdna_end)


    def format_dna_allele(self) -> str:
        """
        The function `format_dna_allele` generates an HGVS DNA allele based on the mutation type and
        alleles provided.
        :return: The function `format_dna_allele` returns a string representing the HGVS DNA allele. The
        specific format of the returned string depends on the value of the `mutation_type` attribute of
        the object. The possible return values are:
        """
        
        if self.mutation_type == '=':
            # No change.
            # example: 101A=
            return self.ref_allele + '='

        if self.mutation_type == '>':
            # SNP.
            # example: 101A>C
            return self.ref_allele + '>' + self.alt_allele

        elif self.mutation_type == 'delins':
            # Indel.
            # example: 112_117delAGGTCAinsTG, 112_117delinsTG
            return 'del' + self.ref_allele + 'ins' + self.alt_allele

        elif self.mutation_type in ('del', 'dup'):
            # Delete, duplication.
            # example: 1000_1003delATG, 1000_1003dupATG
            return self.mutation_type + self.ref_allele

        elif self.mutation_type == 'ins':
            # Insert.
            # example: 1000_1001insATG
            return self.mutation_type + self.alt_allele

        elif self.mutation_type == 'inv':
            return self.mutation_type

        else:
            raise AssertionError(
                "unknown mutation type: '%s'" % self.mutation_type)


    def format_cdna(self) -> str:
        """
        The function "format_cdna" generates an HGVS cDNA allele by combining the cDNA coordinates and
        the DNA allele.

        Some examples include:
          Substitution: 101A>C,
          Indel: 3428delCinsTA, 1000_1003delATG, 1000_1001insATG

        :return: a string that represents the HGVS cDNA allele.
        """
        
        return self.format_cdna_coords() + self.format_dna_allele()


    def format_protein(self) -> str:
        """
        The `format_protein` function generates an HGVS protein name based on different scenarios such
        as no change, change, frameshift, and range change.

        Some examples include:
          No change: Glu1161=
          Change: Glu1161Ser
          Frameshift: Glu1161_Ser1164?fs

        :return: The method `format_protein` returns a string representing the HGVS protein name.
        """
        
        if (self.start == self.end and
                self.ref_allele == self.ref2_allele == self.alt_allele):
            # Match.
            # Example: Glu1161=
            pep_extra = self.pep_extra if self.pep_extra else '='
            return self.ref_allele + str(self.start) + pep_extra

        elif (self.start == self.end and
              self.ref_allele == self.ref2_allele and
              self.ref_allele != self.alt_allele):
            # Change.
            # Example: Glu1161Ser
            return (self.ref_allele + str(self.start) +
                    self.alt_allele + self.pep_extra)

        elif self.start != self.end:
            # Range change.
            # Example: Glu1161_Ser1164?fs
            return (self.ref_allele + str(self.start) + '_' +
                    self.ref2_allele + str(self.end) +
                    self.pep_extra)
        elif self.pep_extra:
            return self.pep_extra

        else:
            return None


    def format_coords(self) -> str:
        """
        The function `format_coords` generates a string representation of HGVS cDNA coordinates.
        :return: a string that represents the HGVS cDNA coordinates. If the start and end coordinates
        are the same, it returns just the start coordinate. Otherwise, it returns a string in the format
        "start_end".
        """
        
        # Format coordinates.
        if self.start == self.end:
            return str(self.start)
        else:
            return "%s_%s" % (self.start, self.end)


    def format_genome(self) -> str:
        """
        The function "format_genome" generates an HGVS genomic allele by combining the formatted
        coordinates and DNA allele.

        Some examples include:
          Substitution: 1000100A>T
          Indel: 1000100_1000102delATG

        :return: a string that represents the HGVS genomic allele.
        """
        
        return self.format_coords() + self.format_dna_allele()


    def get_raw_coords(self, transcript:object = None) -> tuple:
        """
        The function `get_raw_coords` returns the genomic coordinates based on the given transcript or
        the provided chromosomal coordinates.
        
        :param transcript: The `transcript` parameter is an object that represents a transcript. It is
        used to retrieve genomic coordinates based on the type of HGVS name (`self.kind`). The
        `transcript` object should have the following attributes and methods:
        :type transcript: object
        :return: a tuple containing the genomic coordinates. The tuple consists of three elements: the
        chromosome, the start position, and the end position.
        """
        
        if self.kind in ('c', 'n'):
            chrom = transcript.tx_position.chrom
            start = transcript.cdna_to_genomic_coord(self.cdna_start)
            end = transcript.cdna_to_genomic_coord(self.cdna_end)

            if not transcript.tx_position.is_forward_strand:
                start, end = end, start

            if start > end:
                raise AssertionError(
                    "cdna_start cannot be greater than cdna_end")
        elif self.kind in ('g', 'm'):
            chrom = self.chrom
            start = self.start
            end = self.end
        else:
            raise NotImplementedError(
                'Coordinates are not available for this kind of HGVS name "%s"'
                % self.kind)

        # Check coordinate span is equal to reference bases
        if self.ref_allele:
            coordinate_span = end - start + 1  # Ref will always be >=1 base
            ref_length = len(self.ref_allele)
            if coordinate_span != ref_length:
                raise InvalidHGVSName("Coordinate span (%d) not equal to ref length %d" % (coordinate_span, ref_length))

        return chrom, start, end


    def get_ref_coords(self, transcript:object = None) -> tuple:
        """
        The function "get_ref_coords" returns the genomic coordinates of the reference allele, taking
        into account different mutation types.
        
        :param transcript: The `transcript` parameter is an optional object that represents a transcript
        or gene. It is used to retrieve the genomic coordinates of the reference allele
        :type transcript: object
        :return: a tuple containing the genomic coordinates of the reference allele. The tuple consists
        of three elements: the chromosome, the start position, and the end position.
        """
        
        chrom, start, end = self.get_raw_coords(transcript)

        if self.mutation_type == "ins":
            # Inserts have empty interval.
            if start < end:
                start += 1
                end -= 1
            else:
                end = start - 1

        elif self.mutation_type == "dup":
            end = start - 1
        return chrom, start, end


    def get_vcf_coords(self, transcript:object = None) -> tuple:
        """
        The function "get_vcf_coords" returns the genomic coordinates of the reference allele in
        VCF-style, with left-padding for indels.
        
        :param transcript: The `transcript` parameter is an object that represents a transcript or gene.
        It is used to retrieve the genomic coordinates of the reference allele
        :type transcript: object
        :return: a tuple containing the genomic coordinates of the reference allele in VCF-style. The
        tuple consists of three elements: the chromosome, the start position, and the end position.
        """
        
        chrom, start, end = self.get_ref_coords(transcript)

        # Inserts and deletes require left-padding by 1 base
        if self.mutation_type in ("=", ">"):
            pass
        elif self.mutation_type in ("del", "ins", "dup", "delins"):
            # Indels have left-padding.
            start -= 1
        else:
            raise NotImplementedError("Unknown mutation_type '%s'" %
                                      self.mutation_type)
        return chrom, start, end
    

    def get_ref_alt(self, is_forward_strand:bool = True, raw_dup_alleles:bool = False) -> tuple:
        """
        The function `get_ref_alt` returns the reference and alternate alleles, with an option to modify
        duplications to look like inserts.
        
        :param is_forward_strand: The parameter `is_forward_strand` is a boolean flag that indicates
        whether the alleles should be returned for the forward strand or the reverse complement strand.
        If `is_forward_strand` is `True`, the alleles will be returned as is. If `is_forward_strand` is
        `False`,, defaults to True
        :type is_forward_strand: bool (optional)
        :param raw_dup_alleles: The `raw_dup_alleles` parameter is a boolean flag that determines
        whether the raw values of duplicated alleles should be returned. By default, it is set to
        `False`, which means that if the mutation type is a duplication (`dup`), the reference allele
        will be represented as an empty string, defaults to False
        :type raw_dup_alleles: bool (optional)
        :return: The function `get_ref_alt` returns a tuple containing the reference and alternate
        alleles.
        """

        if self.kind == 'p':
            raise NotImplementedError(
                'get_ref_alt is not implemented for protein HGVS names')
        alleles = [self.ref_allele, self.alt_allele]

        # Represent duplications are inserts.
        if not raw_dup_alleles and self.mutation_type == "dup":
            alleles[0] = ""
            alleles[1] = alleles[1][:len(alleles[1]) // 2]

        if is_forward_strand:
            return alleles
        else:
            return tuple(map(revcomp, alleles))
