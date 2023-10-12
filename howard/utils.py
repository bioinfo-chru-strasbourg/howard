"""
Helper functions.
"""

from __future__ import absolute_import
from __future__ import unicode_literals

import operator
import re
import time
import polars as pl

from howard.objects.transcript import Transcript, CDNA_Match, Exon
from howard.objects.cdna import CDNACoord, CDNA_START_CODON, CDNA_STOP_CODON
from howard.objects.genome import GenomeSubset
from howard.objects.hgvs import HGVSName, InvalidHGVSName, CODON_1, CODON_3, CODON_FULL
from howard.objects.variant import Position, justify_indel, normalize_variant, revcomp


def read_refgene(infile):
    """
    The function "read_refgene" reads a genePred file with an extra column at the front and returns the
    genePred data.

    refGene = genePred with extra column at front (and ignored ones after)
    
    :param infile: The input file containing the refGene data
    :return: the result of calling the function `read_genepred` with the argument `infile` and the
    keyword argument `skip_first_column` set to `True`.
    """

    return read_genepred(infile, skip_first_column=True)


def read_genepred(infile, skip_first_column=False):
    """
    The function `read_genepred` reads a file in GenePred extension format and yields a dictionary for
    each line, containing information about a gene.
    
    :param infile: The `infile` parameter is the input file object that contains the gene annotation
    data in the GenePred format. It is used to read the lines of the file and extract the necessary
    information
    :param skip_first_column: The `skip_first_column` parameter is a boolean flag that determines
    whether to skip the first column of the input file when parsing the genePred format. By default, it
    is set to `False`, which means the first column (usually the transcript ID) will be included in the
    output. If you, defaults to False (optional)
    
    GenePred extension format:
    http://genome.ucsc.edu/FAQ/FAQformat.html#GenePredExt

    Column definitions:
    0. string name;                 "Name of gene (usually transcript_id from GTF)"
    1. string chrom;                "Chromosome name"
    2. char[1] strand;              "+ or - for strand"
    3. uint txStart;                "Transcription start position"
    4. uint txEnd;                  "Transcription end position"
    5. uint cdsStart;               "Coding region start"
    6. uint cdsEnd;                 "Coding region end"
    7. uint exonCount;              "Number of exons"
    8. uint[exonCount] exonStarts;  "Exon start positions"
    9. uint[exonCount] exonEnds;    "Exon end positions"
    10. uint id;                    "Unique identifier"
    11. string name2;               "Alternate name (e.g. gene_id from GTF)"
    """
    for line in infile:
        # Skip comments.
        if line.startswith('#') or line.startswith('bin'):
            continue
        row = line.rstrip('\n').split('\t')
        if skip_first_column:
            row = row[1:]

        # Skip trailing ,
        exon_starts = list(map(int, row[8].split(',')[:-1]))
        exon_ends = list(map(int, row[9].split(',')[:-1]))
        exons = list(zip(exon_starts, exon_ends))

        yield {
            'chrom': row[1],
            'start': int(row[3]),
            'end': int(row[4]),
            'id': row[0],
            'strand': row[2],
            'cds_start': int(row[5]),
            'cds_end': int(row[6]),
            'gene_name': row[11],
            'exons': exons,
        }


def make_transcript(transcript_json):
    """
    The function `make_transcript` takes a JSON object representing a transcript and creates a
    Transcript object from it.
    
    :param transcript_json: The `transcript_json` parameter is a JSON object that contains information
    about a transcript. It should have the following keys:
    :return: a Transcript object.
    """

    transcript_name = transcript_json['id']
    if '.' in transcript_name:
        name, version = transcript_name.split('.')
    else:
        name, version = transcript_name, None

    transcript = Transcript(
        name=name,
        version=int(version) if version is not None else None,
        gene=transcript_json['gene_name'],
        tx_position=Position(
            transcript_json['chrom'],
            transcript_json['start'],
            transcript_json['end'],
            transcript_json['strand'] == '+'),
        cds_position=Position(
            transcript_json['chrom'],
            transcript_json['cds_start'],
            transcript_json['cds_end'],
            transcript_json['strand'] == '+'),
        start_codon_transcript_pos=transcript_json.get("start_codon_transcript_pos"),
        stop_codon_transcript_pos=transcript_json.get("stop_codon_transcript_pos"),
    )

    exons = transcript_json['exons']
    exons.sort(key=operator.itemgetter(0))
    cdna_match = transcript_json.get('cdna_match', [])
    cdna_match.sort(key=operator.itemgetter(0))

    if not transcript.tx_position.is_forward_strand:
        exons.reverse()
        cdna_match.reverse()

    # We don't use exons, but run everything through cDNA match so there's just 1 path
    # exons are treated as a perfect cDNA match
    if not cdna_match:
        cdna_match = json_perfect_exons_to_cdna_match(exons)

    for number, (exon_start, exon_end, cdna_start, cdna_end, gap) in enumerate(cdna_match, 1):
        transcript.cdna_match.append(CDNA_Match(transcript=transcript,
                                                tx_position=Position(
                                                    transcript_json['chrom'],
                                                    exon_start,
                                                    exon_end,
                                                    transcript_json['strand'] == '+'),
                                                cdna_start=cdna_start,
                                                cdna_end=cdna_end,
                                                gap=gap,
                                                number=number))

    return transcript


def json_perfect_exons_to_cdna_match(ordered_exons, single=False):
    """
    The function `json_perfect_exons_to_cdna_match` converts a list of ordered exons into a list of cDNA
    matches, where each match consists of the start and end positions of the exon, the start and end
    positions of the corresponding cDNA sequence, and an optional gap list.

    Perfectly matched exons are basically a no-gap case of cDNA match
    single - use a single cDNA match (deletions for introns) - this is currently broken do not use
    
    :param ordered_exons: A list of tuples representing the start and end positions of exons in a gene
    sequence. The exons should be ordered based on their position in the gene
    :param single: The `single` parameter is a boolean flag that determines whether to use a single cDNA
    match or not. If `single` is set to `True`, the function will create a single cDNA match by
    considering deletions for introns. If `single` is set to `False` (, defaults to False (optional)
    :return: a list of lists, where each inner list represents a cDNA match. Each inner list contains
    the start and end positions of the exon, the start and end positions of the corresponding cDNA
    match, and a string representing any gaps (intron lengths) between exons.
    """
    
    cdna_match = []
    if single:
        ordered_exons = list(ordered_exons)
        start = ordered_exons[0][0]
        end = ordered_exons[-1][1]
        last_exon_end = None
        gap_list = []
        cdna_length = 0
        for (exon_start, exon_end) in ordered_exons:
            # end up looking like "M D M D (M=exon, D=intron length)"
            if last_exon_end:
                intron_length = abs(exon_start - last_exon_end)
                gap_list.append("D%d" % intron_length)
            exon_length = exon_end - exon_start
            cdna_length += exon_length
            gap_list.append("M%d" % exon_length)
            last_exon_end = exon_end
        cdna_match = [[start, end, 1, cdna_length, " ".join(gap_list)]]
    else:
        cdna_start = 1
        for (exon_start, exon_end) in ordered_exons:
            exon_length = exon_end - exon_start
            cdna_end = cdna_start + exon_length - 1
            cdna_match.append([exon_start, exon_end, cdna_start, cdna_end, None])
            cdna_start = cdna_end + 1
    return cdna_match


def read_transcripts(refgene_file):
    """
    The function `read_transcripts` reads all transcripts in a RefGene file and returns them as a
    dictionary.
    
    :param refgene_file: The `refgene_file` parameter is the file path to a RefGene file. This file
    contains information about gene transcripts, such as their names, full names, and other relevant
    details. The `read_transcripts` function reads this file and returns a dictionary of transcripts,
    where the keys are the
    :return: a dictionary of transcripts.
    """
    
    transcripts = {}
    for trans in (make_transcript(record)
        for record in read_refgene(refgene_file)):
            transcripts[trans.name] = trans
            transcripts[trans.full_name] = trans

    return transcripts


def get_genomic_sequence(genome, chrom, start, end):
    """
    The function `get_genomic_sequence` returns a sequence for a given genomic region.
    
    :param genome: A dictionary containing genomic sequences for different chromosomes. The keys of the
    dictionary are chromosome names (e.g., 'chr1', 'chr2', etc.), and the values are the corresponding
    genomic sequences
    :param chrom: The chrom parameter represents the chromosome or genomic region from which you want to
    extract the sequence
    :param start: The start parameter is the 1-based coordinate of the beginning of the genomic region
    :param end: The `end` parameter is the end coordinate of the genomic region. It is a 1-based,
    end-inclusive coordinate, meaning that the base at the `end` position is included in the returned
    sequence
    :return: a sequence for the specified genomic region.
    """
    
    if start > end:
        return ''
    else:
        return str(genome[str(chrom)][start - 1:end]).upper()


# def get_allele(hgvs, genome, transcript=None):
#     """
#     The function "get_allele" takes a HGVSName, a genome, and an optional transcript, and returns the
#     chromosome, start and end positions, reference sequence, and alternate sequence of the allele.
    
#     :param hgvs: An object that provides methods for working with HGVS names
#     :param genome: The genome parameter is the genomic sequence from which the allele will be extracted
#     :param transcript: A transcript is a specific version of a gene that is transcribed into RNA. It
#     contains the coding sequence as well as the non-coding regions such as introns and untranslated
#     regions (UTRs). In this context, the transcript parameter is used to specify which transcript is
#     being referred to when retrieving the
#     :return: the chromosome, start position, end position, reference sequence, and alternate sequence
#     for a given HGVS name, genome, and transcript.
#     """
    
#     chrom, start, end = hgvs.get_ref_coords(transcript)
#     _, alt = hgvs.get_ref_alt(
#         transcript.tx_position.is_forward_strand if transcript else True)
#     ref = get_genomic_sequence(genome, chrom, start, end)
#     return chrom, start, end, ref, alt


_indel_mutation_types = set(['ins', 'del', 'dup', 'delins'])


def get_vcf_allele(hgvs, genome, transcript=None):
    """
    The function `get_vcf_allele` takes a HGVS name, a genome, and an optional transcript, and returns a
    VCF-style allele.
    
    :param hgvs: The `hgvs` parameter is an object of type `HGVSName`. It likely contains information
    about a genetic variant, such as the chromosome, start and end positions, and the type of mutation
    (e.g., substitution, deletion, insertion, etc.)
    :param genome: The `genome` parameter is the genomic sequence from which the allele will be
    extracted. It is a string representing the entire genome sequence
    :param transcript: The `transcript` parameter is an optional argument that represents a transcript.
    It is used to retrieve the VCF-style allele from the given HGVSName and genome. If a transcript is
    provided, the function will use it to get the VCF coordinates and the reference and alternate
    alleles. If no
    :return: the chromosome, start position, end position, reference allele, and alternate allele.
    """
    
    chrom, start, end = hgvs.get_vcf_coords(transcript)
    _, alt = hgvs.get_ref_alt(
        transcript.tx_position.is_forward_strand if transcript else True)
    ref = get_genomic_sequence(genome, chrom, start, end)

    # Sometimes we need to retrieve alt from reference
    # Eg NC_000001.11:g.169549811=
    if hgvs.mutation_type == "=":
        alt = ref

    if hgvs.mutation_type in _indel_mutation_types:
        if hgvs.mutation_type == 'dup':
            # No alt supplied: NM_000492.3:c.1155_1156dup
            # Number used:     NM_004119.2(FLT3):c.1794_1811dup18
            # We *know* what the sequence is for "dup18", but not for "ins18"
            if not hgvs.alt_allele or re.match("^N+$", hgvs.alt_allele):
                alt = get_alt_from_sequence(hgvs, genome, transcript)

        # Left-pad alternate allele.
        alt = ref[0] + alt
    return chrom, start, end, ref, alt


def get_alt_from_sequence(hgvs, genome, transcript):
    """
    The function "get_alt_from_sequence" returns a genomic sequence from a given HGVS notation, genome,
    and transcript.
    
    :param hgvs: The `hgvs` parameter is an object that provides methods for working with Human Genome
    Variation Society (HGVS) nomenclature. It likely has a method called `get_raw_coords()` that takes a
    transcript as input and returns the chromosome, start position, and end position of the
    corresponding genomic sequence
    :param genome: The genome parameter refers to the genomic sequence from which the alternative allele
    will be extracted
    :param transcript: The transcript parameter is a string that represents the transcript ID or name
    :return: the genomic sequence from the specified region in the genome.
    """

    chrom, start, end = hgvs.get_raw_coords(transcript)
    return get_genomic_sequence(genome, chrom, start, end)


def matches_ref_allele(hgvs, genome, transcript=None):
    """
    The function `matches_ref_allele` checks if the reference allele in a given HGVS notation matches
    the corresponding genomic sequence.
    
    :param hgvs: The `hgvs` parameter is an object that represents a variant in the Human Genome
    Variation Society (HGVS) format. It contains information about the variant's reference allele,
    alternative allele, and genomic coordinates
    :param genome: The `genome` parameter is the genomic sequence from which the reference allele is
    extracted
    :param transcript: The transcript parameter is an object that represents a transcript. It has a
    property called tx_position which provides information about the position of the transcript on the
    genome, including whether it is on the forward or reverse strand
    :return: True if the reference allele matches the genomic sequence, and False otherwise.
    """
    
    is_forward_strand = transcript.tx_position.is_forward_strand if transcript else True
    ref, _ = hgvs.get_ref_alt(is_forward_strand, raw_dup_alleles=True)  # get raw values so dup isn't always True
    chrom, start, end = hgvs.get_raw_coords(transcript)
    genome_ref = get_genomic_sequence(genome, chrom, start, end)
    return genome_ref == ref


def hgvs_justify_dup(chrom, offset, ref, alt, genome):
    """
    The function `hgvs_justify_dup` determines if an allele is a duplication and justifies it by
    returning the duplicated region if applicable.
    
    :param chrom: The chromosome name where the allele is located
    :param offset: The offset parameter is the 1-index genomic coordinate, which represents the position
    of the variant on the chromosome
    :param ref: The "ref" parameter represents the reference allele, which is the allele that is present
    in the reference genome at the given genomic coordinate
    :param alt: The `alt` parameter represents the alternate allele, which is the allele that differs
    from the reference allele at a specific genomic position
    :param genome: The `genome` parameter is a pygr compatible genome object. It is an object that
    represents a reference genome and provides methods to access genomic sequences
    :return: a tuple containing the chromosome name, offset, reference allele, alternate allele, and
    mutation type.
    """

    if len(ref) == len(alt) == 0:
        # it's a SNP, just return.
        return chrom, offset, ref, alt, '>'

    if len(ref) > 0 and len(alt) > 0:
        # complex indel, don't know how to dup check
        return chrom, offset, ref, alt, 'delins'

    if len(ref) > len(alt):
        # deletion -- don't dup check
        return chrom, offset, ref, alt, 'del'

    indel_seq = alt
    indel_length = len(indel_seq)

    # Convert offset to 0-index.
    offset -= 1

    # Get genomic sequence around the lesion.
    prev_seq = str(
        genome[str(chrom)][offset - indel_length:offset]).upper()
    next_seq = str(
        genome[str(chrom)][offset:offset + indel_length]).upper()

    # Convert offset back to 1-index.
    offset += 1

    if prev_seq == indel_seq:
        offset = offset - indel_length
        mutation_type = 'dup'
        ref = indel_seq
        alt = indel_seq * 2
    elif next_seq == indel_seq:
        mutation_type = 'dup'
        ref = indel_seq
        alt = indel_seq * 2
    else:
        mutation_type = 'ins'

    return chrom, offset, ref, alt, mutation_type


def hgvs_justify_indel(chrom, offset, ref, alt, strand, genome):
    """
    The function `hgvs_justify_indel` justifies an indel (insertion or deletion) according to the HGVS
    standard by determining the genomic sequence around the lesion, identifying the actual lesion
    sequence, and 3' justifying the offset.
    
    :param chrom: The chromosome where the indel is located
    :param offset: The offset parameter represents the position of the indel (insertion or deletion)
    within the chromosome or genomic sequence
    :param ref: The `ref` parameter represents the reference allele of the variant. It is a string that
    contains the nucleotide sequence of the reference allele
    :param alt: The `alt` parameter in the `hgvs_justify_indel` function represents the alternate allele
    sequence for an indel variant. It is the sequence that replaces the reference allele sequence
    (`ref`) at the specified `offset` position on the `chrom` chromosome
    :param strand: The parameter "strand" represents the orientation of the DNA strand where the indel
    is located. It can have two possible values: "+" or "-". The "+" strand refers to the forward
    strand, while the "-" strand refers to the reverse complement strand
    :param genome: The `genome` parameter is a dictionary that contains the genomic sequence for each
    chromosome. The keys of the dictionary are the chromosome names (e.g., "chr1", "chr2", etc.), and
    the values are the corresponding genomic sequences
    :return: the variables chrom, offset, ref, and alt.
    """

    if len(ref) == len(alt) == 0:
        # It's a SNP, just return.
        return chrom, offset, ref, alt

    if len(ref) > 0 and len(alt) > 0:
        # Complex indel, don't know how to justify.
        return chrom, offset, ref, alt

    # Get genomic sequence around the lesion.
    window_size = 100
    size = window_size + max(len(ref), len(alt))
    start = max(offset - size, 0)
    end = offset + size
    seq = str(genome[str(chrom)][start - 1:end]).upper()
    cds_offset = offset - start

    # indel -- strip off the ref base to get the actual lesion sequence
    is_insert = len(alt) > 0
    if is_insert:
        indel_seq = alt
        cds_offset_end = cds_offset
    else:
        indel_seq = ref
        cds_offset_end = cds_offset + len(indel_seq)

    # Now 3' justify (vs. cDNA not genome) the offset
    justify = 'right' if strand == '+' else 'left'
    offset, _, indel_seq = justify_indel(
        cds_offset, cds_offset_end, indel_seq, seq, justify)
    offset += start

    if is_insert:
        alt = indel_seq
    else:
        ref = indel_seq

    return chrom, offset, ref, alt


def hgvs_normalize_variant(chrom, offset, ref, alt, genome, transcript=None):
    """
    The function `hgvs_normalize_variant` converts a variant in VCF-style to HGVS-style by adjusting the
    offset, reference and alternate alleles, and determining the mutation type.
    
    :param chrom: The chromosome where the variant is located
    :param offset: The offset parameter represents the position of the variant within the chromosome. It
    is an integer value
    :param ref: The `ref` parameter represents the reference allele in a variant
    :param alt: The `alt` parameter represents the alternate allele in a variant. It is a string that
    represents the alternative nucleotide(s) or sequence(s) at a specific position in the genome
    :param genome: The `genome` parameter is the reference genome sequence. It is used to perform
    certain operations on the variant, such as justifying indels and representing duplications
    :param transcript: The `transcript` parameter is an optional argument that represents the transcript
    or gene in which the variant occurs. It is used to determine the strand of the gene and to perform
    certain operations on the variant. If no transcript is provided, the default value is `None`
    :return: the following values: chrom, offset, ref, alt, and mutation_type.
    """

    if len(ref) == len(alt) == 1:
        if ref == alt:
            mutation_type = '='
        else:
            mutation_type = '>'
    else:
        # Remove 1bp padding
        offset += 1
        ref = ref[1:]
        alt = alt[1:]

        # 3' justify allele.
        strand = transcript.strand if transcript else '+'
        chrom, offset, ref, alt = hgvs_justify_indel(
            chrom, offset, ref, alt, strand, genome)

        # Represent as duplication if possible.
        chrom, offset, ref, alt, mutation_type = hgvs_justify_dup(
            chrom, offset, ref, alt, genome)
    return chrom, offset, ref, alt, mutation_type


def parse_hgvs_name(hgvs_name, genome, transcript=None, get_transcript=lambda name: None, flank_length=30, normalize=True, lazy=False, indels_start_with_same_base=True):
    """
    The function `parse_hgvs_name` takes an HGVS name, a genome object, and optional parameters, and
    returns the chromosome, start position, reference allele, and alternate allele of the variant
    described by the HGVS name.
    
    :param hgvs_name: The HGVS name to parse
    :param genome: A pygr compatible genome object. This object represents the reference genome and
    provides methods to access genomic sequences and annotations
    :param transcript: The transcript parameter is an optional argument that represents the transcript
    corresponding to the HGVS name. It is used to determine the reference sequence for the variant. If
    not provided, the get_transcript function is used to retrieve the transcript based on the HGVS name.
    If neither transcript nor get_transcript is
    :param get_transcript: A function that takes a transcript name as input and returns the
    corresponding transcript object. If not provided, the default behavior is to return None
    :param flank_length: The `flank_length` parameter is an integer that specifies the length of the
    flanking sequence to include when normalizing the variant allele. This is used in the
    `normalize_variant` function to determine the reference allele and normalize the variant allele
    according to the VCF standard, defaults to 30 (optional)
    :param normalize: A boolean parameter that determines whether the allele should be normalized
    according to the VCF standard. If set to True, the allele will be normalized; if set to False, the
    allele will not be normalized, defaults to True (optional)
    :param lazy: The `lazy` parameter is a boolean flag that determines whether or not to discard
    version information from the incoming transcript or gene. If `lazy` is set to `True`, the version
    information will be discarded. If `lazy` is set to `False`, the version information will be included
    in the, defaults to False (optional)
    :param indels_start_with_same_base: The parameter "indels_start_with_same_base" is a boolean flag
    that determines whether or not to strip the common prefix from indels when normalizing alleles. If
    set to True, the common prefix will not be stripped, defaults to True (optional)
    :return: The function `parse_hgvs_name` returns a tuple containing the chromosome, start position,
    reference allele, and alternate allele of the parsed HGVS name.
    """

    hgvs = HGVSName(hgvs_name)

    # Determine transcript.
    if hgvs.kind in ('c', 'n') and not transcript:
        if '.' in hgvs.transcript and lazy:
            hgvs.transcript, _ = hgvs.transcript.split('.')
        elif '.' in hgvs.gene and lazy:
            hgvs.gene, _ = hgvs.gene.split('.')
        if get_transcript:
            if hgvs.transcript:
                transcript = get_transcript(hgvs.transcript)
            elif hgvs.gene:
                transcript = get_transcript(hgvs.gene)
        if not transcript:
            raise ValueError('transcript is required')

    if transcript and hgvs.transcript in genome:
        # Reference sequence is directly known, use it.
        genome = GenomeSubset(genome, transcript.tx_position.chrom,
                              transcript.tx_position.chrom_start,
                              transcript.tx_position.chrom_stop,
                              hgvs.transcript)

    chrom, start, _, ref, alt = get_vcf_allele(hgvs, genome, transcript)
    if normalize:
        nv = normalize_variant(chrom, start, ref, [alt], genome,
                               flank_length=flank_length,
                               indels_start_with_same_base=indels_start_with_same_base)
        chrom, start, ref, [alt] = nv.variant
    return (chrom, start, ref, alt)



def cdna_to_protein(hgvs, offset, genome, chrom, transcript, ref, alt, mutation_type, codon_type:str = "3"):
    """
    The function `cdna_to_protein` takes in various parameters related to a genetic mutation and returns
    an updated HGVS object with additional protein information.
    
    :param hgvs: The parameter `hgvs` is an object that represents a variant in the Human Genome
    Variation Society (HGVS) format. It contains information about the variant, such as the cDNA start
    and end positions
    :param offset: The offset is a numerical value that represents the starting position of the genomic
    sequence in the reference genome. It is used to calculate the genomic position of the mutation
    :param genome: The `genome` parameter is a dictionary that represents the genomic sequence. It
    contains the chromosome as the key and the corresponding DNA sequence as the value
    :param chrom: The `chrom` parameter represents the chromosome on which the mutation occurs
    :param transcript: The `transcript` parameter is a string that represents the transcript ID or name.
    It is used to identify the specific transcript in the genome
    :param ref: The parameter "ref" is a string that represents the reference nucleotide sequence. It is
    used to determine the codons in the DNA sequence
    :param alt: The `alt` parameter in the `cdna_to_protein` function is a string that represents the
    alternate nucleotide sequence for a mutation
    :param mutation_type: The `mutation_type` parameter is a string that represents the type of
    mutation. It can have the following values:
    :param codon_type: The `codon_type` parameter is a string that specifies the type of codon
    translation to be used. It can have one of the following values:, defaults to 3
    :type codon_type: str (optional)
    :return: the updated `hgvs` object.
    """

    if hgvs.cdna_start.offset == 0 and not hgvs.pep_extra:

        # InDel ? MNV ?
        is_indel = mutation_type not in [">"]
        is_mnv = mutation_type in ["delins"] and len(ref) == len(alt)

        # gap
        gap_offset = -1
        gap_cdna = -1
        gap_cdna_hgvs = +1
        if mutation_type in [">"]:
            gap_offset = -1
            gap_cdna = -1
            gap_cdna_hgvs = +1
        if mutation_type in ["ins"]:
            gap_offset = -1
            gap_cdna = 0
            gap_cdna_hgvs = +1
        if mutation_type in ["del"]:
            gap_offset = -1
            gap_cdna = -1
            gap_cdna_hgvs = +1

        genomic_position = offset + gap_offset
        cdna_position_start = hgvs.cdna_start.coord + gap_cdna
        cdna_position_end = hgvs.cdna_end.coord + gap_cdna

        # offsets / positions
        offset_protein_mod = int(int(cdna_position_start) % 3)
        offset_protein = int(int(cdna_position_start) / 3)
        offset_protein_hgvs = offset_protein + gap_cdna_hgvs
        offset_protein_end = int(int(cdna_position_end) / 3)
        offset_genomic_codon_start = genomic_position - offset_protein_mod #- gap_ins
        offset_genomic_codon_end = offset_genomic_codon_start + 3 + ( (offset_protein_end - offset_protein) * 3 )

        # ref sequence
        seq_ref = str(genome[str(chrom)][offset_genomic_codon_start:offset_genomic_codon_end])
        
        # alt sequence
        if is_indel and not is_mnv:
            seq_alt = ""
        else:
            seq_alt_split = []
            seq_alt_split.extend(seq_ref)
            alt_split = []
            alt_split.extend(alt)
            i = 0
            while i < len(ref):
                nuleotide = alt_split[i].upper()
                seq_alt_split[offset_protein_mod+i] = nuleotide
                i += 1
            seq_alt = "".join(seq_alt_split)

        # transcript strand
        if not transcript.tx_position.is_forward_strand:
            seq_ref = revcomp(seq_ref)
            seq_alt = revcomp(seq_alt)

        # Split codons
        seq_ref_split_codon = re.findall('...',seq_ref)
        seq_alt_split_codon = re.findall('...',seq_alt)
        if not is_mnv:
            if len(seq_ref_split_codon):
                seq_ref_split_codon = [seq_ref_split_codon[0]]
                if len(seq_alt_split_codon):
                    seq_alt_split_codon = [seq_alt_split_codon[0]]

        # codon type
        if codon_type == "1":
            codon_translate = CODON_1
        elif codon_type == "3":
            codon_translate = CODON_3
        elif codon_type == "FULL":
            codon_translate = CODON_FULL
        else:
            codon_translate = CODON_3

        # codons
        codon3_ref = ""
        codon3_alt = ""
        for codon in seq_ref_split_codon:
            codon3_ref += str(codon_translate.get(codon.upper()))
        for codon in seq_alt_split_codon:
            codon3_alt += str(codon_translate.get(codon.upper()))

        # indel alt
        if is_indel and not is_mnv:
            codon3_alt = "fs"
        
        # add protein infos
        hgvs.pep_extra = f"{codon3_ref}{offset_protein_hgvs}{codon3_alt}"

    # return
    return hgvs


def variant_to_hgvs_name(chrom, offset, ref, alt, genome, transcript, transcript_protein=None, exon=None, max_allele_length=4, use_counsyl=False, codon_type:str = "3"):
    """
    The function `variant_to_hgvs_name` takes in genomic coordinates, alleles, and other parameters, and
    returns a HGVS-style name for the variant.
    
    :param chrom: The chromosome name where the variant is located
    :param offset: The `offset` parameter represents the genomic offset of the allele. It is the
    position of the variant on the chromosome
    :param ref: The reference allele at the given genomic coordinate
    :param alt: The `alt` parameter is the alternate allele. In genetics, a variant or mutation can
    occur at a specific position in the genome, and the `alt` allele represents the alternative
    nucleotide or sequence at that position compared to the reference genome
    :param genome: A pygr compatible genome object, which represents the reference genome sequence
    :param transcript: The `transcript` parameter is the transcript corresponding to the allele. It is
    used to determine the type of coordinates to use in the HGVS name (either genomic coordinates or
    cDNA coordinates). If the transcript is not available, the function will use genomic coordinates
    :param transcript_protein: The `transcript_protein` parameter is an optional argument that
    represents the protein sequence corresponding to the transcript. It is used to populate the
    `transcript_protein` attribute of the `HGVSName` object
    :param exon: The `exon` parameter is an optional argument that represents the exon number or
    identifier associated with the variant. It is used to populate the `exon` attribute of the
    `HGVSName` object. If provided, it will be included in the final HGVS name generated by the function
    :param max_allele_length: The `max_allele_length` parameter is used to determine whether to
    represent the alleles as their actual sequence or as the length of the sequence. If the length of
    the reference allele or alternate allele is greater than `max_allele_length`, then the length of the
    allele is used instead of the actual, defaults to 4 (optional)
    :param use_counsyl: A boolean flag indicating whether to use Counsyl-specific rules for single-base
    indels, defaults to False (optional)
    :param codon_type: The parameter `codon_type` is a string that specifies the type of codon numbering
    to be used in the HGVS name. It is used in the `cdna_to_protein` function to determine the type of
    codon numbering to be used in the protein-level HGVS name. The, defaults to 3
    :type codon_type: str (optional)
    :return: an object of type HGVSName.
    """

    # Convert VCF-style variant to HGVS-style.
    chrom, offset, ref, [alt] = normalize_variant(
        chrom, offset, ref, [alt], genome).variant
    chrom, offset, ref, alt, mutation_type = hgvs_normalize_variant(
        chrom, offset, ref, alt, genome, transcript)

    # Populate HGVSName parse tree.
    hgvs = HGVSName()

    # Populate coordinates.
    if mutation_type == 'ins':
        # Insert uses coordinates around the insert point.
        offset_start = offset - 1
        offset_end = offset
    else:
        offset_start = offset
        offset_end = offset + len(ref) - 1

    if not transcript:
        # Use genomic coordinate when no transcript is available.
        hgvs.kind = 'g'
        hgvs.start = offset_start
        hgvs.end = offset_end
    else:
        # Use cDNA coordinates.
        if transcript.is_coding:
            hgvs.kind = 'c'
        else:
            hgvs.kind = 'n'

        is_single_base_indel = (
            (mutation_type == 'ins' and len(alt) == 1) or
            (mutation_type in ('del', 'delins', 'dup') and len(ref) == 1))

        if mutation_type == '>' or (use_counsyl and is_single_base_indel):
            # Use a single coordinate.
            hgvs.cdna_start = transcript.genomic_to_cdna_coord(offset)
            hgvs.cdna_end = hgvs.cdna_start
        else:
            if transcript.strand == '-':
                offset_start, offset_end = offset_end, offset_start
            hgvs.cdna_start = transcript.genomic_to_cdna_coord(offset_start)
            hgvs.cdna_end = transcript.genomic_to_cdna_coord(offset_end)

        # pep extra
        hgvs = cdna_to_protein(hgvs, offset, genome, chrom, transcript, ref, alt, mutation_type, codon_type=codon_type)
        

    # Populate prefix.
    if transcript:
        hgvs.transcript = transcript.full_name
        hgvs.gene = transcript.gene.name

    # Add transcript protein
    if transcript_protein:
        hgvs.transcript_protein = transcript_protein

    # Add transcript protein
    if exon:
        hgvs.exon = exon

    # Convert alleles to transcript strand.
    if transcript and transcript.strand == '-':
        ref = revcomp(ref)
        alt = revcomp(alt)

    # Convert to allele length if alleles are long.
    ref_len = len(ref)
    alt_len = len(alt)
    if ((mutation_type == 'dup' and ref_len > max_allele_length) or
            (mutation_type != 'dup' and
             (ref_len > max_allele_length or alt_len > max_allele_length))):
        ref = str(ref_len)
        alt = str(alt_len)

    # Populate alleles.
    hgvs.mutation_type = mutation_type
    hgvs.ref_allele = ref
    hgvs.alt_allele = alt

    return hgvs


def format_hgvs_name(chrom, offset, ref, alt, genome, transcript, transcript_protein=None, exon=None, use_prefix=True, use_gene=True, use_protein=False, use_counsyl=False, max_allele_length=4, full_format=False, use_version=False, codon_type:str = "3"):
    """
    The `format_hgvs_name` function generates a HGVS name from a genomic coordinate.
    
    :param chrom: The `chrom` parameter represents the chromosome name. It is a string that specifies
    the chromosome on which the variant occurs
    :param offset: The `offset` parameter represents the genomic offset of the allele, which is the
    position of the variant on the chromosome. It is used to generate the HGVS name based on the genomic
    coordinate
    :param ref: The `ref` parameter represents the reference allele. In genetics, a variant or mutation
    can occur at a specific position in the genome, resulting in a change from the reference allele to
    an alternate allele. The `ref` parameter specifies the sequence of the reference allele at that
    position
    :param alt: The `alt` parameter represents the alternate allele. In genetics, a variant or mutation
    can occur at a specific position in the genome, resulting in a change from the reference allele to
    an alternate allele. The `alt` parameter specifies the sequence of the alternate allele at that
    position
    :param genome: A pygr compatible genome object, which is used to retrieve genomic sequences and
    annotations. It provides methods to access genomic information such as chromosome names, sequences,
    and gene annotations
    :param transcript: The `transcript` parameter is the transcript corresponding to the allele. It is
    used to generate the HGVS name based on the genomic coordinate
    :param transcript_protein: The `transcript_protein` parameter is an optional argument that
    represents the protein transcript corresponding to the cDNA transcript. It is used to generate the
    protein HGVS name if it exists
    :param exon: The `exon` parameter is used to specify the exon number in the HGVS name. It is an
    optional parameter and is used to generate a more specific HGVS name when needed
    :param use_prefix: A boolean indicating whether to include a transcript/gene/chromosome prefix in
    the HGVS name. If set to True, the prefix will be included; if set to False, the prefix will be
    excluded, defaults to True (optional)
    :param use_gene: A boolean parameter that determines whether to include the gene name in the HGVS
    prefix. If set to True, the gene name will be included; if set to False, the gene name will be
    excluded, defaults to True (optional)
    :param use_protein: A boolean parameter that determines whether to include protein HGVS notation in
    the generated HGVS name. If set to True, the protein HGVS notation will be included if it exists. If
    set to False, only the genomic and transcript HGVS notation will be included, defaults to False
    (optional)
    :param use_counsyl: The `use_counsyl` parameter is a boolean parameter that determines whether to
    use Counsyl-specific formatting for the HGVS name. If set to True, the HGVS name will be formatted
    according to Counsyl's specific guidelines. If set to False, the HGVS name will be, defaults to
    False (optional)
    :param max_allele_length: The `max_allele_length` parameter is used to determine the maximum length
    of the allele. If the length of the allele is greater than the specified `max_allele_length`, then
    the allele length will be used in the HGVS name instead of the actual allele sequence. By default,
    the `, defaults to 4 (optional)
    :param full_format: A boolean parameter that determines whether to use the full HGVS format or not.
    If set to True, the HGVS name will include the gene name, transcript name, exon number (if
    provided), and the amino acid change (if protein information is available). If set to False, the
    HGVS, defaults to False (optional)
    :param use_version: A boolean parameter that determines whether to include the version number of the
    transcript in the HGVS name. If set to True, the version number will be included; if set to False,
    the version number will be excluded, defaults to False (optional)
    :param codon_type: The `codon_type` parameter is a string that specifies the type of codon numbering
    to be used in the HGVS name. It can have one of the following values:, defaults to 3
    :type codon_type: str (optional)
    :return: a formatted HGVS name generated from a genomic coordinate.
    """
    
    hgvs = variant_to_hgvs_name(chrom, offset, ref, alt, genome, transcript, transcript_protein=transcript_protein, exon=exon,
                                max_allele_length=max_allele_length,
                                use_counsyl=use_counsyl, codon_type=codon_type)
    return hgvs.format(use_prefix=use_prefix, use_gene=use_gene, use_exon=exon, use_protein=use_protein, full_format=full_format, use_version=use_version)


def create_refseq_table(conn, refseq_table:str = "refseq", refseq_file:str = None) -> str:
    """
    The function `create_refseq_table` creates a table in a database with the specified name and
    structure, either using a file or without a file.
    
    :param conn: The `conn` parameter is a connection object that represents a connection to a database.
    It is used to execute SQL queries and interact with the database
    :param refseq_table: The `refseq_table` parameter is a string that specifies the name of the table
    that will be created in the database to store the RefGene data, defaults to refseq
    :type refseq_table: str (optional)
    :param refseq_file: The `refseq_file` parameter is a string that specifies the path to a file
    containing the data for the refGene table. If this parameter is provided, the function will create
    the refGene table in the database using the data from the file. If this parameter is not provided,
    the function will
    :type refseq_file: str
    :return: the name of the refseq table that was created or used.
    """

    # RefGene in database
    if refseq_table == "refseq":
        refseq_structure = {
            "bin": "INTEGER",
            "name": "STRING",
            "chrom": "STRING",
            "strand": "STRING",
            "txStart": "INTEGER",
            "txEnd": "INTEGER",
            "cdsStart": "INTEGER",
            "cdsEnd": "INTEGER",
            "exonCount": "INTEGER",
            "exonStarts": "STRING",
            "exonEnds": "STRING",
            "score": "INTEGER",
            "name2": "STRING",
            "cdsStartStat": "STRING",
            "cdsEndStat": "STRING",
            "exonFrames": "STRING"
        }
        query = f"CREATE TABLE {refseq_table} AS SELECT * FROM read_csv_auto('{refseq_file}',HEADER=False,columns={refseq_structure},skip=1)"
    elif refseq_table == "refseqlink":
        refseq_structure = {
            "id": "STRING",
            "status": "STRING",
            "name": "STRING",
            "product": "STRING",
            "mrnaAcc": "STRING",
            "protAcc": "STRING",
            "locusLinkId": "STRING",
            "omimId": "STRING",
            "hgnc": "STRING",
            "genbank" :"STRING",
            "pseudo": "STRING",
            "gbkey": "STRING",
            "source": "STRING",
            "gene_biotype": "STRING",
            "gene_synonym": "STRING",
            "ncrna_class": "STRING",
            "note": "STRING",
            "description": "STRING",
            "externalId": "STRING",
            }
        query = rf"CREATE TABLE {refseq_table} AS SELECT *, regexp_extract(mrnaAcc, '(.*)\..*', 1) AS 'mrnaAcc_without_ver', regexp_extract(protAcc, '(.*)\..*', 1) AS 'protAcc_without_ver', mrnaAcc AS 'mrnaAcc_with_ver', protAcc AS 'protAcc_with_ver' FROM read_csv_auto('{refseq_file}',HEADER=False,columns={refseq_structure},skip=1)"

    # Create tabel with file
    if refseq_file:
        #print(f"Create table refGene '{refseq_table}' with file '{refseq_file}'")
        conn.query(query)
    # Create table without file
    else:
        #print(f"Create table refGene '{refseq_table}' without file")
        sql_structure_list = []
        for col in refseq_structure:
            col_format = refseq_structure.get(col,"VARCHAR").replace("STRING","VARCHAR")
            sql_structure_list.append(f" {col} {col_format}")
        sql_structure = ",".join(sql_structure_list)
        conn.query(f"CREATE TABLE {refseq_table}({sql_structure})")

    return refseq_table


def get_refseq_table(conn, refseq_table:str = "refseq", refseq_file:str = None) -> str:
    """
    The function `get_refseq_table` checks if a table named `refseq` exists in a database, and if not,
    creates it using the `create_refseq_table` function.
    
    :param conn: The parameter `conn` is expected to be a connection object that allows you to interact
    with a database. It could be an instance of a database connector class, such as `pymysql.connect()`
    for MySQL or `psycopg2.connect()` for PostgreSQL
    :param refseq_table: The parameter "refseq_table" is a string that specifies the name of the table
    in the database where the refGene data will be stored. If this table already exists in the database,
    the function will return the name of the existing table. If the table does not exist, the function
    will create, defaults to refseq
    :type refseq_table: str (optional)
    :param refseq_file: The `refseq_file` parameter is the name or path of the file that contains the
    refGene data. This file is used to populate the refGene table in the database
    :type refseq_file: str
    :return: the name of the refseq_table.
    """

    # Existing tables
    existing_tables = conn.query(f"SHOW TABLES").df()

    # refGene already table exists
    if refseq_table in list(existing_tables["name"]):
        return refseq_table
    # Create refGene table
    else:
        refseq_table = create_refseq_table(conn, refseq_table=refseq_table, refseq_file=refseq_file)
        return refseq_table


def get_transcript(transcripts:dict, transcript_name:str) -> Transcript:
    """
    The function `get_transcript` takes a dictionary of transcripts and a name as input, and returns the
    transcript associated with that name.
    
    :param transcripts: A dictionary containing transcripts as values, with names as keys
    :type transcripts: dict
    :param name: The name parameter is a string that represents the name of the transcript that you want
    to retrieve from the transcripts dictionary
    :type name: str
    :return: the value associated with the given name key in the transcripts dictionary.
    """

    return transcripts.get(transcript_name)

