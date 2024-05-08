"""
Models for representing genomic elements.
"""

from __future__ import unicode_literals

from collections import namedtuple

from lazy import lazy

from howard.objects.cdna import CDNA_START_CODON, CDNA_STOP_CODON, CDNACoord


# The Gene class is a basic template for representing genes in Python.
class Gene(object):
    def __init__(self, name):
        self.name = name


# The Transcript class is a blueprint for creating objects that represent a transcript.
# A gene may have multiple transcripts with different combinations of exons.
# We need both exons and cdna_match as need to know exact exon boundaries to work out flanking
class Transcript(object):

    def __init__(
        self,
        name: str,
        version: str,
        gene: str,
        tx_position: int,
        cds_position: int,
        is_default: bool = False,
        cdna_match: list = None,
        start_codon_transcript_pos: int = None,
        stop_codon_transcript_pos: int = None,
    ) -> None:
        """
        The function initializes an object with various attributes related to a gene and its transcript.

        :param name: A string representing the name of the coding sequence
        :type name: str
        :param version: The `version` parameter is a string that represents the version of the object.
        It is used to track changes or updates to the object over time
        :type version: str
        :param gene: The `gene` parameter is a string that represents the gene associated with the
        coding sequence
        :type gene: str
        :param tx_position: The `tx_position` parameter represents the position of the transcript. It is
        an integer value that indicates the position of the transcript in the genome
        :type tx_position: int
        :param cds_position: The `cds_position` parameter represents the position of the coding sequence
        (CDS) within the transcript. It is an integer value that indicates the starting position of the
        CDS within the transcript sequence
        :type cds_position: int
        :param is_default: The `is_default` parameter is a boolean flag that indicates whether the
        instance of the class is the default version of the gene. It is set to `False` by default, but
        can be set to `True` if the instance is the default version, defaults to False
        :type is_default: bool (optional)
        :param cdna_match: The `cdna_match` parameter is a list that contains the positions of the
        matching cDNA sequences. It is an optional parameter and if not provided, it defaults to an
        empty list
        :type cdna_match: list
        :param start_codon_transcript_pos: The parameter "start_codon_transcript_pos" is an optional
        parameter that represents the transcript position of the start codon. It is used to store the
        pre-calculated transcript position of the start codon for a specific gene
        :type start_codon_transcript_pos: int
        :param stop_codon_transcript_pos: The parameter `stop_codon_transcript_pos` is an optional
        integer that represents the transcript position of the stop codon. It is used to store the
        pre-calculated transcript coordinate of the stop codon. If not provided, it will be set to
        `None`
        :type stop_codon_transcript_pos: int
        """

        self.name = name
        self.version = version
        self.gene = Gene(gene)
        self.tx_position = tx_position
        self.cds_position = cds_position
        self.is_default = is_default
        self.cdna_match = cdna_match or []
        # Optional pre-calculated transcript coordinates
        self._start_codon_transcript_pos = start_codon_transcript_pos
        self._stop_codon_transcript_pos = stop_codon_transcript_pos

    @property
    def full_name(self) -> str:
        """
        The function `full_name` returns the full name of an object, including its version if it exists.
        :return: a string. If the `version` attribute of the object is not `None`, it returns a string
        in the format `name.version`. Otherwise, it returns just the `name` attribute.
        """

        if self.version is not None:
            return "%s.%d" % (self.name, self.version)
        else:
            return self.name

    @property
    def is_coding(self) -> bool:
        """
        The function checks if a coding transcript has a non-zero length coding sequence.
        :return: a boolean value indicating whether the coding transcript has a coding sequence (CDS)
        with a non-zero length.
        """

        return self.cds_position.chrom_stop - self.cds_position.chrom_start > 0

    @property
    def strand(self) -> str:
        """
        The function returns a string '+' if the tx_position is on the forward strand, and '-' if it is
        on the reverse strand.
        :return: a string that represents the strand of the given `self.tx_position`. If
        `self.tx_position.is_forward_strand` is `True`, then the string returned is '+'. Otherwise, the
        string returned is '-'.
        """

        return "+" if self.tx_position.is_forward_strand else "-"

    @lazy
    def ordered_cdna_match(self) -> str:
        """
        The function "ordered_cdna_match" sorts a list of cdna_match objects based on their
        tx_position.chrom_start attribute and returns the sorted list.
        :return: a sorted list of cdna_match objects.
        """

        transcript_strand = self.tx_position.is_forward_strand
        cdna_match = list(self.cdna_match)
        cdna_match.sort(key=lambda cm: cm.tx_position.chrom_start)
        if not transcript_strand:
            cdna_match.reverse()
        return cdna_match

    def get_cds_start_stop(self) -> tuple:
        """
        The function "get_cds_start_stop" returns the start and stop positions of a coding sequence,
        taking into account the direction of the strand.
        :return: a tuple containing the start and stop positions of the coding sequence (CDS).
        """

        (start, stop) = (self.cds_position.chrom_start, self.cds_position.chrom_stop)
        if not self.tx_position.is_forward_strand:
            (stop, start) = (start, stop)
        return start, stop

    @lazy
    def start_codon(self) -> int:
        """
        The function returns the transcript position of the start codon.
        :return: the transcript position of the start codon.
        """

        if self._start_codon_transcript_pos is not None:
            return self._start_codon_transcript_pos

        start_genomic_coordinate, _ = self.get_cds_start_stop()
        return self._get_transcript_position(start_genomic_coordinate)

    @lazy
    def stop_codon(self) -> int:
        """
        The function returns the transcript position of the stop codon.
        :return: The method `stop_codon` returns an integer, which represents the transcript position of
        the stop codon.
        """

        if self._stop_codon_transcript_pos is not None:
            return self._stop_codon_transcript_pos

        _, stop_genomic_coordinate = self.get_cds_start_stop()
        return self._get_transcript_position(stop_genomic_coordinate)

    def _get_transcript_position(
        self, genomic_coordinate: int, label: str = None
    ) -> int:
        """
        The function `_get_transcript_position` calculates the position of a genomic coordinate within a
        transcript.

        :param genomic_coordinate: The `genomic_coordinate` parameter is an integer representing a
        specific position on the genome. It is the coordinate for which we want to determine the
        corresponding position in the transcript
        :type genomic_coordinate: int
        :param label: The `label` parameter is an optional string that represents a label or identifier
        for the genomic coordinate. It is used in the error message if the genomic coordinate is not
        found in any of the exons. If no label is provided, a default label is generated using the value
        of the `genomic
        :type label: str
        :return: an integer representing the position of the genomic coordinate within the transcript.
        """

        transcript_strand = self.tx_position.is_forward_strand

        cdna_offset = 0
        for cdna_match in self.ordered_cdna_match:
            start = cdna_match.tx_position.chrom_start
            stop = cdna_match.tx_position.chrom_stop

            if start <= genomic_coordinate <= stop:
                # We're inside this match
                if transcript_strand:
                    position = genomic_coordinate - start
                else:
                    position = stop - genomic_coordinate
                return cdna_offset + position + cdna_match.get_offset(position)
            else:
                cdna_offset += cdna_match.length
        if label is None:
            label = "Genomic coordinate: %d" % genomic_coordinate
        raise ValueError("%s is not in any of the exons" % label)

    def cdna_to_genomic_coord(self, coord: object) -> int:
        """
        The function `cdna_to_genomic_coord` converts a HGVS cDNA coordinate to a genomic coordinate.

        :param coord: The parameter `coord` is an object that represents a cDNA coordinate. It is used
        to specify a position along a cDNA sequence
        :type coord: object
        :return: an integer value, which represents the genomic coordinate corresponding to the given
        cDNA coordinate.
        """

        transcript_strand = self.tx_position.is_forward_strand

        # compute starting position along spliced transcript.
        if coord.landmark == CDNA_START_CODON:
            utr5p = self.start_codon if self.is_coding else 0

            if coord.coord > 0:
                cdna_pos = utr5p + coord.coord
            else:
                cdna_pos = utr5p + coord.coord + 1
        elif coord.landmark == CDNA_STOP_CODON:
            if coord.coord < 0:
                raise ValueError(
                    "CDNACoord cannot have a negative coord and "
                    "landmark CDNA_STOP_CODON"
                )
            cdna_pos = self.stop_codon + coord.coord
        else:
            raise ValueError('unknown CDNACoord landmark "%s"' % coord.landmark)

        # 5' flanking sequence (no need to account for gaps)
        if cdna_pos < 1:
            if transcript_strand:
                return self.tx_position.chrom_start + cdna_pos
            else:
                return self.tx_position.chrom_stop - cdna_pos + 1

        # Walk along transcript until we find an exon that contains cdna_pos.
        for cdna_match in self.ordered_cdna_match:
            if cdna_match.cdna_start <= cdna_pos <= cdna_match.cdna_end:
                match_pos = cdna_pos - cdna_match.cdna_start
                match_pos -= cdna_match.get_offset(match_pos)
                # Compute genomic coordinate using offset.
                if transcript_strand:
                    # Plus strand.
                    start = cdna_match.tx_position.chrom_start + 1
                    return start + match_pos + coord.offset
                else:
                    # Minus strand.
                    end = cdna_match.tx_position.chrom_stop
                    return end - match_pos - coord.offset
        else:
            # 3' flanking sequence (no need to account for gaps)
            if transcript_strand:
                return self.cds_position.chrom_stop + coord.coord
            else:
                return self.cds_position.chrom_start + 1 - coord.coord

    def genomic_to_cdna_coord(self, genomic_coord: int) -> object:
        """
        The function `genomic_to_cdna_coord` converts a genomic coordinate to a cDNA coordinate and
        offset, taking into account exons, strand, and coding transcript information.

        :param genomic_coord: The `genomic_coord` parameter is an integer representing a genomic
        coordinate
        :type genomic_coord: int
        :return: an object of type `CDNACoord`.
        """

        exons = [exon.get_as_interval() for exon in self.ordered_cdna_match]

        if len(exons) == 0:
            return None

        strand = self.strand
        distances = [exon.distance(genomic_coord) for exon in exons]
        min_distance_to_exon = min(map(abs, distances))
        if min_distance_to_exon:
            # We're outside of exon - so need to find closest point
            for exon in exons:
                distance = exon.distance(genomic_coord)
                if abs(distance) == min_distance_to_exon:

                    # Outside the exon.
                    if distance > 0:
                        genomic_nearest_exon = exon.chrom_start + 1
                    else:
                        genomic_nearest_exon = exon.chrom_end

                    if strand == "+":
                        distance *= -1

                    nearest_exon_coord = self._exon_genomic_to_cdna_coord(
                        genomic_nearest_exon
                    )

                    # If outside transcript, don't use offset.
                    if (
                        genomic_coord < self.tx_position.chrom_start + 1
                        or genomic_coord > self.tx_position.chrom_stop
                    ):
                        nearest_exon_coord += distance
                        distance = 0
                    cdna_coord = CDNACoord(nearest_exon_coord, distance)
                    break
            else:
                raise ValueError("Could not find closest exon!")  # Should never happen
        else:
            coord = self._exon_genomic_to_cdna_coord(genomic_coord)
            cdna_coord = CDNACoord(coord, 0)

        # Adjust coordinates for coding transcript.
        if self.is_coding:
            # Detect if position before start codon.
            utr5p = self.start_codon
            cdna_coord.coord -= utr5p
            if cdna_coord.coord <= 0:
                cdna_coord.coord -= 1
            else:
                # Detect if position is after stop_codon.
                stop_codon = self.stop_codon
                stop_codon -= utr5p
                if (
                    cdna_coord.coord > stop_codon
                    or cdna_coord.coord == stop_codon
                    and cdna_coord.offset > 0
                ):
                    cdna_coord.coord -= stop_codon
                    cdna_coord.landmark = CDNA_STOP_CODON

        return cdna_coord

    def _exon_genomic_to_cdna_coord(self, genomic_coord: int) -> int:
        """
        The function `_exon_genomic_to_cdna_coord` converts a genomic coordinate to a cDNA coordinate
        within an exon.

        :param genomic_coord: The `genomic_coord` parameter represents the genomic coordinate that you
        want to convert to cDNA coordinate. It is an integer value that specifies the position on the
        genome
        :type genomic_coord: int
        :return: the cDNA coordinate corresponding to the given genomic coordinate.
        """

        cdna_offset = 0
        for cdna_match in self.ordered_cdna_match:
            # Inside the exon.
            if (
                cdna_match.tx_position.chrom_start
                <= genomic_coord
                <= cdna_match.tx_position.chrom_stop
            ):
                if self.strand == "+":
                    position = genomic_coord - (cdna_match.tx_position.chrom_start + 1)
                else:
                    position = cdna_match.tx_position.chrom_stop - genomic_coord
                offset = cdna_match.get_offset(position, validate=False)
                return cdna_offset + position + offset + 1
            cdna_offset += cdna_match.length

        raise ValueError(f"Couldn't find {genomic_coord=}!")

    def find_exon_number(self, offset: int) -> int:
        """
        The function `find_exon_number` returns the exon number for a given position.

        :param offset: The offset parameter represents a position in the genome. It is an integer value
        that indicates the position of interest within the genome
        :type offset: int
        :return: an integer value, which represents the exon number for a given position.
        """

        exon_number = 1
        for cdna_match in self.ordered_cdna_match:
            if (
                cdna_match.tx_position.chrom_start
                <= offset
                <= cdna_match.tx_position.chrom_stop
            ):
                return exon_number
            exon_number += 1
        return None


BED6Interval_base = namedtuple(
    "BED6Interval_base",
    ("chrom", "chrom_start", "chrom_end", "name", "score", "strand"),
)


# The BED6Interval class is a subclass of BED6Interval_base.
class BED6Interval(BED6Interval_base):

    def distance(self, offset: int) -> int:
        """
        The `distance` function calculates the distance between an offset and an interval, returning
        zero if the offset is inside the interval, a positive value if the interval comes after the
        offset, and a negative value if the interval comes before the offset.
        if offset is inside the exon, distance is zero.
        otherwise, distance is the distance to the nearest edge.
        distance is positive if the exon comes after the offset.
        distance is negative if the exon comes before the offset.

        :param offset: The offset parameter represents a position or point in the genome. It is an
        integer value that indicates the position within the genome sequence
        :type offset: int
        :return: an integer value, which represents the distance to the interval.
        """

        start = self.chrom_start + 1
        end = self.chrom_end

        if start <= offset <= end:
            return 0

        start_distance = start - offset
        end_distance = offset - end

        if abs(start_distance) < abs(end_distance):
            return start_distance
        else:
            return -end_distance


# The Exon class is a blueprint for creating objects that represent exons in genetic sequences.
# We still need exons to work out the flanking boundaries
class Exon(object):

    def __init__(self, transcript: Transcript, tx_position: int, number: int) -> None:
        """
        The function initializes an object with a transcript, a position in the transcript, and a
        number.

        :param transcript: The `transcript` parameter is of type `Transcript`. It represents a
        transcript object that contains information about a conversation or dialogue
        :type transcript: Transcript
        :param tx_position: The `tx_position` parameter represents the position of the transcript in a
        list or array. It is an integer value that indicates the index of the transcript in the list or
        array
        :type tx_position: int
        :param number: The "number" parameter is an integer that represents a specific number. It is
        used as a parameter in the constructor of a class
        :type number: int
        """

        self.transcript = transcript
        self.tx_position = tx_position
        self.number = number

    @property
    def name(self) -> str:
        """
        The function returns a string that combines the name of the transcript and a number.
        :return: a string that combines the name of the transcript with the number. The format of the
        string is "{transcript name}.{number}".
        """

        return "%s.%d" % (self.transcript.name, self.number)

    def get_as_interval(self, coding_only: bool = False) -> object:
        """
        The function `get_as_interval` returns the coding region for an exon as a BED6Interval object.
        This function returns a BED6Interval objects containing  position
        information for this exon. This may be used as input for
        pybedtools.create_interval_from_list() after casting chrom_start
        and chrom_end as strings.

        :param coding_only: The `coding_only` parameter is a boolean flag that determines whether to
        include only exons in the coding region. If `coding_only` is set to `True`, the function will
        check if the exon is completely outside the coding region defined by the transcript's CDS
        (coding sequence) position, defaults to False
        :type coding_only: bool (optional)
        :return: a BED6Interval object.
        """

        exon_start = self.tx_position.chrom_start
        exon_stop = self.tx_position.chrom_stop

        # Get only exon coding region if requested
        if coding_only:
            if (
                exon_stop <= self.transcript.cds_position.chrom_start
                or exon_start >= self.transcript.cds_position.chrom_stop
            ):
                return None
            exon_start = max(exon_start, self.transcript.cds_position.chrom_start)
            exon_stop = min(
                max(exon_stop, self.transcript.cds_position.chrom_start),
                self.transcript.cds_position.chrom_stop,
            )

        return BED6Interval(
            self.tx_position.chrom,
            exon_start,
            exon_stop,
            self.name,
            ".",
            self.strand,
        )

    @property
    def strand(self):
        strand = "+" if self.tx_position.is_forward_strand else "-"
        return strand


# The CDNA_Match class is a subclass of the Exon class.
class CDNA_Match(Exon):

    def __init__(
        self,
        transcript: Transcript,
        tx_position: int,
        cdna_start: int,
        cdna_end: int,
        gap: int,
        number: int,
    ) -> None:
        """
        The function initializes a CDNA_Match object with specified attributes.

        :param transcript: The `transcript` parameter is an instance of the `Transcript` class. It
        represents the transcript that the CDNA match belongs to
        :type transcript: Transcript
        :param tx_position: The `tx_position` parameter represents the position of the transcript in the
        genome. It is an integer value
        :type tx_position: int
        :param cdna_start: The `cdna_start` parameter represents the starting position of the cDNA
        match. It is an integer value
        :type cdna_start: int
        :param cdna_end: The `cdna_end` parameter represents the end position of the cDNA match. It is
        an integer value that indicates the position of the last nucleotide in the cDNA sequence that
        matches the transcript
        :type cdna_end: int
        :param gap: The "gap" parameter represents the number of nucleotides that are missing or
        inserted in the cDNA sequence compared to the reference transcript sequence. It indicates the
        presence of gaps or insertions in the alignment between the cDNA and the reference transcript
        :type gap: int
        :param number: The `number` parameter represents the number of the CDNA match. It is used to
        uniquely identify each CDNA match object
        :type number: int
        """

        super(CDNA_Match, self).__init__(transcript, tx_position, number)
        self.cdna_start = cdna_start
        self.cdna_end = cdna_end
        self.gap = gap

    @property
    def length(self) -> int:
        """
        The function calculates the length of a sequence by subtracting the start position from the end
        position and adding 1.
        :return: The length of the sequence, calculated by subtracting the cdna_start from the cdna_end
        and adding 1.
        """

        return self.cdna_end - self.cdna_start + 1

    def get_offset(self, position: int, validate: bool = True) -> int:
        """
        The `get_offset` function calculates the offset for a given position in a cDNA sequence based on
        the GAP attribute.
        cdna_match GAP attribute looks like: 'M185 I3 M250' which is code/length
        @see https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md#the-gap-attribute
        codes operation
        M 	match
        I 	insert a gap into the reference sequence
        D 	insert a gap into the target (delete from reference)
        If you want the whole exon, then pass the end

        :param position: The `position` parameter is an integer that represents the position in the
        sequence. It is used to calculate the offset based on the GAP attribute of the cDNA match
        :type position: int
        :param validate: The `validate` parameter is a boolean flag that determines whether to perform
        validation checks during the calculation of the offset. If `validate` is set to `True`, the
        function will raise a `ValueError` if the given position falls within an insertion or deletion
        gap. If `validate` is set, defaults to True
        :type validate: bool (optional)
        :return: an integer value representing the offset for a given position in the cDNA sequence.
        """

        if not self.gap:
            return 0

        position_1_based = position + 1
        cdna_match_index = 1
        offset = 0
        for gap_op in self.gap.split():
            code = gap_op[0]
            length = int(gap_op[1:])
            if code == "M":
                cdna_match_index += length
            elif code == "I":
                if validate and position_1_based < cdna_match_index + length:
                    raise ValueError(
                        "Coordinate (%d) inside insertion (%s) - no mapping possible!"
                        % (position_1_based, gap_op)
                    )
                offset += length
            elif code == "D":
                if validate and position < cdna_match_index + length:
                    raise ValueError(
                        "Coordinate (%d) inside deletion (%s) - no mapping possible!"
                        % (position_1_based, gap_op)
                    )
                offset -= length
            else:
                raise ValueError(f"Unknown code in cDNA GAP: {gap_op}")

            if cdna_match_index > position_1_based:
                break

        return offset
