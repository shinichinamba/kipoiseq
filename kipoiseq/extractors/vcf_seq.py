import abc
import logging
from typing import Union

from pyfaidx import Sequence, complement
from kipoiseq.dataclasses import Interval
from kipoiseq.extractors import (
    BaseExtractor,
    FastaStringExtractor,
    MultiSampleVCF,
)

from kipoiseq import __version__
from deprecation import deprecated

__all__ = [
    'VariantSeqExtractor',
    'SingleVariantVCFSeqExtractor',
    'SingleSeqVCFSeqExtractor'
]


class IntervalSeqBuilder(list):
    """
    String builder for `pyfaidx.Sequence` and `Interval` objects.
    """

    def restore(self, sequence: Sequence, fixed_len=False):
        """
        Args:
          sequence: `pyfaidx.Sequence` which convert all interval inside
            to `Seqeunce` objects.
          fixed_len: if True, the return sequence will have the same length
            as the `interval` by N-padding
        """
        for i, interval in enumerate(self):
            # interval.end can be bigger than interval.start
            interval_len = max(0, interval.end - interval.start)

            if type(self[i]) == Interval:
                start = interval.start - sequence.start
                end = start + interval_len
                if fixed_len:
                    start_overflow = min(-start, interval_len)
                    sequence_len = sequence.end - sequence.start
                    end_overflow = min(end - sequence_len, interval_len)
                    if (start_overflow > 0):
                        pad_start = start_overflow * 'N'
                        start = max(0, start)
                    else:
                        pad_start = ''
                        start = min(start, sequence_len)
                    if (end_overflow > 0):
                        pad_end = end_overflow * 'N'
                        end = min(end, sequence_len)
                    else:
                        pad_end = ''
                        end = max(0, end)
                    self[i] = Sequence(
                        seq = pad_start + sequence[start: end].seq + pad_end,
                        start = interval.start,
                        end = interval.start + interval_len)
                else:
                    self[i] = sequence[start: end]

    def _concat(self):
        for sequence in self:
            if type(sequence) != Sequence:
                raise TypeError('Intervals should be restored with `restore`'
                                ' method before calling concat method!')
            yield sequence.seq

    def concat(self):
        """
        Build the string from sequence objects.

        Returns:
          str: the final sequence.
        """
        return ''.join(self._concat())


class VariantSeqExtractor(BaseExtractor):

    def __init__(self, fasta_file: str = None, reference_sequence: BaseExtractor = None, use_strand=True):
        """
        Sequence extractor which allows to obtain the alternative sequence,
        given some interval and variants inside this interval.

        Args:
            fasta_file: path to the fasta file (can be gzipped)
            reference_sequence: extractor returning the reference sequence given some interval
            use_strand (bool): if True, the extracted sequence
                is reverse complemented in case interval.strand == "-"
        """
        self._use_strand = use_strand

        if fasta_file is not None:
            if reference_sequence is not None:
                raise ValueError(
                    "either fasta_file or ref_seq_extractor have to be specified")
            self._ref_seq_extractor = FastaStringExtractor(
                fasta_file, use_strand=False)
        else:
            if reference_sequence is None:
                raise ValueError(
                    "either fasta_file or ref_seq_extractor have to be specified")
            self._ref_seq_extractor = reference_sequence

    @property
    @deprecated(deprecated_in="1.0",
                # removed_in="2.0",
                current_version=__version__,
                details="Use `ref_seq_extractor` instead")
    def fasta(self):
        return self._ref_seq_extractor

    @property
    def ref_seq_extractor(self) -> BaseExtractor:
        """

        Returns:
            The reference sequence extractor of this object
        """
        return self._ref_seq_extractor

    def extract(self, interval, variants, anchor, fixed_len=True, use_strand=None, **kwargs):
        """

        Args:
            interval: pybedtools.Interval Region of interest from
                which to query the sequence. 0-based
            variants: List[cyvcf2.Variant]: variants overlapping the `interval`.
                can also be indels. 1-based
            anchor: absolution position w.r.t. the interval start. (0-based).
                E.g. for an interval of `chr1:10-20` the anchor of 10 denotes
                the point chr1:10 in the 0-based coordinate system.
            fixed_len: if True, the return sequence will have the same length
                as the `interval` (e.g. `interval.end - interval.start`)
            use_strand (bool, optional): if True, the extracted sequence
                is reverse complemented in case interval.strand == "-".
                Overrides `self.use_strand`

        Returns:
            A single sequence (`str`) with all the variants applied.
        """
        first_half_results = self._extract_first_half(interval, variants, anchor, fixed_len)
        return self._extract_second_half(interval, first_half_results, fixed_len, use_strand)
    
    def _extract_first_half(self, interval, variants, anchor, fixed_len):
        """
        The first half of the extract() function.
        The extract() function was splitted for the ad-hoc speed up of SingleSeqVCFSeqExtractor.extract()
        """
        # Preprocessing
        anchor = max(min(anchor, interval.end), interval.start)
        variant_pairs = self._variant_to_sequence(variants)

        # 1. Split variants overlapping with anchor
        # and interval start end if not fixed_len
        variant_pairs = self._split_overlapping(variant_pairs, anchor)

        if not fixed_len:
            variant_pairs = self._split_overlapping(
                variant_pairs, interval.start, which='right')
            variant_pairs = self._split_overlapping(
                variant_pairs, interval.end, which='left')

        variant_pairs = list(variant_pairs)

        # 2. split the variants into upstream and downstream
        # and sort the variants in each interval
        upstream_variants = sorted(
            filter(lambda x: x[0].start >= anchor, variant_pairs),
            key=lambda x: x[0].start
        )

        downstream_variants = sorted(
            filter(lambda x: x[0].start < anchor, variant_pairs),
            key=lambda x: x[0].start,
            reverse=True
        )

        # 3. Extend start and end position for deletions
        if fixed_len:
            istart, iend = self._updated_interval(
                interval, upstream_variants, downstream_variants)
        else:
            istart, iend = interval.start, interval.end
        
        return (anchor, istart, iend, upstream_variants, downstream_variants)


    def _extract_second_half(self, interval, first_half_results, fixed_len, use_strand=None):
        """
        The second half of the extract() function
        """
        anchor, istart, iend, upstream_variants, downstream_variants = first_half_results
        # 4. Iterate from the anchor point outwards. At each
        # register the interval from which to take the reference sequence
        # as well as the interval for the variant
        down_sb = self._downstream_builder(
            downstream_variants, interval, anchor, istart)

        up_sb = self._upstream_builder(
            upstream_variants, interval, anchor, iend)

        # 5. fetch the sequence and restore intervals in builder
        seq = self._fetch(interval, istart, iend, error_if_invalid = not fixed_len)
        up_sb.restore(seq, fixed_len=fixed_len)
        down_sb.restore(seq, fixed_len=fixed_len)

        # 6. Concate sequences from the upstream and downstream splits. Concat
        # upstream and downstream sequence. Cut to fix the length.
        down_str = down_sb.concat()
        up_str = up_sb.concat()

        if fixed_len:
            down_str, up_str = self._cut_to_fix_len(
                down_str, up_str, interval, anchor)

        seq = down_str + up_str

        if use_strand is None:
            use_strand = self.use_strand
        if use_strand and interval.strand == '-':
            # reverse-complement
            seq = complement(seq)[::-1]

        return seq

    @staticmethod
    def _variant_to_sequence(variants):
        """
        Convert `cyvcf2.Variant` objects to `pyfaidx.Seqeunce` objects
        for reference and variants.
        """
        for v in variants:
            ref = Sequence(name=v.chrom, seq=v.ref,
                           start=v.start, end=v.start + len(v.ref))
            alt = Sequence(name=v.chrom, seq=v.alt,
                           start=v.start, end=v.start + len(v.alt))
            yield ref, alt

    @staticmethod
    def _split_overlapping(variant_pairs, anchor, which='both'):
        """
        Split the variants hitting the anchor into two
        """
        for ref, alt in variant_pairs:
            if ref.start < anchor < ref.end:
                mid = anchor - ref.start
                if which == 'left' or which == 'both':
                    yield ref[:mid], alt[:mid]
                if which == 'right' or which == 'both':
                    yield ref[mid:], alt[mid:]
            else:
                yield ref, alt

    @staticmethod
    def _updated_interval(interval, up_variants, down_variants):
        istart = interval.start
        iend = interval.end

        for ref, alt in up_variants:
            if iend <= ref.start:
                break
            diff_len = len(alt) - len(ref)
            if diff_len < 0:
                iend -= diff_len

        for ref, alt in down_variants:
            if istart >= ref.end:
                break
            diff_len = len(alt) - len(ref)
            if diff_len < 0:
                istart += diff_len

        return istart, iend

    @staticmethod
    def _downstream_builder(down_variants, interval, anchor, istart):
        down_sb = IntervalSeqBuilder()

        prev = anchor
        for ref, alt in down_variants:
            if ref.end <= istart:
                break
            down_sb.append(Interval(interval.chrom, ref.end, prev))
            down_sb.append(alt)
            prev = ref.start
        down_sb.append(Interval(interval.chrom, istart, prev))
        down_sb.reverse()

        return down_sb

    @staticmethod
    def _upstream_builder(up_variants, interval, anchor, iend):
        up_sb = IntervalSeqBuilder()

        prev = anchor
        for ref, alt in up_variants:
            if ref.start >= iend:
                break
            up_sb.append(Interval(interval.chrom, prev, ref.start))
            up_sb.append(alt)
            prev = ref.end
        up_sb.append(Interval(interval.chrom, prev, iend))

        return up_sb

    def _fetch(self, interval, istart, iend, error_if_invalid = True):
        # fetch interval, ignore strand
        if istart < 0:
            if error_if_invalid:
                raise ValueError(
                    'Requested start coordinate was negative and set to 0: %s' % str(istart))
            else:
                logging.warning(
                    'Requested start coordinate was negative and set to 0: %s' % str(istart))

            istart = 0
        seq = self.ref_seq_extractor.extract(
            Interval(interval.chrom, istart, iend))
        iend_fetched = istart+len(seq)
        if iend_fetched != iend:
            if error_if_invalid:
                raise ValueError(
                    'Requested end coordinate was larger than that of Sequence and the former was set to the latter: %s to %s'
                    % (str(iend), str(iend_fetched)))
            else:
                logging.warning(
                    'Requested end coordinate was larger than that of Sequence and the former was set to the latter: %s to %s'
                    % (str(iend), str(iend_fetched)))
        seq = Sequence(name=interval.chrom, seq=seq, start=istart, end=iend_fetched)
        return seq

    @staticmethod
    def _cut_to_fix_len(down_str, up_str, interval, anchor):
        down_len = anchor - interval.start
        up_len = interval.end - anchor
        down_str = down_str[-down_len:] if down_len else ''
        up_str = up_str[: up_len] if up_len else ''
        return down_str, up_str


class _BaseVCFSeqExtractor(BaseExtractor, metaclass=abc.ABCMeta):
    """
    Base class to fetch sequence in which variants applied based
    on given vcf file.
    """

    def __init__(self, fasta_file, vcf_file):
        """
        Args:
          fasta_file: path to the fasta file (can be gzipped)
          vcf_file: path to the fasta file (need be bgzipped and indexed)
        """
        self.fasta_file = fasta_file
        self.vcf_file = vcf_file
        self.variant_extractor = VariantSeqExtractor(fasta_file)
        self.vcf = MultiSampleVCF(vcf_file)

    @abc.abstractmethod
    def extract(self, interval: Interval, *args, **kwargs) -> str:
        raise NotImplementedError()


class SingleVariantVCFSeqExtractor(_BaseVCFSeqExtractor):
    """
    Fetch list of sequence in which each variant applied based
    on given vcf file.
    """

    def extract(self, interval, anchor=None, sample_id=None, phase=None, fixed_len=True):
        for variant in self.vcf.fetch_variants(interval, sample_id, phase):
            yield self.variant_extractor.extract(
                interval,
                variants=[variant],
                anchor=anchor,
                fixed_len=fixed_len
            )


class SingleSeqVCFSeqExtractor(_BaseVCFSeqExtractor):
    """
    Fetch sequence in which all variant applied based on given vcf file.
    Use _extract_first_half() and _extract_second_half() directly for the ad-hoc speed up
    """

    def extract(self, interval, anchor=None, sample_id=None, phase=None, fixed_len=True):
        variants=self.vcf.fetch_variants(interval, sample_id, phase)
        if not fixed_len:
            seq = self.variant_extractor.extract(
                interval,
                variants=variants,
                anchor=anchor,
                fixed_len=fixed_len)
        else:
            # first, we try the sequence extraction with variants in the interval
            anchor, istart, iend, upstream_variants, downstream_variants = self.variant_extractor._extract_first_half(
                interval, variants, anchor, fixed_len)
            
            # We iteratively extend the interval if the interval is insufficient to get the fixed-length sequence due to deletions
            downstream_additional_width = interval.start - istart
            if downstream_additional_width > 0:
                istart = interval.start
            while downstream_additional_width > 0:
                downstream_variants, istart, downstream_additional_width = self._add_downstream_variants(
                    sample_id, phase, downstream_variants, interval.chrom, istart, anchor, downstream_additional_width)
            
            upstream_additional_width = iend - interval.end
            if upstream_additional_width > 0:
                iend = interval.end
            while upstream_additional_width > 0:
                upstream_variants, iend, upstream_additional_width = self._add_upstream_variants(
                    sample_id, phase, upstream_variants, interval.chrom, iend, upstream_additional_width)
            
            # Getting the sequence
            seq = self.variant_extractor._extract_second_half(
                interval, (anchor, istart, iend, upstream_variants, downstream_variants), fixed_len)
        return seq

    def _add_downstream_variants(self, sample_id, phase, downstream_variants, chrom, istart, anchor, additional_width):
        # Get additional downstream variants
        # We don't have to calculate the istart/iend precisely, 
        # rather we'd like to update istart so that (anchor - istart) got equal to or larger than the required width
        new_istart = istart - additional_width
        add_down_variants = self.vcf.fetch_variants(Interval(chrom, max(0, new_istart), max(0, istart)), sample_id, phase)
        add_down_variants = list(self.variant_extractor._variant_to_sequence(add_down_variants))
        add_down_variants.sort(key=lambda x: x[0].start, reverse=True)
        downstream_variants.extend(add_down_variants)
        new_add_width = 0
        for ref, alt in add_down_variants:
            new_add_width -= len(alt) - len(ref)
            if ref.end > anchor:
                # This function is called only when there are some downstream deletions in the previous interval.
                # Therefore, if the end position is larger than the anchor position, the variants must be overlapped with
                # at least one deletions.
                raise ValueError('There are at least two overlapped variants (Note that the consistency of variants is not routinely checked in this function)')
        return (downstream_variants, new_istart, new_add_width)

    def _add_upstream_variants(self, sample_id, phase, upstream_variants, chrom, iend, additional_width):
        # Get additional upstream variants
        new_iend = iend + additional_width
        add_up_variants = self.vcf.fetch_variants(Interval(chrom, iend, new_iend), sample_id, phase)
        add_up_variants = list(self.variant_extractor._variant_to_sequence(add_up_variants))
        add_up_variants.sort(key=lambda x: x[0].start)
        upstream_variants.extend(add_up_variants)
        new_add_width = 0
        for ref, alt in add_up_variants:
            new_add_width -= len(alt) - len(ref)
        return (upstream_variants, new_iend, new_add_width)
