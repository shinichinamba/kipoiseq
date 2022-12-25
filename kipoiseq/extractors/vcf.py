import abc
import logging
from typing import Tuple, Iterable, List, Union, Iterator
from itertools import islice
from collections import defaultdict
from kipoiseq.dataclasses import Variant, Interval
from kipoiseq.variant_source import VariantFetcher
from kipoiseq.extractors.vcf_query import VariantIntervalQueryable
from kipoiseq.utils import batch_iter

try:
    from cyvcf2 import VCF
except ImportError:
    VCF = object

__all__ = [
    'MultiSampleVCF'
]


class MultiSampleVCF(VariantFetcher, VCF):

    def __init__(self, *args, **kwargs):
        from cyvcf2 import VCF
        VCF.__init__(self, *args, **kwargs, strict_gt=True)
        self.sample_mapping = dict(zip(self.samples, range(len(self.samples))))

    def fetch_variants(self, interval, sample_id=None, phase=None):
        region = self._region(interval)
        for cy_variant in self(region):
            for variant in self._variants_from_cyvcf2(cy_variant):
                if sample_id is None:
                    has_variant = True
                else:
                    if phase is None:
                        has_variant = self.has_variant(variant, sample_id)
                    else:
                        has_variant = self.has_variant_phased(variant, sample_id, phase)
                if has_variant:
                    yield variant

    def _variants_from_cyvcf2(self, cy_variant):
        # in case deletion is present
        ALTs = cy_variant.ALT or ['']
        # single REF can have multiple ALT
        for alt in ALTs:
            v = Variant.from_cyvcf_and_given_alt(cy_variant, alt)
            if 'N' in alt or '*' in alt:
                logging.warning(
                    'Undefined variant %s are not supported: Skip' % str(v))
                continue
            yield v

    @staticmethod
    def _region(interval: Interval):
        return '%s:%d-%d' % (interval.chrom, interval.start + 1, interval.end)

    def has_variant(self, variant: Variant, sample_id: str) -> bool:
        gt_type = variant.source.gt_types[self.sample_mapping[sample_id]]
        return self._has_variant_gt(gt_type)

    def has_variant_phased(self, variant: Variant, sample_id: str, phase: int) -> bool:
        genotypes = variant.source.genotypes[self.sample_mapping[sample_id]]
        return self._has_variant_phased(genotypes, phase)

    @staticmethod
    def _has_variant_gt(gt_type: int) -> bool:
        return gt_type != 0 and gt_type != 2
    
    @staticmethod
    def _has_variant_phased(genotypes: list, phase: int) -> bool:
        """
        # Arguments
            phase: phasing index, either int(0) or int(1).
        """
        if phase == 0 or phase == 1:
            if genotypes[2]:
                return genotypes[phase] == 1
            else:
                logging.warning(
                    'Unphased genotype was ignored')
                return False
        else:
            raise ValueError("phase must be int(0) or int(1)")

    def __next__(self):
        return Variant.from_cyvcf(super().__next__())

    def __iter__(self) -> Iterable[Variant]:
        while True:
            try:
                cy_variant = super().__next__()
            except StopIteration:
                break
            yield from self._variants_from_cyvcf2(cy_variant)

    def batch_iter(self, batch_size=10000) -> Iterable[Iterable[Variant]]:
        """Iterates variants in vcf file.

        # Arguments
            vcf_file: path of vcf file.
            batch_size: size of each batch.
        """
        variants = iter(self)
        yield from batch_iter(variants, batch_size=batch_size)

    def query_variants(self, intervals: List[Interval], sample_id=None, phase=None,
                       progress=False) -> VariantIntervalQueryable:
        """Fetch variants for given multi-intervals from vcf file
        for sample if sample id is given.
        To consider phasing, specify phase as either int(0) or int(1).
        If phase is given, unphased variants will be ignored.

        # Arguments
            intervals (List[pybedtools.Interval]): list of Interval objects
            sample_id (str, optional): sample id in vcf file.
            phase (int, optional): phase of the variant to be considered

        # Returns
          VCFQueryable: queryable object whihc allow you to query the
            fetched variatns.

        # Example
            To fetch variants if quality more than 10 and there is
            a variant in interval
            ```
              >>> MultiSampleVCF(vcf_path) \
                    .query_variants(intervals) \
                    .filter(lambda variant: variant.qual > 10) \
                    .filter_range(NumberVariantQuery(max_num=1)) \
                    .to_vcf(output_path)
            ```
        """
        pairs = [(self.fetch_variants(i, sample_id=sample_id, phase=phase), i)
                 for i in intervals]
        return VariantIntervalQueryable(self, pairs, progress=progress)

    def query_all(self, progress=False) -> VariantIntervalQueryable:
        """Convert all variants into queryable object without interval so
        interval filters will not work.

        # Example
            To fetch variants if quality more than 10
            ```
              >>> MultiSampleVCF(vcf_path) \
                    .query_all() \
                    .filter(lambda variant: variant.qual > 10) \
                    .to_vcf(output_path)
            ```
        """
        pairs = [(iter(self), None)]
        return VariantIntervalQueryable(self, pairs, progress=progress)

    def get_variant(self, variant: Union[Variant, str]) -> Variant:
        """Returns variant from vcf file. Lets you use vcf file as dict.

        # Arguments:
            variant: variant object or variant id as string.

        # Returns
            Variant object.

        # Example
            ```python
              >>> MultiSampleVCF(vcf_path).get_variant("chr1:4:T:['C']")
            ```
        """
        if type(variant) == str:
            variant = Variant.from_str(variant)

        variants = self.fetch_variants(
            Interval(variant.chrom, variant.pos - 1, variant.pos))
        for v in variants:
            if v.ref == variant.ref and v.alt == variant.alt:
                return v
        raise KeyError('Variant %s not found in vcf file.' % str(variant))

    def get_variants(self, variants: Iterable[Union[str, Variant]],
                     regions=None, variant_gap=150) -> List[Variant]:
        """Returns list of variants from vcf file. Lets you use vcf file as dict.

        # Arguments:
            variants: list of variants
            regions: list of regions to seek for variants.
              Automatically generated from variants if not given.
            strategy: strategy if there is not variant in region.

        # Returns
           List of variants
        """
        variants = [
            Variant.from_str(v) if type(v) == str else v
            for v in variants
        ]
        regions = regions or self._regions_from_variants(
            variants, variant_gap=variant_gap)
        variant_map = dict()

        for r in regions:
            r_variants = self.fetch_variants(r)
            for v in r_variants:
                variant_map[v] = v

        return [variant_map.get(v) for v in variants]

    def _regions_from_variants(self, variants: List[Variant], variant_gap=150):
        regions = list()

        for chrom, vs in self._group_variants_by_chrom(variants).items():
            starts = sorted(v.pos for v in vs)

            start_i = starts[0]
            prev_i = starts[0]
            for i in starts[1:]:
                if prev_i + 150 < i:
                    regions.append(Interval(chrom, start_i - 1, prev_i))
                    start_i = i
                prev_i = i

            regions.append(Interval(chrom, start_i - 1, prev_i))

        return regions

    def _group_variants_by_chrom(self, variants: List[Variant]):
        chroms = defaultdict(set)

        for v in variants:
            chroms[v.chrom].add(v)

        return dict(chroms)

    def get_samples(self, variant: Variant) -> str:
        """Fetchs sample names which have given variants

        # Arguments
            variant: variant object.

        # Returns
            Dict[str, int]: Dict of sample which have variant and gt as value.
        """
        return dict(filter(lambda x: self._has_variant_gt(x[1]),
                           zip(self.samples, variant.source.gt_types)))
