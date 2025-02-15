import pytest
from conftest import vcf_file, phased_vcf_file, sample_5kb_fasta_file
from cyvcf2 import VCF
from pyfaidx import Sequence
from kipoiseq.dataclasses import Variant, Interval
from kipoiseq.extractors.vcf_seq import IntervalSeqBuilder, \
    VariantSeqExtractor, SingleSeqVCFSeqExtractor, SingleVariantVCFSeqExtractor
from kipoiseq.extractors import FastaStringExtractor

fasta_file = sample_5kb_fasta_file

intervals = [
    Interval('chr1', 4, 10),
    Interval('chr1', 5, 30),
    Interval('chr1', 20, 30)
]


@pytest.fixture
def interval_seq_builder():
    return IntervalSeqBuilder([
        Interval('chr1', 10, 13),
        Interval('chr1', 13, 14),
        Sequence(seq='TAGC', start=14, end=18),
        Interval('chr1', 18, 20)
    ])


def test_interval_seq_builder_restore(interval_seq_builder):
    sequence = Sequence(seq='CCCCATCGTT', start=10, end=20)
    interval_seq_builder.restore(sequence)
    assert interval_seq_builder[0].seq == 'CCC'
    assert interval_seq_builder[1].seq == 'C'
    assert interval_seq_builder[2].seq == 'TAGC'
    assert interval_seq_builder[3].seq == 'TT'

    interval_seq_builder.append(Interval('chr1', 5, 10))
    interval_seq_builder.restore(sequence)
    assert interval_seq_builder[4].seq == ''

    interval_seq_builder.append(Interval('chr1', 20, 25))
    interval_seq_builder.restore(sequence)
    assert interval_seq_builder[5].seq == ''

    interval_seq_builder.append(Interval('chr1', 10, 5))
    interval_seq_builder.restore(sequence)
    assert interval_seq_builder[6].seq == ''

    interval_seq_builder.append(Interval('chr1', 25, 20))
    interval_seq_builder.restore(sequence)
    assert interval_seq_builder[7].seq == ''

    interval_seq_builder.append(Interval('chr1', 5, 7))
    interval_seq_builder.restore(sequence, fixed_len=True)
    assert interval_seq_builder[8].seq == 'NN'

    interval_seq_builder.append(Interval('chr1', 9, 11))
    interval_seq_builder.restore(sequence, fixed_len=True)
    assert interval_seq_builder[9].seq == 'NC'

    interval_seq_builder.append(Interval('chr1', 19, 21))
    interval_seq_builder.restore(sequence, fixed_len=True) 
    assert interval_seq_builder[10].seq == 'TN'

    interval_seq_builder.append(Interval('chr1', 23, 25))
    interval_seq_builder.restore(sequence, fixed_len=True) 
    assert interval_seq_builder[11].seq == 'NN'

    interval_seq_builder.append(Interval('chr1', 5, 25))
    interval_seq_builder.restore(sequence, fixed_len=True) 
    assert interval_seq_builder[12].seq == 'NNNNNCCCCATCGTTNNNNN'

    interval_seq_builder.append(Interval('chr1', 11, 9))
    interval_seq_builder.restore(sequence, fixed_len=True)
    assert interval_seq_builder[13].seq == ''

    interval_seq_builder.append(Interval('chr1', 21, 19))
    interval_seq_builder.restore(sequence, fixed_len=True)
    assert interval_seq_builder[14].seq == '' 
    

def test_interval_seq_builder_concat(interval_seq_builder):
    with pytest.raises(TypeError):
        interval_seq_builder.concat()

    sequence = Sequence(seq='CCCCATCGNN', start=10, end=20)
    interval_seq_builder.restore(sequence)
    assert interval_seq_builder.concat() == 'CCCCTAGCNN'


@pytest.fixture
def variant_seq_extractor():
    return VariantSeqExtractor(fasta_file)


def test__split_overlapping(variant_seq_extractor):
    pair = (Sequence(seq='AAA', start=3, end=6),
            Sequence(seq='T', start=3, end=4))
    splited_pairs = list(variant_seq_extractor._split_overlapping([pair], 5))

    assert splited_pairs[0][0].seq == 'AA'
    assert splited_pairs[0][1].seq == 'T'
    assert splited_pairs[1][0].seq == 'A'
    assert splited_pairs[1][1].seq == ''

    pair = (Sequence(seq='TT', start=3, end=5),
            Sequence(seq='AAA', start=3, end=6))
    splited_pairs = list(variant_seq_extractor._split_overlapping([pair], 4))

    assert splited_pairs[0][0].seq == 'T'
    assert splited_pairs[0][1].seq == 'A'
    assert splited_pairs[1][0].seq == 'T'
    assert splited_pairs[1][1].seq == 'AA'


def test_extract(variant_seq_extractor):
    variants = [Variant.from_cyvcf(v) for v in VCF(vcf_file)]

    interval = Interval('chr1', 2, 9)

    seq = variant_seq_extractor.extract(interval, variants, anchor=5)
    assert len(seq) == interval.end - interval.start
    assert seq == 'CGAACGT'

    interval = Interval('chr1', 2, 9, strand='-')
    seq = variant_seq_extractor.extract(interval, variants, anchor=5)
    assert len(seq) == interval.end - interval.start
    assert seq == 'ACGTTCG'

    interval = Interval('chr1', 4, 14)
    seq = variant_seq_extractor.extract(interval, variants, anchor=7)
    assert len(seq) == interval.end - interval.start
    assert seq == 'AACGTAACGT'

    interval = Interval('chr1', 4, 14)
    seq = variant_seq_extractor.extract(interval, variants, anchor=4)
    assert len(seq) == interval.end - interval.start
    assert seq == 'GAACGTAACG'

    interval = Interval('chr1', 2, 5)
    seq = variant_seq_extractor.extract(interval, variants, anchor=3)
    assert len(seq) == interval.end - interval.start
    assert seq == 'GCG'

    interval = Interval('chr1', 24, 34)
    seq = variant_seq_extractor.extract(interval, variants, anchor=27)
    assert len(seq) == interval.end - interval.start
    assert seq == 'TGATAACGTA'

    interval = Interval('chr1', 25, 35)
    seq = variant_seq_extractor.extract(interval, variants, anchor=34)
    assert len(seq) == interval.end - interval.start
    assert seq == 'TGATAACGTA'

    interval = Interval('chr1', 34, 44)
    seq = variant_seq_extractor.extract(interval, variants, anchor=37)
    assert len(seq) == interval.end - interval.start
    assert seq == 'AACGTAACGT'

    interval = Interval('chr1', 34, 44)
    seq = variant_seq_extractor.extract(interval, variants, anchor=100)
    assert len(seq) == interval.end - interval.start
    assert seq == 'AACGTAACGT'

    interval = Interval('chr1', 5, 11, strand='+')
    seq = variant_seq_extractor.extract(
        interval, variants, anchor=10, fixed_len=False)
    assert seq == 'ACGTAA'

    interval = Interval('chr1', 0, 3, strand='+')
    seq = variant_seq_extractor.extract(
        interval, variants, anchor=10, fixed_len=False)
    assert seq == 'ACG'

    interval = Interval('chr1', 0, 3, strand='+')
    ref_seq_extractor = FastaStringExtractor(fasta_file, use_strand=True)
    seq = VariantSeqExtractor(reference_sequence=ref_seq_extractor).extract(
        interval, variants, anchor=10, fixed_len=False)
    assert seq == 'ACG'

    interval = Interval('chr1', 4999, 5001, strand='+')
    with pytest.raises(ValueError):
        variant_seq_extractor.extract(
            interval, variants, anchor=10, fixed_len=False)
    seq = variant_seq_extractor.extract(
        interval, variants, anchor=10, fixed_len=True)
    assert seq == 'AN'

    interval = Interval('chr1', -1, 2, strand='+')
    with pytest.raises(ValueError):
        variant_seq_extractor.extract(
            interval, variants, anchor=10, fixed_len=False)
    seq = variant_seq_extractor.extract(
        interval, variants, anchor=10, fixed_len=True)
    assert seq == 'NAC'

    interval2 = Interval('chr1', 27, 30, strand='+')
    seq = variant_seq_extractor.extract(
        interval2, variants, anchor=30, fixed_len=False)
    assert seq == 'TA'
    seq = variant_seq_extractor.extract(
        interval2, variants, anchor=30, fixed_len=True)
    assert seq == 'ATA'
    
    interval3 = Interval('chr1', -1, 2, strand='-')
    seq_list = variant_seq_extractor.batch_extract(
        [interval, interval2, interval3], variants, anchors=[10, 30, 10], fixed_len=True)
    assert seq_list == ['NAC', 'ATA', 'GTN']

@pytest.fixture
def single_variant_vcf_seq_extractor():
    return SingleVariantVCFSeqExtractor(fasta_file, vcf_file)


def test_single_variant_vcf_seq_extract(single_variant_vcf_seq_extractor):
    interval = Interval('chr1', 2, 9)
    seqs = single_variant_vcf_seq_extractor.extract(interval, anchor=3)
    assert next(seqs) == 'GCAACGT'
    assert next(seqs) == 'GTGAACG'


@pytest.fixture
def single_seq_vcf_seq_extractor():
    return SingleSeqVCFSeqExtractor(fasta_file, vcf_file)


def test_single_seq_vcf_seq_extract(single_seq_vcf_seq_extractor):
    interval = Interval('chr1', 2, 9)
    seq = single_seq_vcf_seq_extractor.extract(interval, anchor=3)
    assert seq == 'GCGAACG'

@pytest.fixture
def single_seq_vcf_seq_extractor_phased():
    return SingleSeqVCFSeqExtractor(fasta_file, phased_vcf_file)


def test_single_seq_vcf_seq_extract_phased(single_seq_vcf_seq_extractor_phased):
    interval = Interval('chr1', 24, 29)
    seq = single_seq_vcf_seq_extractor_phased.extract(interval, anchor=24, sample_id='NA00002', phase=1, fixed_len=True)
    assert seq == 'GATGC'
    
    interval2 = Interval('chr1', 45, 50)
    seq = single_seq_vcf_seq_extractor_phased.extract(interval2, anchor=50, sample_id='NA00002', phase=1, fixed_len=True)
    assert seq == 'TTGTA'
    
    seq_list = single_seq_vcf_seq_extractor_phased.batch_extract([interval, interval2], anchors=[24, 50], sample_id='NA00002', phase=1, fixed_len=True)
    assert seq_list == ['GATGC', 'TTGTA']
