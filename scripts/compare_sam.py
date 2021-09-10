'''
Compares two SAM/BAM files and report a summary.

Nae-Chyun Chen
Johns Hopkins University
2021
'''
import argparse
import pysam
import re
import sys

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-q', '--input_query',
        help='Path to a query SAM/BAM file.'
    )
    parser.add_argument(
        '-b', '--input_baseline',
        help='Path to a baseline SAM/BAM file.'
    )
    parser.add_argument(
        '-o', '--out', default='',
        help='Path to the ouput report. ['': print to sys.stdout]'
    )
    parser.add_argument(
        '-c', '--num_err_printed', default=10, type=int,
        help='Number of printed errors. Set to -1 to print all. [10]'
    )
    parser.add_argument(
        '-cp', '--categories_printed', default='pos_idy',
        help=('Caterogies that erroneous records are printed. '
              'Options: "pos", "idy", "pos_idy", "cigar", "tlen". '
              'Set to "none" to turn this off. '
              'Split by commas, e.g. `--categories_printed pos, idy`. '
              '["pos_idy"]'
        )
    )
    parser.add_argument(
        '-p', '--allowed_posdiff', default=1, type=int,
        help='Allowed bases of positional difference. [1]'
    )
    parser.add_argument(
        '-i', '--identity_cutoff', default=0.8, type=float,
        help=(
            'Identify cutoff. An alignment with an identity >= this value '
            'is considered as correct. [0.8]')
    )
    parser.add_argument(
        '-m', '--min_mapq', default=0, type=int,
        help='Min MAPQ to consider. Alignments with lower MAPQ are considered as invalid. [0]'
    )
    parser.add_argument(
        '-mr', '--max_posdiff_reported', default=sys.maxsize, type=int,
        help=(
            'Only report records within this value of positional difference. '
            'Set to a negative value to turn off. [sys.maxsize]')
    )
    args = parser.parse_args()
    return args


class SamUtils:
    ''' Convert a typical CIGAR string to the expanded form
    Example:
        s = CompareSamSummary(0)
        print(s._expand_cigar('3M2D5M'))
        # MMMDDMMMMM
    '''
    @staticmethod
    def _expand_cigar(cigarstring: str, consumes: str) -> str:
        re_cigar = re.compile('[MIDS]+')
        ops = re_cigar.findall(cigarstring)
        lens = re_cigar.split(cigarstring)
        ecigar = ''
        for i, op in enumerate(ops):
            if consumes == 'ref' and not (op in ['M', 'D']):
                continue
            elif consumes == 'query' and not (op in ['M', 'I', 'S']):
                continue
            ecigar += int(lens[i]) * op
        return ecigar

    ''' Return true if a cigar operator consumes QUERY 
    Op table:
        0 M BAM_CMATCH
        1 I BAM_CINS
        2 D BAM_CDEL
        3 N BAM_CREF_SKIP
        4 S BAM_CSOFT_CLIP
        5 H BAM_CHARD_CLIP
        6 P BAM_CPAD
        7 = BAM_CEQUAL
        8 X BAM_CDIFF
        9 B BAM_CBACK
    '''
    @staticmethod
    def cigar_op_consumes_query(op: int) -> bool:
        return op in [0, 1, 4, 7, 8]

    ''' Return true if a cigar operator consumes REFERENCE '''
    @staticmethod
    def cigar_op_consumes_ref(op: int) -> bool:
        return op in [0, 2, 3, 7, 8]


class CigarSegment:
    def __init__(self, cigartuple, q_offset, r_offset) -> None:
        self.op = cigartuple[0]
        self.size = cigartuple[1]
        self.r_offset = r_offset
        self.q_offset = q_offset
        self.intercept = r_offset - q_offset

    def __repr__(self) -> str:
        return f'{self.op}\t{self.size}\t{self.q_offset:<10d}\t{self.intercept}'
    
    def truncate_front(self, y) -> int:
        assert y.q_offset > self.q_offset
        diff = y.q_offset - self.q_offset
        self.size = self.size - diff
        self.q_offset += diff
        self.r_offset += diff

    def overlap(self, y) -> int:
        assert (self.q_offset == y.q_offset)
        if self.op == y.op and self.intercept == y.intercept:
            return min(self.size, y.size)
        else:
            return 0


class CigarSegments:
    def __init__(self, aln) -> None:
        stats = []
        r_offset = aln.reference_start
        q_offset = 0
        for i, (op, n) in enumerate(aln.cigartuples):
            # cigartuples: e.g [(0, 30), (1, 5), (0, 20)] -> 30M5I20M
            stats.append(CigarSegment((op, n), q_offset, r_offset))
            if SamUtils.cigar_op_consumes_query(op):
                q_offset += n
            if SamUtils.cigar_op_consumes_ref(op):
                r_offset += n
        self.cigar = stats

    def __repr__(self) -> str:
        out = 'OP\tSIZE\tQ_OFFSET\tINTETCEPT\n'
        for s in self.cigar:
            out += (s.__repr__() + '\n')
        return out

    def __getitem__(self, key):
        return self.cigar[key]

    def __len__(self):
        return len(self.cigar)

    def pop(self, key) -> None:
        self.cigar.pop(key)


class CompareSamSummary():
    '''
    Inputs:
        - allowed_posdiff: allowed bases of positional difference
        - num_err_printed: number of printed errors
        - max_posdiff_reported: only report records within this value of positional difference
        - fn_out: path to the ouput report
        - aln_filter: alignment filtering criteria
    '''
    def __init__(
        self, allowed_posdiff: int=1, num_err_printed: int=10,
        max_posdiff_reported: int=sys.maxsize, identity_cutoff: float=0.8,
        fn_out: str='', aln_filter: dict={'MAPQ': 0}) -> None:
        self.posdiff = []
        self.mapq_diff = []
        self.cigar_diff = []
        self.records = []
        # Unaligned or invalid alignments
        self.unmapped_records = [[], []]
        self.invalid_records = [[], []] # MAPQ = 255
        self.identity = []
        self.tlendiff = []

        self.allowed_posdiff = allowed_posdiff
        self.identity_cutoff = identity_cutoff
        self.num_err_printed = num_err_printed
        self.max_posdiff_reported = max_posdiff_reported
        self.fn_out = fn_out
        self.aln_filter = aln_filter

    ''' Caculate CIGAR identity between two sequences
    Core algorithmic idea:
        In each iteration, compare the query_offsets between two sequences.
        Truncate the sequence with a smaller query_offset (there can not be
        a match in this scenario).
        When the query_offsets are equal, pop the smaller CIGAR segment. 
        We also check if the CIGAR op and intercept are equal. If equal, 
        increment `idy` with the size of the popped segment.

    This algorithm can be executed in O(len(rx)+len(ry)) time, where rx and ry
    are the number of CIGAR runs for the input sequences respectively.
    '''
    def _calc_identity(
        self, query: pysam.AlignedSegment, baseline: pysam.AlignedSegment
    ) -> float:
        idy = 0
        posdiff = self._calc_posdiff(query, baseline)
        if posdiff < 0:
            return 0
        elif posdiff == 0 and query.cigarstring == baseline.cigarstring:
            # Quick return for nicely-matched query-baseline pairs
            return 1

        query_cigar_offsets = CigarSegments(query)
        baseline_cigar_offsets = CigarSegments(baseline)
        while len(query_cigar_offsets) and len(baseline_cigar_offsets):
            b = baseline_cigar_offsets[0]
            q = query_cigar_offsets[0]
            if q.q_offset == b.q_offset:
                idy += q.overlap(b)
                if q.size == b.size:
                    baseline_cigar_offsets.pop(0)
                    query_cigar_offsets.pop(0)
                elif q.size < b.size:
                    query_cigar_offsets.pop(0)
                else:
                    baseline_cigar_offsets.pop(0)
            elif q.q_offset < b.q_offset:
                q.truncate_front(b)
            else:
                b.truncate_front(q)
        
        return idy / query.infer_read_length()

    '''
    Calculate the difference of the leftmost aligned positions between query and baseline.
    If query and baseline are aligned to different contigs, return -1.
    Otherwise results are always >= 0.
    '''
    def _calc_posdiff(self, query: pysam.AlignedSegment, baseline: pysam.AlignedSegment) -> int:
        if query.reference_name == baseline.reference_name:
            return abs(query.reference_start - baseline.reference_start)
        else:
            return -1
    
    '''
    Calculate the difference of TLEN between query and baseline.
    If query and baseline are aligned to different contigs, return -1.
    Otherwise results are always >= 0.
    '''
    def _calc_tlendiff(self, query: pysam.AlignedSegment, baseline: pysam.AlignedSegment) -> int:
        if query.reference_name == baseline.reference_name:
            return abs(query.template_length - baseline.template_length)
        else:
            return -1

    '''
    Check if either of query or baseline is unmapped. For an unmapped alignment, add it to 
    `self.unmapped_records`. If any alignment in the pair is unmapped, return false; return 
    true otherwise.
    '''
    def _check_unmap(self, query: pysam.AlignedSegment, baseline: pysam.AlignedSegment) -> bool:
        is_unmapped = False
        if query.is_unmapped:
            self.unmapped_records[0].append(
                query.query_name + ('_1' if query.is_read1 else '_2'))
            is_unmapped = True
        if baseline.is_unmapped:
            self.unmapped_records[1].append(
                baseline.query_name + ('_1' if baseline.is_read1 else '_2'))
            is_unmapped = True
        return is_unmapped

    '''
    Check if either of query or baseline has low MAPQ. For an alignment not passed the filter,
    add it to `self.invalid_records`. If any alignment in the pair does not pass, return false;
    return true otherwise.
    '''
    def _check_low_qual(
        self, query: pysam.AlignedSegment, baseline: pysam.AlignedSegment
    ) -> bool:
        is_low_qual = False
        if query.mapping_quality == 255 or query.mapping_quality < self.aln_filter['MAPQ']:
            self.invalid_records[0].append(
                query.query_name + ('_1' if query.is_read1 else '_2'))
            is_low_qual = True
        if baseline.mapping_quality == 255 or baseline.mapping_quality < self.aln_filter['MAPQ']:
            self.invalid_records[1].append(
                baseline.query_name + ('_1' if baseline.is_read1 else '_2'))
            is_low_qual = True
        return is_low_qual

    '''
    Read one record from query and baseline and update the summary
    
    Inputs:
      - query/baseline (pysam alignment record)
    '''
    def update(
        self, query: pysam.AlignedSegment, baseline: pysam.AlignedSegment
    ) -> None:
        if (not query) or (not baseline):
            return

        # Whether to filter the baseline dataset is a bit tricky - on one hand, we should
        # not filter the "truth"; on the other, there can be failed heuristics, especially for
        # unmapped and/or low-mapping-quality alignments.
        # Exlude unmapped alignments
        if self._check_unmap(query, baseline):
            return
        # Mapped but invalid
        if self._check_low_qual(query, baseline):
            return
        
        self.posdiff.append(self._calc_posdiff(query, baseline))
        self.records.append([query, baseline])
        self.mapq_diff.append(query.mapping_quality - baseline.mapping_quality)
        self.cigar_diff.append(query.cigarstring == baseline.cigarstring)
        self.identity.append(self._calc_identity(query, baseline))
        self.tlendiff.append(self._calc_tlendiff(query, baseline))

    ''' Report summary results '''
    def _print_records(self, f_out, by='pos') -> None:
        cnt = 0
        show = False
        for i, rec in enumerate(self.records):
            if cnt >= self.num_err_printed:
                break
            if by == 'pos':
                show = (self.posdiff[i] > self.allowed_posdiff and 
                        self.posdiff[i] <= self.max_posdiff_reported)
            elif by == 'idy':
                show = (self.identity[i] < self.identity_cutoff)
            elif by == 'pos_idy':
                show = (self.posdiff[i] > self.allowed_posdiff and 
                        self.posdiff[i] <= self.max_posdiff_reported)
                show &= (self.identity[i] < self.identity_cutoff)
            elif by == 'cigar':
                show = not self.cigar_diff[i]
            elif by == 'tlen':
                show = (self.tlendiff[i] >= 50 and self.tlendiff[i] >= 0)

            if show:
                query = rec[0]
                baseline = rec[1]
                msg_query = (
                    '    '
                    f'{query.flag:5d}\t{query.reference_name:6s}\t'
                    f'{query.reference_start+1:10d}\t'
                    f'{query.mapping_quality:3d}\t{query.cigarstring}\t'
                    f'{query.template_length}')
                msg_baseline = (
                    '    '
                    f'{baseline.flag:5d}\t{baseline.reference_name:6s}\t'
                    f'{baseline.reference_start+1:10d}\t'
                    f'{baseline.mapping_quality:3d}\t{baseline.cigarstring}\t'
                    f'{baseline.template_length}')
                print(f'{query.query_name}', file=f_out)
                print(f'  p_diff = {self.posdiff[i]:<10d}\t{msg_query}', file=f_out)
                print(f'  idy    = {self.identity[i]:.4f}\t{msg_baseline}', file=f_out)
                cnt += 1
        return

    def report(self, cat_printed: list=[]) -> None:
        if self.fn_out == '':
            f_out = sys.stdout
        else:
            f_out = open(self.fn_out, 'w')

        if len(self.posdiff) == 0:
            print('Zero matched records. Exit.', file=f_out)
            exit(1)

        print('## Position', file=f_out)
        num_pos_match = sum([i >= 0 and i < self.allowed_posdiff for i in self.posdiff])
        print(f'{num_pos_match / len(self.posdiff):.6f} '
              f'({num_pos_match}/{len(self.posdiff)})',
              file=f_out)
        if self.num_err_printed > 1 and 'pos' in cat_printed:
            self._print_records(f_out, by='pos')

        print('## Identity', file=f_out)
        num_idy = sum([i >= self.identity_cutoff for i in self.identity])
        print(f'{num_idy / len(self.identity):.6f} ({num_idy}/{len(self.identity)})',
              file=f_out)
        if self.num_err_printed > 1 and 'idy' in cat_printed:
            self._print_records(f_out, by='idy')
        print('## Average Identity', file=f_out)
        avg_idy = sum([i for i in self.identity]) / len(self.identity)
        print(f'{avg_idy:.6f}', file=f_out)

        print('## Position || Identity', file=f_out)
        num_pos_idy = sum([d >= self.identity_cutoff or (self.posdiff[i] >= 0 and self.posdiff[i] < self.allowed_posdiff) for i, d in enumerate(self.identity)])
        print(f'{num_pos_idy / len(self.identity):.6f} ({num_pos_idy}/{len(self.identity)})',
              file=f_out)
        if self.num_err_printed > 1 and 'pos_idy' in cat_printed:
            self._print_records(f_out, by='pos_idy')

        print('## MAPQ', file=f_out)
        print((f'{self.mapq_diff.count(0) / len(self.mapq_diff):.6f} '
               f'({self.mapq_diff.count(0)}/{len(self.mapq_diff)})'), file=f_out)

        print('## CIGAR', file=f_out)
        print((f'{self.cigar_diff.count(True) / len(self.cigar_diff):.6f} '
               f'({self.cigar_diff.count(True)}/{len(self.cigar_diff)})'), file=f_out)
        if self.num_err_printed > 1 and 'cigar' in cat_printed:
            self._print_records(f_out, by='cigar')
        
        print('## TLEN', file=f_out)
        num_tlen_match = sum([i <= 50 and i >= 0 for i in self.tlendiff])
        print(f'{num_tlen_match / len(self.tlendiff):.6f} '
              f'({num_tlen_match}/{len(self.tlendiff)})',
              file=f_out)
        if self.num_err_printed > 1 and 'tlen' in cat_printed:
            self._print_records(f_out, by='tlen')

        print('## Unaligned', file=f_out)
        set_unaligned = set(self.unmapped_records[0] + self.unmapped_records[1])
        print((f'{len(set_unaligned)} (query={len(self.unmapped_records[0])}, '
               f'baseline={len(self.unmapped_records[1])})'), file=f_out)

        print('## Invalid (MAPQ=255)', file=f_out)
        set_invalid = set(self.invalid_records[0] + self.invalid_records[1])
        print((f'{len(set_invalid)} (query={len(self.invalid_records[0])}, '
               f'baseline={len(self.invalid_records[1])})'), file=f_out)


''' Read a SAM/BAM file as a dictionary (key: query name; value: pysam alignment record) '''
def read_sam_as_dict(fn: str) -> dict:
    dict_reads = {}
    f = pysam.AlignmentFile(fn, 'r')
    for read in f:
        if not read.is_paired:
            segment_idx = 0
        else:
            if read.is_read1:
                segment_idx = 0
            else:
                segment_idx = 1

        if dict_reads.setdefault(read.query_name, [[], []]):
            dict_reads[read.query_name][segment_idx] = read
    return dict_reads


def compare_sam(args):
    aln_filter = {'MAPQ': args.min_mapq}
    summary = CompareSamSummary(
        allowed_posdiff=args.allowed_posdiff,
        identity_cutoff=args.identity_cutoff,
        num_err_printed=args.num_err_printed,
        max_posdiff_reported=args.max_posdiff_reported,
        fn_out=args.out, aln_filter=aln_filter)
    dict_query = read_sam_as_dict(args.input_query)
    dict_baseline = read_sam_as_dict(args.input_baseline)
    for _, [name, [first_seg, second_seg]] in enumerate(dict_query.items()):
        if dict_baseline.get(name):
            summary.update(query=first_seg, baseline=dict_baseline[name][0])
            summary.update(query=second_seg, baseline=dict_baseline[name][1])
    if args.categories_printed == 'none':
        cat_printed = []
    else:
        cat_printed = args.categories_printed.split(',')
        for c in cat_printed:
            assert c in ['pos', 'idy', 'pos_idy', 'cigar', 'tlen']
    summary.report(cat_printed)


if __name__ == '__main__':
    args = parse_args()
    compare_sam(args)

