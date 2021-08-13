"""Compares two SAM/BAM files and report a summary.
"""
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
        '-c', '--num_err_printed', default=20, type=int,
        help='Number of printed errors. Set to -1 to print all. [20]'
    )
    parser.add_argument(
        '-p', '--allowed_pos_diff', default=1, type=int,
        help='Allowed bases of positional difference. [1]'
    )
    parser.add_argument(
        '-m', '--min_mapq', default=0, type=int,
        help='Min MAPQ to consider. Alignments with lower MAPQ are considered as invalid. [0]'
    )
    parser.add_argument(
        '--max_posdiff_reported', default=-1, type=int,
        help=(
            'Only report records within a positional difference. '
            'Set to a negative value to turn off. [-1]')
    )
    args = parser.parse_args()
    return args


class CompareSamSummary():
    def __init__(self, allowed_pos_diff):
        self.pos_diff = []
        self.mapq_diff = []
        self.cigar_diff = []
        self.records = []
        self.allowed_pos_diff = allowed_pos_diff
        # Unaligned or invalid alignments
        self.unmapped_records = [[], []]
        self.invalid_records = [[], []] # MAPQ = 255
        self.identity = []

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


    #TODO this is currently incorrect
    ''' Caculate CIGAR identity between two sequences

    Note that I didn't implment it with a fast algorithm. There should be some way to speed this
    up.
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

        # query_ecigar = self._expand_cigar(cigarstring=query.cigarstring, consumes='query')
        # baseline_ecigar = self._expand_cigar(cigarstring=baseline.cigarstring, consumes='query')
        # query_ecigar = self._expand_cigar(cigarstring=query.cigarstring, consumes='ref')
        # baseline_ecigar = self._expand_cigar(cigarstring=baseline.cigarstring, consumes='ref')
        query_ecigar = self._expand_cigar(cigarstring=query.cigarstring, consumes=None)
        baseline_ecigar = self._expand_cigar(cigarstring=baseline.cigarstring, consumes=None)

        def _get_cigar_offsets(aln):
            stats = aln.cigartuples
            # e.g [(0, 30), (1, 5), (0, 20)] -> 30M5I20M
            r_offset = aln.reference_start
            q_offset = 0
            for i, (op, n) in enumerate(stats):
                stats[i] = (op, n, q_offset, r_offset, r_offset - q_offset)
                if op in [0, 1, 4, 7, 8]:
                    q_offset += n
                if op in [0, 2, 3, 7, 8]:
                    r_offset += n
            return stats
        query_cigar_offsets = _get_cigar_offsets(query)
        baseline_cigar_offsets = _get_cigar_offsets(baseline)
        while query_cigar_offsets and baseline_cigar_offsets:
            for i
            if baseline_cigar_offsets[0][3] < query_cigar_offsets[0][3]:
                
        print(query)
        print(query_cigar_offsets)
        print(baseline)
        print(baseline_cigar_offsets)
        exit(0)


        # diff = query.reference_start - baseline.reference_start
        # if diff > 0:
        #     baseline_ecigar = baseline_ecigar[diff:]
        # elif diff < 0:
        #     query_ecigar = query_ecigar[diff:]
        # idxb = 0
        # idxq = 0
        # for i, q in enumerate(query_ecigar):
        #     if i >= len(baseline_ecigar):
        #         break
        #     if q not in ['M', 'I', 'S']:
        #         idxq -= 1
        #     if baseline_ecigar[i] not in ['M', 'I', 'S']:
        #         idxb -= 1
        #     idxq += 1
        #     idxb += 1
        #     if q == baseline_ecigar[idxb]:
        #         idy += 1
        
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
        self, query: pysam.AlignedSegment,
        baseline: pysam.AlignedSegment, aln_filter: dict
    ) -> bool:
        is_low_qual = False
        if query.mapping_quality == 255 or query.mapping_quality < aln_filter['MAPQ']:
            self.invalid_records[0].append(
                query.query_name + ('_1' if query.is_read1 else '_2'))
            is_low_qual = True
        if baseline.mapping_quality == 255 or baseline.mapping_quality < aln_filter['MAPQ']:
            self.invalid_records[1].append(
                baseline.query_name + ('_1' if baseline.is_read1 else '_2'))
            is_low_qual = True
        return is_low_qual

    '''
    Read one record from query and baseline and update the summary
    
    Inputs:
      - query/baseline (pysam alignment record)
      - aln_filter (dict): alignment filtering criteria
    '''
    def update(
        self, query: pysam.AlignedSegment,
        baseline: pysam.AlignedSegment, aln_filter: dict
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
        if self._check_low_qual(query, baseline, aln_filter):
            return
        
        self.pos_diff.append(self._calc_posdiff(query, baseline))
        self.records.append([query, baseline])
        self.mapq_diff.append(query.mapping_quality - baseline.mapping_quality)
        self.cigar_diff.append(query.cigarstring == baseline.cigarstring)
        self.identity.append(self._calc_identity(query, baseline))

    ''' Report summary results '''
    def report(self, fn_out: str, num_err_printed: int, max_posdiff_reported: int) -> None:
        if fn_out == '':
            f_out = sys.stdout
        else:
            f_out = open(fn_out, 'w')

        if len(self.pos_diff) == 0:
            print('Zero matched records. Exit.', file=f_out)
            exit(1)

        print('## Position', file=f_out)
        num_pos_match = sum([i >= 0 and i < self.allowed_pos_diff for i in self.pos_diff])
        print(
            f'{num_pos_match / len(self.pos_diff)} ({num_pos_match}/{len(self.pos_diff)})',
            file=f_out)
        cnt = 0
        for i, rec in enumerate(self.records):
            if self.pos_diff[i] > self.allowed_pos_diff:
                if max_posdiff_reported >= 0 and self.pos_diff[i] > max_posdiff_reported:
                    continue
                query = rec[0]
                baseline = rec[1]
                msg_query = (
                    '    '
                    f'{query.flag:5d}\t{query.reference_name:6s}\t'
                    f'{query.reference_start+1:10d}\t'
                    f'{query.mapping_quality:3d}\t{query.cigarstring}')
                msg_baseline = (
                    '    '
                    f'{baseline.flag:5d}\t{baseline.reference_name:6s}\t'
                    f'{baseline.reference_start+1:10d}\t'
                    f'{baseline.mapping_quality:3d}\t{baseline.cigarstring}')
                print(f'{query.query_name}', file=f_out)
                print(f'  p_diff = {self.pos_diff[i]:<10d}\t{msg_query}', file=f_out)
                # print(f'{" ":17s}\t{msg_baseline}', file=f_out)
                print(f'  idy    = {self.identity[i]:.4f}\t{msg_baseline}', file=f_out)
                cnt += 1
            if cnt >= num_err_printed and num_err_printed >= 0:
                break

        print('## Identity', file=f_out)
        num_idy = sum([i >= 0.8 for i in self.identity])
        print(
            f'{num_idy / len(self.identity)} ({num_idy}/{len(self.identity)})',
            file=f_out)

        print('## MAPQ', file=f_out)
        print((f'{self.mapq_diff.count(0) / len(self.mapq_diff)} '
               f'({self.mapq_diff.count(0)}/{len(self.mapq_diff)})'), file=f_out)

        print('## CIGAR', file=f_out)
        print((f'{self.cigar_diff.count(True) / len(self.cigar_diff)} '
               f'({self.cigar_diff.count(True)}/{len(self.cigar_diff)})'), file=f_out)

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
            segmend_idx = 0
        else:
            if read.is_read1:
                segment_idx = 0
            else:
                segment_idx = 1

        if dict_reads.setdefault(read.query_name, [[], []]):
            dict_reads[read.query_name][segment_idx] = read
    return dict_reads


def compare_sam(args):
    summary = CompareSamSummary(args.allowed_pos_diff)
    dict_query = read_sam_as_dict(args.input_query)
    dict_baseline = read_sam_as_dict(args.input_baseline)
    aln_filter = {'MAPQ': args.min_mapq}
    for i_q, [name, [first_seg, second_seg]] in enumerate(dict_query.items()):
        if dict_baseline.get(name):
            summary.update(
                query=first_seg, baseline=dict_baseline[name][0],
                aln_filter=aln_filter)
            summary.update(
                query=second_seg, baseline=dict_baseline[name][1],
                aln_filter=aln_filter)
    summary.report(
        fn_out=args.out, num_err_printed=args.num_err_printed,
        max_posdiff_reported=args.max_posdiff_reported)


if __name__ == '__main__':
    args = parse_args()
    compare_sam(args)

