"""Compares two SAM/BAM files and report a summary.
"""
import argparse
import pysam
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
    args = parser.parse_args()
    return args


class Summary():
    def __init__(self, allowed_pos_diff):
        self.pos_diff = []
        self.mapq_diff = []
        self.cigar_diff = []
        self.records = []
        self.allowed_pos_diff = allowed_pos_diff
        # Unaligned or invalid alignments
        self.unmapped_records = [[], []]
        self.invalid_records = [[], []] # MAPQ = 255

    def update(self, query, baseline, aln_filter):
        try:
            is_unmapped_or_invalid = False
            if query.is_unmapped:
                self.unmapped_records[0].append(query.query_name + ('_1' if query.is_read1 else '_2'))
                is_unmapped_or_invalid = True
            if baseline.is_unmapped:
                self.unmapped_records[1].append(baseline.query_name + ('_1' if baseline.is_read1 else '_2'))
                is_unmapped_or_invalid = True
            if is_unmapped_or_invalid:
                return
            # Mapped but invalid
            if query.mapping_quality == 255 or query.mapping_quality < aln_filter['MAPQ']:
                self.invalid_records[0].append(query.query_name + ('_1' if query.is_read1 else '_2'))
                is_unmapped_or_invalid = True
            if baseline.mapping_quality == 255 or baseline.mapping_quality < aln_filter['MAPQ']:
                self.invalid_records[1].append(baseline.query_name + ('_1' if baseline.is_read1 else '_2'))
                is_unmapped_or_invalid = True
            if is_unmapped_or_invalid:
                return

            if query.reference_name == baseline.reference_name:
                self.pos_diff.append(abs(query.reference_start - baseline.reference_start))
            else:
                self.pos_diff.append(-1)
            self.records.append([query, baseline])
            self.mapq_diff.append(query.mapping_quality - baseline.mapping_quality)
            self.cigar_diff.append(query.cigarstring == baseline.cigarstring)
        except:
            print('query=', query.query_name, query.reference_name, query.reference_start)
            print('baseline=', baseline.query_name, baseline.reference_name, baseline.reference_start)

    def report(self, fn_out, num_err_printed):
        if fn_out == '':
            f_out = sys.stdout
        else:
            f_out = open(fn_out, 'w')
        print('## Position', file=f_out)
        num_pos_match = sum([i >= 0 and i < self.allowed_pos_diff for i in self.pos_diff])
        print(f'{num_pos_match / len(self.pos_diff)} ({num_pos_match}/{len(self.pos_diff)})', file=f_out)
        cnt = 0
        for i in range(len(self.records)):
            # if self.pos_diff[i] < 100 and self.pos_diff[i] > self.allowed_pos_diff:
            if self.pos_diff[i] > self.allowed_pos_diff:
                query = self.records[i][0]
                baseline = self.records[i][1]
                msg_query = (
                    '    '
                    f'{query.reference_name}\t{query.reference_start:10d}\t'
                    f'{query.mapping_quality:3d}\t{query.cigarstring}')
                msg_baseline = (
                    '    '
                    f'{baseline.reference_name}\t{baseline.reference_start:10d}\t'
                    f'{baseline.mapping_quality:3d}\t{baseline.cigarstring}')
                print(f'{self.pos_diff[i]:10d}\t{msg_query}')
                print(f'{" ":10s}\t{msg_baseline}')
                cnt += 1
            if cnt >= num_err_printed and num_err_printed >= 0:
                break

        print('## MAPQ', file=f_out)
        print(f'{self.mapq_diff.count(0) / len(self.mapq_diff)} ({self.mapq_diff.count(0)}/{len(self.mapq_diff)})', file=f_out)

        print('## CIGAR', file=f_out)
        print(f'{self.cigar_diff.count(True) / len(self.cigar_diff)} ({self.cigar_diff.count(True)}/{len(self.cigar_diff)})', file=f_out)

        print('## Unaligned', file=f_out)
        set_unaligned = set(self.unmapped_records[0] + self.unmapped_records[1])
        print((f'{len(set_unaligned)} (query={len(self.unmapped_records[0])}, '
               f'baseline={len(self.unmapped_records[1])})'), file=f_out)

        print('## Invalid (MAPQ=255)', file=f_out)
        set_invalid = set(self.invalid_records[0] + self.invalid_records[1])
        print((f'{len(set_invalid)} (query={len(self.invalid_records[0])}, '
               f'baseline={len(self.invalid_records[1])})'), file=f_out)

def read_sam_as_dict(fn):
    dict_reads = {}
    f = pysam.AlignmentFile(fn, 'r')
    for read in f.fetch():
        if (not read.is_paired) or (read.is_paired and read.is_read1):
            segment_idx = 0
        elif read.is_paired and read.is_read2:
            segment_idx = 1

        if dict_reads.setdefault(read.query_name, [[], []]):
            #dict_reads[read.query_name][segment_idx] = read_info
            dict_reads[read.query_name][segment_idx] = read
    return dict_reads


def compare_sam(args):
    summary = Summary(args.allowed_pos_diff)
    dict_query = read_sam_as_dict(args.input_query)
    dict_baseline = read_sam_as_dict(args.input_baseline)
    aln_filter = {'MAPQ': args.min_mapq}
    for i_q, [name, [first_seg, second_seg]] in enumerate(dict_query.items()):
        if dict_baseline.get(name):
            summary.update(
                query=first_seg, baseline=dict_baseline[name][0], aln_filter=aln_filter)
            summary.update(
                query=second_seg, baseline=dict_baseline[name][1], aln_filter=aln_filter)
    summary.report(args.out, args.num_err_printed)


if __name__ == '__main__':
    args = parse_args()
    compare_sam(args)


