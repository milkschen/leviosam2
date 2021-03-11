"""Compares two SAM files and report a summary.
"""
import argparse
import sys

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-q', '--input_query',
        help='Path to a query SAM file.'
    )
    parser.add_argument(
        '-b', '--input_baseline',
        help='Path to a baseline SAM file.'
    )
    parser.add_argument(
        '-o', '--out', default=sys.stdout,
        help='Path to the ouput report. [stdout]'
    )
    parser.add_argument(
        '-c', '--num_err_printed', default=20, type=int,
        help='Number of printed errors. Set to -1 to print all. [20]'
    )
    parser.add_argument(
        '-p', '--allowed_pos_diff', default=1, type=int,
        help='Allowed bases of positional difference. [1]'
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

    def update(self, query, baseline):
        if len(query) == 0 or len(baseline) == 0:
            return
        if query[0] == baseline[0]: # contig
            self.pos_diff.append(abs(query[1] - baseline[1]))
        else:
            self.pos_diff.append(-1)
        self.records.append([query, baseline])
        self.mapq_diff.append(query[2] - baseline[2])
        self.cigar_diff.append(query[3] == baseline[3])

    def report(self, f_out, num_err_printed):
        print('## Position', file=f_out)
        num_pos_match = sum([i >= 0 and i < self.allowed_pos_diff for i in self.pos_diff])
        print(f'{num_pos_match / len(self.pos_diff)} ({num_pos_match}/{len(self.pos_diff)})', file=f_out)
        #print(f'{self.pos_diff.count(0) / len(self.pos_diff)} ({self.pos_diff.count(0)}/{len(self.pos_diff)})', file=f_out)
        cnt = 0
        for i in range(len(self.records)):
            if self.pos_diff[i] < 100 and self.pos_diff[i] > self.allowed_pos_diff:
                print(self.pos_diff[i], self.records[i], file=f_out)
                cnt += 1
            if cnt >= num_err_printed and num_err_printed >= 0:
                break

        print('## MAPQ', file=f_out)
        print(f'{self.mapq_diff.count(0) / len(self.mapq_diff)} ({self.mapq_diff.count(0)}/{len(self.mapq_diff)})', file=f_out)

        print('## CIGAR', file=f_out)
        print(f'{self.cigar_diff.count(True) / len(self.cigar_diff)} ({self.cigar_diff.count(True)}/{len(self.cigar_diff)})', file=f_out)


def process_sam_line(line, dict_sam):
    line = line.split()
    name = line[0]
    flag = int(line[1])
    contig = line[2]
    pos = int(line[3])
    mapq = int(line[4])
    cigar = line[5]

    # Ignore unmapped reads.
    if flag & 4:
        return

    if (flag & 1 and flag & 64) or (not (flag & 1)):
        if dict_sam.setdefault(name, [[], []]):
            dict_sam[name][0] = [contig, pos, mapq, cigar]
    elif (flag & 1 and flag & 128):
        if dict_sam.setdefault(name, [[], []]):
            dict_sam[name][1] = [contig, pos, mapq, cigar]
    else:
        print(line)


def read_sam_as_dict(fn):
    dict_sam = {}
    with open(fn, 'r') as f:
        for line in f:
            if line[0] != '@':
                process_sam_line(line, dict_sam)
    return dict_sam


def compare_sam(args):
    summary = Summary(args.allowed_pos_diff)
    dict_query = read_sam_as_dict(args.input_query)
    dict_baseline = read_sam_as_dict(args.input_baseline)
    for i_q, [name, [first_seg, second_seg]] in enumerate(dict_query.items()):
        if dict_baseline.get(name):
            summary.update(first_seg, dict_baseline[name][0])
            summary.update(second_seg, dict_baseline[name][1])
    summary.report(args.out, args.num_err_printed)


if __name__ == '__main__':
    args = parse_args()
    compare_sam(args)


