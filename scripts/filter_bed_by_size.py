'''
Filter a BED file by segment size.

Nae-Chyun Chen
Johns Hopkins University
2021
'''
import argparse
import pysam
import sys

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-b', '--input_bed', required=True,
        help='Path to the input mappability BED file. [required]'
    )
    parser.add_argument(
        '-o', '--out', default='',
        help='Path to the ouput report. ['': print to sys.stdout]'
    )
    parser.add_argument(
        '-s', '--size', type=int,
        help='Min allowed size of a BED record (smaller ones are removed). [required]'
    )
    args = parser.parse_args()
    return args


def filter_bed_by_size(args):
    fb = open(args.input_bed, 'r')
    if args.out == '':
        fo = sys.stdout
    else:
        fo = open(args.out, 'w')

    for line in fb:
        line = line.split()
        start = int(line[1])
        end = int(line[2])
        if end - start >= args.size:
            print(f'{line[0]}\t{start}\t{end}', file=fo)


if __name__ == '__main__':
    args = parse_args()
    filter_bed_by_size(args)
