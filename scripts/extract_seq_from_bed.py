'''
Extract subsequences from a FASTA file using info from a BED

Nae-Chyun Chen
Johns Hopkins University
2021
'''
import argparse
import pysam
import sys
import leviosam_utils

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-f', '--fasta', required=True,
        help='Path to the input FASTA file. [required]'
    )
    parser.add_argument(
        '-b', '--bed', required=True,
        help='Path to the input BED file. [required]'
    )
    parser.add_argument(
        '-o', '--out', default='',
        help='Path to the extracted subsequences (FASTA file). ['': print to sys.stdout]'
    )
    args = parser.parse_args()
    return args


def extract_seq_from_bed(args):
    ref = leviosam_utils.read_fasta(args.fasta)
    if args.out == '':
        fo = sys.stdout
    else:
        fo = open(args.out, 'w')

    fb = open(args.bed, 'r')
    for line in fb:
        line = line.split('\t')
        contig = line[0]
        start = int(line[1])
        end = int(line[2])
        # Note that FASTA coords are usually represented in 1-based;
        # other systems here are 0-based
        print(f'>{contig}:{start+1}-{end+1}', file=fo)
        print(f'{ref[contig][start:end]}', file=fo)


if __name__ == '__main__':
    args = parse_args()
    extract_seq_from_bed(args)


