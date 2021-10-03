'''
Reads a mappability BED file and returns highly-mappable regions.

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
        '-c', '--min_mappability', default=1.0, type=float,
        help='Min mappability to be considered as highly-mappable. [1.0]'
    )
    parser.add_argument(
        '-k', '--kmer_size', required=True, type=int,
        help='Kmer size for each BED record. [required]'
    )
    args = parser.parse_args()
    return args


def get_mappable_regions(args):
    fb = open(args.input_bed, 'r')
    if args.out == '':
        fo = sys.stdout
    else:
        fo = open(args.out, 'w')

    chrom = ''
    high_map_start = 0
    high_map_end = 0
    for line in fb:
        line = line.split()
        # Reset if sees a new chromosome
        if line[0] != chrom:
            chrom = line[0]
            high_map_start = 0
            high_map_end = 0
        if line[3] == 'inf' or float(line[3]) >= args.min_mappability:
            high_map_end = int(line[1]) + args.kmer_size
        else:
            low_map_start = int(line[1])
            if high_map_start < low_map_start:
                print(f'{chrom}\t{high_map_start}\t{low_map_start}', file=fo)
                # print(f'{chrom}\t{high_map_start}\t{high_map_end}', file=fo)
            high_map_start = low_map_start + args.kmer_size



if __name__ == '__main__':
    args = parse_args()
    get_mappable_regions(args)
