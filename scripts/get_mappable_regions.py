'''
Reads a mappability BED file and returns highly-mappable regions

Nae-Chyun Chen
Johns Hopkins University
2021
'''
import argparse
import typing
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


def print_to_bed(
    chrom: chr, high_map_start: int, low_map_start: int, fo: typing.TextIO
) -> None:
    if high_map_start < low_map_start:
        print(f'{chrom}\t{high_map_start}\t{low_map_start}', file=fo)


def is_mappable(fields: list, min_mappability: float) -> bool:
    if fields[3] == 'inf' or float(fields[3]) < min_mappability:
        return False
    return True


def get_mappable_regions(args):
    fb = open(args.input_bed, 'r')
    if args.out == '':
        fo = sys.stdout
    else:
        fo = open(args.out, 'w')

    chrom = ''
    high_map_start = 0
    low_map_start = 0
    state = None
    for line in fb:
        line = line.split()
        # Reset when seeing a new chromosome
        if line[0] != chrom:
            chrom = line[0]
            high_map_start = 0
            low_map_start = 0
            state = None

        if not is_mappable(line, args.min_mappability):
            # mappable to unmappable
            if state != 'unmappable':
                low_map_start = int(line[1])
                print_to_bed(chrom, high_map_start, low_map_start, fo)
                state = 'unmappable'
            high_map_start = low_map_start + args.kmer_size
        else:
            if state != 'mappable':
                if int(line[1]) > high_map_start:
                    high_map_start = int(line[1])
                state = 'mappable'

    if state == 'mappable':
        print_to_bed(chrom, high_map_start, int(line[1])+args.kmer_size, fo)


if __name__ == '__main__':
    args = parse_args()
    get_mappable_regions(args)
