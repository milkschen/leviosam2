'''
Mask sequences in a FASTA with annotations in a BED file

Nae-Chyun Chen
Johns Hopkins University
2021
'''
import argparse
import leviosam_utils
import pysam
import sys

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
        help='Path to the ouput FASTA file. ['': print to sys.stdout]'
    )
    args = parser.parse_args()
    return args


def mask_fasta_with_bed(
    fasta: str, bed: str
) -> dict:
    ref = leviosam_utils.read_fasta(fasta)
    fb = open(bed, 'r')

    for line in fb:
        line = line.split()
        contig = line[0]
        start = int(line[1])
        end = int(line[2])
        size = end - start
        assert size >= 0
        ref[contig] = ref[contig][:start] + 'N' * size + ref[contig][end:]

    fb.close()
    return ref


def print_fasta(ref: dict, out: str) -> None:
    if out == '':
        fo = sys.stdout
    else:
        fo = open(out, 'w')

    for i, (k, v) in enumerate(ref.items()):
        print(f'>{k}', file=fo)
        for i in range(0, len(v), 60):
            print(f'{v[i:i+60]}', file=fo)
    fo.close()


if __name__ == '__main__':
    args = parse_args()
    ref = mask_fasta_with_bed(fasta=args.fasta, bed=args.bed)
    print_fasta(ref=ref, out=args.out)


