'''
@A00744:46:HV3C3DSXX:2:2551:23493:1063
GACAGCGGGCACTCAGTTCCACTAACGAATGATGCCCGTGTGGACAGACAGAATGATGGACAGGCAGATGAATGCGTGGGCTTTATGTGAAACAGGTCCTCTTGGTTGTTGACAAGATACTGTTTTAAAGTTCCATTTTGCCATACTTCGA
+
FF,F::::,FFF::F:FFFFFF,F:F:F:F:FFF,FFFFFFFFFFFFFFFFFFFFF,::FFF::FFFFF:FFFFFF:F:FF,FF:F:FFFFF,FF:FFF:F:FFF:FFFF:FFFF::F:FFFFFF:F::FF:F,FFFFF::FFFFFFFF,F
'''

import argparse
import gzip
import random
import sys

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-r1', '--read1',
        help='Path to R1 FASTQ'
    )
    parser.add_argument(
        '-r2', '--read2',
        help='Path to R2 FASTQ'
    )
    parser.add_argument(
        '-op', '--out_prefix',
        help='Output FASTQ prefix'
    )
    parser.add_argument(
        '-s', '--sample_rate', type=float, default=0,
        help=('Sample rate (e.g. -s 0.05 samples ~5% of the reads). '
              'Should be within 0 and 1. [0].')
    )
    args = parser.parse_args()
    return args


def sample_fq(args):
    f1 = gzip.open(args.read1, 'rb')
    f2 = gzip.open(args.read2, 'rb')
    fo1 = gzip.open(args.out_prefix + '-R1.fq.gz', 'wb')
    fo2 = gzip.open(args.out_prefix + '-R2.fq.gz', 'wb')
    
    rd = 0
    for l1 in f1:
        # l2 = f2.readline().decode('ascii').rstrip()
        l2 = f2.readline()
        # l1 = l1.decode('ascii').rstrip()
        # if l1.startswith('@'):
        if l1.startswith(b'@'):
            rd = random.random()
        if rd <= 1 - args.sample_rate:
            # print(l1, file=sys.stderr)
            continue
        fo1.write(l1)
        fo2.write(l2)


if __name__ == '__main__':
    args = parse_args()
    if args.read1 is None:
        exit(1)
    sample_fq(args)

