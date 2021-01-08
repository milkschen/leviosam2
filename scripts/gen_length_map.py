import argparse
import sys

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-f', '--faidx',
        help='Path to the .fa.idx file.'
    )
    parser.add_argument(
        '-o', '--out',
        help='Path to the output length map.'
    )
    args = parser.parse_args()
    return args


def gen_length_map(args):
    f = open(args.faidx, 'r')
    if not args.out:
        f_out = sys.stdout
    else:
        f_out = open(args.out, 'w')
    for line in f:
        line = line.split()
        print(f'{line[0]}\t{line[1]}', file=f_out)


if __name__ == '__main__':
    args = parse_args()
    gen_length_map(args)

