
import argparse
import sys
import re


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-c', '--chain', required=True,
        help='Path to the input chain file.'
    )
    parser.add_argument(
        '-o', '--out', default='',
        help='Path to the output verbose chain file. [empty string]'
    )
    args = parser.parse_args()
    return args


def change_sides(in_fn: str, out_fn: str) -> None:
    f = open(in_fn, 'r')
    fo = open(out_fn, 'w')
    for line in f:
        if not line:
            continue
        line = line.rstrip()
        fields = re.split(r'[\s\t]+', line)
        if fields[0] == 'chain':
            source = fields[2]
            source_len = int(fields[3])
            s_start = int(fields[5])
            s_end = int(fields[6])
            dest = fields[7]
            dest_len = int(fields[8])
            strand = fields[9]
            d_start = int(fields[10])
            d_end = int(fields[11])
            changed = fields[:2] + fields[7:12] + fields[2:7] + [fields[12]]
            if strand == '-':
                changed[5] = str(dest_len - d_end)
                changed[6] = str(dest_len - d_start)
                changed[10] = str(source_len - s_end)
                changed[11] = str(source_len - s_start)
                changed[4] = '+'
                changed[9] = '-'
            changed_str = ' '.join(changed)
            print(f'{changed_str}', file=fo)
        elif len(fields) == 3:
            l = int(fields[0])
            ds = int(fields[1])
            dd = int(fields[2])
            changed = [fields[0], fields[2], fields[1]]
            changed_str = ' '.join(changed)
            print(f'{changed_str}', file=fo)
        elif len(fields) == 1 and fields[0] != '':
            l = int(fields[0])
            print(fields[0], file=fo)
            print('', file=fo)


if __name__ == '__main__':
    args = parse_args()

    print('Input chain:', args.chain, file=sys.stderr)
    print('Output chain (verbose):', args.out, file=sys.stderr)
    change_sides(in_fn=args.chain, out_fn=args.out)
