''' Add info for a chain file to facilitate debugging.
'''
import argparse
import re
import sys

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-c', '--chain',
        help='Path to the input chain file.'
    )
    parser.add_argument(
        '-o', '--out',
        help='Path to the output verbose chain file.'
    )
    args = parser.parse_args()
    return args


def verbosify_chain(fn_chain, fn_out=None):
    f = open(fn_chain, 'r')
    if not fn_out:
        fo = sys.stderr
    else:
        fo = open(fn_out, 'w')
    for line in f:
        if not line:
            continue
        line = line.rstrip()
        fields = re.split(r'[\s\t]+', line)
        if fields[0] == 'chain':
            print(line, file=fo)
            source = fields[2]
            s_start = int(fields[5])
            s_end = int(fields[6])
            dest = fields[7]
            dest_len = int(fields[8])
            strand = fields[9]
            if strand == '+':
                d_start = int(fields[10])
                d_end = int(fields[11])
            else:
                d_start = dest_len - int(fields[10])
                d_end = dest_len - int(fields[11])
        elif len(fields) == 3:
            l = int(fields[0])
            ds = int(fields[1])
            dd = int(fields[2])
            if strand == '+':
                msg = f'\t{source}:{s_start}-{s_start+l}=>{dest}:{d_start}-{d_start+l} ({d_start-s_start})'
                print(line + msg, file=fo)
                s_start += (l + ds)
                d_start += (l + dd)
            else:
                msg = f'\t{source}:{s_start}-{s_start+l}=>{dest}:{d_start}-{d_start-l} ({d_start-s_start})'
                print(line + msg, file=fo)
                s_start += (l + ds)
                d_start -= (l + dd)
        elif len(fields) == 1 and fields[0] != '':
            l = int(fields[0])
            if strand == '+':
                msg = f'\t\t\t{source}:{s_start}-{s_start+l}=>{dest}:{d_start}-{d_start+l}'
            else:
                msg = f'\t\t\t{source}:{s_start}-{s_start+l}=>{dest}:{d_start}-{d_start-l}'
            print(line + msg + '\n', file=fo)



if __name__ == '__main__':
    args = parse_args()
    print(args.chain)
    print(args.out)
    verbosify_chain(fn_chain=args.chain, fn_out=args.out)
