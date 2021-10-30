'''
Add info for a chain file to facilitate debugging

Nae-Chyun Chen
Johns Hopkins University
2021
'''
import argparse
import leviosam_utils
import pysam
import re
import sys

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
    parser.add_argument(
        '-b', '--bed_prefix', default='',
        help='Prefix to the BED files that include all the chain segments. [empty string]'
    )
    parser.add_argument(
        '-s', '--summary', default='',
        help='Path to the output summary. Leave empty for no output. [empty string]'
    )
    parser.add_argument(
        '-f1', '--ref1', default='',
        help='Path to the source reference (optional).'
    )
    parser.add_argument(
        '-f2', '--ref2', default='',
        help='Path to the destination reference (optional).'
    )
    args = parser.parse_args()
    return args


def reverse_complement(seq):
    d = {'A': 'T', 'a': 'T', 'C': 'G', 'c': 'G',
         'G': 'C', 'g': 'G', 'T': 'A', 't': 'A',
         'N': 'N'}
    rc = ''
    for s in seq:
        if s in d:
            rc += d[s]
        else:
            print(f'Base "{s}" is not a known nucleotide and is converted to N',
                  file=sys.stderr)
            rc += 'N'
    return rc[::-1]


def compute_hamming_dist(
    forward,
    ref1, contig1, start1, end1,
    ref2, contig2, start2, end2
):
    if contig1 in ref1:
        s1 = ref1[contig1][start1: end1]
    else:
        print(f'Warning : {contig1} not in ref1', file=sys.stderr)
        return 0
    if contig2 in ref2:
        s2 = ref2[contig2][start2: end2]
        if not forward:
            s2 = reverse_complement(s2)
    else:
        print(f'Warning : {contig2} not in ref2', file=sys.stderr)
        return 0

    try:
        assert (len(s1) == len(s2))
    except:
        print('Error: lengths do not match', file=sys.stderr)
        print(len(s1), len(s2), file=sys.stderr)
        print(f'{contig1}:{start1}-{end1}', file=sys.stderr)
        print(f'{contig2}:{start2}-{end2}', file=sys.stderr)
        exit(1)

    if len(s1) == 0:
        return 0
    idy = 0
    for i_s, s in enumerate(s1):
        if s == s2[i_s]:
            idy += 1
    idy /= len(s1)
    return idy


''' Write a chain record in the BED format '''
def write_to_summary(fs_fn, fs, strand, l, hd, source, s_start, dest, d_start):
    if fs_fn:
        if strand == '+':
            print((f'{l}\t{hd:.6f}\t{source}\t{s_start}\t{s_start+l}\t+'
                   f'\t{dest}\t{d_start}\t{d_start+l}'), file=fs)
        else:
            print((f'{l}\t{hd:.6f}\t{source}\t{s_start}\t{s_start+l}\t+'
                   f'\t{dest}\t{d_start-l}\t{d_start}'), file=fs)


def verbosify_chain(args):
    f = open(args.chain, 'r')
    if args.out == '':
        fo = sys.stderr
    else:
        fo = open(args.out, 'w')
    if args.bed_prefix != '':
        f_sbed = open(args.bed_prefix + '-source.bed', 'w')
        f_dbed = open(args.bed_prefix + '-dest.bed', 'w')

    ref1 = leviosam_utils.read_fasta(args.ref1)
    ref2 = leviosam_utils.read_fasta(args.ref2)

    check_hdist = True if (ref1 != {} and ref2 != {}) else False
    if not check_hdist:
        hd = None
    if args.summary:
        assert check_hdist == True
        fs = open(args.summary, 'w')
        # Write header
        print(f'SIZE\tHDIST\tSOURCE\tS_START\tS_END\tSTRAND\tDEST\tD_START\tD_END', file=fs)
    else:
        fs = None

    total_bases = 0
    total_bases_idy = 0
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
            total_bases += l
            ds = int(fields[1])
            dd = int(fields[2])
            if strand == '+':
                msg = f'\t{source}:{s_start}-{s_start+l}=>{dest}:{d_start}-{d_start+l} ({d_start-s_start})'
                if args.bed_prefix != '':
                    f_sbed.write(f'{source}\t{s_start}\t{s_start+l}\n')
                    f_dbed.write(f'{dest}\t{d_start}\t{d_start+l}\n')
                if check_hdist:
                    hd = compute_hamming_dist(
                        True,
                        ref1, source, s_start, s_start+l,
                        ref2, dest, d_start, d_start+l)
                    msg += f'\t{hd}'
            else:
                msg = f'\t{source}:{s_start}-{s_start+l}=>{dest}:{d_start}-{d_start-l} ({d_start-s_start})'
                if args.bed_prefix != '':
                    f_sbed.write(f'{source}\t{s_start}\t{s_start+l}\n')
                    f_dbed.write(f'{dest}\t{d_start-l}\t{d_start}\n')
                if check_hdist:
                    hd = compute_hamming_dist(
                        False,
                        ref1, source, s_start, s_start+l,
                        ref2, dest, d_start-l, d_start)
                    msg += f'\t{hd}'
            print(line + msg, file=fo)
            write_to_summary(args.summary, fs, strand, l, hd, source, s_start, dest, d_start)
            if strand == '+':
                s_start += (l + ds)
                d_start += (l + dd)
            else:
                s_start += (l + ds)
                d_start -= (l + dd)
            # write_to_summary(args.summary, fs, strand, l, hd, source, s_start, dest, d_start)
            if check_hdist:
                total_bases_idy += (l * hd)
        elif len(fields) == 1 and fields[0] != '':
            l = int(fields[0])
            total_bases += l
            if strand == '+':
                msg = f'\t\t\t{source}:{s_start}-{s_start+l}=>{dest}:{d_start}-{d_start+l}'
                if args.bed_prefix != '':
                    f_sbed.write(f'{source}\t{s_start}\t{s_start+l}\n')
                    f_dbed.write(f'{dest}\t{d_start}\t{d_start+l}\n')
                if check_hdist:
                    hd = compute_hamming_dist(
                        True,
                        ref1, source, s_start, s_start+l,
                        ref2, dest, d_start, d_start+l)
                    msg += f'\t{hd}'
            else:
                msg = f'\t\t\t{source}:{s_start}-{s_start+l}=>{dest}:{d_start}-{d_start-l}'
                if args.bed_prefix != '':
                    f_sbed.write(f'{source}\t{s_start}\t{s_start+l}\n')
                    f_dbed.write(f'{dest}\t{d_start-l}\t{d_start}\n')
                if check_hdist:
                    hd = compute_hamming_dist(
                        False,
                        ref1, source, s_start, s_start+l,
                        ref2, dest, d_start-l, d_start)
                    msg += f'\t{hd}'
            write_to_summary(args.summary, fs, strand, l, hd, source, s_start, dest, d_start)
            if check_hdist:
                total_bases_idy += (l * hd)
            print(line + msg + '\n', file=fo)

    print(f'Total number of gapless aligned bases = {total_bases}', file=sys.stderr)
    print(f'Total number of gapless matched bases = {total_bases_idy}', file=sys.stderr)


if __name__ == '__main__':
    args = parse_args()

    print('Input chain:', args.chain, file=sys.stderr)
    print('Output chain (verbose):', args.out, file=sys.stderr)

    assert (args.ref1 != '' and args.ref2 != '') or (args.ref1 == '' and args.ref2 == '')
    print('Ref1:', args.ref1, file=sys.stderr)
    print('Ref2:', args.ref2, file=sys.stderr)

    verbosify_chain(args)
