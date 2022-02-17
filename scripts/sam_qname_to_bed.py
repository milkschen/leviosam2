'''
Convert SAM QNAME to BED records

This only works for QNAMEs like `chr1:1-100, which would be converted
to `chr1\t0\t100`

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
        '-s', '--sam', required=True,
        help='Path to the input SAM file. [required]'
    )
    parser.add_argument(
        '-f', '--flag_include', default=None, type=int, action='append',
        help='Only process records with specified flags. [None]'
    )
    parser.add_argument(
        '-F', '--flag_exclude', default=None, type=int, action='append',
        help='Exclude records with specified flags. [None]'
    )
    parser.add_argument(
        '-m', '--exclude_aln_num_range', default=None, type=str,
        help='Only print records outside a range of alignments, e.g. setting `-m 1-10` prints qnames aligned zero times or >10 times. [None]'
    )
    parser.add_argument(
        '-o', '--out', default='',
        help='Path to the ouput BED file. ['': print to sys.stdout]'
    )
    args = parser.parse_args()
    return args


def print_record(qname, fo) -> None:
    qname = qname.split(':')
    contig = qname[0]
    start = int(qname[1].split('-')[0]) - 1
    end = qname[1].split('-')[1]
    print(f'{contig}\t{start}\t{end}', file=fo)


def sam_qname_to_bed(args) -> None:
    f = pysam.AlignmentFile(args.sam)
    if args.out == '':
        fo = sys.stdout
    else:
        fo = open(args.out, 'w')
    
    aln_dict = {}
    for r in f:
        qname = r.query_name
        if qname not in aln_dict:
            aln_dict[qname] = 0
        if args.flag_include != None:
            is_valid = True
            for f_inc in args.flag_include:
                if not r.flag & f_inc:
                    is_valid = False
                    break
            if not is_valid:
                continue
        if args.flag_exclude != None:
            is_valid = True
            for f_exc in args.flag_exclude:
                if r.flag & f_exc:
                    is_valid = False
                    break
            if not is_valid:
                continue
        # qname = r.query_name.split(':')[0]
        aln_dict[qname] += 1

    for _, (qname, times) in enumerate(aln_dict.items()):
        if args.exclude_aln_num_range is None:
            print_record(qname, fo)
        else:
            min_exc_aln_num = int(args.exclude_aln_num_range.split('-')[0])
            max_exc_aln_num = int(args.exclude_aln_num_range.split('-')[1])
            if times < min_exc_aln_num or times > max_exc_aln_num:
                print_record(qname, fo)


if __name__ == '__main__':
    args = parse_args()
    sam_qname_to_bed(args)

