'''
Extract the mate of reads in a single-end FASTQ from a BAM file

This script is designed for the pipeline where
    (a) a BAM file is split by an indicator
    (b) reads in a sub BAM file are converted to FASTQ,
        where some are properly paired but some aren't
    (c) we use this script to find the mate for those that aren't paired

Nae-Chyun Chen
Johns Hopkins University
2021
'''
import argparse
import os
import pysam
import sys

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-i', '--input', required=True,
        help='Path to a SAM/BAM. [required]'
    )
    parser.add_argument(
        '-s', '--reads', required=True,
        help='Path to a single-end FASTQ. [required]'
    )
    parser.add_argument(
        '-op', '--out_prefix', required=True,
        help='Prefix of outputs (a properly-paired BAM and a pair of FASTQs. [required]'
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


def extract_unpaired_reads(
    fn_reads: str, out_prefix: str, fn_input: str):
    dict_reads = {}
    read_name = ''
    with open(fn_reads, 'r') as fr:
        for i, line in enumerate(fr):
            line = line.rstrip()
            if i % 4 == 0:
                read_name = line[1:]
                assert not dict_reads.get(read_name)
                dict_reads[read_name] = [None, None]
            elif i % 4 == 1:
                dict_reads[read_name][0] = line
            elif i % 4 == 2:
                continue
            elif i % 4 == 3:
                dict_reads[read_name][1] = line
                read_name = None
    print(f'Read {i/4} records from {fn_reads}', file=sys.stderr)

    fn_paired_output = out_prefix + '-paired.bam'
    fn_fq1 = out_prefix + '-singleton-R1.fq'
    fn_fq2 = out_prefix + '-singleton-R2.fq'
    for fn in [fn_paired_output, fn_fq1, fn_fq2]:
        if os.path.exists(fn):
            print((f'Error: file {fn} exists. '
                   f'Please remove it if you want to proceed'),
                   file=sys.stderr)
            exit(1)
    fa = pysam.AlignmentFile(fn_input, 'r')
    fo_paired = pysam.AlignmentFile(fn_paired_output, 'wb', template=fa)
    fo_fq1 = open(fn_fq1, 'w')
    fo_fq2 = open(fn_fq2, 'w')
    print(f'Read from {fn_input}')
    print(f'Write filtered BAM to {fn_paired_output}')
    print(f'Write extracted R1 to {fn_fq1}')
    print(f'Write extracted R2 to {fn_fq2}')

    for read in fa:
        if read.query_name in dict_reads:
            seq = read.get_forward_sequence()
            qual = ''.join(map(lambda x: chr( x+33 ), read.query_qualities))
            # See https://github.com/pysam-developers/pysam/issues/839#issuecomment-560130300
            # if read.is_reverse:
            #     qual = read.query_qualities[::-1]
            #     seq = reverse_complement(read.query_sequence)
            # else:
            #     qual = read.query_qualities
            #     seq = read.query_sequence
            singleton = dict_reads[read.query_name]
            if read.is_read1 and not read.is_secondary and not read.is_supplementary:
                fo_fq1.write(f'@{read.query_name}\n{seq}\n+\n{qual}\n')
                fo_fq2.write(f'@{read.query_name}\n{singleton[0]}\n+\n{singleton[1]}\n')
            elif read.is_read2 and not read.is_secondary and not read.is_supplementary:
                fo_fq2.write(f'@{read.query_name}\n{seq}\n+\n{qual}\n')
                fo_fq1.write(f'@{read.query_name}\n{singleton[0]}\n+\n{singleton[1]}\n')
        else:
            fo_paired.write(read)


if __name__ == '__main__':
    args = parse_args()
    extract_unpaired_reads(
        fn_reads=args.reads, out_prefix=args.out_prefix, fn_input=args.input)
