'''
Compare two sets of FASTQ files

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
        '-g1', '--gold_r1', required=True,
        help='Path to the truth FASTQ-R1 file. [required]'
    )
    parser.add_argument(
        '-g2', '--gold_r2', required=True,
        help='Path to the truth FASTQ-R2 file. [required]'
    )
    parser.add_argument(
        '-r1', '--input_r1', required=True,
        help='Path to the truth FASTQ-R1 file. [required]'
    )
    parser.add_argument(
        '-r2', '--input_r2', required=True,
        help='Path to the truth FASTQ-R2 file. [required]'
    )
    parser.add_argument(
        '-k', '--sampled_len', default=30, type=int,
        help='Only store the first `k`-bp of a seq to save memory. [30]'
    )
    parser.add_argument(
        '--num_printed_err', default=3, type=int,
        help='Num of printed errors during processing. [3]'
    )
    args = parser.parse_args()
    return args


def compare_fastq(
    gold_fn1: str, gold_fn2: str,
    input_fn1: str, input_fn2: str,
    k: int, max_cnt: int
) -> None:
    gold_r1 = {}
    with pysam.FastxFile(gold_fn1) as f:
        for r in f:
            gold_r1[r.name] = [r.sequence[:k]]
    gold_r2 = {}
    with pysam.FastxFile(gold_fn2) as f:
        for r in f:
            gold_r2[r.name] = [r.sequence[:k]]

    r1 = {}
    with pysam.FastxFile(input_fn1) as f:
        for r in f:
            r1[r.name] = [r.sequence[:k]]
    r2 = {}
    with pysam.FastxFile(input_fn2) as f:
        for r in f:
            r2[r.name] = [r.sequence[:k]]

    err = {'only_r1': 0, 'only_r2': 0, 'both': 0, 'correct': 0}
    printed_err = {'only_r1': 0, 'only_r2': 0, 'both': 0, 'correct': 0}
    for i, (n, v1) in enumerate(r1.items()):
        r1_diff = (gold_r1[n] != v1)
        r2_diff = (gold_r2[n] != r2[n])
        if r1_diff and r2_diff:
            err['both'] += 1
        elif r1_diff:
            err['only_r1'] += 1
            if printed_err['only_r1'] < max_cnt:
                printed_err['only_r1'] += 1
                print(f'[only_r1] {n}')
        elif r2_diff:
            err['only_r2'] += 1
            if printed_err['only_r2'] < max_cnt:
                printed_err['only_r2'] += 1
                print(f'[only_r2] {n}')
        else:
            err['correct'] += 1

    print(err)


if __name__ == '__main__':
    args = parse_args()

    print(f'gold_fn1 = {args.gold_r1}')
    print(f'gold_fn2 = {args.gold_r2}')
    print(f'fn1 = {args.input_r1}')
    print(f'fn2 = {args.input_r2}')

    compare_fastq(
        gold_fn1=args.gold_r1,
        gold_fn2=args.gold_r2,
        input_fn1=args.input_r1,
        input_fn2=args.input_r2,
        k=args.sampled_len,
        max_cnt=args.num_printed_err
    )

