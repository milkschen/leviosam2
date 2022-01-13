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
    parser.add_argument(
        '--small_gold', action='store_true',
        help='Activate to save gold reads in dicts. [off]'
    )
    args = parser.parse_args()
    return args


def compare_fastq_large_gold(
    gold_fn1: str, gold_fn2: str,
    input_fn1: str, input_fn2: str,
    k: int, max_cnt: int
) -> None:
    print(f'Reading r1...', file=sys.stderr)
    r1 = {}
    with pysam.FastxFile(input_fn1) as f:
        for r in f:
            r1[r.name] = r.sequence[:k]
    print(f'R1: {len(r1)}')
    print(f'Reading r2...', file=sys.stderr)
    r2 = {}
    with pysam.FastxFile(input_fn2) as f:
        for r in f:
            r2[r.name] = r.sequence[:k]
    print(f'R2: {len(r2)}')
    if len(r1) != len(r2):
        print(f'Warning: lengths of R1 and R2 differ')

    err_dict = {}
    print(f'Reading g1...', file=sys.stderr)
    with pysam.FastxFile(gold_fn1) as f:
        for r in f:
            name = r.name.split()[0]
            if name in r1:
                # if seq not equal
                if r.sequence[:k] != r1[name][:k]:
                    print(name)
                    print(r)
                    print('  ' + r.sequence[:k])
                    print('  ' + r1[name][:k])
                    print()
                    err_dict[name] = [1, 0]
    print(len(err_dict))
    print(f'Reading g2...', file=sys.stderr)
    with pysam.FastxFile(gold_fn2) as f:
        for r in f:
            if r.name in r2:
                name = r.name.split()[0]
                # if seq not equal
                if r.sequence[:k] != r2[name][:k]:
                    if err_dict.get(name):
                        err_dict[name] = [err_dict[name][0], 1]
                    else:
                        err_dict[name] = [0, 1]

    err = {'only_r1': 0, 'only_r2': 0, 'both': 0, 'correct': 0}
    printed_err = {'only_r1': 0, 'only_r2': 0, 'both': 0, 'correct': 0}
    for i, (n, result) in enumerate(err_dict.items()):
        if result == [1, 1]:
            err['both'] += 1
            if printed_err['both'] < max_cnt:
                printed_err['both'] += 1
                print(f'[both] {n}')
        elif result == [0, 1]:
            err['only_r2'] += 1
            if printed_err['only_r2'] < max_cnt:
                printed_err['only_r2'] += 1
                print(f'[only_r2] {n}')
        elif result == [1, 0]:
            err['only_r1'] += 1
            if printed_err['only_r1'] < max_cnt:
                printed_err['only_r1'] += 1
                print(f'[only_r1] {n}')
    err['correct'] = len(r1) - err['only_r1'] - err['only_r2'] - err['both']
    print(err)


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

    if args.small_gold:
        compare_fastq(
            gold_fn1=args.gold_r1,
            gold_fn2=args.gold_r2,
            input_fn1=args.input_r1,
            input_fn2=args.input_r2,
            k=args.sampled_len,
            max_cnt=args.num_printed_err
        )
    else:
        compare_fastq_large_gold(
            gold_fn1=args.gold_r1,
            gold_fn2=args.gold_r2,
            input_fn1=args.input_r1,
            input_fn2=args.input_r2,
            k=args.sampled_len,
            max_cnt=args.num_printed_err
        )

