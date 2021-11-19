'''
Summarize an alignment tag for a BAM file and show quantiles

Nae-Chyun Chen
Johns Hopkins University
2021
'''
import argparse
import numpy as np
import pysam
import sys

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-f', '--input', required=True,
        help='Path to the input SAM/BAM file. [required]'
    )
    parser.add_argument(
        '-t', '--tag', default=None, type=str,
        help='Alignment tag to summarize, e.g. AS. [None]'
    )
    parser.add_argument(
        '-q', '--mapq', action='store_true',
        help='Set to summarize MAPQ. [False]'
    )
    parser.add_argument(
        '-T', '--isize', action='store_true',
        help='Set to summarize template length (absolute value). [False]'
    )
    parser.add_argument(
        '-c', '--clipped_fraction', action='store_true',
        help='Set to summarzie clipped fraction. [False]'
    )
    args = parser.parse_args()
    return args


def summarize_aln_features(
    fn_input: str, tag: str=None,
    mapq: bool=False, isize: int=None, clipped_fraction: float=None
) -> None:
    results = []
    failed_cnt = 0
    other_cnt = 0
    total = 0
    with pysam.AlignmentFile(fn_input) as f:
        for read in f:
            if read.is_secondary or read.is_supplementary or read.is_unmapped:
                continue
            total += 1
            try:
                if tag:
                    results.append(int(read.get_tag(tag)))
                elif mapq:
                    results.append(read.mapping_quality)
                elif isize:
                    # If the alignment not considered paired properly, use isize=0
                    if read.flag & 2:
                        results.append(abs(read.template_length))
                    else:
                        other_cnt += 1
                elif clipped_fraction:
                    results.append(read.query_alignment_length / read.query_length)
            except:
                failed_cnt += 1

    print(f'Summarized {total} primary and aligned reads')
    if failed_cnt > 0:
        print(f'{failed_cnt} reads failed')

    results = np.array(sorted(results))
    l = len(results)
    m = np.mean(results)
    v = np.std(results)

    if tag:
        print(f'{tag} distribution:')
    elif mapq:
        print('MAPQ distribution:')
    elif isize:
        print('Template length distribution:')
    elif clipped_fraction:
        print('Clipped fraction distribution:')
    print(f'  mean = {m:.2f}')
    print(f'   std = {v:.2f}')
    print(f'-3 std = {m - 3 * v:.2f}')
    print(f'-2 std = {m - 2 * v:.2f}')
    print(f'-1 std = {m -     v:.2f}')
    print(f'+1 std = {m +     v:.2f}')
    print(f'+2 std = {m + 2 * v:.2f}')
    print(f'+3 std = {m + 3 * v:.2f}')
    print()

    print(f'   min = {results[0]}')
    print(f'    1% = {results[: int(l * 0.01)][-1]}')
    print(f'    3% = {results[: int(l * 0.03)][-1]}')
    print(f'    5% = {results[: int(l * 0.05)][-1]}')
    print(f'   10% = {results[: int(l * 0.10)][-1]}')
    print(f'   20% = {results[: int(l * 0.20)][-1]}')
    print(f'   50% = {results[: int(l * 0.50)][-1]}')
    print(f'   80% = {results[: int(l * 0.80)][-1]}')
    print(f'   90% = {results[: int(l * 0.90)][-1]}')
    print(f'   95% = {results[: int(l * 0.95)][-1]}')
    print(f'   97% = {results[: int(l * 0.97)][-1]}')
    print(f'   99% = {results[: int(l * 0.99)][-1]}')
    print(f'   max = {results[-1]}')
    if isize:
        print(f'not properly paird = {other_cnt} ({int(other_cnt / total * 100)}%)')


if __name__ == '__main__':
    args = parse_args()
    fn_input = args.input
    tag = args.tag
    mapq = args.mapq
    isize = args.isize
    clipped_fraction = args.clipped_fraction

    check_options = [tag != None, mapq != False, isize != False, clipped_fraction != False]
    if sum(check_options) != 1:
        print('Error: Only one out of the tag/mapq/isize/clipped_fraction options can be set')
        print('[tag, mapq, isize, clipped_fraction]')
        print(check_options)
        exit(1)

    summarize_aln_features(fn_input=fn_input, tag=tag, mapq=mapq,
                      isize=isize, clipped_fraction=clipped_fraction)
