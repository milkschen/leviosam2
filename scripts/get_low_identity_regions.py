"""
Extract low identity regions from a chain summary

Nae-Chyun Chen
Johns Hopkins University
2021
"""
import argparse
import sys

import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-s",
        "--summary",
        required=True,
        help="Path to the input chain summary file. [required]",
    )
    parser.add_argument(
        "-o",
        "--out",
        default="",
        help="Path to the ouput BED file. [" ": print to sys.stdout]",
    )
    parser.add_argument(
        "-l",
        "--identity_cutoff",
        required=True,
        type=float,
        help="Identity cutoff. Values lower than this will be excluded. [required]",
    )
    parser.add_argument(
        "--source",
        action="store_true",
        help="Set to write wrt the source reference; otherwise write wrt the dest reference.",
    )
    args = parser.parse_args()
    return args


def get_low_identity_regions(
    summary: str, out: str, identity_cutoff: float, source: bool
) -> None:
    if source:
        print("Output with respect to the source reference", file=sys.stderr)
    else:
        print("Output with respect to the dest reference", file=sys.stderr)

    df = pd.read_csv(summary, sep="\t")
    fo = open(out, "w")
    for i in range(df.shape[0]):
        if df.iloc[i][1] < identity_cutoff:
            rec = df.iloc[i]
            if source:
                print(f"{rec[2]}\t{rec[3]}\t{rec[4]}", file=fo)
            else:
                print(f"{rec[6]}\t{rec[7]}\t{rec[8]}", file=fo)
    fo.close()


if __name__ == "__main__":
    args = parse_args()
    get_low_identity_regions(
        summary=args.summary,
        out=args.out,
        identity_cutoff=args.identity_cutoff,
        source=args.source,
    )
