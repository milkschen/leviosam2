"""
Convert a FAI file to a BED file

Nae-Chyun Chen
Johns Hopkins University
2021
"""
import argparse
import sys


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-f",
        "--fai",
        required=True,
        help="Path to the input FAI file. [required]",
    )
    parser.add_argument(
        "-o",
        "--out",
        default="",
        help="Path to the ouput BED file. [" ": print to sys.stdout]",
    )
    args = parser.parse_args()
    return args


def fai_to_bed(args):
    fb = open(args.fai, "r")
    if args.out == "":
        fo = sys.stdout
    else:
        fo = open(args.out, "w")

    for line in fb:
        line = line.split()
        print(f"{line[0]}\t0\t{line[1]}", file=fo)


if __name__ == "__main__":
    args = parse_args()
    fai_to_bed(args)
