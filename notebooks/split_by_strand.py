#!/usr/bin/env python3

import argparse
from pathlib import Path


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", required=True, help="Input bed file")

    args = parser.parse_args()
    path = Path(args.input)
    minus_filepath = path.parent / (path.stem + ".minus.bed")
    plus_filepath = path.parent / (path.stem + ".plus.bed")
    nostrand_filepath = path.parent / (path.stem + ".none.bed")
    with open(args.input, "r") as bed_file, open(
        minus_filepath, "w"
    ) as minus_file, open(plus_filepath, "w") as plus_file, open(
        nostrand_filepath, "w"
    ) as none_file:
        header = next(bed_file)
        print(header, file=minus_file, end="")
        print(header, file=plus_file, end="")
        for line in bed_file:
            strand = line.rstrip().split("\t")[5]
            if strand == "+":
                output = plus_file
            elif strand == "-":
                output = minus_file
            else:
                output = none_file
            print(line, file=output, end="")


if __name__ == "__main__":
    main()
