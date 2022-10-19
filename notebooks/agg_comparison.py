#!/usr/bin/env python3

import argparse
import pyBigWig as pbw
import matplotlib.pyplot as plt


def parse_agg_block(input, chrom, start, end):
    acc = []
    for line in input:
        arr = [0] * (end - start)
        line = line.rstrip().split('\t')
        position = int(line[1])
        if (chrom == line[0]) and (start <= position <= end):
            frac = float(line[-1])
            rel_position = position - start
            acc[rel_position] = frac
            acc.append(arr)
    return acc


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", required=True, help="Input bed file, usually from agg-block")
    parser.add_argument("-c", "--comparison", required=True, help="Data to compare to, in bigwig format")
    parser.add_argument("-o", "--output", required=True, help="Path to output graphic")
    parser.add_argument("--chrom", required=True, help="Chromosome of data to compare with")
    parser.add_argument("--start", type=int, required=True, help="Start genomic coordinate of data to compare with")
    parser.add_argument("--end", type=int, required=True, help="End genomic coordinate of data to compare with")
    parser.add_argument("--suptitle", default="Aggregate plot", help="Plot title")
    args = parser.parse_args()

    comparison = pbw.open(args.comparison)
    values = comparison.values(args.chrom, args.start, args.end)

    with open(args.input, 'r') as bed_file:
        agg_blocks = parse_agg_block(args.input, args.chrom, args.start, args.end)

    fig, ax = plt.subplots(1, 1, figsize=(6, 14))
    ax.plot(values)
    fig.suptitle(args.suptitle)
    fig.supxlabel("Genomic position")
    fig.supylabel("Score")

if __name__ == "__main__":
    main()