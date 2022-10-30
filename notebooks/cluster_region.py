#!/usr/bin/env python3
"""
This script takes single molecule analysis results in a bed format, usually
from cawlr sma, and performs k-means clustering to determine different
nucleosome configurations in a given region.
"""

import argparse
from pathlib import Path
import sys
from collections import defaultdict
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans
import numpy as np


def bed_line_cluster_array(line, cluster_start, cluster_stop):
    """
    Parse line from a bed file, converting it into a list which
    each position representing a basepair in the genomic position,
    1 indicates nucleosome, 0 indicates linker, None indicates
    the read ended before that position.
    """
    line = line.rstrip().split("\t")
    read_start = int(line[1])
    read_stop = int(line[2])
    blocks = [int(x) for x in line[10].split(",")]
    starts = [int(x) for x in line[11].split(",")]

    cluster_idxs = {pos: None for pos in range(cluster_start, cluster_stop + 1)}
    for pos in range(read_start, read_stop + 1):
        if pos in cluster_idxs:
            cluster_idxs[pos] = 0.0

    for st, bl in zip(starts, blocks):
        for x in range(read_start + st, read_start + st + bl + 1):
            if x in cluster_idxs:
                cluster_idxs[x] = 1
    cluster_array = [a[1] for a in sorted(cluster_idxs.items(), key=lambda x: x[0])]
    return cluster_array


def pct_full(arr):
    """Count what percent of the read overlaps with the read, None in
    the array means that there was no data for it."""
    n_none = sum(1 for x in arr if x is None)
    return 1 - float(n_none) / len(arr)


def convert_nones(arr):
    """Convert Nones to something else so KMeans is able to cluster on it"""
    return [x if x is not None else 0.5 for x in arr]


def split_clusters(cresults, carrays):
    """Build dictionary mapping KMeans cluster label to the read
    array containing information about positions of nucleosomes"""
    labels = defaultdict(list)
    for idx, res in enumerate(carrays):
        label = cresults[idx]
        labels[label].append(res)
    return labels


def parse_highlights(s: str, region_start, region_end):
    xs = s.split(":")
    strand = xs[-1]
    xs = xs[0].split("-")

    start = int(xs[0])
    start = start if start > region_start else region_start

    end = int(xs[1])
    end = end if end < region_end else region_end

    return (start, end, strand)


def strand_to_color(hstrand: str):
    if hstrand == "+":
        fp_color = "green"
        tp_color = "red"
    elif hstrand == "-":
        fp_color = "red"
        tp_color = "green"
    else:
        fp_color = "black"
        tp_color = "black"
    return fp_color, tp_color


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-i", "--input", required=True, help="Input bed file, usually from cawlr sma"
    )
    parser.add_argument(
        "-s", "--start", type=int, required=True, help="Start of region"
    )
    parser.add_argument("-e", "--end", type=int, required=True, help="End of region")
    parser.add_argument(
        "-p",
        "--pct",
        type=float,
        default=0.9,
        help="Perecent of the region that should be covered for a read to be valid",
    )
    parser.add_argument(
        "-n", "--n-clusters", type=int, default=3, help="Number of clusters"
    )
    parser.add_argument("--suptitle", default="Region name", help="Figure title")

    parser.add_argument(
        "--highlight",
        nargs="*",
        help="Highlight particular regions, usually gene bodies, etc., format is usually {start}-{end}:{strand}",
        default=[],
    )

    args = parser.parse_args()

    input_path = Path(args.input)
    output = input_path.parent / (input_path.stem + ".cluster.png")

    highlights = [parse_highlights(h, args.start, args.end) for h in args.highlight]

    acc = []
    with open(args.input, "r") as bedfile:
        next(bedfile)  # Skip header
        for line in bedfile:
            arr = bed_line_cluster_array(
                line, cluster_start=args.start, cluster_stop=args.end
            )
            if pct_full(arr) > args.pct:
                acc.append(convert_nones(arr))
    if not acc:
        print("Empty results, input maybe empty or lower % threshold")
        sys.exit(1)

    kmeans = KMeans(n_clusters=args.n_clusters)
    results = kmeans.fit_predict(acc)

    label_to_arr = split_clusters(results, acc)

    fig, axs = plt.subplots(
        nrows=args.n_clusters, ncols=1, sharex=True, figsize=(15, 6)
    )
    for idx, arrs in label_to_arr.items():
        axs[idx].imshow(arrs, aspect="auto", interpolation="none")
        for (hstart, hend, hstrand) in highlights:
            fp_color, tp_color = strand_to_color(hstrand)
            axs[idx].axvspan(
                hstart - args.start, hend - args.start, alpha=0.5, color="grey"
            )
            axs[idx].axvline(
                hstart - args.start,
                color=fp_color,
                alpha=0.5,
                linewidth=2,
                linestyle=":",
            )
            axs[idx].axvline(
                hend - args.start, color=tp_color, alpha=0.5, linewidth=2, linestyle=":"
            )
    fig.add_subplot(111, frameon=False)
    plt.tick_params(
        labelcolor="none",
        which="both",
        top=False,
        bottom=False,
        left=False,
        right=False,
    )

    axs[-1].set_xticks(np.linspace(0, (args.end - args.start), 20, dtype=int))
    axs[-1].set_xticklabels(
        np.linspace(args.start, args.end, 20, dtype=int), rotation=45
    )

    plt.subplots_adjust(top=0.95, bottom=0.15)
    plt.xlabel("Genomic coordinate", labelpad=30)
    plt.ylabel("Reads")
    plt.title(args.suptitle)

    plt.savefig(output, dpi=500)
    print("Output file to: ", output)


if __name__ == "__main__":
    main()
