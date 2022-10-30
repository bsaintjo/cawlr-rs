#!/usr/bin/env python3

import argparse
import itertools
import pyBigWig as pbw
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from sklearn.preprocessing import StandardScaler

import cluster_region

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-i", "--input", required=True, help="Input bed file, usually from agg-block"
    )
    parser.add_argument(
        "-c",
        "--comparison",
        nargs="?",
        required=True,
        help="Data to compare to, in bigwig format",
    )
    parser.add_argument("-o", "--output", required=True, help="Path to output graphic")
    parser.add_argument(
        "--chrom", required=True, help="Chromosome of data to compare with"
    )
    parser.add_argument(
        "--start",
        type=int,
        required=True,
        help="Start genomic coordinate of data to compare with",
    )
    parser.add_argument(
        "--end",
        type=int,
        required=True,
        help="End genomic coordinate of data to compare with",
    )
    parser.add_argument(
        "--highlight", nargs="*", help="Region to highlight", default=[]
    )
    parser.add_argument("--suptitle", default="Aggregate plot", help="Plot title")
    args = parser.parse_args()

    columns = ["pos", "score", "group"]
    bw_dfs = []
    for bw_filepath in args.comparison:
        bw = pbw.open(bw_filepath)
        xs = bw.values(args.chrom, args.start, args.end)
        df = pd.DataFrame(zip(args.start, args.end + 1), xs, itertools.repeat("ortho"))
        bw_dfs.append(df)
    bw_df = pd.concat(bw_dfs, ignore_index=True, sort=False)
    bw_df.columns = columns

    agg_blocks = pd.read_csv(args.input, delimiter="\t", header=None)
    agg_blocks = agg_blocks.sort_values(by=1)
    agg_blocks = agg_blocks[(args.start <= agg_blocks[1]) & (agg_blocks[1] <= args.end)]
    agg_blocks[5] = agg_blocks[2].rolling(75, center=True).sum()
    agg_blocks["rolling_total"] = agg_blocks[3].rolling(75, center=True).sum()
    agg_blocks["rolling_avg"] = agg_blocks[5] / agg_blocks["rolling_total"]
    agg_blocks[7] = "sample"
    only_needed = agg_blocks[[1, 6, 7]]
    only_needed.columns = columns
    combined = pd.concat([only_needed, agg_blocks], ignore_index=True)

    scaler = StandardScaler()
    agg_blocks[6] = scaler.fit_transform(agg_blocks["rolling_avg"])

    _, ax = plt.subplots(1, 1, figsize=(10, 4))
    sns.lineplot(data=combined, x="pos", y="score", hue="group")
    for hl in args.highlights:
        (hstart, hend, hstrand) = cluster_region.parse_highlights(hl, args.start, args.end)
        fp_color, tp_color = cluster_region.strand_to_color(hstrand)
        ax.axvspan(hstart, hend, alpha=0.5, color="grey")
        ax.axvline(
            hstart,
            color=fp_color,
            alpha=0.5,
            linewidth=2,
            linestyle=":",
        )
        ax.axvline(hend, color=tp_color, alpha=0.5, linewidth=2, linestyle=":")
    plt.title(args.suptitle)
    plt.xlabel("Genomic position")
    plt.ylabel("Normalized Score")
    plt.savefig(args.input + ".agg.png", dpi=200)


if __name__ == "__main__":
    main()
