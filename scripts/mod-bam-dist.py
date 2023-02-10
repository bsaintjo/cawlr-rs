#!/usr/bin/env python3

import pysam
import argparse
import matplotlib.pyplot as plt
from dataclasses import dataclass
import random


@dataclass
class ScoreDist:
    name: str
    scores: list[float]

    def __init__(self, name: str, filepath: str, cutoff: float = 0.0):
        self.name = name
        self.scores = bam_to_scores(filepath, cutoff)

    def plot_hist(self, ax, color: str):
        ax.hist(
            self.scores,
            bins=25,
            range=(0.0, 1.0),
            histtype="step",
            density=True,
            alpha=0.8,
            color=color,
            label=self.name,
        )


def bam_to_scores(filepath: str, cutoff: float) -> list[float]:
    bam_file = pysam.AlignmentFile(filepath)
    all_scores = list()
    for alignment in bam_file:
        if alignment.has_tag("Ml"):
            for x in alignment.get_tag("Ml"):
                if random.random() > cutoff:
                    all_scores.append(x / 255.0)
    print(f"{filepath} with {len(all_scores)} scores")
    return all_scores


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-n",
        "--neg-ctrl",
        required=True,
        help="Input BAM file with MM/MD tags from negative control",
    )
    parser.add_argument(
        "-p",
        "--pos-ctrl",
        required=True,
        help="Input BAM file with MM/MD tags from positive control",
    )
    parser.add_argument(
        "-s",
        "--sample",
        help="Input BAM file with MM/MD tags from sample control",
    )
    args = parser.parse_args()
    neg_scores = ScoreDist("(-) ctrl", args.neg_ctrl, cutoff=0.95)
    pos_scores = ScoreDist("(+) ctrl", args.pos_ctrl)

    fig, ax = plt.subplots(nrows=1, ncols=1)
    neg_scores.plot_hist(ax, "purple")
    pos_scores.plot_hist(ax, "orange")
    if args.sample is not None:
        sample_scores = ScoreDist("sample data", args.sample)
        sample_scores.plot_hist(ax, "green")
    fig.legend()
    fig.savefig("test.png")


if __name__ == "__main__":
    main()
