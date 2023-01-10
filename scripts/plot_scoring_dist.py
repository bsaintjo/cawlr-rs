#!/usr/bin/env python3

import argparse
import matplotlib.pyplot as plt
import pickle
from dataclasses import dataclass
import itertools


@dataclass
class ScoreDist:
    filename: str
    dist: list[float]
    color: str


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-i", "--input", nargs="+", help="Score files to plot distribution"
    )
    parser.add_argument("-o", "--output", required=True, help="Output filename")

    args = parser.parse_args()

    all_sdists = []
    colors = itertools.cycle(["red", "blue", "green", "orange", "yellow"])
    for ms_filepath, color in zip(args.input, colors):
        with open(ms_filepath, "rb") as ms_file:
            data = pickle.load(ms_file)
            dist = data["bins"]
            sdist = ScoreDist(ms_filepath, dist, color)
            all_sdists.append(sdist)
    for sdist in all_sdists:
        n_bins = len(sdist.dist)
        dist_sum = sum(sdist.dist)
        bar_heights = [x / dist_sum for x in sdist.dist]
        plt.bar(
            x=[a / n_bins for a in range(0, n_bins)],
            height=bar_heights,
            width=1 / n_bins,
            color=sdist.color,
            alpha=0.3,
            label=sdist.filename,
        )
    plt.legend()
    plt.xlabel("Score")
    plt.ylabel("Density")
    plt.savefig(args.output)


if __name__ == "__main__":
    main()
