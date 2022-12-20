#!/usr/bin/env python3

import argparse
from typing import List
import pyarrow.feather as feather
import matplotlib.pyplot as plt

def table_to_gmm_scores(table) -> List[float]:
    scores = []
    for chunk in table:
        for read in chunk:
            for score in read["scores"]:
                signal_score = score["signal_score"].as_py()
                kmer = score["kmer"]
                if signal_score is not None:
                    scores.append(signal_score)
    return scores


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", nargs="+", help="Score files to plot distribution")
    parser.add_argument("-o", "--output", required=True, help="Output filename")

    args = parser.parse_args()

    score_file = args.input[0]
    table = feather.read_table(score_file)
    print(table)


if __name__ == "__main__":
    main()
