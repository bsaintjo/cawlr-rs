#!/usr/bin/env python3
"""
This script provides an example of reading a cawlr score file using pyarrow.

In this example, we load an Arrow file passed via a command line. We change
the name of the first read via its metadata, and save it to a new file. We
then create a new file from scratch and then overwrite our previous example.
"""

import argparse
from pyarrow import feather, Table


def from_scratch():
    """This function gives an example of constructing a dictionary for
    creating a scored Arrow file from scratch."""

    # metadata gives information about the read, chromosome, position, length,
    # etc.
    metadata = {
        "chrom": "chrI",
        "length": 100,
        "name": "example",
        # seq will always be an empty string
        "seq": "",
        "start": 100,
        # True is positive strand, False is negative strand
        "strand": True,
    }

    # scores is a list of dictionaries, each item represents a position within
    # the read and the score it has.
    scores = [
        {
            # left most position represent the position
            # ie T is 101, A is 102, T is 103, T is 104, etc.
            "kmer": "TATTGA",
            "pos": 101,
            # For your case, "score" and "signal_score" should be the same
            # This is where you will put in the scores output by the method
            # you are working with.
            "score": 0.123,
            "signal_score": 0.123,
            # For your case, skip_score is always 0.0
            "skip_score": 0.0,
            # For your case, skipped will always be set to False
            "skipped": False,
        },
        # Second position
        {
            # If space/time are becoming more of an issue, only write scores
            # for kmers that contain modification motif.
            # ie for GpC methyltransferase, only score if kmer contains a GC
            # ie for Addseq, only score if kmer contains an AT or TA
            "kmer": "ATTGAC",
            "pos": 102,
            "score": 0.456,
            "signal_score": 0.456,
            "skip_score": 0.0,
            "skipped": False,
        },
    ]

    # Each sequencing read is a dictionary with the metadata dictionary and
    # list of score dictionaries.
    cawlr_read = {"metadata": metadata, "scores": scores}

    scored = {"scored": [cawlr_read]}
    return scored


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", required=True, help="Scored Arrow file")
    parser.add_argument(
        "-o", "--output", required=True, help="Output Scored Arrow file name"
    )
    args = parser.parse_args()

    table = feather.read_table(args.input)
    # fdict is a dictionary with one key "scored",
    # whose value is a list of dictionaries
    fdict = table.to_pydict()

    fdict["scored"][0]["metadata"]["name"] = "new_read_name"
    # ......^ Access the list of dictionaries
    # ..............^ Each item in the list is a read, lets just take the
    #                 first one in this example.
    # ....................^ Each read is a dictionary with two keys,
    #                       "metadata" and "scores"
    # ...............................^ This returns a dictionary with several
    #                                  several pieces of info, one of which is
    #                                  the read name.

    output = Table.from_pydict(fdict)
    feather.write_feather(output, args.output)

    # Another example, creating it from scratch and overwriting the last output
    # Comment out this portion if you want to do something with the modified
    # version
    new_fdict = from_scratch()
    # Note: Currently this line will fail because pyarrow doesn't support
    # dense_union. Currently fixing to resolve this issue, for now, don't pass
    # table.schema
    # output = Table.from_pydict(new_fdict, schema=table.schema)
    output = Table.from_pydict(new_fdict)
    feather.write_feather(output, args.output)


if __name__ == "__main__":
    main()
