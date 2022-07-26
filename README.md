# cawlr 𓅨

[![License](https://img.shields.io/badge/license-BSD_3--Clause-informational)](./LICENSE)

**C**hromatin **a**ccesibility **w**ith **l**ong **r**eads (`cawlr`) is a tool for detecting accessible regions of chromatin on single molecules. The tool works with nanopore data via [`nanopolish`](https://github.com/jts/nanopolish).

<!-- TODO or PacBio data via [`kineticsTools`](https://github.com/PacificBiosciences/kineticsTools). -->

`cawlr` is flexible and can work with data that uses various modifications for chromatin accessibility detection. Outputs of all `cawlr` modules are in Apache Arrow format can be manipulated via any programming language that has an [Apache Arrow](https://arrow.apache.org/install/) library. Furthermore, `cawlr` includes various scripts for plotting and evaluating the results.

## Table of Contents

- [cawlr 𓅨](#cawlr-𓅨)
  - [Table of Contents](#table-of-contents)
  - [Installation](#installation)
    - [Installing rust via rustup](#installing-rust-via-rustup)
    - [Installing system dependencies](#installing-system-dependencies)
      - [Installing on Ubuntu/Debian](#installing-on-ubuntudebian)
      - [Installing on CentOS 7](#installing-on-centos-7)
      - [Installing on MacOS](#installing-on-macos)
      - [Installing on Windows](#installing-on-windows)
    - [Installing cawlr](#installing-cawlr)
      - [Latest from git](#latest-from-git)
      - [(Optional) Run tests](#optional-run-tests)
    - [Python dependencies](#python-dependencies)
  - [Nanopore data preparation](#nanopore-data-preparation)
  - [Usage](#usage)
    - [`cawlr collapse`](#cawlr-collapse)
    - [`cawlr train`](#cawlr-train)
    - [`cawlr rank`](#cawlr-rank)
    - [`cawlr score`](#cawlr-score)
    - [`cawlr sma`](#cawlr-sma)
  - [Models](#models)
  - [Example `cawlr` vignette](#example-cawlr-vignette)
  - [QC Scripts](#qc-scripts)
    - [`plot_gmm_models.py`](#plot_gmm_modelspy)
    - [`plot_read_length_vs_mod_rate.py`](#plot_read_length_vs_mod_ratepy)
    - [`plot_scores_with_inference.py`](#plot_scores_with_inferencepy)
  - [Citations](#citations)

## Installation

### Installing rust via rustup

It is recommended to install the Rust compiler is through [rustup.rs](https://rustup.rs/).

### Installing system dependencies

Ensure you have these installed on your system before installing.

- make
- openblas
- perl
- gcc
- gfortran

#### Installing on Ubuntu/Debian

```bash
sudo apt update
sudo apt install cmake libopenblas-dev
```

#### Installing on CentOS 7

```bash
yum install epel-release
yum install perl make gcc gcc-gfortran openblas-devel
```

#### Installing on MacOS

TODO

#### Installing on Windows

TODO

### Installing cawlr

#### Latest from git

```bash
git clone https://github.com/bsaintjo/cawlr-rs.git
cd cawlr-rs
cargo install --path .
```

While the tool was developed to be as portable as possible, you can sacrifice portability for small speed increases by installing with compiling options shown below. The `fast` feature requires GCC version >= 4.9

```bash
RUSTFLAGS="-C target-cpu=native" cargo install --path . --features fast
```

#### (Optional) Run tests

```bash
cargo test
```

### Python dependencies

## Nanopore data preparation

In order to prepare data for `cawlr` you need to install:

- [`samtools`](http://www.htslib.org/)
- [`nanopolish`](https://github.com/jts/nanopolish)
- [`minimap2`](https://github.com/lh3/minimap2)

Example command of running nanopolish to set up inputs

```bash
$ nanopolish index -d path/to/fast5 reads.fasta
$ minimap2 -ax map-ont --sam-hit-only --secondary=no genome.fasta reads.fasta | \
    samtools sort -o aln.sorted.bam -T reads.tmp
$ samtools index aln.sorted.bam
$ nanopolish eventalign --reads reads.fasta \
    --bam aln.sorted.bam \
    --genome genome.fasta \
    --scale-events --samples \
    --print-read-names >eventalign.txt
```

<!-- TODO Confirm that cawlr collapse will fail without `--sam-hit-only` -->

For training the models, its good practice to avoid using multiple alignments for the same read.
To avoid this, you can filter the bam files with the command below before running `nanopolish`. The command filters with the `-f 20` to filter out reads with a MAPQ below 20, and `-F 0x900` removes supplementary and secondary reads as well.

```bash
# After minimap2 command
$ samtools view -b -q 20 -F 0x900 aln.sorted.bam >filtered.aln.sorted.bam
```

Once completed you can confirm that alignments have been filtered correctly with:

```bash
samtools flagstats filtered.aln.sorted.bam
```

## Usage

### `cawlr collapse`

```bash
cawlr-collapse 

USAGE:
    cawlr collapse [OPTIONS] --input <INPUT> --bam <BAM> --output <OUTPUT>

OPTIONS:
    -b, --bam <BAM>              Path to BAM alignment file used in nanopolish eventalign
    -c, --capacity <CAPACITY>    Number of eventalign records to hold in memory [default: 2048]
    -h, --help                   Print help information
    -i, --input <INPUT>          Path to nanopolish eventalign output with samples column
    -o, --output <OUTPUT>        Path to output file in Apache Arrow format, defaults to stdout if
                                 no argument provided
    -q, --quiet                  Less output per occurrence
    -v, --verbose                More output per occurrence
```

### `cawlr train`

```bash
$ cawlr help train
cawlr-train 
For each kmer, train a two-component gaussian mixture model and save models to a file

USAGE:
    cawlr train [OPTIONS] --input <INPUT> --output <OUTPUT> --genome <GENOME>

OPTIONS:
    -g, --genome <GENOME>      Path to genome fasta file
    -h, --help                 Print help information
    -i, --input <INPUT>        Positive or negative control output from cawlr collapse
    -o, --output <OUTPUT>      Path to resulting pickle file
    -s, --samples <SAMPLES>    Number of samples per kmer to allow [default: 50000]
```

### `cawlr rank`

```bash
$ cawlr help rank
cawlr-rank 
Rank each kmer by the Kulback-Leibler Divergence and between the trained models

USAGE:
    cawlr rank [OPTIONS] --pos-ctrl <POS_CTRL> --neg-ctrl <NEG_CTRL> --output <OUTPUT>

OPTIONS:
    -h, --help                   Print help information
        --neg-ctrl <NEG_CTRL>    Negative control output from cawlr train
    -o, --output <OUTPUT>        Path to output file
        --pos-ctrl <POS_CTRL>    Positive control output from cawlr train
        --samples <SAMPLES>      Ranks are estimated via sampling, higher value for samples means it
                                 takes longer for cawlr rank to run but the ranks will be more
                                 accurate [default: 100000]
        --seed <SEED>            Ranks are estimated via sampling, so to keep values consistent
                                 between subsequent runs a seed value is used [default: 2456]
```

### `cawlr score`

```text
$ cawlr help score
cawlr-score 
Score each kmer with likelihood based on positive and negative controls

USAGE:
    cawlr score [OPTIONS] --input <INPUT> --output <OUTPUT> --pos-ctrl <POS_CTRL> --neg-ctrl <NEG_CTRL> --ranks <RANKS> --genome <GENOME>

OPTIONS:
        --cutoff <CUTOFF>        [default: 10]
    -g, --genome <GENOME>        Path to fasta file for organisms genome, must have a .fai file from
                                 samtools faidx
    -h, --help                   Print help information
    -i, --input <INPUT>          Path to Apache Arrow file from cawlr collapse
    -m, --motif <MOTIF>          
        --neg-ctrl <NEG_CTRL>    Negative control file from cawlr train
    -o, --output <OUTPUT>        Path to output file
        --pos-ctrl <POS_CTRL>    Positive control file from cawlr train
    -r, --ranks <RANKS>          Path to rank file from cawlr rank
```

### `cawlr sma`

## Models

TODO: Point out the models that are provided by `cawlr`

## Example `cawlr` vignette

TODO

## QC Scripts

### `plot_gmm_models.py`

### `plot_read_length_vs_mod_rate.py`

### `plot_scores_with_inference.py`

## Citations
