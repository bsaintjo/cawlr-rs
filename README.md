# cawlr 𓅨

[![License](https://img.shields.io/badge/license-BSD_3--Clause-informational)](./LICENSE)

**C**hromatin **a**ccesibility **w**ith **l**ong **r**eads (`cawlr`) is a tool for detecting accessible regions of chromatin on single molecules. The tool works with nanopore data via [`nanopolish`](https://github.com/jts/nanopolish) or PacBio data via [`kineticsTools`](https://github.com/PacificBiosciences/kineticsTools).

`cawlr` is flexible and can work with data that uses various modifications for chromatin accessibility detection. Outputs of all `cawlr` modules are in parquet format can be manipulated via any programming language that has an [Apache Arrow](https://arrow.apache.org/install/) library. Furthermore, `cawlr` includes various scripts for plotting and evaluating the results.

- [cawlr 𓅨](#cawlr-𓅨)
  - [Installation](#installation)
    - [Installing rust via rustup](#installing-rust-via-rustup)
    - [Installing cawlr](#installing-cawlr)
      - [via `cargo install` (Not yet uploaded)](#via-cargo-install-not-yet-uploaded)
      - [Latest from git](#latest-from-git)
      - [(Optional) Run tests](#optional-run-tests)
    - [Python dependencies](#python-dependencies)
  - [Running nanopolish](#running-nanopolish)
  - [Usage](#usage)
    - [`cawlr preprocess`](#cawlr-preprocess)
    - [`cawlr train`](#cawlr-train)
    - [`cawlr rank`](#cawlr-rank)
    - [`cawlr score`](#cawlr-score)
    - [`cawlr sma`](#cawlr-sma)
  - [QC Scripts](#qc-scripts)
    - [`plot_gmm_models.py`](#plot_gmm_modelspy)
    - [`plot_read_length_vs_mod_rate.py`](#plot_read_length_vs_mod_ratepy)
    - [`plot_scores_with_inference.py`](#plot_scores_with_inferencepy)
  - [Citations](#citations)

## Installation

### Installing rust via rustup

It is recommended to use rustup to install the tools for compiling rust code at [rustup.rs](https://rustup.rs/).

### Installing cawlr

#### via `cargo install` (Not yet uploaded)

```bash
cargo install cawlr
```

#### Latest from git

```bash
git clone https://github.com/bsaintjo/cawlr-rs.git
cd cawlr-rs
cargo install --path .
```

#### (Optional) Run tests

```bash
cargo test
```

### Python dependencies

## Running nanopolish

Example command of running nanopolish to set up inputs

```bash
nanopolish index -d path/to/fast5s reads.fasta
minimap2
```

## Usage

### `cawlr preprocess`

```bash
$ cawlr help preprocess
cawlr-preprocess 
Calculates mean per-read per-position and optionally filters data based on a given region

USAGE:
    cawlr preprocess [OPTIONS] --input <INPUT> --bam <BAM> --genome <GENOME> --output <OUTPUT>

OPTIONS:
    -b, --bam <BAM>          path to nanopolish eventalign output with samples column
    -c, --chrom <CHROM>      output only includes data from this chromosome
    -g, --genome <GENOME>    path to genome file
    -h, --help               Print help information
    -i, --input <INPUT>      path to nanopolish eventalign output with samples column
    -o, --output <OUTPUT>    path to output file in parquet format
        --start <START>      output only includes data that aligns at or after this position, should
                             be set with --chrom TODO: Throw error if set without --chrom
        --stop <STOP>        output only includes data that aligns at or before this position,
                             should be set with --chrom TODO: Throw error if set without --chrom

```

### `cawlr train`

```bash
$ cawlr help train
cawlr-train 
For each kmer, train a two-component gaussian mixture model and save models to a file

USAGE:
    cawlr train --input <INPUT> --output <OUTPUT>

OPTIONS:
    -h, --help               Print help information
    -i, --input <INPUT>      Parquet file of positive or negative control from cawlr preprocess
    -o, --output <OUTPUT>    
```

### `cawlr rank`

```bash
$ cawlr help rank
cawlr-rank 
Rank each kmer by the symmetrical Kulback-Leibler Divergence and output results

USAGE:
    cawlr rank [OPTIONS] --pos-ctrl <POS_CTRL> --neg-ctrl <NEG_CTRL> --output <OUTPUT>

OPTIONS:
    -h, --help                   Print help information
        --neg-ctrl <NEG_CTRL>    
    -o, --output <OUTPUT>        
        --pos-ctrl <POS_CTRL>    
        --samples <SAMPLES>      [default: 10000]
        --seed <SEED>            [default: 2456]
```

### `cawlr score`

```text
$ cawlr help score
cawlr-score 

USAGE:
    cawlr score --input <INPUT> --output <OUTPUT> --pos-ctrl <POS_CTRL> --neg-ctrl <NEG_CTRL> --ranks <RANKS> --genome <GENOME> --cutoff <CUTOFF>

OPTIONS:
        --cutoff <CUTOFF>        
    -g, --genome <GENOME>        
    -h, --help                   Print help information
    -i, --input <INPUT>          
        --neg-ctrl <NEG_CTRL>    
    -o, --output <OUTPUT>        
        --pos-ctrl <POS_CTRL>    
    -r, --ranks <RANKS>          

```

### `cawlr sma`

## QC Scripts

### `plot_gmm_models.py`

### `plot_read_length_vs_mod_rate.py`

### `plot_scores_with_inference.py`

## Citations
