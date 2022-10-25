use std::{fs, path::PathBuf, process::Command};

use cawlr::utils;
use clap::Parser;

#[derive(Parser)]
struct Args {
    #[clap(short, long)]
    genome: PathBuf,

    #[clap(short, long)]
    pos_fast5s: PathBuf,

    #[clap(short, long)]
    pos_reads: PathBuf,

    #[clap(short, long)]
    neg_fast5s: PathBuf,

    #[clap(short, long)]
    neg_reads: PathBuf,

    #[clap(short, long)]
    output_dir: PathBuf,

    #[clap(long)]
    nanopolish_path: Option<PathBuf>,

    #[clap(long)]
    minimap2_path: Option<PathBuf>,

    #[clap(long)]
    samtools_path: Option<PathBuf>,
}

fn aln_reads(
    minimap2: PathBuf,
    genome: PathBuf,
    reads: PathBuf,
    output: PathBuf,
) -> eyre::Result<()> {
    Command::new(minimap2)
        .arg("-ax")
        .arg("map-ont")
        .arg("--sam-hit-only")
        .arg("--secondary=no")
        .args(&["-t", "4"])
        .arg(genome)
        .arg(reads)
        .arg("-o")
        .arg(output);
    Ok(())
}

fn main() -> eyre::Result<()> {
    let args = Args::parse();
    let nanopolish = utils::find_binary("nanopolish", &args.nanopolish_path)?;
    let minimap2 = utils::find_binary("minimap2", &args.minimap2_path)?;
    let samtools = utils::find_binary("samtools", &args.samtools_path)?;

    fs::create_dir_all(&args.output_dir)?;

    Ok(())
}
