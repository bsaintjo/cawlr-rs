use std::{
    fs,
    path::PathBuf,
    process::{Command, Stdio},
};

use cawlr::utils;
use clap::Parser;
use eyre::Result;

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

fn np_index(nanopolish: PathBuf, fast5s: PathBuf, reads: PathBuf) -> Result<()> {
    let mut cmd = Command::new(nanopolish);
    cmd.arg("index").arg("-d").arg(fast5s).arg(reads);
    log::info!("{cmd:?}");
    cmd.output()?;
    Ok(())
}

fn aln_reads(
    minimap2: PathBuf,
    samtools: PathBuf,
    genome: PathBuf,
    reads: PathBuf,
    output: PathBuf,
) -> eyre::Result<()> {
    let map_cmd = Command::new(minimap2)
        .arg("-ax")
        .arg("map-ont")
        .arg("--sam-hit-only")
        .arg("--secondary=no")
        .args(&["-t", "4"])
        .arg(genome)
        .arg(reads)
        .stdout(Stdio::piped())
        .spawn()?;
    let sam_cmd = Command::new(samtools)
        .arg("sort")
        .arg("--write-index")
        .arg("-T")
        .arg("reads.tmp")
        .arg("-o")
        .arg(output)
        .stdin(map_cmd.stdout.unwrap())
        .output()?;
    Ok(())
}

fn eventalign_collapse(nanopolish: PathBuf, reads: PathBuf, bam: PathBuf) -> Result<()> {
    let cmd = Command::new(nanopolish)
        .arg("eventalign")
        .arg("-r")
        .arg(reads)
        .arg("-b")
        .arg(bam)
        .arg("-t")
        .arg("4")
        .arg("--print-read-names")
        .arg("--samples")
        .stdout(Stdio::piped())
        .spawn()?;
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
