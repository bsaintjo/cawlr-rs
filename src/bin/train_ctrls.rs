use std::{
    path::{Path, PathBuf},
    process::{Command, Stdio},
};

use cawlr::utils::find_binary;
use clap::Parser;

fn collapse_piped(np_bin: &Path, reads: &Path, bam: &Path, genome: &Path) -> eyre::Result<()> {
    let nanopolish = Command::new(np_bin)
        .arg("eventalign")
        .arg("--reads")
        .arg(reads)
        .arg("--bam")
        .arg(bam)
        .arg("--genome")
        .arg(genome)
        .args(&["--scale-events", "--print-read-names"])
        .stdout(Stdio::piped())
        .spawn()?;
    let _cawlr = Command::new("cawlr")
        .arg("collapse")
        .arg("-b")
        .arg(bam)
        .stdin(nanopolish.stdout.unwrap())
        .output()?;
    todo!()
}

#[derive(Parser)]
struct Args {
    #[clap(short, long)]
    name: String,

    #[clap(long)]
    pos_fast5s: PathBuf,

    #[clap(long)]
    pos_reads: PathBuf,

    #[clap(long)]
    neg_fast5s: PathBuf,

    #[clap(long)]
    neg_reads: PathBuf,

    #[clap(long)]
    samtools_filepath: Option<PathBuf>,

    #[clap(long)]
    minimap2_filepath: Option<PathBuf>,

    #[clap(long)]
    nanopolish_filepath: Option<PathBuf>,

    #[clap(short, long)]
    output_dir: PathBuf,
}

fn main() -> eyre::Result<()> {
    let args = Args::parse();

    let samtools = find_binary("samtools", &args.samtools_filepath)?;
    let nanopolish = find_binary("nanopolish", &args.nanopolish_filepath)?;
    let minimap2 = find_binary("minimap2", &args.minimap2_filepath)?;
    todo!()
}
