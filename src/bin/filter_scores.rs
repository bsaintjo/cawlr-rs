use std::{fs::File, path::PathBuf};

use anyhow::Result;
use cawlr::arrow::{load_apply, save, wrap_writer, MetadataExt, ScoredRead};
use clap::Parser;

#[derive(Parser)]
struct Args {
    /// Arrow input file from cawlr score
    #[clap(short, long, required = true)]
    input: PathBuf,

    /// Arrow output file path
    #[clap(short, long)]
    output: PathBuf,

    #[clap(long)]
    chrom: Option<String>,

    #[clap(long)]
    start: Option<usize>,

    #[clap(long)]
    stop: Option<usize>,

    // Threshold for score to be modified
    #[clap(long, default_value_t = 0.0)]
    modification_threshold: f64,

    // Minimum level of modification allowed
    #[clap(long, default_value_t = 0.0)]
    percent_modified_min: f64,

    // Maximum level of modification allowed
    #[clap(long, default_value_t = 100.0)]
    percent_modified_max: f64,

    // Minimum read length allowed
    #[clap(long, default_value_t = 0)]
    read_length_min: u64,

    // Maximum read length allowed
    #[clap(long, default_value_t=u64::MAX)]
    read_length_max: u64,
}

fn percent_mod(read: &ScoredRead, threshold: f64) -> f64 {
    let above = read
        .scores()
        .iter()
        .filter(|s| s.score() >= threshold)
        .count() as f64;
    above / read.length() as f64
}

fn filter_by(args: &Args, read: &ScoredRead) -> bool {
    let pmod = percent_mod(read, args.modification_threshold);
    read.length() >= args.read_length_min
        && read.length() < args.read_length_max
        && pmod >= args.percent_modified_min
        && pmod <= args.percent_modified_max
}

fn main() -> Result<()> {
    let args = Args::parse();

    let reader = File::open(&args.input)?;
    let writer = File::create(&args.output)?;
    let schema = ScoredRead::schema();
    let mut writer = wrap_writer(writer, &schema)?;

    load_apply(reader, |reads: Vec<ScoredRead>| {
        let reads = reads
            .into_iter()
            .filter(|read| filter_by(&args, read))
            .collect::<Vec<_>>();
        save(&mut writer, &reads)?;
        Ok(())
    })?;

    Ok(())
}
