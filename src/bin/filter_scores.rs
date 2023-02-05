use std::{fs::File, path::PathBuf};

use libcawlr::arrow::{
    arrow_utils::{load_read_write, wrap_writer},
    metadata::MetadataExt,
    scored_read::{Score, ScoredRead},
};
use clap::Parser;
use eyre::Result;

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

    // Score must be greater than or equal to this value to count as modified.
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

impl Default for Args {
    fn default() -> Self {
        Args {
            input: PathBuf::default(),
            output: PathBuf::default(),
            chrom: None,
            start: None,
            stop: None,
            modification_threshold: 0.0,
            percent_modified_min: 0.0,
            percent_modified_max: 100.0,
            read_length_min: 0,
            read_length_max: u64::MAX,
        }
    }
}

fn percent_mod(scores: &[Score], threshold: f64) -> f64 {
    if scores.is_empty() {
        0.0
    } else {
        let above = scores.iter().filter(|s| s.score >= threshold).count() as f64;
        above / scores.len() as f64
    }
}

fn filter_by(args: &Args, read: &ScoredRead) -> bool {
    let pmod = percent_mod(read.scores(), args.modification_threshold);
    read.seq_length() >= args.read_length_min
        && read.seq_length() < args.read_length_max
        && pmod >= args.percent_modified_min
        && pmod <= args.percent_modified_max
}

fn main() -> Result<()> {
    let args = Args::parse();

    let reader = File::open(&args.input)?;
    let writer = File::create(&args.output)?;
    let schema = ScoredRead::schema();
    let writer = wrap_writer(writer, &schema)?;

    load_read_write(reader, writer, |reads: Vec<ScoredRead>| {
        let reads = reads
            .into_iter()
            .filter(|read| filter_by(&args, read))
            .collect();
        Ok(reads)
    })?;

    Ok(())
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_percent_mod_empty() {
        let read = ScoredRead::default();
        let pmod = percent_mod(read.scores(), 0.0);
        assert_eq!(pmod, 0.0f64)
    }

    #[test]
    fn test_percent_mod() {
        let scores = [10.0, 20.0, 30.0]
            .into_iter()
            .map(|s| Score {
                score: s,
                ..Default::default()
            })
            .collect::<Vec<_>>();
        let pmod = percent_mod(&scores, 25.0);
        assert_eq!(pmod, (1. / 3.));

        // let scores = [10.0, 20.0, 30.0].into_iter().map(|s|
        // ScoreBuilder::default().score(s).build()).collect::<Vec<_>>();
        let pmod = percent_mod(&scores, 40.0);
        assert_eq!(pmod, 0.0f64);

        let pmod = percent_mod(&scores, 0.0);
        assert_eq!(pmod, 1.0f64);

        let scores = [0.0, 100.0]
            .into_iter()
            .map(|s| Score {
                score: s,
                ..Default::default()
            })
            .collect::<Vec<_>>();

        let pmod = percent_mod(&scores, 0.0);
        assert_eq!(pmod, 1.0f64);

        let pmod = percent_mod(&scores, 100.0);
        assert_eq!(pmod, 0.5f64);
    }
}
