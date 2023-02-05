use std::{fs::File, io::BufReader, path::PathBuf};

use clap::Parser;
use libcawlr::{motif::Motif, npsmlr};

#[derive(Parser, Debug)]
pub struct ScoreCmd {
    /// Input arrow file, usually from cawlr collapse
    #[clap(short, long)]
    input: PathBuf,

    /// Path to positive control model, usually from cawlr train
    #[clap(short, long)]
    pos_ctrl: PathBuf,

    /// Path to negative control model, usually from cawlr train
    #[clap(short, long)]
    neg_ctrl: PathBuf,

    /// Path to ranks file, usually from cawlr rank
    #[clap(short, long)]
    ranks: PathBuf,

    /// Output arrow file of scored reads
    #[clap(short, long)]
    output: PathBuf,

    /// Motifs to score on, at least 1 motif must be provided
    #[clap(short, long, required=true, num_args=1.., value_delimiter=',')]
    motif: Vec<Motif>,

    /// Values less than -cutoff for the positive and negative control will
    /// be filtered
    #[clap(short, long, default_value_t = 10.0)]
    cutoff: f64,

    /// If an events has more than freq_thresh samples, it will be filtered
    #[clap(short, long, default_value_t = 10)]
    freq_thresh: usize,
}

impl ScoreCmd {
    pub fn run(self) -> eyre::Result<()> {
        let reader = BufReader::new(File::open(self.input)?);
        let writer = File::create(self.output)?;
        let mut score_options =
            npsmlr::ScoreOptions::load(self.pos_ctrl, self.neg_ctrl, self.ranks)?;
        score_options
            .freq_thresh(self.freq_thresh)
            .cutoff(self.cutoff)
            .motifs(self.motif)
            .run(reader, writer)
    }
}
