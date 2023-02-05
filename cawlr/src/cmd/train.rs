use std::{fs::File, io::BufReader, path::PathBuf};

use clap::Parser;
use libcawlr::{
    motif::{all_bases, Motif},
    npsmlr::train::TrainOptions,
};

#[derive(Debug, Parser)]
pub struct TrainCmd {
    /// Input arrow file, usually from cawlr collapse
    #[clap(short, long)]
    pub input: PathBuf,

    /// Pickle file containing model parameters
    #[clap(short, long)]
    pub output: PathBuf,

    /// Number of samples to use to train GMM
    #[clap(long, default_value_t = 50000)]
    pub samples: usize,

    /// Train a single component GMM (ie fit a single Gaussian)
    #[clap(long)]
    pub single: bool,

    /// Filter outliers with DBSCAN algorithm
    #[clap(long)]
    pub dbscan: bool,

    /// Path to SQLite database used for storing training data,
    /// otherwise created in temporary file and removed after completion
    #[clap(long)]
    pub db_path: Option<PathBuf>,

    /// Only train on kmers containing these motifs, can speed up training
    /// time
    #[clap(short, long, value_delimiter = ',')]
    pub motif: Vec<Motif>,
}

impl TrainCmd {
    pub fn run(mut self) -> eyre::Result<()> {
        log::info!("Train command");
        let reader = BufReader::new(File::open(self.input)?);
        let writer = File::create(self.output)?;
        if self.motif.is_empty() {
            log::info!("No motifs found, will train on all motifs");
            self.motif = all_bases();
        }
        TrainOptions::default()
            .n_samples(self.samples)
            .db_path(self.db_path)
            .single(self.single)
            .dbscan(self.dbscan)
            .motifs(self.motif)
            .run(reader, writer)?;
        Ok(())
    }
}
