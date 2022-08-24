use std::{
    fs::File,
    io::{self, Read},
    path::{Path, PathBuf},
};

use anyhow::Result;
use clap::{IntoApp, Parser, Subcommand};
use clap_verbosity_flag::Verbosity;
use human_panic::setup_panic;
#[cfg(feature = "mimalloc")]
use mimalloc::MiMalloc;

mod arrow;
mod bkde;
mod collapse;
mod context;
mod index;
mod motif;
mod plus_strand_map;
mod rank;
mod score;
mod sma;
mod train;
mod utils;

use motif::Motif;
use train::Model;
use utils::CawlrIO;

#[cfg(feature = "mimalloc")]
#[global_allocator]
static GLOBAL: MiMalloc = MiMalloc;

#[derive(Parser, Debug)]
#[clap(author, version, about, long_about=None)]
/// Chromatin accessibility with long reads.
struct Args {
    #[clap(short, long)]
    debug: bool,

    #[clap(flatten)]
    verbose: Verbosity,

    #[clap(subcommand)]
    command: Commands,
}

#[derive(Subcommand, Debug)]
enum Commands {
    Collapse {
        /// Path to nanopolish eventalign output with samples column, or stdin
        /// if not provided.
        #[clap(short, long)]
        input: Option<PathBuf>,

        /// Path to BAM alignment file used in nanopolish eventalign
        #[clap(short, long)]
        bam: PathBuf,

        #[clap(short, long)]
        /// Path to output file in Apache Arrow format, defaults to stdout if no
        /// argument provided.
        output: String,

        #[clap(short, long, default_value_t = 2048)]
        /// Number of eventalign records to hold in memory.
        capacity: usize,
    },

    Index {
        #[clap(short, long)]
        input: String,
    },

    /// For each kmer, train a two-component gaussian mixture model and save
    /// models to a file
    Train {
        #[clap(short, long)]
        /// Positive or negative control output from cawlr collapse
        input: String,

        #[clap(short, long)]
        /// Path to resulting pickle file
        output: String,

        #[clap(short, long)]
        /// Path to genome fasta file
        genome: String,

        #[clap(short, long, default_value_t = 50_000)]
        /// Number of samples per kmer to allow
        samples: usize,
    },

    /// Rank each kmer by the Kulback-Leibler Divergence and between the trained
    /// models
    Rank {
        #[clap(long)]
        /// Positive control output from cawlr train
        pos_ctrl: String,

        #[clap(long)]
        /// Negative control output from cawlr train
        neg_ctrl: String,

        #[clap(short, long)]
        /// Path to output file
        output: String,

        #[clap(long, default_value_t = 2456)]
        /// Ranks are estimated via sampling, so to keep values consistent
        /// between subsequent runs a seed value is used
        seed: u64,

        /// Ranks are estimated via sampling, higher value for samples means it
        /// takes longer for cawlr rank to run but the ranks will be more
        /// accurate
        #[clap(long, default_value_t = 100_000_usize)]
        samples: usize,
    },

    /// Score each kmer with likelihood based on positive and negative controls
    Score {
        /// Path to Apache Arrow file from cawlr collapse
        #[clap(short, long)]
        input: String,

        /// Path to output file
        #[clap(short, long)]
        output: String,

        /// Positive control file from cawlr train
        #[clap(long)]
        pos_ctrl: String,

        /// Negative control file from cawlr train
        #[clap(long)]
        neg_ctrl: String,

        /// Path to rank file from cawlr rank
        #[clap(short, long)]
        ranks: String,

        /// Path to fasta file for organisms genome, must have a .fai file from
        /// samtools faidx
        #[clap(short, long)]
        genome: String,

        /// Threshold for current value to be considered reasonable
        #[clap(long, default_value_t = 10.0)]
        cutoff: f64,

        /// Threshold for kmer model to be used
        #[clap(long, default_value_t = 0.05)]
        p_value_threshold: f64,

        /// Only score in kmers that contain this motif, by default will score
        /// all kmers
        #[clap(short, long)]
        motif: Option<Vec<Motif>>,
    },
    Sma {
        /// Path to scored data from cawlr score
        #[clap(short, long)]
        input: String,

        /// Path to output file
        #[clap(short, long)]
        output: Option<String>,

        /// Path to score from positive control dataset
        #[clap(long)]
        pos_ctrl_scores: String,

        /// Path to score from negative control dataset
        #[clap(long)]
        neg_ctrl_scores: String,

        /// Only that contain this motif will be used to perform single molecule
        /// analysis, by default will use all kmers
        #[clap(short, long)]
        motif: Option<Vec<Motif>>,

        /// Number of scores sampled to create kernel density estimate
        #[clap(long, default_value_t = 10_000_usize)]
        kde_samples: usize,

        /// Set seed to have reproducible sampling
        #[clap(long, default_value_t = 2456_u64)]
        seed: u64,
    },
}

fn main() -> Result<()> {
    setup_panic!();
    let args = Args::parse();
    env_logger::Builder::new()
        .filter_level(args.verbose.log_level_filter())
        .init();
    match args.command {
        Commands::Collapse {
            input,
            bam,
            output,
            capacity,
        } => {
            if capacity == 0 {
                let mut cmd = Args::command();
                cmd.error(
                    clap::ErrorKind::InvalidValue,
                    "Capacity must be greater than 0",
                )
                .exit();
            }
            let final_input: Box<dyn Read> = {
                if let Some(path) = input {
                    Box::new(File::open(path)?)
                } else {
                    let stdin = io::stdin().lock();
                    Box::new(stdin)
                }
            };

            let mut collapse = collapse::CollapseOptions::try_new(&bam, &output)?;
            collapse.capacity(capacity).progress(true);
            collapse.run(final_input)?;
        }
        Commands::Index { input } => {
            index::index(input)?;
        }
        Commands::Train {
            input,
            output,
            genome,
            samples,
        } => {
            log::info!("Train command");
            let train = train::Train::try_new(&input, &genome, samples)?;
            let model = train.run()?;
            model.save(output)?;
        }

        Commands::Rank {
            pos_ctrl,
            neg_ctrl,
            output,
            seed,
            samples,
        } => {
            let pos_ctrl_db = Model::load(pos_ctrl)?;
            let neg_ctrl_db = Model::load(neg_ctrl)?;
            let kmer_ranks = rank::RankOptions::new(seed, samples).rank(&pos_ctrl_db, &neg_ctrl_db);
            kmer_ranks.save(output)?;
        }

        Commands::Score {
            input,
            output,
            pos_ctrl,
            neg_ctrl,
            ranks,
            genome,
            cutoff,
            p_value_threshold,
            motif,
        } => {
            let fai_file = format!("{}.fai", genome);
            let fai_file_exists = Path::new(&fai_file).exists();
            if !fai_file_exists {
                let mut cmd = Args::command();
                cmd.error(
                    clap::ErrorKind::MissingRequiredArgument,
                    "Missing .fai index file, run samtools faidx on genome file.",
                )
                .exit();
            }

            motif.iter().for_each(|ms| {
                ms.iter().for_each(|m| {
                    if m.len_motif() > 6 {
                        let mut cmd = Args::command();
                        cmd.error(
                            clap::ErrorKind::InvalidValue,
                            "Length of motif must be less than 6 (size of kmer)",
                        )
                        .exit();
                    }
                })
            });

            log::debug!("Motifs parsed: {motif:?}");
            let scoring = score::ScoreOptions::try_new(
                &pos_ctrl,
                &neg_ctrl,
                &genome,
                &ranks,
                &output,
                cutoff,
                p_value_threshold,
                motif,
            )?;
            scoring.run(input)?;
        }

        Commands::Sma {
            input,
            output,
            pos_ctrl_scores: pos_control_scores,
            neg_ctrl_scores: neg_control_scores,
            motif,
            kde_samples,
            seed,
        } => {
            let mut builder = sma::Builder::new(pos_control_scores, neg_control_scores);
            builder
                .kde_samples(kde_samples)
                .seed(seed)
                .try_motifs(motif)
                .try_output_file(output);
            let sma = builder.build()?;
            sma.run(input)?;
        }
    }
    Ok(())
}
