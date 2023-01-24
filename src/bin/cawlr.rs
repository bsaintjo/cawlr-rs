use std::{
    fs::File,
    io::{self, BufReader, BufWriter, Read},
    path::{Path, PathBuf},
};

use cawlr::{
    arrow::{
        arrow_utils::{load_apply2, load_read_write_arrow},
        eventalign::Eventalign,
        scored_read::ScoredRead,
    },
    bkde::BinnedKde,
    collapse::CollapseOptions,
    filter::{FilterOptions, Region},
    index,
    motif::{all_bases, Motif},
    npsmlr::{self, train::TrainOptions},
    rank::RankOptions,
    score::ScoreOptions,
    score_model,
    sma::SmaOptions,
    train::{self, Model, Train, TrainStrategy},
    utils::{self, CawlrIO},
};
use clap::{error::ErrorKind, CommandFactory, Parser, Subcommand};
use clap_verbosity_flag::Verbosity;
use eyre::Result;
use human_panic::setup_panic;
#[cfg(feature = "mimalloc")]
use mimalloc::MiMalloc;

#[cfg(feature = "mimalloc")]
#[global_allocator]
static GLOBAL: MiMalloc = MiMalloc;

fn parse_strategy(src: &str) -> Result<TrainStrategy, String> {
    match src {
        "all" => Ok(TrainStrategy::AllSamples),
        "avg" => Ok(TrainStrategy::AvgSample),
        _ => Err(String::from("Invalid strategy: either 'avg' or 'all'")),
    }
}

#[derive(Debug, Subcommand)]
enum QCCmd {
    Score {
        #[clap(short, long)]
        input: PathBuf,
    },

    Eventalign {
        #[clap(short, long)]
        input: PathBuf,
    },
}

#[derive(Debug, Subcommand)]
enum FilterCmd {
    Score {
        #[clap(short, long)]
        input: PathBuf,

        /// Arrow file output
        #[clap(short, long)]
        output: PathBuf,

        #[clap(short, long, num_args = 1..)]
        region: Vec<Region>,
    },

    Eventalign {
        #[clap(short, long)]
        input: PathBuf,

        /// Arrow file output
        #[clap(short, long)]
        output: PathBuf,

        #[clap(short, long, num_args = 1..)]
        region: Vec<Region>,
    },
}

#[derive(Debug, Subcommand)]
enum NpsmlrCmd {
    /// Train using algorithm adapted from NP-SMLR
    Train {
        /// Input arrow file, usually from cawlr collapse
        #[clap(short, long)]
        input: PathBuf,

        /// Pickle file containing model parameters
        #[clap(short, long)]
        output: PathBuf,

        /// Only train on kmers containing these motifs, can speed up training
        /// time
        #[clap(short, long)]
        motif: Vec<Motif>,

        /// Number of samples to use to train GMM
        #[clap(long, default_value_t = 50000)]
        samples: usize,

        /// Train a single component GMM (ie fit a single Gaussian)
        #[clap(long)]
        single: bool,

        /// Filter outliers with DBSCAN algorithm
        #[clap(long)]
        dbscan: bool,

        #[clap(long)]
        db_path: Option<PathBuf>,
    },

    /// Score using algorithm adapted from NP-SMLR
    Score {
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
        #[clap(short, long, num_args(1..))]
        motif: Vec<Motif>,

        /// Values less than -cutoff for the positive and negative control will
        /// be filtered
        #[clap(short, long, default_value_t = 10.0)]
        cutoff: f64,

        /// If an events has more than freq_thresh samples, it will be filtered
        #[clap(short, long, default_value_t = 10)]
        freq_thresh: usize,
    },
}

#[derive(Parser, Debug)]
#[clap(author, version, about, long_about=None)]
/// Chromatin accessibility with long reads.
struct Args {
    #[clap(flatten)]
    verbose: Verbosity,

    #[clap(subcommand)]
    command: Commands,
}

#[derive(Subcommand, Debug)]
enum Commands {
    #[clap(subcommand)]
    QC(QCCmd),

    #[clap(subcommand)]
    Npsmlr(NpsmlrCmd),

    /// Preprocess nanopolish eventalign output
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
        output: Option<PathBuf>,

        #[clap(short, long, default_value_t = 2048)]
        /// Number of eventalign records to hold in memory.
        capacity: usize,
    },

    /// Create bed file of the reads in the Arrow file
    ///
    /// Output file will be named {input}.idx.bed
    Index {
        /// Arrow file from collapse or score
        #[clap(short, long)]
        input: PathBuf,
    },

    /// Filter Arrow output file based on genomic coordinates
    #[clap(subcommand)]
    Filter(FilterCmd),

    /// For each kmer, train a two-component gaussian mixture model and save
    /// models to a file
    Train {
        /// Positive or negative control output from cawlr collapse
        #[clap(short, long)]
        input: PathBuf,

        /// Path to resulting pickle file
        #[clap(short, long)]
        output: PathBuf,

        /// Path to genome fasta file
        #[clap(short, long)]
        genome: PathBuf,

        /// Number of samples per kmer to allow
        #[clap(short, long, default_value_t = 50_000)]
        samples: usize,

        /// Number of threads to use for training, by default num cpus
        #[clap(short = 'j', long)]
        num_threads: Option<usize>,

        /// Pick what data is used to train models
        ///
        /// Either train on individual samples using "all" or just the average
        /// using "avg"
        #[clap(long, default_value_t = TrainStrategy::AllSamples, value_parser=parse_strategy)]
        strategy: train::TrainStrategy,
    },

    /// Rank each kmer by the Kulback-Leibler Divergence and between the trained
    /// models
    Rank {
        /// Positive control output from cawlr train
        #[clap(long)]
        pos_ctrl: PathBuf,

        /// Negative control output from cawlr train
        #[clap(long)]
        neg_ctrl: PathBuf,

        /// Path to output file
        #[clap(short, long)]
        output: PathBuf,

        /// Ranks are estimated via sampling, so to keep values consistent
        /// between subsequent runs a seed value is used
        #[clap(long, default_value_t = 2456)]
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
        input: PathBuf,

        /// Path to output file
        #[clap(short, long)]
        output: PathBuf,

        /// Positive control file from cawlr train
        #[clap(long)]
        pos_ctrl: PathBuf,

        /// Negative control file from cawlr train
        #[clap(long)]
        neg_ctrl: PathBuf,

        /// Path to rank file from cawlr rank
        #[clap(short, long)]
        ranks: PathBuf,

        /// Path to fasta file for organisms genome, must have a .fai file from
        /// samtools faidx
        #[clap(short, long)]
        genome: PathBuf,

        /// Threshold for current value to be considered reasonable
        #[clap(long, default_value_t = 10.0)]
        cutoff: f64,

        /// Threshold for kmer model to be used
        #[clap(long, default_value_t = 0.05)]
        p_value_threshold: f64,

        /// Only score in kmers that contain this motif, by default will score
        /// all kmers. Format = "{position of modified base}:{motif}", ie "2:GC"
        /// if the C in GC is the modified base.
        #[clap(short, long)]
        motif: Option<Vec<Motif>>,
    },
    /// Compute kernel density estimate of control score data
    ModelScores {
        /// Arrow output from cawlr score
        #[clap(short, long)]
        input: PathBuf,

        /// Pickle file containing estimated kernel density estimate values
        #[clap(short, long)]
        output: PathBuf,

        /// Number of bins used to estimate the kernel density estimate
        #[clap(short, long, default_value_t = 10_000)]
        bins: u32,

        /// Number of scores sampled from the input to compute the kernel
        /// density estimate
        #[clap(short, long, default_value_t = 10_000)]
        samples: usize,
    },
    /// Infer nucleosome positions on single molecules
    Sma {
        /// Path to scored data from cawlr score
        #[clap(short, long)]
        input: PathBuf,

        /// Path to output file
        #[clap(short, long)]
        output: Option<PathBuf>,

        /// Output from cawlr model-scores for treated control sample
        #[clap(long)]
        pos_ctrl_scores: PathBuf,

        /// Output from cawlr model-scores for untreated control sample
        #[clap(long)]
        neg_ctrl_scores: PathBuf,

        /// Only that contain this motif will be used to perform single molecule
        /// analysis, by default will use all kmers
        #[clap(short, long)]
        motif: Option<Vec<Motif>>,
    },
}

fn main() -> Result<()> {
    setup_panic!();
    jane_eyre::install()?;

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
                cmd.error(ErrorKind::InvalidValue, "Capacity must be greater than 0")
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

            let final_output = utils::stdout_or_file(output.as_ref())?;
            let final_output = BufWriter::new(final_output);

            let mut collapse = CollapseOptions::from_writer(final_output, &bam)?;
            collapse.capacity(capacity).progress(true);
            collapse.run(final_input)?;
        }
        Commands::Index { input } => {
            index::index(input)?;
        }
        Commands::Filter(FilterCmd::Eventalign {
            input,
            output,
            region,
        }) => {
            let filters = FilterOptions::new(region);
            let reader = File::open(input)?;
            let writer = File::create(output)?;
            load_read_write_arrow(reader, writer, |xs: Vec<Eventalign>| {
                Ok(xs.into_iter().filter(|x| filters.any_valid(x)).collect())
            })?;
        }

        Commands::Filter(FilterCmd::Score {
            input,
            output,
            region,
        }) => {
            let filters = FilterOptions::new(region);
            let reader = File::open(input)?;
            let writer = File::create(output)?;
            load_read_write_arrow(reader, writer, |xs: Vec<ScoredRead>| {
                Ok(xs.into_iter().filter(|x| filters.any_valid(x)).collect())
            })?;
        }

        Commands::Train {
            input,
            output,
            genome,
            samples,
            strategy,
            num_threads,
        } => {
            log::info!("Train command");
            let mut n_logical_cores = num_cpus::get();

            if let Some(n) = num_threads {
                rayon::ThreadPoolBuilder::new()
                    .num_threads(n)
                    .build_global()?;
                n_logical_cores = n;
            }

            log::info!("Using {n_logical_cores} logical cores");
            log::info!("Using strategy: {strategy}");
            let train = Train::try_new(input, genome, samples, strategy)?;
            let model = train.run()?;
            model.save_as(output)?;
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
            let kmer_ranks = RankOptions::new(seed, samples).rank(&pos_ctrl_db, &neg_ctrl_db);
            kmer_ranks.save_as(output)?;
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
            let fai_file = format!("{}.fai", genome.display());
            let fai_file = Path::new(&fai_file);
            log::debug!("fasta index file filename: {fai_file:?}");
            if !fai_file.exists() {
                let mut cmd = Args::command();
                cmd.error(
                    ErrorKind::MissingRequiredArgument,
                    "Missing .fai index file, run samtools faidx on genome file.",
                )
                .exit();
            }

            motif.iter().for_each(|ms| {
                ms.iter().for_each(|m| {
                    if m.len_motif() > 6 {
                        let mut cmd = Args::command();
                        cmd.error(
                            ErrorKind::InvalidValue,
                            "Length of motif must be less than 6 (size of kmer)",
                        )
                        .exit();
                    }
                })
            });

            log::debug!("Motifs parsed: {motif:?}");
            let mut scoring =
                ScoreOptions::try_new(&pos_ctrl, &neg_ctrl, &genome, &ranks, &output)?;
            scoring.cutoff(cutoff).p_value_threshold(p_value_threshold);
            if let Some(motifs) = motif {
                scoring.motifs(motifs);
            }
            scoring.run(input)?;
        }

        Commands::ModelScores {
            input,
            output,
            bins,
            samples,
        } => {
            let file = File::open(input)?;
            let bkde = score_model::Options::default()
                .bins(bins)
                .samples(samples)
                .run(file)?;
            bkde.save_as(output)?;
        }

        Commands::Sma {
            input,
            output,
            pos_ctrl_scores,
            neg_ctrl_scores,
            motif,
        } => {
            let pos_bkde = BinnedKde::load(pos_ctrl_scores)?;
            let neg_bkde = BinnedKde::load(neg_ctrl_scores)?;
            let writer = utils::stdout_or_file(output.as_ref())?;
            let motifs = {
                if let Some(motif) = motif {
                    motif
                } else {
                    all_bases()
                }
            };
            let mut sma = SmaOptions::new(pos_bkde, neg_bkde, motifs, writer);
            if let Some(output_filename) = output {
                let track_name = output_filename
                    .file_name()
                    .ok_or_else(|| eyre::eyre!("Not a filename"))?
                    .to_str()
                    .unwrap();
                sma.track_name(track_name);
            }
            sma.run(input)?;
        }
        Commands::QC(cmd) => match cmd {
            QCCmd::Score { input } => {
                let reader = BufReader::new(File::open(input)?);
                load_apply2(reader, |_xs: ScoredRead| Ok(()))?;
            }
            QCCmd::Eventalign { input } => {
                let reader = BufReader::with_capacity(1024 * 32, File::open(input)?);
                load_apply2(reader, |_xs: Eventalign| Ok(()))?;
            }
        },

        Commands::Npsmlr(cmd) => match cmd {
            NpsmlrCmd::Train {
                input,
                output,
                motif,
                samples,
                single,
                dbscan,
                db_path,
            } => {
                let reader = BufReader::new(File::open(input)?);
                let writer = File::create(output)?;
                TrainOptions::default()
                    .n_samples(samples)
                    .db_path(db_path)
                    .single(single)
                    .dbscan(dbscan)
                    .motifs(motif)
                    .run(reader, writer)?;
            }
            NpsmlrCmd::Score {
                input,
                pos_ctrl,
                neg_ctrl,
                ranks,
                output,
                cutoff,
                freq_thresh,
                motif,
            } => {
                let reader = BufReader::new(File::open(input)?);
                let writer = File::create(output)?;
                let mut score_options = npsmlr::ScoreOptions::load(pos_ctrl, neg_ctrl, ranks)?;
                score_options
                    .freq_thresh(freq_thresh)
                    .cutoff(cutoff)
                    .motifs(motif)
                    .run(reader, writer)?;
            }
        },
    }
    Ok(())
}
