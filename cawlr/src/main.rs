mod cmd;
mod file;
mod pipeline;

use std::{
    fs::File,
    io::BufReader,
    path::{Path, PathBuf},
};

use clap::{error::ErrorKind, CommandFactory, Parser, Subcommand};
use clap_verbosity_flag::Verbosity;
use eyre::Result;
use file::ValidPathBuf;
use human_panic::setup_panic;
use libcawlr::{
    arrow::{
        arrow_utils::{load_apply2, load_read_write_arrow},
        eventalign::Eventalign,
        io::ModFile,
        scored_read::ScoredRead,
    },
    bkde::BinnedKde,
    filter::FilterOptions,
    index,
    motif::{all_bases, Motif},
    rank::RankOptions,
    region::Region,
    score::ScoreOptions,
    score_model,
    sma::SmaOptions,
    train::{self, Model, Train, TrainStrategy},
    utils::{self, CawlrIO},
};
#[cfg(feature = "mimalloc")]
use mimalloc::MiMalloc;
use pipeline::PipelineCmds;

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
    Train(cmd::train::TrainCmd),

    /// Score using algorithm adapted from NP-SMLR
    Score(cmd::score::ScoreCmd),
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

    /// Pipelines for running multiple commands at once
    #[clap(subcommand)]
    Pipeline(PipelineCmds),

    /// Preprocess nanopolish eventalign output
    Collapse(cmd::collapse::CollapseCmd),

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
        pos_ctrl: ValidPathBuf,

        /// Negative control output from cawlr train
        #[clap(long)]
        neg_ctrl: ValidPathBuf,

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
        input: ValidPathBuf,

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

        /// Bam tag to use for modification detection. This is only used if the
        /// input is a BAM file, usually as input from another tool. This is on
        /// the MM tag in the bam file with typical format such as C+m
        /// for methylation on the top strand. For more information, see
        /// section 1.7 of the Sequence Alignment/Map Optional Fields
        /// Specification link: https://samtools.github.io/hts-specs/SAMtags.pdf
        #[clap(short, long)]
        tag: Option<Vec<u8>>,
    },
    /// Infer nucleosome positions on single molecules
    Sma {
        /// Path to scored data from cawlr score
        #[clap(short, long)]
        input: ValidPathBuf,

        /// Path to output file
        #[clap(short, long)]
        output: Option<PathBuf>,

        /// Output from cawlr model-scores for treated control sample
        #[clap(long)]
        pos_ctrl_scores: ValidPathBuf,

        /// Output from cawlr model-scores for untreated control sample
        #[clap(long)]
        neg_ctrl_scores: ValidPathBuf,

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
    let log_level_filter = args.verbose.log_level_filter();

    match args.command {
        Commands::Collapse(cmd) => cmd.run()?,
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
            tag,
        } => {
            let mod_file = ModFile::open_path(input, tag)?;
            let bkde = score_model::Options::default()
                .bins(bins)
                .samples(samples)
                .run_modfile(mod_file)?;
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
            NpsmlrCmd::Train(cmd) => cmd.run()?,
            NpsmlrCmd::Score(cmd) => cmd.run()?,
        },
        Commands::Pipeline(plcmd) => plcmd.run(log_level_filter)?,
    }
    Ok(())
}
