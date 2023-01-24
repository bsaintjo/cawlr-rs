use std::{
    fs::{self, File},
    io::BufReader,
    path::{Path, PathBuf},
    process::{Command, Stdio},
};

use cawlr::{
    collapse::CollapseOptions,
    motif::Motif,
    npsmlr::{train::TrainOptions, ScoreOptions},
    rank::RankOptions,
    score_model::Options,
    train::Model,
    utils::{self, wrap_cmd, wrap_cmd_output, CawlrIO},
};
use clap::Parser;
use eyre::Result;
use fnv::FnvHashMap;
use log::LevelFilter;

#[derive(Parser)]
struct Args {
    /// Path to genome fasta file
    #[clap(short, long)]
    genome: PathBuf,

    /// Directory containing fast5s for positive control
    #[clap(long)]
    pos_fast5s: PathBuf,

    /// Path to single fasta/q file or directory of fasta/q of reads from the positive control
    #[clap(long)]
    pos_reads: PathBuf,

    /// Optional path to sequencing_summary.txt file for positive control,
    /// speeds up nanopolish indexing
    #[clap(long)]
    pos_summary: Option<PathBuf>,

    /// Directory containing fast5s for negative control
    #[clap(long)]
    neg_fast5s: PathBuf,

    /// Path to single fasta/q file or directory of fasta/q of reads from the negative control
    #[clap(long)]
    neg_reads: PathBuf,

    /// Optional path to sequencing_summary.txt file for negative control,
    /// speeds up nanopolish indexing
    #[clap(long)]
    neg_summary: Option<PathBuf>,

    /// Output directory for pipeline
    #[clap(short, long)]
    output_dir: PathBuf,

    /// Path to nanopolish tool, optional if in docker container or in PATH
    #[clap(long)]
    nanopolish_path: Option<PathBuf>,

    /// Path to minimap2 tool, optional if in docker container or in PATH
    #[clap(long)]
    minimap2_path: Option<PathBuf>,

    /// Path to samtools tool, optional if in docker container or in PATH
    #[clap(long)]
    samtools_path: Option<PathBuf>,

    /// Number of threads to use
    #[clap(short = 'j', long, default_value_t = 4)]
    n_threads: usize,

    /// Motifs of modification to filter on, format is "{position}:{motif}"
    /// ie for GpC motif, motif is "2:GC"
    #[clap(short, long, num_args=1..)]
    motifs: Vec<Motif>,
}

fn np_index(
    nanopolish: &Path,
    fast5s: &Path,
    reads: &Path,
    summary: &Option<PathBuf>,
) -> Result<()> {
    let mut cmd = Command::new(nanopolish);
    cmd.arg("index").arg("-d").arg(fast5s);
    if let Some(summary) = summary {
        cmd.arg("-s").arg(summary);
    }
    cmd.arg(reads);
    log::info!("{cmd:?}");
    cmd.output()?;
    Ok(())
}

fn aln_reads(
    minimap2: &Path,
    samtools: &Path,
    genome: &Path,
    reads: &Path,
    output: &Path,
) -> eyre::Result<()> {
    let mut map_cmd = Command::new(minimap2);
    map_cmd
        .arg("-ax")
        .arg("map-ont")
        .arg("--sam-hit-only")
        .arg("--secondary=no")
        .args(["-t", "4"])
        .arg(genome)
        .arg(reads)
        .stdout(Stdio::piped());
    log::info!("{map_cmd:?}");
    let map_output = map_cmd.spawn()?;

    let mut sam_cmd = Command::new(samtools);
    sam_cmd
        .arg("sort")
        .arg("--write-index")
        .arg("-T")
        .arg("reads.tmp")
        .arg("-o")
        .arg(output)
        .stdin(map_output.stdout.unwrap());
    log::info!("{sam_cmd:?}");
    sam_cmd.output()?;
    Ok(())
}

fn eventalign_collapse(
    nanopolish: &Path,
    reads: &Path,
    bam: &Path,
    genome: &Path,
    output: &Path,
) -> Result<()> {
    let cmd = Command::new(nanopolish)
        .arg("eventalign")
        .arg("-r")
        .arg(reads)
        .arg("-b")
        .arg(bam)
        .arg("-g")
        .arg(genome)
        .arg("-t")
        .arg("4")
        .arg("--scale-events")
        .arg("--print-read-names")
        .arg("--samples")
        .stdout(Stdio::piped())
        .spawn()?
        .stdout
        .ok_or_else(|| eyre::eyre!("Could not capture stdout"))?;
    let reader = BufReader::new(cmd);
    let mut collapse = CollapseOptions::try_new(bam, output)?;
    collapse.run(reader)?;
    Ok(())
}

fn train_npsmlr(collapse_file: &Path, single: bool) -> Result<Model> {
    let train_opts = TrainOptions::default().dbscan(true).single(single);
    let reader = File::open(collapse_file)?;
    let model = train_opts.run_model(reader)?;
    Ok(model)
}

fn rank_models(
    rank_output: &Path,
    pos_model: &Model,
    neg_model: &Model,
) -> Result<FnvHashMap<String, f64>> {
    let mut rank_opts = RankOptions::default();
    let ranks = rank_opts.rank(pos_model, neg_model);
    ranks.save_as(rank_output)?;
    Ok(ranks)
}

// Takes a path reads and checks if it is a directory. If its a directory, find all the fastqs
// and concatenate them all into a single file.
fn reads_to_single_reads(reads: &Path, name: &str, output_dir: &Path) -> Result<PathBuf> {
    if reads.is_dir() {
        log::info!("detected directory for neg-reads, cat-ing into a single file");
        let output_filepath = output_dir.join(name);
        let output_file = File::create(&output_filepath)?;
        let glob = format!(
            "{}/*fastq",
            reads.as_os_str().to_str().ok_or(eyre::eyre!(
                "Failed to convert path into str, unicdoe issue?"
            ))?
        );
        Command::new("cat").arg(glob).stdout(output_file);
        Ok(output_filepath)
    } else {
        Ok(reads.to_path_buf())
    }
}

fn main() -> eyre::Result<()> {
    let args = Args::parse();

    let nanopolish = utils::find_binary("nanopolish", &args.nanopolish_path)?;
    let minimap2 = utils::find_binary("minimap2", &args.minimap2_path)?;
    let samtools = utils::find_binary("samtools", &args.samtools_path)?;

    fs::create_dir_all(&args.output_dir)?;

    let log_file = args.output_dir.join("log.txt");
    simple_logging::log_to_file(log_file, LevelFilter::Info)?;

    let neg_reads = reads_to_single_reads(&args.neg_reads, "neg_reads.fastq", &args.output_dir)?;
    let pos_reads = reads_to_single_reads(&args.pos_reads, "pos_reads.fastq", &args.output_dir)?;

    wrap_cmd("nanopolish index for (+) ctrl", || {
        np_index(&nanopolish, &args.pos_fast5s, &pos_reads, &args.pos_summary)
    })?;
    wrap_cmd("nanopolish index for (-) ctrl", || {
        np_index(&nanopolish, &args.neg_fast5s, &neg_reads, &args.neg_summary)
    })?;

    let pos_aln = args.output_dir.join("pos.bam");
    wrap_cmd("align (+) ctrl reads", || {
        aln_reads(&minimap2, &samtools, &args.genome, &pos_reads, &pos_aln)
    })?;
    let neg_aln = args.output_dir.join("neg.bam");
    wrap_cmd("align (-) ctrl reads", || {
        aln_reads(&minimap2, &samtools, &args.genome, &neg_reads, &neg_aln)
    })?;

    let pos_collapse = args.output_dir.join("pos_collapse.arrow");
    wrap_cmd("nanopolish eventalign (+) ctrl | cawlr collapse", || {
        eventalign_collapse(
            &nanopolish,
            &pos_reads,
            &pos_aln,
            &args.genome,
            &pos_collapse,
        )
    })?;

    let neg_collapse = args.output_dir.join("neg_collapse.arrow");
    wrap_cmd("nanopolish eventalign (-) ctrl | cawlr collapse", || {
        eventalign_collapse(
            &nanopolish,
            &neg_reads,
            &neg_aln,
            &args.genome,
            &neg_collapse,
        )
    })?;

    let pos_train = args.output_dir.join("pos_train.pickle");
    let neg_train = args.output_dir.join("neg_train.pickle");

    let pos_model = wrap_cmd_output("Train (+) ctrl", || train_npsmlr(&pos_collapse, false))?;
    pos_model.save_as(pos_train)?;
    let neg_model = wrap_cmd_output("Train (-) ctrl", || train_npsmlr(&neg_collapse, true))?;
    neg_model.save_as(neg_train)?;

    let rank_output = args.output_dir.join("ranks.pickle");
    let ranks = wrap_cmd_output("ranking model kmers", || {
        rank_models(&rank_output, &pos_model, &neg_model)
    })?;

    let score_opts = ScoreOptions::new(pos_model, neg_model, ranks, 10, 10.0, args.motifs.clone());

    let pos_scores_path = args.output_dir.join("pos_scored.arrow");
    wrap_cmd("Scoring (+) ctrl", || {
        let pos_collapse = File::open(&pos_collapse)?;
        let pos_scores = File::create(&pos_scores_path)?;
        score_opts.run(pos_collapse, &pos_scores)
    })?;

    let neg_scores_path = args.output_dir.join("neg_scored.arrow");
    wrap_cmd("Scoring (-) ctrl", || {
        let neg_collapse = File::open(&neg_collapse)?;
        let neg_scores = File::create(&neg_scores_path)?;
        score_opts.run(neg_collapse, neg_scores)
    })?;

    wrap_cmd("(+) model score dist", || {
        let pos_scores = File::open(&pos_scores_path)?;
        let pos_bkde_path = args.output_dir.join("pos_model_scores.pickle");
        let pos_bkde = Options::default().run(pos_scores)?;
        pos_bkde.save_as(pos_bkde_path)
    })?;

    wrap_cmd("(-) model score dist", || {
        let neg_scores = File::open(&neg_scores_path)?;
        let neg_bkde_path = args.output_dir.join("neg_model_scores.pickle");
        let neg_bkde = Options::default().run(neg_scores)?;
        neg_bkde.save_as(neg_bkde_path)
    })?;

    Ok(())
}
