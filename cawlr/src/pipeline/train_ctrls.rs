use std::{
    fs::{self, File},
    io::{BufRead, BufReader, BufWriter, Write},
    path::{Path, PathBuf},
    process::{Command, Stdio},
};

use clap::Parser;
use eyre::{Context, Result};
use fnv::FnvHashMap;
use libcawlr::{
    collapse::CollapseOptions,
    motif::Motif,
    npsmlr::{train::TrainOptions, ScoreOptions},
    rank::RankOptions,
    score_model::Options,
    train::Model,
    utils::{self, check_if_failed, wrap_cmd, wrap_cmd_output, CawlrIO},
};
use log::LevelFilter;

use crate::file::ValidPathBuf;

#[derive(Parser, Debug)]
pub struct TrainCtrlPipelineCmd {
    /// Path to genome fasta file
    #[clap(short, long)]
    genome: ValidPathBuf,

    /// Directory containing fast5s for positive control
    #[clap(long)]
    pos_fast5: PathBuf,

    /// Path to single fasta/q file or directory of fasta/q of reads from the
    /// positive control
    #[clap(long)]
    pos_reads: PathBuf,

    /// Optional path to sequencing_summary.txt file for positive control,
    /// speeds up nanopolish indexing
    #[clap(long)]
    pos_summary: Option<PathBuf>,

    /// Directory containing fast5s for negative control
    #[clap(long)]
    neg_fast5: PathBuf,

    /// Path to single fasta/q file or directory of fasta/q of reads from the
    /// negative control
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

    /// Motifs of modification to filter on, separated by commas, format is
    /// "{position}:{motif}" ie for GpC and CpG motif , motif is "2:GC,1:CG"
    #[clap(short, long, required=true, num_args=1.., value_delimiter=',')]
    motifs: Vec<Motif>,
}

fn np_index(
    nanopolish: &Path,
    fast5s: &Path,
    reads: &Path,
    summary: &Option<PathBuf>,
    log_file: File,
) -> Result<()> {
    let mut cmd = Command::new(nanopolish);
    cmd.arg("index").arg("-d").arg(fast5s);
    if let Some(summary) = summary {
        cmd.arg("-s").arg(summary);
    }
    cmd.arg(reads);
    cmd.stderr(log_file);
    log::info!("{cmd:?}");
    let output = cmd.output()?;
    check_if_failed(output).wrap_err("nanopolish index failed")
}

fn aln_reads(
    minimap2: &Path,
    samtools: &Path,
    genome: &ValidPathBuf,
    reads: &Path,
    output: &Path,
    output_dir: &Path,
    log_file: File,
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
        .stdout(Stdio::piped())
        .stderr(log_file.try_clone()?);
    log::info!("{map_cmd:?}");
    let map_output = map_cmd.spawn()?;

    let mut sam_cmd = Command::new(samtools);
    sam_cmd
        .arg("sort")
        .arg("--write-index")
        .arg("-T")
        .arg(output_dir)
        .arg("-o")
        .arg(output)
        .stderr(log_file)
        .stdin(map_output.stdout.unwrap());
    log::info!("{sam_cmd:?}");
    let output = sam_cmd.output()?;
    check_if_failed(output).wrap_err("minimap2 | samtools failed")
}

fn eventalign_collapse(
    nanopolish: &Path,
    reads: &Path,
    bam: &Path,
    genome: &ValidPathBuf,
    output: &Path,
    log_file: File,
) -> Result<()> {
    let mut cmd = Command::new(nanopolish);
    cmd.arg("eventalign")
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
        .arg("--samples");
    log::info!("nanopolish cmd: {cmd:?}");
    let mut cmd = cmd.stdout(Stdio::piped()).stderr(log_file).spawn()?;
    let stdout = cmd
        .stdout
        .take()
        .ok_or_else(|| eyre::eyre!("Could not capture stdout"))?;
    let reader = BufReader::new(stdout);
    let mut collapse = CollapseOptions::try_new(bam, output)?;
    collapse.run(reader)?;
    Ok(())
}

fn train_npsmlr(collapse_file: &Path, db_file: &Path, single: bool, motifs: &[Motif]) -> Result<Model> {
    let train_opts = TrainOptions::default()
        .dbscan(true)
        .single(single)
        .db_path(Some(db_file.to_path_buf()))
        .motifs(motifs.to_vec());
    let reader = BufReader::new(File::open(collapse_file)?);
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

// Takes a path reads and checks if it is a directory. If its a directory, find
// all the fastqs and concatenate them all into a single file.
fn reads_to_single_reads(reads: &Path, name: &str, output_dir: &Path) -> Result<PathBuf> {
    if reads.is_dir() {
        log::info!("Detected directory, concatenating into a single fastq file.");
        let output_filepath = output_dir.join(name);
        let mut output_file = BufWriter::new(File::create(&output_filepath)?);
        let fastq_matcher = format!(
            "{}/**/*fastq",
            reads.as_os_str().to_str().ok_or(eyre::eyre!(
                "Failed to convert path into str, unicode issue?"
            ))?
        );
        let mut n_fastq_files = 0;
        for fastq in glob::glob(&fastq_matcher)? {
            let fastq = fastq?;
            n_fastq_files += 1;
            log::info!("Found fastq: {}", fastq.display());
            let mut fastq_file = BufReader::new(File::open(fastq)?);
            loop {
                let buf_len = {
                    let buf = fastq_file.fill_buf()?;
                    if buf.is_empty() {
                        break;
                    }
                    output_file.write_all(buf)?;
                    buf.len()
                };
                fastq_file.consume(buf_len);
            }
        }
        output_file.flush()?;

        if n_fastq_files == 0 {
            return Err(eyre::eyre!(
                "No fastq files processed, check if directory contained files ending with .fastq"
            ));
        } else {
            log::info!("Processed {n_fastq_files} fastq files");
        }

        Ok(output_filepath)
    } else {
        Ok(reads.to_path_buf())
    }
}

pub fn run(args: TrainCtrlPipelineCmd) -> eyre::Result<()> {
    log::info!("{args:?}");
    let nanopolish = utils::find_binary("nanopolish", &args.nanopolish_path)?;
    let minimap2 = utils::find_binary("minimap2", &args.minimap2_path)?;
    let samtools = utils::find_binary("samtools", &args.samtools_path)?;

    fs::create_dir_all(&args.output_dir)?;

    let log_file_path = args.output_dir.join("log.txt");
    let log_file = File::create(log_file_path)?;
    simple_logging::log_to(log_file.try_clone()?, LevelFilter::Info);

    let neg_reads = reads_to_single_reads(&args.neg_reads, "neg_reads.fastq", &args.output_dir)?;
    let pos_reads = reads_to_single_reads(&args.pos_reads, "pos_reads.fastq", &args.output_dir)?;

    wrap_cmd("nanopolish index for (+) ctrl", || {
        np_index(
            &nanopolish,
            &args.pos_fast5,
            &pos_reads,
            &args.pos_summary,
            log_file.try_clone()?,
        )
    })?;
    wrap_cmd("nanopolish index for (-) ctrl", || {
        np_index(
            &nanopolish,
            &args.neg_fast5,
            &neg_reads,
            &args.neg_summary,
            log_file.try_clone()?,
        )
    })?;

    let pos_aln = args.output_dir.join("pos.bam");
    wrap_cmd("align (+) ctrl reads", || {
        aln_reads(
            &minimap2,
            &samtools,
            &args.genome,
            &pos_reads,
            &pos_aln,
            &args.output_dir,
            log_file.try_clone()?,
        )
    })?;
    let neg_aln = args.output_dir.join("neg.bam");
    wrap_cmd("align (-) ctrl reads", || {
        aln_reads(
            &minimap2,
            &samtools,
            &args.genome,
            &neg_reads,
            &neg_aln,
            &args.output_dir,
            log_file.try_clone()?,
        )
    })?;

    let pos_collapse = args.output_dir.join("pos_collapse.arrow");
    wrap_cmd("nanopolish eventalign (+) ctrl | cawlr collapse", || {
        eventalign_collapse(
            &nanopolish,
            &pos_reads,
            &pos_aln,
            &args.genome,
            &pos_collapse,
            log_file.try_clone()?,
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
            log_file.try_clone()?,
        )
    })?;

    let pos_train = args.output_dir.join("pos_train.pickle");
    let neg_train = args.output_dir.join("neg_train.pickle");

    let pos_db_file = args.output_dir.join("pos.db.sqlite3");
    let neg_db_file = args.output_dir.join("neg.db.sqlite3");

    let pos_model = wrap_cmd_output("Train (+) ctrl", || {
        log::info!("Starting  + training");
        train_npsmlr(&pos_collapse, &pos_db_file, false, &args.motifs)
    })?;
    pos_model.save_as(pos_train)?;
    let neg_model = wrap_cmd_output("Train (-) ctrl", || {
        log::info!("Starting - training");
        train_npsmlr(&neg_collapse, &neg_db_file, true, &args.motifs)
    })?;
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
        score_opts.run(pos_collapse, &pos_scores)?;
        log::info!("Finished scoring positive control");
        Ok(())
    })?;

    let neg_scores_path = args.output_dir.join("neg_scored.arrow");
    wrap_cmd("Scoring (-) ctrl", || {
        let neg_collapse = File::open(&neg_collapse)?;
        let neg_scores = File::create(&neg_scores_path)?;
        score_opts.run(neg_collapse, neg_scores)?;
        log::info!("Finished scoring positive control");
        Ok(())
    })?;

    wrap_cmd("(+) model score dist", || {
        let pos_scores = File::open(&pos_scores_path)?;
        let pos_bkde_path = args.output_dir.join("pos_model_scores.pickle");
        let pos_bkde = Options::default().run(pos_scores)?;
        pos_bkde.save_as(pos_bkde_path)?;
        log::info!("Completed BKDE for (+) control");
        Ok(())
    })?;

    wrap_cmd("(-) model score dist", || {
        let neg_scores = File::open(&neg_scores_path)?;
        let neg_bkde_path = args.output_dir.join("neg_model_scores.pickle");
        let neg_bkde = Options::default().run(neg_scores)?;
        neg_bkde.save_as(neg_bkde_path)?;
        log::info!("Completed BKDE for (-) control");
        Ok(())
    })?;

    Ok(())
}
