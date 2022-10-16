use std::{
    fs::{self, File},
    path::PathBuf,
    process::{Command, Stdio},
    time::Duration,
};

use cawlr::{
    collapse::CollapseOptions, filter::Region, motif::{Motif, all_bases}, score::ScoreOptions, sma::SmaOptions,
    utils,
};
use clap::Parser;
use indicatif::{ProgressBar, ProgressFinish};

#[derive(Parser)]
struct Args {
    /// Used to name output bed file
    #[clap(short, long)]
    name: String,

    /// Region of interested {chromosome}:{start}-{stop}
    #[clap(short, long)]
    locus: Region,

    /// Where to output results
    #[clap(short, long)]
    output_dir: PathBuf,

    /// Path to bam file to filter on the locus
    #[clap(short, long)]
    bam: PathBuf,

    /// Path to full fastq, doesn't need to be filtered
    #[clap(short, long)]
    fastq: PathBuf,

    /// Path to genome
    #[clap(short, long)]
    genome: PathBuf,

    /// Path to postive control model, from cawlr train
    #[clap(long)]
    pos_model: PathBuf,

    /// Path to postive control scores, from cawlr model-scores
    #[clap(long)]
    pos_scores: PathBuf,

    /// Path to negative control model, from cawlr train
    #[clap(long)]
    neg_model: PathBuf,

    /// Path to negative control scores, from cawlr model-scores
    #[clap(long)]
    neg_scores: PathBuf,

    /// Path to ranks file, from cawlr ranks
    #[clap(short, long)]
    ranks: PathBuf,

    /// Motifs to analyze, formatted "2:GC" if second base C is modified
    /// Can have more than one
    #[clap(short, long)]
    motifs: Option<Vec<Motif>>,

    /// Path to nanopolish binary, if not specified will look in $PATH
    #[clap(long)]
    nanopolish_path: Option<PathBuf>,

    /// Path to samtools binary, if not specified will look in $PATH
    #[clap(long)]
    samtools_path: Option<PathBuf>,
}

fn pb(msg: &'static str) -> ProgressBar {
    let p = ProgressBar::new_spinner()
        .with_message(msg)
        .with_finish(ProgressFinish::AndLeave);
    p.enable_steady_tick(Duration::from_secs(1));
    p
}

fn wrap_cmd<F>(msg: &'static str, f: F) -> eyre::Result<()>
where
    F: Fn() -> eyre::Result<()>,
{
    let p = pb(msg);
    f()?;
    p.finish();
    Ok(())
}

fn main() -> eyre::Result<()> {
    let args = Args::parse();
    let motifs = args.motifs.ok_or(eyre::eyre!("Need atleast 1 motif"))?;
    let nanopolish = utils::find_binary("nanopolish", &args.nanopolish_path)?;

    fs::create_dir_all(&args.output_dir)?;
    let filtered_bam = args.output_dir.join("filtered.bam");

    wrap_cmd("Running samtools", || {
        let samtools = utils::find_binary("samtools", &args.samtools_path)?;
        Command::new(samtools)
            .arg("view")
            .arg("-hb")
            .arg("--write-index")
            .arg(&args.bam)
            .arg(format!("{}", args.locus))
            .arg("-o")
            .arg(&filtered_bam)
            .output()?;
        Ok(())
    })?;

    let p = pb("Running nanopolish eventalign");
    let eventalign_path = args.output_dir.join("eventalign.tsv");
    let eventalign = File::create(&eventalign_path)?;
    let eventalign_stdout = Stdio::from(eventalign.try_clone()?);

    Command::new(nanopolish)
        .arg("eventalign")
        .arg("--reads")
        .arg(args.fastq)
        .arg("--bam")
        .arg(&filtered_bam)
        .arg("--genome")
        .arg(&args.genome)
        .arg("--scale-events")
        .arg("--print-read-names")
        .arg("--samples")
        .args(&["-t", "4"])
        .stdout(eventalign_stdout)
        .output()?;
    p.finish();

    let p = pb("Running cawlr collapse");
    let collapse = args.output_dir.join("collapse.arrow");
    let eventalign = File::open(&eventalign_path)?;
    CollapseOptions::try_new(args.bam, &collapse)?.run(eventalign)?;
    p.finish();

    let p = pb("Running cawlr score");
    let scored = args.output_dir.join("score.arrow");
    let mut scoring = ScoreOptions::try_new(
        args.pos_model,
        args.neg_model,
        args.genome,
        args.ranks,
        scored.clone(),
    )?;
    scoring.motifs(motifs);
    scoring.run(collapse)?;
    p.finish();

    let p = pb("Running cawlr sma");
    let sma = args.output_dir.join(format!("{}.sma.bed", args.name));
    SmaOptions::try_new(args.pos_scores, args.neg_scores, all_bases(), sma)?.run(scored)?;
    p.finish();

    Ok(())
}
