mod cmd;

use std::{
    fs::{self, File},
    path::{Path, },
    process::{Command, Stdio}, ffi::OsStr,
};

use cawlr::{
    collapse::CollapseOptions,
    utils::{self, wrap_cmd, }, region::Region, agg_blocks, sma::SmaOptions, motif::all_bases,
};
pub use cmd::AnalyzeCmd;
use log::LevelFilter;

pub fn parse_name_from_output_dir<P: AsRef<Path>>(path: P) -> eyre::Result<String> {
    let name = path
        .as_ref()
        .file_name()
        .ok_or(eyre::eyre!("Invalid input directory"))?
        .to_str()
        .ok_or(eyre::eyre!("Invalid path name"))?;
    Ok(name.to_string())
}

fn cluster_region_cmd<S: AsRef<OsStr>>(
    region: &Region,
    pct: f64,
    n_clusters: usize,
    name: &str,
    highlights: &[String],
    sma_path: S,
) -> Command {
    let mut cmd = Command::new("cluster_region.py");
    cmd.arg("-p")
        .arg(pct.to_string())
        .arg("-s")
        .arg(region.start().to_string())
        .arg("-e")
        .arg(region.end().to_string())
        .arg("--suptitle")
        .arg(name)
        .arg("-n")
        .arg(n_clusters.to_string())
        .arg("-i")
        .arg(&sma_path);

    if !highlights.is_empty() {
        cmd.arg("--highlight");
        cmd.args(highlights);
    }
    cmd
}

pub fn run(args: AnalyzeCmd) -> eyre::Result<()> {

    if args.overwrite && args.output_dir.exists() {
        fs::remove_dir_all(&args.output_dir)?;
    }
    fs::create_dir_all(&args.output_dir)?;

    let log_file = args.output_dir.join("log.txt");
    simple_logging::log_to_file(log_file, LevelFilter::Info)?;
    log::info!("{args:?}");

    let name = parse_name_from_output_dir(&args.output_dir)?;
    let motifs = args.motifs.clone().ok_or(eyre::eyre!("Need atleast 1 motif"))?;
    let nanopolish = utils::find_binary("nanopolish", &args.nanopolish_path)?;

    let filtered_bam = args.output_dir.join("filtered.bam");
    wrap_cmd("Running samtools", || {
        let samtools = utils::find_binary("samtools", &args.samtools_path)?;
        let mut cmd = Command::new(samtools);
        cmd.arg("view")
            .arg("-hb")
            .arg("--write-index")
            .arg(&args.bam)
            .arg(format!("{}", args.locus))
            .arg("-o")
            .arg(&filtered_bam);
        log::info!("{cmd:?}");
        log::info!("Output file: {}", filtered_bam.display());
        cmd.output()?;
        Ok(())
    })?;

    let eventalign_path = args.output_dir.join("eventalign.tsv");
    wrap_cmd("nanopolish eventalign", || {
        let eventalign = File::create(&eventalign_path)?;
        let eventalign_stdout = Stdio::from(eventalign.try_clone()?);

        let mut cmd = Command::new(&nanopolish);
        cmd.arg("eventalign")
            .arg("--reads")
            .arg(&args.reads)
            .arg("--bam")
            .arg(&filtered_bam)
            .arg("--genome")
            .arg(&args.genome)
            .arg("--scale-events")
            .arg("--print-read-names")
            .arg("--samples")
            .arg("-t")
            .arg(args.n_threads.to_string())
            .stdout(eventalign_stdout);
        log::info!("{cmd:?} >{}", eventalign_path.display());
        cmd.output()?;
        Ok(())
    })?;

    let collapse = args.output_dir.join("collapse.arrow");
    wrap_cmd("cawlr collapse", || {
        let eventalign = File::open(&eventalign_path)?;
        CollapseOptions::try_new(&args.bam, &collapse)?
            .progress(false)
            .run(eventalign)
    })?;

    let scored = args.output_dir.join("score.arrow");
    wrap_cmd("cawlr score", || {
        let mut scoring =
            cawlr::npsmlr::ScoreOptions::load(&args.pos_model, &args.neg_model, &args.ranks)?;
        scoring.motifs(motifs.clone());
        let collapse_file = File::open(&collapse)?;
        let score_file = File::create(&scored)?;
        scoring.run(collapse_file, score_file)
    })?;

    let track_name = format!("{name}.cawlr.sma");
    let sma = args.output_dir.join(format!("{track_name}.bed"));
    wrap_cmd("cawlr sma", || {
        let mut sma_opts =
            SmaOptions::try_new(&args.pos_scores, &args.neg_scores, all_bases(), &sma)?;
        sma_opts.track_name(&track_name);
        sma_opts.run(&scored)
    })?;

    let agg_output = args.output_dir.join(format!("{track_name}.tsv"));
    wrap_cmd("Aggregating blocks", || {
        agg_blocks::run(&sma, Some(&agg_output))
    })?;

    wrap_cmd("Splitting by strand", || {
        let mut cmd = Command::new("split_by_strand.py");
        cmd.arg("-i").arg(&sma);
        log::info!("{cmd:?}");
        cmd.output()?;
        Ok(())
    })?;

    let minus_filepath: &Path = sma.file_stem().unwrap().as_ref();
    let minus_filepath = sma
        .parent()
        .unwrap()
        .join(format!("{}.minus.bed", minus_filepath.display()));

    let plus_filepath: &Path = sma.file_stem().unwrap().as_ref();
    let plus_filepath = sma
        .parent()
        .unwrap()
        .join(format!("{}.plus.bed", plus_filepath.display()));

    wrap_cmd("Clustering all reads", || {
        let mut cmd = cluster_region_cmd(
            &args.locus,
            args.pct,
            args.n_clusters,
            &format!("{name} {} all", args.locus),
            &args.highlights,
            &sma,
        );
        log::info!("{cmd:?}");
        let output = cmd.output()?;
        log::info!("Exit code: {}", output.status);
        Ok(())
    })?;

    wrap_cmd("Clustering (+) reads", || {
        let mut cmd = cluster_region_cmd(
            &args.locus,
            args.pct,
            args.n_clusters,
            &format!("{name} {} plus", args.locus),
            &args.highlights,
            &plus_filepath,
        );
        log::info!("{cmd:?}");
        let output = cmd.output()?;
        log::info!("Exit code: {}", output.status);
        Ok(())
    })?;

    wrap_cmd("Clustering (-) reads", || {
        let mut cmd = cluster_region_cmd(
            &args.locus,
            args.pct,
            args.n_clusters,
            &format!("{name} {} minus", args.locus),
            &args.highlights,
            &minus_filepath,
        );
        log::info!("{cmd:?}");
        let output = cmd.output()?;
        log::info!("Exit code: {}", output.status);
        Ok(())
    })?;

    Ok(())
}
