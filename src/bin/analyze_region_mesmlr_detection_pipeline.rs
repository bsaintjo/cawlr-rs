use std::{
    fs::{self, File},
    io::{BufRead, BufReader, BufWriter, Write},
    path::PathBuf,
    process::Command,
};

mod convert_detection;

use libcawlr::utils::parse_name_from_output_dir;
use libcawlr::{agg_blocks, motif::all_bases, region::Region, sma::SmaOptions, utils::wrap_cmd};
use clap::Parser;
use log::LevelFilter;

#[derive(Debug, Parser)]
struct Args {
    /// Detection file from NP-SMLR
    #[clap(short, long)]
    detection: PathBuf,

    #[clap(short, long)]
    bam: Option<PathBuf>,

    /// Locus to filter on, in form of {chrom}:{start}-{end}
    #[clap(short, long)]
    locus: Region,

    /// Path to positive control scores
    #[clap(short, long)]
    pos_scores: PathBuf,

    /// Path to negative control scores
    #[clap(short, long)]
    neg_scores: PathBuf,

    /// Path to output directory
    #[clap(short, long)]
    output_dir: PathBuf,

    // #[clap(long)]
    // name: String,
    /// Number of clusters to use for clustering script
    #[clap(long, default_value_t = 3)]
    n_clusters: usize,

    /// Percent of read that should overlap region to be clustered
    #[clap(long)]
    pct: f64,

    /// Regions to highlight during clustering
    #[clap(long)]
    highlights: Vec<String>,

    #[clap(long, default_value_t = false)]
    overwrite: bool,
}

fn convert_chrom(s: &str) -> eyre::Result<&str> {
    let conv = match s {
        "ref|NC_001133|" => "chrI",
        "ref|NC_001134|" => "chrII",
        "ref|NC_001135|" => "chrIII",
        "ref|NC_001136|" => "chrIV",
        "ref|NC_001137|" => "chrV",
        "ref|NC_001138|" => "chrVI",
        "ref|NC_001139|" => "chrVII",
        "ref|NC_001140|" => "chrVIII",
        "ref|NC_001141|" => "chrIX",
        "ref|NC_001142|" => "chrX",
        "ref|NC_001143|" => "chrXI",
        "ref|NC_001144|" => "chrXII",
        "ref|NC_001145|" => "chrXIII",
        "ref|NC_001146|" => "chrXIV",
        "ref|NC_001147|" => "chrXV",
        "ref|NC_001148|" => "chrXVI",
        "ref|NC_001224|" => "chrM",
        _ => return Err(eyre::eyre!("Not valid convertable chrom: {s}")),
    };
    Ok(conv)
}

fn main() -> eyre::Result<()> {
    let args = Args::parse();

    if args.overwrite && args.output_dir.exists() {
        fs::remove_dir_all(&args.output_dir)?;
    }
    fs::create_dir_all(&args.output_dir)?;

    let log_file = args.output_dir.join("log.txt");
    simple_logging::log_to_file(log_file, LevelFilter::Info)?;
    log::info!("{args:?}");

    let name = parse_name_from_output_dir(&args.output_dir)?;

    let filtered_output_path = args.output_dir.join("filtered_detection.txt");
    let mut writer = BufWriter::new(File::create(&filtered_output_path)?);

    wrap_cmd("Filtering detection.txt", || {
        let detection = BufReader::new(File::open(&args.detection)?);
        for line in detection.lines() {
            let line = line?;
            let mut splitted = line.split('\t').collect::<Vec<_>>();
            let chrom = convert_chrom(splitted[0])?;
            let pos = splitted[1].parse::<u64>()?;
            if (chrom == args.locus.chrom())
                && (args.locus.start() <= pos)
                && (pos <= args.locus.end())
            {
                splitted[0] = chrom;
                let result = splitted.join("\t");
                writeln!(&mut writer, "{}", result)?;
            }
        }
        writer.flush()?;
        Ok(())
    })?;

    let converted_output = args.output_dir.join("converted.arrow");
    wrap_cmd("Convert to score.arrow", || {
        convert_detection::run(&filtered_output_path, &args.bam, &converted_output)
    })?;

    let track_name = format!("{}.sma", name);
    let sma = args.output_dir.join(format!("{}.bed", track_name));
    wrap_cmd("cawlr sma", || {
        let mut sma_opts =
            SmaOptions::try_new(&args.pos_scores, &args.neg_scores, all_bases(), &sma)?;
        sma_opts.track_name(&track_name);
        sma_opts.run(&converted_output)
    })?;

    let agg_output = args.output_dir.join(format!("{}.tsv", track_name));
    wrap_cmd("Aggregating", || agg_blocks::run(&sma, Some(&agg_output)))?;

    wrap_cmd("Clustering reads", || {
        let mut cmd = Command::new("cluster_region.py");
        cmd.arg("-p")
            .arg(args.pct.to_string())
            .arg("-s")
            .arg(args.locus.start().to_string())
            .arg("-e")
            .arg(args.locus.end().to_string())
            .arg("--suptitle")
            .arg(format!("{name} {}", args.locus))
            .arg("-n")
            .arg(args.n_clusters.to_string())
            .arg("-i")
            .arg(&sma);
        if !args.highlights.is_empty() {
            cmd.arg("--highlight");
            cmd.args(&args.highlights);
        }
        log::info!("{cmd:?}");
        let output = cmd.output()?;
        log::info!("{}", output.status);
        Ok(())
    })?;
    Ok(())
}
