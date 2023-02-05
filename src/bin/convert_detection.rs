use std::{
    error::Error,
    fs::File,
    io::BufWriter,
    path::{Path, PathBuf},
};

use libcawlr::{
    arrow::{
        arrow_utils::{save, wrap_writer},
        metadata::{Metadata, MetadataExt, Strand},
        scored_read::{Score, ScoredRead},
    },
    plus_strand_map::PlusStrandMap,
};
use clap::Parser;
use indicatif::{ProgressBar, ProgressStyle};
use serde::Deserialize;

#[derive(Parser)]
struct Args {
    /// detection.txt output from NP-SMLR tool
    #[clap(short, long)]
    input: PathBuf,

    /// bam file for adding strand information
    #[clap(short, long)]
    bam: Option<PathBuf>,

    /// Arrow file for use in cawlr sma
    #[clap(short, long)]
    output: PathBuf,
}

#[derive(Deserialize)]
struct DetectionLine {
    chrom: String,
    pos: u64,
    kmer: String,
    read_name: String,
    _pos_log_prob: f64,
    _neg_log_prob: f64,
    score: f64,
}

impl DetectionLine {
    fn read_name(&self) -> &str {
        &self.read_name
    }
}

fn convert_to_read(dlines: &[DetectionLine]) -> ScoredRead {
    let chrom = dlines[0].chrom.clone();
    let read_name = dlines[0].read_name.clone();
    let start = dlines.iter().map(|dline| dline.pos).min().unwrap();
    let end = dlines.iter().map(|dline| dline.pos).max().unwrap();
    let meta = Metadata::new(
        read_name,
        chrom,
        start,
        end - start + 1,
        Strand::unknown(),
        String::new(),
    );
    let scores: Vec<Score> = dlines
        .iter()
        .map(|dline| {
            Score::new(
                dline.pos,
                dline.kmer.clone(),
                false,
                Some(dline.score),
                0.0,
                dline.score,
            )
        })
        .collect();
    ScoredRead::new(meta, scores)
}

pub fn run(input: &Path, bam: &Option<PathBuf>, output: &Path) -> eyre::Result<()> {
    let strand_map = {
        if let Some(bam_file) = bam {
            PlusStrandMap::from_bam_file(bam_file)?
        } else {
            PlusStrandMap::default()
        }
    };

    let style = ProgressStyle::with_template("[{elapsed_precise}] [{binary_bytes_per_sec}] {msg}")?;
    let pb = ProgressBar::new_spinner().with_style(style);

    let reader = File::open(input)?;
    let writer = File::create(output)?;
    let schema = ScoredRead::schema();
    let mut writer = wrap_writer(BufWriter::new(writer), &schema)?;
    let reader = pb.wrap_read(reader);
    let mut builder = csv::ReaderBuilder::new()
        .has_headers(false)
        .delimiter(b'\t')
        .from_reader(reader);
    let mut iter = builder.deserialize::<DetectionLine>().flatten();
    let mut acc = vec![iter.next().unwrap()];
    let mut curr_read = acc[0].read_name().to_owned();

    for dline in iter {
        if dline.read_name() == curr_read {
            acc.push(dline);
        } else {
            let mut read = convert_to_read(&acc);
            if let Some(plus_stranded) = strand_map.get(read.name()) {
                read.metadata.strand = if plus_stranded {
                    Strand::plus()
                } else {
                    Strand::minus()
                };
            }
            save(&mut writer, &[read])?;
            curr_read = dline.read_name().to_owned();
            acc = vec![dline];
        }
    }
    if !acc.is_empty() {
        let read = convert_to_read(&acc);
        save(&mut writer, &[read])?;
    }
    writer.finish()?;
    Ok(())
}

#[allow(dead_code)]
fn main() -> Result<(), Box<dyn Error>> {
    let args = Args::parse();

    run(&args.input, &args.bam, &args.output)?;
    Ok(())
}
