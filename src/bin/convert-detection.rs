use std::{error::Error, fs::File, io::BufWriter, path::PathBuf};

use cawlr::{save, wrap_writer, Metadata, Score, ScoredRead, Strand};
use clap::Parser;
use serde::Deserialize;

#[derive(Parser)]
struct Args {
    #[clap(short, long)]
    input: PathBuf,

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
        Strand::Unknown,
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

fn main() -> Result<(), Box<dyn Error>> {
    let args = Args::parse();

    let reader = File::open(&args.input)?;
    let writer = File::create(&args.output)?;
    let schema = ScoredRead::schema();
    let mut writer = wrap_writer(BufWriter::new(writer), &schema)?;
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
            let read = convert_to_read(&acc);
            save(&mut writer, &[read])?;
            curr_read = dline.read_name().to_owned();
            acc = vec![dline];
        }
    }
    writer.finish()?;
    Ok(())
}
