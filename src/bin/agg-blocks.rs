use std::{
    fs::File,
    io::{LineWriter, Write, BufReader, BufRead},
    path::PathBuf,
};

use clap::Parser;
use fnv::{FnvHashMap, FnvHashSet};
use serde::{de::IgnoredAny, Deserialize};
use serde_with::{formats::CommaSeparator, serde_as, StringWithSeparator};

#[derive(Parser)]
struct Args {
    /// Bed file, usually from cawlr sma
    #[clap(short, long)]
    input: PathBuf,

    /// Output tsv containing chromosome, position, frac overlapped
    #[clap(short, long)]
    output: PathBuf,
}

#[serde_as]
#[derive(Deserialize)]
struct Bed {
    chrom: String,
    start: u64,
    stop: u64,
    _extra: IgnoredAny,
    _score: IgnoredAny,
    _strand: IgnoredAny,
    _thick_start: IgnoredAny,
    _thick_end: IgnoredAny,
    _item_rgb: IgnoredAny,
    _bcount: IgnoredAny,
    #[serde_as(as = "StringWithSeparator::<CommaSeparator, u64>")]
    bsizes: Vec<u64>,
    #[serde_as(as = "StringWithSeparator::<CommaSeparator, u64>")]
    bstarts: Vec<u64>,
}

impl Bed {
    fn iter_counts(self) -> impl Iterator<Item = Position> {
        self.bsizes
            .into_iter()
            .zip(self.bstarts.into_iter())
            .map(move |(a, b)| Position::new(self.chrom.clone(), self.start + a + b))
    }

    fn overlaps(self) -> FnvHashSet<Position> {
        self.iter_counts().collect()
    }
}

#[derive(Default)]
struct Count {
    count: u64,
    total: u64,
}

impl Count {
    fn both(&mut self) {
        self.count += 1;
        self.total += 1;
    }

    fn total(&mut self) {
        self.total += 1;
    }

    fn frac(&self) -> f64 {
        (self.count as f64) / (self.total as f64)
    }
}

#[derive(Eq, Hash, PartialEq, Clone)]
struct Position {
    chrom: String,
    pos: u64,
}

impl Position {
    fn new(chrom: String, pos: u64) -> Self {
        Self { chrom, pos }
    }
}

fn main() -> eyre::Result<()> {
    let args = Args::parse();
    let mut input = BufReader::new(File::open(args.input)?);
    // Skip header
    input.read_line(&mut String::new())?;  
    let mut reader = csv::ReaderBuilder::new()
        .has_headers(true)
        .delimiter(b'\t')
        .from_reader(input.into_inner());
    let reader = reader.deserialize::<Bed>();
    let mut counts: FnvHashMap<Position, Count> = FnvHashMap::default();
    for line in reader {
        let line = line?;
        let chrom = line.chrom.clone();
        let start = line.start;
        let stop = line.stop;
        let overlapped = line.overlaps();
        (start..stop).for_each(|pos| {
            let pos = Position::new(chrom.clone(), pos);
            let e = counts.entry(pos.clone()).or_default();
            if overlapped.contains(&pos) {
                e.both();
            } else {
                e.total();
            }
        });
    }

    let mut output = LineWriter::new(File::open(args.output)?);
    for (p, c) in counts.into_iter() {
        writeln!(
            &mut output,
            "{}\t{}\t{}\t{}",
            p.chrom,
            p.pos,
            c.total,
            c.frac()
        )?;
    }
    Ok(())
}
