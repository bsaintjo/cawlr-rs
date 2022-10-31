use std::{
    fs::File,
    io::{BufRead, BufReader},
    path::{Path, PathBuf},
};

use csv::StringRecord;
use fnv::{FnvHashMap, FnvHashSet};
use serde::{de::IgnoredAny, Deserialize};
use serde_with::{formats::CommaSeparator, serde_as, StringWithSeparator};

use crate::utils::stdout_or_file;

#[derive(Eq, Hash, PartialEq, Clone)]
struct Position {
    chrom: String,
    pos: u64,
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

impl Position {
    fn new(chrom: String, pos: u64) -> Self {
        Self { chrom, pos }
    }
}

#[serde_as]
#[derive(Deserialize)]
pub struct Bed {
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

    pub fn bstarts(&self) -> &[u64] {
        &self.bstarts
    }

    pub fn bsizes(&self) -> &[u64] {
        &self.bsizes
    }
}

pub fn run(input: &Path, output: Option<&PathBuf>) -> eyre::Result<()> {
    let input = BufReader::new(File::open(input)?);
    // Skip header

    let mut counts: FnvHashMap<Position, Count> = FnvHashMap::default();
    for rec in input.lines().skip(1) {
        let rec = rec?;
        let line: Vec<&str> = rec.split('\t').collect();
        let line = StringRecord::from(line);
        let line = line.deserialize::<Bed>(None)?;
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

    let mut output = stdout_or_file(output)?;
    for (p, c) in counts.into_iter() {
        writeln!(
            &mut output,
            "{}\t{}\t{}\t{}\t{}",
            p.chrom,
            p.pos,
            c.count,
            c.total,
            c.frac()
        )?;
    }
    Ok(())
}
