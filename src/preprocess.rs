use std::{
    collections::{hash_map::Entry, HashMap},
    fs::File,
    io::Read,
    path::Path,
};

use anyhow::Result;
use bio_types::{genome::AbstractInterval, sequence::SequenceRead};
use indicatif::{ProgressBar, ProgressFinish, ProgressStyle};
use rstats::Stats;
use rust_htslib::bam::{ext::BamRecordExtensions, IndexedReader, Read as BamRead, Record};
use serde::{Deserialize, Serialize};
use serde_with::{serde_as, CommaSeparator, StringWithSeparator};

use crate::reads::{LData, PreprocessRead};

pub(crate) type ReadToPreprocessRead = HashMap<Vec<u8>, PreprocessRead>;
pub(crate) type ReadPosToSamples = HashMap<(Vec<u8>, u64), Data>;

fn process_bam<P>(filename: P) -> Result<ReadToPreprocessRead>
where
    P: AsRef<Path>,
{
    let mut acc = HashMap::new();
    let mut bam = IndexedReader::from_path(filename)?;
    let mut record = Record::new();
    // TODO map only unique reads by filtering on mapq or flags
    // TODO should i grab the sequence and store it for later?
    while let Some(result) = bam.read(&mut record) {
        if result.is_ok() {
            let name = record.name().to_owned();
            let key_name = name.clone();
            let chrom = record.contig().to_owned();
            let start = record.reference_start() as usize;
            let len = record.seq_len() as usize;
            let data = Vec::new();
            let read = PreprocessRead::new(name, chrom, start, len, data);
            if let Some(ld) = acc.insert(key_name, read) {
                log::warn!("Duplicate read of encountered in bam with data {:?}", ld);
            }
        }
    }
    Ok(acc)
}

#[serde_as]
#[derive(Clone, Debug, Serialize, Deserialize)]
struct Npr {
    contig: String,

    position: u64,

    #[serde(skip)]
    _reference_kmer: String,

    read_name: Vec<u8>,

    #[serde(skip)]
    _strand: String,

    #[serde(skip)]
    _event_index: u32,

    #[serde(skip)]
    _event_level_mean: f32,

    #[serde(skip)]
    _event_stdv: f32,

    event_length: f64,

    model_kmer: String,

    #[serde(skip)]
    _model_mean: f32,

    #[serde(skip)]
    _model_stdv: f32,

    #[serde(skip)]
    _standardized_level: f32,

    #[serde_as(as = "StringWithSeparator::<CommaSeparator, f64>")]
    samples: Vec<f64>,
}

pub(crate) struct Process {
    chrom: Option<String>,
    start: Option<u64>,
    stop: Option<u64>,
}

// TODO: Store ProgressBar in process and add indicatif support across
/// Reads a nanopolish output file, collapses samples across multiple lines,
/// finds the mean for all the samples, and produces an Output for
/// each collapsed line.
impl Process {
    pub(crate) fn new() -> Self {
        Self {
            chrom: None,
            start: None,
            stop: None,
        }
    }

    pub(crate) fn chrom(mut self, chrom: Option<String>) -> Self {
        self.chrom = chrom;
        self
    }

    pub(crate) fn start(mut self, start: Option<u64>) -> Self {
        self.start = start;
        self
    }

    pub(crate) fn stop(mut self, stop: Option<u64>) -> Self {
        self.stop = stop;
        self
    }

    pub(crate) fn run<P, R>(&self, filename: P, bam_filename: R) -> Result<Vec<PreprocessRead>>
    where
        P: AsRef<Path>,
        R: AsRef<Path>,
    {
        let rp_to_samples = self.with_file(filename)?;
        let mut bam_to_pr = process_bam(bam_filename)?;
        for (key, datum) in rp_to_samples.into_iter() {
            if let Some(pr) = bam_to_pr.get_mut(&key.0) {
                let mean = datum.samples.amean()?;
                let data = LData::new(datum.pos, datum.kmer, mean, datum.time);
                pr.get_mut_data().push(data);
            }
        }
        Ok(bam_to_pr
            .into_values()
            .filter(|lr| lr.is_empty())
            .collect())
    }

    pub(crate) fn with_file<P>(&self, filename: P) -> Result<ReadPosToSamples>
    where
        P: AsRef<Path>,
    {
        let file = File::open(filename)?;
        let style = ProgressStyle::default_spinner()
            .template("{spinner} [{elapsed_precise}] {binary_bytes_per_sec} {msg}")
            .on_finish(ProgressFinish::AndLeave);
        let file = ProgressBar::new_spinner()
            .with_message("Processing data")
            .with_style(style)
            .wrap_read(file);
        Ok(self.with_reader(file))
    }

    pub(crate) fn with_reader<R>(&self, input: R) -> ReadPosToSamples
    where
        R: Read,
    {
        let mut acc: ReadPosToSamples = HashMap::new();
        let mut reader = csv::ReaderBuilder::new()
            .delimiter(b'\t')
            .from_reader(input);
        for line in reader.deserialize().flatten() {
            log::debug!("{:?}", line);
            let line: Npr = line;

            // Filter records that are not on chromosome
            if let Some(ref chrom) = self.chrom {
                if chrom == &line.contig {
                    continue;
                }
            }
            if let Some(start) = self.start {
                if line.position < start {
                    continue;
                }
            }
            if let Some(stop) = self.stop {
                if line.position > stop {
                    continue;
                }
            }

            let key = (line.read_name.clone(), line.position);
            match acc.entry(key) {
                Entry::Occupied(mut oe) => {
                    log::debug!("Occupied");
                    oe.get_mut().update_from_npr(line);
                }
                Entry::Vacant(ve) => {
                    log::debug!("Vacant");
                    let d = Data::new_from_npr(line);
                    ve.insert(d);
                }
            }
        }
        acc
    }
}

#[derive(Debug)]
pub(crate) struct Data {
    pos: u64,
    kmer: String,
    samples: Vec<f64>,
    time: f64,
}

impl Data {
    fn new(pos: u64, kmer: String, samples: Vec<f64>, time: f64) -> Self {
        Self {
            pos,
            kmer,
            samples,
            time,
        }
    }

    fn new_from_npr(npr: Npr) -> Self {
        let samples = npr.samples;
        let pos = npr.position;
        let kmer = npr.model_kmer;
        let time = npr.event_length;
        Self::new(pos, kmer, samples, time)
    }

    fn update_from_npr(&mut self, mut npr: Npr) {
        self.samples.append(&mut npr.samples);
        self.time += npr.event_length;
    }
}

#[cfg(test)]
mod tests {
    use std::io::Cursor;

    use super::*;

    #[test]
    fn test_single_line() {
        let line = b"chrI\t1\tCACACC\t6246746d-97e5-4ace-9362-60d2faffdb93\tt\t529\t96.94\t1.722\t0.00150\tCACACC\t96.60\t1.74\t0.16\t95.5522,99.8873,94.8586,96.2459,97.8065,97.2863";
        let cursor = Cursor::new(line);
        let p = Process::new().with_reader(cursor);
        eprintln!("{:?}", p);
    }
}
