use std::{
    collections::{hash_map::Entry, HashMap},
    fs::File,
    io::Read,
    path::Path,
};

use anyhow::Result;
use indicatif::{ProgressBar, ProgressFinish, ProgressStyle};
use rstats::Stats;
use serde::{Deserialize, Serialize};
use serde_with::{serde_as, CommaSeparator, StringWithSeparator};

#[serde_as]
#[derive(Clone, Debug, Serialize, Deserialize)]
struct Npr {
    contig: String,

    position: u64,

    #[serde(skip)]
    _reference_kmer: String,

    read_name: String,

    #[serde(skip)]
    _strand: String,

    #[serde(skip)]
    _event_index: u32,

    #[serde(skip)]
    _event_level_mean: f32,

    #[serde(skip)]
    _event_stdv: f32,

    event_length: f32,

    model_kmer: String,

    #[serde(skip)]
    _model_mean: f32,

    #[serde(skip)]
    _model_stdv: f32,

    #[serde(skip)]
    _standardized_level: f32,

    #[serde_as(as = "StringWithSeparator::<CommaSeparator, f32>")]
    samples: Vec<f32>,
}

pub(crate) struct Process {
    data: HashMap<(String, u64), Data>,
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
            data: HashMap::new(),
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

    pub(crate) fn with_file<P>(self, filename: P) -> Result<Vec<Output>>
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

    pub(crate) fn with_reader<R>(mut self, input: R) -> Vec<Output>
    where
        R: Read,
    {
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
            match self.data.entry(key) {
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
        self.data
            .into_values()
            .flat_map(|d| d.into_output())
            .collect()
    }
}

#[derive(Debug, Serialize, Deserialize)]
pub(crate) struct Output {
    mean: f64,
    time: f32,
    #[serde(flatten)]
    npr: NanopolishRecord,
}

impl Output {
    pub(crate) fn new(mean: f64, time: f32, npr: NanopolishRecord) -> Self {
        Self { mean, time, npr }
    }

    /// Get a reference to the output's npr.
    pub(crate) fn contig(&self) -> &str {
        &self.npr.contig
    }

    pub(crate) fn pos(&self) -> u64 {
        self.npr.position
    }

    pub(crate) fn mean(&self) -> f64 {
        self.mean
    }

    pub(crate) fn metadata(self) -> NanopolishRecord {
        self.npr
    }
}

#[derive(Serialize, Deserialize)]
pub(crate) struct Data {
    samples: Vec<f32>,
    time: f32,
    #[serde(flatten)]
    npr: NanopolishRecord,
}

impl Data {
    fn new(samples: Vec<f32>, time: f32, npr: NanopolishRecord) -> Self {
        Self { samples, time, npr }
    }

    fn new_from_npr(npr: Npr) -> Self {
        let samples = npr.samples;
        let contig = npr.contig;
        let position = npr.position;
        let read_name = npr.read_name;
        let kmer = npr.model_kmer;
        let time = npr.event_length;
        let metadata = NanopolishRecord::new(contig, position, read_name, kmer);
        Self::new(samples, time, metadata)
    }

    fn update_from_npr(&mut self, mut npr: Npr) {
        self.samples.append(&mut npr.samples);
        self.time += npr.event_length;
    }

    // TODO: serde_arrow doesn't support lists or sequences, so we need to calculate
    // final statistics and use that in the output
    fn into_output(self) -> Result<Output> {
        let mean = self.samples.amean()?;
        Ok(Output::new(mean, self.time, self.npr))
    }
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub(crate) struct NanopolishRecord {
    contig: String,
    position: u64,
    read_name: String,
    kmer: String,
}

impl NanopolishRecord {
    fn new(contig: String, position: u64, read_name: String, kmer: String) -> Self {
        NanopolishRecord {
            contig,
            position,
            read_name,
            kmer,
        }
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
