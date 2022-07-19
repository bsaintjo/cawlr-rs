use std::{
    fs::File,
    io::{BufWriter, Read, Write},
    path::{Path, PathBuf},
};

use anyhow::{Context, Result};
use arrow2::io::ipc::write::FileWriter;
use bio::alphabets::dna::revcomp;
use indicatif::{ProgressBar, ProgressBarIter, ProgressFinish, ProgressStyle};
use rstats::Stats;
use serde::Deserialize;
use serde_with::{serde_as, CommaSeparator, StringWithSeparator};

use crate::arrow::{self, save, Eventalign, Signal};

fn empty_from_npr(npr: &Npr) -> Eventalign {
    let name = npr.read_name.clone();
    let chrom = npr.contig.clone();
    let start = npr.position;
    let length = 0;
    let seq = String::new();
    Eventalign::empty(name, chrom, start, length, seq)
}

/// Takes a vector of nanpolish records and converts them into a Eventalign.
/// Since strand info isn't immediately available, it is inferred by comparing
/// the reference kmer to the model kmer. If they are the same then it is a plus
/// strand read, if it isn't
/// Differences between model kmer and reference kmer can be due to,
/// - Mismatch
/// - NNNNNN, ie nanopolish failed to choose any 6-mer model
/// - Opposite strand +? above
/// TODO: Strand inference will almost always work but needs to handle edge
/// cases, ie mismatch in low complexity region at end of the read can lead to
/// strand "switching"
/// TODO Strand inference by difference in event_index
fn nprs_to_eventalign(nprs: Vec<Npr>) -> Result<Option<Eventalign>> {
    let mut eventalign = nprs
        .get(0)
        .ok_or(anyhow::anyhow!("Empty nprs"))
        .map(empty_from_npr)?;
    let mut stop = eventalign.start_0b();
    for npr in nprs.into_iter() {
        stop = npr.position;
        let position = npr.position;
        let kmer = npr.model_kmer;
        let ref_kmer = npr.reference_kmer;
        if kmer == ref_kmer {
            *eventalign.strand_mut() = arrow::Strand::plus();
        }
        let rev_ref_kmer = revcomp(ref_kmer.as_bytes());
        if kmer.as_bytes() == rev_ref_kmer {
            *eventalign.strand_mut() = arrow::Strand::minus();
        }
        let mean = npr
            .samples
            .amean()
            .context("No signal sample values to take mean of, malformed input?")?;
        let time = npr.event_length;
        let signal = Signal::new(position, ref_kmer, mean, time);
        eventalign.signal_data_mut().push(signal);
    }
    log::debug!("Read stop {stop}");

    // Handle last edge case with multi-mapped reads, throwing away the read if
    // length calculation leads to overflow
    if let Some(len) = stop.checked_sub(eventalign.start_0b()) {
        *eventalign.length_mut() = len + 1;
    } else {
        return Ok(None);
    }

    // Unable to infer read strand so we remove the read
    if eventalign.strand().is_unknown_strand() {
        return Ok(None);
    }

    if eventalign.strand().is_minus_strand() {
        for signal in eventalign.signal_data_mut().iter_mut() {
            let rev_kmer = revcomp(signal.kmer().as_bytes());
            let rev_kmer = String::from_utf8(rev_kmer)?;
            signal.set_kmer(rev_kmer)
        }
    }
    Ok(Some(eventalign))
}

/// Create spinner that wraps an IO read iterator
fn spin_iter<I: Read>(iter: I, show_progress: bool) -> ProgressBarIter<I> {
    let pb = if show_progress {
        ProgressBar::new_spinner()
    } else {
        ProgressBar::hidden()
    };
    let style = ProgressStyle::default_spinner()
        .template("{spinner} [{elapsed_precise}] {binary_bytes_per_sec} {msg}")
        .on_finish(ProgressFinish::AndLeave);
    pb.with_message("Processing eventalign data")
        .with_style(style)
        .wrap_read(iter)
}

pub struct CollapseOptions<W: Write> {
    input: PathBuf,
    writer: FileWriter<W>,
    capacity: usize,
    progress: bool,
}

impl CollapseOptions<BufWriter<File>> {
    pub fn try_new<P, Q>(input: P, output: Q) -> Result<Self>
    where
        P: AsRef<Path>,
        Q: AsRef<Path>,
    {
        // let writer = File::create(output)?;
        let writer = File::create(output)?;
        let writer = BufWriter::new(writer);
        CollapseOptions::from_writer(input.as_ref().to_owned(), writer)
    }
}

impl<W: Write> CollapseOptions<W> {
    pub fn capacity(&mut self, capacity: usize) -> &mut Self {
        self.capacity = capacity;
        self
    }

    pub fn progress(&mut self, progress: bool) -> &mut Self {
        self.progress = progress;
        self
    }

    pub(crate) fn from_writer(input: PathBuf, writer: W) -> Result<Self> {
        let schema = arrow::Eventalign::schema();
        let writer = arrow::wrap_writer(writer, &schema)?;
        Ok(CollapseOptions {
            input,
            writer,
            capacity: 2048,
            progress: false,
        })
    }

    fn save_eventalign(&mut self, eventaligns: &[Eventalign]) -> Result<()> {
        save(&mut self.writer, eventaligns)
    }

    fn close(mut self) -> Result<()> {
        self.writer.finish()?;
        Ok(())
    }

    pub fn run(mut self) -> Result<()> {
        let file = File::open(&self.input)?;
        let file = spin_iter(file, self.progress);
        let mut builder = csv::ReaderBuilder::new().delimiter(b'\t').from_reader(file);
        let mut npr_iter = builder.deserialize().peekable();

        let mut flats = Vec::with_capacity(self.capacity);

        // First get a single line to intialize with, if None we have reached the end of
        // the iterator
        while let Some(line) = npr_iter.next() {
            let mut acc: Vec<Npr> = Vec::new();
            let mut line: Npr = line?;
            log::debug!("line: {:?}", line);

            // Peek at next lines and advance the iterator only if they belong to the same
            // read
            // Since multiple mapping reads can potentially be next to each other, we need
            // to use the event index to differentiate between the reads. This can
            // potentially cause subtle bugs when the event index of the first multi-mapped
            // read is one minus the event index of the next read.
            // eg. chrom position read_name event_index
            // chr1 100 .. ABC .. 82
            // chr1 10 .. ABC .. 83  <-- New alignment of same read but event_index happens
            // to show the read as being contiguous
            // In order to handle this last
            // case if the length calculation fails in creating the Eventalign,
            // the read will get thrown out.
            while let Some(read_line) = npr_iter.next_if(|new_npr: &Result<Npr, _>| {
                let new_npr = new_npr.as_ref().expect("Parsing failed");
                (new_npr.read_name == line.read_name)
                    && (new_npr.event_index.abs_diff(line.event_index) == 1)
            }) {
                let mut read_line = read_line?;
                log::debug!("same read: {:?}", read_line);
                // Data from same position, split across two lines
                if read_line.position == line.position {
                    log::debug!("Same position");
                    line.samples.append(&mut read_line.samples);
                    line.event_length += read_line.event_length;
                    line.event_index = read_line.event_index
                } else {
                    log::debug!("New position, same read, pushing");
                    // Data from next position in read
                    acc.push(line);
                    // Update line to look for possible repeating lines after it
                    line = read_line;
                }
            }

            // End of the iterator, handle the last value
            if npr_iter.peek().is_none() {
                log::debug!("Iterator finished");
                acc.push(line)
            }
            // Add converted reads to buffer and if buffer is full
            // write to disc and clear buffer
            if let Some(eventalign) = nprs_to_eventalign(acc)? {
                flats.push(eventalign);
            }
            if flats.len() >= self.capacity {
                self.save_eventalign(&flats)?;
                flats.clear();
            }
        }

        // If reads are left in the buffer, save those
        if !flats.is_empty() {
            self.save_eventalign(&flats)?;
        }
        self.close()
    }
}

#[serde_as]
#[derive(Default, Clone, Debug, Deserialize)]
struct Npr {
    contig: String,

    position: u64,

    reference_kmer: String,

    read_name: String,

    #[serde(skip)]
    _strand: String,

    event_index: i64,

    #[serde(skip)]
    _event_level_mean: f64,

    #[serde(skip)]
    _event_stdv: f64,

    event_length: f64,

    model_kmer: String,

    #[serde(skip)]
    _model_mean: f64,

    #[serde(skip)]
    _model_stdv: f64,

    #[serde(skip)]
    _standardized_level: f64,

    #[serde_as(as = "StringWithSeparator::<CommaSeparator, f64>")]
    samples: Vec<f64>,
}

#[cfg(test)]
mod test {

    use assert_fs::TempDir;

    use super::*;
    use crate::arrow::{load_apply, load_iter, MetadataExt};

    #[test]
    fn test_collapse() -> Result<()> {
        let temp_dir = TempDir::new()?;
        let filepath = "extra/single_read.eventalign.txt";
        let output = temp_dir.path().join("test");
        let collapse = CollapseOptions::try_new(filepath, &output)?;
        collapse.run()?;

        let output = File::open(output)?;
        let x = load_iter(output).next().unwrap().unwrap();
        assert_eq!(x.len(), 1);
        let read = &x[0];
        assert_eq!(read.strand(), arrow::Strand::plus());
        assert_eq!(read.chrom(), "chrXIII");
        assert_eq!(read.start_0b(), 182504);
        assert_eq!(read.end_1b_excl(), 182682);

        assert_eq!(read.seq_stop_1b_excl(), 182687);
        assert_eq!(read.seq_length(), 183);

        Ok(())
    }

    #[test]
    fn test_collapse_big() -> Result<()> {
        let temp_dir = TempDir::new()?;
        let filepath = "extra/neg_control.eventalign.txt";
        let output = temp_dir.path().join("test");
        let collapse = CollapseOptions::try_new(filepath, &output)?;
        collapse.run()?;

        let output = File::open(output)?;
        let mut loads = 0;
        let mut acc = Vec::new();
        load_apply(output, |eventaligns: Vec<Eventalign>| {
            loads += 1;
            for eventalign in eventaligns.into_iter() {
                acc.push(eventalign);
            }
            Ok(())
        })?;
        assert_eq!(loads, 1);
        assert_eq!(acc.len(), 98);
        Ok(())
    }
}
