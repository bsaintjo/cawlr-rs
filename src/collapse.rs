use std::{fs::File, io::Read, path::Path};

use anyhow::Result;
use indicatif::{ProgressBar, ProgressBarIter, ProgressFinish, ProgressStyle};
use parquet::arrow::ArrowWriter;
use rstats::Stats;
use serde::Deserialize;
use serde_arrow::{to_record_batch, trace_schema, Schema};
use serde_with::{serde_as, CommaSeparator, StringWithSeparator};

use crate::reads::FlatLReadLData;

fn nprs_to_flatreads(nprs: &[Npr]) -> Result<Vec<FlatLReadLData>> {
    log::debug!("nprs length: {}", nprs.len());
    let read_start = nprs.get(0).ok_or(anyhow::anyhow!("Empty nprs"))?.position;
    log::debug!("Read start {read_start}");
    let mut stop = 0;
    let mut acc = Vec::with_capacity(nprs.len());
    for npr in nprs.iter() {
        let name = npr.read_name.as_bytes();
        let chrom = &npr.contig;
        stop = npr.position;
        let start = read_start as usize;
        let length = 0;
        let seq: &[u8] = &[];
        let position = npr.position;
        let kmer = &npr.model_kmer;
        let mean = npr.samples.amean().expect("No values");
        let time = npr.event_length;
        acc.push(FlatLReadLData::new(
            name, chrom, start, length, seq, position, kmer, mean, time,
        ));
    }
    log::debug!("Read stop {stop}");
    acc.iter_mut().for_each(|fnpr| {
        *fnpr.length_mut() = (stop - read_start) as usize;
    });
    Ok(acc)
}

/// Create spinner that wraps an IO read iterator
fn spin_iter<I: Read>(iter: I) -> ProgressBarIter<I> {
    let pb = ProgressBar::new_spinner();
    let style = ProgressStyle::default_spinner()
        .template("{spinner} [{elapsed_precise}] {binary_bytes_per_sec} {msg}")
        .on_finish(ProgressFinish::AndLeave);
    pb.with_message("Processing eventalign data")
        .with_style(style)
        .wrap_read(iter)
}

pub(crate) struct CollapseOptions {
    input: String,
    writer: ArrowWriter<File>,
    schema: Schema,
    capacity: usize,
}

impl CollapseOptions {
    pub(crate) fn try_new<P>(input: &str, output: P, capacity: usize) -> Result<Self>
    where
        P: AsRef<Path>,
    {
        let example = vec![FlatLReadLData::default()];
        let schema = trace_schema(&example)?;
        let batches = to_record_batch(&example, &schema)?;
        let writer = File::create(output)?;
        let writer = ArrowWriter::try_new(writer, batches.schema(), None)?;
        Ok(CollapseOptions {
            input: input.to_owned(),
            writer,
            schema,
            capacity,
        })
    }

    fn save_flatreads(&mut self, flat_reads: &[FlatLReadLData]) -> Result<()> {
        let batches = to_record_batch(&flat_reads, &self.schema)?;
        self.writer.write(&batches)?;
        Ok(())
    }

    pub(crate) fn close(mut self) -> Result<()> {
        self.writer.close()?;
        Ok(())
    }

    pub(crate) fn run(&mut self) -> Result<()> {
        let file = File::open(&self.input)?;
        let file = spin_iter(file);
        let mut builder = csv::ReaderBuilder::new().delimiter(b'\t').from_reader(file);
        let mut npr_iter = builder.deserialize().peekable();

        let mut acc: Vec<Npr> = Vec::new();
        let mut flats: Vec<FlatLReadLData> = Vec::with_capacity(self.capacity);

        // First get a single line to intialize with, if None we have reached the end of
        // the iterator
        while let Some(line) = npr_iter.next() {
            let mut line: Npr = line?;
            log::debug!("line: {:?}", line);

            // Peek at next lines and advance the iterator only if they belong to the same
            // read
            while let Some(read_line) = npr_iter.next_if(|new_npr: &Result<Npr, _>| {
                let new_read_name = &new_npr.as_ref().expect("Parsing failed").read_name;
                new_read_name == &line.read_name
            }) {
                let mut read_line = read_line?;
                log::debug!("same read: {:?}", read_line);
                // Data from same position, split across two lines
                if read_line.position == line.position {
                    log::debug!("Same position");
                    line.samples.append(&mut read_line.samples);
                    line.event_length += read_line.event_length;
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
            let mut flat_reads = nprs_to_flatreads(&acc)?;
            flats.append(&mut flat_reads);
            if flats.len() >= self.capacity {
                self.save_flatreads(&flats)?;
                flats.clear();
            }

            acc.clear();
        }

        // If reads are left in the buffer, save those
        if !flats.is_empty() {
            self.save_flatreads(&flats)?;
        }
        Ok(())
    }
}

#[serde_as]
#[derive(Clone, Debug, Deserialize)]
struct Npr {
    contig: String,

    position: u64,

    #[serde(skip)]
    _reference_kmer: String,

    read_name: String,

    #[serde(skip)]
    _strand: String,

    #[serde(skip)]
    _event_index: u64,

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

    use crate::{utils::CawlrIO, reads::{LRead, LData}};

    use super::*;

    #[test]
    fn test_collapse() -> Result<()> {
        let temp_dir = TempDir::new()?;
        let filepath = "extra/single_read.eventalign.txt";
        let output = temp_dir.path().join("test");
        let mut collapse = CollapseOptions::try_new(filepath, &output, 2048)?;
        collapse.run()?;
        collapse.close()?;

        let xs: Vec<LRead<LData>> = CawlrIO::load(&output)?;
        assert_eq!(xs.len(), 1);
        assert_eq!(xs[0].data()[0].pos(), 182504);
        let data_len = xs[0].data().len();
        assert_eq!(xs[0].data()[data_len - 1].pos(), 182681);
        Ok(())
    }
}
