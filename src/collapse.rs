use std::{fs::File, path::Path};

use anyhow::Result;
use parquet::arrow::ArrowWriter;
use rstats::Stats;
use serde::Deserialize;
use serde_arrow::{to_record_batch, trace_schema, Schema};
use serde_with::{serde_as, CommaSeparator, StringWithSeparator};

use crate::reads::FlatLReadLData;

pub(crate) struct CollapseOptions {
    writer: ArrowWriter<File>,
    schema: Schema,
}

impl CollapseOptions {
    pub(crate) fn try_new<P>(output: P) -> Result<Self>
    where
        P: AsRef<Path>,
    {
        let example = vec![FlatLReadLData::default()];
        let schema = trace_schema(&example)?;
        let batches = to_record_batch(&example, &schema)?;
        let writer = File::create(output)?;
        let writer = ArrowWriter::try_new(writer, batches.schema(), None)?;
        Ok(CollapseOptions { writer, schema })
    }

    fn save_nprs(&mut self, nprs: &[Npr]) -> Result<()> {
        let read_start = nprs.get(0).ok_or(anyhow::anyhow!("Empty nprs"))?.position;
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
        acc.iter_mut().for_each(|fnpr| {
            *fnpr.length_mut() = (stop - read_start) as usize;
        });
        let batches = to_record_batch(&acc, &self.schema)?;
        self.writer.write(&batches)?;
        Ok(())
    }

    pub(crate) fn close(mut self) -> Result<()> {
        self.writer.close()?;
        Ok(())
    }

    pub(crate) fn run<P>(&mut self, filepath: P) -> Result<()>
    where
        P: AsRef<Path>,
    {
        let file = File::open(filepath)?;
        let mut builder = csv::ReaderBuilder::new().delimiter(b'\t').from_reader(file);
        let mut npr_iter = builder.deserialize().peekable();

        // First get a single line to intialize with, if None we have reached the end of
        // the iterator
        while let Some(line) = npr_iter.next() {
            let mut acc: Vec<Npr> = Vec::new();
            let mut line: Npr = line?;

            // Peek at next lines and advance the iterator only if they belong to the same
            // read
            while let Some(read_line) = npr_iter.next_if(|new_npr: &Result<Npr, _>| {
                new_npr.as_ref().expect("Parsing failed").contig == line.contig
            }) {
                let mut read_line = read_line?;
                // Data from same position, split across two lines
                if read_line.position == line.position {
                    line.samples.append(&mut read_line.samples);
                    line.event_length += read_line.event_length;
                } else {
                    // Data from next position in read
                    acc.push(line);
                    // Update line to look for possible repeating lines after it
                    line = read_line;
                }
            }
            // We are starting a new read so save the results before continuing
            // with the next read
            self.save_nprs(&acc)?;
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

    use super::*;

    #[test]
    fn test_collapse() -> Result<()> {
        let temp_dir = TempDir::new()?;
        let filepath = "extra/single_read.eventalign.txt";
        let output = temp_dir.path().join("test");
        let mut collapse = CollapseOptions::try_new(output)?;
        collapse.run(filepath)?;
        collapse.close()?;
        Ok(())
    }
}
