use std::{
    collections::{hash_map::Entry, HashMap},
    fs::File,
};

use anyhow::Result;
use indicatif::{ProgressBar, ProgressIterator, ProgressStyle};
use parquet::{arrow::ArrowWriter, file::properties::WriterProperties};
use serde::{Deserialize, Serialize};

#[derive(Clone, Debug, Serialize, Deserialize)]
pub(crate) struct NanopolishRecord {
    contig: String,
    position: u32,
    read_name: String,
    kmer: String,
    event_mean: f32,
}

impl NanopolishRecord {
    fn new(
        contig: String,
        position: u32,
        read_name: String,
        kmer: String,
        event_mean: f32,
    ) -> Self {
        NanopolishRecord {
            contig,
            position,
            read_name,
            kmer,
            event_mean,
        }
    }

    pub(crate) fn contig(&self) -> &str {
        &self.contig
    }

    pub(crate) fn pos_as_u64(&self) -> u64 {
        self.position as u64
    }

    pub(crate) fn event_mean(&self) -> f32 {
        self.event_mean
    }

    fn empty(contig: &str, position: u32, read_name: &str, kmer: &str) -> Self {
        NanopolishRecord::new(
            contig.to_owned(),
            position,
            read_name.to_owned(),
            kmer.to_owned(),
            f32::default(),
        )
    }
}

/// Reads a nanopolish output file, collapses samples across multiple lines,
/// finds the mean for all the samples, and produces a NanopolishRecord for
/// each collapsed line.
pub(crate) fn preprocess(
    input: &str,
    chrom: &Option<String>,
    start: &Option<u32>,
    stop: &Option<u32>,
) -> Result<Vec<NanopolishRecord>> {
    let source = File::open(input)?;
    let bar_style =
        ProgressStyle::default_spinner().template("[{elapsed_precise}] {spinner} {msg}");
    let bar = ProgressBar::new_spinner()
        .with_message("Processing nanopolish file.")
        .with_style(bar_style);
    let input = bar.wrap_read(source);

    let mut position_to_events: HashMap<(String, u32), Vec<f32>> = HashMap::new();
    let mut position_to_record: HashMap<(String, u32), NanopolishRecord> = HashMap::new();
    let mut reader = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .from_reader(input);
    // TODO Use serde_with::StringDelimiter
    for result in reader.records() {
        let record = result?;
        let contig = record.get(0).unwrap();

        // Filter records that are not on chromosome
        if let Some(chrom) = chrom {
            if chrom == contig {
                continue;
            }
        }
        let position = record.get(1).unwrap().parse::<u32>()?;
        if let &Some(start) = start {
            if position < start {
                continue;
            }
        }
        if let &Some(stop) = stop {
            if position > stop {
                continue;
            }
        }

        let read_name = record.get(3).unwrap();
        let kmer = record.get(9).unwrap();

        let mut event_lvls: Vec<f32> = record
            .get(13)
            .unwrap()
            .split(',')
            .map(|x| x.parse::<f32>())
            .collect::<Result<_, _>>()?;
        let key = (read_name.to_owned(), position);
        match position_to_events.entry(key.clone()) {
            Entry::Occupied(mut oe) => {
                oe.get_mut().append(&mut event_lvls);
            }
            Entry::Vacant(ve) => {
                ve.insert(event_lvls);
            }
        }

        position_to_record
            .entry(key)
            .or_insert_with(|| NanopolishRecord::empty(contig, position, read_name, kmer));
    }
    bar.finish();

    let bar_style = ProgressStyle::default_bar()
        .on_finish(indicatif::ProgressFinish::Abandon)
        .template("[{elapsed_precise}] {bar} {pos:>7}/{len:7} {msg}");
    let bar = ProgressBar::new(position_to_events.len() as u64)
        .with_message("Getting means.")
        .with_style(bar_style);
    position_to_events
        .into_iter()
        .progress_with(bar)
        .for_each(|(k, v)| {
            let n = v.len() as f32;
            let x: f32 = v.iter().sum();
            let mean = x / n;
            if let Some(np) = position_to_record.get_mut(&k) {
                np.event_mean = mean;
            } else {
                log::error!("Event in position_to_events but not in position_to_record");
            }
        });
    Ok(position_to_record.values().cloned().collect())
}

/// Writes out a Vec<NanopolishRecord> to a .parquet file
pub(crate) fn write_records_to_parquet(output: &str, records: Vec<NanopolishRecord>) -> Result<()> {
    eprintln!("Writing to .parquet file: {}", output);
    let schema = serde_arrow::Schema::from_records(&records)?;
    let batches = serde_arrow::to_record_batch(&records, &schema)?;
    let out = File::create(output)?;
    let props = WriterProperties::builder().build();
    let mut writer = ArrowWriter::try_new(out, batches.schema(), Some(props))?;
    writer.write(&batches)?;
    writer.close()?;
    Ok(())
}
