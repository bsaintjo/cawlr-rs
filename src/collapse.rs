use std::{
    fs::File,
    io::{BufWriter, Read, Write},
    path::Path,
};

use anyhow::Result;
use arrow2::io::ipc::write::FileWriter;
use bio::alphabets::dna::revcomp;
use indicatif::{ProgressBar, ProgressBarIter, ProgressFinish, ProgressStyle};
use serde::Deserialize;
use serde_with::{serde_as, StringWithSeparator};
use serde_with::formats::CommaSeparator;
use statrs::statistics::Statistics;

use crate::{
    arrow::{self, save, Eventalign, Metadata, Signal, Strand},
    plus_strand_map::PlusStrandMap,
};

fn empty_from_npr(npr: Npr) -> Eventalign {
    let name = npr.read_name().to_string();
    let chrom = npr.contig().to_string();
    let start = npr.position;
    let length = 1;
    let seq = String::new();
    let metadata = Metadata::new(name, chrom, start, length, Strand::unknown(), seq);
    let signal_data = vec![Signal::new(
        npr.position,
        npr.reference_kmer().to_string(),
        npr.samples().mean(),
        npr.event_length,
        npr.samples,
    )];

    Eventalign::new(metadata, signal_data)
}

/// Takes a vector of nanpolish records and converts them into a Eventalign.
fn nprs_to_eventalign(
    mut nprs: impl Iterator<Item = Npr>,
    strand_map: &PlusStrandMap,
) -> Result<Option<Eventalign>> {
    let mut eventalign = nprs
        .next()
        .ok_or(anyhow::anyhow!("Empty nprs"))
        .map(empty_from_npr)?;
    let mut stop = eventalign.start_0b();
    for npr in nprs {
        stop = npr.position;
        let position = npr.position;
        let ref_kmer = npr.reference_kmer().to_string();
        let mean = npr.samples().mean();

        if mean.is_nan() {
            return Err(anyhow::anyhow!("No signal samples values, malformed input"));
        }

        let time = npr.event_length;
        let signal = Signal::new(position, ref_kmer, mean, time, npr.samples);
        eventalign.signal_data_mut().push(signal);
    }

    // Update strand from bam file results
    let strand = strand_map.get(eventalign.name());
    if let Some(b) = strand {
        let strand_ptr = eventalign.strand_mut();
        *strand_ptr = if b {
            arrow::Strand::plus()
        } else {
            arrow::Strand::minus()
        }
    } else {
        log::warn!("Read {} could not find strand", eventalign.name())
    }

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

    // Reverse kmer
    if eventalign.strand().is_minus_strand() {
        for signal in eventalign.signal_data_mut().iter_mut() {
            let rev_kmer = revcomp(signal.kmer().as_bytes());
            let rev_kmer = String::from_utf8(rev_kmer)?;
            signal.set_kmer(rev_kmer)
        }
    }
    log::debug!("Parsed Eventalign: {eventalign:.2?} ");
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
        .unwrap();
    pb.with_message("Processing eventalign data")
        .with_style(style)
        .with_finish(ProgressFinish::AndLeave)
        .wrap_read(iter)
}

pub struct CollapseOptions<W: Write> {
    writer: FileWriter<W>,
    strand_db: PlusStrandMap,
    capacity: usize,
    progress: bool,
}

impl CollapseOptions<BufWriter<File>> {
    pub fn try_new<Q, R>(bam_file: Q, output: R) -> Result<Self>
    where
        Q: AsRef<Path>,
        R: AsRef<Path>,
    {
        let writer = File::create(output)?;
        let writer = BufWriter::new(writer);
        CollapseOptions::from_writer(writer, bam_file)
    }
}

impl<W: Write> CollapseOptions<W> {
    fn new(writer: FileWriter<W>, strand_db: PlusStrandMap) -> Self {
        Self {
            writer,
            strand_db,
            capacity: 2048,
            progress: false,
        }
    }

    pub fn capacity(&mut self, capacity: usize) -> &mut Self {
        self.capacity = capacity;
        self
    }

    pub fn progress(&mut self, progress: bool) -> &mut Self {
        self.progress = progress;
        self
    }

    pub(crate) fn from_writer<R>(writer: W, bam_file: R) -> Result<Self>
    where
        R: AsRef<Path>,
    {
        let strand_db = PlusStrandMap::from_bam_file(bam_file)?;
        let schema = arrow::Eventalign::schema();
        let writer = arrow::wrap_writer(writer, &schema)?;
        Ok(CollapseOptions::new(writer, strand_db))
    }

    fn save_eventalign(&mut self, eventaligns: &[Eventalign]) -> Result<()> {
        save(&mut self.writer, eventaligns)
    }

    fn close(&mut self) -> Result<()> {
        self.writer.finish()?;
        Ok(())
    }

    pub fn run<R>(&mut self, input: R) -> Result<()>
    where
        R: Read,
    {
        let file = spin_iter(input, self.progress);
        let mut builder = csv::ReaderBuilder::new().delimiter(b'\t').from_reader(file);
        let mut npr_iter = builder.deserialize();

        let mut idx_diff = 1;
        let npr: Npr = npr_iter.next().expect(
            "No data, check if eventalign has data; nanopolish eventalign may have failed",
        )?;
        let mut position = npr.position;

        let mut acc = vec![npr];
        let mut flats = Vec::with_capacity(self.capacity);

        for line in npr_iter {
            if let Ok(mut next_npr) = line {
                let last = acc.last().unwrap();
                let read_name = last.read_name();
                let event_idx = last.event_index();
                if (next_npr.read_name() == read_name)
                    && (next_npr.event_index().abs_diff(event_idx) == idx_diff)
                {
                    // Same read, possibly new kmer or same
                    if next_npr.position == position {
                        // Same read, same kmer
                        let npr_mut = acc.last_mut().unwrap();
                        npr_mut.samples.append(&mut next_npr.samples);
                        npr_mut.event_length += next_npr.event_length;
                        npr_mut.event_index = next_npr.event_index;
                    } else {
                        // Same read, different kmer
                        position = next_npr.position;
                        acc.push(next_npr);
                    }
                } else {
                    // New read, write data and move forward
                    if let Some(eventalign) = nprs_to_eventalign(acc.drain(..), &self.strand_db)? {
                        flats.push(eventalign);
                    }

                    if flats.len() >= self.capacity {
                        self.save_eventalign(&flats)?;
                        flats.clear();
                    }
                    acc.push(next_npr);
                }
                idx_diff = 1;
            } else {
                log::warn!("Parsing failed: {line:?}");
                idx_diff += 1;
            }
        }

        if !acc.is_empty() {
            if let Some(eventalign) = nprs_to_eventalign(acc.drain(..), &self.strand_db)? {
                flats.push(eventalign);
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
#[derive(Default, Clone, Debug, Deserialize, PartialEq)]
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

    #[serde(skip)]
    _model_kmer: String,

    #[serde(skip)]
    _model_mean: f64,

    #[serde(skip)]
    _model_stdv: f64,

    #[serde(skip)]
    _standardized_level: f64,

    #[serde_as(as = "StringWithSeparator::<CommaSeparator, f64>")]
    samples: Vec<f64>,
}

impl Npr {
    fn contig(&self) -> &str {
        &self.contig
    }

    fn read_name(&self) -> &str {
        &self.read_name
    }

    fn event_index(&self) -> i64 {
        self.event_index
    }

    fn samples(&self) -> &[f64] {
        &self.samples
    }

    fn reference_kmer(&self) -> &str {
        &self.reference_kmer
    }
}

#[cfg(test)]
mod test {

    use std::io::Cursor;

    use assert_fs::TempDir;

    use super::*;
    use crate::arrow::{load_apply, load_iter, Metadata, MetadataExt, Strand};

    #[test]
    fn test_collapse() -> Result<()> {
        let temp_dir = TempDir::new()?;
        let filepath = "extra/single_read.eventalign.txt";
        let input = File::open(filepath)?;
        let bam_file = "extra/single_read.bam";
        let output = temp_dir.path().join("test");
        let mut collapse = CollapseOptions::try_new(bam_file, &output)?;
        collapse.run(input)?;

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
        let input = File::open(filepath)?;
        let bam_file = "extra/neg_control.bam";
        let output = temp_dir.path().join("test");
        let mut collapse = CollapseOptions::try_new(bam_file, &output)?;
        collapse.run(input)?;

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

    #[test]
    fn test_malformed() {
        let lines: &[u8] = b"contig	position	reference_kmer	read_name	strand	event_index	event_level_mean	event_stdv	event_length	model_kmer	model_mean	model_stdv	standardized_level	samples
chr1	199403040	ATATAA	c25d27a8-0eec-4e7d-96f9-b8e730a25832	t	3919	86.81	0.500	0.00100	TTATAT	87.94	1.88	-0.59	87.1186,87.4749,86.406,86.2279
chr1	199403040	ATATAA	c25d27a8-0eec-4e7d-96f9-b8e730a25832	t	3918	87.01		72.4013,75.9601,78.395,77.6458
";
        // let file = File::open("extra/malformed.eventalign.txt").unwrap();
        let mut builder = csv::ReaderBuilder::new()
            .has_headers(true)
            .delimiter(b'\t')
            .from_reader(lines);
        let mut iter = builder.deserialize::<Npr>();

        let npr = Npr {
            contig: "chr1".to_string(),
            position: 199403040,
            reference_kmer: "ATATAA".to_string(),
            read_name: "c25d27a8-0eec-4e7d-96f9-b8e730a25832".to_string(),
            samples: vec![87.1186, 87.4749, 86.406, 86.2279],
            event_index: 3919,
            event_length: 0.00100,
            ..Default::default()
        };

        let next = iter.next().unwrap();
        assert!(next.is_ok());

        assert_eq!(next.unwrap(), npr);
        assert!(iter.next().unwrap().is_err());

        let mut strand_db = PlusStrandMap::default();
        strand_db.insert(b"c25d27a8-0eec-4e7d-96f9-b8e730a25832" as &[u8], true);

        let schema = arrow::Eventalign::schema();
        let writer = arrow::wrap_writer(Vec::new(), &schema).unwrap();
        let mut opts = CollapseOptions::new(writer, strand_db);
        let res = opts.run(lines);
        assert!(res.is_ok());

        let reader = Cursor::new(opts.writer.into_inner());
        let x = load_iter(reader).next().unwrap().unwrap();

        let target = Eventalign::new(
            Metadata::new(
                "c25d27a8-0eec-4e7d-96f9-b8e730a25832".to_string(),
                "chr1".to_string(),
                199403040,
                1,
                Strand::plus(),
                String::new(),
            ),
            vec![Signal::new(
                199403040,
                "ATATAA".to_string(),
                npr.samples().mean(),
                0.00100,
                vec![87.1186, 87.4749, 86.406, 86.2279],
            )],
        );

        pretty_assertions::assert_eq!(x[0], target);
    }

    #[test]
    fn test_diff_idx() {
        let lines: &[u8] = b"contig	position	reference_kmer	read_name	strand	event_index	event_level_mean	event_stdv	event_length	model_kmer	model_mean	model_stdv	standardized_level	samples
chr1	199403040	ATATAA	c25d27a8-0eec-4e7d-96f9-b8e730a25832	t	3919	86.81	0.500	0.00100	TTATAT	87.94	1.88	-0.59	87.1186,87.4749,86.406,86.2279
chr1	199403040	ATATAA	c25d27a8-0eec-4e7d-96f9-b8e730a25832	t	3918	87.01		72.4013,75.9601,78.395,77.6458
chr1	199403041	GATATA	c25d27a8-0eec-4e7d-96f9-b8e730a25832	t	3917	106.85	4.255	0.00100	TATATC	107.52	3.75	-0.18	99.4103,108.674,110.277,109.03
";
        let mut strand_db = PlusStrandMap::default();
        strand_db.insert(b"c25d27a8-0eec-4e7d-96f9-b8e730a25832" as &[u8], true);

        let schema = arrow::Eventalign::schema();
        let writer = arrow::wrap_writer(Vec::new(), &schema).unwrap();
        let mut opts = CollapseOptions::new(writer, strand_db);
        let res = opts.run(lines);
        assert!(res.is_ok());

        let reader = Cursor::new(opts.writer.into_inner());
        let x = load_iter(reader).next().unwrap().unwrap();

        let target = Eventalign::new(
            Metadata::new(
                "c25d27a8-0eec-4e7d-96f9-b8e730a25832".to_string(),
                "chr1".to_string(),
                199403040,
                2,
                Strand::plus(),
                String::new(),
            ),
            vec![
                Signal::new(
                    199403040,
                    "ATATAA".to_string(),
                    [87.1186, 87.4749, 86.406, 86.2279].mean(),
                    0.00100,
                    vec![87.1186, 87.4749, 86.406, 86.2279],
                ),
                Signal::new(
                    199403041,
                    "GATATA".to_string(),
                    [99.4103, 108.674, 110.277, 109.03].mean(),
                    0.00100,
                    vec![99.4103, 108.674, 110.277, 109.03],
                ),
            ],
        );
        pretty_assertions::assert_eq!(x[0], target);
    }
}
