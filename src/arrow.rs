use std::{
    io::{Read, Seek, Write},
    ops::Index,
    slice::SliceIndex,
    sync::Arc,
};

use anyhow::Result;
use arrow2::{
    array::Array,
    chunk::Chunk,
    datatypes::{Field, Schema},
    error::Error,
    io::ipc::{
        read::{read_file_metadata, FileReader},
        write::{FileWriter, WriteOptions},
    },
};
use arrow2_convert::{
    deserialize::{ArrowDeserialize, TryIntoCollection},
    field::ArrowField,
    serialize::{ArrowSerialize, TryIntoArrow},
    ArrowField,
};

#[derive(Debug, Clone, ArrowField, Default)]
pub(crate) struct Metadata {
    name: String,
    chrom: String,
    start: u64,
    length: u64,
    strand: Strand,
    seq: String,
}

impl Metadata {
    fn new(
        name: String,
        chrom: String,
        start: u64,
        length: u64,
        strand: Strand,
        seq: String,
    ) -> Self {
        Self {
            name,
            chrom,
            start,
            length,
            strand,
            seq,
        }
    }

    /// Get a reference to the metadata's name.
    #[must_use]
    pub(crate) fn name(&self) -> &str {
        self.name.as_ref()
    }

    /// Get a reference to the metadata's chrom.
    #[must_use]
    pub(crate) fn chrom(&self) -> &str {
        self.chrom.as_ref()
    }

    /// Get the metadata's start.
    #[must_use]
    pub(crate) fn start(&self) -> u64 {
        self.start
    }

    /// Get the metadata's length.
    #[must_use]
    pub(crate) fn length(&self) -> u64 {
        self.length
    }

    /// Get the metadata's strand.
    #[must_use]
    pub(crate) fn strand(&self) -> Strand {
        self.strand
    }
}

#[derive(Debug, Clone, ArrowField)]
pub(crate) struct Signal {
    pos: u64,
    kmer: String,
    signal_mean: f64,
    signal_time: f64,
}

impl Signal {
    pub(crate) fn new(pos: u64, kmer: String, signal_mean: f64, signal_time: f64) -> Self {
        Self {
            pos,
            kmer,
            signal_mean,
            signal_time,
        }
    }

    /// Get the signal's signal mean.
    #[must_use]
    pub(crate) fn mean(&self) -> f64 {
        self.signal_mean
    }

    /// Get a reference to the signal's kmer.
    #[must_use]
    pub(crate) fn kmer(&self) -> &str {
        self.kmer.as_ref()
    }

    /// Get the signal's pos.
    #[must_use]
    pub(crate) fn pos(&self) -> u64 {
        self.pos
    }

    /// Set the signal's kmer.
    pub(crate) fn set_kmer(&mut self, kmer: String) {
        self.kmer = kmer;
    }
}

#[derive(Default, Debug, Copy, Clone, ArrowField, PartialEq, Eq)]
// Currently arrow2 doesn't support newtypes or enums so for now it is a struct
pub(crate) struct Strand {
    strand: i8,
}

impl Strand {
    pub(crate) fn plus() -> Self {
        Strand { strand: 1i8 }
    }

    pub(crate) fn minus() -> Self {
        Strand { strand: -1i8 }
    }

    pub(crate) fn unknown() -> Self {
        Strand { strand: 0i8 }
    }

    pub(crate) fn is_plus_strand(&self) -> bool {
        self == &Strand::plus()
    }

    pub(crate) fn is_minus_strand(&self) -> bool {
        self == &Strand::minus()
    }

    pub(crate) fn is_unknown_strand(&self) -> bool {
        self == &Strand::unknown()
    }

    pub(crate) fn as_str(&self) -> &str {
        if self.is_plus_strand() {
            "+"
        } else if self.is_minus_strand() {
            "-"
        } else {
            "."
        }
    }
}

#[derive(Debug, Clone, ArrowField)]
pub(crate) struct Eventalign {
    metadata: Metadata,
    signal_data: Vec<Signal>,
}

impl Eventalign {
    pub(crate) fn new(metadata: Metadata, signal_data: Vec<Signal>) -> Self {
        Self {
            metadata,
            signal_data,
        }
    }

    pub(crate) fn name(&self) -> &str {
        self.metadata.name.as_str()
    }

    /// 1-based indexing
    pub(crate) fn start_ob(&self) -> u64 {
        self.metadata.start + 1
    }

    /// 0-based indexing
    pub(crate) fn start_zb(&self) -> u64 {
        self.metadata.start
    }

    pub(crate) fn length_mut(&mut self) -> &mut u64 {
        &mut self.metadata.length
    }

    pub(crate) fn length(&self) -> u64 {
        self.metadata.length
    }

    /// 1-based indexing inclusive stop
    pub(crate) fn stop_ob(&self) -> u64 {
        self.start_ob() + self.length() - 1
    }

    /// 0-based indexing inclusive stop
    pub(crate) fn stop_zb(&self) -> u64 {
        self.start_zb() + self.length() - 1
    }

    pub(crate) fn seq_stop_zb(&self) -> u64 {
        self.stop_zb() + 5
    }

    pub(crate) fn chrom(&self) -> &str {
        self.metadata.chrom.as_str()
    }

    pub(crate) fn strand(&self) -> Strand {
        self.metadata.strand
    }

    pub(crate) fn strand_mut(&mut self) -> &mut Strand {
        &mut self.metadata.strand
    }

    pub(crate) fn empty(name: String, chrom: String, start: u64, length: u64, seq: String) -> Self {
        let strand = Strand::unknown();
        let metadata = Metadata::new(name, chrom, start, length, strand, seq);
        Eventalign::new(metadata, Vec::new())
    }

    pub(crate) fn schema() -> Schema {
        let data_type = Self::data_type();
        Schema::from(vec![Field::new("eventalign", data_type, false)])
    }

    /// Get a mutable reference to the eventalign's signal data.
    #[must_use]
    pub(crate) fn signal_data_mut(&mut self) -> &mut Vec<Signal> {
        &mut self.signal_data
    }

    pub(crate) fn signal_iter(&self) -> impl Iterator<Item = &Signal> {
        self.signal_data.iter()
    }

    pub(crate) fn metadata(&self) -> &Metadata {
        &self.metadata
    }
}

#[derive(Default, Debug, Clone, ArrowField)]
pub(crate) struct Score {
    pos: u64,
    kmer: String,
    skipped: bool,
    signal_score: Option<f64>,
    skip_score: f64,
    score: f64,
}

impl Score {
    pub(crate) fn new(
        pos: u64,
        kmer: String,
        skipped: bool,
        signal_score: Option<f64>,
        skip_score: f64,
        score: f64,
    ) -> Self {
        Self {
            pos,
            kmer,
            skipped,
            signal_score,
            skip_score,
            score,
        }
    }

    /// Get the score's score.
    #[must_use]
    pub(crate) fn score(&self) -> f64 {
        self.score
    }

    /// Get a reference to the score's kmer.
    #[must_use]
    pub(crate) fn kmer(&self) -> &str {
        self.kmer.as_ref()
    }
}

#[derive(Debug, Clone, ArrowField)]
pub struct ScoredRead {
    metadata: Metadata,
    scores: Vec<Score>,
}

impl ScoredRead {
    fn new(metadata: Metadata, scores: Vec<Score>) -> Self {
        ScoredRead { metadata, scores }
    }

    pub(crate) fn metadata(&self) -> &Metadata {
        &self.metadata
    }

    pub(crate) fn length(&self) -> u64 {
        self.metadata.length
    }

    pub(crate) fn from_read_with_scores(eventalign: Eventalign, scores: Vec<Score>) -> Self {
        let metadata = eventalign.metadata;
        ScoredRead::new(metadata, scores)
    }

    pub(crate) fn schema() -> Schema {
        let data_type = Self::data_type();
        Schema::from(vec![Field::new("scored", data_type, false)])
    }

    pub(crate) fn scores(&self) -> &[Score] {
        &self.scores
    }

    pub(crate) fn to_expanded_scores(&self) -> ExpandedScores {
        let n = self.length() as usize;
        let mut acc = vec![None; n];
        for score in self.scores.iter() {
            let rel_pos = score.pos - self.metadata.start;
            if rel_pos >= self.length() {
                log::warn!(
                    "Read contains data outside length, {}, start+length is {}+{} but score at {}",
                    self.metadata.name,
                    self.metadata.start,
                    self.metadata.length,
                    score.pos
                );
                continue;
            } else {
                acc[rel_pos as usize] = Some(score);
            }
        }
        ExpandedScores(acc)
    }
}

pub(crate) struct ExpandedScores<'a>(Vec<Option<&'a Score>>);

impl<'a, I: SliceIndex<[Option<&'a Score>]>> Index<I> for ExpandedScores<'a> {
    type Output = I::Output;

    fn index(&self, index: I) -> &Self::Output {
        self.0.index(index)
    }
}

pub(crate) fn wrap_writer<W>(writer: W, schema: &Schema) -> Result<FileWriter<W>>
where
    W: Write,
{
    let options = WriteOptions { compression: None };
    let fw = FileWriter::try_new(writer, schema, None, options)?;
    Ok(fw)
}

pub(crate) fn save<W, T>(writer: &mut FileWriter<W>, x: &[T]) -> Result<()>
where
    T: ArrowField<Type = T> + ArrowSerialize + 'static,
    W: Write,
{
    let arrow_array: Arc<dyn Array> = x.try_into_arrow()?;
    let chunks = Chunk::new(vec![arrow_array]);
    writer.write(&chunks, None)?;
    Ok(())
}

pub(crate) fn load<R>(mut reader: R) -> Result<FileReader<R>>
where
    R: Read + Seek,
{
    let metadata = read_file_metadata(&mut reader)?;
    let reader = FileReader::new(reader, metadata, None);
    Ok(reader)
}

/// Apply a function to chunks of data loaded from an Arrow Feather File.
pub(crate) fn load_apply<R, F, T>(reader: R, mut func: F) -> Result<()>
where
    R: Read + Seek,
    F: FnMut(Vec<T>) -> anyhow::Result<()>,
    T: ArrowField<Type = T> + ArrowDeserialize + 'static,
    for<'a> &'a <T as ArrowDeserialize>::ArrayType: IntoIterator,
{
    let feather = load(reader)?;
    for read in feather {
        if let Ok(chunk) = read {
            for arr in chunk.into_arrays().into_iter() {
                let eventaligns: Vec<T> = arr.try_into_collection()?;
                func(eventaligns)?;
            }
        } else {
            log::warn!("Failed to load arrow chunk")
        }
    }
    Ok(())
}

// TODO Refactor multiple maps
pub(crate) fn load_iter<R>(
    mut reader: R,
) -> impl Iterator<Item = Result<Vec<Eventalign>, Error>>
where
    R: Read + Seek,
{
    let metadata = read_file_metadata(&mut reader).unwrap();
    let reader = FileReader::new(reader, metadata, None);
    reader
        .map(|x| x.map(|c| c.into_arrays().into_iter()))
        .map(|q| {
            q.map(|b| {
                b.flat_map(|r| {
                    let v: Vec<Eventalign> = r.try_into_collection().unwrap();
                    v
                })
                .collect::<Vec<_>>()
            })
        })
}

#[cfg(test)]
mod test {
    use std::io::Cursor;

    use arrow2_convert::deserialize::TryIntoCollection;

    use super::*;

    #[test]
    fn test_single_read() {}

    #[test]
    fn test_round_trip() {
        let metadata = Metadata::new(
            "abc".to_string(),
            "chrI".to_string(),
            0u64,
            100u64,
            Strand::plus(),
            String::new(),
        );

        let signal = Signal::new(1u64, "AAAAAA".to_string(), 80.0f64, 0.01f64);

        let eventalign = Eventalign::new(metadata, vec![signal]);
        let x = [eventalign.clone(), eventalign];

        let schema = Eventalign::schema();

        let file = vec![];
        let mut writer = wrap_writer(file, &schema).unwrap();
        save(&mut writer, &x).unwrap();
        writer.finish().unwrap();

        let reader = writer.into_inner();
        let reader = Cursor::new(reader);

        let filereader = load(reader).unwrap();
        for chunks in filereader.flatten() {
            for row in chunks.into_arrays().into_iter() {
                let _: Vec<Eventalign> = row.try_into_collection().unwrap();
            }
        }
    }

    #[test]
    fn test_expanded_scores() {
        let metadata = Metadata {
            start: 100,
            length: 10,
            ..Default::default()
        };
        let scores = vec![
            Score {
                pos: 101,
                ..Default::default()
            },
            Score {
                pos: 105,
                ..Default::default()
            },
        ];

        let read = ScoredRead::new(metadata, scores);
        let expanded = read.to_expanded_scores();
        assert_eq!(expanded.0.len(), 10);
        assert!(expanded.0[0].is_none());
        assert!(expanded.0[1].is_some());
        assert!(expanded.0[5].is_some());
        assert!(expanded.0.get(10).is_none());
    }

    #[test]
    fn test_expanded_scores_outside() {
        let metadata = Metadata {
            start: 100,
            length: 10,
            ..Default::default()
        };
        let scores = vec![
            Score {
                pos: 100,
                ..Default::default()
            },
            Score {
                pos: 110,
                ..Default::default()
            },
        ];

        let read = ScoredRead::new(metadata, scores);
        let expanded = read.to_expanded_scores();
        assert_eq!(expanded.0.len(), 10);
        assert!(expanded.0[0].is_some());
        assert!(expanded.0[1].is_none());
        assert!(expanded.0[9].is_none());
    }
}
