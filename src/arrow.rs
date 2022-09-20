use std::{
    io::{Read, Seek, Write},
    ops::Index,
    slice::SliceIndex,
    borrow::Borrow, fmt::Display,
};

use anyhow::Result;
use arrow2::{
    array::Array,
    chunk::Chunk,
    datatypes::{Field, Schema},
    error::Error,
    io::ipc::{
        read::{read_file_metadata, FileReader},
        write::{FileWriter, WriteOptions, Compression},
    },
};
use arrow2_convert::{
    deserialize::{ArrowDeserialize, TryIntoCollection, arrow_array_deserialize_iterator},
    field::ArrowField,
    serialize::{ArrowSerialize, TryIntoArrow},
    ArrowField,
};
use itertools::Itertools;

pub trait MetadataExt {
    fn metadata(&self) -> &Metadata;

    fn name(&self) -> &str {
        self.metadata().name.as_ref()
    }

    /// Get a reference to the metadata's chrom.
    fn chrom(&self) -> &str {
        self.metadata().chrom.as_ref()
    }

    /// Get the metadata's start.
    fn start_0b(&self) -> u64 {
        self.metadata().start
    }

    /// Get the metadata's start.
    fn start_1b(&self) -> u64 {
        self.metadata().start + 1
    }

    /// Get the metadata's length.
    fn np_length(&self) -> u64 {
        self.metadata().length
    }

    /// Get the metadata's strand.
    fn strand(&self) -> Strand {
        self.metadata().strand
    }

    /// stop)
    fn seq_stop_1b_excl(&self) -> u64 {
        self.metadata().start + self.seq_length()
    }

    /// stop]
    fn seq_length(&self) -> u64 {
        self.metadata().length + 5
    }

    fn end_1b_excl(&self) -> u64 {
        self.seq_stop_1b_excl() - 5
    }
}

impl MetadataExt for Metadata {
    fn metadata(&self) -> &Metadata {
        self
    }
}

impl MetadataExt for ScoredRead {
    fn metadata(&self) -> &Metadata {
        &self.metadata
    }
}

impl MetadataExt for Eventalign {
    fn metadata(&self) -> &Metadata {
        &self.metadata
    }
}

/// Represents the genomic coordinates and other information about a sequencing
/// read.
///
/// Note: All coordinate data will be zero-based for the start and one based
/// (zero-based not inclusive) for the end
#[derive(Debug, Clone, ArrowField, Default, PartialEq, Eq)]
pub struct Metadata {
    name: String,
    chrom: String,
    start: u64,
    length: u64,
    strand: Strand,
    seq: String,
}

impl Metadata {
    pub(crate) fn new(
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

    // /// Get a reference to the metadata's name.
    // pub(crate) fn name(&self) -> &str {
    //     self.name.as_ref()
    // }

    // /// Get a reference to the metadata's chrom.
    // pub(crate) fn chrom(&self) -> &str {
    //     self.chrom.as_ref()
    // }

    // /// Get the metadata's start.
    // pub(crate) fn start(&self) -> u64 {
    //     self.start
    // }

    // /// Get the metadata's length.
    // pub(crate) fn length(&self) -> u64 {
    //     self.length
    // }

    // /// Get the metadata's strand.
    // pub(crate) fn strand(&self) -> Strand {
    //     self.strand
    // }
}

#[derive(Debug, Clone, ArrowField, Default, PartialEq)]
pub(crate) struct Signal {
    pos: u64,
    kmer: String,
    signal_mean: f64,
    signal_time: f64,
    samples: Vec<f64>,
}

impl Signal {
    pub(crate) fn new(
        pos: u64,
        kmer: String,
        signal_mean: f64,
        signal_time: f64,
        samples: Vec<f64>,
    ) -> Self {
        Self {
            pos,
            kmer,
            signal_mean,
            signal_time,
            samples,
        }
    }

    /// Get the signal's signal mean.
    pub(crate) fn mean(&self) -> f64 {
        self.signal_mean
    }

    /// Get a reference to the signal's kmer.
    pub(crate) fn kmer(&self) -> &str {
        self.kmer.as_ref()
    }

    /// Get the signal's pos.
    pub(crate) fn pos(&self) -> u64 {
        self.pos
    }

    /// Set the signal's kmer.
    pub(crate) fn set_kmer(&mut self, kmer: String) {
        self.kmer = kmer;
    }
}

#[derive(Debug, Copy, Clone, ArrowField, PartialEq, Eq)]
// Currently arrow2 doesn't support newtypes or enums so for now it is a struct
#[arrow_field(type = "dense")]
pub enum Strand {
    Plus,
    Minus,
    Unknown,
}

impl Default for Strand {
    fn default() -> Self {
        Strand::Unknown
    }
}

impl Display for Strand {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.as_str())
    }
}

impl Strand {
    pub(crate) fn plus() -> Self {
        // Strand { strand: 1i8 }
        Strand::Plus
    }

    pub(crate) fn minus() -> Self {
        Strand::Minus
    }

    pub(crate) fn unknown() -> Self {
        Strand::Unknown
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

    pub(crate) fn as_str(&self) -> &'static str {
        match self {
            Self::Plus => "+",
            Self::Minus => "-",
            Self::Unknown => "."
        }
    }
}

#[derive(Debug, Clone, ArrowField, Default, PartialEq)]
pub struct Eventalign {
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

    pub fn name(&self) -> &str {
        self.metadata.name.as_str()
    }

    /// 1-based indexing
    pub(crate) fn start_1b(&self) -> u64 {
        self.metadata.start + 1
    }

    /// 0-based indexing
    pub(crate) fn start_0b(&self) -> u64 {
        self.metadata.start
    }

    pub(crate) fn length_mut(&mut self) -> &mut u64 {
        &mut self.metadata.length
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

#[derive(Default)]
pub struct ScoreBuilder {
    score: Score,
}

impl ScoreBuilder {
    pub fn new(score: Score) -> Self {
        Self { score }
    }

    pub fn score(mut self, score: f64) -> Self {
        self.score.score = score;
        self
    }

    pub fn build(self) -> Score {
        self.score
    }
}

#[derive(Default, Debug, Clone, ArrowField)]
pub struct Score {
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
    pub fn score(&self) -> f64 {
        self.score
    }

    pub(crate) fn signal_score(&self) -> &Option<f64> {
        &self.signal_score
    }

    /// Get a reference to the score's kmer.
    pub fn kmer(&self) -> &str {
        self.kmer.as_ref()
    }
}

#[derive(Debug, Clone, ArrowField, Default)]
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

    pub fn schema() -> Schema {
        let data_type = Self::data_type();
        Schema::from(vec![Field::new("scored", data_type, false)])
    }

    pub fn scores(&self) -> &[Score] {
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

pub fn wrap_writer<W>(writer: W, schema: &Schema) -> Result<FileWriter<W>>
where
    W: Write,
{
    let options = WriteOptions { compression: Some(Compression::LZ4) };
    let fw = FileWriter::try_new(writer, schema, None, options)?;
    Ok(fw)
}

pub fn save<W, T>(writer: &mut FileWriter<W>, x: &[T]) -> Result<()>
where
    T: ArrowField<Type = T> + ArrowSerialize + 'static,
    W: Write,
{
    let arrow_array: Chunk<Box<dyn Array>> = x.try_into_arrow()?;
    writer.write(&arrow_array, None)?;
    Ok(())
}

pub(crate) fn load<R>(mut reader: R) -> Result<FileReader<R>>
where
    R: Read + Seek,
{
    let metadata = read_file_metadata(&mut reader)?;
    let reader = FileReader::new(reader, metadata, None, None);
    Ok(reader)
}

/// Apply a function to chunks of data loaded from an Arrow Feather File.
///
/// # Example
/// ```rust, ignore
/// # use std::fs::File;
/// # use std::error::Error;
/// # use cawlr::arrow::Eventalign;
/// # use cawlr::arrow::load_apply;
/// # fn main() -> Result<(), Box<dyn Error>> {
/// load_apply(file, |eventalign: Vec<Eventalign>| {
///     // Do stuff with each chunk
/// # Ok(())
/// })
/// # }
/// ```
///
/// TODO Fix example with correct file
pub fn load_apply<R, F, T>(reader: R, mut func: F) -> Result<()>
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

pub fn load_apply_indy<R, F, T>(reader: R, mut func: F) -> Result<()>
where
    R: Read + Seek,
    F: FnMut(T) -> anyhow::Result<()>,
    T: ArrowField<Type = T> + ArrowDeserialize + 'static,
    for<'a> &'a <T as ArrowDeserialize>::ArrayType: IntoIterator,
{
    let feather = load(reader)?;
    for read in feather {
        if let Ok(chunk) = read {
            for arr in chunk.into_arrays().into_iter() {
                let iter = arrow_array_deserialize_iterator(arr.borrow())?;
                for x in iter {
                    func(x)?;
                }
            }
        } else {
            log::warn!("Failed to load arrow chunk")
        }
    }
    Ok(())
}

/// Loops over every chunk in an arrow file, applies the closure to the chunk,
/// and writes the results to the writer.
///
/// TODO make F: ... Result<&[U]> instead to avoid allocation?
pub fn load_read_write<R, W, F, T, U>(
    reader: R,
    mut writer: FileWriter<W>,
    mut func: F,
) -> Result<()>
where
    R: Read + Seek,
    W: Write,
    F: FnMut(Vec<T>) -> anyhow::Result<Vec<U>>,
    T: ArrowField<Type = T> + ArrowDeserialize + 'static,
    U: ArrowField<Type = U> + ArrowSerialize + 'static,
    for<'a> &'a <T as ArrowDeserialize>::ArrayType: IntoIterator,
{
    let feather = load(reader)?;
    for read in feather {
        if let Ok(chunk) = read {
            for arr in chunk.into_arrays().into_iter() {
                let eventaligns: Vec<T> = arr.try_into_collection()?;
                let res = func(eventaligns)?;
                save(&mut writer, &res)?;
            }
        } else {
            log::warn!("Failed to load arrow chunk")
        }
    }
    writer.finish()?;
    Ok(())
}

// TODO Refactor multiple maps
pub(crate) fn load_iter<R>(mut reader: R) -> impl Iterator<Item = Result<Vec<Eventalign>, Error>>
where
    R: Read + Seek,
{
    let metadata = read_file_metadata(&mut reader).unwrap();
    let reader = FileReader::new(reader, metadata, None, None);
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

pub(crate) fn load_iter2<R>(
    mut reader: R,
) -> impl Iterator<Item = Result<impl Iterator<Item = Result<Vec<Eventalign>>>, arrow2::error::Error>>
where
    R: Read + Seek,
{
    let metadata = read_file_metadata(&mut reader).unwrap();
    let reader = FileReader::new(reader, metadata, None, None);
    reader.map_ok(|c| {
        c.into_arrays().into_iter().map(|a| {
            let x: Vec<Eventalign> = a.try_into_collection_as_type::<Eventalign>()?;
            Ok(x)
        })
    })
}

#[cfg(test)]
mod test {
    use std::io::Cursor;

    use arrow2_convert::deserialize::TryIntoCollection;
    use bio::io::fasta::IndexedReader;

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

        let signal = Signal::new(1u64, "AAAAAA".to_string(), 80.0f64, 0.01f64, Vec::new());

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

    #[allow(clippy::read_zero_byte_vec)]
    #[test]
    fn test_fasta_reader_start() {
        let genome = "extra/sacCer3.fa";
        let mut genome = IndexedReader::from_file(&genome).unwrap();
        let chrom = "chrI";
        let start = 71071;
        let stop = start + 6;
        genome.fetch(chrom, start, stop).unwrap();
        let mut seq = Vec::new();
        genome.read(&mut seq).unwrap();

        assert_eq!(b"GCAAGC", seq.as_slice());
    }
}
