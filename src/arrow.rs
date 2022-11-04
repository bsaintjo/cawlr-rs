use std::{
    borrow::Borrow,
    fmt::Display,
    io::{Read, Seek, Write},
    ops::Index,
    slice::SliceIndex,
};

use arrow2::{
    array::Array,
    chunk::Chunk,
    datatypes::{Field, Schema},
    io::ipc::{
        read::{read_file_metadata, FileReader},
        write::{Compression, FileWriter, WriteOptions},
    },
};
use arrow2_convert::{
    deserialize::{arrow_array_deserialize_iterator, ArrowDeserialize, TryIntoCollection},
    field::ArrowField,
    serialize::{ArrowSerialize, TryIntoArrow},
    ArrowField,
};
use eyre::Result;
use itertools::Itertools;
use rv::traits::ContinuousDistr;

/// Trait for getting read information
pub trait MetadataExt {
    fn metadata(&self) -> &Metadata;

    fn name(&self) -> &str {
        self.metadata().name.as_ref()
    }

    /// Chromosome read aligned to
    fn chrom(&self) -> &str {
        self.metadata().chrom.as_ref()
    }

    /// Zero-based start of read alignment (bed-like compatible)
    fn start_0b(&self) -> u64 {
        self.metadata().start
    }

    /// One-based start of read alignment (just start_0b + 1)
    fn start_1b(&self) -> u64 {
        self.metadata().start + 1
    }

    /// Read length where nanopolish output ends
    ///
    /// See [MetadataExt::seq_length] for more detail.
    fn np_length(&self) -> u64 {
        self.metadata().length
    }

    /// Read strand
    fn strand(&self) -> Strand {
        self.metadata().strand
    }

    fn seq_stop_1b_excl(&self) -> u64 {
        self.metadata().start + self.seq_length()
    }

    /// One-based exclusive position, useful for bed-like outputs
    /// stop)
    fn end_1b_excl(&self) -> u64 {
        self.seq_stop_1b_excl() - 5
    }

    /// Length of the entire read
    ///
    /// nanopolish outputs data in 6-mers only, and positions for only the
    /// beginning of the kmer.
    ///
    /// This means the true length of the sequence of the read is 5 + the end of
    /// this output, which this method provides.
    fn seq_length(&self) -> u64 {
        self.metadata().length + 5
    }
}

impl MetadataExt for Metadata {
    fn metadata(&self) -> &Metadata {
        self
    }
}
impl MetadataMutExt for Metadata {
    fn metadata_mut(&mut self) -> &mut Metadata {
        self
    }
}

impl MetadataExt for ScoredRead {
    fn metadata(&self) -> &Metadata {
        &self.metadata
    }
}

impl MetadataMutExt for ScoredRead {
    fn metadata_mut(&mut self) -> &mut Metadata {
        &mut self.metadata
    }
}

pub trait MetadataMutExt {
    fn metadata_mut(&mut self) -> &mut Metadata;

    fn strand_mut(&mut self) -> &mut Strand {
        &mut self.metadata_mut().strand
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
    pub fn new(
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

    pub(crate) fn samples(&self) -> &[f64] {
        &self.samples
    }

    pub(crate) fn score_lnsum<M, N>(&self, pm: &M, nm: &N) -> (f64, f64)
    where
        M: ContinuousDistr<f64>,
        N: ContinuousDistr<f64>,
    {
        self.samples()
            .iter()
            .flat_map(|x| {
                let likelihood_neg = nm.ln_pdf(x);
                let likelihood_pos = pm.ln_pdf(x);
                if likelihood_neg > -10.0 && likelihood_pos > -10.0 {
                    Some((likelihood_pos, likelihood_neg))
                } else {
                    None
                }
            })
            .fold((0.0, 0.0), |acc, elem| (acc.0 + elem.0, acc.1 + elem.1))
    }
}

/// Read orientation relative to a genome
#[derive(Debug, Copy, Clone, ArrowField, PartialEq, Eq)]
#[arrow_field(type = "dense")]
pub enum Strand {
    /// (+) strand
    Plus,
    /// (-) strand
    Minus,
    /// Unable to determine strand
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
    pub fn rgb_str(&self) -> &str {
        match self {
            Strand::Plus => "255,0,0",
            Strand::Minus => "0,0,255",
            Strand::Unknown => "0,0,0",
        }
    }

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

    pub(crate) fn is_minus_strand(&self) -> bool {
        matches!(self, Strand::Minus)
    }

    pub(crate) fn is_unknown_strand(&self) -> bool {
        matches!(self, Strand::Unknown)
    }

    pub(crate) fn as_str(&self) -> &'static str {
        match self {
            Self::Plus => "+",
            Self::Minus => "-",
            Self::Unknown => ".",
        }
    }
}

/// Output representing a single read from nanopolish eventalign
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

    pub fn schema() -> Schema {
        let data_type = Self::data_type();
        Schema::from(vec![Field::new("eventalign", data_type, false)])
    }

    /// Get a mutable reference to the eventalign's signal data.
    pub(crate) fn signal_data_mut(&mut self) -> &mut Vec<Signal> {
        &mut self.signal_data
    }

    /// Iterate over Signal data
    pub(crate) fn signal_iter(&self) -> impl Iterator<Item = &Signal> {
        self.signal_data.iter()
    }

    pub(crate) fn metadata(&self) -> &Metadata {
        &self.metadata
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
    pub fn new(
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

    pub fn with_score(mut self, score: f64) -> Self {
        self.score = score;
        self
    }

    pub(crate) fn signal_score(&self) -> &Option<f64> {
        &self.signal_score
    }

    /// Get a reference to the score's kmer.
    pub fn kmer(&self) -> &str {
        self.kmer.as_ref()
    }

    pub fn pos(&self) -> u64 {
        self.pos
    }
}

/// Represents a single read scored by cawlr score
#[derive(Debug, Clone, ArrowField, Default)]
pub struct ScoredRead {
    metadata: Metadata,
    scores: Vec<Score>,
}

impl ScoredRead {
    pub fn new(metadata: Metadata, scores: Vec<Score>) -> Self {
        ScoredRead { metadata, scores }
    }

    /// Creates new ScoredRead using metadata from Eventalign output
    pub(crate) fn from_read_with_scores(eventalign: Eventalign, scores: Vec<Score>) -> Self {
        let metadata = eventalign.metadata;
        ScoredRead::new(metadata, scores)
    }

    /// Schema used for outputing into Arrow file
    pub fn schema() -> Schema {
        let data_type = Self::data_type();
        Schema::from(vec![Field::new("scored", data_type, false)])
    }

    pub fn scores(&self) -> &[Score] {
        &self.scores
    }

    pub(crate) fn to_expanded_scores(&self) -> ExpandedScores {
        let n = self.np_length() as usize;
        let mut acc = vec![None; n];
        for score in self.scores.iter() {
            let rel_pos = score.pos - self.metadata.start;
            if rel_pos >= self.np_length() {
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

#[cfg(test)]
mod test {
    use std::io::Cursor;

    use arrow2_convert::deserialize::TryIntoCollection;
    use bio::io::fasta::IndexedReader;

    use crate::{wrap_writer, save, arrow_utils::load};

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
