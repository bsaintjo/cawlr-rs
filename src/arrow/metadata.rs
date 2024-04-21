use std::fmt::Display;

use arrow2_convert::{ArrowDeserialize, ArrowField, ArrowSerialize};
use serde::{Deserialize, Serialize};

/// Represents the genomic coordinates and other information about a sequencing
/// read.
///
/// Note: All coordinate data will be zero-based for the start and one based
/// (zero-based not inclusive) for the end
#[derive(
    Debug,
    Clone,
    ArrowField,
    ArrowSerialize,
    ArrowDeserialize,
    Serialize,
    Deserialize,
    Default,
    PartialEq,
    Eq,
)]
pub struct Metadata {
    pub name: String,
    pub chrom: String,
    pub start: u64,
    pub length: u64,
    pub strand: Strand,
    pub seq: String,
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
}

pub trait MetadataExt {
    fn metadata(&self) -> &Metadata;

    fn is_unaligned(&self) -> bool {
        self.metadata().chrom.is_empty()
    }

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

pub trait MetadataMutExt {
    fn metadata_mut(&mut self) -> &mut Metadata;

    fn strand_mut(&mut self) -> &mut Strand {
        &mut self.metadata_mut().strand
    }
}

/// Read orientation relative to a genome
#[derive(
    Debug,
    Copy,
    Clone,
    ArrowField,
    ArrowSerialize,
    ArrowDeserialize,
    Serialize,
    Deserialize,
    PartialEq,
    Eq,
)]
pub struct Strand {
    strand: i8,
}

impl Default for Strand {
    fn default() -> Self {
        Strand::unknown()
    }
}

impl Display for Strand {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.as_str())
    }
}

impl Strand {
    pub const fn new(strand: i8) -> Self {
        Self { strand }
    }

    pub const fn rgb_str(&self) -> &str {
        match self.strand {
            s if s > 0 => "255,0,0",
            s if s < 0 => "0,0,255",
            0 => "0,0,0",
            _ => panic!(),
        }
    }

    pub const fn plus() -> Self {
        Strand::new(1)
    }

    pub const fn minus() -> Self {
        Strand::new(-1)
    }

    pub const fn unknown() -> Self {
        Strand::new(0)
    }

    pub const fn is_minus_strand(&self) -> bool {
        self.strand < 0
    }

    pub const fn is_unknown_strand(&self) -> bool {
        self.strand == 0
    }

    pub const fn as_str(&self) -> &'static str {
        match self.strand {
            s if s > 0 => "+",
            s if s < 0 => "-",
            0 => ".",
            _ => panic!("Strand enum pattern"),
        }
    }
}
