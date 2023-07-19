use std::{collections::HashMap, hash::Hash};

struct Read {
    metadata: Metadata,
    modifications: HashMap<Tag, Vec<usize>>,
    likelihoods: Vec<f64>,
}

impl Read {
    fn tag_scores(&self, tag: &Tag) -> HashMap<u64, f64> {
        todo!()
    }
}

struct Metadata {
    seq: String,
}

impl Metadata {
    fn top_positions(&self) -> HashMap<Nucleotide, Vec<usize>> {
        todo!()
    }

    fn bottom_positions(&self) -> HashMap<Nucleotide, Vec<usize>> {
        todo!()
    }
}

struct Tag {
    n_before: usize,
    base: Nucleotide,
    strand: TagStrand,
    modification: String,
}

impl Hash for Tag {
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        self.base.hash(state);
        self.strand.hash(state);
        self.modification.hash(state);
    }
}

#[derive(PartialEq, Eq, Hash)]
enum Nucleotide {
    A,
    C,
    G,
    T,
    U,
    N,
}

impl Nucleotide {
    fn dna_complement(&self) -> Self {
        match self {
            Nucleotide::A => Nucleotide::T,
            Nucleotide::C => Nucleotide::G,
            Nucleotide::G => Nucleotide::C,
            Nucleotide::T => Nucleotide::A,
            Nucleotide::U => Nucleotide::A,
            Nucleotide::N => Nucleotide::N,
        }
    }

    fn rna_complement(&self) -> Self {
        match self {
            Nucleotide::A => Nucleotide::U,
            _ => self.dna_complement(),
        }
    }
}

#[derive(PartialEq, Eq, Hash)]
enum TagStrand {
    Top,
    Bottom,
}