use core::fmt;
use std::io::{Read, Seek};

use anyhow::Result;
use bio::{alphabets::dna, io::fasta::IndexedReader};
use fnv::FnvHashMap;

use crate::{
    arrow::{Eventalign, MetadataExt},
    motif::Motif,
};

/// Contains the genomic bases for a given position including additional
/// metadata to handle positions near the end of the genome.
/// Represents the genomic sequence for a read.
pub(crate) struct Context {
    context: Vec<u8>,
    read_start: u64,
    start_slop: u64,
    end_slop: u64,
}

impl fmt::Debug for Context {
    fn fmt(&self, fmt: &mut fmt::Formatter) -> fmt::Result {
        fmt.debug_struct("Context")
            .field("context", &std::str::from_utf8(&self.context).unwrap())
            .field("read_start", &self.read_start)
            .field("start_slop", &self.start_slop)
            .finish()
    }
}

impl Context {
    pub(crate) fn new(context: Vec<u8>, read_start: u64, start_slop: u64, end_slop: u64) -> Self {
        Self {
            context,
            read_start,
            start_slop,
            end_slop,
        }
    }

    /// Genome fasta reader method makes clippy think its wrong but it still
    /// works correctly.
    #[allow(clippy::read_zero_byte_vec)]
    pub(crate) fn from_read<R>(
        genome: &mut IndexedReader<R>,
        chrom_lens: &FnvHashMap<String, u64>,
        read: &Eventalign,
    ) -> Result<Self>
    where
        R: Read + Seek,
    {
        let chrom = read.chrom();
        let chrom_len = *chrom_lens
            .get(chrom)
            .expect("chromosome missing in chrom_lens, different genome used?");
        let start_slop = read.start_0b().min(5);

        let start = if read.start_0b() < 5 {
            0
        } else {
            read.start_0b() - 5
        };

        let stop = read.seq_stop_1b_excl();
        // let end_slop = if (stop + 1) > chrom_len {
        //     0
        // } else {
        //     5.min(chrom_len - (stop + 1))
        // };
        // let stop = if (stop + 1) > chrom_len {
        //     chrom_len
        // } else {
        //     stop + 1
        // };
        genome.fetch(chrom, start, stop)?;
        let mut seq = Vec::new();

        genome.read(&mut seq)?;

        if read.strand().is_minus_strand() {
            log::debug!("Read is on negative");
            seq = seq.into_iter().map(dna::complement).collect();
        }

        Ok(Context::new(seq, read.start_0b(), start_slop, 0u64))
    }

    pub(crate) fn surrounding(&self, pos: u64, motif: &Motif) -> Vec<&[u8]> {
        let true_pos = (pos - self.read_start) + self.start_slop + motif.position_0b() as u64;

        let true_start = if true_pos < 5 {
            0
        } else {
            true_pos - 5
        };

        let mut acc = Vec::new();
        let ctxt_len = self.context.len() as u64;
        for base_pos in true_start..=true_pos {
            if (base_pos + 5) < ctxt_len {
                let base_pos = base_pos as usize;
                acc.push(&self.context[base_pos..=base_pos + 5]);
            }
        }
        acc
    }

    /// Returns None if the position is near the end of the chromosome and it
    /// would return a position with a kmer size less than six
    pub(crate) fn sixmer_at(&self, pos: u64) -> Option<&[u8]> {
        let true_pos = (pos - self.read_start) + self.start_slop;
        let true_pos = true_pos as usize;
        self.context.get(true_pos..=true_pos + 5)
    }

    pub(crate) fn start_slop(&self) -> u64 {
        self.start_slop
    }

    pub(crate) fn end_slop(&self) -> u64 {
        self.end_slop
    }
}

#[cfg(test)]
mod test {
    // use std::io::Cursor;

    // use super::*;
    // use crate::{
    //     arrow::{MetadataExt, Strand},
    //     utils::chrom_lens,
    // };

    // #[test]
    // fn test_context() -> Result<(), anyhow::Error> {
    //     u
    // }
}
