use core::fmt;
use std::io::{Read, Seek};

use anyhow::Result;
use bio::{alphabets::dna, io::fasta::IndexedReader};
use fnv::FnvHashMap;

use crate::arrow::Eventalign;

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
        let start_slop = read.start_zb().min(5);

        let start = if read.start_zb() < 5 {
            0
        } else {
            read.start_zb() - 5
        };

        let stop = read.seq_stop_zb();
        let end_slop = if (stop + 1) > chrom_len {
            0
        } else {
            5.min(chrom_len - (stop + 1))
        };
        let stop = if (stop + 1) > chrom_len {
            chrom_len
        } else {
            stop + 1
        };
        genome.fetch(chrom, start, stop)?;
        let mut seq = Vec::new();

        genome.read(&mut seq)?;

        if read.strand().is_minus_strand() {
            log::debug!("Read is on negative");
            seq = seq.into_iter().map(dna::complement).collect();
        }

        Ok(Context::new(seq, read.start_zb(), start_slop, end_slop))
    }

    pub(crate) fn surrounding(&self, pos: u64) -> Vec<&[u8]> {
        let mut acc = Vec::new();
        let true_pos = (pos - self.read_start) + self.start_slop;

        let true_start = if true_pos < 5 { 0 } else { true_pos - 5 };

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
    use std::io::Cursor;

    use super::*;
    use crate::{arrow::Strand, utils::chrom_lens};

    #[test]
    fn test_context() -> Result<(), anyhow::Error> {
        const FASTA_FILE: &[u8] =
            b">one\nATGCATGCATGCATGCATGCATGCATGCAT\nGCATGCATGCATGCATGCATGCATGCATGC\nATGCAT";
        const FAI_FILE: &[u8] = b"one\t66\t5\t30\t31";

        let mut genome = IndexedReader::new(Cursor::new(FASTA_FILE), FAI_FILE)?;
        let chrom_lens = chrom_lens(&genome);
        assert_eq!(chrom_lens["one"], 66);

        // Truncated start
        let read = Eventalign::empty(String::new(), "one".to_string(), 0, 10, String::new());
        assert_eq!(read.length(), 10);
        assert_eq!(read.start_zb(), 0);
        assert_eq!(read.stop_zb(), 9);
        assert_eq!(read.start_ob(), 1);
        assert_eq!(read.stop_ob(), 10);
        let ctxt = Context::from_read(&mut genome, &chrom_lens, &read)?;
        assert_eq!(ctxt.context.len(), 15);
        // assert_eq!(
        //     std::str::from_utf8(&ctxt.context).unwrap(),
        //     "ATGCATGCATGCATGCATGC"
        // );
        assert_eq!(ctxt.start_slop, 0);
        assert_eq!(ctxt.end_slop, 5);
        assert_eq!(
            ctxt.surrounding(1)
                .into_iter()
                .flat_map(std::str::from_utf8)
                .collect::<Vec<_>>(),
            vec!["ATGCAT", "TGCATG"]
        );

        // Same but minus strand
        let mut read = Eventalign::empty(String::new(), "one".to_string(), 0, 10, String::new());
        *read.strand_mut() = Strand::minus();
        let ctxt = Context::from_read(&mut genome, &chrom_lens, &read)?;
        assert_eq!(ctxt.context.len(), 15);
        assert_eq!(
            std::str::from_utf8(&ctxt.context).unwrap(),
            "TACGTACGTACGTAC"
        );
        assert_eq!(ctxt.start_slop, 0);
        assert_eq!(ctxt.end_slop, 5);
        // assert_eq!(
        //     ctxt.surrounding(1)
        //         .into_iter()
        //         .flat_map(std::str::from_utf8)
        //         .collect::<Vec<_>>(),
        //     vec!["ATGCAT", "TGCATG", "GCATGC", "CATGCA", "ATGCAT", "TGCATG"]
        // );

        // Partial start
        let mut read = Eventalign::empty(String::new(), "one".to_string(), 2, 10, String::new());
        *read.strand_mut() = Strand::minus();
        let ctxt = Context::from_read(&mut genome, &chrom_lens, &read)?;
        assert_eq!(ctxt.context.len(), 17);
        assert_eq!(ctxt.start_slop, 2);
        assert_eq!(ctxt.end_slop, 5);

        // Truncated end
        // Read not possible due to left-adjusted positioning
        let read = Eventalign::empty(String::new(), "one".to_string(), 60, 6, String::new());
        let ctxt = Context::from_read(&mut genome, &chrom_lens, &read)?;
        assert_eq!(ctxt.context.len(), 11);
        assert_eq!(std::str::from_utf8(&ctxt.context).unwrap(), "CATGCATGCAT");

        assert_eq!(ctxt.start_slop, 5);
        assert_eq!(ctxt.end_slop, 0);
        assert_eq!(
            ctxt.surrounding(61)
                .into_iter()
                .flat_map(std::str::from_utf8)
                .collect::<Vec<_>>(),
            vec!["ATGCAT", "TGCATG", "GCATGC", "CATGCA", "ATGCAT"]
        );

        // Partial end
        let read = Eventalign::empty(String::new(), "one".to_string(), 53, 6, String::new());
        let ctxt = Context::from_read(&mut genome, &chrom_lens, &read)?;
        assert_eq!(ctxt.context.len(), 16);
        assert_eq!(ctxt.start_slop, 5);
        assert_eq!(ctxt.end_slop, 2);

        // Middle
        let read = Eventalign::empty(String::new(), "one".to_string(), 30, 10, String::new());
        let ctxt = Context::from_read(&mut genome, &chrom_lens, &read)?;
        assert_eq!(ctxt.context.len(), 20);
        assert_eq!(ctxt.start_slop, 5);
        assert_eq!(ctxt.end_slop, 5);
        Ok(())
    }
}
