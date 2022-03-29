use std::{
    collections::HashMap,
    io::{Read, Seek},
    str::from_utf8,
};

use anyhow::Result;
use bio::io::fasta::IndexedReader;
use rv::{
    prelude::{Gaussian, Mixture},
    traits::Rv,
};

use crate::{
    reads::{PreprocessRead, Score, ScoredRead},
    train::Model,
};

// TODO: test for positions after chr_len
// TODO: Deal with +1 issue at chromosome ends
fn get_genomic_context<R>(
    chrom_lens: &HashMap<String, u64>,
    genome: &mut IndexedReader<R>,
    chr: &str,
    pos: u64,
) -> Result<Vec<u8>>
where
    R: Read + Seek,
{
    let chr_len = chrom_lens.get(chr).unwrap();
    let start = if pos < 5 { 0 } else { pos - 5 };
    let stop = if (pos + 6) > *chr_len {
        *chr_len
    } else {
        pos + 6
    };
    genome.fetch(chr, start, stop)?;
    let mut seq = Vec::new();
    genome.read_iter()?.flatten().for_each(|bp| seq.push(bp));
    Ok(seq)
}

fn get_genome_chrom_lens<R>(genome: &IndexedReader<R>) -> HashMap<String, u64>
where
    R: Read + Seek,
{
    let mut chrom_lens = HashMap::new();
    genome.index.sequences().into_iter().for_each(|sequence| {
        chrom_lens.insert(sequence.name, sequence.len);
    });
    chrom_lens
}

// TODO Remove unwrap to allow for using kmers that are above a threshold and
// output is Option/Result
// TODO Deal with possible failures from f64::partial_cmp
fn choose_best_kmer<'a>(kmer_ranks: &HashMap<String, f64>, context: &'a [u8]) -> &'a [u8] {
    context
        .windows(6)
        .max_by(|&x, &y| {
            let x = from_utf8(x).unwrap();
            let y = from_utf8(y).unwrap();

            let x_rank = kmer_ranks.get(x).unwrap();
            let y_rank = kmer_ranks.get(y).unwrap();

            x_rank.partial_cmp(y_rank).unwrap()
        })
        .expect("Genomic context is empty.")
}

/// Score given signal based on GMM from a positive and negative control.
/// Scoring function based on:
///  Wang, Y. et al. Single-molecule long-read sequencing reveals the chromatin
/// basis of gene expression. Genome Res. 29, 1329–1342 (2019).
/// We don't take the ln(score) for now, only after the probability from the Kde
/// later in cawlr sma
///
/// TODO explore use log-likelihood ratio instead of this scoring
fn score_signal(signal: f64, pos_mix: &Mixture<Gaussian>, neg_mix: &Mixture<Gaussian>) -> Option<f64> {
    let pos_log_proba = pos_mix.f(&signal);
    let neg_log_proba = neg_mix.f(&signal);
    let score = pos_log_proba / (pos_log_proba + neg_log_proba);

    Some(score)
}

fn score_skip(kmer: String, pos_model: &Model, neg_model: &Model) -> Option<f64> {
    let pos_frac = pos_model.skips().get(&kmer);
    let neg_frac = neg_model.skips().get(&kmer);
    match (pos_frac, neg_frac) {
        (Some(&p), Some(&n)) => {
            Some(p / (p + n))
        }
        _ => None,
    }
}

// TODO use rayon for parallel scoring
// TODO: add scoring of skipped events
pub(crate) fn score<R>(
    nprs: Vec<PreprocessRead>,
    pos_models: Model,
    neg_models: Model,
    kmer_ranks: HashMap<String, f64>,
    mut genome: IndexedReader<R>,
) -> Vec<ScoredRead>
where
    R: Read + Seek,
{
    let chrom_lens = get_genome_chrom_lens(&genome);
    nprs.into_iter()
        .map(|npr| {
            let chrom = npr.chrom().to_owned();
            npr.map_data(|ld| {
                let ctxt = get_genomic_context(&chrom_lens, &mut genome, &chrom, ld.pos())
                    .expect("Failed to read genome fasta.");
                let best_kmer = choose_best_kmer(&kmer_ranks, &ctxt);
                let best_kmer = from_utf8(best_kmer).unwrap();
                let pos_model = pos_models.gmms().get(best_kmer).unwrap();
                let neg_model = neg_models.gmms().get(best_kmer).unwrap();
                let signal_ll = score_signal(ld.mean(), pos_model, neg_model).unwrap();
                Score::new(ld.pos(), signal_ll)
            })
        })
        .collect()
}

#[cfg(test)]
mod test {
    use std::io::Cursor;

    use super::*;

    #[test]
    fn test_get_genomic_context() {
        const FASTA_FILE: &[u8] = b">chr1\nGTAGGCTGAAAA\nCCCC";
        const FAI_FILE: &[u8] = b"chr1\t16\t6\t12\t13";

        let chrom = "chr1";
        let pos = 7;
        let mut genome = IndexedReader::new(Cursor::new(FASTA_FILE), FAI_FILE).unwrap();
        let chrom_lens = get_genome_chrom_lens(&genome);
        let seq = get_genomic_context(&chrom_lens, &mut genome, chrom, pos).unwrap();
        assert_eq!(seq, b"AGGCTGAAAAC");

        let pos = 2;
        let seq = get_genomic_context(&chrom_lens, &mut genome, chrom, pos).unwrap();
        assert_eq!(seq, b"GTAGGCTG");

        // eprintln!("{:?}", genome.index.sequences());
        let pos = 14;
        let seq = get_genomic_context(&chrom_lens, &mut genome, chrom, pos).unwrap();
        assert_eq!(seq, b"AAACCCC");
    }
}
