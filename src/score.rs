use std::{
    collections::HashMap,
    fs::File,
    io::{Read, Seek},
    path::Path,
    str::from_utf8,
};

use anyhow::Result;
use bio::io::fasta::IndexedReader;
use parquet::arrow::ArrowWriter;
use rv::{
    prelude::{Gaussian, Mixture},
    traits::Rv,
};
use serde_arrow::{to_record_batch, trace_schema, Schema};

use crate::{
    reads::{FlatLReadScore, LData, LRead, PreprocessRead},
    train::Model,
    utils::CawlrIO,
};

pub(crate) struct ScoreOptions {
    input: String,
    pos_ctrl: String,
    neg_ctrl: String,
    genome: String,
    rank: String,
    writer: ArrowWriter<File>,
    schema: Schema,
}

impl ScoreOptions {
    pub(crate) fn try_new<P>(
        input: &str,
        pos_ctrl: &str,
        neg_ctrl: &str,
        genome: &str,
        rank: &str,
        output: P,
    ) -> Result<Self>
    where
        P: AsRef<Path>,
    {
        let example = vec![FlatLReadScore::default()];
        let schema = trace_schema(&example)?;
        let batches = to_record_batch(&example, &schema)?;
        let writer = File::create(output)?;
        let writer = ArrowWriter::try_new(writer, batches.schema(), None)?;
        Ok(ScoreOptions {
            input: input.to_owned(),
            pos_ctrl: pos_ctrl.to_owned(),
            neg_ctrl: neg_ctrl.to_owned(),
            genome: genome.to_owned(),
            rank: rank.to_owned(),
            writer,
            schema,
        })
    }

    pub(crate) fn close(mut self) -> Result<()> {
        self.writer.close()?;
        Ok(())
    }

    pub(crate) fn run(&mut self) -> Result<()> {
        let nprs: Vec<LRead<LData>> = CawlrIO::load(&self.input)?;
        let pos_ctrl_db = CawlrIO::load(&self.pos_ctrl)?;
        let neg_ctrl_db = CawlrIO::load(&self.neg_ctrl)?;
        let kmer_ranks = CawlrIO::load(&self.rank)?;
        let genome = IndexedReader::from_file(&self.genome)?;
        self.score(nprs, pos_ctrl_db, neg_ctrl_db, kmer_ranks, genome)?;
        Ok(())
    }

    pub(crate) fn score<R>(
        &mut self,
        nprs: Vec<PreprocessRead>,
        pos_models: Model,
        neg_models: Model,
        kmer_ranks: HashMap<String, f64>,
        mut genome: IndexedReader<R>,
    ) -> Result<()>
    where
        R: Read + Seek,
    {
        let chrom_lens = get_genome_chrom_lens(&genome);
        log::debug!("Len nprs {}", nprs.len());
        for npr in nprs.into_iter() {
            let chrom = npr.chrom().to_owned();
            let mut acc = Vec::new();
            for data in npr.data().iter() {
                let ctxt = get_genomic_context(&chrom_lens, &mut genome, &chrom, data.pos())?;
                let best_kmer = choose_best_kmer(&kmer_ranks, &ctxt);
                let best_kmer = from_utf8(best_kmer)?;
                log::debug!("best_kmer: {best_kmer}");
                let pos_model = pos_models.gmms().get(best_kmer).unwrap();
                let neg_model = neg_models.gmms().get(best_kmer).unwrap();
                let signal_score = score_signal(data.mean(), pos_model, neg_model);
                let skip_score = score_skip(data.kmer().to_string(), &pos_models, &neg_models);
                log::debug!("signal score: {signal_score:?}");
                log::debug!("skip score: {skip_score:?}");
                let final_score = signal_score.or(skip_score);
                log::debug!("Final score: {final_score:?}");
                if let Some(score) = final_score {
                    let snpr = FlatLReadScore::new(
                        npr.name(),
                        npr.chrom(),
                        npr.start(),
                        npr.length(),
                        &[],
                        data.pos(),
                        score,
                    );
                    acc.push(snpr);
                }
            }
            if !acc.is_empty() {
                self.save_flatscores(&acc)?;
                acc.clear();
            }
        }
        Ok(())
    }

    fn save_flatscores(&mut self, flats: &[FlatLReadScore]) -> Result<()> {
        let batches = to_record_batch(flats, &self.schema)?;
        self.writer.write(&batches)?;
        Ok(())
    }
}

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
        .map(|x| {
            let x_str = from_utf8(x).unwrap();
            (x, kmer_ranks.get(x_str))
        })
        .filter(|x| x.1.is_some())
        .reduce(|a, b| match (a.1, b.1) {
            (None, _) => b,
            (_, None) => a,
            (Some(x), Some(y)) => {
                if x > y {
                    a
                } else {
                    b
                }
            }
        })
        .expect("Genomic context is empty.")
        .0
}

/// Score given signal based on GMM from a positive and negative control.
/// Scoring function based on:
///  Wang, Y. et al. Single-molecule long-read sequencing reveals the chromatin
/// basis of gene expression. Genome Res. 29, 1329–1342 (2019).
/// We don't take the ln(score) for now, only after the probability from the Kde
/// later in cawlr sma
fn score_signal(
    signal: f64,
    pos_mix: &Mixture<Gaussian>,
    neg_mix: &Mixture<Gaussian>,
) -> Option<f64> {
    let pos_log_proba = pos_mix.f(&signal);
    let neg_log_proba = neg_mix.f(&signal);
    let score = pos_log_proba / (pos_log_proba + neg_log_proba);

    if (pos_mix.ln_f(&signal) < -20.) && (neg_mix.ln_f(&signal) < -20.) {
        None
    } else {
        Some(score)
    }
}

fn score_skip(kmer: String, pos_model: &Model, neg_model: &Model) -> Option<f64> {
    let pos_frac = pos_model.skips().get(&kmer);
    let neg_frac = neg_model.skips().get(&kmer);
    match (pos_frac, neg_frac) {
        (Some(&p), Some(&n)) => Some(p / (p + n)),
        _ => None,
    }
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
