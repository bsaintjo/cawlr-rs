use std::{
    collections::HashMap,
    fs::File,
    io::{Read, Seek},
    path::Path,
    str::from_utf8,
};

use anyhow::Result;
use bio::io::fasta::IndexedReader;
use indicatif::ProgressIterator;
use parquet::arrow::ArrowWriter;
use rv::{
    prelude::{Gaussian, Mixture},
    traits::{Rv, KlDivergence},
};
use serde_arrow::{to_record_batch, trace_schema, Schema};

use crate::{
    reads::{FlatLReadScore, LData, LRead, PreprocessRead},
    train::Model,
    utils::CawlrIO,
};

pub(crate) struct ScoreOptions {
    input: String,
    pos_ctrl: Model,
    neg_ctrl: Model,
    genome: IndexedReader<File>,
    chrom_lens: HashMap<String, u64>,
    rank: HashMap<String, f64>,
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
        let kmer_ranks = CawlrIO::load(rank)?;
        let genome = IndexedReader::from_file(&genome)?;
        let chrom_lens = get_genome_chrom_lens(&genome);
        let pos_ctrl_db = CawlrIO::load(&pos_ctrl)?;
        let neg_ctrl_db = CawlrIO::load(&neg_ctrl)?;
        Ok(ScoreOptions {
            input: input.to_owned(),
            pos_ctrl: pos_ctrl_db,
            neg_ctrl: neg_ctrl_db,
            genome,
            chrom_lens,
            rank: kmer_ranks,
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
        self.score(nprs)?;
        Ok(())
    }

    fn score_data(
        &mut self,
        chrom: &str,
        ld: &LData,
        read_seq: &HashMap<u64, Vec<u8>>,
    ) -> Result<Option<(f64, String)>> {
        // let positions = sixmer_postions(&self.chrom_lens, chrom, ld.pos());
        // log::debug!("positions: {positions:?}");
        // let kmers: Vec<&[u8]> = (positions.0..positions.1)
        //     .flat_map(|x| read_seq.get(&x))
        //     .map(|x| x.as_ref())
        //     .collect();
        // let best_kmer = choose_best_kmer(&self.rank, &kmers);
        // let best_kmer = from_utf8(best_kmer)?.to_owned();
        // log::debug!("best_kmer: {best_kmer}");
        let best_kmer = ld.kmer().to_owned();
        let pos_model = self.pos_ctrl.gmms().get(&best_kmer);
        let neg_model = self.neg_ctrl.gmms().get(&best_kmer);
        match (pos_model, neg_model) {
            (Some(pos_gmm), Some(neg_gmm)) => {
                let signal_score = score_signal(ld.mean(), pos_gmm, neg_gmm);
                let skip_score = score_present(ld.kmer().to_string(), &self.pos_ctrl, &self.neg_ctrl);
                log::debug!("signal score: {signal_score:?}");
                log::debug!("skip score: {skip_score:?}");
                let final_score = signal_score.or(skip_score).map(|x| (x, best_kmer));
                Ok(final_score)
            }
            _ => Ok(None)
        }
    }

    fn score_skipped(
        &mut self,
        pos: u64,
        read_seq: &HashMap<u64, Vec<u8>>,
        strand: &Strand,
    ) -> Result<Option<(f64, String)>> {
        let kmer = read_seq[&pos].clone();

        let kmer = match strand {
            Strand::Plus => kmer,
            Strand::Minus => bio::alphabets::dna::revcomp(kmer),
        };

        let kmer = from_utf8(&kmer)?.to_owned();

        let pos_skip = self.pos_ctrl.skips().get(&kmer).map(|x| 1. - x);
        let neg_skip = self.neg_ctrl.skips().get(&kmer).map(|x| 1. - x);
        match (pos_skip, neg_skip) {
            (Some(p), Some(n)) => Ok(Some((p / (p + n), kmer))),
            _ => Ok(None),
        }
    }

    pub(crate) fn score(&mut self, nprs: Vec<PreprocessRead>) -> Result<()> {
        log::debug!("Len nprs {}", nprs.len());
        for npr in nprs.into_iter().progress() {
            let strand = infer_strand(&npr, &mut self.genome)?;
            let chrom = npr.chrom().to_owned();
            let mut acc = Vec::new();
            log::debug!("Read start: {}", npr.start());
            log::debug!("Read stop: {}", npr.stop());
            let read_seq = context_pos(&self.chrom_lens, &mut self.genome, &npr)?;
            let data_pos = pos_with_data(&npr);
            for pos in npr.start()..=npr.stop() {
                let final_score = {
                    if let Some(ld) = data_pos.get(&pos) {
                        self.score_data(&chrom, ld, &read_seq)?
                    } else {
                        self.score_skipped(pos as u64, &read_seq, &strand)?
                    }
                };
                if let Some((score, kmer)) = final_score {
                    let snpr = FlatLReadScore::new(
                        npr.name(),
                        npr.chrom(),
                        npr.start(),
                        npr.length(),
                        &[],
                        pos as u64,
                        score,
                        &kmer,
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

fn pos_with_data(npr: &LRead<LData>) -> HashMap<usize, &LData> {
    let mut avail_pos = HashMap::new();
    for ld in npr.data().iter() {
        avail_pos.insert(ld.pos() as usize, ld);
    }
    avail_pos
}

fn sixmer_postions(chrom_lens: &HashMap<String, u64>, chrom: &str, pos: u64) -> (u64, u64) {
    let chrom_len = chrom_lens.get(chrom).unwrap();
    let start_pos = if pos < 5 { 0 } else { pos - 5 };
    let stop = if (pos + 6) > *chrom_len {
        *chrom_len
    } else {
        pos + 6
    };
    (start_pos, stop + 1)
}

fn context_pos<R>(
    chrom_lens: &HashMap<String, u64>,
    genome: &mut IndexedReader<R>,
    read: &LRead<LData>,
) -> Result<HashMap<u64, Vec<u8>>>
where
    R: Read + Seek,
{
    let chrom_len = chrom_lens.get(read.chrom()).unwrap();
    log::debug!("chrom length: {chrom_len}");
    let start_pos = if read.start() < 5 {
        0
    } else {
        read.start() - 5
    } as u64;
    let stop = read.stop() as u64;
    let stop = if (stop + 13) > *chrom_len {
        *chrom_len
    } else {
        // Need to extend it again so we can get windows past the end of the read
        // (6 for last overlapping kmer) + (6 include sequence of last kmer) + (1 fetch
        // is exclusive for ending)
        stop + 13
    };
    log::debug!("Calc start: {}", start_pos);
    log::debug!("Calc stop: {}", stop);
    genome.fetch(read.chrom(), start_pos, stop)?;
    let mut seq = Vec::new();
    genome.read(&mut seq)?;

    // Will not iterate if (stop - 6) > start_pos
    if (stop - 13) < start_pos {
        log::warn!(
            "Stop is greater than start, no sequence will be produced {:?}",
            read.name()
        );
    }

    let mut seq_map: HashMap<u64, Vec<u8>> = HashMap::new();
    for (kmer, pos) in seq.windows(6).zip(start_pos..) {
        seq_map.insert(pos, kmer.to_owned());
    }
    Ok(seq_map)
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

fn choose_best_kmer<'a>(kmer_ranks: &HashMap<String, f64>, kmers: &[&'a [u8]]) -> &'a [u8] {
    kmers
        .iter()
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

pub(crate) fn choose_model(neg_mix: &Mixture<Gaussian>) -> &Gaussian {
    let true_neg = rv::misc::argmax(neg_mix.weights());
    let true_neg = true_neg[0];
    let true_neg = &neg_mix.components()[true_neg];
    true_neg
}

pub(crate) fn choose_pos_model<'a>(neg_comp: &Gaussian, pos_mix: &'a Mixture<Gaussian>) -> &'a Gaussian {
    pos_mix.components().iter().reduce(|g1, g2| {
        let g1_kl = g1.kl(neg_comp);
        let g2_kl = g2.kl(neg_comp);
        if g1_kl > g2_kl {
            g1
        } else {
            g2
        }
    }).unwrap()
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
    let neg_mix = choose_model(neg_mix);
    let pos_log_proba = pos_mix.f(&signal);
    let neg_log_proba = neg_mix.f(&signal);
    let score = pos_log_proba / (pos_log_proba + neg_log_proba);

    if (pos_mix.ln_f(&signal) < -20.) && (neg_mix.ln_f(&signal) < -20.) {
        None
    } else {
        Some(score)
    }
}

fn score_present(kmer: String, pos_model: &Model, neg_model: &Model) -> Option<f64> {
    let pos_frac = pos_model.skips().get(&kmer);
    let neg_frac = neg_model.skips().get(&kmer);
    match (pos_frac, neg_frac) {
        (Some(p), Some(n)) => Some(p / (p + n)),
        _ => None,
    }
}

pub(crate) enum Strand {
    Plus,
    Minus,
}

impl Strand {
    pub(crate) fn is_plus_strand(&self) -> bool {
        matches!(self, Strand::Plus)
    }
}

// Nanopolish and everything else is zero-based
pub(crate) fn infer_strand<R>(lread: &LRead<LData>, genome: &mut IndexedReader<R>) -> Result<Strand>
where
    R: Read + Seek,
{
    let ld = &lread.data()[0];
    let chrom = lread.chrom();
    let pos = ld.pos();
    let kmer = ld.kmer();

    genome.fetch(chrom, pos, pos + 6)?;
    let mut seq = Vec::new();
    genome.read(&mut seq)?;

    let seq = from_utf8(&seq)?;
    if seq == kmer {
        Ok(Strand::Plus)
    } else {
        Ok(Strand::Minus)
    }
}

#[cfg(test)]
mod test {
    use std::io::Cursor;

    use super::*;

    #[test]
    fn test_infer_strand() {
        const FASTA_FILE: &[u8] = b">chr1\nGTAGGCTGAAAA\nCCCC";
        const FAI_FILE: &[u8] = b"chr1\t16\t6\t12\t13";

        let chrom = "chr1";
        let mut genome = IndexedReader::new(Cursor::new(FASTA_FILE), FAI_FILE).unwrap();
        let ld = LData::new(6, "TGAAAA".to_owned(), 0.0, 0.0);
        let lread: LRead<LData> =
            LRead::new(Vec::new(), chrom.to_owned(), 1, 20, Vec::new(), vec![ld]);
        let ld = &lread.data()[0];
        let chrom = lread.chrom();
        let pos = ld.pos();
        let kmer = ld.kmer();

        genome.fetch(chrom, pos, pos + 6).unwrap();
        let mut seq = Vec::new();
        genome.read(&mut seq).unwrap();

        let seq = from_utf8(&seq).unwrap();
        assert_eq!(kmer, "TGAAAA");
        assert_eq!(seq, "TGAAAA");
    }
}
