use std::{
    collections::HashMap, fmt::Debug, fs::File, hash::BuildHasher, ops::RangeInclusive, path::Path,
};

use arrow2::io::ipc::write::FileWriter;
use bio::io::fasta::IndexedReader;
use eyre::Result;
use fnv::FnvHashMap;
use rv::{
    prelude::{Gaussian, Mixture},
    traits::{Cdf, KlDivergence, Rv},
};
use statrs::statistics::Statistics;

use crate::{
    arrow::{Eventalign, MetadataExt, Score, ScoredRead, Signal},
    context, load_apply,
    motif::{all_bases, Motif},
    save,
    train::{Model, ModelDB},
    utils::{chrom_lens, CawlrIO},
    wrap_writer,
};

pub struct ScoreOptions {
    pos_ctrl: Model,
    neg_ctrl: Model,
    genome: IndexedReader<File>,
    chrom_lens: FnvHashMap<String, u64>,
    rank: FnvHashMap<String, f64>,
    writer: FileWriter<File>,
    cutoff: f64,
    p_value_threshold: f64,
    motifs: Vec<Motif>,
}

impl ScoreOptions {
    pub fn try_new<P>(
        pos_ctrl_filepath: P,
        neg_ctrl_filepath: P,
        genome_filepath: P,
        rank_filepath: P,
        output: P,
    ) -> Result<Self>
    where
        P: AsRef<Path> + Debug,
    {
        let schema = ScoredRead::schema();
        let writer = File::create(output)?;
        let writer = wrap_writer(writer, &schema)?;
        let kmer_ranks = FnvHashMap::load(rank_filepath)?;
        let genome = IndexedReader::from_file(&genome_filepath)
            .map_err(|_| eyre::eyre!("Failed to read genome file"))?;
        let chrom_lens = chrom_lens(&genome);
        let pos_ctrl_db = Model::load(&pos_ctrl_filepath)?;
        let neg_ctrl_db = Model::load(&neg_ctrl_filepath)?;
        Ok(ScoreOptions {
            pos_ctrl: pos_ctrl_db,
            neg_ctrl: neg_ctrl_db,
            genome,
            chrom_lens,
            rank: kmer_ranks,
            writer,
            cutoff: 10.0,
            p_value_threshold: 0.05,
            motifs: all_bases(),
        })
    }

    pub fn cutoff(&mut self, cutoff: f64) -> &mut Self {
        self.cutoff = cutoff;
        self
    }

    pub fn p_value_threshold(&mut self, p_value_threshold: f64) -> &mut Self {
        self.p_value_threshold = p_value_threshold;
        self
    }

    pub fn motifs<V: Into<Vec<Motif>>>(&mut self, motifs: V) -> &mut Self {
        self.motifs = motifs.into();
        self
    }

    fn close(mut self) -> Result<()> {
        self.writer.finish()?;
        Ok(())
    }

    /// For every read in the input file, try to calculate scores for each base
    /// position and write to file.
    pub fn run<P>(mut self, input: P) -> Result<()>
    where
        P: AsRef<Path>,
    {
        let file = File::open(input)?;
        load_apply(file, |eventaligns| {
            let scored = eventaligns
                .into_iter()
                .flat_map(|e| self.score_eventalign(e))
                .collect();
            self.save(scored)
        })?;
        self.close()
    }

    /// Write batch of scored reads to the writer.
    pub(crate) fn save(&mut self, scored: Vec<ScoredRead>) -> Result<()> {
        save(&mut self.writer, &scored)
    }

    /// Scores a single Eventalign read. For each read, loop over each base pair
    /// position, and if the kmer at the position matches the motif attempt to
    /// score it.
    fn score_eventalign(&mut self, read: Eventalign) -> Result<ScoredRead> {
        let mut acc = Vec::new();
        let context = context::Context::from_read(&mut self.genome, &self.chrom_lens, &read)?;

        log::debug!("{:?}", read.metadata());
        log::debug!("{context:.3?}");

        let data_pos = pos_with_data(&read);
        for pos in read.start_1b()..read.end_1b_excl() {
            // Get kmer and check if kmer matches the motifs, if there are any supplied
            let pos_kmer: Option<(&[u8], &Motif)> = context.sixmer_at(pos).and_then(|k| {
                self.motifs
                    .iter()
                    .find(|m| {
                        let m = m.motif().as_bytes();
                        k.starts_with(m)
                    })
                    .map(|m| (k, m))
            });

            if let Some((kmer, motif)) = pos_kmer {
                let kmer = std::str::from_utf8(kmer).unwrap().to_string();
                log::debug!("Position {pos} kmer: {kmer}");

                let signal_score = self.calc_signal_score(pos, &data_pos);
                let skipping_score = self.calc_skipping_score(pos, &data_pos, &context, motif)?;
                let final_score = signal_score.map_or(skipping_score, |x| x.max(skipping_score));
                let score = Score::new(
                    pos,
                    kmer,
                    signal_score.is_none(),
                    signal_score,
                    skipping_score,
                    final_score,
                );
                log::debug!("final score: {score:.3?}");
                acc.push(score)
            }
        }
        let scored_read = ScoredRead::from_read_with_scores(read, acc);
        Ok(scored_read)
    }

    fn calc_skipping_score(
        &self,
        pos: u64,
        data_pos: &FnvHashMap<u64, &Signal>,
        context: &context::Context,
        motif: &Motif,
    ) -> Result<f64> {
        let sur_kmers = context.surrounding(pos, motif);
        let sur_has_data = surround_has_data(pos, data_pos);
        let skipping_scores = sur_kmers
            .into_iter()
            .zip(sur_has_data.into_iter())
            .flat_map(|(kmer, has_data)| {
                let kmer = std::str::from_utf8(kmer).expect("Invalid kmer");
                let pos_presence = self.pos_ctrl.skips().get(kmer);
                let neg_presence = self.neg_ctrl.skips().get(kmer);
                match (pos_presence, neg_presence) {
                    (Some(&pos_presence), Some(&neg_presence)) => {
                        if has_data {
                            Some(pos_presence / (pos_presence + neg_presence))
                        } else {
                            let pos_absent = 1. - pos_presence;
                            let neg_absent = 1. - neg_presence;
                            Some(pos_absent / (pos_absent + neg_absent))
                        }
                    }
                    _ => None,
                }
            })
            .collect::<Vec<_>>();

        // TODO: Switch to median when it can be correctly handled
        let skip_score = skipping_scores.mean();
        if skip_score.is_nan() {
            Err(eyre::eyre!("No data for calculating median"))
        } else {
            Ok(skip_score)
        }
    }

    /// For a given position, get the values for the position and surrounding
    /// kmers. Filter for the best kmer model, if there is confidence in the
    /// model, otherwise return None.
    fn calc_signal_score(&self, pos: u64, data_pos: &FnvHashMap<u64, &Signal>) -> Option<f64> {
        log::debug!("Calculating signal score");
        let sur_signals = surrounding_signal(pos, data_pos);
        log::debug!("surrounding signals: {sur_signals:.3?}");
        let best_signal = best_surrounding_signal(
            sur_signals,
            &self.rank,
            self.pos_ctrl.gmms(),
            self.neg_ctrl.gmms(),
            self.p_value_threshold,
        );

        log::debug!("Best signal: {best_signal:.3?}");

        best_signal.and_then(|sig| {
            let mean = sig.mean();
            let kmer = sig.kmer();
            let pos_mix = self.pos_ctrl.gmms().get(kmer);
            let neg_mix = self.neg_ctrl.gmms().get(kmer);
            match (pos_mix, neg_mix) {
                (Some(pos_gmm), Some(neg_gmm)) => {
                    let neg_mix = neg_gmm.mixture();
                    let pos_mix = pos_gmm.mixture();
                    score_signal(mean, &pos_mix, &neg_mix, self.cutoff)
                }
                _ => {
                    log::debug!("Missing kmer, unable to score signal.");
                    None
                }
            }
        })
    }
}

fn surrounding_pos(pos: u64) -> RangeInclusive<u64> {
    let start = if pos < 5 { 0 } else { pos - 5 };
    start..=pos
}

/// Return list of kmer positions around a given position pos contain signal
/// current data
fn surround_has_data<S>(pos: u64, signal_map: &HashMap<u64, &Signal, S>) -> Vec<bool>
where
    S: BuildHasher,
{
    let positions = surrounding_pos(pos);
    positions.map(|p| signal_map.get(&p).is_some()).collect()
}

/// Returns None if none of the positions around a genomic position have signal
/// measurements.
fn surrounding_signal<'a, S>(
    pos: u64,
    signal_map: &HashMap<u64, &'a Signal, S>,
) -> Option<Vec<&'a Signal>>
where
    S: BuildHasher,
{
    let positions = surrounding_pos(pos);
    let acc = positions
        .flat_map(|p| signal_map.get(&p))
        .cloned()
        .collect::<Vec<_>>();
    if acc.is_empty() {
        None
    } else {
        Some(acc)
    }
}

/// Return mu and sigma from a Gaussian distribution.
fn extract_components(gauss: &Gaussian) -> (f64, f64) {
    let mu = gauss.mu();
    let sigma = gauss.sigma();
    (mu, sigma)
}

/// Use z-test to calculate the p-value between two Gaussians
fn gauss_to_pvalue(pos_model: &Gaussian, neg_model: &Gaussian) -> f64 {
    let (pos_mu, pos_sigma) = extract_components(pos_model);
    let (neg_mu, neg_sigma) = extract_components(neg_model);

    let zscore = (pos_mu - neg_mu) / ((neg_sigma.powi(2)) + pos_sigma.powi(2)).sqrt();
    zscore_to_tt_pvalue(zscore)

    // 2. * Gaussian::standard().sf(&zscore)
}

fn zscore_to_tt_pvalue(zscore: f64) -> f64 {
    2. * Gaussian::standard().sf(&zscore.abs())
}

/// Filters out surrounding signal for best signal to use for scoring.
/// Will return None if one of the signal's kmers have a z-test p-value less
/// than 0.05.
fn best_surrounding_signal<'a, S>(
    surrounding: Option<Vec<&'a Signal>>,
    ranks: &HashMap<String, f64, S>,
    pos_gmms: &ModelDB,
    neg_gmms: &ModelDB,
    p_value_threshold: f64,
) -> Option<&'a Signal>
where
    S: BuildHasher,
{
    log::debug!("Determine best surrounding signal");
    surrounding.and_then(|signals| {
        signals
            .into_iter()
            // Only use kmers with z-test p-values less than 0.05
            .filter(|&s| {
                log::debug!("Signal: {s:.3?}");
                let kmer = s.kmer();
                if !neg_gmms.contains_key(kmer) || !pos_gmms.contains_key(kmer) {
                    false
                } else {
                    let neg_mix = neg_gmms[kmer].mixture();
                    let pos_mix = pos_gmms[kmer].mixture();
                    let neg_model = choose_model(&neg_mix);
                    let pos_model = choose_pos_model(neg_model, &pos_mix);
                    let pvalue = gauss_to_pvalue(pos_model, neg_model);
                    log::debug!("p-value: {pvalue:.3?}");
                    pvalue < p_value_threshold
                }
            })
            // Of the ones the best, choose the one with the best ranking
            .reduce(|x, y| {
                let x_rank = ranks.get(x.kmer());
                let y_rank = ranks.get(y.kmer());
                match (x_rank, y_rank) {
                    (None, _) => y,
                    (_, None) => x,
                    (Some(a), Some(b)) => {
                        if a > b {
                            x
                        } else {
                            y
                        }
                    }
                }
            })
    })
}

/// Returns HashMap mapping positions as u64 to the respective signal data
/// Useful for iterating through each base pair position and computing results
/// based on if there is data or not
fn pos_with_data(read: &Eventalign) -> FnvHashMap<u64, &Signal> {
    let mut avail_pos = FnvHashMap::default();
    for signal in read.signal_iter() {
        avail_pos.insert(signal.pos(), signal);
    }
    avail_pos
}

/// Return the Gaussian with the highest component weight. This is a heuristic
/// that expects that the highest weight component in the negative control
/// should represent the data from the true negative control distribution.
pub(crate) fn choose_model(neg_mix: &Mixture<Gaussian>) -> &Gaussian {
    let true_neg = rv::misc::argmax(neg_mix.weights());
    let true_neg = true_neg[0];
    let true_neg = &neg_mix.components()[true_neg];
    true_neg
}

/// Given a Gaussian, and a mixture model containing two gaussians, find the
/// Gaussian in the mixture model that is most disimilar based on KL divergence
/// and return it.
/// Should not fail because train should always produce a Mixture model
/// containing two gaussians
pub(crate) fn choose_pos_model<'a>(
    neg_comp: &Gaussian,
    pos_mix: &'a Mixture<Gaussian>,
) -> &'a Gaussian {
    pos_mix
        .components()
        .iter()
        .reduce(|g1, g2| {
            let g1_kl = g1.kl(neg_comp);
            let g2_kl = g2.kl(neg_comp);
            if g1_kl > g2_kl {
                g1
            } else {
                g2
            }
        })
        .unwrap()
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
    cutoff: f64,
) -> Option<f64> {
    log::debug!("Scoring signal: {signal}");
    let neg_mix = choose_model(neg_mix);
    let pos_mix = choose_pos_model(neg_mix, pos_mix);
    let pos_proba = pos_mix.f(&signal);
    let neg_proba = neg_mix.f(&signal);
    let score = pos_proba / (pos_proba + neg_proba);
    log::debug!("Score: {score:.3}");

    let pos_log_proba = pos_mix.ln_f(&signal);
    let neg_log_proba = neg_mix.ln_f(&signal);

    log::debug!("+ Gaussian log proba: {pos_log_proba}");
    log::debug!("- Gaussian log proba: {neg_log_proba}");

    if (pos_log_proba > -cutoff) || (neg_log_proba > -cutoff) {
        log::debug!("Valid score");
        Some(score)
    } else {
        log::debug!("Below cutoff, not scoring.");
        None
    }
}

#[cfg(test)]
mod test {
    use assert_fs::TempDir;
    use float_eq::assert_float_eq;

    use super::*;
    use crate::{collapse::CollapseOptions, motif::Motif, arrow_utils::load_iter};

    #[test]
    fn test_score_signal() {
        let signal = 80.0;
        let cutoff = 10.0;

        let neg_mix = Mixture::new(
            vec![0.9, 0.1],
            vec![
                Gaussian::new(100.0, 1.0).unwrap(),
                Gaussian::new(100.0, 1.0).unwrap(),
            ],
        )
        .unwrap();
        let pos_mix = Mixture::new(
            vec![0.9, 0.1],
            vec![
                Gaussian::new(80.0, 1.0).unwrap(),
                Gaussian::new(100.0, 1.0).unwrap(),
            ],
        )
        .unwrap();

        let result = score_signal(signal, &pos_mix, &neg_mix, cutoff);
        assert!(result.is_some());

        let result = score_signal(1000.0, &pos_mix, &neg_mix, cutoff);
        assert!(result.is_none());
    }

    #[test]
    fn test_zscore_to_tt_pvalue() {
        assert_float_eq!(zscore_to_tt_pvalue(2.9), 0.003_732, abs <= 0.000_001);
        assert_float_eq!(zscore_to_tt_pvalue(0.1), 0.920_344, abs <= 0.000_001);
        assert_float_eq!(zscore_to_tt_pvalue(-0.3), 0.764_177, abs <= 0.000_001);
        assert_float_eq!(zscore_to_tt_pvalue(-3.2), 0.001_374, abs <= 0.000_001);
        assert_float_eq!(zscore_to_tt_pvalue(0.0), 1.0, abs <= 0.000_001);

        // These won't panic
        zscore_to_tt_pvalue(f64::NEG_INFINITY);
        zscore_to_tt_pvalue(f64::NAN);
        zscore_to_tt_pvalue(f64::INFINITY);
    }

    #[test]
    fn test_single_read() -> Result<()> {
        let temp_dir = TempDir::new()?;
        let filepath = "extra/single_read.eventalign.txt";
        let input = File::open(filepath)?;
        let bam_file = "extra/single_read.bam";
        let output = temp_dir.path().join("test");
        let mut collapse = CollapseOptions::try_new(bam_file, &output)?;
        collapse.run(input)?;

        let output = File::open(output)?;
        let reads = load_iter(output).next().unwrap().unwrap();
        let read = &reads[0];

        let genome_file = "extra/sacCer3.fa";
        let mut genome = IndexedReader::from_file(&genome_file)
            .map_err(|_| eyre::eyre!("Failed to read genome file."))?;

        let chrom_lens = chrom_lens(&genome);

        let context = context::Context::from_read(&mut genome, &chrom_lens, read)?;
        assert_eq!(context.start_slop(), 5);
        // assert_eq!(context.end_slop(), 5);

        let m = Motif::new("AT", 2);
        assert_eq!(m.position_0b(), 1);
        assert_eq!(
            context
                .surrounding(182522, &m)
                .into_iter()
                .flat_map(std::str::from_utf8)
                .collect::<Vec<_>>(),
            vec!["ACATAT", "CATATT", "ATATTC", "TATTCA", "ATTCAA", "TTCAAT"]
        );

        Ok(())
    }
}
