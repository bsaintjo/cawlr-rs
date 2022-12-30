use std::{
    io::{Read, Seek, Write},
    path::Path,
};

use eyre::Result;
use fnv::FnvHashMap;

use crate::{
    arrow::Signal,
    arrow_utils::load_read_write_arrow,
    motif::{all_bases, Motif},
    train::Model,
    utils::CawlrIO,
    Eventalign, Score, ScoredRead,
};

pub struct ScoreOptions {
    pos_model: Model,
    neg_model: Model,
    ranks: FnvHashMap<String, f64>,
    freq_thresh: usize,
    cutoff: f64,
    motifs: Vec<Motif>,
}

impl std::fmt::Debug for ScoreOptions {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("ScoreOptions")
            .field("freq_thresh", &self.freq_thresh)
            .field("cutoff", &self.cutoff)
            .field("motifs", &self.motifs)
            .finish_non_exhaustive()
    }
}

fn count_motif_in_kmer(kmer: &str, motif: &Motif) -> usize {
    kmer.matches(motif.motif()).count()
}

#[derive(Debug)]
struct SignalScore<'a> {
    signal: &'a Signal,
    pos_sum: f64,
    neg_sum: f64,
}

impl<'a> SignalScore<'a> {
    fn new(signal: &'a Signal, pos_sum: f64, neg_sum: f64) -> Self {
        Self {
            signal,
            pos_sum,
            neg_sum,
        }
    }
}

impl ScoreOptions {
    pub fn new(
        pos_model: Model,
        neg_model: Model,
        ranks: FnvHashMap<String, f64>,
        freq_thresh: usize,
        cutoff: f64,
        motifs: Vec<Motif>,
    ) -> Self {
        Self {
            pos_model,
            neg_model,
            ranks,
            freq_thresh,
            cutoff,
            motifs,
        }
    }

    pub fn load<P>(pos_model_filepath: P, neg_model_filepath: P, ranks_filepath: P) -> Result<Self>
    where
        P: AsRef<Path>,
    {
        let pos_model = Model::load(pos_model_filepath)?;
        let neg_model = Model::load(neg_model_filepath)?;
        let ranks = FnvHashMap::load(ranks_filepath)?;
        let score_options = ScoreOptions::new(pos_model, neg_model, ranks, 10, 10.0, all_bases());
        log::debug!("Score Options: {score_options:?}");
        Ok(score_options)
    }

    pub fn freq_thresh(&mut self, freq_thresh: usize) -> &mut Self {
        self.freq_thresh = freq_thresh;
        self
    }

    pub fn cutoff(&mut self, cutoff: f64) -> &mut Self {
        self.cutoff = cutoff;
        self
    }

    pub fn motifs(&mut self, motifs: Vec<Motif>) -> &mut Self {
        self.motifs = motifs;
        self
    }

    pub fn run<R, W>(&self, reader: R, writer: W) -> Result<()>
    where
        R: Read + Seek,
        W: Write,
    {
        load_read_write_arrow(reader, writer, |eventaligns: Vec<Eventalign>| {
            let mut scored_reads = Vec::new();
            for eventalign in eventaligns {
                log::debug!("eventalign: {:?}", eventalign.metadata());
                let mut scores = Vec::new();
                let data_map = eventalign
                    .signal_iter()
                    .map(|s| (s.pos(), s))
                    .collect::<FnvHashMap<_, _>>();
                for signal in eventalign.signal_iter() {
                    log::debug!("signal {signal:?}");
                    let kmer = signal.kmer();
                    if let Some(m) = self.motifs.iter().find(|m| kmer.starts_with(m.motif())) {
                        log::debug!("Kmer motif matches {m:?}");
                        let mut kmers = Vec::new();
                        let surrounding = m.surrounding_idxs(signal.pos());
                        for surr in surrounding {
                            log::debug!("Surrounding idx {surr}");
                            if let Some(&s) = data_map.get(&surr) {
                                log::debug!("Surrounding signal: {s:?}");
                                if signal.samples().len() > self.freq_thresh {
                                    log::debug!(
                                        "n samples greater than frequency threshold, skipping"
                                    );
                                    continue;
                                }

                                let kmer = s.kmer();
                                if count_motif_in_kmer(kmer, m) > 1 {
                                    log::debug!("Count of motifs in kmer greater than 1, skipping");
                                    continue;
                                }
                                let pm = self.pos_model.gmms().get(kmer);
                                let nm = self.neg_model.gmms().get(kmer);
                                if let (Some(pm), Some(nm)) = (pm, nm) {
                                    let pos_model = pm.mixture();
                                    let neg_model = nm.single();

                                    if let Some((pos_sum, neg_sum)) = s.score_lnsum(&pos_model, &neg_model) {
                                        kmers.push(SignalScore::new(s, pos_sum, neg_sum));
                                    }
                                }
                            }
                        }
                        let mut best_signal = None;
                        let mut diff = f64::NEG_INFINITY;
                        for ss in kmers.into_iter() {
                            if let Some(&rank) = self.ranks.get(ss.signal.kmer()) {
                                log::debug!("signal score: {ss:?}");
                                if rank > diff {
                                    diff = rank;
                                    best_signal = Some(ss);
                                }
                            }
                        }

                        if let Some(best_signal) = best_signal {
                            log::debug!("Best signal: {best_signal:?}");

                            let exp_me = best_signal.pos_sum.exp();
                            let exp_un = best_signal.neg_sum.exp();

                            let rate = exp_me / (exp_me + exp_un);

                            log::debug!("exp_me: {exp_me}");
                            log::debug!("exp_un: {exp_un}");
                            log::debug!("rate: {rate}");

                            let score = Score::new(
                                signal.pos(),
                                signal.kmer().to_string(),
                                false,
                                Some(rate),
                                0.0,
                                rate,
                            );
                            scores.push(score);
                        }
                    }
                }
                let scored = ScoredRead::from_read_with_scores(eventalign, scores);
                scored_reads.push(scored);
            }
            Ok(scored_reads)
        })?;
        Ok(())
    }
}
