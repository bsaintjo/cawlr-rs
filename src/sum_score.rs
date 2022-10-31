use std::{
    io::{Read, Seek, Write},
    path::Path,
};

use eyre::Result;
use fnv::FnvHashMap;

use crate::{
    arrow::{Signal},
    motif::Motif,
    train::Model,
    Eventalign, Score, ScoredRead, arrow_utils::{load_read_write_arrow, SchemaExt},
};

#[derive(Debug)]
struct ScoreOptions {
    pos_model: Model,
    neg_model: Model,
    ranks: FnvHashMap<String, f64>,
    freq_thresh: usize,
    motifs: Vec<Motif>,
}

fn count_motif_in_kmer(kmer: &str, motif: &Motif) -> usize {
    kmer.matches(motif.motif()).count()
}

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
    fn new(
        pos_model: Model,
        neg_model: Model,
        ranks: FnvHashMap<String, f64>,
        freq_thresh: usize,
        motifs: Vec<Motif>,
    ) -> Self {
        Self {
            pos_model,
            neg_model,
            ranks,
            freq_thresh,
            motifs,
        }
    }

    fn load<P>(pos_model_filepath: P, neg_model_filepath: P, ranks_filepath: P) -> Result<Self>
    where
        P: AsRef<Path>,
    {
        todo!()
    }

    fn freq_thresh(&mut self, freq_thresh: usize) -> &mut Self {
        self.freq_thresh = freq_thresh;
        self
    }

    fn motifs(&mut self, motifs: Vec<Motif>) -> &mut Self {
        self.motifs = motifs;
        self
    }

    fn run<R: Read + Seek, W: Write>(&self, reader: R, writer: W) -> Result<()> {
        let writer = ScoredRead::wrap_writer(writer)?;
        load_read_write_arrow(reader, writer, |eventaligns: Vec<Eventalign>| {
            let mut scored_reads = Vec::new();
            for eventalign in eventaligns {
                let mut scores = Vec::new();
                let data_map = eventalign
                    .signal_iter()
                    .map(|s| (s.pos(), s))
                    .collect::<FnvHashMap<_, _>>();
                for signal in eventalign.signal_iter() {
                    let kmer = signal.kmer();
                    if let Some(m) = self.motifs.iter().find(|m| kmer.starts_with(m.motif())) {
                        let mut kmers = Vec::new();
                        let surround_idx = signal.pos() + m.position_0b() as u64;
                        let surrounding =
                            surround_idx - (5 + m.position_0b() as u64)..=surround_idx;
                        for surr in surrounding {
                            if let Some(&s) = data_map.get(&surr) {
                                if signal.samples().len() > self.freq_thresh {
                                    continue;
                                }

                                let kmer = s.kmer();
                                if count_motif_in_kmer(kmer, m) > 1 {
                                    continue;
                                }
                                let pm = self.pos_model.gmms().get(kmer);
                                let nm = self.neg_model.gmms().get(kmer);
                                if let (Some(pm), Some(nm)) = (pm, nm) {
                                    let pos_model = pm.mixture();
                                    let neg_model = nm.single();

                                    let (pos_sum, neg_sum) = s.score_lnsum(&pos_model, &neg_model);
                                    kmers.push(SignalScore::new(s, pos_sum, neg_sum));
                                }
                            }
                        }
                        let mut best_signal = None;
                        let mut diff = f64::NEG_INFINITY;
                        for ss in kmers.into_iter() {
                            if let Some(&rank) = self.ranks.get(ss.signal.kmer()) {
                                if rank > diff {
                                    diff = rank;
                                    best_signal = Some(ss);
                                }
                            }
                        }

                        if let Some(best_signal) = best_signal {
                            let exp_me = best_signal.pos_sum.exp();
                            let exp_un = best_signal.neg_sum.exp();
                            let rate = exp_me / (exp_me + exp_un);
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
