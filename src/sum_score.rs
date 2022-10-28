use std::{
    io::{Read, Seek},
};

use eyre::Result;
use fnv::FnvHashMap;

use crate::{arrow::Signal, load_apply, train::Model, Eventalign};

struct ScoreOptions {
    pos_model: Model,
    neg_model: Model,
}

const MOTIF: &str = "GC";
const MOD_START: u64 = 1;

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
    fn run<R: Read + Seek>(&self, reader: R) -> Result<()> {
        load_apply(reader, |eventaligns: Vec<Eventalign>| {
            for eventalign in eventaligns {
                let data_map = eventalign
                    .signal_iter()
                    .map(|s| (s.pos(), s))
                    .collect::<FnvHashMap<_, _>>();
                for signal in eventalign.signal_iter() {
                    if signal.kmer().starts_with(MOTIF) {
                        let mut kmers = Vec::new();
                        let surround_idx = signal.pos() + MOD_START;
                        let surrounding = surround_idx - (5 + MOD_START)..=surround_idx;
                        for surr in surrounding {
                            if let Some(&s) = data_map.get(&surr) {
                                let kmer = s.kmer();
                                if self.pos_model.gmms().contains_key(kmer)
                                    && self.neg_model.gmms().contains_key(kmer)
                                {
                                    let pos_model = self.pos_model.gmms()[kmer].mixture();
                                    let neg_model = self.neg_model.gmms()[kmer].single();

                                    let pos_sum = s.score_lnsum(&pos_model);
                                    let neg_sum = s.score_lnsum(&neg_model);
                                    kmers.push(SignalScore::new(s, pos_sum, neg_sum));
                                }
                            }
                        }
                    }
                }
            }
            Ok(())
        })?;
        Ok(())
    }
}
