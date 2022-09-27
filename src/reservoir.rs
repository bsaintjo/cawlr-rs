//! Implements reservoir sampling for sampling scores from samples. This allows
//! cawlr train/model-scores to get a fairer representation of values.
//!
//! Partly necessary because arrow2 FileMetadata blocks value is private in
//! version 0.13, which is needed to determine exactly how many chunks are in an
//! Arrow file.
//!
//! This crate aims to implement both the L and R implementations based
//! on <https://en.wikipedia.org/wiki/Reservoir_sampling>.

use std::collections::HashMap;

use anyhow::Result;
use rand::{rngs::SmallRng, Rng, SeedableRng};
use serde::{Deserialize, Serialize};
use typed_sled::Tree;

use crate::{arrow::Signal, Score};

// TOTRY: move count into separate Hashmap on the Reservoir so it
#[derive(Serialize, Deserialize, Default)]
struct ScoreReservoir {
    scores: Vec<f64>,
}

impl ScoreReservoir {
    fn fill_scores(&mut self, scores: &mut [f64], capacity: usize) {
        todo!()
    }

    fn replace(&mut self, score: f64) {
        todo!()
    }
}

struct Reservoir {
    samples: usize,
    rng: SmallRng,
    counts: HashMap<String, usize>,
    tree: Tree<String, ScoreReservoir>,
}

impl Reservoir {
    fn new(samples: usize, tree: Tree<String, ScoreReservoir>) -> Self {
        let rng = SmallRng::seed_from_u64(2456);
        let counts = HashMap::new();
        Reservoir {
            samples,
            rng,
            counts,
            tree,
        }
    }

    fn add_samples_r(&mut self, score: &Signal) -> Result<()> {
        let kmer = score.kmer().to_string();
        log::debug!("Adding samples for kmer {kmer}");
        let mut scores = score.samples().to_owned();
        let kcount = self.counts.entry(kmer.clone()).or_default();
        if *kcount >= self.samples {
            log::debug!("Kmer full, replacing reservoir");
            let mut acc = Vec::new();
            for x in scores {
                let chance = self.rng.gen_range(0..*kcount);
                log::debug!("count: {}, chance: {chance}", *kcount);
                if chance < self.samples {
                    acc.push((chance, x));
                }
                *kcount += 1;
            }
            if !acc.is_empty() {
                self.tree.fetch_and_update(&kmer, |sr| {
                    let mut sr = sr.unwrap();
                    for &(idx, x) in acc.iter() {
                        sr.scores[idx] = x;
                    }
                    Some(sr)
                })?;
            }
        } else {
            log::debug!("Filling values for kmer");
            *kcount += scores.len();
            log::debug!("Kmer {kmer} Reservoir count: {}", *kcount);
            self.tree.fetch_and_update(&kmer, |sr| {
                let mut sr = sr.unwrap_or_default();
                sr.scores.append(&mut scores);
                Some(sr)
            })?;
            // s.scores.append(&mut scores);
        }
        Ok(())
    }

    // fn add_samples_l(&mut self, signal: &Signal) -> Result<()> {
    //     let kmer = signal.kmer().to_string();
    //     let mut w: f64 = (self.rng.gen::<f64>() / self.samples as
    // f64).ln().exp();     if let Some(mut s) = self.tree.get(&kmer)? {
    //         let mut scores = signal.samples().to_owned();
    //         if s.count >= self.samples {
    //             for x in scores {
    //                 s.count += 1; // Essentially an index of the number of times
    // seen                 let chance =
    //                     s.count + ((self.rng.gen::<f64>().ln()) / (1. -
    // w).ln()).floor() as usize;                 if chance <= self.samples {
    //                     let rand_idx = self.rng.gen_range(0..self.samples);
    //                     s.scores[rand_idx] = x;
    //                     w *= (self.rng.gen::<f64>() / self.samples as
    // f64).ln().exp();                 }
    //             }
    //         } else {
    //             s.count += scores.len();
    //             s.scores.append(&mut scores);
    //         }
    //         self.tree.insert(&kmer, &s)?;
    //     }
    //     Ok(())
    // }
}

#[cfg(test)]
mod test {
    use assert_fs::{prelude::PathChild, TempDir};
    use sled::Config;

    use super::*;

    #[test_log::test]
    fn test_reservoir_r() {
        let tmp_dir = TempDir::new().unwrap();
        let db_path = tmp_dir.child("db");
        let db = Config::new().path(db_path).temporary(true).open().unwrap();
        let tree = typed_sled::Tree::open(&db, "id");
        let mut reservoir = Reservoir::new(100, tree);
        for x in 0..1000 {
            let x = x as f64;
            let signal = Signal::new(1, String::from("AAAAAA"), 0.0, 0.0, vec![x]);
            reservoir.add_samples_r(&signal).unwrap();
        }

        pretty_assertions::assert_eq!(reservoir.counts["AAAAAA"], 1000)
    }
}
