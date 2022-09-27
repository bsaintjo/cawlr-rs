//! Implements reservoir sampling for sampling scores from samples. This allows
//! cawlr train/model-scores to get a fairer representation of values.
//!
//! Partly necessary because arrow2 FileMetadata blocks value is private in
//! version 0.13, which is needed to determine exactly how many chunks are in an
//! Arrow file.
//!
//! This crate aims to implement both the L and R implementations based
//! on <https://en.wikipedia.org/wiki/Reservoir_sampling>.

use anyhow::Result;
use rand::{rngs::SmallRng, Rng};
use serde::{Deserialize, Serialize};
use typed_sled::Tree;

use crate::arrow::Signal;

// TOTRY: move count into separate Hashmap on the Reservoir so it 
#[derive(Serialize, Deserialize)]
struct ScoreReservoir {
    count: usize,
    scores: Vec<f64>,
}

struct Reservoir {
    samples: usize,
    rng: SmallRng,
    tree: Tree<String, ScoreReservoir>,
}

impl Reservoir {
    fn add_samples_r(&mut self, score: &Signal) -> Result<()> {
        let kmer = score.kmer().to_string();
        if let Some(mut s) = self.tree.get(&kmer)? {
            let mut scores = score.samples().to_owned();
            if s.count >= self.samples {
                for x in scores {
                    let chance = self.rng.gen_range(0..s.count);
                    if chance <= self.samples {
                        s.scores[chance] = x;
                    }
                    s.count += 1;
                }
            } else {
                s.count += scores.len();
                s.scores.append(&mut scores);
            }
            self.tree.insert(&kmer, &s)?;
        }
        Ok(())
    }

    fn add_samples_l(&mut self, signal: &Signal) -> Result<()> {
        let kmer = signal.kmer().to_string();
        let mut w: f64 = (self.rng.gen::<f64>() / self.samples as f64).ln().exp();
        if let Some(mut s) = self.tree.get(&kmer)? {
            let mut scores = signal.samples().to_owned();
            if s.count >= self.samples {
                for x in scores {
                    s.count += 1; // Essentially an index of the number of times seen
                    let chance =
                        s.count + ((self.rng.gen::<f64>().ln()) / (1. - w).ln()).floor() as usize;
                    if chance <= self.samples {
                        let rand_idx = self.rng.gen_range(0..self.samples);
                        s.scores[rand_idx] = x;
                        w *= (self.rng.gen::<f64>() / self.samples as f64).ln().exp();
                    }
                }
            } else {
                s.count += scores.len();
                s.scores.append(&mut scores);
            }
            self.tree.insert(&kmer, &s)?;
        }
        Ok(())
    }
}
