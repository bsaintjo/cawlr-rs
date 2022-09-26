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

#[derive(Serialize, Deserialize)]
struct ScoreReservoir {
    count: usize,
    len: usize,
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
                    if chance <= s.len {
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
        todo!()
    }
}
