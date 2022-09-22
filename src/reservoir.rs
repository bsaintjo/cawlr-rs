//! Implements reservoir sampling for sampling scores from samples. This allows
//! cawlr train/model-scores to get a fairer representation of values.
//!
//! Partly necessary because arrow2 FileMetadata blocks value is private in
//! version 0.13, which is needed to determine exactly how many chunks are in an
//! Arrow file.
//!
//! This crate aims to implement both the L and R implementations based
//! on https://en.wikipedia.org/wiki/Reservoir_sampling.

use std::collections::HashMap;

use rand::{rngs::SmallRng, Rng};

use crate::arrow::Score;

struct ScoreReservoir {
    count: usize,
    scores: Vec<Score>,
}

struct Reservoir {
    samples: usize,
    rng: SmallRng,
    kmers: HashMap<String, ScoreReservoir>,
}

impl Reservoir {
    fn is_full(&self, kmer: &str) -> bool {
        self.kmers[kmer].count >= self.samples
    }

    fn add_sample(&mut self, score: &Score) {
        let kmer = score.kmer();
        if !self.is_full(kmer) {
            let sr = self.kmers.get_mut(kmer).unwrap();
            sr.scores.push(score.clone());
            sr.count += 1;
        } else {
            let count = self.kmers[kmer].count;
            let j = self.rng.gen_range(1..count);
            if j < self.samples {
                self.kmers.get_mut(kmer).unwrap().scores[j] = score.clone();
            }
        }
    }
}
