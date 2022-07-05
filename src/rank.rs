use fnv::{FnvHashMap, FnvHashSet};
use rand::{prelude::SmallRng, SeedableRng};
use rv::{
    prelude::{Gaussian, Mixture},
    traits::Rv,
};
use serde::{Deserialize, Serialize};

use crate::{
    score::{choose_model, choose_pos_model},
    train::Model,
};

pub type Ranks = FnvHashMap<String, f64>;

pub struct RankOptions {
    rng: SmallRng,
    n_samples: usize,
}

impl Default for RankOptions {
    fn default() -> Self {
        let rng = SmallRng::seed_from_u64(2456);
        RankOptions {
            rng,
            n_samples: 10_000,
        }
    }
}

#[derive(Serialize, Deserialize)]
pub(crate) struct Rankings {
    ranks: FnvHashMap<String, f64>,
    hi_rank_pos_idx: FnvHashMap<String, usize>,
    p_values: FnvHashMap<String, f64>,
}

impl Rankings {
    fn new(
        ranks: FnvHashMap<String, f64>,
        hi_rank_pos_idx: FnvHashMap<String, usize>,
        p_values: FnvHashMap<String, f64>,
    ) -> Self {
        Self {
            ranks,
            hi_rank_pos_idx,
            p_values,
        }
    }

    pub(crate) fn ranks(&self) -> &FnvHashMap<String, f64> {
        &self.ranks
    }

    pub(crate) fn hi_rank_pos_idx(&self) -> &FnvHashMap<String, usize> {
        &self.hi_rank_pos_idx
    }

    pub(crate) fn p_values(&self) -> &FnvHashMap<String, f64> {
        &self.p_values
    }
}

impl RankOptions {
    pub(crate) fn new(seed: u64, n_samples: usize) -> Self {
        let rng = SmallRng::seed_from_u64(seed);
        RankOptions { rng, n_samples }
    }

    // Approximate the Kulback-Leibler Divergence for the two GMMs as mentioned in
    // J. R. Hershey and P. A. Olsen, "Approximating the Kullback Leibler Divergence
    // Between Gaussian Mixture Models," 2007 IEEE International Conference on
    // Acoustics, Speech and Signal Processing - ICASSP '07, 2007, pp.
    // IV-317-IV-320, doi: 10.1109/ICASSP.2007.366913.
    //
    // TODO: Check if some normalization is required for this
    fn kl_approx(&mut self, pos_ctrl: &Mixture<Gaussian>, neg_ctrl: &Mixture<Gaussian>) -> f64 {
        let neg_ctrl = choose_model(neg_ctrl);
        let pos_ctrl = choose_pos_model(neg_ctrl, pos_ctrl);
        let samples: Vec<f32> = pos_ctrl.sample(self.n_samples, &mut self.rng);
        let total: f64 = samples
            .into_iter()
            .map(|sample| {
                let p = pos_ctrl.ln_f(&sample);
                let n = neg_ctrl.ln_f(&sample);
                p - n
            })
            .sum();
        total / self.count()
    }

    fn count(&self) -> f64 {
        self.n_samples as f64
    }

    pub fn rank(&mut self, pos_ctrl: &Model, neg_ctrl: &Model) -> Ranks {
        let mut kmer_ranks = FnvHashMap::default();
        let pos_ctrl_kmers = pos_ctrl.gmms().keys().collect::<FnvHashSet<&String>>();
        let neg_ctrl_kmers = neg_ctrl.gmms().keys().collect::<FnvHashSet<&String>>();
        let kmers = pos_ctrl_kmers.intersection(&neg_ctrl_kmers);
        for &kmer in kmers {
            let pos_ctrl_model = &pos_ctrl.gmms()[kmer];
            let neg_ctrl_model = &neg_ctrl.gmms()[kmer];
            let kl = self.kl_approx(pos_ctrl_model, neg_ctrl_model);
            kmer_ranks.insert(kmer.clone(), kl);
        }
        kmer_ranks
    }
}
