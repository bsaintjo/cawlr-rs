use std::{
    collections::{HashMap, HashSet},
    fs::File,
};

use anyhow::Result;
use rand::{SeedableRng, prelude::SmallRng};
use rv::{
    prelude::{Gaussian, Mixture},
    traits::Rv,
};
use serde_pickle::from_reader;

use crate::train::ModelDB;

// TODO switch to parquet
pub(crate) fn load_models(input: &str) -> Result<ModelDB> {
    let file = File::open(input)?;
    let model_db = from_reader(file, Default::default())?;
    Ok(model_db)
}

// TODO switch to parquet
pub(crate) fn save_kmer_ranks(output: &str, kmer_ranks: HashMap<String, f64>) -> Result<()> {
    let mut file = File::create(output)?;
    serde_pickle::to_writer(&mut file, &kmer_ranks, Default::default())?;
    Ok(())
}


pub struct RankOptions {
    rng: SmallRng,
    n_samples: usize
}

impl RankOptions {
    pub(crate) fn new(seed: u64, n_samples: usize) -> Self {
        let rng = SmallRng::seed_from_u64(seed);
        RankOptions {
            rng,
            n_samples
        }
    }

    // Approximate the Kulback-Leibler Divergence for the two GMMs as mentioned in
    // J. R. Hershey and P. A. Olsen, "Approximating the Kullback Leibler Divergence
    // Between Gaussian Mixture Models," 2007 IEEE International Conference on
    // Acoustics, Speech and Signal Processing - ICASSP '07, 2007, pp.
    // IV-317-IV-320, doi: 10.1109/ICASSP.2007.366913.
    fn kl_approx(&mut self, pos_ctrl: &Mixture<Gaussian>, neg_ctrl: &Mixture<Gaussian>) -> f64 {
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

    pub(crate) fn rank(&mut self, pos_ctrl: ModelDB, neg_ctrl: ModelDB) -> HashMap<String, f64> {
        let mut kmer_ranks = HashMap::new();
        let pos_ctrl_kmers = pos_ctrl.keys().collect::<HashSet<&String>>();
        let neg_ctrl_kmers = neg_ctrl.keys().collect::<HashSet<&String>>();
        let kmers = pos_ctrl_kmers.intersection(&neg_ctrl_kmers);
        for &kmer in kmers {
            let pos_ctrl_model = &pos_ctrl[kmer];
            let neg_ctrl_model = &neg_ctrl[kmer];
            let kl = self.kl_approx(pos_ctrl_model, neg_ctrl_model);
            kmer_ranks.insert(kmer.clone(), kl);
        }
        kmer_ranks
    }
}