use std::collections::HashMap;

use anyhow::Result;
use linfa::{traits::Fit, DatasetBase, ParamGuard};
use linfa_clustering::GaussianMixtureModel;
use ndarray::Array;
use rv::prelude::{Gaussian, Mixture};

use crate::reads::PreprocessRead;

pub(crate) type ModelDB = HashMap<String, Mixture<Gaussian>>;
type KmerMeans = HashMap<String, Vec<f64>>;

pub(crate) struct Train {
    acc: KmerMeans,
}

impl Train {
    pub(crate) fn new() -> Self {
        Self {
            acc: HashMap::new(),
        }
    }

    pub(crate) fn run(mut self, reads: Vec<PreprocessRead>) -> Result<ModelDB> {
        let mut gmms = HashMap::new();
        self.reads_to_kmer_means(reads);
        for (kmer, kmer_mean) in self.acc.into_iter() {
            let gmm = train_gmm(kmer_mean)?;
            gmms.insert(kmer, gmm);
        }

        Ok(gmms)
    }

    fn read_to_kmer_means(&mut self, read: PreprocessRead) {
        for ld in read.into_iter() {
            let mean = ld.mean();
            let kmer = ld.into_kmer();
            self.acc.entry(kmer).or_insert(Vec::new()).push(mean);
        }
    }

    fn reads_to_kmer_means(&mut self, reads: Vec<PreprocessRead>) {
        reads.into_iter().for_each(|lr| self.read_to_kmer_means(lr))
    }
}

fn train_gmm(means: Vec<f64>) -> Result<Mixture<Gaussian>> {
    let len = means.len();
    let shape = (len, 1);
    let means = Array::from_shape_vec(shape, means)?;
    let data = DatasetBase::from(means);

    let n_clusters = 2;
    let n_runs = 10;
    let tolerance = 1e-4f64;
    let gmm = GaussianMixtureModel::params(n_clusters)
        .n_runs(n_runs)
        .tolerance(tolerance)
        .check()?
        .fit(&data)?;
    let weights = gmm.weights().iter().cloned().collect::<Vec<f64>>();
    let means = gmm.means().iter();
    let covs = gmm.covariances().iter();
    let gausses = means
        .zip(covs)
        .map(|(&mean, &sigma)| Gaussian::new(mean, sigma).unwrap())
        .collect::<Vec<Gaussian>>();
    let mm = Mixture::new_unchecked(weights, gausses);

    Ok(mm)
}
