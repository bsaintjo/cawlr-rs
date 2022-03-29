use std::collections::{HashMap, HashSet};

use anyhow::Result;
use linfa::{traits::Fit, DatasetBase, ParamGuard};
use linfa_clustering::GaussianMixtureModel;
use ndarray::Array;
use rv::prelude::{Gaussian, Mixture};
use serde::{Deserialize, Serialize};

use crate::reads::PreprocessRead;

pub(crate) type ModelDB = HashMap<String, Mixture<Gaussian>>;
type KmerMeans = HashMap<String, Vec<f64>>;

#[derive(Serialize, Deserialize)]
pub(crate) struct Model {
    gmms: ModelDB,
    skips: HashMap<String, f64>,
}

impl Model {
    /// Get a reference to the model's gmms.
    pub(crate) fn gmms(&self) -> &ModelDB {
        &self.gmms
    }

    /// Get a reference to the model's skips.
    pub(crate) fn skips(&self) -> &HashMap<String, f64> {
        &self.skips
    }
}

impl Model {
    pub(crate) fn new(gmms: ModelDB, skips: HashMap<String, f64>) -> Self {
        Self { gmms, skips }
    }
}

struct Skips {
    count: usize,
    total: usize,
}

impl Skips {
    fn new(count: usize, total: usize) -> Self {
        Self { count, total }
    }

    fn plus_both(&mut self) {
        self.total += 1;
        self.count += 1;
    }

    fn plus_count(&mut self) {
        self.count += 1;
    }

    fn had_score(&mut self, is_score: bool) {
        if is_score {
            self.plus_both();
        } else {
            self.plus_count();
        }
    }
}

impl Default for Skips {
    fn default() -> Self {
        Skips::new(0, 0)
    }
}

struct KmerSkips(HashMap<Vec<u8>, Skips>);

impl KmerSkips {
    fn new() -> Self {
        Self(HashMap::new())
    }
}

pub(crate) struct Train {
    acc: KmerMeans,
    skips: KmerSkips,
}

impl Train {
    pub(crate) fn new() -> Self {
        Self {
            acc: HashMap::new(),
            skips: KmerSkips::new(),
        }
    }

    pub(crate) fn run(mut self, reads: Vec<PreprocessRead>) -> Result<Model> {
        for read in reads.into_iter() {
            self.read_to_kmer_means(&read);
            self.read_to_skip_counts(&read);
        }

        let mut gmms = HashMap::new();
        for (kmer, kmer_mean) in self.acc.into_iter() {
            if kmer_mean.len() > 1 {
                let gmm = train_gmm(kmer_mean)?;
                gmms.insert(kmer, gmm);
            }
        }

        let mut ratios = HashMap::new();
        for (kmer, skips) in self.skips.0.into_iter() {
            let kmer = String::from_utf8(kmer)?;
            let ratio = (skips.count as f64) / (skips.total as f64);
            ratios.insert(kmer, ratio);
        }

        let model = Model::new(gmms, ratios);

        Ok(model)
    }

    fn read_to_kmer_means(&mut self, read: &PreprocessRead) {
        for ld in read.iter() {
            let mean = ld.mean();
            let kmer = ld.kmer().to_owned();
            self.acc.entry(kmer).or_default().push(mean);
        }
    }

    fn read_to_skip_counts(&mut self, read: &PreprocessRead) {
        let mut pos_scores = HashSet::new();
        for ld in read.iter() {
            pos_scores.insert(ld.pos() as usize);
        }
        for (kmer, pos) in read.seq().windows(6).zip(read.start()..) {
            let has_score = pos_scores.contains(&pos);
            let kskip = self.skips.0.entry(kmer.to_owned()).or_default();
            kskip.had_score(has_score);
        }
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
