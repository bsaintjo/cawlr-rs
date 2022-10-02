use std::{
    borrow::Borrow,
    collections::HashMap,
    fmt::Debug,
    fs::File,
    path::{Path, PathBuf},
};

use bio::io::fasta::IndexedReader;
use eyre::Result;
use fnv::{FnvHashMap, FnvHashSet};
use linfa::{
    traits::{Fit, Transformer},
    DatasetBase, ParamGuard,
};
use linfa_clustering::{Dbscan, GaussianMixtureModel};
use ndarray::Array;
use rayon::prelude::{IntoParallelIterator, ParallelIterator};
use rv::prelude::{Gaussian, Mixture};
use serde::{Deserialize, Serialize};

use crate::arrow::{load_apply, Eventalign, MetadataExt};

pub(crate) type ModelDB = FnvHashMap<String, ModelParams>;
type KmerMeans = FnvHashMap<String, Vec<f64>>;

#[derive(Debug, Serialize, Deserialize, PartialEq)]
pub struct ModelParams {
    // weight for a
    weight: f64,
    mu_a: f64,
    sigma_a: f64,

    // weight is 1 - weight
    mu_b: f64,
    sigma_b: f64,
}

impl ModelParams {
    pub fn new(weight: f64, mu_a: f64, sigma_a: f64, mu_b: f64, sigma_b: f64) -> Self {
        Self {
            weight,
            mu_a,
            sigma_a,
            mu_b,
            sigma_b,
        }
    }

    fn weight_a(&self) -> f64 {
        self.weight
    }

    fn weight_b(&self) -> f64 {
        1. - self.weight
    }

    pub fn mixture(&self) -> Mixture<Gaussian> {
        let g1 = Gaussian::new_unchecked(self.mu_a, self.sigma_a);
        let g2 = Gaussian::new_unchecked(self.mu_b, self.sigma_b);
        let components = vec![g1, g2];
        let weights = vec![self.weight_a(), self.weight_b()];
        Mixture::new_unchecked(weights, components)
    }
}

impl<T: Borrow<Mixture<Gaussian>>> From<T> for ModelParams {
    fn from(mix: T) -> Self {
        let mix: &Mixture<Gaussian> = mix.borrow();
        let weight = mix.weights()[0];
        let components = mix.components();
        let mu_a = components[0].mu();
        let sigma_a = components[0].sigma();

        let mu_b = components[1].mu();
        let sigma_b = components[1].sigma();
        ModelParams::new(weight, mu_a, sigma_a, mu_b, sigma_b)
    }
}

#[derive(Default, Serialize, Deserialize)]
pub struct Model {
    gmms: ModelDB,
    skips: FnvHashMap<String, f64>,
}

impl Model {
    /// Get a reference to the model's gmms.
    pub(crate) fn gmms(&self) -> &ModelDB {
        &self.gmms
    }

    /// Get a reference to the model's skips.
    pub(crate) fn skips(&self) -> &FnvHashMap<String, f64> {
        &self.skips
    }
}

impl Model {
    pub(crate) fn new(gmms: ModelDB, skips: FnvHashMap<String, f64>) -> Self {
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

    fn plus_total(&mut self) {
        self.total += 1;
    }

    fn had_score(&mut self, is_score: bool) {
        if is_score {
            self.plus_both();
        } else {
            self.plus_total();
        }
    }
}

impl Default for Skips {
    fn default() -> Self {
        Skips::new(0, 0)
    }
}

struct KmerSkips(FnvHashMap<Vec<u8>, Skips>);

impl KmerSkips {
    fn new() -> Self {
        Self(FnvHashMap::default())
    }
}

pub struct Train {
    acc: KmerMeans,
    skips: KmerSkips,
    genome: IndexedReader<File>,
    feather: PathBuf,
    samples: usize,
}

impl Train {
    pub fn try_new<P, Q>(filename: P, genome: Q, samples: usize) -> Result<Self, eyre::Error>
    where
        P: AsRef<Path>,
        Q: AsRef<Path> + Debug,
    {
        let genome =
            IndexedReader::from_file(&genome).map_err(|_| eyre::eyre!("Failed to read genome."))?;
        let feather = filename.as_ref().to_owned();
        Ok(Self {
            acc: FnvHashMap::default(),
            skips: KmerSkips::new(),
            genome,
            feather,
            samples,
        })
    }

    fn kmer_means_insufficient(&self) -> bool {
        self.acc.is_empty() || insufficient(&self.acc, self.samples)
    }

    fn kmer_skips_insufficient(&self) -> bool {
        self.skips.0.is_empty() || self.skips.0.values().any(|x| x.total < self.samples)
    }

    pub fn run(mut self) -> Result<Model> {
        let file = File::open(&self.feather)?;
        load_apply(file, |eventaligns| {
            for eventalign in eventaligns.into_iter() {
                if self.kmer_means_insufficient() || self.kmer_skips_insufficient() {
                    self.read_to_kmer_means(&eventalign);
                    self.read_to_skip_counts(&eventalign)?;
                }
            }
            Ok(())
        })?;

        // let mut gmms = self.acc;
        let gmms = self
            .acc
            .into_par_iter()
            .filter_map(|item| {
                if let Ok(Some(gmm)) = train_gmm(item.1) {
                    Some((item.0, ModelParams::from(gmm)))
                } else {
                    None
                }
            })
            .collect();

        // for (kmer, kmer_mean) in x {
        //     if kmer_mean.len() > 1 {
        //         let gmm = train_gmm(kmer_mean);
        //         if let Ok(Some(gmm)) = gmm {
        //             gmms.insert(kmer, gmm);
        //         } else {
        //             log::warn!("Failed to train for kmer {kmer}.");
        //         }
        //     }
        // }

        let mut ratios = FnvHashMap::default();
        for (kmer, skips) in self.skips.0.into_iter() {
            let kmer = String::from_utf8(kmer)?;
            let ratio = (skips.count as f64) / (skips.total as f64);
            ratios.insert(kmer, ratio);
        }

        let model = Model::new(gmms, ratios);

        Ok(model)
    }

    fn read_to_kmer_means(&mut self, read: &Eventalign) {
        for signal in read.signal_iter() {
            let kmer = signal.kmer().to_owned();
            let entry = self.acc.entry(kmer).or_default();
            if entry.len() > self.samples {
                continue;
            }
            let samples = signal.samples();
            entry.extend_from_slice(samples);
        }
    }

    fn read_to_skip_counts(&mut self, read: &Eventalign) -> Result<()> {
        let mut pos_scores = FnvHashSet::default();
        for signal in read.signal_iter() {
            pos_scores.insert(signal.pos());
        }
        let read_seq = self.get_read_seq(read)?;
        for (kmer, pos) in read_seq.windows(6).zip(read.start_0b()..) {
            let has_score = pos_scores.contains(&pos);
            let kskip = self.skips.0.entry(kmer.to_owned()).or_default();
            kskip.had_score(has_score);
        }
        Ok(())
    }

    /// Get a mutable reference to the train's genome.
    pub(crate) fn genome_mut(&mut self) -> &mut IndexedReader<File> {
        &mut self.genome
    }

    // TODO: Use Context instead
    // Genome fasta reader method makes clippy think its wrong but it is actually
    // correct.
    #[allow(clippy::read_zero_byte_vec)]
    fn get_read_seq(&mut self, read: &Eventalign) -> Result<Vec<u8>> {
        let strand = read.strand();
        let chrom = read.chrom();
        let start = read.start_0b();
        self.genome_mut()
            .fetch(chrom, start, read.seq_stop_1b_excl())?;
        let mut seq = Vec::new();
        self.genome_mut().read(&mut seq)?;
        let seq = if strand == crate::arrow::Strand::plus() {
            seq
        } else {
            bio::alphabets::dna::revcomp(seq)
        };
        Ok(seq)
    }
}

fn train_gmm(means: Vec<f64>) -> Result<Option<Mixture<Gaussian>>> {
    let len = means.len();
    let shape = (len, 1);
    let means = Array::from_shape_vec(shape, means)?;
    let data = DatasetBase::from(means);
    let min_points = 3;
    let data = Dbscan::params(min_points)
        .tolerance(1e-3)
        .check()?
        .transform(data);

    // Filter the noise from the Dbscan output
    // TODO: Should be a faster way to do this but for now convert to vec and back
    // before training

    let recs = data
        .records()
        .as_slice()
        .expect("Getting records failed after DBSCAN");
    let targets = data
        .targets()
        .as_slice()
        .expect("Getting targets failed after DBSCAN");

    let obs: Vec<f64> = recs
        .iter()
        .zip(targets.iter())
        .filter_map(
            |(&x, cluster)| {
                if cluster.is_some() {
                    Some(x)
                } else {
                    None
                }
            },
        )
        .collect();
    if obs.len() < 2 {
        log::warn!("Not enough values left in observations");
        return Ok(None);
    }

    let len = obs.len();
    let shape = (len, 1);
    let means = Array::from_shape_vec(shape, obs)?;
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
        .map(|(&mean, &sigma)| Gaussian::new_unchecked(mean, sigma.sqrt()))
        .collect::<Vec<Gaussian>>();
    let mm = Mixture::new_unchecked(weights, gausses);

    Ok(Some(mm))
}

fn insufficient<K, V, S>(dict: &HashMap<K, Vec<V>, S>, n: usize) -> bool {
    dict.values().any(|f| f.len() < n)
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_insufficient() {
        let n = 5;
        let x = vec![0; 10];
        let y = vec![1; 20];
        let mut dict = HashMap::new();
        dict.insert('x', x);
        dict.insert('y', y);
        assert!(!insufficient(&dict, n));

        let n = 40;
        assert!(insufficient(&dict, n))
    }

    #[test]
    fn test_model_params() {
        let g1 = Gaussian::new(1., 2.).unwrap();
        let g2 = Gaussian::new(3., 4.).unwrap();
        let components = vec![g1, g2];
        let weights = vec![0.7, 0.3];

        let mix = Mixture::new_unchecked(weights, components);
        let params = ModelParams::from(&mix);

        let answer = ModelParams::new(0.7, 1., 2., 3., 4.);

        pretty_assertions::assert_eq!(params, answer);
    }
}
