use std::io::{Read, Seek, Write};

use eyre::Result;
use linfa::{
    traits::{Fit, Transformer},
    DatasetBase, ParamGuard,
};
use linfa_clustering::{Dbscan, GaussianMixtureModel};
use ndarray::Array;
use rv::prelude::{Gaussian, Mixture};
use typed_sled::Tree;

use crate::{
    arrow_utils::load_read_arrow,
    motif::{all_bases, Motif},
    train::{mix_to_mix, Model},
    utils::CawlrIO,
    Eventalign,
};

pub struct TrainOptions {
    n_samples: usize,
    single: bool,
    dbscan: bool,
    motifs: Vec<Motif>,
}

fn extend_merge(
    _key: String,
    old_value: Option<Vec<f64>>,
    mut new_value: Vec<f64>,
) -> Option<Vec<f64>> {
    if let Some(mut ov) = old_value {
        ov.append(&mut new_value);
        Some(ov)
    } else {
        Some(new_value)
    }
}

impl Default for TrainOptions {
    fn default() -> Self {
        TrainOptions {
            n_samples: 50000,
            single: false,
            dbscan: false,
            motifs: all_bases(),
        }
    }
}

impl TrainOptions {
    pub fn n_samples(mut self, n_samples: usize) -> Self {
        self.n_samples = n_samples;
        self
    }

    pub fn single(mut self, single: bool) -> Self {
        self.single = single;
        self
    }

    pub fn dbscan(mut self, dbscan: bool) -> Self {
        self.dbscan = dbscan;
        self
    }

    pub fn motifs(mut self, motifs: Vec<Motif>) -> Self {
        self.motifs = motifs;
        self
    }

    pub fn run<R, W>(self, input: R, mut writer: W) -> Result<()>
    where
        R: Read + Seek,
        W: Write,
    {
        let model = self.run_model(input)?;
        model.save(&mut writer)?;
        Ok(())
    }

    pub fn run_model<R>(self, input: R) -> Result<Model>
    where
        R: Read + Seek,
    {
        let db = sled::Config::new().temporary(true).open()?;
        let tree = Tree::open(&db, "npsmlr_train");
        tree.set_merge_operator(extend_merge);
        load_read_arrow(input, |eventaligns: Vec<Eventalign>| {
            for eventalign in eventaligns.into_iter() {
                for signal in eventalign.signal_iter() {
                    let kmer = signal.kmer();
                    let samples = signal.samples();
                    tree.merge(&kmer.to_string(), &samples.to_vec())?;
                }
            }
            Ok(())
        })?;

        self.train_gmms(tree)
    }

    fn train_gmms(&self, tree: Tree<String, Vec<f64>>) -> Result<Model> {
        let mut model = Model::default();
        for item in tree.iter() {
            let (kmer, samples) = item?;
            if let Some(gmm) = self.train_gmm(samples) {
                model.insert_gmm(kmer, gmm);
            }
        }
        Ok(model)
    }

    fn train_gmm(&self, samples: Vec<f64>) -> Option<Mixture<Gaussian>> {
        let len = samples.len();
        let shape = (len, 1);
        let means = Array::from_shape_vec(shape, samples).unwrap();
        let mut data = DatasetBase::from(means);
        if self.dbscan {
            let min_points = 3;
            let dataset = Dbscan::params(min_points)
                .tolerance(1e-3)
                .check()
                .unwrap()
                .transform(data);
            let recs = dataset
                .records()
                .as_slice()
                .expect("Getting records failed after DBSCAN");
            let targets = dataset
                .targets()
                .as_slice()
                .expect("Getting targets failed after DBSCAN");

            let filtered: Vec<f64> = recs
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
            if filtered.len() < 2 {
                log::warn!("Not enough values left in observations");
                return None;
            }

            let len = filtered.len();
            let shape = (len, 1);
            let filtered_results = Array::from_shape_vec(shape, filtered).unwrap();
            data = DatasetBase::from(filtered_results);
        }

        let n_clusters = if self.single { 1 } else { 2 };
        let n_runs = 10;
        let tolerance = 1e-4f64;
        let gmm = GaussianMixtureModel::params(n_clusters)
            .n_runs(n_runs)
            .tolerance(tolerance)
            .check()
            .unwrap()
            .fit(&data)
            .unwrap();
        let mm = mix_to_mix(&gmm);
        Some(mm)
    }
}
