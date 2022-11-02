use std::io::{Read, Seek};

use eyre::Result;
use linfa::{
    traits::{Fit, Transformer},
    DatasetBase, ParamGuard,
};
use linfa_clustering::{Dbscan, GaussianMixtureModel};
use ndarray::Array;
use rv::prelude::{Gaussian, Mixture};

use crate::{arrow_utils::load_read_arrow, train::Model, Eventalign};

pub struct TrainOptions {
    n_samples: usize,
    single: bool,
    dbscan: bool,
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

impl TrainOptions {
    pub fn run<R>(mut self, input: R) -> Result<Model>
    where
        R: Read + Seek,
    {
        let db = sled::Config::new().temporary(true).open()?;
        let tree = typed_sled::Tree::<String, Vec<f64>>::open(&db, "npsmlr_train");
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
        let data = DatasetBase::from(means);
        let min_points = 3;
        let data = Dbscan::params(min_points)
            .tolerance(1e-3)
            .check()
            .unwrap()
            .transform(data);
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
            return None;
        }

        let len = obs.len();
        let shape = (len, 1);
        let means = Array::from_shape_vec(shape, obs).unwrap();
        let data = DatasetBase::from(means);

        let n_clusters = 2;
        let n_runs = 10;
        let tolerance = 1e-4f64;
        let gmm = GaussianMixtureModel::params(n_clusters)
            .n_runs(n_runs)
            .tolerance(tolerance)
            .check()
            .unwrap()
            .fit(&data)
            .unwrap();
        todo!()
    }
}
