use std::{
    collections::{HashMap, HashSet},
    io::{Read, Seek, Write},
};

use eyre::Result;
// use itertools::Itertools;
use linfa::{
    traits::{Fit, Transformer},
    DatasetBase, ParamGuard,
};
use linfa_clustering::{Dbscan, GaussianMixtureModel};
use ndarray::Array;
use rusqlite::Connection;
use rv::prelude::{Gaussian, Mixture};

use crate::{
    arrow_utils::load_read_arrow,
    motif::{all_bases, Motif},
    train::{mix_to_mix, Model},
    utils::CawlrIO,
    Eventalign,
};

// fn all_kmers() -> Vec<String> {
//     ["A", "C", "G", "T"]
//         .iter()
//         .permutations(6)
//         .map(|xs| xs.iter().join(""))
//         .collect()
// }

pub struct TrainOptions {
    n_samples: usize,
    single: bool,
    dbscan: bool,
    motifs: Vec<Motif>,
}

// fn extend_merge(
//     _key: String,
//     old_value: Option<Vec<f64>>,
//     mut new_value: Vec<f64>,
// ) -> Option<Vec<f64>> {
//     if let Some(mut ov) = old_value {
//         ov.append(&mut new_value);
//         Some(ov)
//     } else {
//         Some(new_value)
//     }
// }

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
        let db_path = std::env::temp_dir().join("npsmlr.db");
        let db = Db::open(db_path.as_os_str().to_str().unwrap())?;
        let mut tree = HashMap::new();
        let mut finished: HashSet<String> = HashSet::new();
        // let db = sled::Config::new().path(tmp_dir).temporary(true).open()?;
        // let tree = Tree::open(&db, "npsmlr_train");
        // tree.set_merge_operator(extend_merge);
        load_read_arrow(input, |eventaligns: Vec<Eventalign>| {
            for eventalign in eventaligns.into_iter() {
                if finished.len() > 4096 {
                    break;
                }

                for signal in eventalign.signal_iter() {
                    let kmer = signal.kmer();
                    if finished.contains(kmer) {
                        continue;
                    }
                    let mut samples = signal.samples().to_vec();
                    let kmer_samples: &mut Vec<f64> = tree.entry(kmer.to_string()).or_default();
                    kmer_samples.append(&mut samples);
                    if kmer_samples.len() > self.n_samples {
                        finished.insert(kmer.to_string());
                    }
                }
            }
            Ok(())
        })?;

        self.train_gmms(tree)
    }

    fn train_gmms(&self, tree: HashMap<String, Vec<f64>>) -> Result<Model> {
        let mut model = Model::default();
        for (kmer, samples) in tree.into_iter() {
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

struct Db(Connection);

impl Db {
    fn open(path: &str) -> eyre::Result<Self> {
        Ok(Db(Connection::open(path)?))
    }
    fn init(&self) -> eyre::Result<()> {
        self.0.execute(
            "CREATE TABLE data (
                id      INTEGER PRIMARY KEY,
                kmer    TEXT NOT NULL,
                sample  REAL NOT NULL
            );",
            (),
        )?;
        Ok(())
    }

    fn add_read(&mut self, eventalign: Eventalign) -> eyre::Result<()> {
        let tx = self.0.transaction()?;
        for signal in eventalign.signal_iter() {
            let kmer = signal.kmer();
            for sample in signal.samples() {
                tx.execute(
                    "INSERT INTO data (kmer, sample) VALUES (?1, ?2)",
                    (kmer, sample),
                )?;
            }
        }

        tx.commit()?;
        Ok(())
    }
}

// #[cfg(test)]
// mod test {
//     use super::*;

//     #[test]
//     fn test_all_kmers() {
//         let kmers = all_kmers();
//         assert_eq!(4096, kmers.len());
//     }
// }
