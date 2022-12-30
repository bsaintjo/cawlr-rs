use std::{
    io::{Read, Seek, Write},
    path::Path,
};

use eyre::Result;
use itertools::Itertools;
use linfa::{
    traits::{Fit, Transformer},
    DatasetBase, ParamGuard,
};
use linfa_clustering::{Dbscan, GaussianMixtureModel, GmmError, GmmInitMethod};
use ndarray::Array;
use rusqlite::{named_params, Connection};
use rv::prelude::{Gaussian, Mixture};

use crate::{
    arrow_utils::load_read_arrow_measured,
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

fn all_kmers() -> Vec<String> {
    let mut kmers: Vec<String> = vec![String::new()];
    let bases = ["A", "C", "G", "T"];
    for _ in 0..6 {
        let mut acc = Vec::new();
        for base in bases {
            for s in kmers.iter() {
                let mut xs = s.clone();
                xs.push_str(base);
                acc.push(xs);
            }
        }
        kmers = acc;
    }
    kmers
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
        let mut db = Db::open(db_path)?;
        load_read_arrow_measured(input, |eventaligns: Vec<Eventalign>| {
            db.add_reads(eventaligns)?;
            Ok(())
        })?;

        self.train_gmms(db)
    }

    fn train_gmms(&self, db: Db) -> Result<Model> {
        let mut model = Model::default();
        for kmer in all_kmers() {
            log::info!("Training on kmer {kmer}");
            let samples = db.get_kmer_samples(&kmer, self.n_samples)?;
            log::info!("n samples: {}", samples.len());
            if !samples.is_empty() {
                if let Some(gmm) = self.train_gmm(samples) {
                    log::info!("Training successful!");
                    model.insert_gmm(kmer, gmm);
                }
            }
        }
        Ok(model)
    }

    fn train_gmm(&self, samples: Vec<f64>) -> Option<Mixture<Gaussian>> {
        let samples = samples.into_iter().filter(|x| x.is_finite()).collect_vec();
        if samples.is_empty() {
            return None;
        }
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
            .fit(&data);
        if let Err(GmmError::MinMaxError(_)) = gmm {
            let data = data.records.into_raw_vec();
            log::error!("Failed with MinMaxError, raw data \n{data:?}");
            return None;
        }
        let mm = mix_to_mix(&gmm.unwrap());
        Some(mm)
    }
}

struct Db(Connection);

impl Db {
    fn open<P: AsRef<Path>>(path: P) -> eyre::Result<Self> {
        let path = path.as_ref();
        if path.exists() {
            std::fs::remove_file(path)?;
        }
        let db = Db(Connection::open(path)?);
        db.init()?;
        db.create_idx()?;
        Ok(db)
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

    fn create_idx(&self) -> eyre::Result<()> {
        self.0.execute("CREATE INDEX kmer_idx on data (kmer)", ())?;
        self.0.pragma_update(None, "journal_mode", "WAL")?;
        self.0.pragma_update(None, "synchronous", "NORMAL")?;
        self.0.pragma_update(None, "cache_size", -64000)?;
        Ok(())
    }

    fn add_reads(&mut self, es: Vec<Eventalign>) -> eyre::Result<()> {
        let tx = self.0.transaction()?;
        let mut stmt = tx.prepare("INSERT INTO data (kmer, sample) VALUES (?1, ?2)")?;
        for eventalign in es.into_iter() {
            log::debug!("Processing {:?}", eventalign.metadata());
            for signal in eventalign.signal_iter() {
                let kmer = signal.kmer();
                for sample in signal.samples() {
                    if !(&40.0..=&170.0).contains(&sample) {
                        log::warn!("Uncharacteristic signal measurement {sample}");
                        continue;
                    }
                    if sample.is_finite() {
                        stmt.execute((kmer, sample))?;
                    }
                }
            }
        }
        stmt.finalize()?;

        tx.commit()?;
        Ok(())
    }

    fn get_kmer_samples(&self, kmer: &str, n_samples: usize) -> eyre::Result<Vec<f64>> {
        let mut stmt = self
            .0
            .prepare("SELECT sample FROM data where kmer = :kmer ORDER BY RANDOM() LIMIT :n")?;
        let rows = stmt.query_map(named_params! {":kmer": kmer, ":n": n_samples}, |row| {
            row.get::<usize, f64>(0)
        })?;
        let mut samples = Vec::new();
        for sample in rows {
            samples.push(sample?)
        }
        Ok(samples)
    }
}

#[cfg(test)]
mod test {
    use assert_fs::TempDir;

    // use quickcheck::quickcheck;
    use super::*;
    use crate::arrow::Signal;

    #[test]
    fn test_all_kmers() {
        let kmers = all_kmers();
        assert_eq!(kmers.len(), 4096);
    }

    #[test]
    fn test_db_no_kmer() {
        let tmp_dir = TempDir::new().unwrap();
        let db_path = tmp_dir.join("test.db");
        let mut db = Db::open(db_path).expect("Failed to open database file");
        let eventalign = Eventalign::default();
        db.add_reads(vec![eventalign]).expect("Unable to add read");
        let samples = db
            .get_kmer_samples("ABCDEF", 5000)
            .expect("Unable to get samples");
        assert!(samples.is_empty());
    }

    #[test]
    fn test_db() {
        let tmp_dir = TempDir::new().unwrap();
        let db_path = tmp_dir.join("test.db");
        let test_cases = vec![
            ("AAAAAA", vec![100.0; 3], true),
            ("GGGGGG", vec![20.0; 4], false),
            ("CCCCCC", vec![300.0; 2], false),
        ];
        let mut db = Db::open(db_path).expect("Failed to open database file");
        let signal_data = test_cases
            .iter()
            .enumerate()
            .map(|(i, (k, xs, _))| Signal::new(i as u64, k.to_string(), 1.0, 0.5, xs.clone()))
            .collect::<Vec<_>>();
        let mut eventalign = Eventalign::default();
        *eventalign.signal_data_mut() = signal_data;
        db.add_reads(vec![eventalign]).expect("Unable to add read");

        for (k, xs, unfiltered) in test_cases.into_iter() {
            let err_msg = format!("Unable to retrieve kmer values for {k}");
            let samples = db.get_kmer_samples(k, 5000).expect(&err_msg);
            if unfiltered {
                assert_eq!(samples, xs);
            } else {
                assert!(samples.is_empty(), "{k}");
            }
        }
    }

    #[test]
    fn test_train() {
        let cases = vec![
            1.0,
            2.0,
            3.0,
            4.0,
            1.2,
            1.3,
            2.3,
            2.2,
            3.3,
            f64::NAN,
            f64::NEG_INFINITY,
            f64::INFINITY,
        ];
        let opts = TrainOptions::default();
        let xs = opts.train_gmm(cases);
        assert!(xs.is_some(), "first");

        let case = Vec::new();
        let xs = opts.train_gmm(case);
        assert!(xs.is_none(), "empty");

        let case = vec![f64::NAN; 10];
        let xs = opts.train_gmm(case);
        assert!(xs.is_none(), "NAN");

        let case = vec![f64::INFINITY; 10];
        let xs = opts.train_gmm(case);
        assert!(xs.is_none(), "INFINITY");

        let case = vec![f64::NEG_INFINITY; 100];
        let xs = opts.train_gmm(case);
        assert!(xs.is_none(), "NEG_INFINITY");

        let case = vec![100.0, 100.0, 0.0, -0.0];
        let xs = opts.train_gmm(case);
        assert!(xs.is_some(), "hundreds");

        let case = vec![1.5028554895297472e233, 0.0];
        let xs = opts.train_gmm(case);
        assert!(xs.is_none(), "overflow");
    }

    // quickcheck! {
    //     fn valid_prop(xs: Vec<f32>) -> bool {
    //         let opts = TrainOptions::default();
    //         let xs = xs.into_iter().map(|x| x as f64).collect();
    //         let res = opts.train_gmm(xs);
    //         !matches!(res, Some(Err(GmmError::MinMaxError(_))))
    //     }
    // }
}
