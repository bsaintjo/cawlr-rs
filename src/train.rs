use std::{collections::HashMap, fs::File};

use linfa::{traits::Fit, DatasetBase, ParamGuard};
use linfa_clustering::GaussianMixtureModel;
use polars::prelude::{DataFrame, ParquetReader, SerReader};

use anyhow::Result;
use rv::prelude::{Gaussian, Mixture};

type ModelDB = HashMap<String, Mixture<Gaussian>>;

pub(crate) fn save_gmm(output: &str, model_db: ModelDB) -> Result<()> {
    let mut file = File::create(output)?;
    serde_pickle::ser::to_writer(&mut file, &model_db, Default::default())?;
    Ok(())
}

fn train_gmm(df: DataFrame) -> Result<Mixture<Gaussian>> {
    let data = df.column("event_mean")?.f64()?.to_ndarray()?;
    let old_dim = data.shape();
    let new_dim = (old_dim[0], 1);
    let data = data.into_shape(new_dim)?;
    let data = DatasetBase::from(data);

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

pub(crate) fn train(input: &str) -> Result<ModelDB> {
    let file = File::open(input)?;
    let df = ParquetReader::new(file).finish()?;
    let group_idxs = df.groupby(["kmer"])?.groups()?;
    let mut series_iter = group_idxs.iter();
    let kmers = series_iter.next().expect("no kmer series");
    let idxs = series_iter.next().expect("no idxs series");
    let kmers_to_idxs = kmers.utf8()?.into_iter().zip(idxs.list()?.into_iter());
    let mut kmer_to_gmm = HashMap::new();
    for (kmer, idxs) in kmers_to_idxs {
        if let Some(kmer) = kmer {
            if let Some(idxs) = idxs {
                let kdf = df.take(idxs.u32()?)?;
                if kdf.height() < 2 {
                    log::warn!(
                        "Warning: kmer {} doesn't have minimum number of values, skipping...",
                        kmer
                    );
                    continue;
                }
                let result = train_gmm(kdf)?;
                kmer_to_gmm.insert(kmer.to_owned(), result);
            }
        }
    }
    Ok(kmer_to_gmm)
}

#[cfg(test)]
mod tests {
    use polars::df;
    use polars::series::Series;
    use polars::prelude::NamedFrom;
    use rv::traits::ContinuousDistr;

    use super::*;

    #[test]
    fn test_train_gmm() -> Result<()> {
        let df = df!("event_mean" => &[0.1, 0.2, 0.3])?;
        assert_eq!(1 + 1, 2);
        let mm = train_gmm(df)?;
        mm.pdf(&0.2);
        Ok(())
    }
}