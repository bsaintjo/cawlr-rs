use std::fs::File;

use criterion_stats::univariate::kde::{kernel::Gaussian, Kde};
use rv::misc::linspace;
use serde::{Deserialize, Serialize};
use serde_pickle::from_reader;

use crate::utils::CawlrIO;

#[derive(Serialize, Deserialize)]
pub struct BinnedKde {
    bins: Vec<f64>,
}

impl BinnedKde {
    fn new(bins: Vec<f64>) -> Self {
        Self { bins }
    }

    pub(crate) fn from_kde(n_bins: i32, kde: &Kde<f64, Gaussian>) -> Self {
        // TODO explore using a different linspace implementation, only want positive
        // values
        let mut bins: Vec<f64> = linspace(0., 1., n_bins)
            .into_iter()
            .map(|x| kde.estimate(x))
            .collect();
        let total: f64 = bins.iter().sum();
        // Normalize so area approximately sums to 1
        bins.iter_mut().for_each(|x| *x /= total);
        BinnedKde::new(bins)
    }

    pub(crate) fn pmf_from_score(&self, x: f64) -> f64 {
        let idx = x * (self.bins.len() - 1) as f64;
        let idx = idx.round() as usize;
        self.bins[idx]
    }
}

impl CawlrIO for BinnedKde {
    fn save<P>(&self, filename: P) -> anyhow::Result<()>
        where
            P: AsRef<std::path::Path>,
            Self: Sized {
        let mut file = File::create(filename)?;
        serde_pickle::to_writer(&mut file, &self, Default::default())?;
        Ok(())
    }

    fn load<P>(filename: P) -> anyhow::Result<Self>
        where
            P: AsRef<std::path::Path>,
            Self: Sized {
        let file = File::open(filename)?;
        let bkde = from_reader(file, Default::default())?;
        Ok(bkde)
    }
}

#[cfg(test)]
mod test {
    use criterion_stats::univariate::{kde::Bandwidth, Sample};
    use float_eq::assert_float_eq;
    use rand::{prelude::SmallRng, SeedableRng};
    use rv::{prelude::Beta, traits::Rv};

    use super::*;

    #[test]
    fn test_bkde() {
        let mut rng = SmallRng::seed_from_u64(1234);
        let beta = Beta::new_unchecked(5.0, 5.0);
        let samples: Vec<f64> = beta.sample(100, &mut rng);
        let samples = Sample::new(&samples);
        let kde = Kde::new(samples, Gaussian, Bandwidth::Silverman);
        for n_bins in [1_000, 10_000, 100_000] {
            let bkde = BinnedKde::from_kde(n_bins, &kde);
            assert_eq!(bkde.bins.len(), n_bins as usize);

            // Testing edges, these should not panic
            bkde.pmf_from_score(0.0);
            bkde.pmf_from_score(1.0);
            bkde.pmf_from_score(0.99999);
            bkde.pmf_from_score(0.00001);

            let total: f64 = linspace(0.0, 1.0, 5000).into_iter().sum();

            for x in linspace(0.0, 1.0, 5000) {
                assert_float_eq!(kde.estimate(x) / total, bkde.pmf_from_score(x), abs <= 0.01);
            }
        }
    }
}
