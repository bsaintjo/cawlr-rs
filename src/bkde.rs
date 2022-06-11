use criterion_stats::univariate::kde::{kernel::Gaussian, Kde};
use rv::misc::linspace;

pub(crate) struct BinnedKde {
    bins: Vec<f64>,
}

impl BinnedKde {
    fn new(bins: Vec<f64>) -> Self {
        Self { bins }
    }

    pub(crate) fn from_kde(n_bins: i32, kde: Kde<f64, Gaussian>) -> Self {
        let bins = linspace(0., 1., n_bins)
            .into_iter()
            .map(|x| kde.estimate(x))
            .collect();
        BinnedKde::new(bins)
    }

    pub(crate) fn pmf_from_score(&self, x: f64) -> f64 {
        let idx = x * (self.bins.len() - 1) as f64;
        let idx = idx.round() as usize;
        self.bins[idx]
    }
}

#[cfg(test)]
mod test {
    use criterion_stats::univariate::{Sample, kde::Bandwidth};
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
        let bkde = BinnedKde::from_kde(1000, kde);
        assert_eq!(bkde.bins.len(), 1000);

        bkde.pmf_from_score(0.0);
        bkde.pmf_from_score(1.0);
        bkde.pmf_from_score(0.99999);
        bkde.pmf_from_score(0.00001);
    }
}
