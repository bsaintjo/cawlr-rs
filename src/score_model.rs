use std::io::{Read, Seek};

use criterion_stats::univariate::{
    kde::{kernel::Gaussian, Bandwidth, Kde},
    Sample,
};
use eyre::Result;
use rand::{rngs::SmallRng, seq::SliceRandom, SeedableRng};

use crate::{
    arrow::{load_apply, ScoredRead},
    bkde::BinnedKde,
};

pub struct Options {
    samples: usize,
    bins: u32,
    rng: SmallRng,
}

impl Default for Options {
    fn default() -> Self {
        let rng = SmallRng::seed_from_u64(2456);
        let n_samples = 10_000;
        let n_bins = 10_000;
        Options::new(n_samples, n_bins, rng)
    }
}

impl Options {
    fn new(n_samples: usize, n_bins: u32, rng: SmallRng) -> Self {
        Self {
            samples: n_samples,
            bins: n_bins,
            rng,
        }
    }

    pub fn bins(&mut self, bins: u32) -> &mut Self {
        self.bins = bins;
        self
    }

    pub fn samples(&mut self, samples: usize) -> &mut Self {
        self.samples = samples;
        self
    }

    pub fn run<R>(&mut self, reader: R) -> Result<BinnedKde>
    where
        R: Read + Seek,
    {
        let scores = extract_samples_from_reader(reader)?;
        let scores: Vec<f64> = scores
            .choose_multiple(&mut self.rng, self.samples)
            .cloned()
            .collect();
        let kde = sample_kde(&scores)?;
        let bkde = BinnedKde::from_kde(self.bins as i32, &kde);
        Ok(bkde)
    }
}

fn sample_kde(samples: &[f64]) -> Result<Kde<f64, Gaussian>> {
    if samples.is_empty() {
        eyre::bail!("Score file does not contain any values.");
    }
    let samples = Sample::new(samples);
    Ok(Kde::new(samples, Gaussian, Bandwidth::Silverman))
}

pub fn extract_samples_from_reader<R>(reader: R) -> Result<Vec<f64>>
where
    R: Read + Seek,
{
    let mut scores = Vec::new();
    load_apply(reader, |reads: Vec<ScoredRead>| {
        let mut samples = extract_samples(&reads);
        scores.append(&mut samples);
        Ok(())
    })?;
    Ok(scores)
}

// TODO Use full score instead of signal score
pub fn extract_samples(reads: &[ScoredRead]) -> Vec<f64> {
    reads
        .iter()
        .flat_map(|lr| {
            lr.scores()
                .iter()
                .flat_map(|score| score.signal_score())
                .filter(|x| !x.is_nan())
                .copied()
                .collect::<Vec<_>>()
        })
        .collect()
}
