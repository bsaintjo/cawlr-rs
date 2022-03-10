use std::collections::HashMap;

use criterion_stats::univariate::{
    kde::{kernel::Gaussian, Bandwidth, Kde},
    Sample,
};

use crate::reads::ScoredRead;

fn run(
    pos_scores: Vec<ScoredRead>,
    neg_scores: Vec<ScoredRead>,
    reads: Vec<ScoredRead>,
    ranks: HashMap<String, f64>,
) {
    let pos_samples = pos_scores
        .into_iter()
        .flat_map(|lr| lr.into_means())
        .collect::<Vec<f64>>();
    let pos_samples = Sample::new(&pos_samples);
    let pos_kde = Kde::new(pos_samples, Gaussian, Bandwidth::Silverman);

    let neg_samples = neg_scores
        .into_iter()
        .flat_map(|lr| lr.into_means())
        .collect::<Vec<f64>>();
    let neg_samples = Sample::new(&neg_samples);
    let neg_kde = Kde::new(neg_samples, Gaussian, Bandwidth::Silverman);
}
