use std::collections::HashMap;

use criterion_stats::univariate::{Sample, kde::{Kde, kernel::Gaussian, Bandwidth}};
use linfa::Float;

use crate::score::ScoredRecord;

fn run(pos_scores: Vec<ScoredRecord>, neg_scores: Vec<ScoredRecord>, reads: HashMap<String, ScoredRecord>) {
    let pos_samples = pos_scores
        .into_iter()
        .map(|snpr| snpr.get_score())
        .collect::<Vec<f64>>();
    let pos_samples = Sample::new(&pos_samples);
    let pos_kde = Kde::new(pos_samples, Gaussian, Bandwidth::Silverman);

    let neg_samples = neg_scores
        .into_iter()
        .map(|snpr| snpr.get_score())
        .collect::<Vec<f64>>();
    let neg_samples = Sample::new(&neg_samples);
    let neg_kde = Kde::new(neg_samples, Gaussian, Bandwidth::Silverman);

    
}
