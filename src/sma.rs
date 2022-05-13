use criterion_stats::univariate::{
    kde::{kernel::Gaussian, Bandwidth, Kde},
    Sample,
};
use fnv::FnvHashMap;
use nalgebra::DMatrix;
use ndarray::{Array, Array2};

use crate::arrow::ScoredRead;

fn run(
    pos_scores: Vec<ScoredRead>,
    neg_scores: Vec<ScoredRead>,
    reads: Vec<ScoredRead>,
    ranks: FnvHashMap<String, f64>,
) {
    let pos_samples = extract_samples(&pos_scores);
    let pos_kde = sample_kde(&pos_samples);

    let neg_samples = extract_samples(&neg_scores);
    let neg_kde = sample_kde(&neg_samples);
}

fn extract_samples(reads: &[ScoredRead]) -> Vec<f64> {
    reads
        .into_iter()
        .flat_map(|lr| {
            lr.scores()
                .iter()
                .map(|score| score.score())
                .collect::<Vec<_>>()
        })
        .collect()
}

fn sample_kde(samples: &[f64]) -> Kde<f64, Gaussian> {
    let samples = Sample::new(samples);
    Kde::new(samples, Gaussian, Bandwidth::Silverman)
}

fn run_read(
    read: ScoredRead,
    pos_kde: &Kde<f64, Gaussian>,
    neg_kde: &Kde<f64, Gaussian>,
    motifs: &[&str],
) -> DMatrix<f64> {
    let mut matrix = init_dmatrix(&read);
    let scores = read.to_expanded_scores();
    for col_idx in 1..=read.length() {
        let col_idx = col_idx as usize;
        // Score will be None if a) No data at that position, or b) Position kmer didn't
        // contain motif of interest
        let score = scores[col_idx - 1].filter(|s| motifs.iter().any(|&m| s.kmer().contains(m)));
        let prev_max = matrix.column(col_idx - 1).max();
        let mut col = matrix.column_mut(col_idx);
        {
            let first: &mut f64 = col.get_mut(0).expect("No values in matrix.");
            let val: f64 = match score {
                Some(score) => prev_max + pos_kde.estimate(score.score()).ln(),
                None => prev_max,
            };
            *first = val;
        }
        for rest in 1..147 {
            let next = col.get_mut(rest).expect("Failed to get column value");
            let val = match score {
                Some(score) => prev_max + neg_kde.estimate(score.score()).ln(),
                None => prev_max,
            };
            *next = val;
        }
    }
    matrix
}

// TODO Make initial value 10/157 for linker, 1/157 for nucleosome positions
fn init_dmatrix(read: &ScoredRead) -> DMatrix<f64> {
    let mut dm = DMatrix::from_element(147usize, read.length() as usize + 1, f64::MIN);
    dm.column_mut(0).fill((1. / 147.0f64).ln());
    dm
}

fn arr_init_matrix(read: ScoredRead) -> Array2<Option<f64>> {
    let len = read.length() as usize;
    let mut matrix = Array::from_elem((147, len + 1), None);
    for x in matrix.column_mut(0) {
        *x = Some(1. / 147.);
    }
    matrix
}
