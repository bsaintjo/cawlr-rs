use std::{collections::{HashMap, HashSet}, fs::File};

use bio::io::fasta::IndexedReader;
use criterion_stats::univariate::{
    kde::{kernel::Gaussian, Bandwidth, Kde},
    Sample,
};
use ndarray::{arr2, Array, Array2};

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

fn init_matrix(read: ScoredRead) -> Vec<Vec<Option<f64>>> {
    unimplemented!()
}

fn arr_init_matrix(read: ScoredRead) -> Array2<Option<f64>> {
    let len = read.length();
    let mut matrix = Array::from_elem((147, len), None);
    for x in matrix.column_mut(0) {
        *x = Some(1. / 147.);
    }
    matrix
}

fn update_linker(idx: usize, score: f64, matrix: &mut Array2<Option<f64>>, kde: Kde<f64, Gaussian>) {
    unimplemented!()
}

fn update_nucleosome(idx: usize, score: f64, matrix: &mut Array2<Option<f64>>, kde: Kde<f64, Gaussian>) {
    unimplemented!()
}

fn is_modifiable(read: ScoredRead, valid: &mut HashSet<usize>, genome: &mut IndexedReader<File>) {
    unimplemented!()
}

fn has_scores(read: ScoredRead, valid: &mut HashSet<usize>) {
    unimplemented!()
}