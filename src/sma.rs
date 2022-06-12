use std::{
    fs::File,
    io::{stdout, Write},
    path::Path,
};

use anyhow::Result;
use criterion_stats::univariate::{
    kde::{kernel::Gaussian, Bandwidth, Kde},
    Sample,
};
use nalgebra::DMatrix;
use rand::{
    prelude::{SliceRandom, SmallRng},
    SeedableRng,
};

use crate::{
    arrow::{load_apply, ScoredRead},
    bkde::BinnedKde,
};

fn score_file_to_bkde(
    kde_samples: usize,
    rng: &mut SmallRng,
    filepath: String,
) -> Result<BinnedKde> {
    let scores_file = File::open(filepath)?;
    let scores = extract_samples_from_file(scores_file)?;
    let scores: Vec<f64> = scores
        .choose_multiple(rng, kde_samples)
        .cloned()
        .collect();
    let kde = sample_kde(&scores);
    let bkde = BinnedKde::from_kde(1000, kde);
    Ok(bkde)
}

pub(crate) struct SmaOptions {
    pos_bkde: BinnedKde,
    neg_bkde: BinnedKde,
    motifs: Vec<String>,
}

impl SmaOptions {
    fn new(
        pos_bkde: BinnedKde,
        neg_bkde: BinnedKde,
        motifs: Vec<String>,
    ) -> Self {
        Self {
            pos_bkde,
            neg_bkde,
            motifs,
        }
    }

    pub(crate) fn try_new(
        pos_scores_filepath: String,
        neg_scores_filepath: String,
        motifs: Option<Vec<String>>,
        kde_samples: usize,
        seed: u64,
    ) -> Result<Self> {
        let mut rng = SmallRng::seed_from_u64(seed);
        let pos_bkde = score_file_to_bkde(kde_samples, &mut rng, pos_scores_filepath)?;
        let neg_bkde = score_file_to_bkde(kde_samples, &mut rng, neg_scores_filepath)?;

        let motifs = motifs.unwrap_or_else(|| {
            vec![
                "A".to_string(),
                "T".to_string(),
                "C".to_string(),
                "G".to_string(),
            ]
        });

        Ok(Self::new(pos_bkde, neg_bkde, motifs))
    }

    pub(crate) fn run<P>(self, scores_filepath: P) -> Result<()>
    where
        P: AsRef<Path>,
    {
        let stdout = stdout();
        let mut handle = stdout.lock();
        let scores_file = File::open(scores_filepath)?;
        load_apply(scores_file, |reads| {
            for read in reads {
                let matrix = run_read(&read, &self.pos_bkde, &self.neg_bkde, self.motifs.as_slice());
                let result = backtrace(matrix);
                let result = states_to_readable(&result);
                for (state, size) in result.into_iter() {
                    let chrom = read.metadata().chrom();
                    let start = read.metadata().start();
                    let stop = start + read.metadata().length();
                    let strand = read.metadata().strand();
                    let strand = strand.as_str();
                    let name = read.metadata().name();

                    writeln!(
                        &mut handle,
                        "{chrom}\t{start}\t{stop}\t{name}\t0\t{strand}\t{state}\t{size}"
                    )?;
                }
            }
            Ok(())
        })
    }
}

fn extract_samples(reads: &[ScoredRead]) -> Vec<f64> {
    reads
        .iter()
        .flat_map(|lr| {
            lr.scores()
                .iter()
                .map(|score| score.score())
                .collect::<Vec<_>>()
        })
        .collect()
}

fn extract_samples_from_file(file: File) -> Result<Vec<f64>> {
    let mut scores = Vec::new();
    load_apply(file, |reads: Vec<ScoredRead>| {
        let mut samples = extract_samples(&reads);
        scores.append(&mut samples);
        Ok(())
    })?;
    Ok(scores)
}

fn sample_kde(samples: &[f64]) -> Kde<f64, Gaussian> {
    let samples = Sample::new(samples);
    Kde::new(samples, Gaussian, Bandwidth::Silverman)
}

fn run_read(
    read: &ScoredRead,
    pos_kde: &BinnedKde,
    neg_kde: &BinnedKde,
    motifs: &[String],
) -> DMatrix<f64> {
    let mut matrix = init_dmatrix(read);
    let scores = read.to_expanded_scores();
    for col_idx in 1..=read.length() {
        let col_idx = col_idx as usize;
        // Score will be None if a) No data at that position, or b) Position kmer didn't
        // contain motif of interest
        let score = scores[col_idx - 1].filter(|s| motifs.iter().any(|m| s.kmer().contains(m)));
        let prev_max = matrix.column(col_idx - 1).max();
        let linker_val = matrix.column(col_idx - 1)[0];
        let nuc_val = matrix.column(col_idx - 1)[146];
        {
            let mut col = matrix.column_mut(col_idx);
            let first: &mut f64 = col.get_mut(0).expect("No values in matrix.");
            let val: f64 = match score {
                Some(score) => prev_max + pos_kde.pmf_from_score(score.score()).ln(),
                None => prev_max,
            };
            *first = val;
        }
        for rest in 1..147 {
            let nuc_prev = matrix.column(col_idx - 1)[rest - 1];
            let mut col = matrix.column_mut(col_idx);
            let next = col.get_mut(rest).expect("Failed to get column value");
            let val = match score {
                Some(score) => prev_max + neg_kde.pmf_from_score(score.score()).ln(),
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

#[derive(Copy, Clone, Debug, PartialEq, Eq)]
enum States {
    Linker,
    Nucleosome,
}

/// Converts a slice of States into a more readable format, to make it easier
/// for downstream parsing in python
fn states_to_readable(states: &[States]) -> Vec<(String, u64)> {
    let init_state = *states
        .get(0)
        .expect("Failed to convert, empty States slice");
    let mut curr = (init_state, 1);
    let mut acc = Vec::new();
    for &state in &states[1..] {
        if state == curr.0 {
            curr.1 += 1;
        } else {
            acc.push(curr);
            curr = (state, 1);
        }
    }
    acc.push(curr);

    acc.into_iter()
        .map(|(x, y)| match x {
            States::Linker => ("linker".to_string(), y),
            States::Nucleosome => ("nucleosome".to_string(), y),
        })
        .collect()
}

fn backtrace(matrix: DMatrix<f64>) -> Vec<States> {
    let mut pos = matrix.ncols() - 1;
    let mut acc = Vec::new();
    while pos > 0 {
        let col = matrix.column(pos).argmax().0;
        if col == 0 {
            acc.push(States::Linker);
            pos -= 1;
        } else {
            let n = if pos > col { col } else { 0 };
            acc.append(&mut vec![States::Nucleosome; n]);
            pos -= col;
        }
    }
    acc.reverse();
    acc
}

#[cfg(test)]
mod test {
    use nalgebra::dmatrix;

    use super::*;

    #[test]
    fn test_backtrace() {
        let matrix = dmatrix![0.1, 0.9, 0.9;
                                                                0.2, 0.3, 0.4;
                                                                0.9, 0.0, 0.0];
        assert_eq!(matrix[(0, 1)], 0.9);

        let answer = vec![States::Linker, States::Linker];
        assert_eq!(backtrace(matrix), answer);
    }

    #[test]
    #[should_panic]
    fn test_states_to_readable_empty() {
        let acc = Vec::new();
        states_to_readable(&acc);
    }

    #[test]
    fn test_states_to_readable() {
        let mut states = Vec::new();

        for (x, y) in [(10, 5), (5, 7)] {
            let mut xs = vec![States::Linker; x];
            let mut ys = vec![States::Nucleosome; y];
            states.append(&mut xs);
            states.append(&mut ys);
        }

        let answer = vec![
            ("linker".to_string(), 10),
            ("nucleosome".to_string(), 5),
            ("linker".to_string(), 5),
            ("nucleosome".to_string(), 7),
        ];
        assert_eq!(states_to_readable(&states), answer);
    }
}
