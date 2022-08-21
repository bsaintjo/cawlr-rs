use std::{collections::VecDeque, fmt::Display, fs::File, io::Write, path::Path};

use anyhow::{Context, Result};
use criterion_stats::univariate::{
    kde::{kernel::Gaussian, Bandwidth, Kde},
    Sample,
};
use itertools::Itertools;
use nalgebra::DMatrix;
use rand::{
    prelude::{SliceRandom, SmallRng},
    SeedableRng,
};

use crate::{
    arrow::{load_apply, Metadata, MetadataExt, ScoredRead},
    bkde::BinnedKde,
    motif,
    motif::Motif,
    utils,
};

pub struct Builder<P> {
    pos_score_file: P,
    neg_score_file: P,
    output_file: Option<P>,
    motifs: Vec<Motif>,
    seed: u64,
    kde_samples: usize,
}

impl<P> Builder<P>
where
    P: AsRef<Path>,
{
    pub fn new(pos_ctrl_file: P, neg_ctrl_file: P) -> Self {
        let motifs = motif::all_bases();
        Builder {
            pos_score_file: pos_ctrl_file,
            neg_score_file: neg_ctrl_file,
            output_file: None,
            motifs,
            seed: 2456,
            kde_samples: 50_000_usize,
        }
    }

    pub fn seed(&mut self, seed: u64) -> &mut Self {
        self.seed = seed;
        self
    }

    pub fn kde_samples(&mut self, kde_samples: usize) -> &mut Self {
        self.kde_samples = kde_samples;
        self
    }

    pub fn motifs(&mut self, motifs: Vec<Motif>) -> &mut Self {
        self.motifs = motifs;
        self
    }

    pub fn try_motifs(&mut self, motifs: Option<Vec<Motif>>) -> &mut Self {
        if let Some(motifs) = motifs {
            self.motifs = motifs;
        }
        self
    }

    pub fn output_file(&mut self, output_file: P) -> &mut Self {
        self.output_file = Some(output_file);
        self
    }

    pub fn try_output_file(&mut self, output_file: Option<P>) -> &mut Self {
        self.output_file = output_file;
        self
    }

    pub fn build(self) -> Result<SmaOptions> {
        let mut rng = SmallRng::seed_from_u64(self.seed);
        let pos_bkde = score_file_to_bkde(self.kde_samples, &mut rng, self.pos_score_file)?;
        let neg_bkde = score_file_to_bkde(self.kde_samples, &mut rng, self.neg_score_file)?;
        let writer = utils::stdout_or_file(self.output_file)?;

        Ok(SmaOptions::new(pos_bkde, neg_bkde, self.motifs, writer))
    }
}

pub fn score_file_to_bkde<P>(
    kde_samples: usize,
    rng: &mut SmallRng,
    filepath: P,
) -> Result<BinnedKde>
where
    P: AsRef<Path>,
{
    let scores_file = File::open(filepath)?;
    let scores = extract_samples_from_file(scores_file)?;
    let scores: Vec<f64> = scores.choose_multiple(rng, kde_samples).cloned().collect();
    let kde = sample_kde(&scores)?;
    let bkde = BinnedKde::from_kde(1000, kde);
    Ok(bkde)
}

pub struct SmaOptions {
    pos_bkde: BinnedKde,
    neg_bkde: BinnedKde,
    motifs: Vec<Motif>,
    writer: Box<dyn Write>,
}

impl SmaOptions {
    fn new(
        pos_bkde: BinnedKde,
        neg_bkde: BinnedKde,
        motifs: Vec<Motif>,
        writer: Box<dyn Write>,
    ) -> Self {
        Self {
            pos_bkde,
            neg_bkde,
            motifs,
            writer,
        }
    }

    pub fn run<P>(mut self, scores_filepath: P) -> Result<()>
    where
        P: AsRef<Path>,
    {
        let scores_file = File::open(scores_filepath)?;
        load_apply(scores_file, |reads| {
            for read in reads {
                let matrix = self.run_read(&read)?;
                let states = matrix.backtrace().into_iter().collect::<Vec<_>>();
                let states_rle = states_to_rle(&states);
                let sma_output = SmaOutput::new(read.metadata(), states_rle);
                writeln!(&mut self.writer, "{}", sma_output)?;
            }
            Ok(())
        })
    }

    pub fn run_read(&self, read: &ScoredRead) -> Result<SmaMatrix> {
        let mut matrix = SmaMatrix::from_read(read);
        let scores = read.to_expanded_scores();
        for col_idx in 1..=read.length() {
            let col_idx = col_idx as usize;
            // Score will be None if a) No data at that position, or b) Position kmer didn't
            // contain motif of interest
            let score =
                scores[col_idx - 1].filter(|s| self.motifs.iter().any(|m| m.within_kmer(s.kmer())));

            // Start nucleosome vs linker
            let linker_val = matrix.probs().column(col_idx - 1)[0];
            let nuc_val = matrix.probs().column(col_idx - 1)[146];
            let (prev_max, prev_max_idx) = {
                if linker_val > nuc_val {
                    (linker_val, 0usize)
                } else {
                    (nuc_val, 146usize)
                }
            };
            matrix.ptrs_mut().column_mut(col_idx)[0] = Some(prev_max_idx);
            let mut col = matrix.probs_mut().column_mut(col_idx);
            let first: &mut f64 = col.get_mut(0).context("No values in matrix.")?;
            let val: f64 = match score {
                Some(score) => prev_max + self.pos_bkde.pmf_from_score(score.score()).ln(),
                None => prev_max,
            };
            *first = val;

            // Within nucleosome
            for rest in 1..147 {
                let prev_idx = rest - 1;
                matrix.ptrs_mut().column_mut(col_idx)[rest] = Some(prev_idx);
                let nuc_prev = matrix.probs().column(col_idx - 1)[prev_idx];
                let mut col = matrix.probs_mut().column_mut(col_idx);
                let next = col.get_mut(rest).context("Failed to get column value")?;
                let val = match score {
                    Some(score) => nuc_prev + self.neg_bkde.pmf_from_score(score.score()).ln(),
                    None => prev_max,
                };
                *next = val;
            }
        }
        Ok(matrix)
    }
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

pub fn extract_samples_from_file(file: File) -> Result<Vec<f64>> {
    let mut scores = Vec::new();
    load_apply(file, |reads: Vec<ScoredRead>| {
        let mut samples = extract_samples(&reads);
        scores.append(&mut samples);
        Ok(())
    })?;
    Ok(scores)
}

fn sample_kde(samples: &[f64]) -> Result<Kde<f64, Gaussian>> {
    if samples.is_empty() {
        return Err(anyhow::anyhow!("Score file does not contain any values."));
    }
    let samples = Sample::new(samples);
    Ok(Kde::new(samples, Gaussian, Bandwidth::Silverman))
}

pub struct SmaMatrix {
    val_matrix: DMatrix<f64>,
    ptr_matrix: DMatrix<Option<usize>>,
}

impl SmaMatrix {
    pub fn new(val_matrix: DMatrix<f64>, ptr_matrix: DMatrix<Option<usize>>) -> Self {
        Self {
            val_matrix,
            ptr_matrix,
        }
    }

    pub fn probs(&self) -> &DMatrix<f64> {
        &self.val_matrix
    }

    pub fn from_read(read: &ScoredRead) -> Self {
        let val_dm = init_dmatrix(read);

        let ptr_dm = DMatrix::from_element(147, read.length() as usize + 1, None);
        Self::new(val_dm, ptr_dm)
    }

    pub fn ptrs(&self) -> &DMatrix<Option<usize>> {
        &self.ptr_matrix
    }

    pub fn probs_mut(&mut self) -> &mut DMatrix<f64> {
        &mut self.val_matrix
    }

    pub fn ptrs_mut(&mut self) -> &mut DMatrix<Option<usize>> {
        &mut self.ptr_matrix
    }

    pub fn backtrace(&self) -> VecDeque<States> {
        unimplemented!()
    }
}

// TODO Make initial value 10/157 for linker, 1/157 for nucleosome positions
pub fn init_dmatrix(read: &ScoredRead) -> DMatrix<f64> {
    let mut dm = DMatrix::from_element(147usize, read.length() as usize + 1, f64::MIN);
    dm.column_mut(0).fill((1. / 147.0f64).ln());
    dm
}

#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub enum States {
    Linker,
    Nucleosome,
}

/// Converts a list of States into a Run-length encoding vector of states
///
/// # Examples
/// ```rust, ignore
/// # fn main() -> anyhow::Result<()> {
/// use cawlr::sma::States::*;
/// let states = vec![Linker, Linker, Linker, Nucleosome, Nucleosome, Linker];
/// let rle = states_to_rle(&states);
/// assert_eq!(rle, [(Linker, 3), (Nucleosome, 2), (Linker, 1)]);
/// # Ok(())
/// # }
/// ```
///
/// # Panics
/// Will panic if states is empty
fn states_to_rle(states: &[States]) -> Vec<(States, u64)> {
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
    acc
}

/// Converts a slice of States into a more readable format, to make it easier
/// for downstream parsing in python
pub fn states_to_readable(states: &[States]) -> Result<Vec<(String, u64)>> {
    if states.is_empty() {
        return Err(anyhow::anyhow!("Empty list of states"));
    }

    let acc = states_to_rle(states);
    let readable = acc
        .into_iter()
        .map(|(x, y)| match x {
            States::Linker => ("linker".to_string(), y),
            States::Nucleosome => ("nucleosome".to_string(), y),
        })
        .collect();
    Ok(readable)
}

pub fn backtrace(matrix: DMatrix<f64>) -> Vec<States> {
    let mut pos = matrix.ncols() - 1;
    let mut acc = Vec::new();
    let nuc_idx = matrix.ncols() - 1;
    while pos > 0 {
        let linker_val = matrix.column(pos)[0];
        let nuc_val = matrix.column(pos)[nuc_idx];
        if linker_val > nuc_val {
            acc.push(States::Linker);
            pos -= 1;
        } else {
            let n = if pos > 147 { 147 } else { pos };
            acc.append(&mut vec![States::Nucleosome; n]);
            pos -= n;
        }
    }
    acc.reverse();
    acc
}

struct SmaOutput<'a> {
    metadata: &'a Metadata,
    states_rle: Vec<(States, u64)>,
}

impl<'a> MetadataExt for SmaOutput<'a> {
    fn metadata(&self) -> &Metadata {
        self.metadata
    }
}

impl<'a> SmaOutput<'a> {
    fn new(metadata: &'a Metadata, states_rle: Vec<(States, u64)>) -> Self {
        Self {
            metadata,
            states_rle,
        }
    }

    // Count number of runs of States::Nucleosomes
    fn num_nuc(&self) -> usize {
        self.states_rle
            .iter()
            .filter(|x| x.0 == States::Nucleosome)
            .count()
    }

    // Converts the RLE States list into two vectors, one containing the starts of
    // the nucleosome positions, the other of their lengths
    fn nuc_starts_lens(&self) -> (Vec<u64>, Vec<u64>) {
        let mut acc = Vec::new();
        let mut curr = self.start_0b();

        for (state, length) in self.states_rle.iter() {
            if *state == States::Nucleosome {
                acc.push((curr, *length));
            }
            curr += length;
        }

        acc.into_iter().unzip()
    }
}

impl<'a> Display for SmaOutput<'a> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let (starts, lengths) = self.nuc_starts_lens();
        let starts = starts.iter().map(|x| x.to_string()).join(",");
        let lengths = lengths.iter().map(|x| x.to_string()).join(",");
        write!(
            f,
            "{}\t{}\t{}\t{}\t0\t.\t{}\t{}\t0,0,0\t{}\t{}\t{}",
            self.chrom(),
            self.start_0b(),
            self.end_1b_excl(),
            self.name(),
            self.start_0b(),
            self.end_1b_excl(),
            self.num_nuc(),
            starts,
            lengths,
        )
    }
}

#[cfg(test)]
mod test {
    use nalgebra::dmatrix;

    use super::*;
    use crate::arrow::Strand;

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
    fn test_states_to_readable_empty() {
        let acc = Vec::new();
        assert!(states_to_readable(&acc).is_err());
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
        assert_eq!(states_to_readable(&states).unwrap(), answer);
    }

    #[test]
    fn test_sma_output() {
        use States::*;
        let metadata = Metadata::new(
            "test".to_string(),
            "chrI".to_string(),
            100,
            50,
            Strand::plus(),
            "".to_string(),
        );

        let states_rle = vec![(Linker, 5), (Nucleosome, 20), (Linker, 5), (Nucleosome, 10)];

        let sma_output = SmaOutput::new(&metadata, states_rle);
        let formatted = format!("{}", sma_output);
        let answer = "chrI\t100\t150\ttest\t0\t.\t100\t150\t0,0,0\t2\t105,130\t20,10";
        assert_eq!(formatted, answer);
    }

    #[test]
    fn test_states_to_rle() {
        use States::*;
        let states = vec![Linker, Linker, Linker, Nucleosome, Nucleosome, Linker];
        let rle = states_to_rle(&states);
        assert_eq!(rle, [(Linker, 3), (Nucleosome, 2), (Linker, 1)]);
    }
}
