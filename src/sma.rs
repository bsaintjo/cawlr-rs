use std::{collections::VecDeque, fmt::Display, fs::File, io::Write, path::Path};

use anyhow::{Context, Result};
use itertools::Itertools;
use nalgebra::DMatrix;

use crate::{
    arrow::{load_apply, Metadata, MetadataExt, ScoredRead, Strand},
    bkde::BinnedKde,
    motif::Motif,
};

pub struct SmaOptions {
    pos_bkde: BinnedKde,
    neg_bkde: BinnedKde,
    motifs: Vec<Motif>,
    writer: Box<dyn Write>,
}

impl SmaOptions {
    pub fn new(
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
        writeln!(
            &mut self.writer,
            "track name=\"cawlr_sma\" itemRgb=\"on\" visibility=2"
        )?;

        let scores_file = File::open(scores_filepath)?;
        load_apply(scores_file, |reads| {
            for read in reads {
                let matrix = self.run_read(&read)?;
                let states = matrix.backtrack();
                let states_rle = states_to_rle(&states);
                let sma_output = SmaOutput::new(&read, states_rle);
                if sma_output.nuc_starts_lens().0.is_empty() {
                    log::warn!("Empty nuc starts: {:?}", sma_output);
                    continue;
                }
                writeln!(&mut self.writer, "{}", sma_output)?;
            }
            Ok(())
        })
    }

    pub fn run_read(&self, read: &ScoredRead) -> Result<SmaMatrix> {
        let mut matrix = SmaMatrix::from_read(read);
        let scores = read.to_expanded_scores();
        for col_idx in 1..=read.np_length() {
            let col_idx = col_idx as usize;
            // Score will be None if a) No data at that position, or b) Position kmer didn't
            // contain motif of interest
            let score =
                scores[col_idx - 1].filter(|s| self.motifs.iter().any(|m| m.within_kmer(s.kmer())));

            // Start nucleosome vs linker
            let linker_val = matrix.probs().column(col_idx - 1)[0];
            let nuc_val = matrix.probs().column(col_idx - 1)[146];
            let (prev_max, prev_max_idx) = {
                if linker_val >= nuc_val {
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

    // TODO Make initial value 10/157 for linker, 1/157 for nucleosome positions
    pub fn from_read(read: &ScoredRead) -> Self {
        let mut val_dm = DMatrix::from_element(147usize, read.np_length() as usize + 1, f64::MIN);
        val_dm.column_mut(0).fill((1. / 147.0f64).ln());

        let ptr_dm = DMatrix::from_element(147, read.np_length() as usize + 1, None);
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

    /// Use the pointer matrix to infer the most likley sequence of states from
    /// the probability matrix.
    ///
    /// First use the probability matrix to find the most likely state it ended
    /// with. Then use that index walk through the pointer matrix and get
    /// each state and what was the state used to update the probability matrix.
    pub fn backtrack(&self) -> VecDeque<States> {
        use States::*;

        let ncols = self.probs().ncols();
        let mut acc = VecDeque::new();
        let mut idx = self.probs().column(ncols - 1).argmax().0;
        let mut pos = self.probs().ncols() - 1;
        // At pos = 0, the column will be the data intialized from 1/147
        while pos > 0 {
            if idx > 0 {
                acc.push_front(Nucleosome);
            } else {
                acc.push_front(Linker);
            }

            if let Some(new_idx) = self.ptrs().column(pos)[idx] {
                idx = new_idx;
            } else {
                break;
            }

            pos -= 1;
        }
        acc
    }
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
fn states_to_rle(states: &VecDeque<States>) -> Vec<(States, u64)> {
    let init_state = *states
        .get(0)
        .expect("Failed to convert, empty States slice");
    let mut curr = (init_state, 1);
    let mut acc = Vec::new();
    for &state in states.iter().skip(1) {
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

#[derive(Debug)]
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
    fn new<M: MetadataExt>(metadata: &'a M, states_rle: Vec<(States, u64)>) -> Self {
        Self {
            metadata: metadata.metadata(),
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
        let (strand_punc, rgb) = match self.strand() {
            Strand::Plus => ("+", "255,0,0"),
            Strand::Minus => ("-", "0,0,255"),
            Strand::Unknown => (".", "0,0,0"),
        };
        // let chrom_start = self.start_0b();
        let (starts, lengths) = self.nuc_starts_lens();
        let bed_start = starts[0];
        let block_count = starts.len();
        let bed_end = starts[block_count - 1] + lengths[block_count - 1];
        let starts = starts.iter().map(|x| (x - bed_start).to_string()).join(",");
        let lengths = lengths.iter().map(|x| x.to_string()).join(",");
        write!(
            f,
            "{0}\t{bed_start}\t{bed_end}\t{1}\t0\t{strand_punc}\t{bed_start}\t{bed_end}\t{rgb}\t{2}\t{lengths}\t{starts}",
            self.chrom(),
            // self.start_0b(),
            self.name(),
            self.num_nuc(),
        )
    }
}

#[cfg(test)]
mod test {
    use nalgebra::dmatrix;

    use super::*;
    use crate::arrow::Strand;

    #[test]
    fn test_backtrack() {
        use States::*;

        let val_matrix = dmatrix![
            -0.5, -0.1, -0.9, -0.9;
            -0.5, -0.2, -0.3, -0.4;
            -0.5, -0.9, -0.0, -0.0;
        ];
        let ptr_matrix = dmatrix![
            None, Some(0), Some(0), Some(1);
            None, Some(0), Some(0), Some(1);
            None, Some(0), Some(0), Some(1);
        ];
        let sma_matrix = SmaMatrix::new(val_matrix, ptr_matrix);
        assert_eq!(sma_matrix.val_matrix[(0, 0)], -0.5);

        let answer: VecDeque<States> = VecDeque::from([Linker, Nucleosome, Linker]);
        assert_eq!(sma_matrix.backtrack(), answer);
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

        let states_rle = vec![
            (Linker, 5),
            (Nucleosome, 20),
            (Linker, 5),
            (Nucleosome, 10),
            (Linker, 100),
        ];

        let sma_output = SmaOutput::new(&metadata, states_rle);
        let formatted = format!("{}", sma_output);
        let answer = "chrI\t105\t140\ttest\t0\t+\t105\t140\t255,0,0\t2\t20,10\t0,25";
        assert_eq!(formatted, answer);
    }

    #[test]
    fn test_states_to_rle() {
        use States::*;
        let states = VecDeque::from([Linker, Linker, Linker, Nucleosome, Nucleosome, Linker]);
        let rle = states_to_rle(&states);
        assert_eq!(rle, [(Linker, 3), (Nucleosome, 2), (Linker, 1)]);
    }
}
