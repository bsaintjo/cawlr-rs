use std::{
    collections::VecDeque,
    fmt::Display,
    fs::File,
    io::{BufWriter, Write},
    path::Path,
};

use eyre::Result;
use itertools::Itertools;
use nalgebra::DMatrix;

use crate::{
    arrow::{load_apply, Metadata, MetadataExt, ScoredRead, Strand},
    bkde::BinnedKde,
    motif::Motif,
    utils::CawlrIO,
};

fn make_scoring_vec(read: &ScoredRead) -> Vec<f64> {
    let mut calling_vec = Vec::new();
    (0..=(read.end_1b_excl() - read.start_0b() + 1)).for_each(|_| calling_vec.push(-1.0));
    (0..read.scores().len()).for_each(|i| {
        let idx = read.scores()[i].pos() - read.start_0b() + 1;
        calling_vec[idx as usize] = read.scores()[i].score();
    });
    calling_vec
}

fn sma<W: Write>(
    writer: &mut W,
    pos_scores: &BinnedKde,
    neg_scores: &BinnedKde,
    read: &ScoredRead,
) -> Result<()> {
    let calling_vec = make_scoring_vec(read);
    let base_num = read.end_1b_excl() - read.start_0b() + 1;

    // Build matrix
    let mut prob_mat = Vec::new();
    (0..base_num + 1).for_each(|_| prob_mat.push([0.0; 148]));

    let mut ptr_mat = Vec::new();
    (0..base_num + 1).for_each(|_| ptr_mat.push([-1isize; 148]));

    // Initialisation
    let initial_rate: f64 = 1. / 148.;
    let log_initial_rate = initial_rate.ln();

    (0..148).for_each(|j| {
        prob_mat[1][j] = log_initial_rate;
        ptr_mat[1][j] = 0;
    });

    // Recursion
    for i in 2..=base_num {
        let i = i as usize;
        let within_linker;
        let mut back_frm_ncls = 0.0;

        if calling_vec[i] == -1. {
            within_linker = prob_mat[i - 1][0];
            if prob_mat[i - 1][147] != 0.0 {
                back_frm_ncls = prob_mat[i - 1][147];
            }
        } else {
            // let k = (calling_vec[i] * 1000.) as usize;
            // within_linker = EMISSION_PGC_ARRAY[k].ln() + prob_mat[i - 1][0];
            within_linker = pos_scores.pmf_from_score(calling_vec[i]).ln() + prob_mat[i - 1][0];

            if prob_mat[i - 1][147] != 0.0 {
                // back_frm_ncls = EMISSION_PGC_ARRAY[k].ln() + prob_mat[i - 1][147];
                back_frm_ncls =
                    pos_scores.pmf_from_score(calling_vec[i]).ln() + prob_mat[i - 1][147];
            }
        }

        if (back_frm_ncls != 0.0) && (back_frm_ncls > within_linker) {
            prob_mat[i][0] = back_frm_ncls;
            ptr_mat[i][0] = 147;
        } else {
            prob_mat[i][0] = within_linker;
            ptr_mat[i][0] = 0;
        }

        if calling_vec[i] == -1. {
            prob_mat[i][1] = prob_mat[i - 1][0];
        } else {
            // let k = (calling_vec[i] * 1000.) as usize;
            // prob_mat[i][1] = EMISSION_NEG_ARRAY[k].ln() + prob_mat[i - 1][0];
            prob_mat[i][1] = neg_scores.pmf_from_score(calling_vec[i]).ln() + prob_mat[i - 1][0];
        }
        ptr_mat[i][1] = 0;

        for j in 2..=147 {
            if calling_vec[i] == -1. && prob_mat[i - 1][j - 1] != 0.0 {
                prob_mat[i][j] = prob_mat[i - 1][j - 1];
            } else {
                // let k = (calling_vec[i] * 1000.) as usize;
                if prob_mat[i - 1][j - 1] != 0. {
                    // prob_mat[i][j] = EMISSION_NEG_ARRAY[k].ln() + prob_mat[i - 1][j - 1];
                    prob_mat[i][j] =
                        neg_scores.pmf_from_score(calling_vec[i]).ln() + prob_mat[i - 1][j - 1];
                }
            }

            if prob_mat[i][j] != 0. {
                ptr_mat[i][j] = (j - 1) as isize;
            }
        }
    }

    let mut max = f64::NEG_INFINITY;
    let mut max_index = -1;
    for j in 0..148 {
        if prob_mat[base_num as usize][j] > max {
            max = prob_mat[base_num as usize][j];
            max_index = j as isize;
        }
    }

    let mut backtrack_vec = Vec::new();
    for i in (1..=base_num).rev() {
        backtrack_vec.push(max_index);
        max_index = ptr_mat[i as usize][max_index as usize] as isize;
    }

    backtrack_vec.reverse();
    let mut ncls_start = 0;
    let mut ncls_end;
    let shift = read.start_0b() - 1;
    let mut in_nucleosome = false;
    let mut nucs = Vec::new();
    for i in 0..backtrack_vec.len() {
        if backtrack_vec[i] > 0 {
            if !in_nucleosome {
                ncls_start = i + 1 + (shift as usize);
                in_nucleosome = true;
            }
        } else if in_nucleosome {
            ncls_end = i + 1 + (shift as usize);
            nucs.push((ncls_start, ncls_end));
            in_nucleosome = false;
        }
    }
    if in_nucleosome {
        nucs.push((ncls_start, read.end_1b_excl() as usize));
    }

    // Add pseudo block at start if read doesn't start with a nucleosome
    if nucs.is_empty() || nucs[0].0 != read.start_0b() as usize {
        nucs.insert(0, (read.start_0b() as usize, read.start_0b() as usize + 1));
    }

    // Add pseduo block at end if read doesn't end with a nucleosome
    let bend = nucs.last().map(|&(_, b)| b).unwrap();
    if bend != read.end_1b_excl() as usize {
        nucs.push((read.end_1b_excl() as usize - 1, read.end_1b_excl() as usize))
    }

    let n_nucs = nucs.len();
    let (starts, blks): (Vec<_>, Vec<_>) = nucs
        .into_iter()
        .map(|(s, e)| (s - read.start_0b() as usize, (e - s)))
        .unzip();
    writeln!(
        writer,
        "{}\t{}\t{}\t{}\t0\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
        read.chrom(),
        read.start_0b(),
        read.end_1b_excl(),
        read.name(),
        read.strand(),
        read.start_0b(),
        read.end_1b_excl(),
        read.strand().rgb_str(),
        n_nucs,
        blks.into_iter().join(","),
        starts.into_iter().join(","),
    )?;
    Ok(())
}

pub struct SmaOptions {
    track_name: Option<String>,
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
            track_name: None,
            pos_bkde,
            neg_bkde,
            motifs,
            writer,
        }
    }

    pub fn try_new<P: AsRef<Path>>(
        pos_scores_path: P,
        neg_scores_path: P,
        motifs: Vec<Motif>,
        output: P,
    ) -> Result<Self> {
        let pos_bkde = BinnedKde::load(pos_scores_path)?;
        let neg_bkde = BinnedKde::load(neg_scores_path)?;
        let writer = BufWriter::new(File::create(output)?);
        let writer = Box::new(writer);
        Ok(SmaOptions::new(pos_bkde, neg_bkde, motifs, writer))
    }

    pub fn track_name<S: Into<String>>(&mut self, track_name: S) -> &mut Self {
        self.track_name = Some(track_name.into());
        self
    }

    pub fn run<P>(mut self, scores_filepath: P) -> Result<()>
    where
        P: AsRef<Path>,
    {
        let track_name = self
            .track_name
            .clone()
            .unwrap_or_else(|| "cawlr_sma".to_string());
        writeln!(
            &mut self.writer,
            "track name=\"{track_name}\" itemRgb=\"on\" visibility=2"
        )?;

        let scores_file = File::open(scores_filepath)?;
        load_apply(scores_file, |reads: Vec<ScoredRead>| {
            for read in reads {
                sma(&mut self.writer, &self.pos_bkde, &self.neg_bkde, &read)?;
                // let matrix = self.run_read(&read)?;
                // let states = matrix.backtrack();
                // let states_rle = states_to_rle(&states);
                // let sma_output = SmaOutput::new(&read, states_rle);
                // if sma_output.nuc_starts_lens().0.is_empty() {
                //     log::warn!("Empty nuc starts: {:?}", sma_output);
                //     continue;
                // }
                // writeln!(&mut self.writer, "{}", sma_output)?;
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
            let first: &mut f64 = col
                .get_mut(0)
                .ok_or_else(|| eyre::eyre!("No values in matrix."))?;
            let val: f64 = match score.and_then(|s| s.signal_score().as_ref()) {
                Some(&signal_score) => prev_max + self.pos_bkde.pmf_from_score(signal_score).ln(),
                None => prev_max,
            };
            *first = val;

            // Within nucleosome
            for rest in 1..147 {
                let prev_idx = rest - 1;
                matrix.ptrs_mut().column_mut(col_idx)[rest] = Some(prev_idx);
                let nuc_prev = matrix.probs().column(col_idx - 1)[prev_idx];
                let mut col = matrix.probs_mut().column_mut(col_idx);
                let next = col
                    .get_mut(rest)
                    .ok_or_else(|| eyre::eyre!("Failed to get column value"))?;
                let val = match score.and_then(|s| s.signal_score().as_ref()) {
                    Some(&signal_score) => {
                        nuc_prev + self.neg_bkde.pmf_from_score(signal_score).ln()
                    }
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
/// # fn main() -> eyre::Result<()> {
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
        let (mut starts, mut lengths) = self.nuc_starts_lens();
        let bed_start = starts[0];
        let bed_end = starts[starts.len() - 1] + lengths[starts.len() - 1];

        if bed_start != 0 {
            starts.insert(0, self.start_0b());
            lengths.insert(0, 1)
        }

        if bed_end != self.end_1b_excl() {
            starts.push(self.end_1b_excl() - 1);
            lengths.push(1);
        }

        let block_count = starts.len();
        let starts = starts
            .iter()
            .map(|x| (x - self.start_0b()).to_string())
            .join(",");
        let lengths = lengths.iter().map(|x| x.to_string()).join(",");
        write!(
            f,
            "{0}\t{1}\t{2}\t{3}\t0\t{strand_punc}\t{bed_start}\t{bed_end}\t{rgb}\t{block_count}\t{lengths}\t{starts}",
            self.chrom(),
            self.start_0b(),
            self.end_1b_excl(),
            self.name(),
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

        let answer: VecDeque<States> = VecDeque::from([Linker, Nucleosome, Nucleosome]);
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
        let answer = "chrI\t100\t150\ttest\t0\t+\t105\t140\t255,0,0\t4\t1,20,10,1\t0,5,30,49";
        pretty_assertions::assert_eq!(formatted, answer);
    }

    #[test]
    fn test_states_to_rle() {
        use States::*;
        let states = VecDeque::from([Linker, Linker, Linker, Nucleosome, Nucleosome, Linker]);
        let rle = states_to_rle(&states);
        assert_eq!(rle, [(Linker, 3), (Nucleosome, 2), (Linker, 1)]);
    }
}
