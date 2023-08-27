use std::{
    fs::File,
    io::{BufWriter, Write},
    path::Path,
};

use eyre::Result;
use itertools::Itertools;

use crate::{
    arrow::{
        arrow_utils::load_apply,
        io::{read_mod_bam_or_arrow, ModFile},
        metadata::MetadataExt,
        scored_read::ScoredRead,
    },
    bkde::BinnedKde,
    motif::Motif,
    utils::CawlrIO,
};

/// Converts all the scores in the read into a vector. Each element is either
/// -1.0 if no value exists, or a score between 0.0 and 1.0.
/// This vector is usually used in the dynamic alignment step later in single
/// molecule analysis.
pub(crate) fn make_scoring_vec(read: &ScoredRead) -> Vec<f64> {
    let mut calling_vec = Vec::new();
    (0..=(read.end_1b_excl() - read.start_0b() + 1)).for_each(|_| calling_vec.push(-1.0));
    let num_scores = read.scores().len();
    log::debug!("N scores: {num_scores}");
    (0..num_scores).for_each(|i| {
        log::debug!("Ith score: {i}");
        log::debug!("Score position: {:?}", read.scores()[i].pos);
        let idx = read.scores()[i].pos - read.start_0b() + 1;
        calling_vec[idx as usize] = read.scores()[i].score;
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
        max_index = ptr_mat[i as usize][max_index as usize];
    }

    backtrack_vec.reverse();
    let mut ncls_start = 0;
    let mut ncls_end;
    let shift = read.start_0b() - 1;
    let mut in_nucleosome = false;
    let mut nucs = Vec::new();
    for (i, bt_idx) in backtrack_vec.into_iter().enumerate() {
        if bt_idx > 0 {
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

/// Loads and stores data used for single molecule analysis.
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

    /// Set bed file track name to track during in the genome browser
    pub fn track_name<S: Into<String>>(&mut self, track_name: S) -> &mut Self {
        self.track_name = Some(track_name.into());
        self
    }

    pub fn run_modfile(mut self, mod_file: ModFile) -> Result<()> {
        let track_name = self
            .track_name
            .clone()
            .unwrap_or_else(|| "cawlr_sma".to_string());
        writeln!(
            &mut self.writer,
            "track name=\"{track_name}\" itemRgb=\"on\" visibility=2"
        )?;

        read_mod_bam_or_arrow(mod_file, |read| {
            if !read.is_unaligned() {
                log::info!("{:?}", read.metadata());
                sma(&mut self.writer, &self.pos_bkde, &self.neg_bkde, &read)?;
            } else {
                log::debug!("Read {} is unaligned, skipping...", read.name())
            }
            Ok(())
        })
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
                log::info!("{:?}", read.metadata());
                sma(&mut self.writer, &self.pos_bkde, &self.neg_bkde, &read)?;
            }
            Ok(())
        })
    }
}
