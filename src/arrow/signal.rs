use arrow2_convert::{ArrowDeserialize, ArrowField, ArrowSerialize};
use rv::traits::ContinuousDistr;
use serde::{Deserialize, Serialize};

#[derive(
    Debug,
    Clone,
    Default,
    PartialEq,
    ArrowField,
    ArrowSerialize,
    ArrowDeserialize,
    Serialize,
    Deserialize,
)]
pub struct Signal {
    pub pos: u64,
    pub kmer: String,
    pub signal_mean: f64,
    pub signal_time: f64,
    pub samples: Vec<f64>,
}

impl Signal {
    pub fn new(
        pos: u64,
        kmer: String,
        signal_mean: f64,
        signal_time: f64,
        samples: Vec<f64>,
    ) -> Self {
        Self {
            pos,
            kmer,
            signal_mean,
            signal_time,
            samples,
        }
    }

    pub fn score_lnsum<M, N>(&self, pm: &M, nm: &N) -> Option<(f64, f64)>
    where
        M: ContinuousDistr<f64>,
        N: ContinuousDistr<f64>,
    {
        let mut samples = self
            .samples
            .iter()
            .filter(|&x| (&40.0..=&170.0).contains(&x))
            .peekable();
        // If iterator is empty, we just return None
        samples.peek()?;
        Some(
            samples
                .flat_map(|x| {
                    let likelihood_neg = nm.ln_pdf(x);
                    let likelihood_pos = pm.ln_pdf(x);
                    if likelihood_neg > -10.0 && likelihood_pos > -10.0 {
                        Some((likelihood_pos, likelihood_neg))
                    } else {
                        None
                    }
                })
                .fold((0.0, 0.0), |acc, elem| (acc.0 + elem.0, acc.1 + elem.1)),
        )
    }
}
