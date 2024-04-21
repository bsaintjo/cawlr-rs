use arrow2::datatypes::{Field, Schema};
use arrow2_convert::{field::ArrowField, ArrowDeserialize, ArrowField, ArrowSerialize};
use serde::{Deserialize, Serialize};

use super::{
    eventalign::Eventalign,
    metadata::{Metadata, MetadataExt},
};

/// Represents a single read scored by cawlr score
#[derive(
    Debug, Clone, ArrowField, ArrowDeserialize, ArrowSerialize, Default, Serialize, Deserialize,
)]
pub struct ScoredRead {
    pub metadata: Metadata,
    pub scores: Vec<Score>,
}

impl ScoredRead {
    pub fn new(metadata: Metadata, scores: Vec<Score>) -> Self {
        ScoredRead { metadata, scores }
    }

    /// Creates new ScoredRead using metadata from Eventalign output
    pub fn from_read_with_scores(eventalign: Eventalign, scores: Vec<Score>) -> Self {
        let metadata = eventalign.metadata;
        ScoredRead::new(metadata, scores)
    }

    /// Schema used for outputing into Arrow file
    pub fn schema() -> Schema {
        let data_type = Self::data_type();
        Schema::from(vec![Field::new("scored", data_type, false)])
    }

    pub fn scores(&self) -> &[Score] {
        &self.scores
    }
}

impl MetadataExt for ScoredRead {
    fn metadata(&self) -> &Metadata {
        &self.metadata
    }
}

#[derive(
    Default, Debug, Clone, ArrowDeserialize, ArrowSerialize, ArrowField, Serialize, Deserialize,
)]
pub struct Score {
    pub pos: u64,
    pub kmer: String,
    pub skipped: bool,
    pub signal_score: Option<f64>,
    // pub skip_score: f64,
    pub score: f64,
}

impl Score {
    pub fn new(
        pos: u64,
        kmer: String,
        skipped: bool,
        signal_score: Option<f64>,
        // skip_score: f64,
        score: f64,
    ) -> Self {
        Self {
            pos,
            kmer,
            skipped,
            signal_score,
            // skip_score,
            score,
        }
    }
}
