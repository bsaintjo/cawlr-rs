use arrow2::datatypes::{Field, Schema};
use arrow2_convert::{field::ArrowField, ArrowField};

use super::{
    metadata::{Metadata, MetadataExt},
    signal::Signal,
};

/// Output representing a single read from nanopolish eventalign
#[derive(Debug, Clone, ArrowField, Default, PartialEq)]
pub struct Eventalign {
    pub metadata: Metadata,
    signal_data: Vec<Signal>,
}

impl Eventalign {
    pub fn new(metadata: Metadata, signal_data: Vec<Signal>) -> Self {
        Self {
            metadata,
            signal_data,
        }
    }

    pub fn schema() -> Schema {
        let data_type = Self::data_type();
        Schema::from(vec![Field::new("eventalign", data_type, false)])
    }

    /// Get a mutable reference to the eventalign's signal data.
    pub fn signal_data_mut(&mut self) -> &mut Vec<Signal> {
        &mut self.signal_data
    }

    /// Iterate over Signal data
    pub fn signal_iter(&self) -> impl Iterator<Item = &Signal> {
        self.signal_data.iter()
    }

    pub fn metadata(&self) -> &Metadata {
        &self.metadata
    }
}

impl MetadataExt for Eventalign {
    fn metadata(&self) -> &Metadata {
        &self.metadata
    }
}
