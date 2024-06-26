use std::{fmt::Display, str::FromStr};

use thiserror::Error;

use crate::arrow::metadata::MetadataExt;

#[derive(Clone, Debug)]
pub struct Region {
    chrom: String,
    start: u64,
    end: u64,
}

impl Region {
    fn new(chrom: String, start: u64, end: u64) -> Self {
        Self { chrom, start, end }
    }

    pub fn chrom(&self) -> &str {
        &self.chrom
    }

    pub fn start(&self) -> u64 {
        self.start
    }

    pub fn end(&self) -> u64 {
        self.end
    }

    pub fn from_bed_line(bed_line: &str) -> Result<Self, FilterError> {
        if bed_line.is_empty() {
            return Err(FilterError::EmptyRegionError);
        }
        let spliter: Vec<_> = bed_line.split('\t').collect();
        let chrom = spliter[0].to_string();
        let start = spliter[1]
            .parse()
            .map_err(|_| FilterError::StartParseError)?;
        let end = spliter[2].parse().map_err(|_| FilterError::EndParseError)?;
        Ok(Region::new(chrom, start, end))
    }

    pub fn valid<M: MetadataExt + ?Sized>(&self, meta: &M) -> bool {
        (meta.chrom() == self.chrom)
            && overlaps(self.start, self.end, meta.start_0b(), meta.end_1b_excl())
    }
}

fn overlaps(a_start: u64, a_end: u64, b_start: u64, b_end: u64) -> bool {
    ((b_start <= a_start) && (a_start <= b_end)) || // End overlaps
        ((b_start <= a_end) && (a_end <= b_end)) || // Other end overlaps
        ((b_start <= a_start) && (a_end <= b_end)) || // a inside b
        ((a_start <= b_start) && (b_end <= a_end))
}

impl Display for Region {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}:{}-{}", self.chrom, self.start, self.end)
    }
}

impl FromStr for Region {
    type Err = FilterError;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        if s.is_empty() {
            return Err(FilterError::EmptyRegionError);
        }

        let spliterr: Vec<_> = s.split(&[':', '-']).collect();
        if spliterr.len() != 3 {
            return Err(FilterError::ParseError);
        }

        let chrom = spliterr[0].to_string();
        let start: u64 = spliterr[1].parse().map_err(|_| FilterError::ParseError)?;
        let end: u64 = spliterr[2].parse().map_err(|_| FilterError::ParseError)?;
        Ok(Region::new(chrom, start, end))
    }
}

#[derive(Error, Debug)]
pub enum FilterError {
    #[error("Failed to parse chromosome")]
    ChromParseError,
    #[error("Failed to parse chromosome")]
    StartParseError,
    #[error("Failed to parse end position")]
    EndParseError,
    #[error("Empty region")]
    EmptyRegionError,
    #[error("Failed to parse correctly")]
    ParseError,
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_overlap() {
        let a = (10, 15);
        let b = (12, 20);
        assert!(overlaps(a.0, a.1, b.0, b.1));

        let c = (20, 30);
        assert!(!overlaps(a.0, a.1, c.0, c.1));

        let inside_a = (12, 14);
        assert!(overlaps(a.0, a.1, inside_a.0, inside_a.1));

        let outside_a = (9, 16);
        assert!(overlaps(a.0, a.1, outside_a.0, outside_a.1));
    }
}
