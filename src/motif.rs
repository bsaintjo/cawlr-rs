use std::{collections::HashSet, str::FromStr};

use thiserror::Error;

#[derive(Error, Debug)]
pub enum MotifError {
    #[error("Invalid format, should be in the form [pos]:[motif]")]
    InvalidFormat,
    #[error("Invalid base, should only be ACGT, uppercase only")]
    InvalidBase,
    #[error("Position should be less than the length of the motif given.")]
    PositionOutsideofMotif,
    #[error("Position is one-based.")]
    PositionOneBased,
    #[error("Position must be positive integer")]
    PositionParseFailed,
    #[error("Additional parts not expected. Invalid format")]
    UnexpectedAdditionalFormat,
}

fn valid_motif_bases(motif: &str) -> bool {
    let bases = HashSet::from(['A', 'C', 'G', 'T']);
    !motif.is_empty() && motif.chars().all(|b| bases.contains(&b))
}

#[derive(Debug)]
pub struct Motif {
    motif: String,
    position: usize,
}

impl Motif {
    pub(crate) fn new<S>(motif: S, position: usize) -> Self
    where
        S: Into<String>,
    {
        Self {
            motif: motif.into(),
            position,
        }
    }

    pub fn parse_from_str<T>(string: T) -> Result<Self, MotifError>
    where
        T: AsRef<str>,
    {
        let string = string.as_ref();
        let mut iter = string.split(':');
        let pos = iter
            .next()
            .ok_or(MotifError::InvalidFormat)?
            .parse::<usize>()
            .map_err(|_| MotifError::PositionParseFailed)?;
        let motif = iter.next().ok_or(MotifError::InvalidFormat)?;
        if !valid_motif_bases(motif) {
            Err(MotifError::InvalidBase)
        } else if pos == 0 {
            Err(MotifError::PositionOneBased)
        } else if pos >= motif.len() {
            Err(MotifError::PositionOutsideofMotif)
        } else if iter.next().is_some() {
            Err(MotifError::UnexpectedAdditionalFormat)
        } else {
            Ok(Motif::new(motif, pos))
        }
    }

    pub fn motif(&self) -> &str {
        self.motif.as_ref()
    }

    pub fn len_motif(&self) -> usize {
        self.motif.len()
    }

    pub fn position_1b(&self) -> usize {
        self.position
    }

    pub fn position_0b(&self) -> usize {
        self.position - 1
    }

    // TODO impl std::str::pattern::Pattern when it stabilizes
    pub fn within_kmer(&self, kmer: &str) -> bool {
        kmer.contains(self.motif())
    }

    pub fn surrounding_idxs(&self, pos: usize) -> Vec<usize> {
        unimplemented!()
    }
}

impl FromStr for Motif {
    type Err = MotifError;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        Motif::parse_from_str(s)
    }
}

pub fn all_bases() -> Vec<Motif> {
    vec![
        Motif::new("A", 1),
        Motif::new("C", 1),
        Motif::new("G", 1),
        Motif::new("T", 1),
    ]
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_motif() {
        let m = Motif::parse_from_str("1:AT");
        assert!(m.is_ok());

        let m = Motif::parse_from_str("1:TA");
        assert!(m.is_ok());

        let m = Motif::parse_from_str("0:TA");
        assert!(m.is_err());

        let m = Motif::parse_from_str("TA:1");
        assert!(m.is_err());

        let m = Motif::parse_from_str("1:ZA");
        assert!(m.is_err());

        let m = Motif::parse_from_str("1:ZAhfd");
        assert!(m.is_err());

        let m = Motif::parse_from_str("3:TA");
        assert!(m.is_err());

        let m = Motif::parse_from_str("");
        assert!(m.is_err());

        let m = Motif::parse_from_str("T");
        assert!(m.is_err());

        let m = Motif::parse_from_str("-1:TG");
        assert!(m.is_err());

        let m = Motif::parse_from_str("2.1:TG");
        assert!(m.is_err());

        let m = Motif::parse_from_str("quack:TG");
        assert!(m.is_err());

        let m = Motif::parse_from_str("1:TA:");
        assert!(m.is_err());
    }
}
