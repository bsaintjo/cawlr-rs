use std::str::FromStr;

use crate::arrow::MetadataExt;

enum FilterError {
    ParseError,
}

struct Region {
    chrom: String,
    start: u64,
    end: u64,
}

impl Region {
    fn new(chrom: String, start: u64, end: u64) -> Self {
        Self { chrom, start, end }
    }

    fn from_bed_line(bed_line: &str) -> Result<Self, FilterError> {
        let spliter: Vec<_> = bed_line.split('\t').collect();
        let chrom = spliter[0].to_string();
        let start = spliter[1].parse().map_err(|_| FilterError::ParseError)?;
        let end = spliter[2].parse().map_err(|_| FilterError::ParseError)?;
        Ok(Region::new(chrom, start, end))
    }

    fn valid<M: MetadataExt>(&self, meta: &M) -> bool {
        (meta.chrom() == self.chrom) && todo!()
    }
}

struct FilterOptions {
    regions: Vec<Region>,
}

impl FilterOptions {
    fn any_valid<M: MetadataExt>(&self, meta: M) -> bool {
        self.regions.iter().any(|r| r.valid(&meta))
    }
}

impl FromStr for Region {
    type Err = FilterError;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let spliterr: Vec<_> = s.split(&[':', '-']).collect();
        let chrom = spliterr[0].to_string();
        let start: u64 = spliterr[1].parse().map_err(|_| FilterError::ParseError)?;
        let end: u64 = spliterr[2].parse().map_err(|_| FilterError::ParseError)?;
        Ok(Region::new(chrom, start, end))
    }
}

fn overlaps(a_start: u64, a_end: u64, b_start: u64, b_end: u64) -> bool {
    (a_start <= b_start) && (a_end >= b_start)
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
    }
}
