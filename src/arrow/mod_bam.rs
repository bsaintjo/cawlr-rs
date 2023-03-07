//! Provides abstraction for dealing with BAM files containing modification data
//!
//! Current uses bam, but should be switched over to rust-htslib or
//! noodles
use std::{
    fs::File,
    io::{self},
    path::Path,
};

use bam::{record::tags::TagValue, BamReader};

use super::{
    metadata::{Metadata, Strand},
    scored_read::{Score, ScoredRead},
};

#[derive(thiserror::Error, Debug)]
pub enum ModBamConversionError {
    #[error("Mm or Ml tag not found")]
    NoTags,

    #[error("No scores found for modification motif")]
    NoScores,
}

pub struct ModBamAlignment<'a> {
    rec: bam::Record,
    base_mod: &'a [u8],
    header: &'a bam::Header,
}

impl<'a> ModBamAlignment<'a> {
    /// Return None if the read is unaligned
    fn from_record(rec: bam::Record, base_mod: &'a [u8], header: &'a bam::Header) -> Option<Self> {
        if rec.start() == -1 {
            None
        } else {
            Some(Self {
                rec,
                base_mod,
                header,
            })
        }
    }

    fn as_metadata(&self) -> Metadata {
        let name = std::str::from_utf8(self.rec.name())
            .expect("bam read name utf8")
            .to_string();
        let start: u64 = self.rec.start().try_into().expect("Start to positive only");
        let length = self.rec.query_len() as u64;
        let strand = if self.rec.flag().is_reverse_strand() {
            Strand::plus()
        } else {
            Strand::minus()
        };
        let ref_id = self.rec.ref_id() as u32;
        let chrom = self
            .header
            .reference_name(ref_id)
            .expect("No reference name")
            .to_string();
        Metadata::new(name, chrom, start, length, strand, String::new())
    }

    fn mod_prob_positions(&self) -> Result<ModProbs, ModBamConversionError> {
        let tags = self.rec.tags();
        let Some(TagValue::String(score_pos, _)) = tags.get(b"Mm") else { return Err(ModBamConversionError::NoTags); };
        let ModPosMM { skipped, positions } = ModPosMM::parse_mm_tag(self.base_mod, score_pos)
            .ok_or(ModBamConversionError::NoTags)?;

        let Some(TagValue::IntArray(score_prob_arr)) = tags.get(b"Ml") else { return Err(ModBamConversionError::NoTags);  };
        let probs = score_prob_arr
            .raw()
            .iter()
            .map(|&x| (x as f64) / 256.)
            .collect::<Vec<_>>();
        let probs = probs[skipped..skipped + positions.len()].to_vec();
        Ok(ModProbs {
            probs,
            positions,
            modbam: self,
        })
    }
}

struct ModProbs<'a> {
    probs: Vec<f64>,
    positions: Vec<u64>,
    modbam: &'a ModBamAlignment<'a>,
}

impl<'a> ModProbs<'a> {
    fn into_scores(self) -> Result<Vec<Score>, ModBamConversionError> {
        let mut pos_acc = 0;
        let mut scores = Vec::with_capacity(self.probs.len());

        let start = u64::try_from(self.modbam.rec.start()).unwrap();
        let mod_base = self.modbam.base_mod[0];
        let kmer = String::from_utf8(vec![mod_base]).unwrap();

        let seq = if self.modbam.rec.flag().is_reverse_strand() {
            self.modbam.rec.sequence().rev_compl(..).collect()
        } else {
            self.modbam.rec.sequence().to_vec()
        };
        let seq_positions = seq
            .into_iter()
            .enumerate()
            .filter_map(|b| if b.1 == mod_base { Some(b.0) } else { None })
            .collect::<Vec<_>>();

        for (prob, pos) in self.probs.into_iter().zip(self.positions.into_iter()) {
            pos_acc += pos;
            let abs_pos: u64 = start + (seq_positions[pos_acc as usize] as u64);
            let score = Score::new(abs_pos, kmer.clone(), false, Some(prob), 0.0, prob);
            scores.push(score);
            pos_acc += 1;
        }

        if scores.is_empty() {
            Err(ModBamConversionError::NoScores)
        } else {
            Ok(scores)
        }
    }
}

impl TryFrom<ModBamAlignment<'_>> for ScoredRead {
    type Error = ModBamConversionError;

    fn try_from(modbam: ModBamAlignment) -> Result<Self, Self::Error> {
        let metadata = modbam.as_metadata();
        let scores = modbam.mod_prob_positions()?.into_scores()?;
        Ok(ScoredRead::new(metadata, scores))
    }
}

struct ModPosMM {
    skipped: usize,
    positions: Vec<u64>,
}

impl ModPosMM {
    fn parse_mm_tag(identifier: &[u8], xs: &[u8]) -> Option<Self> {
        let mut skipped = 0;
        let mut positions = None;

        let mod_base = xs.split(|&b| b == b';');

        for section in mod_base {
            match parse_mod_base(identifier, section) {
                None => return None,
                Some(TagMatches::Matched(pos)) => {
                    positions = Some(pos);
                    break;
                }
                Some(TagMatches::Skipped(n)) => skipped += n,
            }
        }
        positions.map(|ps| ModPosMM {
            skipped,
            positions: ps,
        })
    }
}

/// The Mm tag can contain positions of multiple modifications for a single read
/// If we want to focus on a particular modification, in order to get the
/// modification probabilities, we need to parse the Mm tag so we skip the
/// modifications we are not interested in Example
/// C+m,0,3,5;C+Y,1,3,4
/// If we want to focus on C+Y modification, we need to skip the C+m and count
/// how many positions belong to it So the string would be parsed into
/// [Skipped(3), Matched(vec![1,3,4])]
#[derive(Debug, PartialEq)]
enum TagMatches {
    Skipped(usize),
    Matched(Vec<u64>),
}

// TODO handle unwraps gracefully
fn parse_mod_base(base_identify: &[u8], bs: &[u8]) -> Option<TagMatches> {
    let mut base_and_pos = bs.split(|&b| b == b',');
    if base_and_pos.next()? != base_identify {
        return Some(TagMatches::Skipped(base_and_pos.count()));
    }
    Some(TagMatches::Matched(
        base_and_pos
            .map(|p| std::str::from_utf8(p).unwrap().parse::<u64>().unwrap())
            .collect(),
    ))
}

pub struct BamRecords(BamReader<File>);

impl BamRecords {
    pub(crate) fn from_path<P: AsRef<Path>>(path: P) -> io::Result<Self> {
        Ok(Self(bam::BamReader::from_path(path, 4)?))
    }

    pub(crate) fn from_file(file: File) -> io::Result<Self> {
        Ok(Self(bam::BamReader::from_stream(file, 4)?))
    }
}

pub struct ModBamIter {
    records: BamRecords,
    base_mod: Vec<u8>,
}

impl ModBamIter {
    pub fn new<B>(records: BamRecords, base_mod: B) -> Self
    where
        B: Into<Vec<u8>>,
    {
        let base_mod = base_mod.into();
        Self { records, base_mod }
    }

    pub fn next(&mut self) -> Option<io::Result<Option<ModBamAlignment<'_>>>> {
        let Some(res) = self.records.0.next() else { return None; };
        let Ok(rec) = res else { return Some(Err(res.err().unwrap())); };
        let mba = ModBamAlignment::from_record(rec, &self.base_mod, self.records.0.header());
        Some(Ok(mba))
    }
}

struct ModBaseTag {
    fundamental: u8,
    is_top: bool,
    modified_base: Vec<u8>,
}

#[cfg(test)]
pub(crate) mod test {
    use noodles::sam::record::data::field::{Tag, Value};

    use super::*;

    #[test]
    fn test_not_bam_file() {
        let random_path = "extra/not-real.bam";
        let modbam = BamRecords::from_path(random_path);
        assert!(modbam.is_err())
    }

    #[test]
    fn test_modbam_conversion() -> eyre::Result<()> {
        let example = "extra/modbams/MM-double.bam";
        let base_mod = b"C+m".to_vec();
        let mut modbam = BamRecords::from_path(example)?;
        let header = modbam.0.header().clone();
        let rec = modbam.0.next().unwrap().unwrap();
        let aln = ModBamAlignment {
            rec,
            base_mod: &base_mod,
            header: &header,
        };
        let mod_prob_pos = aln.mod_prob_positions()?;
        assert_eq!(mod_prob_pos.positions, vec![1, 3, 0]);
        assert_eq!(
            mod_prob_pos.probs,
            vec![128. / 256., 153. / 256., 179. / 256.]
        );
        Ok(())
    }

    #[test]
    fn test_modbam_conversion_no_mm_ml() -> eyre::Result<()> {
        let example = "extra/neg_control.bam";
        let base_mod = b"C+m".to_vec();
        let modbam = BamRecords::from_path(example)?;
        let mut rec_iter = ModBamIter::new(modbam, base_mod);
        let aln = rec_iter.next().unwrap().unwrap().unwrap();
        assert!(aln.mod_prob_positions().is_err());
        Ok(())
    }

    #[test]
    fn test_parse_mod_base() {
        let example = b"C+mh,5,12,0";
        let res = parse_mod_base(b"C+mh", example);
        assert_eq!(res, Some(TagMatches::Matched(vec![5, 12, 0])))
    }

    #[test]
    fn test_parse_mod_base_fail() {
        let example = b"C-mh,5,12,0";
        let res = parse_mod_base(b"C+mh", example);
        assert_eq!(res, Some(TagMatches::Skipped(3)))
    }

    // Scratchpad for eventual reimplementation with nooodles
    #[test]
    fn test_noodles() {
        let example = "extra/modbams/MM-double.bam";
        let mut reader = File::open(example).map(noodles::bam::Reader::new).unwrap();
        let header = reader.read_header().unwrap().parse().unwrap();
        reader.read_reference_sequences().unwrap();
        let rec = reader.records(&header).next().unwrap().unwrap();
        let data = rec.data();
        let Value::UInt8Array(ref ml) = data
            .get(Tag::try_from(*b"Ml").unwrap())
            .or(data.get(Tag::BaseModificationProbabilities))
            .unwrap() else { panic!("Not [u8]")};
        let Value::String(mm) = data
            .get(Tag::try_from(*b"Mm").unwrap())
            .or(data.get(Tag::BaseModifications))
            .unwrap() else { panic!("Not str")};
    }
}