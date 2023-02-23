use std::{fs::File, io, path::Path};

use bam::{record::tags::TagValue, BamReader, Header};

use super::{
    metadata::{Metadata, Strand},
    scored_read::ScoredRead,
};

#[derive(thiserror::Error, Debug)]
pub(crate) enum ModBamConversionError {}

pub(crate) struct ModBamAlignment<'a> {
    rec: bam::Record,
    header: &'a bam::Header,
}

impl<'a> ModBamAlignment<'a> {
    fn from_record(rec: bam::Record, header: &'a bam::Header) -> Option<Self> {
        if rec.start() == -1 {
            None
        } else {
            Some(Self { rec, header })
        }
    }
}

impl TryFrom<ModBamAlignment<'_>> for ScoredRead {
    type Error = ModBamConversionError;

    fn try_from(modbam: ModBamAlignment) -> Result<Self, Self::Error> {
        let name = std::str::from_utf8(modbam.rec.name())
            .expect("bam read name utf8")
            .to_string();
        let start: u64 = modbam
            .rec
            .start()
            .try_into()
            .expect("Start to positive only");
        let length = modbam.rec.query_len() as u64;
        let strand = if modbam.rec.flag().is_reverse_strand() {
            Strand::plus()
        } else {
            Strand::minus()
        };
        let ref_id = modbam.rec.ref_id() as u32;
        let chrom = modbam
            .header
            .reference_name(ref_id)
            .expect("No reference name")
            .to_string();
        let metadata = Metadata::new(name, chrom, start, length, strand, String::new());

        let tags = modbam.rec.tags();
        let Some(TagValue::String(score_pos, _)) = tags.get(b"Mm") else { todo!() };
        let Some(TagValue::IntArray(score_prob_arr)) = tags.get(b"Ml") else { todo!() };
        todo!()
    }
}

struct ModPosMM(Vec<u64>);

impl ModPosMM {
    fn parse_mm_tag(xs: &[u8], start: u64) -> Self {
        todo!()
    }
}

pub(crate) struct BamRecords(BamReader<File>);

impl BamRecords {
    pub(crate) fn from_file<P: AsRef<Path>>(path: P) -> io::Result<Self> {
        Ok(Self(bam::BamReader::from_path(path, 4)?))
    }
}

pub struct ModBamIter {
    records: BamRecords,
    header: bam::Header,
}

impl ModBamIter {
    fn new(records: BamRecords) -> Self {
        let header = records.0.header().clone();
        Self { records, header }
    }

    fn next(&mut self) -> Option<io::Result<ModBamAlignment<'_>>> {
        let Some(res) = self.records.0.next() else { return None; };
        let Ok(rec) = res else { return Some(Err(res.err().unwrap())); };
        let mba = ModBamAlignment {
            rec,
            header: self.records.0.header(),
        };
        Some(Ok(mba))
    }
}

#[cfg(test)]
pub(crate) mod test {
    use super::*;

    #[test]
    fn test_modbam_conversion() -> eyre::Result<()> {
        let example = "example.bam";
        let modbam = BamRecords::from_file(example)?;
        let mut rec_iter = ModBamIter::new(modbam);
        while let Some(Ok(aln)) = rec_iter.next() {
            let read: ScoredRead = aln.try_into()?;
        }
        Ok(())
    }
}
