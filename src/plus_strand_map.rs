use std::{collections::hash_map::Entry, path::Path, str::from_utf8};

use anyhow::Result;
use bam::BamReader;
use fnv::FnvHashMap;

#[derive(Default)]
pub(crate) struct PlusStrandMap(FnvHashMap<Vec<u8>, bool>);

impl PlusStrandMap {
    fn new(db: FnvHashMap<Vec<u8>, bool>) -> Self {
        Self(db)
    }

    pub(crate) fn from_bam_file<P: AsRef<Path>>(bam_file: P) -> Result<Self> {
        let mut acc = FnvHashMap::default();
        let reader = BamReader::from_path(bam_file, 2u16)?;
        for record in reader {
            let record = record?;
            let read_name = record.name();

            log::info!("ReadName from bam: {:?}", from_utf8(read_name));

            let plus_stranded = !record.flag().is_reverse_strand();
            match acc.entry(read_name.to_owned()) {
                Entry::Occupied(mut entry) => {
                    let old_stranded = entry.insert(plus_stranded);
                    if old_stranded != plus_stranded {
                        log::warn!("Multimapped read has strand swap");
                    }
                }
                Entry::Vacant(entry) => {
                    entry.insert(plus_stranded);
                }
            }
        }
        Ok(PlusStrandMap::new(acc))
    }

    pub(crate) fn get<B>(&self, read_id: B) -> Option<bool>
    where
        B: AsRef<[u8]>,
    {
        let read_id = read_id.as_ref();
        self.0.get(read_id).cloned()
    }

    pub(crate) fn insert<K: Into<Vec<u8>>>(&mut self, k: K, v: bool) -> Option<bool> {
        self.0.insert(k.into(), v)
    }
}


#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_from_bam_file() {
        let filepath = "extra/single_read.bam";
        let psmap = PlusStrandMap::from_bam_file(filepath).unwrap();
        let read_id: &[u8] = b"20d1aac0-29de-43ae-a0ef-aa8a6766eb70";
        assert!(psmap.0.contains_key(read_id));
        assert_eq!(psmap.get(read_id), Some(true));
    }

    #[test]
    fn test_from_bam_file_neg_strand() {
        let filepath = "extra/pos_control.bam";
        let psmap = PlusStrandMap::from_bam_file(filepath).unwrap();
        let read_id: &[u8] = b"ca10c9e3-61d4-439b-abb3-078767d19f8c";
        assert!(psmap.0.contains_key(read_id));
        assert_eq!(psmap.get(read_id), Some(false));
    }
}
