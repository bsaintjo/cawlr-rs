pub mod arrow_utils;
pub mod eventalign;
pub mod metadata;
pub mod scored_read;
pub mod signal;
mod mod_bam;
pub mod io;

#[cfg(test)]
mod test {
    use std::io::Cursor;

    use arrow2_convert::deserialize::TryIntoCollection;
    use bio::io::fasta::IndexedReader;

    use super::{
        arrow_utils::{load, save, wrap_writer},
        eventalign::Eventalign,
        metadata::{Metadata, Strand},
        signal::Signal,
    };

    #[test]
    fn test_single_read() {}

    #[test]
    fn test_round_trip() {
        let metadata = Metadata::new(
            "abc".to_string(),
            "chrI".to_string(),
            0u64,
            100u64,
            Strand::plus(),
            String::new(),
        );

        let signal = Signal::new(1u64, "AAAAAA".to_string(), 80.0f64, 0.01f64, Vec::new());

        let eventalign = Eventalign::new(metadata, vec![signal]);
        let x = [eventalign.clone(), eventalign];

        let schema = Eventalign::schema();

        let file = vec![];
        let mut writer = wrap_writer(file, &schema).unwrap();
        save(&mut writer, &x).unwrap();
        writer.finish().unwrap();

        let reader = writer.into_inner();
        let reader = Cursor::new(reader);

        let filereader = load(reader).unwrap();
        for chunks in filereader.flatten() {
            for row in chunks.into_arrays().into_iter() {
                let _: Vec<Eventalign> = row.try_into_collection().unwrap();
            }
        }
    }

    #[allow(clippy::read_zero_byte_vec)]
    #[test]
    fn test_fasta_reader_start() {
        let genome = "extra/sacCer3.fa";
        let mut genome = IndexedReader::from_file(&genome).unwrap();
        let chrom = "chrI";
        let start = 71071;
        let stop = start + 6;
        genome.fetch(chrom, start, stop).unwrap();
        let mut seq = Vec::new();
        genome.read(&mut seq).unwrap();

        assert_eq!(b"GCAAGC", seq.as_slice());
    }
}
