use std::{
    fs::File,
    io::{BufWriter, Write},
    path::Path,
};

use anyhow::Result;

use crate::arrow::{load_apply, Eventalign, Metadata};

fn to_bed_line(metadata: &Metadata, chunk_idx: usize, rec_idx: usize) -> String {
    let chrom = metadata.chrom();
    let start = metadata.start();
    let stop = start + metadata.length();
    let read_name = metadata.name();
    let strand = metadata.strand().as_str();
    format!(
        "{chrom}\t{start}\t{stop}\t{read_name}\t0\t{strand}\t{0}\t{1}",
        chunk_idx, rec_idx
    )
}

pub(crate) fn index<P>(filepath: P) -> Result<()>
where
    P: AsRef<Path>,
{
    let file = File::open(&filepath)?;
    let output_filepath = filepath.as_ref().join(".idx.bed");
    let writer = File::create(output_filepath)?;
    let mut writer = BufWriter::new(writer);

    let mut chunk_idx = 0usize;
    load_apply(file, |chunk: Vec<Eventalign>| {
        for (rec_idx, event) in chunk.into_iter().enumerate() {
            let idx_rec = to_bed_line(event.metadata(), chunk_idx, rec_idx);
            writeln!(writer, "{}", idx_rec)?;
        }
        chunk_idx += 1;
        Ok(())
    })?;
    writer.flush()?;
    Ok(())
}
