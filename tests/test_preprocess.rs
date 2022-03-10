use std::rc::Rc;

use anyhow::Result;
use assert_fs::{
    fixture::{ChildPath, PathChild},
    TempDir,
};
use bio_types::genome::AbstractInterval;
use rust_htslib::bam::{
    ext::BamRecordExtensions, Format, Header, HeaderView, Read, Reader, Record, Writer,
};

fn setup(
    bam_filename: &str,
    input_filename: &str,
    output_filename: &str,
) -> Result<(TempDir, ChildPath)> {
    let mut rec = Record::new();
    let seq = vec![b'T'; 20];
    let qual = vec![255 as u8; seq.len()];
    rec.set(b"abc123", None, &seq, &qual);
    rec.set_tid(0i32);
    rec.set_pos(10i64);
    rec.set_insert_size(100i64);
    let hv = HeaderView::from_bytes(b"@SQ\tSN:chrI\tLN:15072423");
    let header = Header::from_template(&hv);
    let hv = Rc::new(hv);
    rec.set_header(hv);

    let temp_dir = TempDir::new()?;
    let bam_fp = temp_dir.child(bam_filename);
    Writer::from_path(&bam_fp, &header, Format::Bam)?.write(&rec)?;
    Ok((temp_dir, bam_fp))
}

#[test]
fn test_bam() -> Result<()> {
    let bam_file = "test_bam";
    let (_temp_dir, bam_fp) = setup(bam_file, "test_input", "test_output")?;
    let mut record = Record::new();
    let mut bam_reader = Reader::from_path(bam_fp)?;
    bam_reader.read(&mut record);
    assert_eq!(record.contig(), "chrI");
    assert_eq!(10, record.reference_start());
    assert_eq!(20, record.seq_len());
    Ok(())
}
