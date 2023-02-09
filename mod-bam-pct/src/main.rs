use std::{fs::File, io::Write, path::PathBuf};

use bio_types::sequence::SequenceRead;
use clap::Parser;
use miette::IntoDiagnostic;
use num_format::{Locale, ToFormattedString};
use rust_htslib::bam::{record::Aux, Read, Reader};

#[derive(Parser)]
struct Args {
    #[clap(short, long)]
    input: PathBuf,

    #[clap(short, long)]
    cutoff: f32,

    #[clap(short, long)]
    output: PathBuf,
}

#[derive(serde::Serialize)]
struct Output<'a> {
    read_id: &'a str,
    modifiable_bases: u128,
    modified_bases: u128,
}

fn main() -> miette::Result<()> {
    let args = Args::parse();

    let mut file = File::create(args.output).into_diagnostic()?;
    writeln!(&mut file, "#cutoff {}", args.cutoff).into_diagnostic()?;
    let mut writer = csv::WriterBuilder::default()
        .delimiter(b'\t')
        .from_writer(file);

    let mut reader = Reader::from_path(args.input).into_diagnostic()?;
    let mut total_likely_modified = 0u128;
    let mut total_bases = 0u128;

    for record in reader.records() {
        let record = record.into_diagnostic()?;
        let read_id = std::str::from_utf8(record.name()).into_diagnostic()?;
        let Ok(Aux::ArrayU8(ml)) = record.aux(b"Ml") else { continue };
        let Ok(Aux::String(mm)) = record.aux(b"Mm") else { continue };
        let mm: Vec<u64> = mm
            .split_terminator(&[',', ';'])
            .skip(1)
            .filter(|x| !x.is_empty())
            .map(|x| x.parse::<u64>())
            .collect::<Result<_, _>>()
            .into_diagnostic()?;
        let modified_bases = ml
            .iter()
            .map(|x| (x as f32 / 256.0))
            .filter(|&x| x > args.cutoff)
            .count() as u128;
        let modifiable_bases = (mm.len() as u128) + (mm.iter().sum::<u64>() as u128);

        total_bases += modifiable_bases;
        total_likely_modified += modified_bases;
        let output = Output {
            read_id,
            modifiable_bases,
            modified_bases,
        };
        writer.serialize(output).into_diagnostic()?;
    }

    println!(
        "Total modifiable bases:       {}",
        total_bases.to_formatted_string(&Locale::en)
    );
    println!(
        "Total bases likely modified:  {}",
        total_likely_modified.to_formatted_string(&Locale::en)
    );
    Ok(())
}
