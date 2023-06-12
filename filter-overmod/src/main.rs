use std::{fs::File, path::PathBuf};

use clap::Parser;
use noodles::sam::record::data::field::{value::Array, Tag, Value};
use rust_htslib::bam::{Header, Read, Reader, Writer};

#[derive(Parser)]
struct Args {
    #[clap(short, long)]
    input: PathBuf,

    #[clap(short, long)]
    output: PathBuf,

    #[clap(long)]
    threshold: f64,

    #[clap(long)]
    min_bases: usize,

    #[clap(long)]
    mod_tag: String,
}

#[derive(Debug)]
struct ModPos {
    label: String,
    positions: Vec<usize>,
    idx: usize,
}

struct ModPosParser {
    s: Vec<String>,
    idx: usize,
}

impl ModPosParser {
    fn new(s: &str) -> Self {
        let s = s.split(';').map(|s| s.to_string()).collect();
        Self { s, idx: 0 }
    }
}

impl Iterator for ModPosParser {
    type Item = ModPos;

    fn next(&mut self) -> Option<Self::Item> {
        if self.s.is_empty() {
            None
        } else {
            let nxt = self.s[0].clone();
            self.s = self.s[1..].to_vec();
            let mut xs = nxt.split(',');
            let label = xs.next().unwrap().to_string();
            let positions = xs.map(|x| x.parse::<usize>().unwrap()).collect::<Vec<_>>();
            let idx = self.idx;
            self.idx += positions.len();
            Some(ModPos {
                label,
                positions,
                idx,
            })
        }
    }
}

struct ModProbs {
    probs: Vec<u8>,
}

impl ModProbs {
    fn new(probs: Vec<u8>) -> Self {
        ModProbs { probs }
    }

    fn pos_probs(modpos: &ModPos) -> &[u8] {
        todo!()
    }
}

fn main2() -> eyre::Result<()> {
    jane_eyre::install()?;

    let args = Args::parse();

    let reader = Reader::from_path(args.input)?;
    let header = Header::from_template(reader.header());
    let writer = Writer::from_path(args.output, &header, rust_htslib::bam::Format::Bam)?;

    println!("Hello, world!");
    Ok(())
}

fn main() -> eyre::Result<()> {
    env_logger::init();
    let args = Args::parse();

    let mut reader = File::open(args.input).map(noodles::bam::Reader::new)?;
    let header = reader.read_header()?;

    let mut writer = File::create(args.output).map(noodles::bam::Writer::new)?;
    writer.write_header(&header)?;

    for result in reader.records(&header) {
        let rec = result?;
        log::debug!("Record: {rec:?}");
        let data = rec.data();
        log::debug!("Data: {data:?}");
        let Some(positions) = data.get(Tag::BaseModifications).or(data.get((*b"Mm").try_into().unwrap())) else {
            log::debug!("No BaseModifications, {:?}", rec.read_name());
            writer.write_record(&header, &rec)?;
            continue;
        };
        log::debug!("positions: {positions:?}");
        let Value::String(pos_str) = positions else { panic!("Not valid")};
        let Some(probs) = data.get(Tag::BaseModificationProbabilities).or(data.get((*b"Ml").try_into().unwrap())) else {
            log::debug!("No BaseModificationProbabilities, {:?}", rec.read_name());
            writer.write_record(&header, &rec)?;
            continue;
        };
        let Value::Array(Array::UInt8(arr)) = probs else { panic!("Invalid datatype for tag")};
        log::debug!("probs: {:?}", arr);

        let modprobs = ModProbs::new(arr.clone());

        let labelled = {
            let mut parser = ModPosParser::new(pos_str.as_ref());
            parser.find(|x| x.label == args.mod_tag)
        };

        log::debug!("modpos: {labelled:?}");

        match labelled {
            Some(mp) => todo!(),
            None => todo!(),
        }

        if labelled.is_none() {
            panic!("No matching label found");
        }
    }
    Ok(())
}
