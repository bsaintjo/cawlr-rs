use std::{
    error::Error,
    fs::File,
    io::{BufRead, BufReader, LineWriter, Write},
    path::PathBuf,
};

use clap::Parser;

#[derive(Parser)]
struct Args {
    #[clap(short, long)]
    input: PathBuf,

    #[clap(short, long)]
    output: PathBuf,

    #[clap(long)]
    chrom: String,

    #[clap(long)]
    start: usize,

    #[clap(long)]
    end: usize,

    #[clap(long)]
    pct: f32,
}

fn main() -> Result<(), Box<dyn Error>> {
    let args = Args::parse();

    let total_length = (args.end - args.start) as f32;

    let input = BufReader::new(File::open(args.input)?);
    let mut output = LineWriter::new(File::create(args.output)?);

    for line in input.lines() {
        let line = line?;
        let components: Vec<&str> = line.split('\t').collect();
        let chrom = components[0];
        let mut start = components[1].parse::<usize>().unwrap();
        let mut end = components[2].parse::<usize>().unwrap();
        if start < args.start {
            start = args.start;
        }

        if end > args.end {
            end = args.end
        }

        let length = (end - start) as f32;
        let pct_overlap = length / total_length;
        if (pct_overlap >= args.pct) && (chrom == args.chrom) {
            let nbytes = output.write(line.as_bytes());
            log::info!("Wrote {nbytes:?} bytes");
        }
    }
    output.flush()?;

    Ok(())
}
