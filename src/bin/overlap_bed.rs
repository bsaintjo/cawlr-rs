//! This tools helps to filter out specific regions from a bed file, usually
//! output from cawlr sma.

use std::{
    error::Error,
    fs::File,
    io::{BufRead, BufReader, LineWriter, Write},
    path::PathBuf,
};

use clap::Parser;

#[derive(Parser)]
struct Args {
    /// Input bed file, usually from cawlr sma
    #[clap(short, long)]
    input: PathBuf,

    /// Output bed file name
    #[clap(short, long)]
    output: PathBuf,

    /// Chromosome to filter on
    #[clap(long)]
    chrom: String,

    /// Start of region to filter on
    #[clap(long)]
    start: usize,

    /// End of region to filter on
    #[clap(long)]
    end: usize,

    /// Percent of the region the bed line should overlap
    #[clap(long)]
    pct: f32,
}

fn filter_line(s: &str, fchrom: &str, fstart: usize, fend: usize, pct: f32) -> bool {
    let components: Vec<&str> = s.split('\t').collect();
    let chrom = components[0];
    let mut start = components[1].parse::<usize>().unwrap();
    let mut end = components[2].parse::<usize>().unwrap();
    if start < fstart {
        start = fstart;
    }

    if end > fend {
        end = fend;
    }

    if end < start {
        return false;
    }
    let length = (end - start) as f32;
    let total_length = (fend - fstart) as f32;
    let pct_overlap = length / total_length;
    (pct_overlap >= pct) && (fchrom == chrom)
}

fn main() -> Result<(), Box<dyn Error>> {
    let args = Args::parse();

    let input = BufReader::new(File::open(args.input)?);
    let mut lines = input.lines();
    // Skip header
    lines.next();
    let mut output = LineWriter::new(File::create(args.output)?);

    for line in lines {
        let line = line?;
        if filter_line(&line, &args.chrom, args.start, args.end, args.pct) {
            let nbytes = output.write(line.as_bytes())?;
            log::info!("Wrote {nbytes:?} bytes");
        }
    }
    output.flush()?;

    Ok(())
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_filter() {
        let line = "chrI\t100\t200\textra\tstuff\n";
        assert!(filter_line(line, "chrI", 90, 150, 0.5));
        assert!(!filter_line(line, "chrI", 90, 95, 0.5));
        assert!(!filter_line(line, "chrI", 300, 400, 0.5));
        assert!(!filter_line(line, "chrII", 90, 150, 0.5));
    }

    #[test]
    fn test_write() {
        let line = "chrI\t100\t200\n";
        let mut buf: LineWriter<Vec<u8>> = LineWriter::new(Vec::new());
        let n = buf.write(line.as_bytes()).unwrap();
        buf.flush().unwrap();
        assert_eq!(n, 13);
    }
}
