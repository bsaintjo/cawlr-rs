use std::{
    fs::File,
    io::{BufRead, BufReader, BufWriter, Write},
    path::PathBuf,
};

use cawlr::filter::Region;
use clap::Parser;
use indicatif::{ProgressBar, ProgressStyle};

#[derive(Parser)]
struct Args {
    #[clap(short, long)]
    input: PathBuf,

    #[clap(short, long)]
    region: Region,

    #[clap(short, long)]
    output: PathBuf,

    #[clap(short, long, default_value_t = true)]
    progress: bool,
}

fn run<R: BufRead, W: Write>(region: &Region, reader: R, writer: &mut W) -> eyre::Result<()> {
    for line in reader.lines() {
        let line = line?;
        let mut splitted = line.split('\t').collect::<Vec<_>>();
        let chrom = convert_chrom(splitted[0])?;
        let pos = splitted[1].parse::<u64>()?;
        if (chrom == region.chrom()) && (region.start() <= pos) && (pos <= region.end()) {
            splitted[0] = chrom;
            let result = splitted.join("\t");
            writeln!(writer, "{}", result)?;
        }
    }
    writer.flush()?;
    Ok(())
}

fn convert_chrom(s: &str) -> eyre::Result<&str> {
    let conv = match s {
        "ref|NC_001133|" => "chrI",
        "ref|NC_001134|" => "chrII",
        "ref|NC_001135|" => "chrIII",
        "ref|NC_001136|" => "chrIV",
        "ref|NC_001137|" => "chrV",
        "ref|NC_001138|" => "chrVI",
        "ref|NC_001139|" => "chrVII",
        "ref|NC_001140|" => "chrVIII",
        "ref|NC_001141|" => "chrIX",
        "ref|NC_001142|" => "chrX",
        "ref|NC_001143|" => "chrXI",
        "ref|NC_001144|" => "chrXII",
        "ref|NC_001145|" => "chrXIII",
        "ref|NC_001146|" => "chrXIV",
        "ref|NC_001147|" => "chrXV",
        "ref|NC_001148|" => "chrXVI",
        "ref|NC_001224|" => "chrM",
        _ => return Err(eyre::eyre!("Not valid convertable chrom: {s}")),
    };
    Ok(conv)
}

fn main() -> eyre::Result<()> {
    let args = Args::parse();
    let style = ProgressStyle::with_template("[{elapsed_precise}] [{binary_bytes_per_sec}] {msg}")?;
    let pb = {
        if args.progress {
            ProgressBar::new_spinner().with_style(style)
        } else {
            ProgressBar::hidden()
        }
    };
    // pb.enable_steady_tick(Duration::from_millis(100));
    let reader = BufReader::new(File::open(args.input)?);
    let reader = pb.wrap_read(reader);
    let mut writer = BufWriter::new(File::create(args.output)?);
    run(&args.region, reader, &mut writer)?;
    pb.finish();
    Ok(())
}
