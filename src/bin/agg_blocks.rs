use std::path::PathBuf;

use cawlr::agg_blocks::run;
use clap::Parser;

#[derive(Parser)]
struct Args {
    /// Bed file, usually from cawlr sma
    #[clap(short, long)]
    input: PathBuf,

    /// Output tsv containing chromosome, position, frac overlapped
    #[clap(short, long)]
    output: Option<PathBuf>,
}

fn main() -> eyre::Result<()> {
    let args = Args::parse();
    run(&args.input, args.output.as_ref())
}
