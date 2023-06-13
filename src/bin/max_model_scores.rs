use clap::Parser;
use libcawlr::{arrow::io::ModFile, score_model::{self, extract_max_samples}, utils::CawlrIO};

#[derive(Parser)]
struct Args {
    #[clap(short, long)]
    input: String,

    #[clap(short, long)]
    output: String,

    #[clap(short, long)]
    tag: String,

    #[clap(short, long)]
    bins: u32,

    #[clap(short, long)]
    samples: usize,
}

fn main() -> eyre::Result<()> {
    let args = Args::parse();
    let mod_file = ModFile::open_path(args.input, Some(args.tag))?;
    let bkde = score_model::Options::default()
        .bins(args.bins)
        .samples(args.samples)
        .run_modfile_with(mod_file, extract_max_samples)?;
    bkde.save_as(args.output)?;
    Ok(())
}
