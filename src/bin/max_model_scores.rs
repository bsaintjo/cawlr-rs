use clap::Parser;
use libcawlr::{
    arrow::io::ModFile,
    score_model::{self, extract_max_samples},
    utils::CawlrIO,
};

#[derive(Parser)]
struct Args {
    #[clap(short, long)]
    input: String,

    #[clap(short, long)]
    output: String,

    #[clap(short, long)]
    tag: String,

    /// Number of bins used to estimate the kernel density estimate
    #[clap(short, long, default_value_t = 10_000)]
    bins: u32,

    /// Number of scores sampled from the input to compute the kernel
    /// density estimate
    #[clap(short, long, default_value_t = 10_000)]
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
