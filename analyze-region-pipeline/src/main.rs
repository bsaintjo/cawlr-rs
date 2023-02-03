use clap::{Parser, Subcommand};

mod analyze;
mod preprocess;
mod file_types;

#[derive(Parser)]
pub struct Args {
    #[clap(subcommand)]
    pub command: Commands,
}

#[derive(Subcommand)]
pub enum Commands {
    Analyze(analyze::AnalyzeCmd),
    Preprocess(preprocess::PreprocessCmd),
}
fn main() -> eyre::Result<()> {
    let args = Args::parse();
    match args.command {
        Commands::Analyze(args) => analyze::run(args),
        Commands::Preprocess(_) => todo!(),
    }
}
