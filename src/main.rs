use std::fs::File;

use clap::{ArgEnum, Parser, Subcommand};
use mimalloc::MiMalloc;

use anyhow::Result;

mod preprocess;
mod score;
mod train;

#[global_allocator]
static GLOBAL: MiMalloc = MiMalloc;

#[derive(Parser, Debug)]
#[clap(author, version, about, long_about=None)]
/// Chromatin accessibility with long reads.
struct Args {
    #[clap(short, long)]
    debug: bool,

    #[clap(subcommand)]
    command: Option<Commands>,
}

#[derive(Subcommand, Debug)]
enum Commands {
    /// Calculates mean per-read per-position and optionally filters data based
    /// on a given region.
    Preprocess {
        #[clap(short, long)]
        /// path to nanopolish eventalign output with samples column
        input: String,

        #[clap(short, long)]
        /// path to output file
        output: String,

        #[clap(short, long)]
        /// output only includes data from this chromosome
        chrom: Option<String>,

        #[clap(long)]
        /// output only includes data that aligns at or after this position, should be set with --chrom
        start: Option<u32>,

        #[clap(long)]
        /// output only includes data that aligns at or before this position, should be set with --chrom
        stop: Option<u32>,

        #[clap(arg_enum, default_value_t=OutputFileType::Parquet)]
        output_filetype: OutputFileType,
    },
    Train {
        #[clap(short, long)]
        /// Parquet file of positive or negative control from cawlr preprocess
        input: String,

        #[clap(short, long)]
        output: String,
    },
    Score {
        #[clap(short, long)]
        input: String,

        #[clap(long)]
        pos_ctrl: String,

        #[clap(long)]
        neg_ctrl: String,
    },
    Sma {
        #[clap(short, long)]
        input: String,

        #[clap(short, long)]
        /// output only includes data from this chromosome
        chrom: String,

        #[clap(long)]
        /// output only includes data that aligns at or after this position, should be set with --chrom
        start: u32,

        #[clap(long)]
        /// output only includes data that aligns at or before this position, should be set with --chrom
        stop: u32,
    },
}

#[derive(Debug, Copy, Clone, PartialEq, Eq, PartialOrd, Ord, ArgEnum)]
enum OutputFileType {
    Arrow,
    Parquet,
}

fn main() -> Result<()> {
    env_logger::init();
    let args = Args::parse();
    match &args.command {
        Some(Commands::Preprocess {
            input,
            output,
            chrom,
            start,
            stop,
            output_filetype,
        }) => {
            log::info!("Preprocess command");
            let nps = preprocess::preprocess(input, chrom, start, stop)?;
            match output_filetype {
                OutputFileType::Parquet => {
                    preprocess::write_records_to_parquet(output, nps)?;
                }
                OutputFileType::Arrow => {
                    preprocess::write_records_to_arrow(output, nps)?;
                }
            }
        }
        Some(Commands::Train { input, output }) => {
            let model_db = train::train(input)?;
            let mut file = File::create(output)?;
            serde_pickle::ser::to_writer(&mut file, &model_db, Default::default())?;
        }

        Some(Commands::Score {
            ..
        }) => {
            unimplemented!()
        }

        Some(Commands::Sma {
            ..
        }) => {
            unimplemented!()
        }

        _ => {
            unimplemented!()
        }
    }
    Ok(())
}
