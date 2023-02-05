use std::{
    fs::File,
    io::{self, BufWriter, Read},
    path::PathBuf,
};

use clap::Parser;
use libcawlr::{collapse::CollapseOptions, utils};

#[derive(Parser, Debug)]
pub struct CollapseCmd {
    /// Path to nanopolish eventalign output with samples column, or stdin
    /// if not provided.
    #[clap(short, long)]
    input: Option<PathBuf>,

    /// Path to BAM alignment file used in nanopolish eventalign
    #[clap(short, long)]
    bam: PathBuf,

    #[clap(short, long)]
    /// Path to output file in Apache Arrow format, defaults to stdout if no
    /// argument provided.
    output: Option<PathBuf>,

    #[clap(short, long, default_value_t = 2048)]
    /// Number of eventalign records to hold in memory.
    capacity: usize,
}

impl CollapseCmd {
    pub fn run(self) -> eyre::Result<()> {
        if self.capacity == 0 {
            return Err(eyre::eyre!("Capacity must be greater than 0"));
        }
        let final_input: Box<dyn Read> = {
            if let Some(path) = self.input {
                Box::new(File::open(path)?)
            } else {
                let stdin = io::stdin().lock();
                Box::new(stdin)
            }
        };

        let final_output = utils::stdout_or_file(self.output.as_ref())?;
        let final_output = BufWriter::new(final_output);

        let mut collapse = CollapseOptions::from_writer(final_output, &self.bam)?;
        collapse.capacity(self.capacity).progress(true);
        collapse.run(final_input)?;
        Ok(())
    }
}
