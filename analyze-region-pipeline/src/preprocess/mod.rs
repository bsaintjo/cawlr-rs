use std::path::PathBuf;

use clap::Parser;


#[derive(Parser)]
pub struct PreprocessCmd {

    #[clap(short, long)]
    pub genome: PathBuf,

    #[clap(long)]
    pub reads: PathBuf,

    #[clap(long)]
    pub fast5: PathBuf,

    #[clap(short, long)]
    pub output_dir: PathBuf,

    /// Path to nanopolish binary, if not specified will look in $PATH
    #[clap(long)]
    pub nanopolish_path: Option<PathBuf>,

    /// Path to samtools binary, if not specified will look in $PATH
    #[clap(long)]
    pub samtools_path: Option<PathBuf>,

    #[clap(long, default_value_t = false)]
    pub overwrite: bool,

    #[clap(short = 'j', long, default_value_t = 4)]
    pub n_threads: usize,
}
