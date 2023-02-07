use std::path::PathBuf;

use clap::Parser;
use libcawlr::{motif::Motif, region::Region};

#[derive(Debug, Parser)]
pub struct AnalyzeCmd {
    /// Region of interested {chromosome}:{start}-{stop}
    #[clap(short, long)]
    pub locus: Region,

    /// Where to output results
    #[clap(short, long)]
    pub output_dir: PathBuf,

    /// Path to bam file to filter on the locus
    #[clap(short, long)]
    pub bam: PathBuf,

    /// Path to full fastq, doesn't need to be filtered
    #[clap(long)]
    pub reads: PathBuf,

    /// Path to genome
    #[clap(short, long)]
    pub genome: PathBuf,

    /// Path to postive control model, from cawlr train
    #[clap(long)]
    pub pos_model: PathBuf,

    /// Path to postive control scores, from cawlr model-scores
    #[clap(long)]
    pub pos_scores: PathBuf,

    /// Path to negative control model, from cawlr train
    #[clap(long)]
    pub neg_model: PathBuf,

    /// Path to negative control scores, from cawlr model-scores
    #[clap(long)]
    pub neg_scores: PathBuf,

    /// Path to ranks file, from cawlr ranks
    #[clap(long)]
    pub ranks: PathBuf,

    /// Number of clusters to use for clustering script
    #[clap(long, default_value_t = 3)]
    pub n_clusters: usize,

    /// Percent of read that should overlap region to be clustered
    #[clap(long)]
    pub pct: f64,

    /// Motifs of modification to filter on, separated by commas, format is
    /// "{position}:{motif}" ie for GpC and CpG motif , motif is "2:GC,1:CG"
    #[clap(short, long, required=true, num_args=1.., value_delimiter=',')]
    pub motifs: Vec<Motif>,

    /// Regions to highlight during clustering
    #[clap(long)]
    pub highlights: Vec<String>,

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
