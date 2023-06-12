mod analyze;
mod preprocess;
mod train_ctrls;
mod external;
mod utils;

use clap::Subcommand;
use log::LevelFilter;

use self::{analyze::AnalyzeCmd, preprocess::PreprocessCmd, train_ctrls::TrainCtrlPipelineCmd};

#[derive(Subcommand, Debug)]
pub enum PipelineCmds {
    /// Train models for a positive and negative control dataset
    TrainCtrls(TrainCtrlPipelineCmd),

    /// Preprocess data, should be done prior to analyze-region
    PreprocessSample(PreprocessCmd),

    /// Analyze a specific locus, producing Genome Browser compatible .bed files
    /// for visualizing nucleosomes on single molecules, and clustering of
    /// nucleosome density
    AnalyzeRegion(AnalyzeCmd),
}

impl PipelineCmds {
    pub fn run(self, log_level_filter: LevelFilter) -> eyre::Result<()> {
        match self {
            PipelineCmds::AnalyzeRegion(args) => analyze::run(args, log_level_filter),
            PipelineCmds::PreprocessSample(cmd) => cmd.run(),
            PipelineCmds::TrainCtrls(cmd) => train_ctrls::run(cmd),
        }
    }
}
