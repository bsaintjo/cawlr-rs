mod analyze;
mod preprocess;
mod train_ctrls_pipeline;

use clap::Subcommand;

use self::{
    analyze::AnalyzeCmd, preprocess::PreprocessCmd, train_ctrls_pipeline::TrainCtrlPipelineCmd,
};

#[derive(Subcommand, Debug)]
pub enum PipelineCmds {
    AnalyzeRegion(AnalyzeCmd),
    PreprocessSample(PreprocessCmd),
    TrainCtrls(TrainCtrlPipelineCmd),
}

impl PipelineCmds {
    pub fn run(self) -> eyre::Result<()> {
        match self {
            PipelineCmds::AnalyzeRegion(args) => analyze::run(args),
            PipelineCmds::PreprocessSample(cmd) => cmd.run(),
            PipelineCmds::TrainCtrls(cmd) => train_ctrls_pipeline::run(cmd),
        }
    }
}
