use std::{
    fs::File,
    path::{Path, PathBuf},
};

use cawlr::{
    collapse::CollapseOptions,
    rank::{RankOptions, Ranks},
    train::{Model, Train, TrainStrategy},
    utils::CawlrIO,
};
use eyre::Result;

const POS_CTRL: &str = "extra/pos_control.eventalign.txt";
const POS_CTRL_BAM: &str = "extra/pos_control.bam";

const NEG_CTRL: &str = "extra/neg_control.eventalign.txt";
const NEG_CTRL_BAM: &str = "extra/neg_control.bam";

const READ: &str = "extra/single_read.eventalign.txt";
const READ_BAM: &str = "extra/single_read.bam";

const GENOME: &str = "extra/sacCer3.fa";

fn train(input: &PathBuf, output: &PathBuf) -> Result<Model> {
    let train = Train::try_new(input, GENOME, 2048, TrainStrategy::AvgSample)?;
    let model = train.run()?;
    model.save_as(output)?;
    Ok(model)
}

fn rank(pos_ctrl: &Model, neg_ctrl: &Model, output: &Path) -> Result<Ranks> {
    let mut rank_opts = RankOptions::default();
    let rankings = rank_opts.rank(pos_ctrl, neg_ctrl);
    rankings.save_as(output)?;
    Ok(rankings)
}

fn main() -> Result<()> {
    let args = std::env::args();
    let output_dir = args
        .into_iter()
        .next()
        .unwrap_or_else(|| "test_data".to_string());
    let output_dir = PathBuf::from(output_dir);

    let pos_output = output_dir.join("pos_collapse");
    let neg_output = output_dir.join("neg_collapse");
    let read_output = output_dir.join("read_collapse");

    let neg_ctrl = File::open(NEG_CTRL)?;
    let pos_ctrl = File::open(POS_CTRL)?;
    let read = File::open(READ)?;

    CollapseOptions::try_new(NEG_CTRL_BAM, &neg_output)?.run(neg_ctrl)?;
    CollapseOptions::try_new(POS_CTRL_BAM, &pos_output)?.run(pos_ctrl)?;
    CollapseOptions::try_new(READ_BAM, read_output)?.run(read)?;

    let pos_model_path = output_dir.join("pos_model");
    let neg_model_path = output_dir.join("neg_model");

    let pos_model = train(&pos_output, &pos_model_path)?;
    let neg_model = train(&neg_output, &neg_model_path)?;

    let rank_path = output_dir.join("ranks");
    let _ranks = rank(&pos_model, &neg_model, &rank_path)?;

    let _read_score_path = output_dir.join("read_scores");

    Ok(())
}
