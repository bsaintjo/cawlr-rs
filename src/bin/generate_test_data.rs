use std::path::{Path, PathBuf};

use anyhow::Result;
use cawlr::{
    collapse::CollapseOptions,
    rank::{RankOptions, Ranks},
    train::{Model, Train},
    utils::CawlrIO,
};

const POS_CTRL: &'static str = "extra/pos_control.eventalign.txt";
const NEG_CTRL: &'static str = "extra/neg_control.eventalign.txt";
const READ: &'static str = "extra/single_read.eventalign.txt";
const GENOME: &'static str = "extra/sacCer3.fa";

fn collapse(input: &str, output: &PathBuf) -> Result<()> {
    let collapse_opts = CollapseOptions::try_new(input, output, 4096)?;
    collapse_opts.run()
}

fn train(input: &PathBuf, output: &PathBuf) -> Result<Model> {
    let train = Train::try_new(input, GENOME, 2048)?;
    let model = train.run()?;
    model.save(output)?;
    Ok(model)
}

fn rank(pos_ctrl: &Model, neg_ctrl: &Model, output: &Path) -> Result<Ranks> {
    let mut rank_opts = RankOptions::default();
    let rankings = rank_opts.rank(pos_ctrl, neg_ctrl);
    rankings.save(output)?;
    Ok(rankings)
}

fn score() {}

fn main() -> Result<()> {
    let args = std::env::args();
    let output_dir = args.into_iter().next().unwrap_or("test_data".to_string());
    let output_dir = PathBuf::from(output_dir);

    let pos_output = output_dir.join("pos_collapse");
    let neg_output = output_dir.join("neg_collapse");
    let read_output = output_dir.join("read_collapse");

    collapse(POS_CTRL, &pos_output)?;
    collapse(NEG_CTRL, &neg_output)?;
    collapse(READ, &read_output)?;

    let pos_model_path = output_dir.join("pos_model");
    let neg_model_path = output_dir.join("neg_model");

    let pos_model = train(&pos_output, &pos_model_path)?;
    let neg_model = train(&neg_output, &neg_model_path)?;

    let rank_path = output_dir.join("ranks");
    let ranks = rank(&pos_model, &neg_model, &rank_path)?;

    let read_score_path = output_dir.join("read_scores");

    Ok(())
}
