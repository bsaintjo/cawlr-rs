use std::{
    fs::File,
    path::{Path, PathBuf},
    str::FromStr,
};

use cawlr::{
    collapse::CollapseOptions,
    motif::Motif,
    npsmlr,
    rank::{RankOptions, Ranks},
    train::Model,
    utils::CawlrIO,
};
use eyre::Result;
use fnv::FnvHashMap;

const POS_CTRL: &str = "extra/pos_control.eventalign.txt";
const POS_CTRL_BAM: &str = "extra/pos_control.bam";

const NEG_CTRL: &str = "extra/neg_control.eventalign.txt";
const NEG_CTRL_BAM: &str = "extra/neg_control.bam";

const READ: &str = "extra/single_read.eventalign.txt";
const READ_BAM: &str = "extra/single_read.bam";

// const GENOME: &str = "extra/sacCer3.fa";

fn train(input: &PathBuf, output: &PathBuf, is_neg_ctrl: bool) -> Result<Model> {
    let input = File::open(input)?;
    let train = npsmlr::train::TrainOptions::default();
    let train = train
        .motifs(vec![
            Motif::from_str("2:AT").unwrap(),
            Motif::from_str("1:TA").unwrap(),
        ])
        .single(is_neg_ctrl);
    let model = train.run_model(input)?;
    model.save_as(output)?;
    Ok(model)
}

fn rank(pos_ctrl: &Model, neg_ctrl: &Model, output: &Path) -> Result<Ranks> {
    let mut rank_opts = RankOptions::default();
    let rankings = rank_opts.rank(pos_ctrl, neg_ctrl);
    rankings.save_as(output)?;
    Ok(rankings)
}

fn score(
    pos_model: &Model,
    neg_model: &Model,
    ranks: &FnvHashMap<String, f64>,
    reader: &Path,
    writer: &Path,
) -> Result<()> {
    let reader = File::open(reader)?;
    let writer = File::create(writer)?;
    let motifs = vec![
        Motif::from_str("2:AT").unwrap(),
        Motif::from_str("1:TA").unwrap(),
    ];
    let score_opts = npsmlr::ScoreOptions::new(
        pos_model.clone(),
        neg_model.clone(),
        ranks.clone(),
        10,
        10.0,
        motifs,
    );
    score_opts.run(reader, writer)?;
    Ok(())
}

fn main() -> Result<()> {
    jane_eyre::install()?;
    let args = std::env::args().collect::<Vec<_>>();
    let output_dir = if args.len() <= 1 {
        "test_data"
    } else {
        &args[1]
    };
    let output_dir = PathBuf::from(output_dir);
    println!("{output_dir:?}");
    if !output_dir.exists() {
        std::fs::create_dir_all(&output_dir)?;
    }

    let pos_output = output_dir.join("pos_collapse");
    let neg_output = output_dir.join("neg_collapse");
    let read_output = output_dir.join("read_collapse");

    let neg_ctrl = File::open(NEG_CTRL)?;
    let pos_ctrl = File::open(POS_CTRL)?;
    let read = File::open(READ)?;

    CollapseOptions::try_new(NEG_CTRL_BAM, &neg_output)?.run(neg_ctrl)?;
    CollapseOptions::try_new(POS_CTRL_BAM, &pos_output)?.run(pos_ctrl)?;
    CollapseOptions::try_new(READ_BAM, &read_output)?.run(read)?;

    let pos_model_path = output_dir.join("pos_model");
    let neg_model_path = output_dir.join("neg_model");

    let pos_model = train(&pos_output, &pos_model_path, false)?;
    let neg_model = train(&neg_output, &neg_model_path, true)?;

    let rank_path = output_dir.join("ranks");
    let ranks = rank(&pos_model, &neg_model, &rank_path)?;

    let read_score_path = output_dir.join("read_scores");
    let pos_score_path = output_dir.join("pos_scores");
    let neg_score_path = output_dir.join("neg_scores");
    score(
        &pos_model,
        &neg_model,
        &ranks,
        &read_output,
        &read_score_path,
    )?;
    score(
        &pos_model,
        &neg_model,
        &ranks,
        &pos_output,
        &pos_score_path,
    )?;
    score(
        &pos_model,
        &neg_model,
        &ranks,
        &neg_output,
        &neg_score_path,
    )?;

    Ok(())
}
