use std::{error::Error, process::Command};

use assert_cmd::prelude::OutputAssertExt;
use assert_fs::TempDir;
use escargot::CargoBuild;
use log::LevelFilter;

#[test]
fn integration_npsmlr() -> Result<(), Box<dyn Error>> {
    env_logger::builder()
        .filter_level(LevelFilter::Info)
        .filter_module("escargot::format", LevelFilter::Off)
        .try_init()?;
    let temp_dir = TempDir::new()?.into_persistent_if(std::env::var("TEST_PERSIST").is_ok());

    log::info!("Building release cawlr");
    let run = CargoBuild::new()
        .bin("cawlr")
        .release()
        .no_default_features()
        .run()?;
    let cawlr = run.path().as_os_str();

    log::info!("Preprocessing positive control");
    let pos_output = temp_dir.path().join("pos_control.output");
    Command::new(cawlr)
        .arg("collapse")
        .arg("-i")
        .arg("extra/pos_control.eventalign.txt")
        .arg("-b")
        .arg("extra/pos_control.bam")
        .arg("-o")
        .arg(&pos_output)
        .env("RUST_BACKTRACE", "full")
        .assert()
        .success();

    log::info!("Preprocessing negative control");
    let neg_output = temp_dir.path().join("neg_control.output");
    Command::new(cawlr)
        .arg("collapse")
        .arg("-i")
        .arg("extra/neg_control.eventalign.txt")
        .arg("-b")
        .arg("extra/neg_control.bam")
        .arg("-o")
        .arg(&neg_output)
        .env("RUST_BACKTRACE", "full")
        .assert()
        .success();

    log::info!("Preprocessing single read.");
    let single_read_output = temp_dir.path().join("single_read.output");
    Command::new(cawlr)
        .arg("collapse")
        .arg("-i")
        .arg("extra/single_read.eventalign.txt")
        .arg("-b")
        .arg("extra/single_read.bam")
        .arg("-o")
        .arg(&single_read_output)
        .env("RUST_BACKTRACE", "full")
        .assert()
        .success();

    Command::new(cawlr)
        .arg("qc")
        .arg("eventalign")
        .arg("-i")
        .arg(&single_read_output)
        .assert()
        .success();

    log::info!("Training on positive control");
    let pos_train = temp_dir.path().join("pos_control.train");
    Command::new(cawlr)
        .arg("npsmlr")
        .arg("train")
        .arg("-i")
        .arg(&pos_output)
        .arg("-o")
        .arg(&pos_train)
        .env("RUST_BACKTRACE", "full")
        .assert()
        .success();

    log::info!("Training on negative control");
    let neg_train = temp_dir.path().join("neg_control.train");
    Command::new(cawlr)
        .arg("npsmlr")
        .arg("train")
        .arg("-i")
        .arg(&neg_output)
        .arg("-o")
        .arg(&neg_train)
        .arg("--single")
        .env("RUST_BACKTRACE", "full")
        .assert()
        .success();

    log::info!("Ranking kmers");
    let ranks = temp_dir.path().join("ranks");
    Command::new(cawlr)
        .arg("rank")
        .arg("--neg-ctrl")
        .arg(&neg_train)
        .arg("--pos-ctrl")
        .arg(&pos_train)
        .arg("-o")
        .arg(&ranks)
        .env("RUST_BACKTRACE", "full")
        .assert()
        .success();

    log::info!("Scoring single read");
    let scores = temp_dir.path().join("single_scores");
    Command::new(cawlr)
        .arg("npsmlr")
        .arg("score")
        .arg("-m")
        .arg("2:GC")
        .arg("--neg-ctrl")
        .arg(&neg_train)
        .arg("--pos-ctrl")
        .arg(&pos_train)
        .arg("-i")
        .arg(&single_read_output)
        .arg("-r")
        .arg(&ranks)
        .arg("-o")
        .arg(&scores)
        .env("RUST_BACKTRACE", "1")
        .assert()
        .success();

    log::info!("Validating score output format");
    Command::new(cawlr)
        .arg("qc")
        .arg("score")
        .arg("-i")
        .arg(&scores)
        .assert()
        .success();

    log::info!("Scoring positive controls");
    let pos_scores = temp_dir.path().join("pos_scores");
    Command::new(cawlr)
        .arg("npsmlr")
        .arg("score")
        .arg("-m")
        .arg("2:GC")
        .arg("--neg-ctrl")
        .arg(&neg_train)
        .arg("--pos-ctrl")
        .arg(&pos_train)
        .arg("-i")
        .arg(&pos_output)
        .arg("-r")
        .arg(&ranks)
        .arg("-o")
        .arg(&pos_scores)
        .env("RUST_BACKTRACE", "1")
        .assert()
        .success();

    log::info!("Scoring negative controls");
    let neg_scores = temp_dir.path().join("neg_scores");
    Command::new(cawlr)
        .arg("npsmlr")
        .arg("score")
        .arg("-m")
        .arg("2:GC")
        .arg("--neg-ctrl")
        .arg(&neg_train)
        .arg("--pos-ctrl")
        .arg(&pos_train)
        .arg("-i")
        .arg(&neg_output)
        .arg("-r")
        .arg(&ranks)
        .arg("-o")
        .arg(&neg_scores)
        .env("RUST_BACKTRACE", "1")
        .assert()
        .success();

    log::info!("Compute pos ctrl kernel density estimate");
    let pos_bkde_model = temp_dir.path().join("pos_bkde_model");
    Command::new(cawlr)
        .arg("model-scores")
        .arg("-i")
        .arg(&pos_scores)
        .arg("--bins")
        .arg("1000")
        .arg("-o")
        .arg(&pos_bkde_model)
        .env("RUST_BACKTRACE", "1")
        .assert()
        .success();

    log::info!("Compute neg ctrl kernel density estimate");
    let neg_bkde_model = temp_dir.path().join("neg_bkde_model");
    Command::new(cawlr)
        .arg("model-scores")
        .arg("-i")
        .arg(&neg_scores)
        .arg("--bins")
        .arg("1000")
        .arg("-o")
        .arg(&neg_bkde_model)
        .env("RUST_BACKTRACE", "1")
        .assert()
        .success();

    log::info!("Single molecule analysis");
    let sma_bed = temp_dir.path().join("sma_bed");
    Command::new(cawlr)
        .arg("sma")
        .arg("--neg-ctrl-scores")
        .arg(&neg_bkde_model)
        .arg("--pos-ctrl-scores")
        .arg(&pos_bkde_model)
        .arg("-i")
        .arg(&scores)
        .arg("-o")
        .arg(&sma_bed)
        .env("RUST_BACKTRACE", "1")
        .assert()
        .success();

    log::info!("Cleaning up");
    temp_dir.close()?;
    Ok(())
}
