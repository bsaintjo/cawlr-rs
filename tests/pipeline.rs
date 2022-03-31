use std::{error::Error, path::Path, process::Command};

use assert_cmd::prelude::OutputAssertExt;
use assert_fs::TempDir;
use escargot::CargoBuild;

#[test]
fn pipeline() -> Result<(), Box<dyn Error>> {
    let temp_dir = TempDir::new()?;

    eprintln!("Building release cawlr");
    let run = CargoBuild::new().bin("cawlr").release().run()?;
    let cawlr = run.path().as_os_str();
    let genome = dunce::realpath("extra/sacCer3.fa")?;

    eprintln!("Preprocessing positive control");
    let pos_output = temp_dir.path().join("pos_control.parquet");
    Command::new(cawlr)
        .arg("preprocess")
        .arg("-b")
        .arg("extra/pos_control.bam")
        .arg("-i")
        .arg("extra/pos_control.eventalign.txt")
        .arg("-g")
        .arg(&genome)
        .arg("-o")
        .arg(&pos_output)
        .env("RUST_BACKTRACE", "full")
        .assert()
        .success();

    eprintln!("Preprocessing negative control");
    let neg_output = temp_dir.path().join("neg_control.parquet");
    Command::new(cawlr)
        .arg("preprocess")
        .arg("-b")
        .arg("extra/neg_control.bam")
        .arg("-i")
        .arg("extra/neg_control.eventalign.txt")
        .arg("-g")
        .arg(&genome)
        .arg("-o")
        .arg(&neg_output)
        .env("RUST_BACKTRACE", "full")
        .assert()
        .success();

    eprintln!("Preprocessing single read.");
    let single_read_output = temp_dir.path().join("single_read.parquet");
    Command::new(cawlr)
        .arg("preprocess")
        .arg("-b")
        .arg("extra/single_read.bam")
        .arg("-i")
        .arg("extra/single_read.eventalign.txt")
        .arg("-g")
        .arg(&genome)
        .arg("-o")
        .arg(&single_read_output)
        .env("RUST_BACKTRACE", "full")
        .assert()
        .success();

    eprintln!("Training on positive control");
    let pos_train = temp_dir.path().join("pos_control.train");
    Command::new(cawlr)
        .arg("train")
        .arg("-i")
        .arg(&pos_output)
        .arg("-o")
        .arg(&pos_train)
        .env("RUST_BACKTRACE", "full")
        .assert()
        .success();

    eprintln!("Training on negative control");
    let neg_train = temp_dir.path().join("neg_control.train");
    Command::new(cawlr)
        .arg("train")
        .arg("-i")
        .arg(&neg_output)
        .arg("-o")
        .arg(&neg_train)
        .env("RUST_BACKTRACE", "full")
        .assert()
        .success();

    eprintln!("Ranking kmers");
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

    eprintln!("Scoring single read");
    let scores = temp_dir.path().join("scores");
    Command::new(cawlr)
        .arg("score")
        .arg("--neg-ctrl")
        .arg(&neg_train)
        .arg("--pos-ctrl")
        .arg(&pos_train)
        .arg("-i")
        .arg(&single_read_output)
        .arg("-r")
        .arg(&ranks)
        .arg("-g")
        .arg(&genome)
        .arg("-o")
        .arg(scores)
        .env("RUST_BACKTRACE", "full")
        .assert()
        .success();

    temp_dir.close()?;
    Ok(())
}
