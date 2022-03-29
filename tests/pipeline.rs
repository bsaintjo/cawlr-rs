use std::{error::Error, process::Command};

use assert_cmd::prelude::OutputAssertExt;
use assert_fs::TempDir;
use escargot::CargoBuild;

// lazy_static::lazy_static! {
//     static ref RUN: CargoRun =
// CargoBuild::new().bin("cawlr").release().run().unwrap();     static ref
// CAWLR: &'static Path = RUN.path(); }

#[test]
fn pipeline() -> Result<(), Box<dyn Error>> {
    let temp_dir = TempDir::new()?;

    eprintln!("Building release cawlr version");
    let run = CargoBuild::new().bin("cawlr").release().run().unwrap();
    let cawlr = run.path().as_os_str();

    eprintln!("Preprocessing positive control");
    let pos_output = temp_dir.path().join("pos_control.pickle");
    Command::new(cawlr)
        .arg("preprocess")
        .arg("-b")
        .arg("extra/pos_control.bam")
        .arg("-i")
        .arg("extra/pos_control.eventalign.txt")
        .arg("-g")
        .arg("hundred/sacCer3.fa")
        .arg("-o")
        .arg(&pos_output)
        .assert()
        .success();

    eprintln!("Preprocessing negative control");
    let neg_output = temp_dir.path().join("neg_control.pickle");
    Command::new(cawlr)
        .arg("preprocess")
        .arg("-b")
        .arg("extra/neg_control.bam")
        .arg("-i")
        .arg("extra/neg_control.eventalign.txt")
        .arg("-g")
        .arg("hundred/sacCer3.fa")
        .arg("-o")
        .arg(&neg_output)
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
        .assert()
        .success();

    temp_dir.close()?;
    Ok(())
}
