use std::error::Error;

use assert_cmd::Command;
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
    log::info!("temp_dir: {}", temp_dir.path().display());

    log::info!("Building release cawlr");
    let run = CargoBuild::new()
        .package("cawlr")
        .release()
        .no_default_features()
        .run()?;
    let cawlr = run.path().as_os_str();
    log::info!("Validate model-scores with bam files");
    let modbam_model_scores = temp_dir.path().join("modbam-model-scores");

    // Should fail but only because no aligned reads are inside
    Command::new(cawlr)
        .arg("model-scores")
        .arg("-i")
        .arg("extra/modbams/MM-double.bam")
        .arg("-t")
        .arg("C+m")
        .arg("--bins")
        .arg("1000")
        .arg("-o")
        .arg(&modbam_model_scores)
        .assert()
        .failure();

    Command::new(cawlr)
        .arg("model-scores")
        .arg("-i")
        .arg("extra/modbams/megalodon-modbam.bam")
        .arg("-t")
        .arg("A+Y")
        .arg("--bins")
        .arg("1000")
        .arg("-o")
        .arg(&modbam_model_scores)
        .assert()
        .success();

    let sma_output = temp_dir.path().join("sma-output");
    log::info!("sma testing");
    Command::new(cawlr)
        .arg("sma")
        .arg("-i")
        .arg("extra/modbams/megalodon-modbam.bam")
        .arg("-t")
        .arg("A+Y")
        .arg("--pos-ctrl-scores")
        .arg(&modbam_model_scores)
        .arg("--neg-ctrl-scores")
        .arg(&modbam_model_scores)
        .arg("-o")
        .arg(sma_output)
        .assert()
        .success();
    Ok(())
}
