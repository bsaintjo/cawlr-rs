use assert_cmd::Command;
use assert_fs::TempDir;
use escargot::CargoBuild;

#[test]
#[ignore = "Only useful to run inside docker"]
fn test_docker_train_ctrl_pipeline() -> eyre::Result<()> {
    let temp_dir = TempDir::new()?.into_persistent_if(std::env::var("TEST_PERSIST").is_ok());
    log::info!("temp_dir: {}", temp_dir.path().display());

    log::info!("Building release cawlr");
    let run = CargoBuild::new()
        .package("cawlr")
        .release()
        .no_default_features()
        .run()?;
    let cawlr = run.path().as_os_str();
    let train_output = temp_dir.join("train_outputs");
    Command::new(cawlr)
        .arg("pipeline")
        .arg("train-ctrls")
        .arg("-g")
        .arg("extra/sacCer3.fa")
        .arg("--pos-fast5")
        .arg("extra/pos-fast5")
        .arg("--pos-reads")
        .arg("extra/pos-fastq")
        .arg("--neg-fast5")
        .arg("extra/neg-fast5")
        .arg("--neg-reads")
        .arg("extra/neg-fastq")
        .arg("--output-dir")
        .arg(&train_output)
        .arg("--motifs")
        .arg("1:TA,2:GC,1:CG")
        .env("RUST_BACKTRACE", "full")
        .assert()
        .success();

    // let reads_file = train_output.join("pos_reads.fastq");
    let preprocess_output = temp_dir.join("preprocess_output");
    Command::new(cawlr)
        .arg("pipeline")
        .arg("preprocess-sample")
        .arg("-g")
        .arg("extra/sacCer3.fa")
        .arg("--reads")
        .arg("extra/pos-fastq")
        .arg("--fast5")
        .arg("extra/pos-fast5")
        .arg("preprocess_output")
        .arg("-o")
        .arg(&preprocess_output)
        .assert()
        .success();
    Ok(())
}
