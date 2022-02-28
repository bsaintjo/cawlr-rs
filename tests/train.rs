use std::{ffi::OsStr, fs::File};

use anyhow::Result;
use assert_cmd::Command;
use assert_fs::{
    fixture::{ChildPath, PathChild},
    NamedTempFile, TempDir,
};
use polars::{
    io::prelude::ParquetWriter,
    prelude::{DataFrame, NamedFrom, Series},
};

struct TestFiles {
    // Keep around TempDir to extend lifetime
    _temp_dir: TempDir,
    input: ChildPath,
    output: ChildPath,
}

impl TestFiles {
    fn new(_temp_dir: TempDir, input: ChildPath, output: ChildPath) -> Self {
        TestFiles {
            _temp_dir,
            input,
            output,
        }
    }

    fn input(&self) -> &OsStr {
        self.input.as_os_str()
    }

    fn output(&self) -> &OsStr {
        self.output.as_os_str()
    }
    fn setup(input_filename: &str, output_filename: &str) -> Result<Self> {
        let temp_dir = TempDir::new()?;
        let input_file_path = temp_dir.child(input_filename);
        let input_file = File::create(&input_file_path)?;
        let df = polars::df!("event_mean" => &[0.1, 0.2, 0.3, 0.4, 0.5],
                                    "kmer" => &["AAAAAA", "AAAAAA", "AAAAAA", "AAAAAA", "AAAAAA"])?;
        ParquetWriter::new(input_file).finish(&df)?;

        let output_file_path = temp_dir.child(output_filename);
        let test_files = TestFiles::new(temp_dir, input_file_path, output_file_path);
        Ok(test_files)
    }

    fn sim_processed_data(n_kmers: usize, samples_per_kmer: usize) -> DataFrame {
        unimplemented!()
    }
}

#[test]
fn test_train() -> Result<()> {
    let test_files = TestFiles::setup("train", "model")?;
    let mut cmd = Command::cargo_bin("cawlr")?;
    cmd.arg("train")
        .arg("-i")
        .arg(test_files.input())
        .arg("-o")
        .arg(test_files.output());
    cmd.assert().try_success()?;
    Ok(())
}

#[test]
fn test_train_to_rank() -> Result<()> {
    let pos_ctrl_files = TestFiles::setup("pos_ctrl.train", "pos_ctrl.models")?;
    let neg_ctrl_files = TestFiles::setup("neg_ctrl.train", "neg_ctrl.models")?;

    let mut pos_cmd = Command::cargo_bin("cawlr")?;
    pos_cmd
        .arg("train")
        .arg("-i")
        .arg(pos_ctrl_files.input())
        .arg("-o")
        .arg(pos_ctrl_files.output());
    pos_cmd.assert().try_success()?;

    let mut neg_cmd = Command::cargo_bin("cawlr")?;
    neg_cmd
        .arg("train")
        .arg("-i")
        .arg(neg_ctrl_files.input())
        .arg("-o")
        .arg(neg_ctrl_files.output());
    neg_cmd.assert().try_success()?;

    let result = NamedTempFile::new("result")?;
    let mut rank_cmd = Command::cargo_bin("cawlr")?;
    rank_cmd
        .arg("rank")
        .arg("--pos-ctrl")
        .arg(pos_ctrl_files.output())
        .arg("--neg-ctrl")
        .arg(neg_ctrl_files.output())
        .arg("-o")
        .arg(result.path());
    rank_cmd.assert().try_success()?;

    Ok(())
}
