use std::ffi::OsStr;
use std::fs::File;

use anyhow::Result;
use assert_cmd::Command;
use assert_fs::fixture::ChildPath;
use assert_fs::{fixture::PathChild, TempDir};
use polars::io::prelude::ParquetWriter;
use polars::prelude::{DataFrame, NamedFrom, Series};

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
}

fn setup() -> Result<TestFiles> {
    let temp_dir = TempDir::new()?;
    let input_file_path = temp_dir.child("test.parquet");
    let input_file = File::create(&input_file_path)?;
    let df = polars::df!("event_mean" => &[0.1, 0.2, 0.3, 0.4, 0.5],
                                  "kmer" => &["AAAAAA", "AAAAAA", "AAAAAA", "AAAAAA", "AAAAAA"])?;
    ParquetWriter::new(input_file).finish(&df)?;

    let output_file_path = temp_dir.child("results");
    let test_files = TestFiles::new(temp_dir, input_file_path, output_file_path);
    Ok(test_files)
}

#[test]
fn test_train() -> Result<()> {
    let test_files = setup()?;
    let mut cmd = Command::cargo_bin("cawlr")?;
    cmd.arg("train")
        .arg("-i")
        .arg(test_files.input())
        .arg("-o")
        .arg(test_files.output());
    cmd.assert().try_success()?;
    Ok(())
}
