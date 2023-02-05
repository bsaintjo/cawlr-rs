pub mod collapse;
pub mod score;
pub mod train;

#[cfg(test)]
mod test {
    use std::path::PathBuf;

    use assert_fs::TempDir;
    use libcawlr::motif::all_bases;

    use super::*;

    #[ignore = "Useful only if integration test is failing"]
    #[test]
    fn test_training_ctrls() -> eyre::Result<()> {
        let temp_dir = TempDir::new()?.into_persistent_if(std::env::var("TEST_PERSIST").is_ok());
        let collapse_output = temp_dir.join("pos_collapse");
        let collapse_cmd = collapse::CollapseCmd {
            input: Some(PathBuf::from("../extra/pos_control.eventalign.txt")),
            bam: PathBuf::from("../extra/pos_control.bam"),
            output: Some(collapse_output.clone()),
            capacity: 2048,
        };
        collapse_cmd.run()?;

        let train_output = temp_dir.join("train_output");
        let train_db_output = temp_dir.join("train_db");
        let train_cmd = train::TrainCmd {
            input: collapse_output,
            output: train_output,
            motif: all_bases(),
            samples: 50000,
            single: false,
            dbscan: true,
            db_path: Some(train_db_output),
        };
        train_cmd.run()?;
        Ok(())
    }
}
