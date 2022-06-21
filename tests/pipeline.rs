use anyhow::Result;
use cawlr::CollapseOptions;

fn test_ctrl(file_path: &str) -> Result<()> {
    let collapse_opts = CollapseOptions::try_new(file_path, "/dev/null", 4096)?;
    collapse_opts.run()
}

#[test_log::test]
fn pipeline() -> Result<()> {
    let pos_ctrl = "extra/pos_control.eventalign.txt";
    let neg_ctrl = "extra/neg_control.eventalign.txt";

    test_ctrl(pos_ctrl)?;
    test_ctrl(neg_ctrl)?;

    Ok(())
}