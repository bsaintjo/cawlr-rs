use anyhow::Result;
use assert_fs::TempDir;
use cawlr::CollapseOptions;

const POS_CTRL: &'static str = "extra/pos_control.eventalign.txt";
const NEG_CTRL: &'static str = "extra/neg_control.eventalign.txt";
const READ: &'static str = "extra/single_read.eventalign.txt";
const DEVNULL: &'static str = "dev/null";

fn test_collapse(input: &str, output: &str) -> Result<()> {
    let collapse_opts = CollapseOptions::try_new(input, output, 4096)?;
    collapse_opts.run()
}

#[test]
#[ignore]
fn test_pos_ctrl() -> Result<()> {
    test_ctrl(POS_CTRL, DEVNULL)
}

#[test]
#[ignore]
fn test_neg_ctrl() -> Result<()> {
    test_ctrl(NEG_CTRL, DEVNULL)
}

#[test]
fn pipeline() -> Result<()> {
    let tmp_dir = TempDir::new()?;

    let pos_output = "pos_ctrl_output";
    let pos_output = tmp_dir.path().join(pos_output);

    let neg_output = "neg_ctrl_output";
    let neg_output = tmp_dir.path().join(neg_output);
    
    let read_output = "read_output";
    let read_output = tmp_dir.path().join(read_output);

    test_collapse(POS_CTRL, pos_output)?;
    test_collapse(NEG_CTRL, neg_output)?;
    test_collapse(READ, read_output)?;

    Ok(())
}