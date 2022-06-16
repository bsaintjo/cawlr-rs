use anyhow::Result;
use cawlr::CollapseOptions;

#[test_log::test]
fn test_pos_ctrl() -> Result<()> {
    let pos_eventalign = "extra/pos_eventalign.txt";
    let collapse_opts = CollapseOptions::try_new(pos_eventalign, "/dev/null", 4096)?;
    collapse_opts.run()?;
    Ok(())
}
