use std::{
    fs::File,
    io::BufReader,
    path::Path,
    process::{Command, Stdio}, ffi::OsStr,
};

use libcawlr::collapse::CollapseOptions;

pub fn eventalign_collapse<P, Q, R, S, T>(
    nanopolish: P,
    reads: Q,
    bam: R,
    genome: S,
    output: T,
    log_file: File,
) -> eyre::Result<()>
where
    P: AsRef<OsStr> + AsRef<Path>,
    Q: AsRef<OsStr> + AsRef<Path>,
    R: AsRef<OsStr> + AsRef<Path>,
    S: AsRef<OsStr> + AsRef<Path>,
    T: AsRef<OsStr> + AsRef<Path>,
{
    let mut cmd = Command::new(nanopolish);
    cmd.arg("eventalign")
        .arg("-r")
        .arg(reads)
        .arg("-b")
        .arg(&bam)
        .arg("-g")
        .arg(genome)
        .arg("-t")
        .arg("4")
        .arg("--scale-events")
        .arg("--print-read-names")
        .arg("--samples");
    log::info!("nanopolish cmd: {cmd:?}");
    let mut cmd = cmd.stdout(Stdio::piped()).stderr(log_file).spawn()?;
    let stdout = cmd
        .stdout
        .take()
        .ok_or_else(|| eyre::eyre!("Could not capture stdout"))?;
    let reader = BufReader::new(stdout);
    let mut collapse = CollapseOptions::try_new(bam, output)?;
    collapse.run(reader)?;
    Ok(())
}
