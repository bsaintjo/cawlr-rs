use std::{fs::File, io, path::Path};

use super::{
    arrow_utils::{is_arrow_file, load_apply_indy},
    mod_bam::{BamRecords, ModBamIter},
    scored_read::ScoredRead,
};

pub enum ModFile {
    Arrow(File),
    ModBam { file: File, mod_tag: Vec<u8> },
}

impl ModFile {
    pub fn open_arrow<P: AsRef<Path>>(path: P) -> io::Result<Self> {
        let file = File::open(path)?;
        Ok(Self::Arrow(file))
    }

    pub fn open_mod_bam<P, B>(path: P, mod_tag: B) -> io::Result<Self>
    where
        P: AsRef<Path>,
        B: Into<Vec<u8>>,
    {
        let file = File::open(path)?;
        Ok(Self::ModBam {
            file,
            mod_tag: mod_tag.into(),
        })
    }

    pub fn open_path<P, B>(path: P, tag: Option<B>) -> eyre::Result<Self>
    where
        P: AsRef<Path>,
        B: Into<Vec<u8>>,
    {
        let tag: Option<Vec<u8>> = tag.map(|t| t.into());
        let mod_file = match (path.as_ref().extension(), tag) {
            (Some(ext), _) if ext == "arrow" => ModFile::open_arrow(&path)?,
            (Some(ext), tag) if ext == "bam" => {
                let Some(tag) = tag else { return Err(eyre::eyre!("Detected bam file but no tag given, please from tag with -t/--tag parameter. See -h/--help for more info"))};
                ModFile::open_mod_bam(&path, tag)?
            }
            (None, tag) if is_bam_file(&path) => {
                let Some(tag) = tag else { return Err(eyre::eyre!("Detected bam file but no tag given, please from tag with -t/--tag parameter. See -h/--help for more info"))};
                ModFile::open_mod_bam(&path, tag)?
            }
            (None, None) if is_arrow_file(&path) => ModFile::open_arrow(&path)?,
            (_, _) => return Err(eyre::eyre!("Failed to detect input as .bam or .arrow file")),
        };
        Ok(mod_file)
    }
}

/// Try to read modification bam data from path, if it fails, try to read as an
/// Arrow file. If both those fail, then error out.
pub fn read_mod_bam_or_arrow<F>(mod_file: ModFile, mut f: F) -> eyre::Result<()>
where
    F: FnMut(ScoredRead) -> eyre::Result<()>,
{
    match mod_file {
        ModFile::Arrow(file) => {
            log::info!("Detected arrow file");
            load_apply_indy(file, f)
        }
        ModFile::ModBam { file, mod_tag } => {
            log::info!("Detected modification bam file");
            let records = BamRecords::from_file(file)?;
            let mut iter = ModBamIter::new(records, mod_tag);
            while let Some(res) = iter.next() {
                if res.is_err() {
                    log::warn!("Failed to convert to modbam to ScoredRead: {res:?}");
                    continue;
                }
                let mba = res.unwrap();

                // TODO Avoid clone by pass it into the error
                let rec = mba.rec.clone();
                match mba.try_into() {
                    Ok(scored_read) => f(scored_read)?,
                    Err(e) => {
                        log::warn!(
                            "{} failed with error {e}",
                            String::from_utf8_lossy(rec.name())
                        );
                    }
                }
            }
            Ok(())
        }
    }
    // match BamRecords::from_file(&path) {
    //     Ok(brs) => {
    //         let mut iter = ModBamIter::new(brs, base_mod);
    //         while let Some(res) = iter.next() {
    //             let mba = res?;
    //             let sr: ScoredRead = mba.try_into()?;
    //             f(sr)?;
    //         }
    //         Ok(())
    //     }
    //     Err(bre) => {
    //         let file = File::open(&path)?;
    //         let err = load_apply_indy(file, f);
    //         err
    //     }
    // }
}

fn is_bam_file<P: AsRef<Path>>(path: P) -> bool {
    BamRecords::from_path(path).is_ok()
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    #[ignore = "Use gen-test-data to create example scored arrow file"]
    fn test_arrow_input() {
        let arrow_file = ModFile::open_arrow("scored-example.arrow").unwrap();
        let res = read_mod_bam_or_arrow(arrow_file, |_| Ok(()));
        assert!(res.is_ok())
    }

    #[test]
    fn test_bam_input() {
        let modbam_file = "extra/modbams/megalodon-modbam.bam";
        let mod_tag = "A+Y";
        let modbam = ModFile::open_mod_bam(modbam_file, mod_tag).unwrap();
        let mut count = 0;
        let res = read_mod_bam_or_arrow(modbam, |_| {
            count += 1;
            Ok(())
        });
        assert!(res.is_ok());
        assert_eq!(count, 1);
    }

    #[test]
    fn test_nonexistent_file() {
        let bad_file = "random_file";
        let arrow_file = ModFile::open_arrow(bad_file);
        assert!(arrow_file.is_err());
        let modbam_file = ModFile::open_mod_bam(bad_file, "C+m");
        assert!(modbam_file.is_err());
    }

    #[test]
    fn test_bam_input_arrow_expected() {
        let modbam_file = "extra/modbams/megalodon-modbam.bam";
        let modbam = ModFile::open_arrow(modbam_file).unwrap();
        let res = read_mod_bam_or_arrow(modbam, |_| Ok(()));
        assert!(res.is_err())
    }

    #[test]
    fn test_not_bam() {
        let path = "extra/neg_control.eventalign.txt";
        let reader = BamRecords::from_path(path);
        assert!(reader.is_err())
    }
}
