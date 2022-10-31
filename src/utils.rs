use std::{
    collections::HashMap,
    fs::File,
    hash::{BuildHasher, Hash},
    io::{stdout, Read, Seek, Write},
    path::{Path, PathBuf}, time::Duration,
};

use bio::io::fasta::IndexedReader;
use eyre::{Context, Result};
use fnv::FnvHashMap;
use indicatif::{ProgressBar, ProgressStyle};
use serde::{de::DeserializeOwned, Serialize};
use serde_pickle::from_reader;
use which::which;

use crate::train::Model;

/// Allows for writing to File or Stdout depending on if a filename is given.
///
/// TODO: Maybe return with the BufWriter wrapping the trait object, like
/// BufWriter<Box<dyn Write>> instead of the how we have now.
pub fn stdout_or_file<P>(filename: Option<&P>) -> Result<Box<dyn Write>>
where
    P: AsRef<Path>,
{
    if let Some(fp) = filename {
        let handle = File::create(fp)?;
        Ok(Box::new(handle))
    } else {
        let handle = stdout().lock();
        Ok(Box::new(handle))
    }
}

pub trait CawlrIO {
    fn save<P>(&self, filename: P) -> Result<()>
    where
        P: AsRef<Path>,
        Self: Sized;
    fn load<P>(filename: P) -> Result<Self>
    where
        P: AsRef<Path>,
        Self: Sized;
}
impl<K, V, S> CawlrIO for HashMap<K, V, S>
where
    K: Eq + Hash + Serialize + DeserializeOwned,
    V: Serialize + DeserializeOwned,
    S: BuildHasher + Default,
{
    fn save<P>(&self, filename: P) -> Result<()>
    where
        P: AsRef<Path>,
    {
        let mut file = File::create(filename)?;
        serde_pickle::to_writer(&mut file, &self, Default::default())?;
        Ok(())
    }

    fn load<P>(filename: P) -> Result<Self>
    where
        P: AsRef<Path>,
    {
        let file = File::open(filename)?;
        let model_db = from_reader(file, Default::default())?;
        Ok(model_db)
    }
}

impl CawlrIO for Model {
    fn save<P>(&self, filename: P) -> Result<()>
    where
        P: AsRef<Path>,
    {
        let mut file = File::create(filename)?;
        serde_pickle::to_writer(&mut file, &self, Default::default())?;
        Ok(())
    }

    fn load<P>(filename: P) -> Result<Self>
    where
        P: AsRef<Path>,
    {
        let file = File::open(filename)?;
        let model_db = from_reader(file, Default::default())?;
        Ok(model_db)
    }
}

/// Get the size of each chromosome in the genome fasta file. Later used if
/// fetching sequences and want to avoid trying to pull sequence past the end of
/// the chromosome.
pub(crate) fn chrom_lens<R>(genome: &IndexedReader<R>) -> FnvHashMap<String, u64>
where
    R: Read + Seek,
{
    let mut chrom_lens = FnvHashMap::default();
    genome.index.sequences().into_iter().for_each(|sequence| {
        chrom_lens.insert(sequence.name, sequence.len);
    });
    chrom_lens
}

pub fn find_binary(name: &'static str, binary_filepath: &Option<PathBuf>) -> eyre::Result<PathBuf> {
    if let Some(p) = binary_filepath {
        Ok(p.to_path_buf())
    } else {
        which(name).wrap_err("Error finding {name}")
    }
}

pub fn wrap_cmd<F>(msg: &'static str, mut f: F) -> eyre::Result<()>
where
    F: FnMut() -> eyre::Result<()>,
{
    let p = ProgressBar::new_spinner()
        .with_style(
            ProgressStyle::with_template("{spinner:.green} [{elapsed_precise}] {msg}").unwrap(),
        )
        .with_message(msg);
    p.enable_steady_tick(Duration::from_millis(100));
    f()?;
    // p.finish_with_message(format!("✅ \"{}\" complete", msg));
    // Ok(())

    if let Ok(()) = f() {
        p.finish_with_message(format!("✅ \"{}\" complete", msg));
        Ok(())
    } else {
        p.finish_with_message(format!("❌ \"{}\" failed", msg));
        Err(eyre::eyre!("Previous command failed, check log.txt"))
    }
}