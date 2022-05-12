use std::{
    collections::HashMap,
    fs::File,
    hash::{BuildHasher, Hash},
    path::Path,
};

use anyhow::Result;
use serde::{de::DeserializeOwned, Serialize};
use serde_pickle::from_reader;

use crate::train::Model;

pub(crate) trait CawlrIO {
    fn save<P>(self, filename: P) -> Result<()>
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
    fn save<P>(self, filename: P) -> Result<()>
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
    fn save<P>(self, filename: P) -> Result<()>
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
