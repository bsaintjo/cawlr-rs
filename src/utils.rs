use std::{collections::HashMap, fs::File, hash::Hash, path::Path, sync::Arc};

use anyhow::Result;
use parquet::{
    arrow::{ArrowReader, ArrowWriter, ParquetFileArrowReader},
    file::serialized_reader::SerializedFileReader,
};
use serde::{de::DeserializeOwned, Serialize};
use serde_arrow::{from_record_batch, to_record_batch, Schema};
use serde_pickle::from_reader;

pub(crate) trait CawlrIO {
    fn save<P>(&self, filename: P) -> Result<()>
    where
        P: AsRef<Path>,
        Self: Sized;
    fn load<P>(filename: P) -> Result<Self>
    where
        P: AsRef<Path>,
        Self: Sized;
}


// todo switch to parquet
impl<T> CawlrIO for T
where
    T: Serialize + DeserializeOwned
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
