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

impl<T> CawlrIO for Vec<T>
where
    T: Serialize + DeserializeOwned + Sized,
{
    fn save<P>(&self, filename: P) -> Result<()>
    where
        P: AsRef<Path>,
    {
        let file = File::create(filename)?;
        let schema = Schema::from_records(&self)?;
        let batches = to_record_batch(&self, &schema)?;
        let mut writer = ArrowWriter::try_new(file, batches.schema(), None)?;
        writer.write(&batches)?;
        writer.close()?;
        Ok(())
    }

    fn load<P>(filename: P) -> Result<Self>
    where
        P: AsRef<Path>,
    {
        let file = File::open(filename)?;
        let file_reader = SerializedFileReader::new(file)?;
        let mut reader = ParquetFileArrowReader::new(Arc::new(file_reader));
        let n = reader.get_metadata().num_row_groups();
        let schema = reader.get_schema()?;
        let schema = Schema::try_from(schema)?;
        let record_reader = reader.get_record_reader(2048)?;
        let mut acc = Vec::with_capacity(n);
        for rb in record_reader.flatten() {
            let mut xs = from_record_batch(&rb, &schema)?;
            acc.append(&mut xs);
        }
        Ok(acc)
    }
}

// todo switch to parquet
impl<K, V> CawlrIO for HashMap<K, V>
where
    K: Serialize + DeserializeOwned + Eq + Hash,
    V: Serialize + DeserializeOwned,
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
