use std::{collections::HashMap, fs::File, hash::{Hash, BuildHasher}, path::Path, sync::Arc};

use anyhow::Result;
use parquet::{
    arrow::{ArrowReader, ArrowWriter, ParquetFileArrowReader},
    file::serialized_reader::SerializedFileReader,
};
use serde::{de::DeserializeOwned, Serialize};
use serde_arrow::{from_record_batch, to_record_batch, trace_schema, Schema};
use serde_pickle::from_reader;
use xxhash_rust::xxh64::Xxh64Builder;

use crate::{
    reads::{FlatLReadLData, Flatten, LData, LRead, Score},
    train::Model,
};

pub(crate) type XHashMap<K, V> = HashMap<K, V, Xxh64Builder>;

pub(crate) fn xxhashmap<K, V>() -> HashMap<K, V, Xxh64Builder> {
    let s = Xxh64Builder::new(2456);
    HashMap::with_hasher(s)
}


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
    S: BuildHasher + Default
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

// todo switch to parquet
impl CawlrIO for Vec<LRead<LData>> {
    fn save<P>(self, filename: P) -> Result<()>
    where
        P: AsRef<Path>,
    {
        let file = File::create(filename)?;
        let results = self.to_flat();
        let schema = trace_schema(&results)?;
        log::debug!("schema {schema:?}");
        let batches = to_record_batch(&results, &schema)?;
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
        let mut acc: Vec<FlatLReadLData> = Vec::with_capacity(n);
        for rb in record_reader.flatten() {
            let mut xs = from_record_batch(&rb, &schema)?;
            acc.append(&mut xs);
        }
        Ok(Flatten::from_flat(acc))
    }
}

impl CawlrIO for Vec<LRead<Score>> {
    fn save<P>(self, filename: P) -> Result<()>
    where
        P: AsRef<Path>,
    {
        let file = File::create(filename)?;
        let results = self.to_flat();
        let schema = Schema::from_records(&results)?;
        let batches = to_record_batch(&results, &schema)?;
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
        let model_db = from_reader(file, Default::default())?;
        Ok(model_db)
    }
}
