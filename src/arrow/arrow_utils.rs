use std::{
    borrow::Borrow,
    io::{Read, Seek, Write},
    marker::PhantomData,
};

use arrow2::{
    array::Array,
    chunk::Chunk,
    datatypes::{Field, Schema},
    io::ipc::{
        read::{read_file_metadata, FileReader},
        write::{Compression, FileWriter, WriteOptions},
    },
};
use arrow2_convert::{
    deserialize::{arrow_array_deserialize_iterator, ArrowDeserialize, TryIntoCollection},
    field::ArrowField,
    serialize::{ArrowSerialize, TryIntoArrow},
};
use eyre::Result;
use indicatif::{style::TemplateError, ProgressBar, ProgressStyle};
use itertools::Itertools;

use super::{eventalign::Eventalign, scored_read::ScoredRead};

// pub struct ArrowWriter<W: Write>(FileWriter<W>);
pub struct ArrowWriter<W: Write, T> {
    inner: FileWriter<W>,
    _type: PhantomData<T>,
}

impl<W: Write, T> ArrowWriter<W, T> {
    pub fn new(inner: FileWriter<W>) -> Self {
        Self {
            inner,
            _type: PhantomData,
        }
    }
}

/// Helper trait to wrap Writers for saving Arrow files. Only needs to implement
/// type_as_str which is used as a tag in the Arrow StructArray.
// TODO: Eventually replace wrap_writer and saving with ArrowWriter and
// SchemaExt
pub trait SchemaExt: ArrowField {
    fn type_as_str() -> &'static str;
    fn wrap_writer<W: Write>(writer: W) -> Result<ArrowWriter<W, Self>>
    where
        Self: Sized,
    {
        let data_type = Self::data_type();
        let str_type = Self::type_as_str();
        let schema = Schema::from(vec![Field::new(str_type, data_type, false)]);
        let options = WriteOptions {
            compression: Some(Compression::LZ4),
        };
        let fw = FileWriter::try_new(writer, &schema, None, options)?;
        Ok(ArrowWriter::new(fw))
    }
}

impl SchemaExt for Eventalign {
    fn type_as_str() -> &'static str {
        "eventalign"
    }
}

impl SchemaExt for ScoredRead {
    fn type_as_str() -> &'static str {
        "scored"
    }
}

/// Wraps writer for use later with [save].
pub fn wrap_writer<W>(writer: W, schema: &Schema) -> Result<FileWriter<W>>
where
    W: Write,
{
    let options = WriteOptions {
        compression: Some(Compression::LZ4),
    };
    let fw = FileWriter::try_new(writer, schema, None, options)?;
    Ok(fw)
}

/// Writes data to Arrow file
pub fn save<W, T>(writer: &mut FileWriter<W>, x: &[T]) -> Result<()>
where
    T: ArrowField<Type = T> + ArrowSerialize + 'static,
    W: Write,
{
    if !x.is_empty() {
        let arrow_array: Chunk<Box<dyn Array>> = x.try_into_arrow()?;
        writer.write(&arrow_array, None)?;
    }
    Ok(())
}

pub fn save_t<W, T>(writer: &mut ArrowWriter<W, T>, x: &[T]) -> Result<()>
where
    T: ArrowField<Type = T> + ArrowSerialize + 'static,
    W: Write,
{
    if !x.is_empty() {
        let arrow_array: Chunk<Box<dyn Array>> = x.try_into_arrow()?;
        writer.inner.write(&arrow_array, None)?;
    }
    Ok(())
}

pub(crate) fn load<R>(mut reader: R) -> Result<FileReader<R>>
where
    R: Read + Seek,
{
    let metadata = read_file_metadata(&mut reader)?;
    let reader = FileReader::new(reader, metadata, None, None);
    Ok(reader)
}

/// Apply a function to chunks of data loaded from an Arrow Feather File.
///
/// # Example
/// ```rust
/// # use std::fs::File;
/// # use std::error::Error;
/// # use std::io::Cursor;
/// # use libcawlr::arrow::eventalign::Eventalign;
/// # use libcawlr::arrow::arrow_utils::load_apply;
/// # use libcawlr::arrow::arrow_utils::wrap_writer;
/// # use libcawlr::arrow::arrow_utils::save;
/// # fn main() -> Result<(), Box<dyn Error>> {
/// #
/// # let e = Eventalign::default();
/// # let file = Vec::new();
/// # let mut writer = wrap_writer(file, &Eventalign::schema())?;
/// # save(&mut writer, &[e])?;
/// # writer.finish()?;
/// # let file = Cursor::new(writer.into_inner());
/// // Need to specify the type, can be Vec<ScoredRead> as well
/// load_apply(file, |eventalign: Vec<Eventalign>| {
///     // Do stuff with each chunk
/// # Ok(())
/// })?;
/// # Ok(())
/// # }
/// ```
pub fn load_apply<R, F, T>(reader: R, mut func: F) -> Result<()>
where
    R: Read + Seek,
    F: FnMut(Vec<T>) -> eyre::Result<()>,
    T: ArrowField<Type = T> + ArrowDeserialize + 'static,
    for<'a> &'a <T as ArrowDeserialize>::ArrayType: IntoIterator,
{
    let feather = load(reader)?;
    for read in feather {
        if let Ok(chunk) = read {
            for arr in chunk.into_arrays().into_iter() {
                let eventaligns: Vec<T> = arr.try_into_collection()?;
                func(eventaligns)?;
            }
        } else {
            log::warn!("Failed to load arrow chunk")
        }
    }
    Ok(())
}

pub fn load_apply2<R, F, T>(reader: R, mut func: F) -> Result<()>
where
    R: Read + Seek,
    F: FnMut(T) -> Result<()>,
    T: ArrowField<Type = T> + ArrowDeserialize + 'static,
    for<'a> &'a <T as ArrowDeserialize>::ArrayType: IntoIterator,
{
    let feather = load(reader)?;
    for chunk in feather {
        for arr in chunk?.arrays() {
            let iter = arrow_array_deserialize_iterator(arr.as_ref())?;
            for x in iter {
                func(x)?;
            }
        }
    }
    Ok(())
}

/// Trying different ways if iterating over files, can be deleted safely
pub fn load_apply_indy<R, F, T>(reader: R, mut func: F) -> Result<()>
where
    R: Read + Seek,
    F: FnMut(T) -> eyre::Result<()>,
    T: ArrowField<Type = T> + ArrowDeserialize + 'static,
    for<'a> &'a <T as ArrowDeserialize>::ArrayType: IntoIterator,
{
    let feather = load(reader)?;
    for read in feather {
        if let Ok(chunk) = read {
            for arr in chunk.into_arrays().into_iter() {
                let iter = arrow_array_deserialize_iterator(arr.borrow())?;
                for x in iter {
                    func(x)?;
                }
            }
        } else {
            log::warn!("Failed to load arrow chunk")
        }
    }
    Ok(())
}

/// Loops over every chunk in an arrow file, applies the closure to the chunk,
/// and writes the results to the writer.
///
/// TODO make F: ... Result<&[U]> instead to avoid allocation?
pub fn load_read_write<R, W, F, T, U>(
    reader: R,
    mut writer: FileWriter<W>,
    mut func: F,
) -> Result<()>
where
    R: Read + Seek,
    W: Write,
    F: FnMut(Vec<T>) -> eyre::Result<Vec<U>>,
    T: ArrowField<Type = T> + ArrowDeserialize + 'static,
    U: ArrowField<Type = U> + ArrowSerialize + 'static,
    for<'a> &'a <T as ArrowDeserialize>::ArrayType: IntoIterator,
{
    let feather = load(reader)?;
    for read in feather {
        if let Ok(chunk) = read {
            for arr in chunk.into_arrays().into_iter() {
                let eventaligns: Vec<T> = arr.try_into_collection()?;
                let res = func(eventaligns)?;
                save(&mut writer, &res)?;
            }
        } else {
            log::warn!("Failed to load arrow chunk")
        }
    }
    writer.finish()?;
    Ok(())
}

/// Takes a ArrowWriter instead of FileWriter to avoid exposing FileWriter
pub fn load_read_write_arrow<R, W, F, T, U>(reader: R, writer: W, mut func: F) -> Result<()>
where
    R: Read + Seek,
    W: Write,
    F: FnMut(Vec<T>) -> eyre::Result<Vec<U>>,
    T: ArrowField<Type = T> + ArrowDeserialize + 'static,
    U: ArrowField<Type = U> + ArrowSerialize + 'static + SchemaExt,
    for<'a> &'a <T as ArrowDeserialize>::ArrayType: IntoIterator,
{
    let feather = load(reader)?;
    let mut writer = U::wrap_writer(writer)?;
    for read in feather {
        if let Ok(chunk) = read {
            for arr in chunk.into_arrays().into_iter() {
                let eventaligns: Vec<T> = arr.try_into_collection()?;
                let res = func(eventaligns)?;
                save_t(&mut writer, &res)?;
            }
        } else {
            log::error!("Failed to load arrow chunk");
            return Err(eyre::eyre!("Failed to load arrow chunk"));
        }
    }
    writer.inner.finish()?;
    Ok(())
}

pub fn load_read_arrow<R, F, T>(reader: R, mut func: F) -> Result<()>
where
    R: Read + Seek,
    F: FnMut(Vec<T>) -> eyre::Result<()>,
    T: ArrowField<Type = T> + ArrowDeserialize + 'static,
    for<'a> &'a <T as ArrowDeserialize>::ArrayType: IntoIterator,
{
    let feather = load(reader)?;
    for read in feather {
        if let Ok(chunk) = read {
            for arr in chunk.into_arrays().into_iter() {
                let eventaligns: Vec<T> = arr.try_into_collection()?;
                func(eventaligns)?;
            }
        } else {
            log::error!("Failed to load arrow chunk");
            return Err(eyre::eyre!("Failed to load arrow chunk"));
        }
    }
    Ok(())
}

fn block_bar(n_blocks: u64) -> Result<ProgressBar, TemplateError> {
    let style = ProgressStyle::default_bar().template("[{elapsed}] - {percent}% - {bar}")?;
    let pb = ProgressBar::new(n_blocks).with_style(style);

    Ok(pb)
}

pub fn load_read_arrow_measured<R, F, T>(reader: R, mut func: F) -> Result<()>
where
    R: Read + Seek,
    F: FnMut(Vec<T>) -> eyre::Result<()>,
    T: ArrowField<Type = T> + ArrowDeserialize + 'static,
    for<'a> &'a <T as ArrowDeserialize>::ArrayType: IntoIterator,
{
    let feather = load(reader)?;
    let n_blocks = feather.metadata().blocks.len();
    let pb = block_bar(n_blocks as u64)?;
    for read in feather {
        if let Ok(chunk) = read {
            for arr in chunk.into_arrays().into_iter() {
                let eventaligns: Vec<T> = arr.try_into_collection()?;
                func(eventaligns)?;
            }
        } else {
            log::error!("Failed to load arrow chunk");
            return Err(eyre::eyre!("Failed to load arrow chunk"));
        }
        pb.tick();
    }
    pb.finish();
    Ok(())
}
// TODO Refactor multiple maps
#[cfg(test)]
pub(crate) fn load_iter<R>(
    mut reader: R,
) -> impl Iterator<Item = Result<Vec<Eventalign>, arrow2::error::Error>>
where
    R: Read + Seek,
{
    let metadata = read_file_metadata(&mut reader).unwrap();
    let reader = FileReader::new(reader, metadata, None, None);
    reader
        .map(|x| x.map(|c| c.into_arrays().into_iter()))
        .map(|q| {
            q.map(|b| {
                b.flat_map(|r| {
                    let v: Vec<Eventalign> = r.try_into_collection().unwrap();
                    v
                })
                .collect::<Vec<_>>()
            })
        })
}

/// Trying different ways if iterating over files, can be deleted safely
#[allow(dead_code)]
pub(crate) fn load_iter2<R>(
    mut reader: R,
) -> impl Iterator<Item = Result<impl Iterator<Item = Result<Vec<Eventalign>>>, arrow2::error::Error>>
where
    R: Read + Seek,
{
    let metadata = read_file_metadata(&mut reader).unwrap();
    let reader = FileReader::new(reader, metadata, None, None);
    reader.map_ok(|c| {
        c.into_arrays().into_iter().map(|a| {
            let x: Vec<Eventalign> = a.try_into_collection_as_type::<Eventalign>()?;
            Ok(x)
        })
    })
}
