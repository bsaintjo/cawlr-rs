use serde::{Deserialize, Serialize};

use crate::utils;

pub(crate) type PreprocessRead = LRead<LData>;
pub(crate) type ScoredRead = LRead<Score>;

trait Position {
    fn position(&self) -> usize;

    fn rel_position(&self, x: usize) -> usize {
        x - self.position()
    }
}

#[derive(Debug, Serialize, Deserialize, Clone)]
pub(crate) struct LData {
    pos: u64,
    kmer: String,
    mean: f64,
    time: f64,
}

impl LData {
    pub(crate) fn new(pos: u64, kmer: String, mean: f64, time: f64) -> Self {
        Self {
            pos,
            kmer,
            mean,
            time,
        }
    }

    /// Get the ldata's pos.
    pub(crate) fn pos(&self) -> u64 {
        self.pos
    }

    /// Get the ldata's mean.
    pub(crate) fn mean(&self) -> f64 {
        self.mean
    }

    pub(crate) fn kmer(&self) -> &str {
        &self.kmer
    }
}

#[derive(Debug, Serialize, Deserialize, Clone)]
pub(crate) struct Score {
    pos: u64,
    score: f64,
}

impl Score {
    pub(crate) fn new(pos: u64, score: f64) -> Self {
        Self { pos, score }
    }
}

impl Position for Score {
    fn position(&self) -> usize {
        self.pos as usize
    }
}

#[derive(Debug, Serialize, Deserialize, Clone)]
pub(crate) struct LRead<T> {
    name: Vec<u8>,
    chrom: String,
    start: usize,
    length: usize,
    seq: Vec<u8>,
    #[serde(flatten)]
    data: Vec<T>,
}

impl<T> LRead<T> {
    pub(crate) fn new(
        name: Vec<u8>,
        chrom: String,
        start: usize,
        length: usize,
        seq: Vec<u8>,
        data: Vec<T>,
    ) -> Self {
        LRead {
            name,
            chrom,
            start,
            length,
            seq,
            data,
        }
    }

    fn empty(name: Vec<u8>, chrom: String, start: usize, length: usize, seq: Vec<u8>) -> Self {
        LRead {
            name,
            chrom,
            start,
            length,
            seq,
            data: Vec::new(),
        }
    }

    pub(crate) fn to_empty_lread<U>(&self) -> LRead<U> {
        LRead::empty(
            self.name.clone(),
            self.chrom.clone(),
            self.start,
            self.length,
            self.seq.clone(),
        )
    }
    pub(crate) fn to_lread_with_data<U>(&self, data: Vec<U>) -> LRead<U> {
        let mut lread = LRead::to_empty_lread(self);
        lread.data = data;
        lread
    }

    pub(crate) fn length(&self) -> usize {
        self.length
    }

    pub(crate) fn is_empty(&self) -> bool {
        self.data.is_empty()
    }

    pub(crate) fn get_mut_data(&mut self) -> &mut Vec<T> {
        &mut self.data
    }

    pub(crate) fn iter(&self) -> impl Iterator<Item = &T> {
        self.data.iter()
    }

    /// Get a reference to the lread's chrom.
    pub(crate) fn chrom(&self) -> &str {
        self.chrom.as_ref()
    }

    /// Get the lread's start.
    pub(crate) fn start(&self) -> usize {
        self.start
    }

    /// Get the lread's stop.
    pub(crate) fn stop(&self) -> usize {
        self.start + self.length
    }

    /// Get a reference to the lread's seq.
    pub(crate) fn seq(&self) -> &[u8] {
        self.seq.as_ref()
    }

    /// Get a reference to the lread's data.
    pub(crate) fn data(&self) -> &[T] {
        self.data.as_ref()
    }
}

impl LRead<Score> {
    fn into_array(self) -> LRead<Option<f64>> {
        let mut arr = vec![None; self.length];
        for score in self.iter() {
            arr[score.position()] = Some(score.score)
        }
        LRead {
            name: self.name,
            chrom: self.chrom,
            start: self.start,
            length: self.length,
            seq: self.seq,
            data: arr,
        }
    }

    pub(crate) fn into_means(self) -> impl Iterator<Item = f64> {
        self.data.into_iter().map(|s| s.score)
    }
}

// TODO: Take self and flat_repr instead of ref
pub trait Flatten {
    type Target;

    fn to_flat(self) -> Self::Target;
    fn from_flat(flat_repr: Self::Target) -> Self;
}

#[derive(Debug, Serialize, Deserialize)]
pub struct FlatLReadScore {
    name: String,
    chrom: String,
    start: usize,
    length: usize,
    seq: String,
    pos: u64,
    score: f64,
}

impl FlatLReadScore {
    fn new(
        name: &str,
        chrom: &str,
        start: usize,
        length: usize,
        seq: &str,
        pos: u64,
        score: f64,
    ) -> Self {
        Self {
            name: name.to_owned(),
            chrom: chrom.to_owned(),
            start,
            length,
            seq: seq.to_owned(),
            pos,
            score,
        }
    }
}

#[derive(Debug, Serialize, Deserialize)]
pub struct FlatLReadLData {
    name: String,
    chrom: String,
    start: usize,
    length: usize,
    seq: String,
    pos: u64,
    kmer: String,
    mean: f64,
    time: f64,
}

impl FlatLReadLData {
    fn new(
        name: &[u8],
        chrom: &str,
        start: usize,
        length: usize,
        seq: &[u8],
        pos: u64,
        kmer: &str,
        mean: f64,
        time: f64,
    ) -> Self {
        Self {
            name: String::from_utf8(name.to_owned()).unwrap(),
            chrom: chrom.to_owned(),
            start,
            length,
            seq: String::from_utf8(seq.to_owned()).unwrap(),
            pos,
            kmer: kmer.to_owned(),
            mean,
            time,
        }
    }
}

#[derive(Hash, PartialEq, Eq)]
pub(crate) struct ReprKey {
    name: Vec<u8>,
    chrom: String,
    start: usize,
}

impl ReprKey {
    pub(crate) fn new(name: &[u8], chrom: &str, start: usize) -> Self {
        ReprKey {
            name: name.to_owned(),
            chrom: chrom.to_owned(),
            start,
        }
    }
}


impl Flatten for Vec<LRead<LData>> {
    type Target = Vec<FlatLReadLData>;

    fn to_flat(self) -> Self::Target {
        self.into_iter()
            .flat_map(|lread| {
                lread.data.into_iter().map(move |ldata| FlatLReadLData {
                    name: String::from_utf8(lread.name.clone()).unwrap(),
                    chrom: lread.chrom.clone(),
                    start: lread.start,
                    length: lread.length,
                    seq: String::from_utf8(lread.seq.clone()).unwrap(),
                    pos: ldata.pos,
                    kmer: ldata.kmer,
                    mean: ldata.mean,
                    time: ldata.time,
                })
            })
            .collect()
    }

    fn from_flat(flat_repr: Self::Target) -> Self {
        let mut acc = utils::xxhashmap();
        for flat in flat_repr.iter() {
            let name = flat.name.to_owned();
            let chrom = flat.chrom.to_owned();
            let start = flat.start;
            let repr_key = ReprKey::new(name.as_bytes(), &chrom, start);
            let pos = flat.pos;
            let mean = flat.mean;
            let time = flat.time;
            let kmer = flat.kmer.clone();
            let ldata = LData::new(pos, kmer, mean, time);

            let val = acc.entry(repr_key).or_insert_with(|| {
                let length = flat.length;
                let seq = flat.seq.to_owned();
                LRead::empty(
                    name.as_bytes().to_owned(),
                    chrom,
                    start,
                    length,
                    seq.as_bytes().to_owned(),
                )
            });
            val.data.push(ldata);
        }
        acc.into_values().collect()
    }
}

impl Flatten for Vec<LRead<Score>> {
    type Target = Vec<FlatLReadScore>;
    fn to_flat(self) -> Self::Target {
        self.into_iter()
            .flat_map(|lread| {
                lread.data.into_iter().map(move |score| FlatLReadScore {
                    name: String::from_utf8(lread.name.clone()).unwrap(),
                    chrom: lread.chrom.clone(),
                    start: lread.start,
                    length: lread.length,
                    seq: String::from_utf8(lread.seq.clone()).unwrap(),
                    pos: score.pos,
                    score: score.score,
                })
            })
            .collect()
    }

    fn from_flat(flat_repr: Self::Target) -> Self {
        let mut acc = utils::xxhashmap();
        for flat in flat_repr.iter() {
            let name = flat.name.to_owned();
            let chrom = flat.chrom.to_owned();
            let start = flat.start;
            let repr_key = ReprKey::new(name.as_bytes(), &chrom, start);
            let pos = flat.pos;
            let likely = flat.score;
            let ldata = Score::new(pos, likely);

            let val = acc.entry(repr_key).or_insert_with(|| {
                let length = flat.length;
                let seq = flat.seq.to_owned();
                LRead::empty(
                    name.as_bytes().to_owned(),
                    chrom,
                    start,
                    length,
                    seq.as_bytes().to_owned(),
                )
            });
            val.data.push(ldata);
        }
        acc.into_values().collect()
    }
}
