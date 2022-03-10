use serde::{Deserialize, Serialize};

pub(crate) type PreprocessRead = LRead<LData>;
pub(crate) type ScoredRead = LRead<Score>;

trait Position {
    fn position(&self) -> usize;

    fn rel_position(&self, x: usize) -> usize {
        x - self.position()
    }
}

#[derive(Debug, Serialize, Deserialize)]
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

    /// Get a reference to the ldata's kmer.
    pub(crate) fn into_kmer(self) -> String {
        self.kmer
    }
}

#[derive(Debug, Serialize, Deserialize)]
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

#[derive(Debug, Serialize, Deserialize)]
pub(crate) struct LRead<T> {
    name: Vec<u8>,
    chrom: String,
    start: usize,
    length: usize,
    #[serde(flatten)]
    data: Vec<T>,
}

impl<T> LRead<T> {
    pub(crate) fn new(
        name: Vec<u8>,
        chrom: String,
        start: usize,
        length: usize,
        data: Vec<T>,
    ) -> Self {
        LRead {
            name,
            chrom,
            start,
            length,
            data,
        }
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

    pub(crate) fn map_data<F, U>(self, f: F) -> LRead<U>
    where
        F: FnMut(T) -> U,
    {
        let new_data = self.data.into_iter().map(f).collect();
        LRead {
            name: self.name,
            chrom: self.chrom,
            start: self.start,
            length: self.length,
            data: new_data,
        }
    }

    pub(crate) fn iter(&self) -> impl Iterator<Item = &T> {
        self.data.iter()
    }

    pub(crate) fn into_iter(self) -> impl Iterator<Item = T> {
        self.data.into_iter()
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
}

impl LRead<Score> {
    fn into_array(self) -> LRead<Option<f64>> {
        let mut arr = vec![None; self.length()];
        for score in self.iter() {
            arr[score.position()] = Some(score.score)
        }
        LRead {
            name: self.name,
            chrom: self.chrom,
            start: self.start,
            length: self.length,
            data: arr,
        }
    }

    pub(crate) fn into_means(self) -> impl Iterator<Item = f64> {
        self.data.into_iter().map(|s| s.score)
    }
}
