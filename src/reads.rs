use serde::{Deserialize, Serialize};

type PreprocessRead = Read<Data>;
type ScoredRead = Read<Score>;

#[derive(Serialize, Deserialize)]
struct Data {
    pos: usize,
    kmer: String,
    mean: f64,
    time: f64,
}

#[derive(Serialize, Deserialize)]
struct Score {
    pos: usize,
    score: f64,
}

#[derive(Serialize, Deserialize)]
struct Read<T> {
    name: String,
    chrom: String,
    start: usize,
    stop: usize,
    #[serde(flatten)]
    data: Vec<T>,
}

impl<T> Read<T> {
    fn new(name: String, chrom: String, start: usize, stop: usize, data: Vec<T>) -> Self {
        Read {
            name,
            chrom,
            start,
            stop,
            data,
        }
    }

    fn map<U>(self, f: fn(T) -> U) -> Read<U> {
        let new_data = self.data.into_iter().map(f).collect();
        Read {
            data: new_data,
            ..self
        }
    }

    fn iter(&self) -> impl Iterator<Item = &T> {
        self.data.iter()
    }
}
