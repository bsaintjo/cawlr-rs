use serde::{Serialize, Deserialize};

type PreprocessRead = Read<Data>;
type ScoredRead = Read<Score>;

#[derive(Serialize, Deserialize)]
struct Read<T> {
    name: String,
    start: usize,
    stop: usize,
    #[serde(flatten)]
    data: Vec<T>
}

#[derive(Serialize, Deserialize)]
struct Data {
    pos: usize,
    kmer: String,
    mean: f64,
    time: f64
}

#[derive(Serialize, Deserialize)]
struct Score {
    pos: usize,
    score: f64
}