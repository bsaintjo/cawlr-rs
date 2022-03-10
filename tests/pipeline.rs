use assert_fs::TempDir;
use rand::{prelude::{SmallRng, IteratorRandom}, SeedableRng};


struct NPTestFiles {
    _temp_dir: TempDir,
    read_size: usize,
    seed: u64,
}

fn random_genomic_sequence(bp: usize, seed: u64) -> Vec<u8> {
    let mut rng = SmallRng::seed_from_u64(seed);
    let mut buf = Vec::with_capacity(bp);
    [b'A', b'C', b'G', b'T'].into_iter().choose_multiple_fill(&mut rng, &mut buf);
    buf
}