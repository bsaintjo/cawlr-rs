struct Motif {
    motif: String,
    position: usize,
}

impl Motif {
    fn new(motif: String, position: usize) -> Self {
        Self { motif, position }
    }

    fn from_string<T>(string: T, position: usize) -> Option<Self>
    where
        T: Into<String>,
    {
        let motif = string.into();
        if position >= motif.len() {
            None
        } else {
            Some(Motif::new(motif, position))
        }
    }

    fn is_valid(&self) -> bool {
        self.position < self.motif.len()
    }
}
