use crate::{arrow::metadata::MetadataExt, region::Region};

pub struct FilterOptions {
    regions: Vec<Region>,
}

impl FilterOptions {
    pub fn new(regions: Vec<Region>) -> Self {
        Self { regions }
    }

    pub fn any_valid<M: MetadataExt + ?Sized>(&self, meta: &M) -> bool {
        self.regions.iter().any(|r| r.valid(meta))
    }
}
