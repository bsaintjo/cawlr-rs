use std::{ffi::OsString, path::PathBuf};

macro_rules! file_newtype {
    ($nt: ident) => {
        #[derive(Clone)]
        pub struct $nt(PathBuf);
        impl From<OsString> for $nt {
            fn from(x: OsString) -> Self {
                $nt(PathBuf::from(x))
            }
        }
    };
}

file_newtype!(GenomeFile);
file_newtype!(BamFile);
file_newtype!(Fast5Directory);