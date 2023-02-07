use std::{path::PathBuf, ffi::OsStr};

use clap::{builder::PathBufValueParser, error::ErrorKind};

#[derive(Clone, Debug)]
pub struct ValidPathBuf(pub PathBuf);

impl AsRef<OsStr> for ValidPathBuf {
    fn as_ref(&self) -> &OsStr {
        self.0.as_ref()
    }
}

impl clap::builder::ValueParserFactory for ValidPathBuf {
    type Parser = ValidPathBufParser;
    fn value_parser() -> Self::Parser {
        ValidPathBufParser
    }
}

#[derive(Clone)]
pub struct ValidPathBufParser;

impl clap::builder::TypedValueParser for ValidPathBufParser {
    type Value = ValidPathBuf;

    fn parse_ref(
        &self,
        cmd: &clap::Command,
        arg: Option<&clap::Arg>,
        value: &std::ffi::OsStr,
    ) -> Result<Self::Value, clap::Error> {
        let val = PathBufValueParser::new().parse_ref(cmd, arg, value)?;
        if !val.exists() {
            let err = clap::Error::raw(
                ErrorKind::ValueValidation,
                format!("Path {value:?} does not exist"),
            )
            .with_cmd(cmd);
            Err(err)
        } else {
            Ok(ValidPathBuf(val))
        }
    }
}
