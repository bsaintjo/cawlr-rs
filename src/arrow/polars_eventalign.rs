use std::path::Path;

#[cfg(test)]
mod test {
    use std::{default, fs::File};

    use polars::lazy::frame::{LazyFileListReader, LazyFrame, OptState};
    use polars::prelude::ListNameSpaceExtension;
    use serde::{Deserialize, Serialize};
    use tempfile::TempDir;

    use super::*;
    use polars::prelude as pl;

    #[derive(Debug, Deserialize, Serialize)]
    struct PolarsEventalign {
        read_name: String,
        contig: String,
        position: Vec<i64>,
        model_kmer: Vec<String>,
        event_length: Vec<f64>,
        samples: Vec<Vec<f64>>,
        n_samples: Vec<u32>,
        sample_mean: Vec<f64>,
        n_positions: u32,
        start: i64,
        end: i64,
        length: i64,
    }

    #[test]
    fn test_polars() {
        std::env::set_var("POLARS_FMT_TABLE_FORMATTING", "UTF8_FULL");
        std::env::set_var("POLARS_FMT_MAX_COLS", "20");
        let df = polars::prelude::LazyCsvReader::new("extra/pos_control.eventalign.txt")
            .with_separator(b'\t')
            .finish()
            .unwrap();
        let tmp_dir = TempDir::new().unwrap();
        let path = tmp_dir.path().join("ipc");
        let df = df
            .with_column(pl::col("samples")
                .str()
                .split(pl::lit(","))
                .list()
                .eval(pl::first().cast(pl::DataType::Float64), false)
                .alias("samples"))
            .group_by([pl::col("read_name"), pl::col("position")])
            .agg([
                pl::col("contig").first(),
                pl::col("model_kmer").first(),
                pl::col("event_length").sum(),
                pl::col("samples").explode(),
                // pl::col("samples").str().split(pl::lit(","))
                // pl::concat_str([pl::col("samples")], ",", true)
                // pl::col("samples").str().split(pl::lit(",")).list().join(pl::lit(","), true)
                                           // pl::concat_list([pl::col("samples")]).unwrap().flatten(),
            ])
            // .with_columns([
            //     pl::col("samples").list().len().alias("n_samples"),
            //     pl::col("samples").list().mean().alias("sample_mean"),
            // ])
            // .sort(["read_name", "position"], Default::default())
            // .group_by_stable([pl::col("read_name"), pl::col("contig")])
            // .agg([pl::all()])
            // .with_columns([
            //     pl::col("position").list().len().alias("n_positions"),
            //     pl::col("position").list().min().alias("start"),
            //     pl::col("position").list().max().alias("end"),
            // ])
            // .with_column((pl::col("end") - pl::col("start")).alias("length"))
            .with_streaming(true);
        let zs = df.clone().explain(true).unwrap();
        println!("{}", zs);
        println!("{:?}", df.clone().collect().unwrap());
        // .sink_parquet(path, Default::default())
        df.sink_ipc(path.clone(), Default::default()).unwrap();
        //     .collect()
        //     .unwrap();
        // let ipc_file = File::open(&path).unwrap();
        // let mut reader = arrow::ipc::reader::FileReader::try_new(ipc_file, None).unwrap();
        // let x = reader.next().unwrap().unwrap();
        // let y: Vec<PolarsEventalign> = serde_arrow::from_record_batch(&x).unwrap();
        // println!("{:?}", y);
        tmp_dir.close().unwrap();
    }
}
