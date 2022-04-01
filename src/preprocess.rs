use std::{collections::hash_map::Entry, fmt::Debug, fs::File, io::Read, path::Path};

use anyhow::Result;
use bio::io::fasta::IndexedReader as GenomeReader;
use bio_types::{genome::AbstractInterval, sequence::SequenceRead};
use indicatif::{ProgressBar, ProgressFinish, ProgressStyle};
use rstats::Stats;
use rust_htslib::bam::{
    ext::BamRecordExtensions,
    record::{Cigar, CigarStringView},
    Read as BamRead, Reader, Record,
};
use serde::{Deserialize, Serialize};
use serde_with::{serde_as, CommaSeparator, StringWithSeparator};

use crate::{
    reads::{LData, LRead, PreprocessRead},
    utils::{self, XHashMap},
};

#[derive(Debug)]
struct QueryLength {
    acc: usize,
}

impl QueryLength {
    fn from_cigar(cigar_string: CigarStringView) -> Self {
        log::debug!("Getting Cigar length");
        let mut acc = 0;

        // let mut citer = cigar_string.iter().peekable();
        // log::debug!("citer length: {}", citer.len());
        // while let Some(cigar) = citer.next_if(|c| matches!(c, Cigar::SoftClip(_))) {
        //     log::debug!("Checking for clips {cigar}");
        //     to_add += cigar.len();
        // }
        cigar_string.iter().for_each(|cigar| {
                log::debug!("cigar {cigar}");
                match cigar {
                    Cigar::SoftClip(_) | Cigar::HardClip(_) | Cigar::Ins(_) => {}
                    _ => acc += cigar.len(),
                }
        });
        // for cigar in citer {
        //     log::debug!("cigar {cigar}");
        //     match cigar {
        //         Cigar::SoftClip(_) | Cigar::HardClip(_) | Cigar::Ins(_) => continue,
        //         _ => acc += cigar.len(),
        //     }
        // }
        QueryLength {
            // to_add: to_add as usize,
            acc: acc as usize,
        }
    }
}

#[derive(PartialEq, Eq, Hash, Debug)]
struct ReadKey(Vec<u8>, String, usize, usize);

type LReadMap = XHashMap<ReadKey, LRead<LData>>;
type ReadKeyMap = XHashMap<(Vec<u8>, String), Vec<(usize, usize)>>;
type ReadChrom = (Vec<u8>, String);
type StartLen = (usize, usize);
type ReadChromPos = (Vec<u8>, String, u64);

fn process_bam<P>(filename: P, faidx: &mut GenomeReader<File>) -> Result<(LReadMap, ReadKeyMap)>
where
    P: AsRef<Path>,
{
    let mut acc: XHashMap<ReadKey, LRead<LData>> = utils::xxhashmap();
    let mut read_keys: XHashMap<ReadChrom, Vec<StartLen>> = utils::xxhashmap();
    let mut bam = Reader::from_path(filename)?;
    let mut record = Record::new();
    // TODO map only unique reads by filtering on mapq or flags
    // TODO should i grab the sequence and store it for later?
    while let Some(result) = bam.read(&mut record) {
        if result.is_ok() {
            let name = record.name().to_owned();
            log::debug!("Read name {}", String::from_utf8_lossy(&name));

            let chrom = record.contig().to_owned();
            log::debug!("chrom {chrom}");

            let qlen = QueryLength::from_cigar(record.cigar());
            log::debug!("QueryLength {qlen:?}");

            let start = record.reference_start() as usize;
            log::debug!("start {start}");

            let length = qlen.acc;
            log::debug!("length {length}");

            faidx.fetch(&chrom, start as u64, (start + length) as u64)?;
            let mut seq = Vec::new();
            faidx.read(&mut seq)?;
            let data = Vec::new();

            let bam_key = ReadKey(name.clone(), chrom.clone(), start, length);
            let read_key = (name.clone(), chrom.clone());
            let read = PreprocessRead::new(name, chrom, start, length, seq, data);

            let keys = read_keys.entry(read_key).or_default();
            keys.push((start, length));

            if let Some(ld) = acc.insert(bam_key, read) {
                log::warn!("Duplicate read of encountered in bam with data {:?}", ld);
            }
        }
    }
    Ok((acc, read_keys))
}

#[serde_as]
#[derive(Clone, Debug, Serialize, Deserialize)]
struct Npr {
    contig: String,

    position: u64,

    #[serde(skip)]
    _reference_kmer: String,

    read_name: String,

    #[serde(skip)]
    _strand: String,

    #[serde(skip)]
    _event_index: u64,

    #[serde(skip)]
    _event_level_mean: f64,

    #[serde(skip)]
    _event_stdv: f64,

    event_length: f64,

    model_kmer: String,

    #[serde(skip)]
    _model_mean: f64,

    #[serde(skip)]
    _model_stdv: f64,

    #[serde(skip)]
    _standardized_level: f64,

    #[serde_as(as = "StringWithSeparator::<CommaSeparator, f64>")]
    samples: Vec<f64>,
}

pub(crate) struct Process {
    chrom: Option<String>,
    start: Option<u64>,
    stop: Option<u64>,
}

// TODO: Store ProgressBar in process and add indicatif support across
/// Reads a nanopolish output file, collapses samples across multiple lines,
/// finds the mean for all the samples, and produces an Output for
/// each collapsed line.
impl Process {
    pub(crate) fn new() -> Self {
        Self {
            chrom: None,
            start: None,
            stop: None,
        }
    }

    pub(crate) fn chrom(mut self, chrom: Option<String>) -> Self {
        self.chrom = chrom;
        self
    }

    pub(crate) fn start(mut self, start: Option<u64>) -> Self {
        self.start = start;
        self
    }

    pub(crate) fn stop(mut self, stop: Option<u64>) -> Self {
        self.stop = stop;
        self
    }

    pub(crate) fn run<P>(
        &self,
        filename: P,
        bam_filename: P,
        genome: P,
    ) -> Result<Vec<PreprocessRead>>
    where
        P: AsRef<Path> + Debug,
    {
        let mut faidx = GenomeReader::from_file(&genome)?;
        let (mut bam_to_pr, read_keys) = process_bam(bam_filename, &mut faidx)?;
        log::debug!("bam map length {}", bam_to_pr.len());
        log::debug!("read map length {}", read_keys.len());
        let rp_to_samples = self.with_file(filename)?;
        for (key, datum) in rp_to_samples.into_iter() {
            if let Some(read_start_lens) = read_keys.get(&(key.0.clone(), key.1.clone())) {
                if !read_start_lens.is_empty() {
                    // Since reads can have multiple primary alignments and there is a potential for
                    // the primary alignments to be on the same strand,
                    let pos = datum.pos as usize;
                    let within_read = read_start_lens.iter().filter(|&&(read_start, read_len)| {
                        let stop = read_start + read_len;
                        log::debug!("{read_start} <= {pos} <= {stop}");
                        (read_start <= pos) && (pos <= stop)
                    });
                    for &(read_start, read_len) in within_read {
                        let key = ReadKey(key.0.clone(), key.1.clone(), read_start, read_len);
                        let mean = datum.samples.amean()?;
                        let data = LData::new(datum.pos, datum.kmer.clone(), mean, datum.time);
                        bam_to_pr
                            .get_mut(&key)
                            .expect("Read key not in bam map {key:?}")
                            .get_mut_data()
                            .push(data);
                    }
                } else {
                    log::warn!("Read found with no starts and lengths!");
                }
            }
        }
        log::debug!("{bam_to_pr:?}");
        Ok(bam_to_pr
            .into_values()
            .filter(|lr| !lr.is_empty())
            .collect())
    }

    pub(crate) fn with_file<P>(&self, filename: P) -> Result<XHashMap<ReadChromPos, Data>>
    where
        P: AsRef<Path>,
    {
        let file = File::open(filename)?;
        let style = ProgressStyle::default_spinner()
            .template("{spinner} [{elapsed_precise}] {binary_bytes_per_sec} {msg}")
            .on_finish(ProgressFinish::AndLeave);
        let file = ProgressBar::new_spinner()
            .with_message("Processing data")
            .with_style(style)
            .wrap_read(file);
        self.with_reader(file)
    }

    pub(crate) fn with_reader<R>(&self, input: R) -> Result<XHashMap<ReadChromPos, Data>>
    where
        R: Read,
    {
        let mut acc: XHashMap<_, Data> = utils::xxhashmap();
        let mut reader = csv::ReaderBuilder::new()
            .delimiter(b'\t')
            .from_reader(input);
        for line in reader.deserialize() {
            let line = line?;
            log::debug!("csv line {:?}", line);
            let line: Npr = line;

            // Filter records that are not on chromosome
            if let Some(ref chrom) = self.chrom {
                if chrom == &line.contig {
                    continue;
                }
            }
            if let Some(start) = self.start {
                if line.position < start {
                    continue;
                }
            }
            if let Some(stop) = self.stop {
                if line.position > stop {
                    continue;
                }
            }

            let key = (
                line.read_name.as_bytes().to_owned(),
                line.contig.clone(),
                line.position,
            );
            match acc.entry(key) {
                Entry::Occupied(mut oe) => {
                    oe.get_mut().update_from_npr(line);
                }
                Entry::Vacant(ve) => {
                    let d = Data::new_from_npr(line);
                    ve.insert(d);
                }
            }
        }
        Ok(acc)
    }
}

#[derive(Debug)]
pub(crate) struct Data {
    pos: u64,
    kmer: String,
    samples: Vec<f64>,
    time: f64,
}

impl Data {
    fn new(pos: u64, kmer: String, samples: Vec<f64>, time: f64) -> Self {
        Self {
            pos,
            kmer,
            samples,
            time,
        }
    }

    fn new_from_npr(npr: Npr) -> Self {
        let samples = npr.samples;
        let pos = npr.position;
        let kmer = npr.model_kmer;
        let time = npr.event_length;
        Self::new(pos, kmer, samples, time)
    }

    fn update_from_npr(&mut self, mut npr: Npr) {
        self.samples.append(&mut npr.samples);
        self.time += npr.event_length;
    }
}

#[cfg(test)]
mod test_super {
    use assert_fs::TempDir;

    use super::*;
    use crate::{preprocess::Process, utils::CawlrIO};

    #[test_log::test]
    fn test_save_load_preprocess() -> Result<()> {
        let temp_dir = TempDir::new()?;
        let output = temp_dir.path().join("output.pickle");

        let input = "extra/single_read.eventalign.txt";
        let bam = "extra/single_read.bam";
        let genome = "./extra/sacCer3.fa";

        let nprs = Process::new().run(input, bam, genome)?;
        nprs.clone().save(output.clone())?;

        let loaded: Vec<LRead<LData>> = CawlrIO::load(output)?;

        log::debug!("save len {}", nprs.len());
        log::debug!("loaded len {}", loaded.len());
        assert_eq!(nprs.len(), loaded.len());
        assert!(!nprs[0].data().is_empty());
        assert!(!loaded[0].data().is_empty());
        Ok(())
    }
}
