use std::{
    fs::{self, File},
    io::{BufRead, BufReader, BufWriter, Write},
    path::{Path, PathBuf},
    process::{Command, Stdio},
};

use clap::Parser;
use eyre::Context;
use libcawlr::utils::{self, check_if_failed};
use log::LevelFilter;

use crate::file::ValidPathBuf;

#[derive(Parser, Debug)]
pub struct PreprocessCmd {
    /// Path to genome fasta file
    #[clap(short, long)]
    pub genome: ValidPathBuf,

    /// Path to fastq or set of fastqs. If a directory is given, the command
    /// will concatenate all .fastq files and write to the output directory.
    /// If a single file is written, it will create a symlink to that file
    /// in the output directory to maintain consistency.
    #[clap(long)]
    pub reads: ValidPathBuf,

    /// Path to a directory of fast5 files
    #[clap(long)]
    pub fast5: ValidPathBuf,

    /// Path to a sequencing_summary.txt file, will speed up nanopolish index
    #[clap(long)]
    pub summary: Option<ValidPathBuf>,

    #[clap(short, long)]
    pub output_dir: PathBuf,

    /// Path to minimap2 binary, if not specified will look in $PATH
    #[clap(long)]
    pub minimap2_path: Option<PathBuf>,

    /// Path to nanopolish binary, if not specified will look in $PATH
    #[clap(long)]
    pub nanopolish_path: Option<PathBuf>,

    /// Path to samtools binary, if not specified will look in $PATH
    #[clap(long)]
    pub samtools_path: Option<PathBuf>,

    #[clap(long, default_value_t = false)]
    pub overwrite: bool,

    #[clap(short = 'j', long, default_value_t = 4)]
    pub n_threads: usize,
}

impl PreprocessCmd {
    pub fn run(self) -> eyre::Result<()> {
        if self.overwrite && self.output_dir.exists() {
            fs::remove_dir_all(&self.output_dir)?;
        }
        fs::create_dir_all(&self.output_dir)?;

        let log_file_path = self.output_dir.join("log.txt");
        let log_file = File::create(log_file_path)?;
        simple_logging::log_to(log_file.try_clone()?, LevelFilter::Info);

        log::info!("{self:?}");
        let reads = self.reads_to_single_reads("reads.fastq")?;
        self.aln_reads(&reads, log_file.try_clone()?)?;
        self.np_index(&reads, log_file.try_clone()?)?;
        Ok(())
    }

    fn np_index(&self, reads: &Path, log_file: File) -> eyre::Result<()> {
        let nanopolish = utils::find_binary("nanopolish", &self.nanopolish_path)?;
        let mut cmd = Command::new(nanopolish);
        cmd.arg("index").arg("-d").arg(&self.fast5);
        if let Some(ref summary) = self.summary {
            cmd.arg("-s").arg(summary);
        }
        cmd.arg(reads);
        cmd.stderr(log_file);
        log::info!("{cmd:?}");
        let output = cmd.output()?;
        check_if_failed(output).wrap_err("nanopolish index failed")
    }

    fn aln_reads(&self, reads: &Path, log_file: File) -> eyre::Result<()> {
        let minimap2 = utils::find_binary("minimap2", &self.minimap2_path)?;
        let samtools = utils::find_binary("samtools", &self.samtools_path)?;
        let mut map_cmd = Command::new(minimap2);
        map_cmd
            .arg("-ax")
            .arg("map-ont")
            .arg("--sam-hit-only")
            .arg("--secondary=no")
            .args(["-t", "4"])
            .arg(&self.genome)
            .arg(reads)
            .stderr(log_file.try_clone()?)
            .stdout(Stdio::piped());
        log::info!("{map_cmd:?}");
        let map_output = map_cmd.spawn()?;

        let mut sam_cmd = Command::new(samtools);
        let aln_bam = self.output_dir.join("aln.bam");
        sam_cmd
            .arg("sort")
            .arg("--write-index")
            .arg("-T")
            .arg(&self.output_dir)
            .arg("-o")
            .arg(aln_bam)
            .stderr(log_file.try_clone()?)
            .stdin(map_output.stdout.unwrap());
        log::info!("{sam_cmd:?}");
        let output = sam_cmd.output()?;

        check_if_failed(output).wrap_err("minimap2 | samtools failed")
    }

    fn reads_to_single_reads(&self, name: &str) -> eyre::Result<PathBuf> {
        let output_filepath = self.output_dir.join(name);
        if self.reads.0.is_dir() {
            log::info!("Detected directory, concatenating into a single fastq file.");
            let mut output_file = BufWriter::new(File::create(&output_filepath)?);
            let fastq_matcher = format!(
                "{}/**/*fastq",
                self.reads.0.as_os_str().to_str().ok_or(eyre::eyre!(
                    "Failed to convert path into str, unicdoe issue?"
                ))?
            );
            let mut n_fastq_files = 0;
            for fastq in glob::glob(&fastq_matcher)? {
                let fastq = fastq?;
                n_fastq_files += 1;
                log::info!("Found fastq: {}", fastq.display());
                let mut fastq_file = BufReader::new(File::open(fastq)?);
                loop {
                    let buf_len = {
                        let buf = fastq_file.fill_buf()?;
                        if buf.is_empty() {
                            break;
                        }
                        output_file.write_all(buf)?;
                        buf.len()
                    };
                    fastq_file.consume(buf_len);
                }
            }

            if n_fastq_files == 0 {
                return Err(eyre::eyre!(
                    "No fastq files processed, check if directory contained files ending with .fastq"
                ));
            } else {
                log::info!("Processed {n_fastq_files} fastq files");
            }
        } else {
            std::os::unix::fs::symlink(&self.reads.0, &output_filepath)?;
        }
        Ok(output_filepath)
    }
}
