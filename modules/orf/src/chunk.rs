//! Core module for detecting open reading frames in a query set of reads
//! Alejandro Gonzales-Irribarren, 2025
//!
//! This module contains the main functions for finding open reading frames (ORFs)
//! in a set of aligned reads.
//!
//! In short, every possible open reading frame (ORF) is detected for every
//! read in the query set. For every potential ORF, learning models and databases
//! are used to determine whether the ORF is a true ORF, a false positive.
//! All the data from each reliable ORF is collected and subjected to another
//! learning model trained with true ORFs and false positives. The process is
//! heavily parallelized to offer fast performance on large datasets.

use flate2::read::MultiGzDecoder;
use genepred::{Bed12, GenePred, Gff, Gtf, Reader, ReaderResult, Strand, Writer, bed::BedFormat};
use rayon::prelude::*;
use twobit::TwoBitFile;

use std::{
    collections::HashMap,
    fmt::Debug,
    fs::{File, create_dir_all},
    io::{BufRead, BufReader, BufWriter, Write},
    path::{Path, PathBuf},
    sync::Arc,
};

use crate::cli::ChunkArgs;

pub fn run_chunk(args: ChunkArgs) {
    let genome = get_sequences(args.sequence);
    let outdir = args.outdir.join("tmp");
    create_dir_all(&outdir).unwrap_or_else(|e| panic!("{}", e));

    let prefix = Arc::new(args.prefix.unwrap_or_default());

    match detect_region_format(&args.regions) {
        Some(RegionFormat::Bed) => process_reader::<Bed12>(
            &args.regions,
            args.chunks,
            &outdir,
            &genome,
            args.upstream_flank,
            args.downstream_flank,
            args.ignore_errors,
            prefix,
        ),
        Some(RegionFormat::Gtf) => process_reader::<Gtf>(
            &args.regions,
            args.chunks,
            &outdir,
            &genome,
            args.upstream_flank,
            args.downstream_flank,
            args.ignore_errors,
            prefix,
        ),
        Some(RegionFormat::Gff) => process_reader::<Gff>(
            &args.regions,
            args.chunks,
            &outdir,
            &genome,
            args.upstream_flank,
            args.downstream_flank,
            args.ignore_errors,
            prefix,
        ),
        None => panic!("ERROR: Unsupported file format"),
    }
}

/// Processes genomic regions in parallel chunks and writes output to FASTA file.
///
/// # Arguments
///
/// - `regions`: Path to the annotation file (BED, GTF, or GFF)
/// - `chunks`: Number of records per parallel processing chunk
/// - `outdir`: Output directory path
/// - `genome`: HashMap of chromosome names to sequences
/// - `upstream_flank`: Bases to extend upstream of first exon
/// - `downstream_flank`: Bases to extend downstream of last exon
/// - `ignore_errors`: Whether to continue on errors
/// - `prefix`: Prefix for output file names
///
/// # Example
///
/// ```rust,ignore
/// use xloci::core::process_reader;
/// use xloci::Feature;
/// use std::collections::HashMap;
///
/// let genome: HashMap<Vec<u8>, Vec<u8>> = HashMap::new();
/// process_reader::<genepred::Gtf>(
///     ...
/// );
/// ```
#[allow(clippy::too_many_arguments)]
fn process_reader<R>(
    regions: &Path,
    chunks: usize,
    outdir: &Path,
    genome: &HashMap<Vec<u8>, Vec<u8>>,
    upstream_flank: usize,
    downstream_flank: usize,
    ignore_errors: bool,
    prefix: Arc<String>,
) where
    R: BedFormat + Into<GenePred> + Send,
{
    Reader::<R>::from_mmap(regions)
        .unwrap_or_else(|e| panic!("{}", e))
        .par_chunks(chunks)
        .unwrap_or_else(|e| panic!("{}", e))
        .for_each(|(idx, chunk)| {
            write_chunk(
                idx,
                chunk,
                genome,
                outdir,
                upstream_flank,
                downstream_flank,
                ignore_errors,
                prefix.clone(),
            )
        });
}

/// Processes a chunk of genomic records and writes extracted sequences to a temporary file.
///
/// # Arguments
///
/// - `idx`: Chunk index for naming output files
/// - `chunk`: Vector of GenePred records to process
/// - `genome`: HashMap of chromosome names to sequences
/// - `upstream_flank`: Bases to extend upstream of first exon
/// - `downstream_flank`: Bases to extend downstream of last exon
/// - `ignore_errors`: Whether to continue on errors
/// - `outdir`: Output directory for chunk files
///
/// # Example
///
/// ```rust,ignore
/// use std::sync::{Arc, Mutex};
/// use xloci::Feature;
///
/// let collector = Some(Arc::new(Mutex::new(Vec::new())));
/// write_chunk(
///     ...
/// );
/// ```
#[allow(clippy::too_many_arguments)]
fn write_chunk(
    idx: usize,
    chunk: Vec<ReaderResult<GenePred>>,
    genome: &HashMap<Vec<u8>, Vec<u8>>,
    outdir: &Path,
    upstream_flank: usize,
    downstream_flank: usize,
    ignore_errors: bool,
    prefix: Arc<String>,
) {
    let tmp = if !prefix.is_empty() {
        outdir.join(format!("{}.tmp_{}.bed", prefix, idx))
    } else {
        outdir.join(format!("tmp_{}.bed", idx))
    };

    let mut writer = BufWriter::new(File::create(&tmp).unwrap_or_else(|e| panic!("{}", e)));

    let mut f_writer =
        BufWriter::new(File::create(tmp.with_extension("fa")).unwrap_or_else(|e| panic!("{}", e)));

    let rescue_dir = outdir.join("rescue");
    let rescue_bed_path = rescue_dir.join(tmp.file_name().unwrap());
    let rescue_fa_path = rescue_bed_path.with_extension("fa");

    // INFO: lazily created only if at least one record needs a flank rescue
    let mut rescue_writer: Option<BufWriter<File>> = None;
    let mut rescue_f_writer: Option<BufWriter<File>> = None;

    chunk
        .into_iter()
        .filter_map(|result| match result {
            Ok(record) => Some(record),
            Err(e) => {
                eprintln!("WARN: Failed to process record: {}", e);
                None
            }
        })
        .for_each(|record| {
            let seq = genome.get(&record.chrom).unwrap_or_else(|| {
                panic!(
                    "ERROR: Chromosome {} not found!",
                    std::str::from_utf8(&record.chrom).unwrap()
                )
            });

            let target = extract_seq(
                &record,
                seq,
                upstream_flank,
                downstream_flank,
                ignore_errors,
            );

            match target {
                Some(target) => write_record(&record, target, &mut writer, &mut f_writer),
                None => {
                    // INFO: rescue: flanks could not be extracted, retry with 0 bp flanks
                    if let Some(target) = extract_seq(&record, seq, 0, 0, ignore_errors) {
                        eprintln!(
                            "WARN: flanks out of bounds, rescuing record {} with 0 bp flanks",
                            record
                        );

                        if rescue_writer.is_none() {
                            create_dir_all(&rescue_dir).unwrap_or_else(|e| panic!("{}", e));
                            rescue_writer = Some(BufWriter::new(
                                File::create(&rescue_bed_path).unwrap_or_else(|e| panic!("{}", e)),
                            ));
                            rescue_f_writer = Some(BufWriter::new(
                                File::create(&rescue_fa_path).unwrap_or_else(|e| panic!("{}", e)),
                            ));
                        }

                        write_record(
                            &record,
                            target,
                            rescue_writer.as_mut().unwrap(),
                            rescue_f_writer.as_mut().unwrap(),
                        );
                    }
                }
            }
        });

    // INFO: flush first so data is written to the file
    writer.flush().unwrap_or_else(|e| panic!("{}", e));
    f_writer.flush().unwrap_or_else(|e| panic!("{}", e));

    // INFO: check actual file size via metadata
    let bed_path = tmp.clone();
    let file_len = std::fs::metadata(&bed_path).map(|m| m.len()).unwrap_or(0);
    if file_len == 0 {
        std::fs::remove_file(&bed_path).unwrap_or_else(|e| panic!("{}", e));
        std::fs::remove_file(bed_path.with_extension("fa")).unwrap_or_else(|e| panic!("{}", e));
    }

    if let Some(mut r_writer) = rescue_writer {
        r_writer.flush().unwrap_or_else(|e| panic!("{}", e));
        if let Some(mut r_f_writer) = rescue_f_writer {
            r_f_writer.flush().unwrap_or_else(|e| panic!("{}", e));
        }

        let rescue_len = std::fs::metadata(&rescue_bed_path)
            .map(|m| m.len())
            .unwrap_or(0);
        if rescue_len == 0 {
            std::fs::remove_file(&rescue_bed_path).unwrap_or_else(|e| panic!("{}", e));
            std::fs::remove_file(&rescue_fa_path).unwrap_or_else(|e| panic!("{}", e));
            std::fs::remove_dir(&rescue_dir).unwrap_or_else(|e| panic!("{}", e));
        }
    }
}

/// Writes a record and its extracted sequence to the BED and FASTA writers,
/// reverse-complementing the sequence for genes on the reverse strand.
fn write_record<W: Write, F: Write>(
    record: &GenePred,
    mut target: Vec<u8>,
    writer: &mut W,
    f_writer: &mut F,
) {
    match &record.strand {
        Some(Strand::Forward) => {}
        Some(Strand::Reverse) => {
            target.reverse();

            for base in target.iter_mut() {
                *base = match *base {
                    b'A' => b'T',
                    b'C' => b'G',
                    b'G' => b'C',
                    b'T' => b'A',
                    b'N' => b'N',
                    b'a' => b't',
                    b'c' => b'g',
                    b'g' => b'c',
                    b't' => b'a',
                    b'n' => b'n',
                    _ => panic!("ERROR: Invalid base"),
                }
            }
        }
        Some(Strand::Unknown) | None => {}
    }

    Writer::<Bed12>::from_record(record, writer)
        .unwrap_or_else(|e| panic!("ERROR: Cannot write record to file: {}", e));

    f_writer.write_all(b">").unwrap_or_else(|e| panic!("{}", e));
    f_writer
        .write_all(record.name().unwrap())
        .unwrap_or_else(|e| panic!("{}", e));
    f_writer
        .write_all(b"\n")
        .unwrap_or_else(|e| panic!("{}", e));
    f_writer
        .write_all(&target)
        .unwrap_or_else(|e| panic!("{}", e));
    f_writer
        .write_all(b"\n")
        .unwrap_or_else(|e| panic!("{}", e));
}

#[derive(Debug)]
#[allow(dead_code)]
enum RangeError {
    Underflow { feature_coord: usize, flank: usize },
}

impl std::fmt::Display for RangeError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            RangeError::Underflow {
                feature_coord,
                flank,
            } => write!(
                f,
                "ERROR: Feature coordinate {} is underflowing by {} bases",
                feature_coord, flank
            ),
        }
    }
}

/// Calculates the slice range for an exon with appropriate flanking regions.
///
/// # Arguments
///
/// - `exon_idx`: Index of the current exon (0-based)
/// - `exon_count`: Total number of exons in the feature
/// - `feature_start`: Start coordinate of the exon
/// - `feature_end`: End coordinate of the exon
/// - `upstream_flank`: Bases to add to first exon start
/// - `downstream_flank`: Bases to add to last exon end
///
/// # Example
///
/// ```rust,ignore
/// // Single exon with flanking
/// let range = slice_range_for_exon(0, 1, 100, 200, 10, 20).unwrap();
/// assert_eq!(range, 90..220);
///
/// // Middle exon (no flanking)
/// let range = slice_range_for_exon(1, 3, 100, 200, 10, 20).unwrap();
/// assert_eq!(range, 100..200);
/// ```
fn slice_range_for_exon(
    exon_idx: usize,
    exon_count: usize,
    exon_start: usize,
    exon_end: usize,
    upstream_flank: usize,
    downstream_flank: usize,
) -> Result<std::ops::Range<usize>, RangeError> {
    let is_first = exon_idx == 0;
    let is_last = exon_idx + 1 == exon_count;

    let start = if is_first {
        exon_start
            .checked_sub(upstream_flank)
            .ok_or(RangeError::Underflow {
                feature_coord: exon_start,
                flank: upstream_flank,
            })?
    } else {
        exon_start
    };

    let end = if is_last {
        exon_end
            .checked_add(downstream_flank)
            .ok_or(RangeError::Underflow {
                feature_coord: exon_end,
                flank: downstream_flank,
            })?
    } else {
        exon_end
    };

    Ok(start..end)
}

/// Extends target sequence with a slice or handles out-of-bounds errors gracefully.
///
/// # Arguments
///
/// - `record`: The GenePred record being processed
/// - `target`: Vector to extend with the slice
/// - `seq`: Full chromosome sequence
/// - `range`: Slice range to extract
/// - `exon_idx`: Index of current exon for error messages
/// - `ignore_errors`: Whether to return false on error instead of panicking
///
/// # Example
///
/// ```rust,ignore
/// let seq = b"ACGTACGTACGT";
/// let mut target = Vec::new();
/// let success = extend_or_handle_oob(&record, &mut target, seq, 0..6, 0, true);
/// assert!(success);
/// assert_eq!(target, b"ACGTAC");
/// ```
fn extend_or_handle_oob(
    record: &GenePred,
    target: &mut Vec<u8>,
    seq: &[u8],
    range: std::ops::Range<usize>,
    exon_idx: usize,
    ignore_errors: bool,
) -> bool {
    if let Some(slice) = seq.get(range.clone()) {
        target.extend_from_slice(slice);
        return true;
    }

    if ignore_errors {
        eprintln!(
            "WARN: out-of-bounds slice for {} exon {}: {:?} (seq_len={})",
            record,
            exon_idx,
            range,
            seq.len()
        );
        false
    } else {
        panic!(
            "ERROR: out-of-bounds slice for {} exon {}: {:?} (seq_len={})",
            record,
            exon_idx,
            range,
            seq.len()
        );
    }
}

/// Extracts genomic sequence for a feature with flanking regions applied.
///
/// # Arguments
///
/// - `record`: GenePred record containing exon coordinates
/// - `seq`: Full chromosome sequence
/// - `upstream_flank`: Bases to extend upstream of first feature
/// - `downstream_flank`: Bases to extend downstream of last feature
/// - `feature_type`: Type of feature to extract (exon, intron, CDS, etc.)
/// - `ignore_errors`: Whether to return None on error instead of panicking
///
/// # Example
///
/// ```rust,ignore
/// let seq = b"ACGTACGTACGT";
/// let extracted = extract_seq(&record, seq, 0, 0, &Feature::Exon, false);
/// assert_eq!(extracted, Some(b"ACGT".to_vec()));
/// ```
fn extract_seq(
    record: &GenePred,
    seq: &[u8],
    upstream_flank: usize,
    downstream_flank: usize,
    ignore_errors: bool,
) -> Option<Vec<u8>> {
    let exons = record.exons();
    let exon_count = exons.len();
    let mut target = Vec::new();

    for (idx, (start, end)) in exons.iter().enumerate() {
        let exon_start = *start as usize;
        let exon_end = *end as usize;

        let range = match slice_range_for_exon(
            idx,
            exon_count,
            exon_start,
            exon_end,
            upstream_flank,
            downstream_flank,
        ) {
            Ok(r) => r,
            Err(err) => {
                if ignore_errors {
                    eprintln!("WARN: {:?} for record {}", err, record);
                    return None;
                } else {
                    panic!("ERROR: {:?} for record {}", err, record);
                }
            }
        };

        if !extend_or_handle_oob(record, &mut target, seq, range, idx, ignore_errors) {
            return None;
        }
    }

    Some(target)
}

/// Supported genomic annotation file formats.
///
/// # Variants
///
/// - `Bed`: BED format (12-column)
/// - `Gtf`: GTF format
/// - `Gff`: GFF format
///
/// # Example
///
/// ```rust,ignore
/// use xloci::core::RegionFormat;
///
/// let format = detect_region_format(Path::new("annotations.gtf"));
/// assert_eq!(format, Some(RegionFormat::Gtf));
/// ```
#[derive(Clone, Copy)]
enum RegionFormat {
    Bed,
    Gtf,
    Gff,
}

/// Detects the genomic annotation format from file extension.
///
/// # Arguments
///
/// - `path`: Path to the annotation file
///
/// # Example
///
/// ```rust,ignore
/// use std::path::Path;
///
/// assert_eq!(detect_region_format(Path::new("file.bed")), Some(RegionFormat::Bed));
/// assert_eq!(detect_region_format(Path::new("file.gtf")), Some(RegionFormat::Gtf));
/// assert_eq!(detect_region_format(Path::new("file.gff")), Some(RegionFormat::Gff));
/// assert_eq!(detect_region_format(Path::new("file.gtf.gz")), Some(RegionFormat::Gtf));
/// assert_eq!(detect_region_format(Path::new("file.txt")), None);
/// ```
fn detect_region_format(path: &Path) -> Option<RegionFormat> {
    match path.extension().and_then(|ext| ext.to_str()) {
        Some("bed") => Some(RegionFormat::Bed),
        Some("gtf") => Some(RegionFormat::Gtf),
        Some("gff") => Some(RegionFormat::Gff),
        Some("gz") => {
            let stem = path.file_stem()?.to_str()?;
            if stem.ends_with(".bed") {
                Some(RegionFormat::Bed)
            } else if stem.ends_with(".gtf") {
                Some(RegionFormat::Gtf)
            } else if stem.ends_with(".gff") {
                Some(RegionFormat::Gff)
            } else {
                None
            }
        }
        _ => None,
    }
}

/// Loads genome sequences from a file (2bit or FASTA format).
///
/// # Arguments
///
/// - `sequence`: Path to the genome file (.fa, .fa.gz, or .2bit)
///
/// # Example
///
/// ```rust,ignore
/// use std::path::PathBuf;
///
/// let genome = get_sequences(PathBuf::from("genome.2bit"));
/// let genome = get_sequences(PathBuf::from("genome.fa"));
/// let genome = get_sequences(PathBuf::from("genome.fa.gz"));
/// ```
pub fn get_sequences(sequence: PathBuf) -> HashMap<Vec<u8>, Vec<u8>> {
    match sequence.extension() {
        Some(ext) => match ext.to_str() {
            Some("2bit") => from_2bit(sequence),
            Some("fa") | Some("fasta") | Some("fna") | Some("gz") => from_fa(sequence),
            _ => panic!("ERROR: Unsupported file format"),
        },
        None => panic!("ERROR: No file extension"),
    }
}

/// Loads genome sequences from a 2bit compressed format file.
///
/// # Arguments
///
/// - `twobit`: Path to the 2bit file
///
/// # Example
///
/// ```rust,ignore
/// use std::path::PathBuf;
///
/// let sequences = from_2bit(PathBuf::from("genome.2bit"));
/// let chr1 = sequences.get(b"chr1");
/// ```
fn from_2bit(twobit: PathBuf) -> HashMap<Vec<u8>, Vec<u8>> {
    let mut genome = TwoBitFile::open_and_read(twobit).expect("ERROR: Cannot open 2bit file");

    let mut sequences = HashMap::new();
    genome.chrom_names().iter().for_each(|chr| {
        let seq = genome
            .read_sequence(chr, ..)
            .unwrap_or_else(|e| panic!("ERROR: {}", e))
            .as_bytes()
            .to_vec();

        sequences.insert(chr.as_bytes().to_vec(), seq);
    });

    sequences
}

/// Loads genome sequences from a FASTA format file (optionally gzipped).
///
/// # Arguments
///
/// - `f`: Path to the FASTA file (.fa or .fa.gz)
///
/// # Example
///
/// ```rust,ignore
/// use std::path::PathBuf;
///
/// let sequences = from_fa(PathBuf::from("genome.fa"));
/// let sequences = from_fa(PathBuf::from("genome.fa.gz"));
/// let chr1 = sequences.get(b"chr1");
/// ```
pub fn from_fa<F: AsRef<Path> + Debug>(f: F) -> HashMap<Vec<u8>, Vec<u8>> {
    let path = f.as_ref();
    let file = File::open(path)
        .unwrap_or_else(|e| panic!("ERROR: cannot open FASTA {}: {}", path.display(), e));

    let mut reader: Box<dyn BufRead> = match path.extension().and_then(|ext| ext.to_str()) {
        Some("gz") => Box::new(BufReader::new(MultiGzDecoder::new(file))),
        _ => Box::new(BufReader::new(file)),
    };

    let mut acc = HashMap::new();
    let mut line = Vec::new();
    let mut header: Option<Vec<u8>> = None;
    let mut seq = Vec::new();

    loop {
        line.clear();
        let bytes_read = reader
            .read_until(b'\n', &mut line)
            .unwrap_or_else(|e| panic!("ERROR: cannot read FASTA {}: {}", path.display(), e));

        if bytes_read == 0 {
            break;
        }

        if line.ends_with(b"\n") {
            line.pop();
        }

        if line.ends_with(b"\r") {
            line.pop();
        }

        if line.is_empty() {
            continue;
        }

        if line[0] == b'>' {
            if let Some(prev_header) = header.replace(line[1..].to_vec()) {
                acc.insert(prev_header, std::mem::take(&mut seq));
            }
        } else {
            seq.extend_from_slice(&line);
        }
    }

    if let Some(last_header) = header {
        acc.insert(last_header, seq);
    }

    acc
}

#[cfg(test)]
mod tests {
    use super::*;

    fn record_from_bed12(data: String) -> GenePred {
        let mut reader: Reader<Bed12> = Reader::from_reader(std::io::Cursor::new(data)).unwrap();
        reader.next().unwrap().unwrap()
    }

    #[test]
    fn test_extract_seq_flank_underflow_rescued_with_zero_flanks() {
        // INFO: transcript close to the chromosome start so the 1000 bp
        // upstream flank underflows and extraction fails
        let record =
            record_from_bed12("chr1\t100\t500\ttx1\t0\t+\t100\t500\t0,0,0\t1\t400\t0".to_string());

        let seq = vec![b'N'; 1000];

        assert!(extract_seq(&record, &seq, 1000, 1000, true).is_none());

        let rescued = extract_seq(&record, &seq, 0, 0, true).unwrap();
        assert_eq!(rescued.len(), 400);
    }

    #[test]
    fn test_extract_seq_flank_out_of_bounds_rescued_with_zero_flanks() {
        // INFO: transcript close to the chromosome end so the 1000 bp
        // downstream flank exceeds the sequence length
        let record =
            record_from_bed12("chr1\t900\t1000\ttx2\t0\t+\t900\t1000\t0,0,0\t1\t100\t0".to_string());

        let seq = vec![b'N'; 1000];

        assert!(extract_seq(&record, &seq, 1000, 1000, true).is_none());

        let rescued = extract_seq(&record, &seq, 0, 0, true).unwrap();
        assert_eq!(rescued.len(), 100);
    }

    #[test]
    fn test_extract_seq_with_flanks() {
        let record =
            record_from_bed12("chr1\t1100\t1500\ttx3\t0\t+\t1100\t1500\t0,0,0\t1\t400\t0".to_string());

        let seq = vec![b'N'; 10000];

        let extracted = extract_seq(&record, &seq, 1000, 1000, true).unwrap();
        assert_eq!(extracted.len(), 2400);
    }
}
