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

use genepred::GenePred;
use hashbrown::HashMap;

use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

use crate::{cli::NetArgs, utils::*};

/// Run nets
///
/// # Arguments
///
///
/// * `args` - CLI arguments
///
/// # Returns
///
/// * `()` - Writes net files to disk
///
/// # Panics
///
/// * If the directory for the net files cannot be created
///
/// # Examples
///
/// ```rust, ignore
/// use orf::cli::NetArgs;
/// use orf::nets::run_nets;
///
/// let args = NetArgs {
///     bed: String::from("tests/data/test.bed"),
///     netstart: String::from("tests/data/test.netstart"),
///     transaid: String::from("tests/data/test.transaid"),
///     outdir: String::from("tests/data/test.outdir"),
/// };
///
/// run_nets(args);
/// ```
pub fn run_nets(args: NetArgs) {
    let dir = args.outdir.join("net");
    std::fs::create_dir_all(&dir)
        .unwrap_or_else(|e| panic!("ERROR: could not create directory -> {e}"));

    let net = dir.join("merged.net");

    let netstart = net_map(args.netstart, NetSource::Netstart);
    let transaid = net_map(args.transaid, NetSource::Transaid);
    let bed = get_bed(&args.bed);

    __join_nets(netstart, transaid, bed, &net);
}

fn __join_nets(
    netstart: HashMap<String, Vec<NetRecord>>,
    mut transaid: HashMap<String, Vec<NetRecord>>,
    bed: HashMap<String, GenePred>,
    outfile: &Path,
) {
    let mut writer = BufWriter::new(
        File::create(outfile).unwrap_or_else(|e| panic!("ERROR: cannot create file -> {e}")),
    );

    netstart.into_iter().for_each(|(id, prediction)| {
        let record = bed.get(&id).unwrap_or_else(|| {
            panic!("ERROR: could not find genepred record for id: {id:?}!");
        });

        for prediction in prediction {
            let prediction = match prediction {
                NetRecord::NetNS(record) => record,
                _ => panic!("ERROR: unexpected record type"),
            };

            let mut found = false;

            let key = format!(
                "{}#{}#{}",
                id, prediction.orf_relative_start, prediction.orf_relative_stop
            );

            let (orf_start, orf_end) = map_absolute_cds(
                record,
                prediction.orf_relative_start as u64,
                prediction.orf_relative_stop as u64,
            );

            // INFO: getting NetTD record for the predicted ORF
            // INFO: format for output is: id start stop strand ns_start td_start td_stop strand
            if let Some(td_records) = transaid.get(&key) {
                found = true;
                for td_record in td_records {
                    let td_record = match td_record {
                        NetRecord::NetTD(record) => record,
                        _ => panic!("ERROR: unexpected record type"),
                    };

                    let line = format!(
                        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                        id,
                        orf_start,
                        orf_end,
                        prediction.orf_relative_start,
                        prediction.orf_relative_stop,
                        record
                            .strand()
                            .unwrap_or_else(|| panic!("ERROR: strand not found for {id:?}!")),
                        prediction.score,
                        td_record.start_score,
                        td_record.stop_score,
                        td_record.integrated_score
                    );

                    writer.write_all(line.as_bytes()).unwrap_or_else(|e| {
                        panic!("ERROR: failed to write record to file -> {e} -> {:?}", line);
                    });
                    writer.write_all(b"\n").unwrap_or_else(|e| {
                        panic!("ERROR: failed to write record to file -> {e} -> {:?}", line);
                    });
                }
            } else {
                let line = format!(
                    "{}\t{}\t{}\t{}\t{}\t{}\t{}\t-1\t-1\t-1",
                    id,
                    orf_start,
                    orf_end,
                    prediction.orf_relative_start,
                    prediction.orf_relative_stop,
                    record
                        .strand()
                        .unwrap_or_else(|| panic!("ERROR: strand not found for {id:?}!")),
                    prediction.score,
                );

                writer.write_all(line.as_bytes()).unwrap_or_else(|e| {
                    panic!("ERROR: failed to write record to file -> {e} -> {:?}", line);
                });
                writer.write_all(b"\n").unwrap_or_else(|e| {
                    panic!("ERROR: failed to write record to file -> {e} -> {:?}", line);
                });
            }

            if found {
                log::debug!("DEBUG: found match for {key:?} -> will remove from transaid");
                transaid.remove(&key);
            }
        }
    });

    if !transaid.is_empty() {
        log::warn!(
            "WARNING: transaid has {} entries remaining that were not matched",
            transaid.len()
        );
        for (key, records) in &transaid {
            log::warn!("WARNING: no match found for {key:?}");

            for record in records {
                let record = match record {
                    NetRecord::NetTD(record) => record,
                    _ => panic!("ERROR: unexpected record type"),
                };

                let gp = bed.get(&record.sequence_id).unwrap_or_else(|| {
                    panic!(
                        "ERROR: could not find genepred record for id: {}",
                        record.sequence_id
                    );
                });

                let (orf_start, orf_end) = map_absolute_cds(
                    gp,
                    record.orf_relative_start as u64,
                    record.orf_relative_stop as u64,
                );

                let line = format!(
                    "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                    record.sequence_id,
                    orf_start,
                    orf_end,
                    record.orf_relative_start,
                    record.orf_relative_stop,
                    gp.strand()
                        .unwrap_or_else(|| panic!("ERROR: strand not found for {gp:?}!")),
                    0.0,
                    record.start_score,
                    record.stop_score,
                    record.integrated_score
                );

                writer.write_all(line.as_bytes()).unwrap_or_else(|e| {
                    panic!("ERROR: failed to write record to file -> {e} -> {:?}", line);
                });
                writer.write_all(b"\n").unwrap_or_else(|e| {
                    panic!("ERROR: failed to write record to file -> {e} -> {:?}", line);
                });
            }
        }
    }
}

/// Map netstart or transaid file to a hashmap
/// with sequence_id as key and a vector of NetRecord
/// as value
///
/// # Arguments
///
///
///* `path` - Path to netstart or transaid file
///* `source` - NetSource enum
///
/// # Returns
///
/// * `HashMap<String, Vec<NetRecord>>` - Hashmap with sequence_id as key and a vector of NetRecord as value
///
/// # Panics
///
/// * If the netstart or transaid file cannot be read
/// * If the netstart or transaid file is invalid
///
/// # Examples
///
/// ```rust, ignore
/// use orf::cli::NetArgs;
/// use orf::nets::{run_nets, net_map};
///
/// let args = NetArgs {
///     bed: String::from("tests/data/test.bed"),
///     netstart: String::from("tests/data/test.netstart"),
///     transaid: String::from("tests/data/test.transaid"),
///     outdir: String::from("tests/data/test.outdir"),
/// };
///
/// run_nets(args);
///
/// let netstart = net_map(args.netstart, NetSource::Netstart);
/// let transaid = net_map(args.transaid, NetSource::Transaid);
/// ```
fn net_map<P: AsRef<Path> + std::fmt::Debug>(
    path: P,
    source: NetSource,
) -> HashMap<String, Vec<NetRecord>> {
    let contents =
        reader(&path).unwrap_or_else(|e| panic!("ERROR: failed to read net file {path:?} -> {e}"));
    let mut accumulator = HashMap::new();

    match source {
        NetSource::Netstart => contents.lines().for_each(|line| {
            if line.starts_with("origin") {
                return;
            }

            // INFO: because the output type a Hash + Vec struct is necessary
            // INFO: since the input is hyper chunked, performance is not affected
            // INFO: key is sequence_id, value is a vector of NetNS records
            let record = NetNS::parse(line);

            if let Ok(record) = record {
                accumulator
                    .entry(record.sequence_id.clone())
                    .or_insert_with(Vec::new)
                    .push(NetRecord::NetNS(record));
            } else {
                eprintln!(
                    "WARN: failed to parse netstart record, skipping -> {:?}",
                    line
                );
            }
        }),
        NetSource::Transaid => contents.lines().for_each(|line| {
            // WARN: this should change to a static pattern with
            // all entries from netstart2
            if line.starts_with("Sequence_ID") {
                return;
            }

            let record = NetTD::parse(line);

            // INFO: header is sequence_id#orf_relative_start#orf_relative_stop
            // for fast lookup from the NetNS side
            let header = format!(
                "{}#{}#{}",
                record.sequence_id, record.orf_relative_start, record.orf_relative_stop
            );

            accumulator
                .entry(header)
                .or_insert_with(Vec::new)
                .push(NetRecord::NetTD(record));
        }),
    }

    accumulator
}

/// NetSource enum
enum NetSource {
    Netstart,
    Transaid,
}

/// NetRecord enum
enum NetRecord {
    NetNS(NetNS),
    NetTD(NetTD),
}

/// Netstart holding struct
pub struct NetNS {
    // Chordata,44,209,55,ENSRNOT00000142028,+,0.145224
    pub orf_relative_start: usize,
    pub orf_relative_stop: usize,
    pub peptide_len: usize,
    pub sequence_id: String,
    pub strand: genepred::Strand,
    pub score: f32,
}

impl NetNS {
    /// Parse netstart line
    pub fn parse(line: &str) -> Result<Self, Box<dyn std::error::Error>> {
        let mut parts = line.split(',');

        // INFO: check if line has 7 parts dropping empty ones
        if parts.clone().filter(|x| !x.is_empty()).count() != 7 {
            return Err(Box::new(std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                format!("ERROR: invalid netstart line -> {:?}", line),
            )));
        }

        // WARN: Note that we are doing -1 to the start and +2 to the stop to match orfipy
        // on the blast side
        let (_, orf_relative_start, orf_relative_stop, peptide_len, sequence_id, strand, score) = (
            parts.next(),
            parts
                .next()
                .unwrap_or_else(|| {
                    panic!("ERROR: failed to parse start from line: {}", line);
                })
                .parse::<usize>()
                .unwrap_or_else(|e| {
                    panic!(
                        "ERROR: failed to parse ORF start from line: {} -> {e}",
                        line
                    );
                })
                - 1,
            parts
                .next()
                .unwrap_or_else(|| {
                    panic!("ERROR: failed to parse ORF end from line: {}", line);
                })
                .parse::<usize>()
                .unwrap_or_else(|e| {
                    panic!("ERROR: failed to parse ORF end from line: {} -> {e}", line);
                })
                + 2,
            parts
                .next()
                .unwrap_or_else(|| {
                    panic!("ERROR: failed to parse peptide length from line: {}", line);
                })
                .parse::<usize>()
                .unwrap_or_else(|e| {
                    panic!(
                        "ERROR: failed to parse peptide length from line: {} -> {e}",
                        line
                    );
                }),
            parts
                .next()
                .unwrap_or_else(|| {
                    panic!("ERROR: failed to sequence_id from line: {}", line);
                })
                .to_string(),
            parts
                .next()
                .unwrap_or_else(|| {
                    panic!("ERROR: failed to parse strand from line: {}", line);
                })
                .chars()
                .next()
                .unwrap_or_else(|| panic!("ERROR: failed to parse strand from line: {}", line)),
            parts
                .next()
                .unwrap_or_else(|| {
                    panic!("ERROR: failed to parse netstart score from line: {}", line);
                })
                .parse::<f32>()
                .unwrap_or_else(|e| {
                    panic!(
                        "ERROR: failed to parse tai netstart score from line: {} -> {e}",
                        line
                    );
                }),
        );

        let strand = match strand {
            '+' => genepred::Strand::Forward,
            '-' => genepred::Strand::Reverse,
            _ => panic!("ERROR: failed to parse strand from line: {}", line),
        };

        Ok(Self {
            orf_relative_start,
            orf_relative_stop,
            peptide_len,
            sequence_id,
            strand,
            score,
        })
    }
}

/// Transaid holding struct
#[allow(unused)]
pub struct NetTD {
    // XM_006253847.5,159,471,0.318,0.527,0.510,0.276,0.555,0.390,105,True,""
    sequence_id: String,
    orf_relative_start: usize, // WARN: needs +1 to match NetNS
    orf_relative_stop: usize,  // WARN: needs +1 to match NetNS
    start_score: f32,
    stop_score: f32,
    kozak_score: f32,
    cai_score: f32,
    gc_score: f32,
    integrated_score: f32,
    peptide_len: usize,
    passed_filter: bool,
}

impl NetTD {
    /// Parse transaid line
    pub fn parse(line: &str) -> Self {
        let mut parts = line.split(',');

        let (
            sequence_id,
            orf_relative_start,
            orf_relative_stop,
            start_score,
            stop_score,
            kozak_score,
            cai_score,
            gc_score,
            integrated_score,
            peptide_len,
            passed_filter,
            _,
        ) = (
            parts
                .next()
                .unwrap_or_else(|| {
                    panic!("ERROR: failed to parse sequence_id from line: {}", line);
                })
                .to_string(),
            parts
                .next()
                .unwrap_or_else(|| {
                    panic!("ERROR: failed to parse ORF start from line: {}", line);
                })
                .parse::<usize>()
                .unwrap_or_else(|e| {
                    panic!(
                        "ERROR: failed to parse ORF start from line: {} -> {e}",
                        line
                    );
                }),
            parts
                .next()
                .unwrap_or_else(|| {
                    panic!("ERROR: failed to parse ORF stop from line: {}", line);
                })
                .parse::<usize>()
                .unwrap_or_else(|e| {
                    panic!("ERROR: failed to parse ORF stop from line: {} -> {e}", line);
                })
                + 3,
            parts
                .next()
                .unwrap_or_else(|| {
                    panic!("ERROR: failed to parse start score from line: {}", line);
                })
                .parse::<f32>()
                .unwrap_or_else(|e| {
                    panic!(
                        "ERROR: failed to parse start score from line: {} -> {e}",
                        line
                    );
                }),
            parts
                .next()
                .unwrap_or_else(|| {
                    panic!("ERROR: failed to parse stop score from line: {}", line);
                })
                .parse::<f32>()
                .unwrap_or_else(|e| {
                    panic!(
                        "ERROR: failed to parse stop score from line: {} -> {e}",
                        line
                    );
                }),
            parts
                .next()
                .unwrap_or_else(|| {
                    panic!("ERROR: failed to parse kozak score from line: {}", line);
                })
                .parse::<f32>()
                .unwrap_or_else(|e| {
                    panic!(
                        "ERROR: failed to parse kozak score from line: {} -> {e}",
                        line
                    );
                }),
            parts
                .next()
                .unwrap_or_else(|| {
                    panic!("ERROR: failed to parse cai score from line: {}", line);
                })
                .parse::<f32>()
                .unwrap_or_else(|e| {
                    panic!(
                        "ERROR: failed to parse cai score from line: {} -> {e}",
                        line
                    );
                }),
            parts
                .next()
                .unwrap_or_else(|| {
                    panic!("ERROR: failed to parse gc score from line: {}", line);
                })
                .parse::<f32>()
                .unwrap_or_else(|e| {
                    panic!("ERROR: failed to parse gc score from line: {} -> {e}", line);
                }),
            parts
                .next()
                .unwrap_or_else(|| {
                    panic!(
                        "ERROR: failed to parse integrated score from line: {}",
                        line
                    );
                })
                .parse::<f32>()
                .unwrap_or_else(|e| {
                    panic!(
                        "ERROR: failed to parse integrated score from line: {} -> {e}",
                        line
                    );
                }),
            parts
                .next()
                .unwrap_or_else(|| {
                    panic!("ERROR: failed to parse peptide length from line: {}", line);
                })
                .parse::<usize>()
                .unwrap_or_else(|e| {
                    panic!(
                        "ERROR: failed to parse peptide length from line: {} -> {e}",
                        line
                    );
                }),
            parts
                .next()
                .unwrap_or_else(|| {
                    panic!("ERROR: failed to parse passed filter from line: {}", line);
                })
                .to_lowercase()
                .parse::<bool>()
                .unwrap_or_else(|e| {
                    panic!(
                        "ERROR: failed to parse passed filter from line: {} -> {e}",
                        line
                    );
                }),
            parts.next(),
        );

        Self {
            sequence_id,
            orf_relative_start,
            orf_relative_stop,
            start_score,
            stop_score,
            kozak_score,
            cai_score,
            gc_score,
            integrated_score,
            peptide_len,
            passed_filter,
        }
    }
}
