use clap::{Parser, Subcommand, ValueEnum};
use std::path::PathBuf;

pub const NMD_DISTANCE: u64 = 55; // 55 bp
pub const WEAK_NMD_DISTANCE: i64 = 80; // 85 bp
pub const ATG_DISTANCE: u64 = 100; // 100 bp
pub const BIG_EXON_DIST_TO_EJ: u64 = 400; // 400 bp

#[derive(Debug, Parser)]
#[command(name = "orf", about = "Open reading frame pipeline wrapper", version = env!("CARGO_PKG_VERSION"))]
pub struct Args {
    #[command(subcommand)]
    pub command: Commands,

    #[arg(
        short = 't',
        long = "threads",
        help = "Number of threads",
        value_name = "THREADS",
        default_value_t = num_cpus::get()
    )]
    pub threads: usize,

    #[arg(
        short = 'L',
        long = "level",
        help = "Logging level",
        value_name = "LEVEL",
        default_value_t = log::Level::Info,
    )]
    pub level: log::Level,
}

#[derive(Debug, Subcommand)]
pub enum Commands {
    /// Run ORF detection using orfipy + diamond + psauron
    Blast(BlastArgs),

    /// Run ORF detection using TranslationAi
    Tai(TaiArgs),

    /// Chunk genomic regions and sequences
    Chunk(ChunkArgs),

    /// Run RNAsamba on a set of reads
    Samba(SambaArgs),

    /// Join nets (NetStart2 [optional] + Transaid)
    Net(NetArgs),
}

#[derive(Debug, Parser)]
pub struct BlastArgs {
    #[arg(
        short = 'f',
        long = "fasta",
        required = true,
        help = "Path to .fa file"
    )]
    pub fasta: PathBuf,

    #[arg(
        short = 'b',
        long = "bed",
        required = true,
        help = "Path to .bed file with candidate regions"
    )]
    pub bed: PathBuf,

    #[arg(
        short = 'o',
        long = "outdir",
        required = false,
        help = "Path to outdir",
        default_value = "."
    )]
    pub outdir: PathBuf,

    #[arg(
        short = 'D',
        long = "database",
        required = true,
        help = "Path to protein database"
    )]
    pub database: PathBuf,

    #[arg(
        short = 't',
        long = "tai",
        required = false,
        help = "Path to translationAi output file"
    )]
    pub tai: Option<PathBuf>,

    #[arg(
        short = 'N',
        long = "net",
        required = false,
        help = "Path to net output file"
    )]
    pub net: Option<PathBuf>,

    #[arg(
        short = 'l',
        long = "orf-min-len",
        required = false,
        help = "Minimum length for ORFs",
        default_value = "50"
    )]
    pub orf_min_len: usize,

    #[arg(
        short = 'p',
        long = "orf-min-percent",
        required = false,
        help = "Minimum percentage of ORF length that is aligned",
        default_value = "0.25"
    )]
    pub orf_min_percent: f32,

    #[arg(
        short = 'u',
        long = "upstream-flank",
        required = false,
        help = "Number of bases upstream of the ORF to include in the context sequence",
        default_value = "0"
    )]
    pub upstream_flank: usize,

    #[arg(
        short = 'd',
        long = "downstream-flank",
        required = false,
        help = "Number of bases downstream of the ORF to include in the context sequence",
        default_value = "0"
    )]
    pub downstream_flank: usize,

    #[arg(
        short = 'n',
        long = "nmd-distance",
        help = "Distance to consider NMD",
        value_name = "DISTANCE",
        default_value_t = NMD_DISTANCE,
    )]
    pub nmd_distance: u64,

    #[arg(
        short = 'w',
        long = "weak-nmd-distance",
        help = "Distance to consider weak NMD",
        value_name = "DISTANCE",
        default_value_t = WEAK_NMD_DISTANCE,
    )]
    pub weak_nmd_distance: i64,

    #[arg(
        short = 'a',
        long = "atg-distance",
        help = "Distance to consider ATG",
        value_name = "DISTANCE",
        default_value_t = ATG_DISTANCE,
    )]
    pub atg_distance: u64,

    #[arg(
        short = 'e',
        long = "big-exon-dist-to-ej",
        help = "Distance to consider big exon to exon junction",
        value_name = "DISTANCE",
        default_value_t = BIG_EXON_DIST_TO_EJ,
    )]
    pub big_exon_dist_to_ej: u64,

    #[arg(
        short = 'P',
        long = "prefix",
        required = false,
        help = "Prefix for output files",
        default_value = "xorf"
    )]
    pub prefix: String,

    #[arg(
        short = 'K',
        long = "keep-temp",
        help = "Keep temporary files",
        action = clap::ArgAction::SetTrue
    )]
    pub keep_temp: bool,

    #[arg(
        long = "engine",
        required = false,
        help = "Engine to use for blast [diamond, mmseqs2]",
        default_value_t = Engine::Diamond,
        value_enum
    )]
    pub engine: Engine,
}

#[derive(Debug, Parser)]
pub struct TaiArgs {
    #[arg(
        short = 'f',
        long = "fasta",
        required = true,
        help = "Path to .fa/.fa.gz"
    )]
    pub fasta: PathBuf,

    #[arg(
        short = 'b',
        long = "bed",
        required = true,
        help = "Path to .bed file with candidate regions"
    )]
    pub bed: PathBuf,

    #[arg(
        short = 'o',
        long = "outdir",
        required = false,
        help = "Path to outdir",
        default_value = "."
    )]
    pub outdir: PathBuf,

    #[arg(
        short = 'u',
        long = "upstream-flank",
        required = false,
        help = "Number of bases upstream of the TSS to include in the chunk",
        default_value = "0"
    )]
    pub upstream_flank: usize,

    #[arg(
        short = 'd',
        long = "downstream-flank",
        required = false,
        help = "Number of bases downstream of the TSS to include in the chunk",
        default_value = "0"
    )]
    pub downstream_flank: usize,
}

#[derive(Debug, Parser)]
pub struct ChunkArgs {
    #[arg(
        short = 's',
        long = "sequence",
        required = true,
        help = "Path to .fa/.fa.gz/.2bit"
    )]
    pub sequence: PathBuf,

    #[arg(
        short = 'r',
        long = "regions",
        required = true,
        help = "Path to .bed file with candidate regions"
    )]
    pub regions: PathBuf,

    #[arg(
        short = 'o',
        long = "outdir",
        required = false,
        help = "Path to outdir",
        default_value = "."
    )]
    pub outdir: PathBuf,

    #[arg(
        short = 'c',
        long = "chunks",
        required = false,
        help = "Number of chunks to split the bed file into",
        default_value = "10000"
    )]
    pub chunks: usize,

    #[arg(
        short = 'u',
        long = "upstream-flank",
        required = false,
        help = "Number of bases upstream of the TSS to include in the chunk",
        default_value = "0"
    )]
    pub upstream_flank: usize,

    #[arg(
        short = 'd',
        long = "downstream-flank",
        required = false,
        help = "Number of bases downstream of the TSS to include in the chunk",
        default_value = "0"
    )]
    pub downstream_flank: usize,

    #[arg(
        short = 'I',
        long = "ignore-errors",
        help = "Ignore slicing errors coming from out-of-bounds coordinates",
        action = clap::ArgAction::SetTrue
    )]
    pub ignore_errors: bool,

    #[arg(
        short = 'p',
        long = "prefix",
        required = false,
        help = "Prefix for output files"
    )]
    pub prefix: Option<String>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, ValueEnum)]
pub enum Engine {
    Diamond,
    Mmseqs2,
}

impl std::fmt::Display for Engine {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Engine::Diamond => write!(f, "diamond"),
            Engine::Mmseqs2 => write!(f, "mmseqs2"),
        }
    }
}

#[derive(Debug, Parser)]
pub struct SambaArgs {
    #[arg(
        short = 'f',
        long = "fasta",
        required = true,
        help = "Path to .fa file"
    )]
    pub fasta: PathBuf,

    #[arg(
        short = 'o',
        long = "outdir",
        default_value = ".",
        help = "Output directory"
    )]
    pub outdir: PathBuf,

    #[arg(
        short = 'w',
        long = "weights",
        required = false,
        help = "Path to model weights file"
    )]
    pub weights: Option<PathBuf>,

    #[arg(
        short = 'u',
        long = "upstream-flank",
        required = false,
        help = "Number of bases upstream of the TSS to include in the chunk",
        default_value = "0"
    )]
    pub upstream_flank: usize,

    #[arg(
        short = 'd',
        long = "downstream-flank",
        required = false,
        help = "Number of bases downstream of the TSS to include in the chunk",
        default_value = "0"
    )]
    pub downstream_flank: usize,
}

#[derive(Debug, Parser)]
pub struct NetArgs {
    #[arg(
        short = 'b',
        long = "bed",
        required = true,
        help = "Path to .bed file with candidate regions"
    )]
    pub bed: PathBuf,

    #[arg(
        short = 'n',
        long = "netstart",
        required = false,
        help = "Path to netstart .csv file"
    )]
    pub netstart: Option<PathBuf>,

    #[arg(
        short = 't',
        long = "transaid",
        required = true,
        help = "Path to transaid .csv file"
    )]
    pub transaid: PathBuf,

    #[arg(
        short = 'o',
        long = "outdir",
        required = false,
        help = "Path to outdir",
        default_value = "."
    )]
    pub outdir: PathBuf,
}
