use clap::{self, Parser};
use log::info;
use simple_logger::init_with_level;

use orf::{
    blast::run_blast,
    chunk::run_chunk,
    cli::{Args, Commands},
    nets::run_nets,
    samba::run_samba,
    tai::run_tai,
};

fn main() {
    let start = std::time::Instant::now();

    let args = Args::parse();

    init_with_level(args.level).unwrap();
    rayon::ThreadPoolBuilder::new()
        .num_threads(args.threads)
        .build()
        .unwrap();

    match args.command {
        Commands::Blast(blast_args) => run_blast(blast_args, args.threads),
        Commands::Tai(args) => run_tai(args),
        Commands::Samba(args) => run_samba(args),
        Commands::Chunk(args) => run_chunk(args),
        Commands::Net(args) => run_nets(args),
    }

    let elapsed = start.elapsed();
    info!("Elapsed time: {:.3?}", elapsed);
}
