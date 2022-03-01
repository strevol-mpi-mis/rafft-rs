use clap::Parser;
use itertools::Itertools;
use std::io::Write;
use std::path::PathBuf;

use rafft::fast_folding::RafftConfig;
use rafft::folding_graph::RafftNodeInfo;
use rafft::{set_global_energy_parameters, set_global_temperature, VIENNA_VERSION};

#[derive(Parser, Debug)]
#[clap(
    name = "rufft",
    version,
    about = "RAFFT implemented in Rust. RNA structure and folding dynamics prediction using fast Fourier transform (https://github.com/strevol-mpi-mis/rafft-rs).",
    after_help(VIENNA_VERSION)
)]
struct Args {
    #[clap(
        parse(from_os_str),
        long = "params",
        short = 'P',
        help = "RNA secondary structure energy parameters."
    )]
    parameters: Option<PathBuf>,
    #[clap(help = "input RNA sequence")]
    sequence: String,
    #[clap(long = "AU", help = "Weight of AU base pairs", default_value = "2.0")]
    au: f64,
    #[clap(long = "GC", help = "Weight of GC base pairs", default_value = "3.0")]
    gc: f64,
    #[clap(long = "GU", help = "Weight of GU base pairs", default_value = "1.0")]
    gu: f64,
    #[clap(
        long = "temperature",
        short = 'T',
        help = "Temperature in °C, passed to ViennaRNA",
        default_value = "37.0"
    )]
    temperature: f64,
    #[clap(
        long = "min-unpaired",
        short = 'u',
        help = "Minimum amount of unpaired positions enclosed by a hairpin loop",
        default_value = "3"
    )]
    min_unpaired: usize,
    #[clap(
        long = "minimum-helix-energy",
        short = 'e',
        help = "Minimum energy [kcal/mol] candidate helix stacks have to contribute; lower (negative) is more stable",
        default_value = "0.0"
    )]
    min_loop_energy: f64,
    #[clap(
        long = "positional-lags",
        short = 'l',
        help = "Number of positional lags to search for possible stems",
        default_value = "100"
    )]
    positional_lags: usize,
    #[clap(
        long = "branch",
        short = 'b',
        help = "Number of branches to explore in construction of the fast-folding graph",
        default_value = "1000"
    )]
    number_of_branches: usize,
    #[clap(
        long = "saved-trajectories",
        short = 's',
        help = "Amount of saved structures per step of the breadth-first search",
        default_value = "1"
    )]
    saved_trajectories: usize,
    #[clap(
        long = "benchmark",
        short = 'B',
        help = "Format output suitable for internal benchmarks"
    )]
    benchmark: bool,
    #[clap(
        long = "compat",
        short = 'c',
        help = "Use an output format compatible to the kinetics scripts of RAFFT. This includes duplicate structures."
    )]
    compat: bool,
    #[clap(
        parse(from_os_str),
        long = "output-edges",
        short = 'o',
        help = "Write edges (pairs of structure indices) to the specified file. The indices correspond to the order of the printed structures."
    )]
    outfile: Option<PathBuf>,
}

fn main() {
    let args = Args::parse();

    if let Some(path) = args.parameters {
        set_global_energy_parameters(path);
    }

    #[allow(clippy::float_cmp)]
    if args.temperature != 37.0 {
        set_global_temperature(args.temperature);
    }

    let rafft_config = RafftConfig::new()
        .maximum_trajectories(args.saved_trajectories)
        .basepair_weights(args.au, args.gc, args.gu)
        .minimum_unpaired_in_hairpins(args.min_unpaired)
        .minimum_loop_energy(args.min_loop_energy)
        .maximum_branches(args.number_of_branches)
        .positional_lags(args.positional_lags);

    let mut ffgraph = rafft_config.folding_graph(&args.sequence);

    ffgraph.construct_trajectories();

    if !args.benchmark {
        if !args.compat {
            ffgraph.iter().for_each(|node| {
                println!(
                    "[{}] {} {:.2}",
                    node.depth,
                    node.structure.to_string(),
                    node.energy as f64 * 0.01
                );
            });
        } else {
            let mut grouped: Vec<(usize, Vec<&RafftNodeInfo>)> = vec![];
            for (depth, nodes) in &ffgraph.iter().group_by(|node| node.depth) {
                let mut nodes = nodes.collect::<Vec<_>>();

                if nodes.len() < args.saved_trajectories {
                    if let Some(previous) = grouped.last() {
                        let mut missing_previous =
                            previous.1[..args.saved_trajectories - nodes.len()].to_vec();

                        nodes.append(&mut missing_previous);
                        assert_eq!(nodes.len(), args.saved_trajectories);

                        nodes.sort_by_key(|node| node.energy);
                    }
                }
                grouped.push((depth, nodes));
            }

            for (depth, nodes) in grouped {
                println!("# ---------{}----------", depth);
                nodes.iter().for_each(|node| {
                    println!(
                        "{} {:.2}",
                        node.structure.to_string(),
                        node.energy as f64 * 0.01
                    );
                });
            }
        }

        if let Some(outfile) = args.outfile {
            if let Ok(mut file) = std::fs::File::create(&outfile) {
                ffgraph.adjacent_indices().for_each(|(i, j)| {
                    writeln!(file, "{} {}", i, j).unwrap();
                })
            }
        }
    } else {
        let mut trajectories: Vec<_> = ffgraph.iter().collect();

        trajectories.sort_by_key(|node| node.energy);

        for node in &trajectories[..args.saved_trajectories.min(trajectories.len())] {
            println!(
                "{} {} {} {:.1} {}",
                args.sequence,
                args.sequence.len(),
                node.structure.to_string(),
                node.energy as f64 * 0.01,
                node.structure.pairs()
            );
        }
    }
}
