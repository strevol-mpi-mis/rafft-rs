use std::path::PathBuf;
use structopt::StructOpt;

use rafft::fast_folding::RafftConfig;
use rafft::{set_global_energy_parameters, set_global_temperature};

#[derive(StructOpt, Debug)]
#[structopt(
    name = "rufft",
    about = "RAFFT implemented in Rust. RNA structure and folding dynamics prediction using fast Fourier transform."
)]
struct Opt {
    #[structopt(
        parse(from_os_str),
        long = "params",
        short = "P",
        help = "RNA secondary structure energy parameters."
    )]
    parameters: Option<PathBuf>,
    #[structopt(help = "input RNA sequence")]
    sequence: String,
    #[structopt(long = "AU", help = "Weight of AU base pairs", default_value = "2.0")]
    au: f64,
    #[structopt(long = "GC", help = "Weight of GC base pairs", default_value = "3.0")]
    gc: f64,
    #[structopt(long = "GU", help = "Weight of GU base pairs", default_value = "1.0")]
    gu: f64,
    #[structopt(
        long = "temperature",
        short = "T",
        help = "Temperature in °C, passed to ViennaRNA",
        default_value = "37.0"
    )]
    temperature: f64,
    #[structopt(
        long = "min-unpaired",
        short = "u",
        help = "Minimum amount of unpaired positions enclosed by a hairpin loop",
        default_value = "3"
    )]
    min_unpaired: usize,
    #[structopt(
        long = "minimum-helix-energy",
        short = "e",
        help = "Minimum energy candidate helix stacks have to contribute; lower (negative) is more stable",
        default_value = "0.0"
    )]
    min_loop_energy: f64,
    #[structopt(
        long = "positional-lags",
        short = "l",
        help = "Minimum amount of unpaired positions enclosed by a hairpin loop",
        default_value = "100"
    )]
    positional_lags: usize,
    #[structopt(
        long = "branch",
        short = "b",
        help = "Number of branches to explore in construction of the fast-folding graph",
        default_value = "1000"
    )]
    number_of_branches: usize,
    #[structopt(
        long = "saved-trajectories",
        short = "s",
        help = "Amount of saved structures per step of the breadth-first search",
        default_value = "1"
    )]
    saved_trajectories: usize,
}

fn main() {
    let opt = Opt::from_args();

    if let Some(path) = opt.parameters {
        set_global_energy_parameters(path);
    }

    #[allow(clippy::float_cmp)]
    if opt.temperature != 37.0 {
        set_global_temperature(opt.temperature);
    }

    let rafft_config = RafftConfig::new()
        .maximum_trajectories(opt.saved_trajectories)
        .basepair_weights(opt.au, opt.gc, opt.gu)
        .minimum_unpaired_in_hairpins(opt.min_unpaired)
        .minimum_loop_energy(opt.min_loop_energy)
        .maximum_branches(opt.number_of_branches)
        .positional_lags(opt.positional_lags);

    let mut ffgraph = rafft_config.folding_graph(&opt.sequence);

    ffgraph.construct_trajectories();

    ffgraph.iter().for_each(|node| {
        println!(
            "[{}] {} {:.2}",
            node.depth,
            node.structure.to_string(),
            node.energy as f64 * 0.01
        );
    });
}
