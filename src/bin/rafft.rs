use std::path::PathBuf;
use structopt::StructOpt;

use rafft::*;

#[derive(StructOpt, Debug)]
#[structopt(
    name = "dmc",
    about = "Diffusion Monte Carlo method to calculate eigenvector centrality in a neutral network."
)]
struct Opt {
    #[structopt(parse(from_os_str), help = "inputfile")]
    filename: PathBuf,
    #[structopt(help = "fraction of nodes to use as population", default_value = "1.0")]
    alpha: f64,
    #[structopt(help = "number of iterations to run", default_value = "100")]
    t: usize,
    #[structopt(
        long = "sync",
        short = "s",
        help = "synchronize every sub_t steps",
        default_value = "1"
    )]
    sub_t: usize,
    #[structopt(help = "number of parallel jobs/tasks", default_value = "4")]
    j: usize,
    #[structopt(
        long = "exact",
        short = "x",
        help = "Use exact method (LOBPCG) to determine largest eigenvalue."
    )]
    exact: bool,
}

fn main() {
    let opt = Opt::from_args();
}
