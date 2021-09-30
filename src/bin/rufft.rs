use std::path::PathBuf;
use structopt::StructOpt;

use rafft::*;

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
    #[structopt(long = "AU", help = "Weight of AU base pairs", default_value = "1.0")]
    au: f64,
    #[structopt(long = "GC", help = "Weight of GC base pairs", default_value = "1.0")]
    gc: f64,
    #[structopt(long = "GU", help = "Weight of GU base pairs", default_value = "1.0")]
    gu: f64,
}

fn main() {
    let opt = Opt::from_args();
}
