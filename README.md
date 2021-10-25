# `rafft-rs`

Rust implementation of [RAFFT](https://github.com/strevol-mpi-mis/RAFFT), an algorithm for [*efficient prediction of RNA folding pathways using the fast Fourier transform*](https://doi.org/10.1101/2021.07.02.450908).

## Prerequisites

- [`Rust`](https://rustup.rs/)
- [`ViennaRNA`](https://www.tbi.univie.ac.at/RNA/#download) including C header files

If `ViennaRNA` is installed in a non-default location, e.g. by using [`bioconda`](https://bioconda.github.io/user/install.html),
it might be necessary to set some [environment variables](https://crates.io/crates/librna-sys) before building:

```sh
export LIBRNA_INCLUDE_DIR=~/bin/miniconda3/envs/yourenv/include/ # adjust as necessary
export LIBRNA_LIB_DIR=~/bin/miniconda3/envs/yourenv/lib/ # adjust as necessary
export CPATH=$LIBRNA_INCLUDE_DIR
```

Alternatively, append `--features librna-sys/auto` to invocations of `cargo` if `pkg-config --exists RNAlib2` returns successfully.

## Building

To build the main CLI executable `rufft`, run

```sh
cargo build --bin rufft --release
```

In addition, to build the dynamic library `librafft.so` containing `python` bindings, run

```sh
cargo build --release --features bindings
```

use `cargo doc --no-deps` to build the documentation.

## Usage

```sh
target/release/rufft -h
```

```python
# It is currently required to have librafft.so in the same directory as your python code
from librafft import FastFoldingGraph, set_temperature, set_energy_parameters

# with optional (named or ordered) parameters:
# number_of_lags = 100,
# number_of_branches = 1000,
# saved_trajectories = 1,
# au = 2.0,
# gc = 3.0,
# gu = 1.0,
# min_unpaired = 3,
# min_loop_energy = 0.0
ffgraph = FastFoldingGraph("GGGUUUGCGGUGUAAGUGCAGCCCGUCUUACACCGUGCGGCACAGGCACUAGUACUGAUGUCGUAUACAGGGCUUUUGACAU")
print(ffgraph.trajectories())
```
