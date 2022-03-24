# `rafft-rs`

Rust implementation of [RAFFT](https://github.com/strevol-mpi-mis/RAFFT), an algorithm for [*efficient prediction of RNA folding pathways using the fast Fourier transform*](https://doi.org/10.1101/2021.07.02.450908).

## Prerequisites

- [`Rust`](https://www.rust-lang.org/tools/install)
- [`ViennaRNA`](https://www.tbi.univie.ac.at/RNA/#download) including C header files

If `ViennaRNA` is installed in a non-default location, e.g. by using [`bioconda`](https://bioconda.github.io/user/install.html),
it might be necessary to set some [environment variables](https://crates.io/crates/librna-sys) before building:

```sh
export LIBRNA_INCLUDE_DIR=~/bin/miniconda3/envs/yourenv/include/ # adjust as necessary
export LIBRNA_LIB_DIR=~/bin/miniconda3/envs/yourenv/lib/ # adjust as necessary
export CPATH=$LIBRNA_INCLUDE_DIR
```

Alternatively, append `--features librna-sys/auto` to invocations of `cargo` if `pkg-config --exists RNAlib2` returns successfully.

### Example Setup

If in doubt, please refer to [`Setup.md`](https://github.com/strevol-mpi-mis/rafft-rs/blob/main/Setup.md) for a reproducible *step-by-step* guide setting up 
a working environment for building `rafft-rs`.

## Building

To build the main CLI executable `rufft`, run

```sh
cargo build --bin rufft --release
```

In addition, to build the dynamic library `librafft.so` containing `python` bindings, run

```sh
cargo build --release --features bindings
```

Use `cargo doc --no-deps` to build the API documentation.

## Usage

### CLI

```sh
target/release/rufft -h
```
```sh
rufft X.Y.Z
RAFFT implemented in Rust. RNA structure and folding dynamics prediction using fast Fourier
transform (https://github.com/strevol-mpi-mis/rafft-rs).

USAGE:
    rufft [OPTIONS] <SEQUENCE>

ARGS:
    <SEQUENCE>    input RNA sequence

OPTIONS:
        --AU <AU>
            Weight of AU base pairs [default: 2.0]

        --GC <GC>
            Weight of GC base pairs [default: 3.0]

        --GU <GU>
            Weight of GU base pairs [default: 1.0]

    -b, --branch <NUMBER_OF_BRANCHES>
            Number of branches to explore in construction of the fast-folding graph [default: 1000]

    -B, --benchmark
            Format output suitable for internal benchmarks

    -c, --compat
            Use an output format compatible to the kinetics scripts of RAFFT. This includes
            duplicate structures.

    -e, --minimum-helix-energy <MIN_LOOP_ENERGY>
            Minimum energy [kcal/mol] candidate helix stacks have to contribute; lower (negative) is
            more stable [default: 0.0]

    -h, --help
            Print help information

    -l, --positional-lags <POSITIONAL_LAGS>
            Number of positional lags to search for possible stems [default: 100]

    -o, --output-edges <OUTFILE>
            Write edges (pairs of structure indices) to the specified file. The indices correspond
            to the order of the printed structures.

    -P, --params <PARAMETERS>
            RNA secondary structure energy parameters.

    -s, --saved-trajectories <SAVED_TRAJECTORIES>
            Amount of saved structures per step of the breadth-first search [default: 1]

    -T, --temperature <TEMPERATURE>
            Temperature in °C, passed to ViennaRNA [default: 37.0]

    -u, --min-unpaired <MIN_UNPAIRED>
            Minimum amount of unpaired positions enclosed by a hairpin loop [default: 3]

    -V, --version
            Print version information

U.V.W
```

#### Example

```sh
target/release/rufft GGGUUUGCGGUGUAAGUGCAGCCCGUCUUACACCGUGCGGCACAGGCACUAGUACUGAUGUCGUAUACAGGGCUUUUGACAU --saved-trajectories 5 --compat
```


### Python Bindings

```python
# It is currently required to have librafft.so in the same directory as your python code
# Alternatively, make sure that librafft.so is in your python path:
# import sys
# sys.path.append("./target/release/")
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

# Afterwards, obtain Metropolis transition rates:
from scipy.sparse import coo_matrix

# with optional parameter beta = 0.61
rates, i_s, j_s = ffgraph.transition_rates()
trmatrix = coo_matrix((rates, (i_s, j_s))).toarray() # or .tocsr()
# see also
#ffgraph.directed_edges()
```

Note that this implementation does not store duplicate structures in the fast-folding graph.
This saves an intermediate processing step where duplicate structures need to be removed 
before any kinetic analysis is pursued.

This does not affect the (unique) structures being discovered.
However, this might change in a future release as it would allow to discover a few more distinct structures
at no asymptotical complexity increase.

A commandline flag `--compat`/`-c` is available to display redundant structures of the fast-folding graph.
Use this flag to generate output that can be used with the existing scripts of 
[RAFFT](https://github.com/strevol-mpi-mis/RAFFT) to produce kinetics.
