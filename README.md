# `rafft-rs`

## Prerequisites

- [`Rust`](https://rustup.rs/)
- [`ViennaRNA`](https://www.tbi.univie.ac.at/RNA/#download) including C header files

If `ViennaRNA` is installed in a non-default location, e.g. [`bioconda`](https://bioconda.github.io/user/install.html),
it might be necessary to set some [environment variables](https://crates.io/crates/librna-sys) before building:

```sh
export LIBRNA_INCLUDE_DIR=~/bin/miniconda3/envs/yourenv/include/ # adjust as necessary
export LIBRNA_LIB_DIR=~/bin/miniconda3/envs/yourenv/lib/ # adjust as necessary
export CPATH=$LIBRNA_INCLUDE_DIR
```

Alternatively, append `--features librna-sys/auto` to invocations of `cargo` if `ViennaRNA` is detectable via `pkg-config`.

## Building

_Note: Both CLI and bindings are not yet working!_

To build the main CLI executable `rufft`, run

```sh
cargo build --release
```

In addition, to build the dynamic library `librafft.so` containing `python` bindings, run

```sh
cargo build --release --features bindings
```

use `cargo doc --no-deps` to build the (incomplete) documentation.

## Usage

```sh
target/release/rufft -h
```

```python
# It is currently required to have librafft.so in the same directory as your python code
import librafft
```

## TODO

- [] license?
