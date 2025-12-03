# crisprme2

[![Crates.io](https://img.shields.io/crates/v/crisprme2.svg)](https://crates.io/crates/crisprme2)
[![Docs.rs](https://docs.rs/crisprme2/badge.svg)](https://docs.rs/crisprme2)
[![CI](https://github.com/z1ko/crisprme2/workflows/CI/badge.svg)](https://github.com/z1ko/crisprme2/actions)

## Installation

Tested with CUDA 12.6, should work with all versions, prefer CUDA 11+

### Cargo

* Install the rust toolchain in order to have cargo installed by following
  [this](https://www.rust-lang.org/tools/install) guide.

### Run

* Build the engine using: `cargo b --release`
* To see all available commands run: `./target/release/crisprme help`
* Example mine:
  ```
  ./target/release/crisprme mine \
    --tgap <target gaps>         \
    --qgap <query gaps>          \
    --mism <mismatches>          \
    <dataset>                    \
    <target sequence length>     \
    <query string>               \
    <output>
  ```

### Preprocessing

The datasets must be firstly preprocess with the following commands:
* If working with a FASTA file: `./target/release/crisprme preprocess <input file> <target sequence length> <delta between windows>`
* If working with a list of ids and sequences: `./target/release/crisprme preprocess-list <input file> <target sequence length>`
