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

* To see all available commands run `cargo run --release -- help`
* Run the example script `test-run-split.sh`

### Docker

* Pull the image `docker pull z1ko/crisprme:cuda12.6.3-ubuntu22.04`
* Run the interactive container `docker run -it z1ko/crisprme:cuda12.6.3-ubuntu22.04`

## License

MIT license([LICENSE](LICENSE))
