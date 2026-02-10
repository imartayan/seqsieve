# SeqSieve

## Build instructions

```sh
git clone https://github.com/imartayan/seqsieve.git
cd seqsieve
RUSTFLAGS="-C target-cpu=native" cargo build --release
cp target/release/seqsieve .
```

## Usage

```
Usage: seqsieve [OPTIONS] --outdir <OUTDIR> --mode <MODE> <PATH>

Arguments:
  <PATH>  Path to the file or directory to filter

Options:
  -o, --outdir <OUTDIR>    Output directory for filtered files
  -m, --mode <MODE>        Mode [possible values: safe, interm, aggressive]
  -t, --threads <THREADS>  Number of threads [default: all]
  -h, --help               Print help
  -V, --version            Print version
```
