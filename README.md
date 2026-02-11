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
Usage: seqsieve [OPTIONS] --outdir <OUTDIR> <PATH>

Arguments:
  <PATH>  Path to the file or directory to filter

Options:
  -o, --outdir <OUTDIR>    Output directory for filtered files
  -s, --stats <STATS>      Output file for statistics [default: stats.tsv]
  -p, --profile <PROFILE>  Profile [default: all profiles] [possible values: safe, interm, aggressive]
  -t, --threads <THREADS>  Number of threads [default: all]
  -h, --help               Print help
  -V, --version            Print version
```
