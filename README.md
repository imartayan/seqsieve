# SeqSieve

This is a small tool that filters all the (possibly compressed) FASTA files in a given directory by discarding sequences based on multiple criteria:
- very short sequences
- invalid nucleotides
- too many Ns
- large homopolymers
- small 1-mer entropy
- small 2-mer entropy

Most of the computation is vectorized using [helicase](https://github.com/imartayan/helicase) and files are processed in parallel.

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
