use clap::{Parser as ClapParser, ValueEnum};
use helicase::config::advanced::COMPUTE_DNA_STRING;
use helicase::dna_format::ColumnarDNA;
use helicase::input::*;
use helicase::*;
use jwalk::WalkDir;
use rayon::prelude::*;

use std::io;
use std::path::{Path, PathBuf};

#[derive(ClapParser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
    /// Path to the file or directory to filter
    path: String,
    /// Output directory for filtered files
    #[arg(short, long)]
    outdir: String,
    /// Mode
    #[arg(short, long, value_enum)]
    mode: Mode,
    /// Number of threads [default: all]
    #[arg(short, long)]
    threads: Option<usize>,
}

#[derive(ValueEnum, Clone, Copy, Debug)]
enum Mode {
    Safe,
    Interm,
    Aggressive,
}

#[derive(Debug)]
struct Params {
    pub min_len: usize,
    pub max_n_frac: f64,
    pub max_homopoly: usize,
    pub min_entropy_1: f64,
    pub min_entropy_2: f64,
    pub strict_nuc: bool,
}

impl Params {
    #[inline(always)]
    const fn new(
        min_len: usize,
        max_n_frac: f64,
        max_homopoly: usize,
        min_entropy_1: f64,
        min_entropy_2: f64,
        strict_nuc: bool,
    ) -> Self {
        Self {
            min_len,
            max_n_frac,
            max_homopoly,
            min_entropy_1,
            min_entropy_2,
            strict_nuc,
        }
    }

    #[inline(always)]
    const fn from_mode(mode: Mode) -> Self {
        match mode {
            Mode::Safe => Self::new(80, 0.05, 60, 0.8, 1.5, true),
            Mode::Interm => Self::new(60, 0.20, 120, 0.2, 0.4, false),
            Mode::Aggressive => Self::new(51, 0.50, 200, 0.0, 0.0, false),
        }
    }
}

#[derive(Debug)]
enum Status {
    Ok,
    TooShort,
    TooManyNs,
    BigHomopoly,
    SmallEntropy1,
    SmallEntropy2,
    InvalidNuc,
}

fn fasta_files<P: AsRef<Path>>(path: P) -> io::Result<Vec<PathBuf>> {
    let path = path.as_ref();
    if !path.exists() {
        return Err(io::Error::other("Path does not exist"));
    }
    if path.is_dir() {
        Ok(fasta_files_in_dir(path))
    } else {
        Ok(vec![path.to_path_buf()])
    }
}

fn fasta_files_in_dir<P: AsRef<Path>>(dir: P) -> Vec<PathBuf> {
    WalkDir::new(dir)
        .min_depth(1)
        .into_iter()
        .filter_map(|entry| {
            let entry = entry.ok()?;
            let name = entry.file_name().to_str()?;
            if is_fasta_file(name) {
                Some(entry.path())
            } else {
                None
            }
        })
        .collect()
}

#[inline(always)]
fn is_fasta_file(name: &str) -> bool {
    name.ends_with(".fa")
        || name.ends_with(".fna")
        || name.ends_with(".fasta")
        || name.ends_with(".fa.gz")
        || name.ends_with(".fna.gz")
        || name.ends_with(".fasta.gz")
        || name.ends_with(".fa.zst")
        || name.ends_with(".fna.zst")
        || name.ends_with(".fasta.zst")
}

fn filter_file(filename: PathBuf, mode: Mode) {
    const CONFIG: Config = ParserOptions::default()
        .dna_columnar()
        .keep_non_actg()
        .config()
        | COMPUTE_DNA_STRING;

    dbg!(&filename);
    // TODO check if output file exists
    let mut parser = FastaParser::<CONFIG, _>::from_file(filename).expect("Failed to parse file");
    let params = Params::from_mode(mode);

    while let Some(_) = parser.next() {
        let seq = parser.get_dna_columnar();
        let status = filter_seq(seq, &params);
        match status {
            Status::Ok => {
                let _ = parser.get_header();
                let _ = parser.get_dna_string();
                // TODO write tmp output file
            }
            Status::TooShort => (),
            Status::TooManyNs => (),
            Status::BigHomopoly => (),
            Status::SmallEntropy1 => (),
            Status::SmallEntropy2 => (),
            Status::InvalidNuc => (),
        }
    }
    // TODO move tmp output file
    // TODO return stats
}

#[inline(always)]
fn filter_seq(seq: &ColumnarDNA, params: &Params) -> Status {
    if seq.len() < params.min_len {
        return Status::TooShort;
    }
    let rem = seq.len() % 64;

    let mut count1 = [0; 4];
    let mut count2 = [0; 16];
    let mut max_hpoly = 1;
    let mut cur_hpoly = 0;

    let (high_bits, hi) = seq.high_bits();
    let (low_bits, lo) = seq.low_bits();
    let (carry_hi, carry_lo) = seq.get(0);
    let mut carry_hi = !carry_hi as u64;
    let mut carry_lo = !carry_lo as u64;

    for (&hi, &lo) in high_bits.iter().zip(low_bits) {
        // TODO check N mask

        let prev_hi = (hi << 1) | carry_hi;
        let prev_lo = (lo << 1) | carry_lo;
        let va = !hi & !lo;
        let vc = !hi & lo;
        let vt = hi & !lo;
        let vg = hi & lo;

        if params.min_entropy_1 > 0.0 {
            count1[0] += va.count_ones() as usize;
            count1[1] += vc.count_ones() as usize;
            count1[2] += vt.count_ones() as usize;
            count1[3] += vg.count_ones() as usize;
        }

        if params.min_entropy_2 > 0.0 {
            let prev_a = !prev_hi & !prev_lo;
            let prev_c = !prev_hi & prev_lo;
            let prev_t = prev_hi & !prev_lo;
            let prev_g = prev_hi & prev_lo;
            count2[0] += (prev_a & va).count_ones() as usize;
            count2[1] += (prev_a & vc).count_ones() as usize;
            count2[2] += (prev_a & vt).count_ones() as usize;
            count2[3] += (prev_a & vg).count_ones() as usize;
            count2[4] += (prev_c & va).count_ones() as usize;
            count2[5] += (prev_c & vc).count_ones() as usize;
            count2[6] += (prev_c & vt).count_ones() as usize;
            count2[7] += (prev_c & vg).count_ones() as usize;
            count2[8] += (prev_t & va).count_ones() as usize;
            count2[9] += (prev_t & vc).count_ones() as usize;
            count2[10] += (prev_t & vt).count_ones() as usize;
            count2[11] += (prev_t & vg).count_ones() as usize;
            count2[12] += (prev_g & va).count_ones() as usize;
            count2[13] += (prev_g & vc).count_ones() as usize;
            count2[14] += (prev_g & vt).count_ones() as usize;
            count2[15] += (prev_g & vg).count_ones() as usize;
        }

        if params.max_homopoly > 0 {
            let eq_hi = !(prev_hi ^ hi);
            let eq_lo = !(prev_lo ^ lo);
            let mut eq = eq_hi & eq_lo;
            if eq == !0 {
                cur_hpoly += 64;
            } else {
                let t = eq.trailing_ones();
                cur_hpoly += t;
                eq >>= t;
                max_hpoly = max_hpoly.max(cur_hpoly);
                while eq != 0 {
                    eq >>= eq.trailing_zeros();
                    let t = eq.trailing_ones();
                    cur_hpoly = t;
                    eq >>= t;
                    max_hpoly = max_hpoly.max(cur_hpoly);
                }
            }
        }

        carry_hi = hi >> 63;
        carry_lo = lo >> 63;
    }

    if rem > 0 {
        let mask = !0 >> (64 - rem);
        // TODO check N mask

        let prev_hi = (hi << 1) | carry_hi;
        let prev_lo = (lo << 1) | carry_lo;
        let va = (!hi & !lo) & mask;
        let vc = (!hi & lo) & mask;
        let vt = (hi & !lo) & mask;
        let vg = (hi & lo) & mask;

        if params.min_entropy_1 > 0.0 {
            count1[0] += va.count_ones() as usize;
            count1[1] += vc.count_ones() as usize;
            count1[2] += vt.count_ones() as usize;
            count1[3] += vg.count_ones() as usize;
        }

        if params.min_entropy_2 > 0.0 {
            let prev_a = (!prev_hi & !prev_lo) & mask;
            let prev_c = (!prev_hi & prev_lo) & mask;
            let prev_t = (prev_hi & !prev_lo) & mask;
            let prev_g = (prev_hi & prev_lo) & mask;
            count2[0] += (prev_a & va).count_ones() as usize;
            count2[1] += (prev_a & vc).count_ones() as usize;
            count2[2] += (prev_a & vt).count_ones() as usize;
            count2[3] += (prev_a & vg).count_ones() as usize;
            count2[4] += (prev_c & va).count_ones() as usize;
            count2[5] += (prev_c & vc).count_ones() as usize;
            count2[6] += (prev_c & vt).count_ones() as usize;
            count2[7] += (prev_c & vg).count_ones() as usize;
            count2[8] += (prev_t & va).count_ones() as usize;
            count2[9] += (prev_t & vc).count_ones() as usize;
            count2[10] += (prev_t & vt).count_ones() as usize;
            count2[11] += (prev_t & vg).count_ones() as usize;
            count2[12] += (prev_g & va).count_ones() as usize;
            count2[13] += (prev_g & vc).count_ones() as usize;
            count2[14] += (prev_g & vt).count_ones() as usize;
            count2[15] += (prev_g & vg).count_ones() as usize;
        }

        if params.max_homopoly > 0 {
            let eq_hi = !(prev_hi ^ hi);
            let eq_lo = !(prev_lo ^ lo);
            let mut eq = (eq_hi & eq_lo) & mask;

            let t = eq.trailing_ones();
            cur_hpoly += t;
            eq >>= t;
            max_hpoly = max_hpoly.max(cur_hpoly);
            while eq != 0 {
                eq >>= eq.trailing_zeros();
                let t = eq.trailing_ones();
                cur_hpoly = t;
                eq >>= t;
                max_hpoly = max_hpoly.max(cur_hpoly);
            }
        }
    }

    let max_hpoly = max_hpoly as usize;
    if max_hpoly > params.max_homopoly {
        return Status::BigHomopoly;
    }

    if params.min_entropy_1 > 0.0 {
        let entropy1 = compute_entropy(&count1);
        if entropy1 < params.min_entropy_1 {
            return Status::SmallEntropy1;
        }
    }

    if params.min_entropy_2 > 0.0 {
        let entropy2 = compute_entropy(&count2);
        if entropy2 < params.min_entropy_2 {
            return Status::SmallEntropy2;
        }
    }

    Status::Ok
}

#[inline(always)]
fn compute_entropy(count: &[usize]) -> f64 {
    let total = count.iter().sum::<usize>() as f64;
    if total == 0.0 {
        return 0.0;
    }
    count
        .iter()
        .map(|&c| {
            let p = c as f64 / total;
            p * p.log2()
        })
        .sum()
}

fn main() {
    let args = Args::parse();
    let _threads = if let Some(t) = args.threads {
        rayon::ThreadPoolBuilder::new()
            .num_threads(t)
            .build_global()
            .unwrap();
        t
    } else {
        rayon::current_num_threads()
    };
    let files = fasta_files(&args.path).expect("Failed to get input files");

    files.into_par_iter().for_each(|filename| {
        filter_file(filename, args.mode);
        // TODO collect stats
    });
}
