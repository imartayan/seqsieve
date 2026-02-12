use arrayvec::ArrayVec;
use clap::{Parser as ClapParser, ValueEnum};
use helicase::config::advanced::COMPUTE_DNA_STRING;
use helicase::dna_format::ColumnarDNA;
use helicase::input::*;
use helicase::*;
use indicatif::ProgressBar;
use jwalk::WalkDir;
use rayon::prelude::*;

use std::fmt::{Display, Formatter};
use std::fs::{self, File};
use std::io::{self, BufWriter, Write};
use std::path::{Path, PathBuf};
use std::sync::Mutex;

#[derive(ClapParser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
    /// Path to the file or directory to filter
    path: PathBuf,
    /// Output directory for filtered files
    #[arg(short, long)]
    outdir: PathBuf,
    /// Output file for statistics
    #[arg(short, long, default_value = "stats.tsv")]
    stats: PathBuf,
    /// Profile [default: all profiles]
    #[arg(short, long, value_enum)]
    profile: Option<Profile>,
    /// Number of threads [default: all]
    #[arg(short, long)]
    threads: Option<usize>,
}

const MAX_NUM_PROFILES: usize = 3;

type SVec<T> = ArrayVec<T, MAX_NUM_PROFILES>;

#[derive(ValueEnum, Clone, Copy, Debug)]
enum Profile {
    Safe,
    Interm,
    Aggressive,
    // Custom,
}

impl Display for Profile {
    #[inline(always)]
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        match self {
            Profile::Safe => write!(f, "safe"),
            Profile::Interm => write!(f, "interm"),
            Profile::Aggressive => write!(f, "aggressive"),
            // Profile::Custom => write!(f, "custom"),
        }
    }
}

#[derive(Debug)]
enum Status {
    Ok,
    TooShort,
    TooManyNs,
    LargeHomopoly,
    SmallEntropy1,
    SmallEntropy2,
    InvalidNuc,
}

#[derive(Debug)]
struct FilterStats {
    pub input: PathBuf,
    pub output: PathBuf,
    pub profile: Profile,
    pub seqs_in: usize,
    pub seqs_out: usize,
    pub bp_in: usize,
    pub bp_out: usize,
    pub drop_short: usize,
    pub drop_n: usize,
    pub drop_badchar: usize,
    pub drop_hpoly: usize,
    pub drop_ent1: usize,
    pub drop_ent2: usize,
}

impl FilterStats {
    #[inline(always)]
    const fn new(input: PathBuf, output: PathBuf, profile: Profile) -> Self {
        Self {
            input,
            output,
            profile,
            seqs_in: 0,
            seqs_out: 0,
            bp_in: 0,
            bp_out: 0,
            drop_short: 0,
            drop_n: 0,
            drop_badchar: 0,
            drop_hpoly: 0,
            drop_ent1: 0,
            drop_ent2: 0,
        }
    }
}

impl Display for FilterStats {
    #[allow(clippy::write_literal)]
    #[inline(always)]
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            self.input.display(),
            self.output.display(),
            self.profile,
            "ok",
            self.seqs_in,
            self.seqs_out,
            self.bp_in,
            self.bp_out,
            self.drop_short,
            self.drop_n,
            self.drop_badchar,
            self.drop_hpoly,
            self.drop_ent1,
            self.drop_ent2,
            0,
            0
        )
    }
}

#[derive(Debug)]
struct SeqThresholds {
    pub min_len: usize,
    pub max_n_frac: f64,
    pub max_homopoly: usize,
    pub min_entropy_1: f64,
    pub min_entropy_2: f64,
    pub strict_nuc: bool,
}

impl SeqThresholds {
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
    const fn from_profile(profile: &Profile) -> Self {
        match profile {
            Profile::Safe => Self::new(80, 0.05, 60, 0.8, 1.5, true),
            Profile::Interm => Self::new(60, 0.20, 120, 0.2, 0.4, false),
            Profile::Aggressive => Self::new(51, 0.50, 200, 0.0, 0.0, false),
            // _ => unimplemented!(),
        }
    }
}

#[derive(Debug)]
struct SeqMetrics {
    pub len: usize,
    pub n_frac: f64,
    pub homopoly: usize,
    pub entropy_1: f64,
    pub entropy_2: f64,
    pub invalid_nuc: bool,
}

impl SeqMetrics {
    #[inline(always)]
    const fn new(
        len: usize,
        n_frac: f64,
        homopoly: usize,
        entropy_1: f64,
        entropy_2: f64,
        invalid_nuc: bool,
    ) -> Self {
        Self {
            len,
            n_frac,
            homopoly,
            entropy_1,
            entropy_2,
            invalid_nuc,
        }
    }

    #[inline(always)]
    const fn verify(&self, t: &SeqThresholds) -> Status {
        if self.len < t.min_len {
            Status::TooShort
        } else if self.n_frac > t.max_n_frac {
            Status::TooManyNs
        } else if self.invalid_nuc && t.strict_nuc {
            Status::InvalidNuc
        } else if self.homopoly > t.max_homopoly {
            Status::LargeHomopoly
        } else if self.entropy_1 < t.min_entropy_1 {
            Status::SmallEntropy1
        } else if self.entropy_2 < t.min_entropy_2 {
            Status::SmallEntropy2
        } else {
            Status::Ok
        }
    }
}

struct OutParam {
    pub profile: Profile,
    pub thresholds: SeqThresholds,
    pub dir: PathBuf,
}

fn fasta_files<P: AsRef<Path>>(path: P) -> io::Result<Vec<PathBuf>> {
    let path = path.as_ref();
    if !path.exists() {
        return Err(io::Error::other("Path does not exist"));
    }
    if path.is_dir() {
        Ok(fasta_files_in_dir(path))
    } else {
        path.to_str()
            .ok_or(io::Error::other("Invalid path encoding"))
            .and_then(|name| {
                if is_fasta_file(name) {
                    Ok(())
                } else {
                    Err(io::Error::other("Invalid path extension"))
                }
            })?;
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

#[inline(always)]
fn basename<P: AsRef<Path>>(filename: &Path, root: P) -> &str {
    let basename = filename
        .strip_prefix(root)
        .expect("Failed to strip root dir from path")
        .to_str()
        .expect("Invalid path encoding");
    basename
        .strip_suffix(".gz")
        .or_else(|| basename.strip_suffix(".zst"))
        .unwrap_or(basename)
}

fn filter_file<P: AsRef<Path>>(
    filename: PathBuf,
    root: P,
    outparams: &[OutParam],
) -> SVec<FilterStats> {
    const CONFIG: Config = ParserOptions::default()
        .dna_columnar()
        .keep_non_actg()
        .config()
        | COMPUTE_DNA_STRING;

    let basename = basename(&filename, root);
    let mut outfiles: SVec<_> = outparams
        .iter()
        .filter_map(|o| {
            let out_filename = o.dir.join(basename);
            if out_filename.is_file() && fs::metadata(&out_filename).unwrap().len() > 0 {
                None
            } else {
                let dir = out_filename.parent().expect("Failed to get parent dir");
                fs::create_dir_all(dir).expect("Failed to create output dir");
                let tmp_filename = out_filename.with_added_extension("tmp");
                let file = File::create(&tmp_filename).expect("Failed to create tmp file");
                let writer = BufWriter::new(file);
                let stats = FilterStats::new(filename.clone(), out_filename, o.profile);
                Some((stats, writer, &o.thresholds))
            }
        })
        .collect();

    if !outfiles.is_empty() {
        let min_len = outfiles.iter().map(|(_, _, t)| t.min_len).min().unwrap();
        let mut parser =
            FastaParser::<CONFIG, _>::from_file(&filename).expect("Failed to parse file");
        while let Some(_) = parser.next() {
            let seq = parser.get_dna_columnar();
            outfiles.iter_mut().for_each(|(stats, _, _)| {
                stats.seqs_in += 1;
                stats.bp_in += seq.len();
            });
            if seq.len() < min_len {
                outfiles.iter_mut().for_each(|(stats, _, _)| {
                    stats.drop_short += 1;
                });
                continue;
            }
            let metrics = seq_metrics(seq);
            outfiles
                .iter_mut()
                .for_each(|(stats, writer, t)| match metrics.verify(t) {
                    Status::Ok => {
                        stats.seqs_out += 1;
                        stats.bp_out += seq.len();
                        writer.write_all(b">").unwrap();
                        writer.write_all(parser.get_header()).unwrap();
                        writer.write_all(b"\n").unwrap();
                        writer.write_all(parser.get_dna_string()).unwrap();
                        writer.write_all(b"\n").unwrap();
                    }
                    Status::TooShort => {
                        stats.drop_short += 1;
                    }
                    Status::TooManyNs => {
                        stats.drop_n += 1;
                    }
                    Status::LargeHomopoly => {
                        stats.drop_hpoly += 1;
                    }
                    Status::SmallEntropy1 => {
                        stats.drop_ent1 += 1;
                    }
                    Status::SmallEntropy2 => {
                        stats.drop_ent2 += 1;
                    }
                    Status::InvalidNuc => {
                        stats.drop_badchar += 1;
                    }
                });
        }
    }

    let stats: SVec<_> = outfiles
        .into_iter()
        .map(|(stats, _, _)| {
            let out = &stats.output;
            let tmp = out.with_added_extension("tmp");
            fs::rename(tmp, out).expect("Failed to move tmp file");
            stats
        })
        .collect();
    stats
}

#[inline(always)]
fn seq_metrics(seq: &ColumnarDNA) -> SeqMetrics {
    let mut count1 = [0; 4];
    let mut count2 = [0; 16];
    let mut n_count = 0;
    let mut max_hpoly = 1;
    let mut cur_hpoly = 0;
    let mut invalid_nuc = false;

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

        count1[0] += va.count_ones() as usize;
        count1[1] += vc.count_ones() as usize;
        count1[2] += vt.count_ones() as usize;
        count1[3] += vg.count_ones() as usize;

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

        carry_hi = hi >> 63;
        carry_lo = lo >> 63;
    }

    let rem = seq.len() % 64;
    if rem > 0 {
        let mask = !0 >> (64 - rem);
        // TODO check N mask

        let prev_hi = (hi << 1) | carry_hi;
        let prev_lo = (lo << 1) | carry_lo;
        let va = (!hi & !lo) & mask;
        let vc = (!hi & lo) & mask;
        let vt = (hi & !lo) & mask;
        let vg = (hi & lo) & mask;

        count1[0] += va.count_ones() as usize;
        count1[1] += vc.count_ones() as usize;
        count1[2] += vt.count_ones() as usize;
        count1[3] += vg.count_ones() as usize;

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

    SeqMetrics::new(
        seq.len(),
        n_count as f64 / seq.len() as f64,
        max_hpoly as usize,
        compute_entropy(&count1),
        compute_entropy(&count2),
        invalid_nuc,
    )
}

#[inline(always)]
fn compute_entropy(count: &[usize]) -> f64 {
    let total = count.iter().sum::<usize>() as f64;
    if total == 0.0 {
        return 0.0;
    }
    count
        .iter()
        .copied()
        .filter(|&c| c > 0)
        .map(|c| {
            let p = c as f64 / total;
            -p * p.log2()
        })
        .sum()
}

fn main() {
    let args = Args::parse();

    if let Some(t) = args.threads {
        rayon::ThreadPoolBuilder::new()
            .num_threads(t)
            .build_global()
            .unwrap();
    }

    let root = if args.path.is_dir() {
        args.path.clone()
    } else {
        args.path
            .parent()
            .expect("Failed to get parent dir")
            .to_path_buf()
    };
    let files = fasta_files(&args.path).expect("Failed to get input files");
    eprintln!("Listed {} files to process", files.len());

    let mut profiles = SVec::new_const();
    if let Some(profile) = args.profile {
        profiles.push(profile);
    } else {
        profiles
            .try_extend_from_slice(&[Profile::Safe, Profile::Interm, Profile::Aggressive])
            .unwrap();
    };
    let outparams: SVec<_> = profiles
        .into_iter()
        .map(|profile| {
            let thresholds = SeqThresholds::from_profile(&profile);
            let dir = args.outdir.join(profile.to_string());
            fs::create_dir_all(&dir).expect("Failed to create output dir");
            OutParam {
                profile,
                thresholds,
                dir,
            }
        })
        .collect();

    if !args.stats.is_file() {
        let mut stats_file = File::create(&args.stats).expect("Failed to create stats file");
        stats_file.write_all(b"input_fasta\toutput_fasta\tprofile\tstatus\tseqs_in\tseqs_out\tbp_in\tbp_out\tdrop_short\tdrop_n\tdrop_badchar\tdrop_hpoly\tdrop_ent1\tdrop_ent2\trc\terr\n").expect("Failed to write header");
    };
    let stats_file = File::options()
        .append(true)
        .open(&args.stats)
        .expect("Failed to open stats file in append mode");
    let stats_writer = Mutex::new(BufWriter::new(stats_file));

    let pb = ProgressBar::new(files.len() as u64);
    files.into_par_iter().for_each(|filename| {
        let stats = filter_file(filename, &root, &outparams);
        let mut writer = stats_writer.lock().unwrap();
        stats.into_iter().for_each(|s| {
            writeln!(writer, "{}", s).expect("Failed to append stats");
        });
        pb.inc(1);
    });
    pb.finish_with_message("Done");
}
