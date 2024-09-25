/*
    Input ASO sequence in 5' -> 3' orientation
        - import input sequences
    Check against a tab separated file of ASOs see
        - [X] Similar ATGC content
        - [X] Levenshtein distance
        - [X] Hamming distance
        - [X] sift3
*/
use std::fs::File;
use std::io;
use std::path::PathBuf;
use clap::{Parser, ValueEnum, ArgAction};
use csv::{ReaderBuilder, StringRecordsIter, Trim};
use log::{debug, info, warn};
use std::error::Error;
use std::rc::Rc;
use distance::{hamming, levenshtein, sift3};

#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
pub struct Cli {
    /// Input ASO sequence. One sequence, in 5' -> 3' orientation
    #[arg(short='a', long="aso-seq")]
    aso_seq: Option<String>,
    /// Process multiple ASO sequences. Conflicts with aso-seq option.
    /// Requires input-aso-file argument
    #[arg(short='m', long= "multiple-aso-seq", conflicts_with = "aso_seq")]
    multiple_aso: bool,
    /// path to input ASO sequences
    /// in csv format, ASO name in column1,
    /// ASO sequences in 5' -> 3' orientation in column2
    /// Any additional information can be entered in lines
    /// starting with #. They won't be read.
    #[arg(long="input-aso-file", requires = "multiple_aso", conflicts_with = "aso_seq")]
    input_aso_file: Option<PathBuf>,
    /// no headers in the input file
    #[arg(long="input-no-header", requires = "input_aso_file",
    action=ArgAction::SetFalse, group = "multi-aso", conflicts_with = "aso_seq")]
    input_header_status: bool,
    /// path to library of existing ASOs
    /// in csv format, ASO name in column1
    /// ASO sequence in 5' -> 3' orientation in column2
    /// Any additional information can be entered in lines
    /// starting with #. They won't be read.
    #[arg(short='l', long="library-aso-file", name="libfile", required = true)]
    library_aso_file: PathBuf,
    /// no headers in the library file
    #[arg(long="library-no-header", name="lib_header", requires = "libfile",
    action=ArgAction::SetFalse)]
    library_header_status: bool,
    /// Display only which distance? Default: Levenshtein
    /// Distance: higher the number, greater the mismatch between sequences
    #[arg(long="list-by", name="List",
    value_enum, ignore_case = true, default_value_t= Dist::Levenshtein)]
    list_by: Dist,
}
#[derive(Debug, PartialEq, Copy, Clone, ValueEnum)]
pub enum Dist {
    Hamming,
    Levenshtein,
    Sift3
}

fn main() {
    env_logger::init(); // Start logging based on the RUST_LOG parameter
    debug!("Parsing commandline arguments");
    let cli = Cli::parse();
    let run_multiple_mode = cli.multiple_aso;
    let library_file_path = cli.library_aso_file.clone();
    info!("Initialising library of ASOs");
    let library_header_status = cli.library_header_status;
    if !library_header_status {
        warn!("Note: Library file has no header. First entry will be processed")
    } else {
        warn!("Note: Library has header, first entry will not be processed.")
    }
    let mut aso_library_reader = ReaderBuilder::new()
        .has_headers(library_header_status)
        .from_path(library_file_path)
        .expect("Unable to open library file. Closing");
    match run_multiple_mode {
        true => {
            debug!("Processing multiple ASO sequences");
            let aso_input_file_path = cli.input_aso_file.clone().unwrap();
            info!("Processing input ASO file {:?}", aso_input_file_path.as_path());
            let input_file_header = cli.input_header_status;
            if !input_file_header {
                warn!("Note: Input ASO file has no header. First entry will be processed")
            } else {
                warn!("Note: Library has header, first entry will not be processed.")
            }
            let mut input_aso_reader = ReaderBuilder::new()
                .has_headers(input_file_header)
                .trim(Trim::All)
                .from_path(aso_input_file_path)
                .expect("Unable to open input ASO file");
            let _ = compute_distance(aso_library_reader.records(), &cli, input_aso_reader.records());
        }
        false => {
            let aso_seq = cli.aso_seq.clone()
                .expect("Enter ASO sequence or provide file path to ASOs");
            debug!("Processing the given input ASO sequence: {}", aso_seq);
            info!("Naming the input ASO {} as testASO_001", aso_seq);
            let aso_input = format!("testASO_001, {}", aso_seq);
            let mut input_aso_reader = ReaderBuilder::new()
                .has_headers(false)
                .trim(Trim::All)
                .from_reader(aso_input.as_bytes());
            let _ = compute_distance(aso_library_reader.records(), &cli, input_aso_reader.records());
        }
    }

}

struct AsoProfile {
    name: String,
    seq: String,
    aso_len: usize,
    atgc: [usize; 4],
    // aso_names: Vec<(String, f32)>
    aso_names: Vec<(Rc<AsoProfile>, f32)>
}

impl AsoProfile {
    fn new(name: String, seq: String) -> Self {
        let name = name;
        let seq = seq;
        let aso_len = seq.len();
        let atgc = atgc_count(&seq);
        AsoProfile {
            name,
            seq,
            aso_len,
            atgc,
            aso_names: vec![],
        }
    }
}

fn atgc_count(seq: &str) -> [usize; 4] {
    let mut count_n = [0; 4];
    count_n[0] += char_windows(seq, 1)
        .filter(|c| c == &"A")
        .count();
    count_n[1] += char_windows(seq, 1)
        .filter(|c| c == &"T")
        .count();
    count_n[2] += char_windows(seq, 1)
        .filter(|c| c == &"G")
        .count();
    count_n[3] += char_windows(seq, 1)
        .filter(|c| c == &"C")
        .count();
    count_n
}

fn char_windows<'a>(src: &'a str, win_size: usize) -> impl Iterator<Item = &'a str> {
    src.char_indices().flat_map(move |(from, _)| {
        src[from..]
            .char_indices()
            .skip(win_size - 1)
            .next()
            .map(|(to, c)| &src[from..from + to + c.len_utf8()])
    })
}

fn compute_distance<R: io::Read>(library: StringRecordsIter<File>, cli: &Cli,
                                 input: StringRecordsIter<R>) -> Result<(), Box<dyn Error>> {
    // compute the ATGC spread of each input source
    // compute the ATGC spread of each library source
    // if ATGC and length match found, calculate all three distances
    let mut input_seq_props: Vec<AsoProfile> = Vec::new();
    let mut library_asos: Vec<Rc<AsoProfile>> = Vec::new();
    let list_method = cli.list_by;
    for input_result in input {
        let record = input_result?;
        if record.len() < 2 {
            panic!("Incomplete file. Name and sequence necessary")
        }
        let name = record.get(0).expect("No name").to_string();
        let seq = record.get(1).expect("No seq").to_string();
        let aso_profile = AsoProfile::new(name, seq);
        input_seq_props.push(aso_profile)
    }
    for library_result in library {
        let record = library_result?;
        if record.len() < 2 {
            panic!("Incomplete file. Name and sequence necessary")
        }
        let seq = record.get(1).expect("No seq").to_string();
        let name = record.get(0).expect("No name").to_string();
        let aso_profile = AsoProfile::new(name, seq);
        let aso_rc = Rc::new(aso_profile);
        library_asos.push(aso_rc.clone());
        let aso_profile = aso_rc;
        input_seq_props.iter_mut().for_each(|in_aso| {
            if in_aso.aso_len == aso_profile.aso_len && in_aso.atgc == aso_profile.atgc
                && in_aso.seq != aso_profile.seq {
                let dist = match list_method {
                    Dist::Hamming => hamming(&in_aso.seq, &aso_profile.seq).unwrap() as f32,
                    Dist::Levenshtein => levenshtein(&in_aso.seq, &aso_profile.seq) as f32,
                    Dist::Sift3 => sift3(&in_aso.seq, &aso_profile.seq)
                };
                in_aso.aso_names.push((aso_profile.clone(), dist))
            }
        })
    }
    println!("{:<10}\t{:<20}\t{:<10}\t{:<20}\t{}", "Input ASO","Seq", "Matching ASO", "Seq", "Distance");
    for aso in input_seq_props.iter_mut() {
        println!("{:<10}\t{:<20}", aso.name, aso.seq);
        aso.aso_names
            .sort_unstable_by(|(_, a), (_, b)|
                a.partial_cmp(b).unwrap());
        for (scramble, distance) in &aso.aso_names {
            println!("{:<10}\t{:<20}\t{:<10}\t{:<20}\t{}", "", "", scramble.name, scramble.seq, distance)
        }
    }
    Ok(())
}

