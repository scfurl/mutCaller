/**

chmod 777 /Users/sfurlan/develop/mutCaller/mutcaller_rust/target/debug/mutcaller_rust
cd /Users/sfurlan/develop/mutCaller/mutcaller_rust
cargo build
/Users/sfurlan/develop/mutCaller/mutcaller_rust/target/debug/mutcaller_rust -t 18 --ibam "/Users/sfurlan/Fred Hutchinson Cancer Research Center/Furlan_Lab - General/experiments/patient_marrows/LKmut/Aligned.out.tagged.sorted_0.05.bam"

/Users/sfurlan/develop/mutCaller/mutcaller_rust/target/debug/fastq -t 8 --ifastq "/Users/sfurlan/Fred Hutchinson Cancer Research Center/Furlan_Lab - General/experiments/patient_marrows/LKmut/sample_filtered.fastq_1M.fastq.gz"

**/


extern crate clap;
extern crate fastq;
extern crate parasailors;

use std::sync::{Arc, Mutex};
use fastq::{parse_path, Record};
// use std::env::args;
use parasailors as align;
use clap::{App, load_yaml};
use std::io::{self, Write};
use std::io::stdout;

#[derive(Clone)]
struct Params {
    ifastq: String,
    umi_len: u8,
    cb_len: u8,
    threads: u8,
}

fn load_params() -> Params {
    let yaml = load_yaml!("params_fastq.yml");
    let params = App::from_yaml(yaml).get_matches();
    let ifastq = params.value_of("ifastq").unwrap();
    let threads = params.value_of("threads").unwrap_or("1");
    let threads = threads.to_string().parse::<u8>().unwrap();
    let umi_len = params.value_of("umi_len").unwrap_or("10");
    let umi_len = umi_len.to_string().parse::<u8>().unwrap() - 1;
    let cb_len = params.value_of("cb_len").unwrap_or("16");
    let cb_len = cb_len.to_string().parse::<u8>().unwrap() - 1;
    Params{
        ifastq: ifastq.to_string(),
        threads: threads as u8,
        umi_len: umi_len as u8,
        cb_len: cb_len as u8,
    }
}

fn main() {
    let params = load_params();
    fastq(&params);
}

fn fastq(params: &Params) {
    let filename: Option<String> = Some(params.ifastq.to_string());
    // Treat "-" as stdin
    let path = match filename.as_ref().map(String::as_ref) {
        None | Some("-") => { None },
        Some(name) => Some(name)
    };


    parse_path(path, |parser| {
        let nthreads = params.threads as usize;
        let results: Vec<usize> = parser.parallel_each(nthreads, |record_sets| {
            // we can initialize thread local variables here.
            let adapter = b"AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT";
            let matrix = align::Matrix::new(align::MatrixType::Identity);
            let profile = align::Profile::new(adapter, &matrix);
            let mut thread_total = 0;
            let handle = std::io::stdout();
            let buf = Mutex::new(BufWriter::new(handle));
            // let mut writer = io::stdout(); // get the global stdout entity
            // let mut handle = io::BufWriter::new(writer); 
            for record_set in record_sets {
                for record in record_set.iter() {
                    let score = align::local_alignment_score(&profile, record.seq(), 5, 1);
                    if score > 8 {
                        thread_total += 1;
                    }
                    record.write(buf.lock().unwrap());
                }
            }

            // The values we return (it can be any type implementing `Send`)
            // are collected from the different threads by
            // `parser.parallel_each` and returned. See doc for a description of
            // error handling.
            thread_total
        }).expect("Invalid fastq file");
        println!("{}", results.iter().sum::<usize>());
    }).expect("Invalid compression");
}