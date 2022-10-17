/*

time ~/develop/mutCaller/mutcaller_rust/target/release/fastq1 -t 1 --fastq1 tests/sequencer_R1.fastq.gz --fastq2 tests/sequencer_R2.fastq.gz


*/

/*
### TO DO ###
measure speed
*/

#[macro_use]
extern crate simple_log;
extern crate clap;
extern crate fastq;
// extern crate parasailors;

use simple_log::LogConfigBuilder;
use bytes::BytesMut;
use fastq::{parse_path, Record, RefRecord, each_zipped};
use clap::{App, load_yaml};
use std::io::{self};

struct Params {
    fastq1: String,
    fastq2: String,
    umi_len: u8,
    cb_len: u8,
    threads: usize,
    max_reads: usize,
    stop: bool,
    debug: bool,
    name_sep: String,
}

fn load_params() -> Params {
    let max_r = 0usize;
    let stop = true;
    let yaml = load_yaml!("params_fastq1.yml");
    let params = App::from_yaml(yaml).get_matches();
    let debug = params.value_of("debug").unwrap_or("false");
    let debug = debug.to_string().parse::<bool>().unwrap();
    let fastq1 = params.value_of("fastq1").unwrap();
    let fastq2 = params.value_of("fastq2").unwrap();
    let threads = params.value_of("threads").unwrap_or("1");
    let threads = threads.to_string().parse::<u8>().unwrap();
    let umi_len = params.value_of("umi_len").unwrap_or("10");
    let umi_len = umi_len.to_string().parse::<u8>().unwrap();
    let cb_len = params.value_of("cb_len").unwrap_or("16");
    let cb_len = cb_len.to_string().parse::<u8>().unwrap();
    let _max_reads = params.value_of("n_reads").unwrap_or("all");
    let name_sep = params.value_of("name_sep").unwrap_or("|BARCODE=");
    if _max_reads == "all"{
        let stop = false;
        let max_r = 0usize;
    }else{
        let stop = true;
        let max_r = _max_reads.to_string().parse::<usize>().unwrap();
    }
    Params{
        fastq1: fastq1.to_string(),
        fastq2: fastq1.to_string(),
        threads: threads as usize,
        umi_len: umi_len as u8,
        cb_len: cb_len as u8,
        max_reads: max_r,
        stop: stop,
        debug: debug,
        name_sep: name_sep.to_string(),
    }
}

fn main() {
    let config = LogConfigBuilder::builder()
        .path("./fastq1.log")
        .size(1 * 100)
        .roll_count(10)
        .time_format("%Y-%m-%d %H:%M:%S.%f") //E.g:%H:%M:%S.%f
        .level("debug")
        .output_file()
        .build();

    let _ = simple_log::new(config);
    info!("starting!");
    let params = load_params();
    eprintln!("Running with {} thread(s)!", &params.threads);
    fastq(&params);
    info!("done!");
}



fn fastq(params: &Params) {
    let mut total: usize = 0;
    let split_at = &params.umi_len + &params.cb_len;
    let sep: Vec::<u8> = "|BARCODE=".as_bytes().to_vec();
    let fastq1 = &params.fastq1;
    let fastq2 = &params.fastq2;
    let mut counts = (0u64, 0u64);
    let mut writer = io::stdout();
    parse_path(Some(fastq1), |parser1| {
        parse_path(Some(fastq2), |parser2| {
            each_zipped(parser1, parser2, |rec1, rec2| {
                // if rec1.is_some() {
                //     counts.0 += 1;
                // }
                // if rec2.is_some() {
                //     counts.1 += 1;
                // }
                if rec1.is_some() & rec2.is_some(){
                    let good = eval_reads(&rec1.unwrap(), &rec2.unwrap());
                }
                // eval_reads(rec1, rec2)
                (true, true)
            })
            .expect("Invalid record.");
        })
        .expect("Unknown format for file 2.");
    })
    .expect("Unknown format for file 1.");

    // println!("Number of reads: ({}, {})", counts.0, counts.1);
}

pub fn eval_reads(rec1: &RefRecord, rec2: &RefRecord) -> bool{
    let mut writer = io::stdout();
    if rec2.seq().contains(&b"N"[0]){
        println!("true");
        false
    }else{
        println!("false");
        true
    }
    // println!("{:?}",rec1.seq())
}

fn fastq_count(params: &Params) {
    let mut total: usize = 0;
    let split_at = params.umi_len + params.cb_len;
    let sep: Vec::<u8> = "|BARCODE=".as_bytes().to_vec();
    let fastq1 = &params.fastq1;
    let fastq2 = &params.fastq2;
    let mut counts = (0u64, 0u64);
    parse_path(Some(fastq1), |parser1| {
        parse_path(Some(fastq2), |parser2| {
            each_zipped(parser1, parser2, |rec1, rec2| {
                if rec1.is_some() {
                    counts.0 += 1;
                }
                if rec2.is_some() {
                    counts.1 += 1;
                }
                (true, true)
            })
            .expect("Invalid record.");
        })
        .expect("Unknown format for file 2.");
    })
    .expect("Unknown format for file 1.");

    println!("Number of reads: ({}, {})", counts.0, counts.1);
}



fn modify_fastq (record: RefRecord, split_at: u8, namesep: &Vec<u8>) -> bool {
    let sep: &mut Vec::<u8> = &mut namesep.clone();
    let mut writer = io::stdout();
    let mut owned_rec = RefRecord::to_owned_record(&record);
    let curr_bytes = BytesMut::from(owned_rec.qual());
    let (_barcode, seq) = &curr_bytes.split_at(split_at.into());
    owned_rec.qual = seq.to_vec();
    let curr_bytes = BytesMut::from(owned_rec.seq());
    let (barcode, seq) = &curr_bytes.split_at(split_at.into());
    owned_rec.seq = seq.to_vec();
    let mut curr_bytes = BytesMut::from(owned_rec.head()).to_vec();
    // let _ = &curr_bytes.push(namesep);
    let _ = &curr_bytes.append(sep);
    let _ = &curr_bytes.append(&mut barcode.to_vec());
    owned_rec.head = curr_bytes;
    owned_rec.sep = Some(vec!['+' as u8]);
    let _ = &owned_rec.write(&mut writer);
    true
}


fn modify_fastq_debug(record: RefRecord, split_at: u8, namesep: &Vec<u8>) -> bool {
    let sep: &mut Vec::<u8> = &mut namesep.clone();
    println!("Before mod");
    println!("{}", String::from_utf8_lossy(record.head()));
    println!("{}", String::from_utf8_lossy(record.seq()));
    println!("{}", String::from_utf8_lossy(record.qual()));
    let mut owned_rec = RefRecord::to_owned_record(&record);
    let curr_bytes = BytesMut::from(owned_rec.qual());
    let (_barcode, seq) = &curr_bytes.split_at(split_at.into());
    owned_rec.qual = seq.to_vec();
    let curr_bytes = BytesMut::from(owned_rec.seq());
    let (barcode, seq) = &curr_bytes.split_at(split_at.into());
    owned_rec.seq = seq.to_vec();
    let mut curr_bytes = BytesMut::from(owned_rec.head()).to_vec();
    let _ = &curr_bytes.append(sep);
    let _ = &curr_bytes.append(&mut barcode.to_vec());
    owned_rec.head = curr_bytes;
    owned_rec.sep = Some(vec!['+' as u8]);
    println!("After mod");
    println!("{}", String::from_utf8_lossy(owned_rec.head()));
    println!("{}", String::from_utf8_lossy(owned_rec.seq()));
    println!("{}", String::from_utf8_lossy(owned_rec.qual()));
    true
}

