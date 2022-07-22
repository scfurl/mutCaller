/*

cd ~/develop/mutCaller/mutcaller_rust
/Users/sfurlan/develop/mutCaller/mutcaller_rust/target/debug/fastq3 --ifastq "/Users/sfurlan/Fred Hutchinson Cancer Research Center/Furlan_Lab - General/experiments/patient_marrows/LKmut/sample_filtered.fastq_1M.fastq.gz"
*/

/*
### TO DO ###
add option for running in debug mode
measure speed
*/

extern crate clap;
extern crate fastq;
extern crate parasailors;

use bytes::BytesMut;
use std::sync::{Arc, Mutex};
use fastq::{parse_path, Record, RefRecord};
use parasailors as align;
use clap::{App, load_yaml};
use std::io::{self, Write};
use std::io::stdout;
use std::vec::Vec;



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
    let mut total: usize = 0;
    let split_at: u8 = 116 - (params.umi_len + params.cb_len);
    let filename: Option<String> = Some(params.ifastq.to_string());
    // Treat "-" as stdin

    let path = match filename.as_ref().map(String::as_ref) {
        None | Some("-") => { None },
        Some(name) => Some(name)
    };
    parse_path(path, |mut parser| {
        
        let stopped = parser.each(|record| {
            // stop parsing if we find a sequnce containing 'N'
            record.validate_dnan();
            modify_fastq(record, &mut total, split_at)
        }).expect("Invalid fastq file");
        if stopped {
            println!("The file contains only sequences with ACTGN");
        } else {
            println!("The file contains invalid sequences");
        }
    }).expect("Invalid compression");
}



fn modify_fastq (record: RefRecord, total: &mut usize, split_at: u8) -> bool {
    let mut writer = io::stdout();
    let mut owned_rec = RefRecord::to_owned_record(&record);
    let mut curr_bytes = BytesMut::from(owned_rec.qual());
    let (seq, barcode) = &curr_bytes.split_at(split_at.into());
    owned_rec.qual = seq.to_vec();
    let mut curr_bytes = BytesMut::from(owned_rec.seq());
    let (seq, mut barcode) = &curr_bytes.split_at(split_at.into());
    owned_rec.seq = seq.to_vec();
    let mut curr_bytes = BytesMut::from(owned_rec.head()).to_vec();
    let _ = &curr_bytes.push(':' as u8);
    let _ = &curr_bytes.append(&mut barcode.to_vec());
    owned_rec.head = curr_bytes;
    owned_rec.sep = Some(vec!['+' as u8]);
    owned_rec.write(&mut writer);
    *total += 1;
    true
}

fn modify_fastq_debug(record: RefRecord, total: &mut usize, split_at: u8) -> bool {
    println!("Before mod");
    println!("{}", String::from_utf8_lossy(record.head()));
    println!("{}", String::from_utf8_lossy(record.seq()));
    println!("{}", String::from_utf8_lossy(record.qual()));
    let mut owned_rec = RefRecord::to_owned_record(&record);
    let mut curr_bytes = BytesMut::from(owned_rec.qual());
    let (seq, barcode) = &curr_bytes.split_at(split_at.into());
    owned_rec.qual = seq.to_vec();
    let mut curr_bytes = BytesMut::from(owned_rec.seq());
    let (seq, mut barcode) = &curr_bytes.split_at(split_at.into());
    owned_rec.seq = seq.to_vec();
    let mut curr_bytes = BytesMut::from(owned_rec.head()).to_vec();
    let _ = &curr_bytes.push(':' as u8);
    let _ = &curr_bytes.append(&mut barcode.to_vec());
    owned_rec.head = curr_bytes;
    owned_rec.sep = Some(vec!['+' as u8]);
    // owned_rec.head = &curr_bytes.append(&mut barcode.to_vec());
    // owned_rec.head = curr_bytes.into_iter().chain(barcode.into_iter()).collect();
    // owned_rec.head = curr_bytes.to_vec();
    println!("After mod");
    println!("{}", String::from_utf8_lossy(owned_rec.head()));
    println!("{}", String::from_utf8_lossy(owned_rec.seq()));
    println!("{}", String::from_utf8_lossy(owned_rec.qual()));
    *total += 1;
    true
}




// use bytes::BytesMut;
// use fastq::{Parser, Record, RefRecord};

// const READS: &str = r#"@read1/ENST00000266263.10;mate1:84-183;mate2:264-363
// GACAGCCAGGGGCCAGCGGGTGGCAGTGCCCAGGACATAGAGAGAGGCAGCACACACGCGGTTGATGGTGAAGCCCGGAATGGCCACAGAGGCTAGAGCC
// +
// IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
// @read2/ENST00000266263.10;mate1:163-262;mate2:283-382
// GATGCCATTGACAAAGGCAAGAAGGCTGGAGAGGTGCCCAGCCCTGAAGCAGGCCGCAGCGCCAGGGTGACTGTGGCTGTGGTGGACACCTTTGTATGGC
// +
// IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
// @read3/ENST00000266263.10;mate1:86-185;mate2:265-364
// GGACAGCCAGGGGCCAGCGGGTGGCAGTGCCCAGGACATAGAGAGAGGCAGCANACACACGGTTGATGGTGAAGCCCGGAATGGCCACAGAGGCTAGAGC
// +
// IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII!IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
// @read4/ENST00000266263.10;mate1:297-396;mate2:401-500
// CAGGAGGAGCTGGGCTTCCCCACTGTTAGGTAGAGCTTGCGCAGGCTGGAGTCCAGGAGGAAATCCACCGACCTGTCAATGGGGTGGATAATGATGGGGA
// +
// IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
// "#;