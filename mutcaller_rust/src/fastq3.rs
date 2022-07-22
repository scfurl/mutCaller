/*

cd ~/develop/mutCaller/mutcaller_rust
/Users/sfurlan/develop/mutCaller/mutcaller_rust/target/debug/fastq3 --ifastq "/Users/sfurlan/Fred Hutchinson Cancer Research Center/Furlan_Lab - General/experiments/patient_marrows/LKmut/sample_filtered.fastq_1M.fastq.gz"
*/


extern crate clap;
extern crate fastq;
extern crate parasailors;

use bytes::BytesMut;
use std::sync::{Arc, Mutex};
use fastq::{parse_path, Record, RefRecord};
// use std::env::args;
use parasailors as align;
use clap::{App, load_yaml};
use std::io::{self, Write};
use std::io::stdout;
use std::vec;
use std::str;

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


//     let parser = Parser::new(READS.as_bytes());
//     parser
//         .each(|record| {
//             modify_qual1(record, &mut total)
//             // println!("Before mod");
//             // println!("{}", String::from_utf8_lossy(record.qual()));
//             // let owned_rec = modify_qual(record);
//             // println!("After mod");
//             // println!("{}", String::from_utf8_lossy(owned_rec.qual()));
//             // total += 1;
//             // true
//         })
//         .expect("Invalid fastq file");
//     println!("{}", total);
// }


// fn modify_qual(record: RefRecord) -> fastq::OwnedRecord {
//     let mut owned_rec = RefRecord::to_owned_record(&record);
//     let mut curr_bytes = BytesMut::from(owned_rec.qual());
//     curr_bytes[0] = b'$';
//     owned_rec.qual = curr_bytes.to_vec();
//     owned_rec
// }

fn modify_fastq(record: RefRecord, total: &mut usize, split_at: u8) -> bool {
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
    owned_rec.seq = b"AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT".to_vec();
    let mut curr_bytes = BytesMut::from(owned_rec.head());
    owned_rec.head = &mut curr_bytes.to_vec().append(&mut barcode.clone().to_vec());
    println!("After mod");
    println!("{}", String::from_utf8_lossy(owned_rec.head()));
    println!("{}", String::from_utf8_lossy(owned_rec.seq()));
    println!("{}", String::from_utf8_lossy(owned_rec.qual()));
    *total += 1;
    true
}