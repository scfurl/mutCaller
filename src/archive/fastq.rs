/*

cd ~/develop/mutCaller/mutcaller_rust
/Users/sfurlan/develop/mutCaller/mutcaller_rust/target/debug/fastq --ifastq "/Users/sfurlan/Fred Hutchinson Cancer Research Center/Furlan_Lab - General/experiments/patient_marrows/LKmut/sample_filtered.fastq_1M.fastq.gz" > test.fastq
cat fastq.log
tail test.fastq
head test.fastq

/Users/sfurlan/develop/mutCaller/mutcaller_rust/target/debug/fastq -t 8 --ifastq "/Users/sfurlan/Fred Hutchinson Cancer Research Center/Furlan_Lab - General/experiments/patient_marrows/LKmut/sample_filtered.fastq_1M.fastq.gz" > test.fastq
cat fastq.log
tail test.fastq
head test.fastq

/Users/sfurlan/develop/mutCaller/mutcaller_rust/target/debug/fastq --ifastq "/Users/sfurlan/Fred Hutchinson Cancer Research Center/Furlan_Lab - General/experiments/patient_marrows/LKmut/sample_filtered.fastq_1M.fastq.gz"


cd ~/develop/mutCaller/mutcaller_rust
/Users/sfurlan/develop/mutCaller/mutcaller_rust/target/debug/fastq --ifastq "/Users/sfurlan/Fred Hutchinson Cancer Research Center/Furlan_Lab - General/experiments/patient_marrows/LKmut/sample_filtered.fastq_1M.fastq.gz" > test.fastq
cat fastq.log
tail test.fastq
head test.fastq

/Users/sfurlan/develop/mutCaller/mutcaller_rust/target/debug/fastq -t 8 --ifastq "/Users/sfurlan/Fred Hutchinson Cancer Research Center/Furlan_Lab - General/experiments/patient_marrows/LKmut/sample_filtered.fastq_1M.fastq.gz" > test.fastq
cat fastq.log
tail test.fastq
head test.fastq

/Users/sfurlan/develop/mutCaller/mutcaller_rust/target/debug/fastq --ifastq "/Users/sfurlan/Fred Hutchinson Cancer Research Center/Furlan_Lab - General/experiments/patient_marrows/LKmut/sample_filtered.fastq_1M.fastq.gz"

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
use fastq::{parse_path, Record, RefRecord};
use clap::{App, load_yaml};
use std::io::{self};

struct Params {
    ifastq: String,
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
    let yaml = load_yaml!("params_fastq.yml");
    let params = App::from_yaml(yaml).get_matches();
    let debug = params.value_of("debug").unwrap_or("false");
    let debug = debug.to_string().parse::<bool>().unwrap();
    let ifastq = params.value_of("ifastq").unwrap();
    let threads = params.value_of("threads").unwrap_or("1");
    let threads = threads.to_string().parse::<u8>().unwrap();
    let umi_len = params.value_of("umi_len").unwrap_or("10");
    let umi_len = umi_len.to_string().parse::<u8>().unwrap();
    let cb_len = params.value_of("cb_len").unwrap_or("16");
    let cb_len = cb_len.to_string().parse::<u8>().unwrap();
    let _max_reads = params.value_of("n_reads").unwrap_or("all");
    let name_sep = params.value_of("name_sep").unwrap_or("|BARCODE=");
    if _max_reads == "all"{
        let _stop = false;
        let _max_r = 0usize;
    }else{
        let _stop = true;
        let _max_r = _max_reads.to_string().parse::<usize>().unwrap();
    }
    Params{
        ifastq: ifastq.to_string(),
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
        .path("./fastq.log")
        .size(1 * 100)
        .roll_count(10)
        .time_format("%Y-%m-%d %H:%M:%S.%f") //E.g:%H:%M:%S.%f
        .level("debug")
        .output_file()
        .build();

    let _ = simple_log::new(config);
    info!("starting!");
    let params = load_params();
    eprintln!("Stopping at {}", params.stop);
    eprintln!("Running with {} thread(s)!", params.threads);
    if params.threads > 1{
        fastq_parallel(params)
    }else{
        fastq(&params)
    }
    info!("done!");
}




fn fastq(params: &Params) {
    let mut total: usize = 0;
    let split_at = params.umi_len + params.cb_len;
    let sep: Vec::<u8> = params.name_sep.as_bytes().to_vec();
    let filename: Option<String> = Some(params.ifastq.to_string());
    // Treat "-" as stdin

    let path = match filename.as_ref().map(String::as_ref) {
        None | Some("-") => { None },
        Some(name) => Some(name)
    };
    parse_path(path, |parser| {
        let _max_reads = params.max_reads;
        let stopped = parser.each( |record| {
            // stop parsing if we find a sequnce containing 'N'
            // let valid = record.validate_dnan();
                if params.debug{
                    modify_fastq_debug(record, split_at, &sep);
                }else{
                    modify_fastq(record, split_at, &sep);
                };
                total += 1;
                true
        }).expect("Invalid fastq file");
        if stopped {
            eprintln!("Completed; {} reads processed!", total);
        } else {
            eprintln!("The file contains invalid sequences");
        }
    }).expect("Invalid compression");
}

fn fastq_parallel(params: Params) {
    //BROKEN
    let filename: Option<String> = Some(params.ifastq.to_string());
    // Treat "-" as stdin
    let path = match filename.as_ref().map(String::as_ref) {
        None | Some("-") => { None },
        Some(name) => Some(name)
    };
    parse_path(path, |parser| {
        let results: Vec::<usize> = parser.parallel_each(params.threads, |record_sets| {
            // let params = params;
            // let split_at: u8 = 116 - (params.umi_len + params.cb_len);
            // let debug: bool = params.debug;
            let mut thread_total = 0;
            for record_set in record_sets {
                for record in record_set.iter() {
                    // let valid = record.validate_dnan(); FIXME LATER
                    let sep: Vec::<u8> = "|BARCODE=".as_bytes().to_vec();
                    let valid = modify_fastq(record, 90, &sep);
                    if valid {
                        thread_total += 1;
                    }
                        // if params.debug{
                        //     let valid = modify_fastq_debug(record, 90);
                        //     if valid {
                        //         thread_total += 1;
                        //     }
                        // }else{
                        //     let valid = modify_fastq(record, 90);
                        //     if valid {
                        //         thread_total += 1;
                        //     }
                        // }
                }
            }
            // The values we return (it can be any type implementing `Send`)
            // are collected from the different threads by
            // `parser.parallel_each` and returned. See doc for a description of
            // error handling.
            thread_total
        }).expect("Invalid fastq file");
        eprintln!("Completed; {} reads processed!", results.iter().sum::<usize>());
    }).expect("Invalid compression");
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