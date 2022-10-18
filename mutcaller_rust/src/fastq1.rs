/*

time ~/develop/mutCaller/mutcaller_rust/target/release/fastq1 -t 1 --fastq1 tests/sequencer_R1.fastq.gz --fastq2 tests/sequencer_R2.fastq.gz | gzip > tests/out1.fq.gz
real    0m1.219s
echo $(zcat < tests/out1.fq.gz | wc -l)/4|bc


time ~/develop/mutCaller/mutcaller_rust/target/release/fastq1 -t 1 --fastq1 tests/sequencer_R1.fastq.gz --fastq2 tests/sequencer_R2.fastq.gz -o tests/out1.fq.gz --barcodes_file /home/sfurlan/develop/mutCaller/data/737K-august-2016.txt.gz

#compare to original mutcaller
cd ~/develop/mutCaller
time ~/develop/mutCaller/mutcaller -u -U 10 -b ~/develop/mutCaller/mutcaller_rust/tests/sequencer_R1.fastq.gz -t ~/develop/mutCaller/mutcaller_rust/tests/sequencer_R2.fastq.gz -l ~/develop/mutCaller/data/737K-august-2016.txt.gz &&
gzip tests/fastq_processed/sample_filtered.fastq
# real  0m23.133s
*/





#[macro_use]
extern crate simple_log;
extern crate clap;
extern crate fastq;
extern crate flate2;
// extern crate parasailors;

use simple_log::LogConfigBuilder;
use bytes::BytesMut;
// use fastq::{parse_path, Record, RefRecord, OwnedRecord, each_zipped};
use fastq::{parse_path, Record, RefRecord, each_zipped};
use clap::{App, load_yaml};
use flate2::write;
use flate2::read;
use flate2::{Compression};
use std::{
    error::Error,
    ffi::OsStr,
    fs::File,
    // io::{prelude::*, BufWriter, Write, BufReader, BufRead},
    io::{self, BufWriter, Write, BufReader, BufRead},
    path::Path
};


// #[derive(Debug)]

fn lines_from_file(filename: &str) -> Vec<String> {
    let path = Path::new(filename);
    let file = match File::open(&path) {
        Err(why) => panic!("couldn't open {}: {}", path.display(), why.description()),
        Ok(file) => file,
    };
    if path.extension() == Some(OsStr::new("gz")){
        let buf = BufReader::new(read::GzDecoder::new(file));
        buf.lines()
            .map(|l| l.expect("Could not parse line"))
            .collect()
    }else{
        let buf = BufReader::new(file);
        buf.lines()
            .map(|l| l.expect("Could not parse line"))
            .collect()
    }
}


// pub fn reader(filename: &str) -> Box<dyn BufRead> {
//     let path = Path::new(filename);
//     let file = match File::open(&path) {
//         Err(why) => panic!("couldn't open {}: {}", path.display(), why.description()),
//         Ok(file) => file,
//     };

//     if path.extension() == Some(OsStr::new("gz")) {
//         Box::new(BufReader::with_capacity(
//             128 * 1024,
//             read::GzDecoder::new(file),
//         ))
//     } else {
//         Box::new(BufReader::with_capacity(128 * 1024, file))
//     }
// }

// pub fn writer(filename: &str) -> BufWriter<W> {
//     let path = Path::new(filename);
//     let file = match File::create(&path) {
//         Err(why) => panic!("couldn't open {}: {}", path.display(), why.description()),
//         Ok(file) => file,
//     };

//     if path.extension() == Some(OsStr::new("gz")) {
//         // Error is here: Created file isn't gzip-compressed
//         BufWriter::new(write::GzEncoder::new(file, Compression::default()))
//     } else {
//         BufWriter::new(file)
//     }
// }


struct Params {
    fastq1: String,
    fastq2: String,
    ofastq: String,
    bcs: String,
    umi_len: u8,
    cb_len: u8,
    threads: usize,
    // max_reads: usize,
    // stop: bool,
    // debug: bool,
    name_sep: String,
}

// pub struct ReadPair<'a> {
//     rec1: RefRecord<'a>,
//     rec2: RefRecord<'a>,
//     pass_qc: Option<bool>,
//     has_n: Option<bool>,
// }

// impl ReadPair <'_> {
//     fn check_n(mut self){
//         if self.rec1.seq().contains(&b"N"[0]) | self.rec2.seq().contains(&b"N"[0]){
//             self.has_n = Some(true);
//         } else{
//             self.has_n = Some(false);
//         }
//     }
//     fn finalcheck(mut self){
//         self.check_n();
//         if self.has_n.unwrap(){
//             self.pass_qc = Some(false);
//         }else{
//             self.pass_qc = Some(true);
//         }
//     }
//     fn process(self){
//         // self.check_n();
//         self.finalcheck();
//         if self.pass_qc.unwrap(){
//             let owned_rec = RefRecord::to_owned_record(&self.rec2);
//             let mut writer = io::stdout();
//             let _ = &owned_rec.write(&mut writer);
//         }
//     }
// }

// pub struct ReadPair{
//     rec1: fastq::OwnedRecord,
//     rec2: fastq::OwnedRecord,
//     pass_qc: Option<bool>,
//     has_n: Option<bool>,
// }

// impl ReadPair {
//     fn check_n(&self){
//         if &self.rec1.seq().contains(&b"N"[0]) | &self.rec2.seq().contains(&b"N"[0]){
//             self.has_n = Some(true);
//         } else{
//             self.has_n = Some(false);
//         }
//     }
//     fn finalcheck(&self){
//         &self.check_n();
//         if self.has_n.unwrap(){
//             self.pass_qc = Some(false);
//         }else{
//             self.pass_qc = Some(true);
//         }
//     }
//     fn process(self){
//         // self.check_n();
//         // self.finalcheck();
//         if self.pass_qc.unwrap(){
//             let mut writer = io::stdout();
//             self.rec2.write(&mut writer);
//         }
//     }
// }

fn load_params() -> Params {
    // let max_r = 0usize;
    // let stop = true;
    let yaml = load_yaml!("params_fastq1.yml");
    let params = App::from_yaml(yaml).get_matches();
    // let debug = params.value_of("debug").unwrap_or("false");
    // let debug = debug.to_string().parse::<bool>().unwrap();
    let fastq1 = params.value_of("fastq1").unwrap();
    let fastq2 = params.value_of("fastq2").unwrap();
    let ofastq = params.value_of("outfastq").unwrap_or("out.fastq.gz");
    let bcs = params.value_of("barcodes_file").unwrap_or("/Users/sfurlan/develop/mutCaller/data/737K-august-2016.txt.gz");
    let threads = params.value_of("threads").unwrap_or("1");
    let threads = threads.to_string().parse::<u8>().unwrap();
    let umi_len = params.value_of("umi_len").unwrap_or("10");
    let umi_len = umi_len.to_string().parse::<u8>().unwrap();
    let cb_len = params.value_of("cb_len").unwrap_or("16");
    let cb_len = cb_len.to_string().parse::<u8>().unwrap();
    // let _max_reads = params.value_of("n_reads").unwrap_or("all");
    let name_sep = params.value_of("name_sep").unwrap_or("|BARCODE=");
    // if _max_reads == "all"{
    //     let _stop = false;
    //     let _max_r = 0usize;
    // }else{
    //     let _stop = true;
    //     let _max_r = _max_reads.to_string().parse::<usize>().unwrap();
    // }
    Params{
        fastq1: fastq1.to_string(),
        fastq2: fastq2.to_string(),
        ofastq: ofastq.to_string(),
        bcs: bcs.to_string(),
        threads: threads as usize,
        umi_len: umi_len as u8,
        cb_len: cb_len as u8,
        // max_reads: max_r,
        // stop: stop,
        // debug: debug,
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
    let mut cbvec = lines_from_file(&params.bcs);
    cbvec.sort_unstable();
    let zip = true;
    let mut total_count: usize = 0;
    let mut nfound_count: usize = 0;
    let mut mmcb_count: usize = 0;
    let split_at = &params.umi_len + &params.cb_len;
    let sep: Vec::<u8> = params.name_sep.as_bytes().to_vec();
    let fastq1 = &params.fastq1;
    let fastq2 = &params.fastq2;
    let _counts = (0u64, 0u64);
    let path = Path::new(&params.ofastq);
    let file = match File::create(&path) {
        Err(why) => panic!("couldn't open {}: {}", path.display(), why.description()),
        Ok(file) => file,
    };
    let mut writer = io::stdout();
    // let mut writer = BufWriter::new(write::GzEncoder::new(file, Compression::default()));
    // let gzext = OsStr::new("gz");
    // let stdout;
    // let filez;
    // let mut writer = match path.extension(){
    //     Some(gzext) => {
    //         filez = BufWriter::new(write::GzEncoder::new(file, Compression::default()));
    //         &filez as &Write
    //     }
    //     _ => {
    //         stdout = io::stdout();
    //         &stdout as &Write
    //     }
    // };
    // let mut writer_file = writer(&params.ofastq);
    parse_path(Some(fastq1), |parser1| {
        parse_path(Some(fastq2), |parser2| {
            each_zipped(parser1, parser2, |rec1, rec2| {
                if rec1.is_some() & rec2.is_some(){
                    let r1 = &rec1.unwrap();
                    let r2 = &rec2.unwrap();
                    if r1.seq().contains(&b"N"[0]) | r2.seq().contains(&b"N"[0]){
                        nfound_count += 1;
                        total_count +=1;
                    }else{
                        total_count +=1;
                        // let barcode = RefRecord::to_owned_record(&r1);
                        let (barcode, _seq) = &r1.seq().split_at(split_at.into());
                        // println!("{:?}", barcode);
                        let (cb, _seq) = barcode.split_at(params.cb_len as usize);
                        match cbvec.binary_search(&std::str::from_utf8(cb).unwrap().to_string()) {
                            Ok(_u) => {
                                let mut readout = RefRecord::to_owned_record(&r2);
                                let mut new_header = BytesMut::from(readout.head()).to_vec();
                                let _ = new_header.append(&mut sep.clone());
                                let _ = new_header.append(&mut barcode.to_vec());
                                readout.head = new_header;
                                // let _ = readout.write(&mut writer(&params.ofastq));
                                let _ = readout.write(&mut writer);
                            }
                            Err(_e) => {
                                mmcb_count +=1;
                            }
                        }
                        // let newhead = &mut sep.clone().append(&mut barcode.to_vec());
                        // let newrecord = <OwnedRecord as Trait>::Record {
                        //     seq: &r2.seq(),
                        //     head: &newhead,
                        //     qual: &r2.qual(),
                        // };
                    }
                }
                (true, true)
            })
            .expect("Invalid record.");
        })
        .expect("Unknown format for file 2.");
    })
    .expect("Unknown format for file 1.");
    eprintln!("Total number of reads processed: {}, {} of these had Ns, {} of these had BC not in whitelist", total_count, nfound_count, mmcb_count);
}

// pub fn eval_reads(_rec1: &RefRecord, rec2: &RefRecord) -> bool{
//     let _writer = io::stdout();
//     if rec2.seq().contains(&b"N"[0]){
//         println!("true");
//         false
//     }else{
//         println!("false");
//         true
//     }
//     // println!("{:?}",rec1.seq())
// }

// fn fastq_count(params: &Params) {
//     let _total: usize = 0;
//     let _split_at = params.umi_len + params.cb_len;
//     let _sep: Vec::<u8> = "|BARCODE=".as_bytes().to_vec();
//     let fastq1 = &params.fastq1;
//     let fastq2 = &params.fastq2;
//     let mut counts = (0u64, 0u64);
//     parse_path(Some(fastq1), |parser1| {
//         parse_path(Some(fastq2), |parser2| {
//             each_zipped(parser1, parser2, |rec1, rec2| {
//                 if rec1.is_some() {
//                     counts.0 += 1;
//                 }
//                 if rec2.is_some() {
//                     counts.1 += 1;
//                 }
//                 (true, true)
//             })
//             .expect("Invalid record.");
//         })
//         .expect("Unknown format for file 2.");
//     })
//     .expect("Unknown format for file 1.");

//     println!("Number of reads: ({}, {})", counts.0, counts.1);
// }



// fn modify_fastq (record: RefRecord, split_at: u8, namesep: &Vec<u8>) -> bool {
//     let sep: &mut Vec::<u8> = &mut namesep.clone();
//     let mut writer = io::stdout();
//     let mut owned_rec = RefRecord::to_owned_record(&record);
//     let curr_bytes = BytesMut::from(owned_rec.qual());
//     let (_barcode, seq) = &curr_bytes.split_at(split_at.into());
//     owned_rec.qual = seq.to_vec();
//     let curr_bytes = BytesMut::from(owned_rec.seq());
//     let (barcode, seq) = &curr_bytes.split_at(split_at.into());
//     owned_rec.seq = seq.to_vec();
//     let mut curr_bytes = BytesMut::from(owned_rec.head()).to_vec();
//     // let _ = &curr_bytes.push(namesep);
//     let _ = &curr_bytes.append(sep);
//     let _ = &curr_bytes.append(&mut barcode.to_vec());
//     owned_rec.head = curr_bytes;
//     owned_rec.sep = Some(vec!['+' as u8]);
//     let _ = &owned_rec.write(&mut writer);
//     true
// }


// fn modify_fastq_debug(record: RefRecord, split_at: u8, namesep: &Vec<u8>) -> bool {
//     let sep: &mut Vec::<u8> = &mut namesep.clone();
//     println!("Before mod");
//     println!("{}", String::from_utf8_lossy(record.head()));
//     println!("{}", String::from_utf8_lossy(record.seq()));
//     println!("{}", String::from_utf8_lossy(record.qual()));
//     let mut owned_rec = RefRecord::to_owned_record(&record);
//     let curr_bytes = BytesMut::from(owned_rec.qual());
//     let (_barcode, seq) = &curr_bytes.split_at(split_at.into());
//     owned_rec.qual = seq.to_vec();
//     let curr_bytes = BytesMut::from(owned_rec.seq());
//     let (barcode, seq) = &curr_bytes.split_at(split_at.into());
//     owned_rec.seq = seq.to_vec();
//     let mut curr_bytes = BytesMut::from(owned_rec.head()).to_vec();
//     let _ = &curr_bytes.append(sep);
//     let _ = &curr_bytes.append(&mut barcode.to_vec());
//     owned_rec.head = curr_bytes;
//     owned_rec.sep = Some(vec!['+' as u8]);
//     println!("After mod");
//     println!("{}", String::from_utf8_lossy(owned_rec.head()));
//     println!("{}", String::from_utf8_lossy(owned_rec.seq()));
//     println!("{}", String::from_utf8_lossy(owned_rec.qual()));
//     true
// }


// use flate2::read;
// use flate2::write;
// use flate2::Compression;
// use std::error::Error;
// use std::ffi::OsStr;
// use std::fs::File;
// use std::io::{self, BufRead, BufReader, BufWriter, Write};
// use std::path::Path;

// /// Read normal or compressed files seamlessly
// /// Uses the presence of a `.gz` extension to decide
// pub fn reader(filename: &str) -> Box<dyn BufRead> {
//     let path = Path::new(filename);
//     let file = match File::open(&path) {
//         Err(why) => panic!("couldn't open {}: {}", path.display(), why.description()),
//         Ok(file) => file,
//     };

//     if path.extension() == Some(OsStr::new("gz")) {
//         Box::new(BufReader::with_capacity(
//             128 * 1024,
//             read::GzDecoder::new(file),
//         ))
//     } else {
//         Box::new(BufReader::with_capacity(128 * 1024, file))
//     }
// }

// Attempting to have a file writer too
// pub fn writer(filename: &str) -> Box<dyn Write> {
//     let path = Path::new(filename);
//     let file = match File::create(&path) {
//         Err(why) => panic!("couldn't open {}: {}", path.display(), why.description()),
//         Ok(file) => file,
//     };

//     if path.extension() == Some(OsStr::new("gz")) {
//         // Error is here: Created file isn't gzip-compressed
//         Box::new(BufWriter::with_capacity(
//             128 * 1024,
//             write::GzEncoder::new(file, Compression::default()),
//         ))
//     } else {
//         Box::new(BufWriter::with_capacity(128 * 1024, file))
//     }
// }

// /// Doing tests
// fn main() -> io::Result<()> {
//     // Test with uncompressed file
//     let filename = "file.txt";
//     println!("Testing reader with uncompressed file: '{}'", filename);
//     let reader_file = reader(filename);
//     for line in reader_file.lines() {
//         println!("{}", line?);
//     }
//     println!();

//     // Test with compressed file
//     let filename = "file.txt.gz";
//     println!("Testing reader with compressed file: '{}'", filename);
//     let reader_file_gz = reader(filename);
//     for line in reader_file_gz.lines() {
//         println!("{}", line?);
//     }
//     println!();

//     // Test writing to uncompressed file
//     let filename = "file.output.txt";
//     println!("Testing writer with compressed file: '{}'", filename);
//     let mut writer_file = writer(filename);
//     for _i in 1..=100 {
//         writer_file.write_all(b"This is the end. Count your chickens.\n")?;
//     }

//     // Test writing to compressed file
//     let filename = "file.output.txt.gz";
//     println!("Testing writer with compressed file: '{}'", filename);
//     let mut writer_file = writer(filename);
//     for _i in 1..=1000 {
//         writer_file.write_all(b"This is the end. Count your chickens.\n")?;
//     }

//     Ok(())
// }

