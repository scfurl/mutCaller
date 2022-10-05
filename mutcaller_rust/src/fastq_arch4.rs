/*


/Users/sfurlan/develop/mutCaller/mutcaller_rust/target/debug/fastq3 --ifastq "/Users/sfurlan/Fred Hutchinson Cancer Research Center/Furlan_Lab - General/experiments/patient_marrows/LKmut/sample_filtered.fastq_1M.fastq.gz"

*/


extern crate clap;
extern crate fastq;
use std::io;
use bio::io::fastq;
use bio::io::fastq::FastqRead;

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


let mut reader = fastq::Reader::new(io::stdin());
let mut writer = fastq::Writer::new(io::stdout());
let mut record = fastq::Record::new();

while let Ok(()) = reader.read(&mut record) {
    if record.is_empty() {
        let check = record.check();
        break;
    }

    let mut sum_qual = record.qual().iter().sum::<u8>() as f64;

    if (sum_qual / record.seq().len() as f64 - 33.0) > 30.0 {
        writer.write_record(&record);
    }
}
}

