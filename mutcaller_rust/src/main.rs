/**

chmod 777 /Users/sfurlan/develop/mutCaller/mutcaller_rust/target/debug/mutcaller_rust
cd /Users/sfurlan/develop/mutCaller/mutcaller_rust
cargo build
/Users/sfurlan/develop/mutCaller/mutcaller_rust/target/debug/mutcaller_rust -t 18 --ibam "/Users/sfurlan/Fred Hutchinson Cancer Research Center/Furlan_Lab - General/experiments/patient_marrows/LKmut/Aligned.out.tagged.sorted_0.05.bam"

#cluster
cd /home/sfurlan/develop/mutCaller/mutcaller_rust
cargo build
/home/sfurlan/develop/mutCaller/mutcaller_rust/target/debug/mutcaller_rust --joiner="_" --no_header --ibam=$OUT/kquant/pseudoalignments.bam | head -n 1000000
**/

extern crate clap;
extern crate bam;

use std::io;
use bam::RecordWriter;
use clap::{App, load_yaml};
use std::str;
// use bam::header::{Header, HeaderEntry};

#[derive(Clone)]

struct Params {
    ibam: String,
    split: String,
    joiner: String,
    threads: usize,
    include_header: bool,
    cb_len: usize, 
    umi_len: usize,
}

fn load_params() -> Params {
    let yaml = load_yaml!("params.yml");
    let params = App::from_yaml(yaml).get_matches();
    let ibam = params.value_of("input_bam").unwrap();
    let joiner = params.value_of("joiner").unwrap_or(":");
    let include_header = params.value_of("no_header").is_some();
    let split = params.value_of("split").unwrap_or("|BARCODE=");
    let cb_len = params.value_of("cb_len").unwrap_or("16");
    let cb_len = cb_len.to_string().parse::<usize>().unwrap() - 1;
    let umi_len = params.value_of("umi_len").unwrap_or("10");
    let umi_len = umi_len.to_string().parse::<usize>().unwrap() - 1;
    let threads = params.value_of("threads").unwrap_or("1");
    let threads = threads.to_string().parse::<usize>().unwrap() - 1;
    Params{
        ibam: ibam.to_string(),
        threads: threads,
        split: split.to_string(),
        joiner: joiner.to_string(),
        include_header: include_header,
        cb_len: cb_len,
        umi_len: umi_len,
    }
}

fn main() {
    let params = load_params();
    mutcaller(&params);
}


fn mutcaller(params: &Params) {
    // println!("Reading bam file : {}",params.ibam.to_string());
    let reader = bam::BamReader::from_path(params.ibam.to_string(), params.threads as u16).unwrap();
    let output = io::BufWriter::new(io::stdout());
    let mut writer = bam::SamWriter::build()
        .write_header(params.include_header)
        .from_stream(output, reader.header().clone()).unwrap();

    for record in reader {
        let record = record.unwrap();
        let seqname = match str::from_utf8(&record.name()) {
            Ok(v) => v,
            Err(e) => panic!("Invalid UTF-8 sequence: {}", e),
        };
        let cbumi= seqname.split(&params.split).nth(1).unwrap().to_string();
        let modified_name = seqname.replace(&params.split, &params.joiner);
        let (cb_umi_s1, cb_umi_s2) = cbumi.split_at(params.cb_len+1);
        let mut newrecord = record.clone();
        newrecord.tags_mut().push_string(b"CB", &cb_umi_s1.as_bytes());
        newrecord.tags_mut().push_string(b"CR", &cb_umi_s1.as_bytes());
        newrecord.tags_mut().push_string(b"UB", &cb_umi_s2.as_bytes());
        newrecord.tags_mut().push_string(b"UR", &cb_umi_s2.as_bytes());
        newrecord.set_name(modified_name.bytes());
        writer.write(&newrecord).unwrap();
    }
}
