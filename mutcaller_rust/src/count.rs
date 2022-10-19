/**

chmod 777 /Users/sfurlan/develop/mutCaller/mutcaller_rust/target/debug/mutcaller_rust
cd /Users/sfurlan/develop/mutCaller/mutcaller_rust
cargo build
/Users/sfurlan/develop/mutCaller/mutcaller_rust/target/debug/mutcaller_rust -t 18 --ibam "/Users/sfurlan/Fred Hutchinson Cancer Research Center/Furlan_Lab - General/experiments/patient_marrows/LKmut/Aligned.out.tagged.sorted_0.05.bam"

#cluster
cd /home/sfurlan/develop/mutCaller/mutcaller_rust
cargo build
~/develop/mutCaller/mutcaller_rust/target/debug/count --ibam=$OUT/kquant/pseudoalignments.bam | head

###test dataset.
##make test
samtools view -hb -s 0.01 $OUT/kquant/pseudoalignments.bam > /home/sfurlan/develop/mutCaller/data/bams/test.bam

###test dataset.
##make test
samtools view -hb -s 0.01 $OUT/kquant/pseudoalignments.bam > /home/sfurlan/develop/mutCaller/data/bams/test.bam

##run test

~/develop/mutCaller/mutcaller_rust/target/release/count --ibam=/home/sfurlan/develop/mutCaller/data/bams/test.bam | head
time ~/develop/mutCaller/addTags.py -u 10 -c 16 /home/sfurlan/develop/mutCaller/data/bams/test.bam | samtools view -hbo kquant/Aligned.out.tagged.bam > test.bam
time ~/develop/mutCaller/mutcaller_rust/target/release/count -t 24 --ibam=/home/sfurlan/develop/mutCaller/data/bams/test.bam

**/

extern crate clap;
extern crate bam;

use std::io;
// use std::io::{self, Write};
// use bam::RecordWriter;
use clap::{App, load_yaml};
use std::str;
// use bam::header::{Header, HeaderEntry};





// #[derive(Clone)]

struct Params {
    ibam: String,
    split: String,
    joiner: String,
    threads: usize,
    cb_len: usize, 
    // umi_len: usize,
    read_len: usize,
}

fn load_params() -> Params {
    let yaml = load_yaml!("params.yml");
    let params = App::from_yaml(yaml).get_matches();
    let ibam = params.value_of("input_bam").unwrap();
    eprintln!("opening: {}", ibam.to_string());
    let joiner = params.value_of("joiner").unwrap_or(":");
    let split = params.value_of("split").unwrap_or("|BARCODE=");
    let cb_len = params.value_of("cb_len").unwrap_or("16");
    let cb_len = cb_len.to_string().parse::<usize>().unwrap() - 1;
    let read_len = params.value_of("read_len").unwrap_or("90");
    let read_len = read_len.to_string().parse::<usize>().unwrap();
    // let umi_len = params.value_of("umi_len").unwrap_or("10");
    // let umi_len = umi_len.to_string().parse::<usize>().unwrap() - 1;
    let threads = params.value_of("threads").unwrap_or("1");
    let threads = threads.to_string().parse::<usize>().unwrap() - 1;
    Params{
        ibam: ibam.to_string(),
        threads: threads,
        split: split.to_string(),
        joiner: joiner.to_string(),
        cb_len: cb_len,
        // umi_len: umi_len,
        read_len: read_len,
    }
}

fn main() {
    let params = load_params();
    count(&params);
}


fn count(params: &Params) {
    let mut total: usize = 0;
    let mut goodreadcount: usize = 0;
    let (_read_threads, _write_threads) = if (*&params.threads as i8) > 2{
        (((*&params.threads/2) -1) as u16, ((*&params.threads/2) -1) as u16)
    } else {
        (0 as u16, 0 as u16)
    };
    let reader = bam::BamReader::from_path(params.ibam.to_string(), 0).unwrap();
    let _output = io::BufWriter::new(io::stdout());
    let header = reader.header().clone();
    let data = header.reference_names();
    let mut seqnames = Vec::new();
    for seq in data {
        // write!(f, "{}\n", seq).expect("unable to write");
        seqnames.push(seq)
    }
    // let array: [T; N]  = data.[into_]iter()
    //     .collect::<Vec<T>>()
    //     .try_into()
    //     .unwrap()
    // let mut oldseqname = "None";
    for record in reader {
        total += 1;
        let newrecord = record.unwrap();
        let seqname = match str::from_utf8(&newrecord.name()) {
            Ok(v) => v,
            Err(e) => panic!("Invalid UTF-8 sequence: {}", e),
        };
        // if &seqname == &oldseqname {
        //     continue;
        // }
        let cbumi= seqname.split(&params.split).nth(1).unwrap().to_string();
        let _modified_name = seqname.replace(&params.split, &params.joiner);
        let (cb_umi_s1, cb_umi_s2) = cbumi.split_at(params.cb_len+1);
        let mut good_read = false;
        let cigarmatch = format!("{}M", *&params.read_len);
        let cigar = newrecord.cigar().to_string();
        if cigar == cigarmatch{
            good_read = true
        }
        if newrecord.mapq() < 255 as u8 {
            good_read = false
        }
        if good_read{
            goodreadcount += 1;
            // oldseqname = seqname;
            println!("{}\t{}\t{}", cb_umi_s1, cb_umi_s2, seqnames[newrecord.ref_id() as usize].to_string());
            // println!("{}\t{}\t{}\t{}", cb_umi_s1, cb_umi_s2, seqnames[newrecord.ref_id() as usize].to_string(), newrecord.start().to_string());
            // println!("{}\t{}\t{}\t1", cb_umi_s1, cb_umi_s2, newrecord.ref_id().to_string());
            // println!("{}\t{}\t{}\t1\t(CIGAR: {})", cb_umi_s1, cb_umi_s2, newrecord.ref_id().to_string(), cigar);
        }
        // let total_str =  &total.to_string().chars();
        // if total_str.count() > 7 as usize{
        //     let str_temp = total_str.rev().take(6).to_string();
        //     if str_temp == "000000" {
        //         eprintln!("Processed {} total reads so far", &total);
        //     }
        // }
    }
    eprintln!("Completed; {} total reads processed!", &total);
    eprintln!("{} good reads counted!", &goodreadcount);
}
