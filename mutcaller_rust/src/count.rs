/**

cd ~/develop/mutCaller/mutcaller_rust/tests
time ~/develop/mutCaller/mutcaller_rust/target/release/count -t 24 --ibam=kquant/pseudoalignments.bam -a kallisto > counts_k.txt


time ~/develop/mutCaller/mutcaller_rust/target/release/count -t 24 --ibam=Aligned.sortedByCoord.out.tagged.bam > counts_s.txt

**/

extern crate clap;
extern crate bam;

use std::io;
use clap::{App, load_yaml};
use std::str;





// #[derive(Clone)]

struct Params {
    ibam: String,
    aligner: String,
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
    let aligner = params.value_of("aligner").unwrap_or("STAR");
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
        aligner: aligner.to_string(),
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
    if params.aligner == "kallisto"{
        count_kallisto(&params);
        return;
    }
    if params.aligner == "STAR"{
        count_star(&params);
    }
    
}


fn count_star(params: &Params) {
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
        seqnames.push(seq)
    }
    for record in reader {
        total += 1;
        let newrecord = record.unwrap();
        let seqname = match str::from_utf8(&newrecord.name()) {
            Ok(v) => v,
            Err(e) => panic!("Invalid UTF-8 sequence: {}", e),
        };
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
        if good_read && ((newrecord.flag().to_string()=="Flag(16)") | (newrecord.flag().to_string()=="Flag(0)")){
            goodreadcount += 1;
            println!("{} {} {} {}", cb_umi_s1, cb_umi_s2, seqnames[newrecord.ref_id() as usize].to_string(), newrecord.start());
        }
    }
    eprintln!("Completed; {} total reads processed!", &total);
    eprintln!("{} good reads counted!", &goodreadcount);
}



fn count_kallisto(params: &Params) {
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
        seqnames.push(seq)
    }
    for record in reader {
        total += 1;
        let newrecord = record.unwrap();
        let seqname = match str::from_utf8(&newrecord.name()) {
            Ok(v) => v,
            Err(e) => panic!("Invalid UTF-8 sequence: {}", e),
        };
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
        if good_read && ((newrecord.flag().to_string()=="Flag(16)") | (newrecord.flag().to_string()=="Flag(0)")){
            goodreadcount += 1;
            println!("{} {} {}", cb_umi_s1, cb_umi_s2, seqnames[newrecord.ref_id() as usize].to_string());
        }
    }
    eprintln!("Completed; {} total alignments processed!", &total);
    eprintln!("{} good alignments counted!", &goodreadcount);
}
