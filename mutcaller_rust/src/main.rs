/**

chmod 777 /Users/sfurlan/develop/mutCaller/mutcaller_rust/target/debug/mutcaller_rust
cd /Users/sfurlan/develop/mutCaller/mutcaller_rust
cargo build
/Users/sfurlan/develop/mutCaller/mutcaller_rust/target/debug/mutcaller_rust -t 18 --ibam "/Users/sfurlan/Fred Hutchinson Cancer Research Center/Furlan_Lab - General/experiments/patient_marrows/LKmut/Aligned.out.tagged.sorted_0.05.bam"
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
    threads: usize,
}

fn load_params() -> Params {
    let yaml = load_yaml!("params.yml");
    let params = App::from_yaml(yaml).get_matches();
    let ibam = params.value_of("input_bam").unwrap();
    let split = params.value_of("split").unwrap_or("|BARCODE=");
    let threads = params.value_of("threads").unwrap_or("1");
    let threads = threads.to_string().parse::<usize>().unwrap() - 1;
    Params{
        ibam: ibam.to_string(),
        threads: threads,
        split: split.to_string(),
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
        .write_header(true)
        .from_stream(output, reader.header().clone()).unwrap();

    for record in reader {
        let record = record.unwrap();
        // writer.write(&record.name).unwrap();
        let seqname = match str::from_utf8(&record.name()) {
            Ok(v) => v,
            Err(e) => panic!("Invalid UTF-8 sequence: {}", e),
        };
        let mut ss= seqname.split(&params.split);
        // println!("{:?}", ss.nth(1));
        let mut newrecord = record.clone();
        newrecord.set_name(ss.nth(1).unwrap().bytes());
        // for s in vec {
        //     println!("{}", ss)
        // }
        writer.write(&newrecord).unwrap();
    }
}


// fn mutcaller(params: &Params) {
//     // Creating a header.
//     let mut header = Header::new();
//     println!("Reading bam file : {}",params.ibam.to_string());
//     // Header line          "@HD  VN:1.6  SO:Coordinate".
//     let mut header_line = HeaderEntry::header_line("1.6".to_string());
//     header_line.push(b"SO", "Coordinate".to_string());
//     header.push_entry(header_line).unwrap();
//     // Reference line       "@SQ  SN:chr1  LN:10000".
//     header.push_entry(HeaderEntry::ref_sequence("chr1".to_string(), 10000)).unwrap();
 
//     // Write SAM to stdout.
//     let output = io::BufWriter::new(io::stdout());
//     let mut writer = bam::SamWriter::from_stream(output, header).unwrap();
 
//     // Create a new record, set its name to "Read_1",
//     // reference id to 0, start to 10 (both 0-based).
//     let mut record = bam::Record::new();
//     record.set_name("Read_1".bytes());
//     record.set_ref_id(0);
//     record.set_start(10);
//     // Set the record to be on reverse strand.
//     record.flag_mut().set_strand(false);
//     // Set sequence and qualities (qualities without +33), and cigar.
//     record.set_seq_qual("ACGT".bytes(), [10_u8, 20, 30, 10].iter().cloned()).unwrap();
//     record.set_cigar("2M1I1M".bytes()).unwrap();
//     // Add NM tag.
//     record.tags_mut().push_num(b"NM", 1);
 
//     writer.write(&record).unwrap();
//     writer.finish().unwrap();
//     // Above code would print the following SAM file:
//     // @HD VN:1.6  SO:Coordinate
//     // @SQ SN:chr1 LN:10000
//     // Read_1  16  chr1    11  0   2M1I1M  *   0   0   ACGT    +5?+    NM:i:1
 
//     println!("Aligned pairs:");
//     for (read_pos, ref_pos) in record.aligned_pairs() {
//         println!("    {:?} {:?}", read_pos, ref_pos);
//     }
//     // Aligned pairs:
//     //     Some(0) Some(10)
//     //     Some(1) Some(11)
//     //     Some(2) None
//     //     Some(3) Some(12)
// }