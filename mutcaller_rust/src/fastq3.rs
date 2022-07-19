
// /Users/sfurlan/develop/mutCaller/mutcaller_rust/target/debug/fastq3 




use bytes::BytesMut;
use fastq::{Parser, Record, RefRecord};

const READS: &str = r#"@read1/ENST00000266263.10;mate1:84-183;mate2:264-363
GACAGCCAGGGGCCAGCGGGTGGCAGTGCCCAGGACATAGAGAGAGGCAGCACACACGCGGTTGATGGTGAAGCCCGGAATGGCCACAGAGGCTAGAGCC
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
@read2/ENST00000266263.10;mate1:163-262;mate2:283-382
GATGCCATTGACAAAGGCAAGAAGGCTGGAGAGGTGCCCAGCCCTGAAGCAGGCCGCAGCGCCAGGGTGACTGTGGCTGTGGTGGACACCTTTGTATGGC
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
@read3/ENST00000266263.10;mate1:86-185;mate2:265-364
GGACAGCCAGGGGCCAGCGGGTGGCAGTGCCCAGGACATAGAGAGAGGCAGCANACACACGGTTGATGGTGAAGCCCGGAATGGCCACAGAGGCTAGAGC
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII!IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
@read4/ENST00000266263.10;mate1:297-396;mate2:401-500
CAGGAGGAGCTGGGCTTCCCCACTGTTAGGTAGAGCTTGCGCAGGCTGGAGTCCAGGAGGAAATCCACCGACCTGTCAATGGGGTGGATAATGATGGGGA
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
"#;

fn main() {
    let mut total: usize = 0;

    let parser = Parser::new(READS.as_bytes());
    parser
        .each(|record| {
            // modify_qual1(record, &mut total)
            println!("Before mod");
            println!("{}", String::from_utf8_lossy(record.qual()));
            let owned_rec = modify_qual(record);
            println!("After mod");
            println!("{}", String::from_utf8_lossy(owned_rec.qual()));
            total += 1;
            true
        })
        .expect("Invalid fastq file");
    println!("{}", total);
}

fn modify_qual(record: RefRecord) -> fastq::OwnedRecord {
    let mut owned_rec = RefRecord::to_owned_record(&record);
    let mut curr_bytes = BytesMut::from(owned_rec.qual());
    curr_bytes[0] = b'$';
    owned_rec.qual = curr_bytes.to_vec();
    owned_rec
}

fn modify_qual1(record: RefRecord, total: &mut usize) -> bool {
    println!("Before mod");
    println!("{}", String::from_utf8_lossy(record.qual()));
    let mut owned_rec = RefRecord::to_owned_record(&record);
    let mut curr_bytes = BytesMut::from(owned_rec.qual());
    curr_bytes[0] = b'$';
    owned_rec.qual = curr_bytes.to_vec();
    println!("After mod");
    println!("{}", String::from_utf8_lossy(owned_rec.qual()));
    *total += 1;
    true
}