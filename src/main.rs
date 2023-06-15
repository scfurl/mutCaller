
use clap::{App, load_yaml};


pub mod mutcaller;
pub mod countbam;
pub mod vcf;

use crate::mutcaller::mutcaller_run;
use crate::countbam::{countbam_run};
use crate::vcf::{read_vcf_compressed, read_vcf_uncompressed, guess_vcf, guess_compression};
use crate::mutcaller::read_csv_str;

fn main() {
    let yaml = load_yaml!("cli.yml");
    let params = App::from_yaml(yaml).get_matches();
    if let Some(_params) = params.subcommand_matches("UNALIGNED") {
    	mutcaller_run();
	}
	if let Some(_params) = params.subcommand_matches("ALIGNED") {
    	countbam_run()
	}
    if let Some(params) = params.subcommand_matches("VCF") {
        let qual = params.value_of("qual").unwrap_or("95.0");
        let mut verbose = true;
        if params.is_present("quiet") {
                verbose = false
        };
        let qual_p = qual.parse::<f64>();
        let vcf_file = params.value_of("variants").unwrap_or("/Users/sfurlan/develop/mutCaller/tests/var.vcf.gz").to_string();
        if verbose {
            eprintln!("Reading variants file: {}\n", &vcf_file);
        }
        let is_vcf = guess_vcf(&vcf_file);
        let is_compressed = guess_compression(&vcf_file);
        let data = {
            if is_vcf.unwrap(){
                if is_compressed.unwrap(){
                    read_vcf_compressed(&vcf_file.to_string(), &qual_p.unwrap(), &verbose).unwrap() 
               } else {
                    read_vcf_uncompressed(&vcf_file.to_string(), &qual_p.unwrap(), &verbose).unwrap() 
               }
            } else {
                read_csv_str(vcf_file, true).unwrap()
            }
        };
        eprintln!("Records:\n");
        for variant in data {
            eprintln!("{}", variant);
        }
    }

}
