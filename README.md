<img width="150" alt="image" src="mutcaller.png">


#                       mutCaller
Pipeline for extracting variants from single cell genomics data


### Overview

mutCaller is a command line tool written in Rust for identifying SNVs in single cell data that contains a cell barcode (CB) and UMI such as 10X genomics data files.  With mutCaller, you can supply variants (VCF file support coming), and obtain a counts file of the number of variants that map to your query and their associated cell-barcodes and umis.  mutCaller using the takes as input either unaligned fastqs (UNALIGNED function) or an aligned BAM file with the CB and UMI as tags (ALIGNED function).

### Installation

mutCaller is written in Rust.  To install the rust compiler go to https://www.rust-lang.org/tools/install.  mutCaller requires two additional tools be available on the command line, minimap2 (https://github.com/lh3/minimap2) and samtools (https://samtools.github.io). 

To install mutCaller:
1. clone the repository by typing `https://github.com/furlan-lab/mutCaller.git` from the location you want to build from
2. enter the cloned repo by typing `cd mutCaller`
3. build by typing `cargo build --release`
4. the build process will create a self contained binary executable file in `targets/release` directory called `mutCaller`
5. move this binary elsewhere if desired (ideally somewhere referenced by your PATH environment variable - e.g. `~/.local/bin`)

### Updates

**version 0.22** - 5/9/23 - first alpha release


### Usage

##### Variants File.

First users will generate a simple, headerless 'variants_file' that lists the variants to be counted (tsv.file). The alleles_file should look something like this:

```sh
seqname\tstart\tref_nt\tquery_nt\tname
```
**More detailed explanation:**
1. seqname - e.g. 'chr1', 'chr2', etc
2. position - 1-indexed position (e.g, '112450407')
3. ref_nt - nucleotide of the reference at this location
4. query_nt - nucleotide of your query
5. name - a string given to the name the variant in the outpull file

**three full lines should look something like this**

```sh
chr12   112450407   A   G   PTPN11_227A>G
chr12   208248389   G   A   IDH1_132G>A
chr17   7674220 C   T   TP53_248C>T
```

##### Barcodes file

A barcode whitelist.  v3 3prime and v2 5prime 10X whitelists are provided in the data folder.


##### Invocation

**To run mutCaller simply type:**

```sh
mutCaller ALIGNED --bam <file.bam> -a <variants.tsv> -o <folder>
mutCaller UNALIGNED --barcodes_file <barcodes_file> --fastq1 <fastq1> --fastq2 <fastq2> --genome <genome> --variants <variants.tsv>
```


##### Help menu

```sh
mutCaller 0.3
copyright Scott Furlan
mutCaller is a command line tool for aligning and counting long-read sequence specific for HLA alleles in single cell
libraries

USAGE:
    mutCaller [FLAGS] [OPTIONS] --alleles <alleles_file> --bam <bam>

FLAGS:
    -h, --help       Prints help information
    -s, --ret_seq    include this flag to return sequence in molecule_info text file
    -V, --version    Prints version information
    -v, --verbose    verbose

OPTIONS:
    -l, --level <align_level>       align to 'genome', 'transcriptome', or 'both'; default is 'both'
    -a, --alleles <alleles_file>    table of hla genes to search for (tsv file)
    -b, --bam <bam>                 input bam
    -c, --cb_tag <cb>               character to parse cell barcode; default = 'CB'
    -o, --out <output_folder>       folder for output; default 'out'
    -t, --threads <threads>         threads
    -u, --umi_tag <umi>             character to parse umi; default = 'XM'
```
 
### Output
 
mutCaller will output some combination of the following files (depending on alignment level):
```sh
Aligned_mm2_sorted_gene.bam          = minimap2 output bam file aligned to genome, sorted by read name 
Aligned_mm2_sorted_mRNA.bam          = minimap2 output bam file aligned to transcriptome, sorted by read name 
Aligned_mm2_sorted_gene.bam.bai      = index for above bam
Aligned_mm2_sorted_mRNA.bam.bai      = index for above bam
counts_gene.txt.gz                   = counts file; columns: CB, UMI, allele, read_count
counts._mRNAtxt.gz                   = counts file; columns: CB, UMI, allele, read_count
align_gene.fa                        = fasta reference file used for minimap2 alignment
align_mRNA.fa                        = fasta reference file used for minimap2 alignment
molecules_info_gene.txt.gz           = a file listing alignment metrics for each molecule; columns: CB, nb_tag (if present), UMI, allele, start_pos, mapq, cigar, NM, AS, s1, de); optionally including sequence in the last column
molecules_info_mRNA.txt.gz           = a file listing alignment metrics for each molecule; columns: CB, nb_tag (if present), UMI, allele, start_pos, mapq, cigar, NM, AS, s1, de); optionally including sequence in the last column
```
See minimap2 manual (https://lh3.github.io/minimap2/minimap2.html) for a discussion of the molecule_info metrics

See Isoseq FAQ for a discussion of the nb tag (https://isoseq.how/umi/isoseq-correct.html)










