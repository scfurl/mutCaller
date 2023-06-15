<img width="200" alt="image" src="mutcaller.png">


#                       mutCaller working examples

### Build mutCaller and install binary

```sh
loc=~/develop/mutCaller # or location where you have cloned the repository
cd $loc  
cargo build --release && cp target/release/mutcaller ~/.local/bin
```

### Run UNALIGNED on short read fastqs using mm2

```sh
loc=~/develop/mutCaller # or location where you have cloned the repository
bc=$loc/data/737K-august-2016.txt.gz  #barcode whitelist
fa=/Users/sfurlan/refs/genome.fa #genome location i.e. GRCh38
#fa=/fh/fast/furlan_s/grp/refs/GRCh38/refdata-gex-GRCh38-2020-A/fasta/genome.fa
mutcaller UNALIGNED -t 8 -g $fa -b $bc -s $loc/tests/variants.tsv -o out_mm2 \
          -i $loc/tests/sequencer_R1.fastq.gz \
          -j $loc/tests/sequencer_R2.fastq.gz
```


### Run UNALIGNED on short read fastqs using STAR (better performance on short reads than mm2 for unclear reasons)

```sh
loc=~/develop/mutCaller # or location where you have cloned the repository
bc=$loc/data/737K-august-2016.txt.gz  #barcode whitelist
fa=/fh/fast/furlan_s/grp/refs/GRCh38/refdata-gex-GRCh38-2020-A/star
ml SAMtools/1.11-GCC-10.2.0 #make sure samtools is accessible
mutcaller UNALIGNED \
                        -t 8 -g $fa -b $bc -s $loc/tests/variants.tsv -a STAR -l /app/software/CellRanger/6.0.1/lib/bin/STAR \
                        -o out_star -i $loc/tests/sequencer_R1.fastq.gz \
                        -j $loc/tests/sequencer_R2.fastq.gz
```


### Compare to cbsniffer

```sh
loc=~/develop/mutCaller # or location where you have cloned the repository
module load SAMtools/1.16.1-GCC-11.2.0 # make samtools accessible on command line
ml R/4.1.2-foss-2021b # make R accessible on command line
cd $loc
git checkout bash
cd mutcaller_rust
git status
cargo build --release
bc=$loc/data/737K-august-2016.txt.gz  #barcode whitelist
out=$loc/cbsniffer
mkdir $out
cd $out
$loc/mutcaller -u -U 10 -b $loc/tests/sequencer_R1.fastq.gz -t $loc/tests/sequencer_R2.fastq.gz -l $bc &&
export fq2=$out/fastq_processed/sample_filtered.fastq &&
export fq3=$out/fastq_processed/sample_filtered_header.fastq &&
$loc/mutcaller_rust/target/release/fastq -t 1 --ifastq ${fq2} > ${fq3}
export transcriptome=/shared/biodata/ngs/Reference/10X/refdata-gex-GRCh38-2020-A #location to cellranger friendly reference
/app/software/CellRanger/6.0.1/lib/bin/STAR --genomeDir $transcriptome/star --readFilesIn ${fq3} --readNameSeparator space --runThreadN 24 --outSAMunmapped Within KeepPairs --outSAMtype BAM SortedByCoordinate &&
#Rscript $loc/scripts/quantReads.R &&
$loc/scripts/addTag.py -u 10 -c 16 Aligned.sortedByCoord.out.bam | samtools view -hbo Aligned.out.tagged.sorted.bam &&
samtools index -@ 24 Aligned.out.tagged.sorted.bam
rm Aligned.sortedByCoord.out.bam &&
rm fastq.log Log.out Log.progress.out SJ.out.tab &&
rm -R fastq_processed &&
rm -R mutcaller &&
rm -R _STARtmp
zcat $bc > $loc/data/737K-august-2016.txt
sed -i 's/$/-1/g' $loc/data/737K-august-2016.txt
cat $loc/data/737K-august-2016.txt | head
python $loc/scripts/cb_sniffer.py Aligned.out.tagged.sorted.bam $loc/tests/variants_cb_sniffer.tsv $loc/data/737K-august-2016.txt test

```

### Check a variants.tsv file to make sure it is compatible (zipped and unzipped ok)
```sh
loc=~/develop/mutCaller 
mutcaller VCF -s $loc/tests/variants.tsv
mutcaller VCF -s $loc/tests/variants.tsv.gz
```

### Look at variants file (VCF format) after filtering out quality scores below 98 (zipped and unzipped ok)
```sh
loc=~/develop/mutCaller 
cargo build --release && cp ~/develop/mutCaller/target/release/mutcaller ~/.local/bin
mutcaller VCF -s $loc/tests/var.vcf.gz -q 98
```


### Run ALIGNED on a bam file using VCF file filtered on 98
```sh
loc=~/develop/mutCaller
mutcaller ALIGNED -b $loc/tests/lr.bam -s $loc/tests/var.vcf.gz -q 98 -t 1 -o out_long
```

# THANKS FOR TRYING mutCaller!!






