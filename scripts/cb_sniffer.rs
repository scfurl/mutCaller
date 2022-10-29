use std::*;
use std::collections::HashMap;

use collections::{defaultdict};
struct GenomicPosition {
chrm: ST0,
start: ST1,
end: ST2,
ref: ST3,
alt: ST4,
gene: ST5,
event: ST6,
}

impl GenomicPosition {
fn __init__<T0>(&self, line: T0)  {
"
		:param line:
		 for each variant holds genomics positions
		 and other informations
		";
let region = line.strip().split("	");
if region.len() < 7 {
logger.debug("Insufficient columns fields please check variant file");
sys.exit(1);
}
self.chrm = region[0];
self.start = i32::from(region[1]);
self.end = i32::from(region[2]);
self.ref = region[3];
self.alt = region[4];
self.gene = region[5];
self.event = region[6];
}
fn classify<RT>(&self) -> RT {
"
		for each variant classifies if its a snp or variant
		:return:
		0 for snp
		int if its insertion or deletion
		";
let v = vec!["A", "T", "G", "C"];
if self.ref&&v.iter().any(|&x| x == self.alt) {
return 0;
} else {
if self.alt == "-" {
return self.ref.len();
} else {
if self.ref == "-" {
return self.alt.len();
} else {
if self.ref.len() == self.alt.len() {
return self.alt.len();
}
}
}
}
}
fn good_barcodes<T0, T1, RT>(cls: T0, f: T1) -> RT {
"

		:param f: good barcodes file
		:return: dict with good barcodes
		";
let barc = HashMap::new();
// with!(open(f, "r") as barcodes) //unsupported
{
for line in barcodes {
let lines = line.strip();
barc[lines] = 1;
}
}
return barc;
}
fn count_barcodes<T0, T1, T2, T3, T4, RT>(&self, bam_fil: T0, bar: T1, indel: T2, mapq: T3, baseq: T4) -> RT {
"

		:param bam_fil: bam_file
		:param bar: good barcodes dict
		:param indel: SNP = 0 indel= int(window size) insertion or deletion
		:param mapq: mapping quality
		:param baseq: base quality
		:return: two dicts 1) list of all barcodes at given var pos
		2) dict with ref and alt barcodes
		";
let barcodes = defaultdict(list);
let barUcodes = defaultdict(list);
let bar_count = [("ref", vec![]), ("alt", vec![])].iter().cloned().collect::<HashMap<_,_>>();
// with!(pysam.AlignmentFile(bam_fil, "rb") as pile) //unsupported
{
logger.info("processing variant: {}	{}	{}	{}	{}	{}	{}".format(self.chrm, self.start, self.end, self.ref, self.alt, self.gene, self.event));
for pileupcolumn in pile.pileup(self.chrm, (self.start - 1), self.start, true, "nofilter", 100000000) {
for read in pileupcolumn.pileups {
if read.alignment.has_tag("CB")&&read.alignment.has_tag("UB") {
if bar.iter().all(|&x| x != read.alignment.get_tag("CB")) {
continue;
}
if indel > 0 {
if read.is_refskip||i32::from(read.alignment.mapping_quality) < i32::from(mapq) {
continue;
}
} else {
if read.is_del||read.is_refskip||i32::from(read.alignment.query_qualities[read.query_position]) < i32::from(baseq)||i32::from(read.alignment.mapping_quality) < i32::from(mapq) {
continue;
}
}
let mut q_name = read.alignment.query_name;
let mut q_tag = read.alignment.get_tag("CB");
let mut qU_tag = read.alignment.get_tag("UB");
let mut qstring = ((q_tag + ":") + qU_tag);
barcodes["{}".format(q_name)].append(q_tag);
barUcodes["{}".format(q_name)].append(qstring);
if indel > 0 {
q_name = read.alignment.query_name;
q_tag = read.alignment.get_tag("CB");
qU_tag = read.alignment.get_tag("UB");
qstring = ((q_tag + ":") + qU_tag);
barcodes["{}".format(q_name)].append(q_tag);
barUcodes["{}".format(q_name)].append(qstring);
if self.alt == "-" {
if read.alignment.has_tag("CB")&&read.query_position == None {
let mut alt_tag = read.alignment.get_tag("CB");
let mut altU_tag = read.alignment.get_tag("UB");
let mut ustringa = ((alt_tag + ":") + altU_tag);
bar_count["alt"].append(ustringa);
} else {
let mut ref_tag = read.alignment.get_tag("CB");
let mut refU_tag = read.alignment.get_tag("UB");
let mut ustring = ((ref_tag + ":") + refU_tag);
bar_count["ref"].append(ustring);
}
} else {
if read.alignment.has_tag("CB")&&read.indel == indel&&read.alignment.query_alignment_sequence[(read.query_position + 1)..((read.query_position + read.indel) + 1)] == self.alt {
let mut alt_tag = read.alignment.get_tag("CB");
let mut altU_tag = read.alignment.get_tag("UB");
let mut ustringa = ((alt_tag + ":") + altU_tag);
bar_count["alt"].append(ustringa);
} else {
let mut ref_tag = read.alignment.get_tag("CB");
let mut refU_tag = read.alignment.get_tag("UB");
let mut ustring = ((ref_tag + ":") + refU_tag);
bar_count["ref"].append(ustring);
}
}
} else {
if read.alignment.query_sequence[read.query_position] != self.alt {
let mut ref_tag = read.alignment.get_tag("CB");
let mut refU_tag = read.alignment.get_tag("UB");
let mut ustring = ((ref_tag + ":") + refU_tag);
bar_count["ref"].append(ustring);
} else {
if read.alignment.query_sequence[read.query_position] == self.alt {
let mut alt_tag = read.alignment.get_tag("CB");
let mut altU_tag = read.alignment.get_tag("UB");
let mut ustringa = ((alt_tag + ":") + altU_tag);
bar_count["alt"].append(ustringa);
}
}
}
}
}
}
}
return (barUcodes, bar_count);
}
fn consensus_calling<T0, T1, T2, T3, T4, RT>(&self, total_barcodes: T0, ub_counts: T1, wt: T2, wu: T3, wc: T4) -> RT {
"

		:param total_barcodes: list of all the barcodes at variant pos
		:param ub_counts: ref/alt barcode counts
		:param wt: file handle to write vaf information for each variant
		:param wu: file handle to write CB:UB information for each variant
		:param wc: file handle to write CB information for each variant
		:return: variant information

		";
let uniqub_barcodes_raw = set();
let mut ub_raw = vec![];
for tupru in total_barcodes.values() {
for sampleru in tupru {
ub_raw.push(sampleru);
uniqub_barcodes_raw.add(sampleru);
}
}
logger.info("Consensus calc variant: {}	{}	{}	{}	{}	{}	{}".format(self.chrm, self.start, self.end, self.ref, self.alt, self.gene, self.event));
let d = HashMap::new();
let mut UBuniq_barcodes_raw = vec![];
for utags in uniqub_barcodes_raw {
if ub_counts["alt"].iter().any(|&x| x == utags)&&ub_counts["ref"].iter().any(|&x| x == utags) {
unref = ub_counts["ref"].count(utags);
unalt = ub_counts["alt"].count(utags);
utotl = (unref + unalt);
un1ref = ub_counts["ref"].count(utags);
un1alt = ub_counts["alt"].count(utags);
let ut1totl = (un1ref + un1alt);
if (unref/utotl) < float(0.75) {
if (unalt/utotl) < float(0.75) {
continue;
} else {
un1alt = 1;
un1ref = 0;
let mut u1totl = (un1ref + un1alt);
}
} else {
if (unalt/utotl) < float(0.75) {
un1alt = 0;
un1ref = 1;
u1totl = (un1ref + un1alt);
} else {
un1alt = 1;
un1ref = 1;
u1totl = (un1ref + un1alt);
}
}
} else {
un1ref = set(ub_counts["ref"]).collect::<Vec<_>>().count(utags);
un1alt = set(ub_counts["alt"]).collect::<Vec<_>>().count(utags);
let mut u1totl = (un1ref + un1alt);
unref = ub_counts["ref"].count(utags);
unalt = ub_counts["alt"].count(utags);
utotl = (unref + unalt);
}
if un1alt {
UBuniq_barcodes_raw.push(utags);
}
let ubtager = "{chrm}	{st}	{en}	{ref}	{alt}	{type}	{gene}	{bar}	{r}	{v}	{tot}
".format(self.chrm, self.start, self.end, self.ref, self.alt, self.event, self.gene, utags, unref, unalt, utotl);
wu.write(ubtager);
if d.iter().all(|&x| x != utags.split(":")[0]) {
d[utags.split(":")[0]] = [("ref", un1ref), ("alt", un1alt)].iter().cloned().collect::<HashMap<_,_>>();
} else {
d[utags.split(":")[0]]["ref"] += un1ref;
d[utags.split(":")[0]]["alt"] += un1alt;
}
}
let mut uni_alt = vec![];
let mut tuni_alt = vec![];
let mut CBuniq_barcodes_raw = vec![];
for (i, v) in d.items() {
if v["alt"] >= 1 {
CBuniq_barcodes_raw.push(i);
}
uni_alt.push(v["alt"]);
let tot1 = (v["alt"] + v["ref"]);
tuni_alt.push(tot1);
let cbtager = "{chrm}	{st}	{en}	{ref}	{alt}	{type}	{gene}	{bar}	{ref1}	{alt1}	{tot}
".format(self.chrm, self.start, self.end, self.ref, self.alt, self.event, self.gene, i, v["alt"], v["ref"], tot1);
wc.write(cbtager);
}
let Uuni_alt_c = uni_alt.iter().sum();
let tUuni_alt_c = tuni_alt.iter().sum();
let UBfin_umi = ",".join(UBuniq_barcodes_raw);
let CBfin_umi = ",".join(CBuniq_barcodes_raw);
let try_dummy = { //unsupported
let mut uvU = round((float(Uuni_alt_c)/float(tUuni_alt_c)), 2);
};
let except!(ZeroDivisionError) = { //unsupported
uvU = float(0.0);
};
let snp = "{chrm}	{st}	{en}	{ref}	{alt}	{type}	{gene}	{Uunidepth}	{Uunc}	{Uuvaf}	{umi}	{Uumi}
".format(self.chrm, self.start, self.end, self.ref, self.alt, self.event, self.gene, tUuni_alt_c, UBfin_umi, uvU, CBfin_umi, Uuni_alt_c);
wt.write(snp);
logger.info("completed processing variant: {}	{}	{}	{}	{}	{}	{}".format(self.chrm, self.start, self.end, self.ref, self.alt, self.gene, self.event));
return snp;
} 
}
fn main() {
let parser = argparse.ArgumentParser("Parse CB barcodes from Single cell rna seq data");
parser.add_argument("bam_file", "BAM file");
parser.add_argument("variant_file", "variants file with header");
parser.add_argument("barcodes", "list of good barcodes file");
parser.add_argument("upn", "upn/sample name: will be used as prefix for out_file");
parser.add_argument("-f", "--filter", int, 0, "number of reads required per barcode default: 0");
parser.add_argument("-mq", "--mapq", int, 0, "Skip read with mapq smaller than default : 0");
parser.add_argument("-bq", "--baseq", int, 1, "Skip bases with base quality less than default : 1");
let args = parser.parse_args();
let bam_file = args.bam_file;
let variants = args.variant_file;
let barcodes_good = args.barcodes;
let outfile = args.upn;
let counts = args.filter;
let bq = args.baseq;
let mq = args.mapq;
logger.basicConfig((outfile + ".log"), "w+", logger.DEBUG, "%(asctime)s %(levelname)s %(message)s");
logger.info("Start process");
let num_lines_variants = (open(variants).iter().map(|line| 1).collect::<Vec<_>>().iter().sum() - 1);
logger.info("Number of variants:	{}".format(num_lines_variants));
let num_lines_barcode = open(barcodes_good).iter().map(|lines| 1).collect::<Vec<_>>().iter().sum();
logger.info("Number of good barcodes:	{}".format(num_lines_barcode));
let bars = GenomicPosition::good_barcodes(barcodes_good);
// with!(open((outfile + "_AllCounts.tsv"), "w+") as var, open((outfile + "_counts_CB.tsv"), "w+") as CB, open((outfile + "_counts_UB.tsv"), "w+") as UB, open(variants, "r") as regions) //unsupported
{
let var_header = vec!["chr", "start", "end", "ref", "alt", "gene", "type", "UB_DEPTH", "UB_ALT", "UB_VAF", "CB_barcodes", "CB:UB_barcodes"];
let var_h = ("	".join(var_header) + "
");
var.write(var_h);
let CB_header = vec!["chrm", "start", "end", "ref", "alt", "type", "gene", "barcode", "ref_count", "alt_count", "total_CB"];
let CB_h = ("	".join(CB_header) + "
");
CB.write(CB_h);
let UB_header = vec!["chrm", "start", "end", "ref", "alt", "type", "gene", "barcode", "ref_count", "alt_count", "total_UB"];
let UB_h = ("	".join(UB_header) + "
");
UB.write(UB_h);
next(regions);
for lines in regions {
let x = GenomicPosition(lines);
println!("{:?} {:?} {:?} {:?} {:?} ",x.chrm, x.start, x.event, x.gene, x.classify());
let barcode_counts = x.count_barcodes(bam_file, bars, x.classify(), mq, bq);
let counters = x.consensus_calling(barcode_counts[0], barcode_counts[1], var, UB, CB);
}
}
}
logger.info("end process");