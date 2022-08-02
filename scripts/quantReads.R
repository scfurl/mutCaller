require("GenomicAlignments")
require("Rsamtools")
require("rtracklayer")
require("BiocParallel")
register(MulticoreParam(workers=36))
mode <- "Union"
gtffile <-file.path(Sys.getenv("transcriptome"), "/genes/genes.gtf")
yield_size <- 8e6
if(!file.exists(gtffile)){stop(paste0("GTF file not found: ", gtffile))}
gff0 <- import(gtffile)
idx <- mcols(gff0)$type == "exon"
genes<- split(gff0[idx], mcols(gff0[idx])[["gene_id"]])
filenames <- "Aligned.sortedByCoord.out.bam"
if(!all(file.exists(filenames))){stop("All Bam files not found")}
bamfiles <- BamFileList(filenames, yieldSize=yield_size)
names(bamfiles) <- filenames
se <- summarizeOverlaps(features=genes, 
                          reads=bamfiles,
                          mode=mode,
                          singleEnd=TRUE,
                          ignore.strand=TRUE)
mcols(se)[["gene_id"]]=mcols(gff0)[["gene_id"]][match(rownames(se), mcols(gff0)[["gene_id"]])]
mcols(se)[['gene_short_name']]<-mcols(gff0)[["gene_name"]][match(rownames(se), mcols(gff0)[["gene_id"]])]
colnames(se)<-filenames
colData(se)<-DataFrame(df)
saveRDS(se, "summarizedExp.RDS")

