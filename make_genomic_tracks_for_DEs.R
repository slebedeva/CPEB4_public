## tracks for DEs

## working dir is where this source file is
basedir <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(basedir)

## important! otherwise Hs.org.db fails ##https://support.bioconductor.org/p/9136239/#9136319
options(connectionObserver = NULL)

mypackages <- c("plyranges"
                ,"tidyverse"
                ,"reshape2"
                ,"magrittr"
                ,"Biostrings"
                ,"GenomicRanges"
                ,"GenomicFeatures"
                ,"rtracklayer"
                ,"Rsamtools"
                ,"rio"
)
suppressPackageStartupMessages(lapply(mypackages, require, character.only=T))

load("")
hg19=("GRCh37.p13.genome.fa.bgz")

## we want to make bedgraph for plus and minus strand for each of the bam file with TC and Tdel
## it has to be separated by strand of the cluster 

## seq info for genomic track export, only take canoncial chr
canonical_chr <- paste0("chr",c(1:22,"X","Y","M"))
myseqinfo=Seqinfo(genome = "hg19")[canonical_chr]


make_DE_bigwigs=function(bam, bed, myname){
  bed=rtracklayer::import(bed)
  ## make label
  bed$which_label <- paste0(seqnames(bed),":",ranges(bed))
  p_param <- PileupParam(max_depth = 1e6) 
  sbp <- ScanBamParam(which = bed)
  mypile<-pileup(file = BamFile(bam), scanBamParam = sbp, pileupParam = p_param)
  
  ## how to handle NA? put to 0 counts 
  mypile[is.na(mypile)] <- 0
  gr <- with(mypile, GRanges(seqnames=seqnames,ranges = IRanges(start=pos,width=1),strand = strand))
  mypile$ref <- as.character(getSeq(x=FaFile(hg19), gr))
  
  ## need strand of the cluster
  bedstr=data.table(which_label=bed$which_label,bedstrand=as.character(strand(bed)))
  mypile$bedstrand=bedstr$bedstrand[match(mypile$which_label,bedstr$which_label)]
  
  #'TC'
  TCs <- subset(mypile, (ref=="T" & nucleotide=="C" & strand=="+" & bedstrand=="+") |
                  (ref=="A" & nucleotide=="G" & strand=="-" & bedstrand=="-") )
  TCgr <- suppressWarnings(
    with(TCs, GRanges(seqnames=seqnames,ranges = IRanges(start=pos,width=1),strand=strand, seqinfo = myseqinfo)))
  TCgr$score=ifelse(TCs$strand=="+",TCs$count,-TCs$count)
  rtracklayer::export.bw(subset(TCgr,strand=="+"), paste0("data/bigwigs/",myname,"_TC.plus.bw"))
  rtracklayer::export.bw(subset(TCgr,strand=="-"), paste0("data/bigwigs/",myname,"_TC.minus.bw"))
  
  # T deletions
  TDs <- subset(mypile, (ref=="T" & nucleotide=="-" & strand=="+" & bedstrand=="+") |
                  (ref=="A" & nucleotide=="-" & strand=="-" & bedstrand=="-") )
  TDgr <- suppressWarnings(
    with(TDs, GRanges(seqnames=seqnames,ranges = IRanges(start=pos,width=1),strand=strand, seqinfo = myseqinfo)))
  TDgr$score=ifelse(TDs$strand=="+",TDs$count,-TDs$count)
  rtracklayer::export.bw(subset(TDgr,strand=="+"), paste0("data/bigwigs/",myname,"_Tdel.plus.bw"))
  rtracklayer::export.bw(subset(TDgr,strand=="-"), paste0("data/bigwigs/",myname,"_Tdel.minus.bw"))
}

## ctrl
CtrlBams=allbams%>%.[grep("DMSO|Ctrl",.)]
mynames=CtrlBams%>%sub("data/bams/","",.)%>%sub(".merged_uniq.bam","",.)
for(i in 1:length(CtrlBams)){
  make_DE_bigwigs(bed="data/beds/CPEB4_Ctrl_hg19.bed"
                  ,bam=CtrlBams[i],myname = mynames[i])
}

## rmd
RMDbams=allbams%>%.[grep("RMD",.)]
mynames=RMDbams%>%sub("data/bams/","",.)%>%sub(".merged_uniq.bam","",.)
for(i in 1:length(RMDbams)){
  make_DE_bigwigs(bed="data/beds/CPEB4_RMD_hg19.bed"
                  ,bam=RMDbams[i],myname = mynames[i])
}

## txt table for target genes (because paper only has excel and it is a no)
load("data/target_genes.RData")
fwrite(target_genes, file = "CPEB4_target_genes.csv")