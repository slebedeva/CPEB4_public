### Analyze omniCLIP peaks for CPEB4

## This script will take raw omniCLIP output,
## merge overlapping peaks, assign them to genes and transcript regions, 
## and add the count of diagnostic events (T-to-C converisons and T deletions).
## It will output viewable bed files and the list of CPEB4 target genes.

## working dir is where this source file is
basedir <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(basedir)


###################################################### load environment ########################################################################################

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

############## load data ##############



## get genome sequence (only once)
if(!(file.exists("GRCh37.p13.genome.fa")|file.exists("GRCh37.p13.genome.fa.gz")|file.exists("GRCh37.p13.genome.fa.bgz"))){
  system("wget http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/GRCh37.p13.genome.fa.gz")
  system("gunzip -c GRCh37.p13.genome.fa.gz | bgzip  > GRCh37.p13.genome.fa.bgz")
}
hg19=("GRCh37.p13.genome.fa.bgz")


####### generate gencode v19 annotation data #################

## get gtf (only once)
if(!file.exists("gencode.v19.annotation.gtf.gz")){system("wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz")}
gtf <- rtracklayer::import(con = "gencode.v19.annotation.gtf.gz")
canonical_chr <- paste0("chr",c(1:22,"X","Y","M"))
## gene2name
gene2name <- unique(data.frame(gene_id = remove_dots(gtf$gene_id), gene_name=gtf$gene_name, gene_status=gtf$gene_status, gene_type=gtf$gene_type, stringsAsFactors = F))
tx2name <- with(subset(gtf, type=="transcript"), unique(data.frame(transcript_id, gene_name, stringsAsFactors = F)))

save(gtf,canonical_chr,gene2name,tx2name,file="data/gtf.RData")

########## generate transcript regions ###################

mytxdb <- makeTxDbFromGFF(file="gencode.v19.annotation.gtf.gz")
saveDb(mytxdb,file = "data/mytxdb.RData")

threeutrs <-threeUTRsByTranscript(mytxdb, use.names=T)
fiveutrs <- fiveUTRsByTranscript(mytxdb, use.names=T)
cdss <- cdsBy(mytxdb, by="tx", use.names=T)
intr <- intronsByTranscript(mytxdb, use.names=T)
ex <- exonsBy(mytxdb, by="tx", use.names=T)

threutrs <-unlist(threeutrs)
threutrs$region<-"3'UTR"
fivutrs <- unlist(fiveutrs)
fivutrs$region <- "5'UTR"
cdsss <- unlist(cdss)
cdsss$region <- "CDS"
## protein coding regions
pcg <- c(fivutrs,threutrs,cdsss)
## non coding 
ncexs <- unlist(ex)[!names(unlist(ex))%in%names(pcg)] ## noncoding=which do not have utrs or cdss
ncexs$region <- "non-coding"
int <- unlist(intr)
int$region <- "intron"
myregions <- c(pcg,ncexs,int)

save(myregions,file="data/myregions.RData")

##############import CLIP peaks (bed and txt) ##############################

## import omniCLIP predictions 
ctrl<-read.table("data/omniCLIP/Ctrl/pred.txt", header=T, stringsAsFactors = F)
rmd<-read.table("data/omniCLIP/RMD/pred.txt", header=T, stringsAsFactors = F)

cbed<-rtracklayer::import("data/omniCLIP/Ctrl/pred.bed")
rbed<-rtracklayer::import("data/omniCLIP/RMD/pred.bed")


############# make reduced binding sites  #################

## This will flatten all overlapping sites keeping information on their original names.
## Also it will fuse together any sites closer than 5nt apart because omniCLIP often splits sites with many deletions. It should reduce amount of too short (1nt and so on) sites a bit.
## I am mainly doing it to prevent having duplicate peaks in further analysis.
## This will be my main peak annotation file to work with.

make_reduced_Bed <- function(bed){
  bed_red <- GenomicRanges::reduce(bed, drop.empty.ranges = T, with.revmap=T, min.gapwidth=5)  ## fuse regions which are too close ## right now primitively with increasing gapwidth 
  bed_red$orig_names <- unlist(lapply(bed_red$revmap, function(x){ return(paste(bed$name[x],collapse=";")) }))
  bed_red$orig_scores <- unlist(lapply(bed_red$revmap, function(x){ return(paste(bed$score[x],collapse=";")) }))
  ## add biological gene names
  bed_red$gene_names <- unlist(lapply(bed_red$revmap, function(x){
    return(paste(
      sort(unique(gene2name$gene_name[match(sub("_.+","",bed$name[x]),gene2name$gene_id)]))
      ,collapse=";"))}))
  ## unique id
  bed_red$name <- make.names(bed_red$gene_names, unique = T)
  return(bed_red)
}

ctrl_red <- make_reduced_Bed(cbed)
rmd_red <- make_reduced_Bed(rbed)

ctrl_red$condition <- "ctrl"
rmd_red$condition <- "rmd"




########################## add TC conversions #################################

## do Rsamtools pileup on merged bam files to count TC conversions and Tdel 
## only canonical chr (in bed) 

## bam file list
allbams <- list.files("data/bams", pattern = ".bam$", full.names = T)

## this function will add number of conversions as metadata to the reduced peak bed file 
### as well as maxTC as thick position (if there are TCs)
add_conv_to_bed <- function(allbams, bed,myname){
  ## make label
  bed$which_label <- paste0(seqnames(bed),":",ranges(bed))
  p_param <- PileupParam(max_depth = 1e6) 
  sbp <- ScanBamParam(which=bed)
  mypile1<-pileup(file = BamFile(allbams[1]), scanBamParam = sbp, pileupParam = p_param)
  mypile2<-pileup(file = BamFile(allbams[2]), scanBamParam = sbp, pileupParam = p_param)
  mypile3<-pileup(file = BamFile(allbams[3]), scanBamParam = sbp, pileupParam = p_param)
  mypile <- dplyr::full_join(mypile1,mypile2,by=c("seqnames","pos","strand","nucleotide","which_label"),suffix=c("_1","_2"))
  mypile <- dplyr::full_join(mypile,mypile3,by=c("seqnames","pos","strand","nucleotide","which_label")) 
  ## how to handle NA? put to 0 counts 
  mypile[is.na(mypile)] <- 0
  gr <- with(mypile, GRanges(seqnames=seqnames,ranges = IRanges(start=pos,width=1),strand = "+")) ## if I want A to G I have to have always Watson strand
  mypile$ref <- as.character(getSeq(x=FaFile(hg19), gr))
  mycoverage <- 
    mypile%>%
    group_by(which_label)%>%
    summarize(
      coverage=sum(count,count_1,count_2, na.rm = T))
  TCs <- subset(mypile, (ref=="T" & nucleotide=="C" & strand=="+") | (ref=="A" & nucleotide=="G" & strand=="-") )
  TDs <- subset(mypile, (ref=="T" & nucleotide=="-" & strand=="+") | (ref=="A" & nucleotide=="-" & strand=="-") )
  # first summarize replicates per position and then max over all replicates per cluster
  TCperCluster <- TCs%>%group_by(which_label,pos)%>%
    summarize(max=max(count,count_1,count_2,na.rm = T),sum=sum(count,count_1,count_2,na.rm = T))%>%
    group_by(which_label)%>%
    summarize(thickTCSum=pos[which.max(sum)], thickTCMax=pos[which.max(max)],
              sumTC=sum(sum,na.rm = T),maxTC=max(max,na.rm = T))
  TDperCluster <- TDs%>%group_by(which_label,pos)%>%
    summarize(max=max(count,count_1,count_2,na.rm = T),sum=sum(count,count_1,count_2,na.rm = T))%>%
    group_by(which_label)%>%
    summarize(thickTDSum=pos[which.max(sum)], thickTDMax=pos[which.max(max)],
              sumTD=sum(sum,na.rm = T),maxTD=max(max,na.rm = T))
  
  ## add to bed and replace NA
  mcols(bed) <- cbind(mcols(bed), 
                      mycoverage[match(bed$which_label,mycoverage$which_label),"coverage"],
                      TDperCluster[match(bed$which_label,TDperCluster$which_label),-1],
                      TCperCluster[match(bed$which_label,TCperCluster$which_label),-1])
  bed$sumDE <- cbind(bed$sumTD,bed$sumTC)%>%rowSums(na.rm=T) ## convert NA to 0
  return(bed)}

ctrl1 <- add_conv_to_bed(allbams%>%.[!grepl("RMD",.)],ctrl_red)
rmd1 <- add_conv_to_bed(allbams%>%.[grepl("RMD",.)],rmd_red)

## this takes a long time, so save intermittent result
save(ctrl, rmd, ctrl1, rmd1, file="data/bed_reduced_tmp.RData")

######################### intersect reduced peaks with transcript regions ###############################

## function to add the region of overlap to peak bed file
assign_regions <- function(bed,regions){
  hits <- findOverlaps(bed,regions) 
  gr.over <- pintersect(bed[queryHits(hits)],regions[subjectHits(hits)]) ## range by range overlaps
  w <- width(gr.over) ## overlap width
  ## this will give overlap for every annotated transcript so most peaks will overlap many,many transcripts in the same region
  peak2regiontx <- data.table(stringsAsFactors = F,
                              OvWidth=w,
                              peakname=gr.over$name,
                              region=regions[subjectHits(hits)]$region,
                              tx=names(regions[subjectHits(hits)]))
  # summarize unique region per peak as well as all regions with respective transcripts and overlap width
  peak2region <- peak2regiontx[,list(region = paste(unique(region), collapse="_")
                                     ,region_tx_ov=paste(region,tx,OvWidth,collapse=";"))
                               ,by=peakname]
  # add metadata to flattened bed file with converisons
  mcols(bed) <- cbind(mcols(bed),peak2region[match(bed$name,peak2region$peakname),c("region","region_tx_ov")])
  return(bed)
}

ctrl2 <- assign_regions(ctrl1,myregions)
rmd2 <- assign_regions(rmd1,myregions)

## apparently one cluster has NA annotation because it is overlapping a gene annotation but no transcript annotation (between 2 transcripts). 
## I would remove it since it does not make any sense.
ctrl2 <- ctrl2[!is.na(ctrl2$region)]
rmd2 <- rmd2[!is.na(rmd2$region)]

########################## group reduced peaks by DE strength ############################

## define cut offs for target groups by diagnostic events
mycuts=c(0,1,10,100,1000,5000)

ctrl2$target_group <- Hmisc::cut2(ctrl2$sumDE, cuts = mycuts)%>%gsub("\\[| ","",.)%>%gsub("10\\)","9",.)%>%gsub("100\\)","99",.)%>%gsub("1000\\)","999",.)%>%gsub("5000\\)","4999",.)%>%gsub("38127\\]","more",.)%>%gsub("\\,","\\-",.)
rmd2$target_group <- Hmisc::cut2(rmd2$sumDE,  cuts = mycuts)%>%gsub("\\[| ","",.)%>%gsub("10\\)","9",.)%>%gsub("100\\)","99",.)%>%gsub("1000\\)","999",.)%>%gsub("5000\\)","4999",.)%>%gsub("9399\\]","more",.)%>%gsub("\\,","\\-",.)

################ make viewable bed to load into igv ###########################

ctrl2$score <- ctrl2$sumDE
ctrl2$name2 <- ctrl2$name
ctrl2$name <- paste(ctrl2$name,ctrl2$condition,sep="@")
## make "thick" point a max of TC or TD and if DEs are absent, just the middle of the cluster
ctrl2$thick <- IRanges(start=ifelse(!is.na(ctrl2$thickTCMax),ctrl2$thickTCMax,
                                    ifelse(!is.na(ctrl2$thickTDMax),ctrl2$thickTDMax,
                                           (start(ranges(ctrl2))+round(width(ranges(ctrl2))/2)))),
                                    width=1)
export.bed(ctrl2, "data/beds/CPEB4_Ctrl_hg19.bed")
rmd2$score <- rmd2$sumDE
rmd2$name2 <- rmd2$name
rmd2$name <- paste(rmd2$name,rmd2$condition,sep="@")
rmd2$thick <- IRanges(start=ifelse(!is.na(rmd2$thickTCMax),rmd2$thickTCMax,
                                    ifelse(!is.na(rmd2$thickTDMax),rmd2$thickTDMax,
                                           (start(ranges(rmd2))+round(width(ranges(rmd2))/2)))),
                       width=1)

export.bed(rmd2, "data/beds/CPEB4_RMD_hg19.bed") 

############ this is our main bed peak file for further analysis: ####################

save(rmd, ctrl,rmd2,ctrl2,file="data/bed_reduced.RData") 

system("rm data/bed_reduced_tmp.RData")

###################### summary table of target genes ####################################


#load("data/bed_reduced.RData")

## About ambigous target genes (overlapping genes on the same strand that share same peak):
## I will keep them in "gene groups" (like protein groups in mass spec)

target_genes <- data.frame(c(rmd2,ctrl2))%>%dplyr::group_by(gene_names)%>%dplyr::summarize(conditionList=str_c(unique(condition), collapse = ";")
                                                                                  , regionList=str_c(unique(region), collapse = ";")
                                                                                  , whichLabel=str_c(which_label, collapse=";")
                                                                                  , orig_names=str_c(orig_names,collapse = ";")
                                                                                  , orig_scores=str_c(orig_scores,collapse = ";")
                                                                                  , sumCovPerGene=sum(coverage, na.rm = T)
                                                                                  , sumTCPerGene=sum(sumTC,na.rm = T)
                                                                                  , sumTDPerGene=sum(sumTD,na.rm = T)
                                                                                  , sumDEPerGene=sum(sumDE,na.rm = T))%>%dplyr::arrange(desc(sumDEPerGene)) 

## add target groups per gene and unique reigons 
target_genes$gene_target_group <- Hmisc::cut2(target_genes$sumDEPerGene,  cuts = mycuts)%>%gsub("\\[| ","",.)%>%gsub("10\\)","9",.)%>%gsub("100\\)","99",.)%>%gsub("1000\\)","999",.)%>%gsub("5000\\)","4999",.)%>%gsub("46027\\]","more",.)%>%gsub("\\,","\\-",.)

## final table of target genes
## this is also the Supplementary table S6
#write.csv(target_genes, file = "table_S6_cpeb4_target_genes.csv", row.names = F) 

save(target_genes,file="data/target_genes.RData")

############### version info for reproducibility #################

sink(file="R_version_info.txt")
sessionInfo()
sink()

