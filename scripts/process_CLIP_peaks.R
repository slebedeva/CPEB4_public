### Analyze omniCLIP peaks for CPEB4

## This script will generate reduced CLIP peak objects with counts of diagnostic events and annotation of transcript regions and genes

## Only use it within manuscript_figures.R otherwise all paths will be broken!
if(!exists("basedir")){message("Please source me within manuscript_figures.R")}
message("Processing omniCLIP peaks... This will take some time.")



###################################################### load environment ########################################################################################

## important! otherwise Hs.org.db fails ##https://support.bioconductor.org/p/9136239/#9136319
options(connectionObserver = NULL)

mypackages <- c("config" ## careful it masks "merge"!!
                ,"plyranges"
                ,"tidyverse"
                ,"reshape2"
                ,"Hmisc"
                ,"magrittr"
                ,"ggplot2"
                ,"ggrepel"
                ,"ggpubr"
                ,"ggunchained"
                ,"lattice"
                ,"DESeq2"
                ,"Biostrings"
                ,"GenomicRanges"
                ,"GenomicFeatures"
                ,"rtracklayer"
                ,"Rsamtools"
                ,"motifStack"
                ,"Gviz"
                ,"ellipse"
                ,"BSgenome.Hsapiens.UCSC.hg19"
                ,"genomation"
                ,"RCAS"
                ,"SRAdb"
                ,"fgsea"
                ,"cowplot"
                ,"gridGraphics"
                ,"magick"
                ,"rio"
                ,"extrafont"
		,"data.table"
)

pl=suppressWarnings(suppressPackageStartupMessages(lapply(mypackages, require, character.only=T)))
if(sum(!unlist(pl))){message("failed to load: \n",paste(mypackages[!unlist(pl)],collapse = "\n"), "\nanyway proceeding...")}





################## load annotation ###################

config=config::get()

message("Loading annotation...")

hg19=config$genome
mytxdb=loadDb(file = file.path(ann_dir,"mytxdb.RData"))
for(x in c("gtf.RData","myregions.RData")){
  message("loading ", file.path(ann_dir,x))
  load(file.path(ann_dir,x))
}

############## import raw omniCLIP peaks (bed and txt) ##############################

## import omniCLIP predictions 
ctrl<-read.table(file.path(datadir,"omniCLIP/Ctrl/pred.txt"), header=T, stringsAsFactors = F)
rmd<-read.table(file.path(datadir,"omniCLIP/RMD/pred.txt"), header=T, stringsAsFactors = F)

cbed<-rtracklayer::import(file.path(datadir,"omniCLIP/Ctrl/pred.bed"))
rbed<-rtracklayer::import(file.path(datadir,"omniCLIP/RMD/pred.bed"))


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


message("Making reduced peaks....")
ctrl_red <- make_reduced_Bed(cbed)
rmd_red <- make_reduced_Bed(rbed)

ctrl_red$condition <- "ctrl"
rmd_red$condition <- "rmd"




########################## add TC conversions #################################

## do Rsamtools pileup on merged bam files to count TC conversions and Tdel 
## only canonical chr (in bed) 

## bam file list
allbams <- list.files(file.path(datadir,"bams"), pattern = ".bam$", full.names = T)

## this function will add number of conversions as metadata to the reduced peak bed file 
### as well as maxTC as thick position (if there are TCs) - otherwise it is just a center of the cluster
### this function uses (R)samtools pileup
add_conv_to_bed <- function(allbams, bed,myname){
  ## make label (important: keep strand information!)
  bed$which_label <- paste0(seqnames(bed),":",ranges(bed),":",strand(bed))
  p_param <- PileupParam(max_depth = 1e6) 
  sbp <- ScanBamParam(which=bed)
  mypile1<-pileup(file = BamFile(allbams[1]), scanBamParam = sbp, pileupParam = p_param)
  gc()
  mypile2<-pileup(file = BamFile(allbams[2]), scanBamParam = sbp, pileupParam = p_param)
  gc()
  mypile3<-pileup(file = BamFile(allbams[3]), scanBamParam = sbp, pileupParam = p_param)
  gc()
  mypile <- dplyr::full_join(mypile1,mypile2,by=c("seqnames","pos","strand","nucleotide","which_label"),suffix=c("_1","_2"))
  mypile <- dplyr::full_join(mypile,mypile3,by=c("seqnames","pos","strand","nucleotide","which_label")) 
  ## how to handle NA? put it to 0 counts 
  mypile[is.na(mypile)] <- 0
  ## get the reference nucleotide
  gr <- with(mypile, GRanges(seqnames=seqnames,ranges = IRanges(start=pos,width=1),strand = "+")) ## if I want A to G I have to have always Watson strand - bam alignment always is in respect to Watson strand
  mypile$ref <- as.character(getSeq(x=FaFile(hg19), gr))
  mycoverage <- 
    mypile%>%
    dplyr::group_by(which_label,strand)%>%
    dplyr::summarize(
      coverage=sum(count,count_1,count_2, na.rm = T))
  mycoverage$which_label=paste(mycoverage$which_label,mycoverage$strand,sep=":")
  mycoverage$strand=NULL # strand is kept in the label, otherwise granges with this metadata column cannot be created
  TCs <- subset(mypile, (ref=="T" & nucleotide=="C" & strand=="+") | (ref=="A" & nucleotide=="G" & strand=="-") )
  TDs <- subset(mypile, (ref=="T" & nucleotide=="-" & strand=="+") | (ref=="A" & nucleotide=="-" & strand=="-") )
  # quickly check that my collection makes sense: there should be much less of TAs and TGs - correct!
  #TAs <- subset(mypile, (ref=="T" & nucleotide=="A" & strand=="+") | (ref=="A" & nucleotide=="T" & strand=="-") )
  #TGs <- subset(mypile, (ref=="T" & nucleotide=="G" & strand=="+") | (ref=="A" & nucleotide=="C" & strand=="-") )
  DEs=subset(mypile, (ref=="T" & nucleotide%in%c("C","-") & strand=="+") | (ref=="A" & nucleotide%in%c("G","-") & strand=="-") )#insert also a summary of all DEs together
  # first summarize replicates per position and then max over all replicates per cluster
  TCperCluster <- TCs%>%dplyr::group_by(which_label,pos,strand)%>%
    dplyr::summarize(max=max(count,count_1,count_2,na.rm = T),sum=sum(count,count_1,count_2,na.rm = T))%>%
    dplyr::group_by(which_label,strand)%>%
    dplyr::summarize(thickTCSum=pos[which.max(sum)], thickTCMax=pos[which.max(max)],
              sumTC=sum(sum,na.rm = T),maxTC=max(max,na.rm = T))
  TDperCluster <- TDs%>%dplyr::group_by(which_label,pos,strand)%>%
  dplyr::summarize(max=max(count,count_1,count_2,na.rm = T),sum=sum(count,count_1,count_2,na.rm = T))%>%
  dplyr::group_by(which_label,strand)%>%
  dplyr::summarize(thickTDSum=pos[which.max(sum)], thickTDMax=pos[which.max(max)],
            sumTD=sum(sum,na.rm = T),maxTD=max(max,na.rm = T))
  ## count all DEs together
  DEperCluster <- DEs%>%dplyr::group_by(which_label,pos,strand)%>%
  dplyr::summarize(max=max(count,count_1,count_2,na.rm = T),sum=sum(count,count_1,count_2,na.rm = T))%>%
  dplyr::group_by(which_label,strand)%>%
  dplyr::summarize(thickDESum=pos[which.max(sum)], thickDEMax=pos[which.max(max)],
            sumDE=sum(sum,na.rm = T),maxDE=max(max,na.rm = T))
  DEperCluster$which_label=paste(DEperCluster$which_label,DEperCluster$strand,sep=":")
  TCperCluster$which_label=paste(TCperCluster$which_label,TCperCluster$strand,sep=":")
  TDperCluster$which_label=paste(TDperCluster$which_label,TDperCluster$strand,sep=":")
  DEperCluster$strand=NULL
  TDperCluster$strand=NULL
  TCperCluster$strand=NULL
  ## add to bed and replace NA
  mcols(bed) <- cbind(mcols(bed) 
                      ,mycoverage[match(bed$which_label,mycoverage$which_label),"coverage"]
                      ,DEperCluster[match(bed$which_label,DEperCluster$which_label),-1]
                      ,TCperCluster[match(bed$which_label,TCperCluster$which_label),-1]
                      ,TDperCluster[match(bed$which_label,TDperCluster$which_label),-1]
                      )
  bed$sumDE = (cbind(bed$sumTD,bed$sumTC)%>%rowSums(na.rm=T)) ## sum of deletions and conversions ## converts NA to 0 - important downstream!
  return(bed)}

message("Counting diagnostic events...")

ctrl1 <- add_conv_to_bed(allbams%>%.[!grepl("RMD",.)],ctrl_red)
rmd1 <- add_conv_to_bed(allbams%>%.[grepl("RMD",.)],rmd_red)

## this takes a long time, so save intermittent result
save(ctrl, rmd, ctrl1, rmd1, file=file.path(resultdir,"bed_reduced_tmp.RData"))

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

load(file.path(resultdir,"bed_reduced_tmp.RData"))

message("Assigning transcript regions...")

ctrl2 <- assign_regions(ctrl1,myregions)
rmd2 <- assign_regions(rmd1,myregions)

## apparently one cluster has NA annotation because it is overlapping a gene annotation but no transcript annotation (this gene has 2 disconnected transcript isoforms). 
## I would remove it since it does not make any sense.
ctrl2 <- ctrl2[!is.na(ctrl2$region)]
rmd2 <- rmd2[!is.na(rmd2$region)]

########################## group reduced peaks by DE strength ############################

## define cut offs for target groups by diagnostic events
mycuts=c(0,1,10,100,1000,5000)

ctrl2$target_group <- Hmisc::cut2(ctrl2$sumDE, cuts = mycuts)%>%gsub("\\[| ","",.)%>%gsub("10\\)","9",.)%>%gsub("100\\)","99",.)%>%gsub("1000\\)","999",.)%>%gsub("5000\\)","4999",.)%>%gsub("5000.+","5000,more",.)%>%gsub("\\,","\\-",.)
rmd2$target_group <- Hmisc::cut2(rmd2$sumDE,  cuts = mycuts)%>%gsub("\\[| ","",.)%>%gsub("10\\)","9",.)%>%gsub("100\\)","99",.)%>%gsub("1000\\)","999",.)%>%gsub("5000\\)","4999",.)%>%gsub("5000.+","5000,more",.)%>%gsub("\\,","\\-",.)

################ make viewable bed to load into igv ###########################

message("Preparing bed files...")
ctrl2$score <- ctrl2$sumDE
ctrl2$score[is.na(ctrl2$score)]<-0
ctrl2$name2 <- paste(ctrl2$name,ctrl2$condition,sep="@")

## for fair comparison of revision, keep thick as before (max)
#ctrl2$thick <- IRanges(start=ifelse(!is.na(ctrl2$thickTCMax),ctrl2$thickTCMax,
#                                    ifelse(!is.na(ctrl2$thickTDMax),ctrl2$thickTDMax,
#                                           (start(ranges(ctrl2))+round(width(ranges(ctrl2))/2)))),
#                       width=1)
## make "thick" point a max of TC (prioritize) or TD and if DEs are absent, just the middle of the cluster
## since the score is sum DE, it would also make sense to have sumDE as the thick point
## if there is no DE we take midpoint but we should put NA actually!
# ctrl2$thick <- IRanges(start=ifelse(!is.na(ctrl2$thickDESum),ctrl2$thickTESum
                        ##, (start(ranges(ctrl2))+round(width(ranges(ctrl2))/2))))
#                                     ,NA  
#                                     ,width=1)
rmd2$score <- rmd2$sumDE
rmd2$score[is.na(rmd2$score)]<-0
rmd2$name2 <- paste(rmd2$name,rmd2$condition,sep="@")
### old way
#rmd2$thick <- IRanges(start=ifelse(!is.na(rmd2$thickTCMax),rmd2$thickTCMax,
#                                   ifelse(!is.na(rmd2$thickTDMax),rmd2$thickTDMax,
#                                          (start(ranges(rmd2))+round(width(ranges(rmd2))/2)))),
#                      width=1)
## new way
# rmd2$thick <- IRanges(start=ifelse(!is.na(rmd2$thickDESum),rmd2$thickTESum
#                                     ,NA  
#                                     ,width=1)

## we cannot use NA in Granges, so we will just take a middle if the peak has no DEs.
ctrl2$thick <- IRanges(start=ifelse(!is.na(ctrl2$thickDESum),ctrl2$thickDESum,
                                    (start(ranges(ctrl2))+round(width(ranges(ctrl2))/2))),
                       width=1)

rmd2$thick <- IRanges(start=ifelse(!is.na(rmd2$thickDESum),rmd2$thickDESum,
                                   (start(ranges(rmd2))+round(width(ranges(rmd2))/2))),
                      width=1)


export.bed(ctrl2, file.path(resultdir, "CPEB4_Ctrl_hg19.bed"))
export.bed(rmd2, file.path(resultdir, "CPEB4_RMD_hg19.bed")) 

############ this is our main bed peak file for further analysis: ####################

save(rmd, ctrl,rmd2,ctrl2,file= file.path(resultdir, "bed_reduced.RData")) 

###################### summary table of target genes ####################################

## About ambigous target genes (overlapping genes on the same strand that share same peak):
## I will keep them in "gene groups" (like protein groups in mass spec)

message("Writing the table of target genes...")

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
target_genes$gene_target_group <- 
  Hmisc::cut2(target_genes$sumDEPerGene,  cuts = mycuts)%>%
  gsub("\\[| ","",.)%>%gsub("10\\)","9",.)%>%gsub("100\\)","99",.)%>%gsub("1000\\)","999",.)%>%gsub("5000\\)","4999",.)%>%
  gsub("[0-9]+\\]","more",.)%>%gsub("\\,","\\-",.)
target_genes

save(target_genes,file=file.path(resultdir, "target_genes.RData"))

write.csv(target_genes, file = file.path(resultdir,"table_S6_cpeb4_target_genes.csv"), row.names = F) 

