#### make figures for CPEB4 manuscript

## working dir is where this source file is
basedir <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(basedir)
message("your working directory is: ", getwd())

## put plots separately in a folder
plotdir <- file.path(basedir,"plots")
if(! dir.exists(plotdir)){dir.create(plotdir)}
## folder for source data
datadir <- file.path(basedir,"data/source_data")
if(! dir.exists(datadir)){dir.create(datadir)}



###################################################### load environment ########################################################################################

## important! otherwise Hs.org.db fails ##https://support.bioconductor.org/p/9136239/#9136319
options(connectionObserver = NULL)

mypackages <- c("plyranges"
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
                )

#### try to install missing packages if any (can easily fail) #####
to_install=mypackages[!mypackages%in%installed.packages()]

if(length(to_install)>0){
  message("missing packages: ", paste(to_install,collapse=", "), "\ntrying to install...")
  BiocManager::install(pkgs = to_install,update = FALSE, ask = FALSE)
  failed=mypackages[!mypackages%in%installed.packages()]
  message("failed to install: ", paste(failed,collapse=", "), "\nanyway proceeding...")
}

pl=suppressWarnings(suppressPackageStartupMessages(lapply(mypackages, require, character.only=T)))
if(sum(!unlist(pl))){message("failed to load: \n",paste(mypackages[!unlist(pl)],collapse = "\n"), "\nanyway proceeding...")}


############### generate data ###################

if(!file.exists("data/bed_reduced.RData")){source("process_CLIP_peaks.R")}
if(!file.exists("data/DEseq.RData")){source("Differential_expression_analysis.R")}


############## load data ##############

## font size for pltos
global_size=14

## custom color
my_palette = c("grey10",RColorBrewer::brewer.pal(n = 9, "Blues")[2:7])

## define cut offs for target groups by diagnostic events
mycuts=c(0,1,10,100,1000,5000)

## get genome sequence (only once)
if(!(file.exists("GRCh37.p13.genome.fa")|file.exists("GRCh37.p13.genome.fa.gz")|file.exists("GRCh37.p13.genome.fa.bgz"))){
  system("wget http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/GRCh37.p13.genome.fa.gz")
  system("gunzip -c GRCh37.p13.genome.fa.gz | bgzip  > GRCh37.p13.genome.fa.bgz")
  }
hg19=("GRCh37.p13.genome.fa.bgz")

## annotation
## get gencode v19 annotation data 
load("data/gtf.RData")

## transcript regions
load("data/myregions.RData")

## make reduced utrs 
ncbi_utrs <- rtracklayer::import("annotation/ncbi_refseq_utrs.bed") ## downloaded from ucsc table browser
red_utrs <- GenomicRanges::reduce(ncbi_utrs,with.revmap=T) ##names needed for mapToTx
canonical_chr <- paste0("chr",c(1:22,"X","Y","M"))
red_utrs <- red_utrs%>%subset(., .@seqnames%in%canonical_chr)
red_utrs$name <- ncbi_utrs$name[unlist(lapply(red_utrs$revmap, function(x) return(x[[1]])))]
names(red_utrs) <- red_utrs$name


## merged CLIP peaks
load("data/bed_reduced.RData") ## these peaks are reduced (overlapping peak calls are merged) and have regions assigned and DE count


## import target genes
load("data/target_genes.RData")

## differential expresison
load("data/DEseq.RData") 


## bam file list
allbams <- list.files("data/bams", pattern = ".bam$", full.names = T)
## index bam (only once)
if(!file.exists("data/bams/RMD2_CLIP.merged_uniq.bam.bai")){lapply(allbams,function(x) Rsamtools::indexBam(x))}

### get ARE and IEG gene databases

## add AREs 
if(!file.exists("Complete 3 UTR Data.xls")){
  system("wget https://brp.kfshrc.edu.sa/ared/Content/Uploads/Complete%203%20UTR%20ARE%20Data.rar") 
  system("unrar e Complete\ 3\ UTR\ ARE\ Data.rar") ## need to install it on a container
}
AREdb=rio::import('Complete 3 UTR Data.xls')

## add fantom IEG 
fantom_ieg=rio::import("https://ndownloader.figshare.com/files/3263948", skip=1)%>%dplyr::slice(1:232)

##clean up 
if(file.exists("wget-log")){system("rm wget-log*")}


########### Figures ################

########### Fig 4 ##########

############ 4B: pie charts transcript categories for sites > 10 DEs #####################

## same color scheme everywhere - put all shared into "multiple" category otherwise too many colors!
ctrl2$tx_region <- ifelse(!grepl("_",ctrl2$region),ctrl2$region,
                          ifelse(grepl("3'UTR",ctrl2$region),"3UTR or other","multiple non-3'UTR"))
rmd2$tx_region <- ifelse(!grepl("_",rmd2$region),rmd2$region,
                         ifelse(grepl("3'UTR",rmd2$region),"3UTR or other","multiple non-3'UTR"))
allreg <- levels(factor(c(ctrl2$tx_region,rmd2$tx_region)))
myCol <- rev(RColorBrewer::brewer.pal(n = length(allreg), "BrBG"))
names(myCol) <- allreg


## only show pie chart for >=10 DEs
mypiec <- mcols(ctrl2)[,c("tx_region","target_group")]%>%as.data.frame%>%dplyr::filter(!target_group%in%(ctrl2$target_group %>% factor %>% levels %>% .[1:2]))%>%group_by(tx_region)%>%tally()%>%mutate(sum=sum(n), perc=n/sum, condition="DMSO")
mypier <- mcols(rmd2)[,c("tx_region","target_group")]%>%as.data.frame%>%dplyr::filter(!target_group%in%(rmd2$target_group %>% factor %>% levels %>% .[1:2]))%>%group_by(tx_region)%>%tally()%>%mutate(sum=sum(n), perc=n/sum, condition="RMD")
mypie <- rbind(mypiec,mypier)
ppie <-
  ggplot(mypie, aes(x="", y=perc, fill=tx_region))+
  geom_bar(stat="identity", color="grey") +
  coord_polar(theta = "y", direction = 1, start=0)+
  theme_void(base_size=global_size)+
  scale_fill_manual(values = myCol, name = "sites with >= 10 DEs")+
  facet_wrap(~condition,dir="v")+
  NULL

write_csv(mypie, file = "data/source_data/Fig_4B.csv")

save_plot("plots/Fig_4B.pdf",ppie)

################# 4D: scatterplot of motif counts #################################

### calculate background: what is kmer count in all UTRs?
## use reduced UTRs to avoid double counting same sequence

## OK my background is 3UTRs and I count kmers in everything... 
## but how would I scale intronic background to proportion of introns? not sure 
## ncbi UTR
myseq <- RNAStringSet(x = getSeq(red_utrs,x=FaFile(hg19)))
## gencode UTR
#myseq <- RNAStringSet(x = getSeq(reduce_ranges(subset(myregions$region=="3'UTR")),x=FaFile(hg19)))
myfreq <- oligonucleotideFrequency(myseq, width=6, step=1) ##6mers
mykmercount <- colSums(myfreq)%>%sort(decreasing=T)
forgg <- mykmercount%>%as.data.frame%>%rownames_to_column(var="kmer")

## can also count in random 50 nt windows in UTRs
## more correctly: in random windows taken by proportion of the regions same as CLIP (by count or length)
## get all reduced regions (exons and introns)
#all_regions_red=reduce_ranges(myregions)


## separate intronic and UTR regions
# kmer_list=list()
# for(region in c("3'UTR","intron")){
#   myseq <- RNAStringSet(x = getSeq(reduce_ranges(subset(myregions,region==region)),x=FaFile(hg19)))
#   myfreq <- oligonucleotideFrequency(myseq, width=6, step=1) ##6mers
#   mykmercount <- colSums(myfreq)%>%sort(decreasing=T)
#   kmer_list[[region]] <- mykmercount%>%as.data.frame%>%rownames_to_column(var="kmer")
# }
# 

###make granges for crosslink centers ("thick" position= highest DE per cluster per rep)
## not filtered out sites with 0 DEs - their center is just middle point ##
myxlc <- GRanges(seqnames = ctrl2@seqnames,
                 ranges = ctrl2$thick,
                 strand = ctrl2@strand,
                 region=ctrl2$region,
                 name=ctrl2$name,
                 score=ctrl2$score)
myxlr <- GRanges(seqnames = rmd2@seqnames,
                 ranges = rmd2$thick,
                 strand = rmd2@strand,
                 region=rmd2$region,
                 name=rmd2$name,
                 score=rmd2$score)


## main figure: kmer count
k=6 # count 6-mers having in mind AAUAAA
mywindow=50
## the window is asymmetrical because the last 6-mer will fit on position 45:50 
count_kmers <- function(myxl, k=6, mywindow=50){
  myxl_cond <- plyranges::mutate(width = mywindow, anchor_center(myxl)) ## make crosslink centered regions of 50nt (xl will be 25)
  myseq <- RNAStringSet(x = getSeq(myxl_cond,x=FaFile(hg19)))
  myfreq <- oligonucleotideFrequency(myseq, width=k, step=1)
  mykmercount <- colSums(myfreq)%>%sort(decreasing=T)
  forgg <- mykmercount%>%as.data.frame%>%rownames_to_column(var="kmer")#%>%mutate(kmer=factor(kmer,levels=kmer))
  return(list(kmer_cnt=forgg,sequences=myseq))
}
forggc <- count_kmers(myxlc,k,mywindow)#Ctrl
forggr <- count_kmers(myxlr,k,mywindow)#RMD

### split by UTR/intron (only pure annotations)
#forggintc=count_kmers(subset(myxlc,region=="intron"),k,mywindow)
#forgg3utrc=count_kmers(subset(myxlc,region=="3'UTR"),k,mywindow)
#forggintr=count_kmers(subset(myxlr,region=="intron"),k,mywindow)
#forgg3utrr=count_kmers(subset(myxlr,region=="3'UTR"),k,mywindow)
forgg5utrc=count_kmers(subset(myxlc,region=="5'UTR"),k,mywindow)
forgg5utrr=count_kmers(subset(myxlr,region=="5'UTR"),k,mywindow)

## unite ctrl, rmd and background
forggall <- forggc$kmer_cnt%>%
  dplyr::full_join(forggr$kmer_cnt,by="kmer")%>%
  dplyr::full_join(forgg,.,by="kmer")%>%set_colnames(c("kmer","UTR","DMSO","RMD"))%>%
  mutate(DMSO_ratio=DMSO/UTR,RMD_ratio=RMD/UTR)%>%arrange(desc(DMSO_ratio),desc(RMD_ratio))

## use the same background count

##intron
forggallint <- forggintc$kmer_cnt%>%
 dplyr::full_join(forggintr$kmer_cnt,by="kmer")%>%
 dplyr::full_join(forgg,.,by="kmer")%>% #kmer_list$intron
  set_colnames(c("kmer","UTR","DMSO","RMD"))%>%
 mutate(DMSO_ratio=DMSO/UTR,RMD_ratio=RMD/UTR)%>%arrange(desc(DMSO_ratio),desc(RMD_ratio))

### 3UTR
forggall3utr <- forgg3utrc$kmer_cnt%>%
 dplyr::full_join(forgg3utrr$kmer_cnt,by="kmer")%>%
 dplyr::full_join(forgg,.,by="kmer")%>% #kmer_list$`3'UTR`
  set_colnames(c("kmer","UTR","DMSO","RMD"))%>%
 mutate(DMSO_ratio=DMSO/UTR,RMD_ratio=RMD/UTR)%>%arrange(desc(DMSO_ratio),desc(RMD_ratio))


##5UTR
forggall5utr <- forgg5utrc$kmer_cnt%>%
  dplyr::full_join(forgg5utrr$kmer_cnt,by="kmer")%>%
  dplyr::full_join(forgg,.,by="kmer")%>% #kmer_list$`3'UTR`
  set_colnames(c("kmer","UTR","DMSO","RMD"))%>%
  mutate(DMSO_ratio=DMSO/UTR,RMD_ratio=RMD/UTR)%>%arrange(desc(DMSO_ratio),desc(RMD_ratio))

## highlight CPE but because we use 6-mers we need all substrings of length 6 that fit into canonical 7-mer UUUUUAU (and others)
allsubstr <- function(x, n) unique(substring(x, 1:(nchar(x) - n + 1), n:nchar(x)))## generate all substrings of length n: stolen https://stackoverflow.com/a/35561930
## CPE consensus: http://genesdev.cshlp.org/content/28/13/1498.full "UUUUA1–2U" or "UUUUUAU" or UUUUAACA https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1817753/ 
##Mendez has even more: TTTTAT’, ‘TTTTAAT’, ‘TTTTAAAT’, ‘TTTTACT’, ‘TTTTCAT’, ‘TTTTAAGT’, ‘TTTTGT #https://www.biorxiv.org/content/10.1101/2021.03.11.434803v1.full
CPEs=c(allsubstr("UUUUUAU",6), allsubstr("UUUUAAU",6), allsubstr("UUUUAAAU",6) ##more canonical
       #,allsubstr("UUUUAACA",6), allsubstr("UUUUACU",6),  allsubstr("UUUUCAU",6), allsubstr("UUUUAAGU",6), allsubstr("UUUUGU",6) ##less canonical with CG
)
forggall$canonical_CPE<-forggall$kmer%in%CPEs
## what to label all CPE plus top nonCPE:
## just plot counts normalized to general UTR abundance
## ranking should be both for dmso and rmd of course
cpes=subset(forggall,canonical_CPE)
topkmers=subset(forggall,!canonical_CPE & (DMSO_ratio > 0.08 & RMD_ratio > 0.04))




plotkmer=function(forggall){
  g=ggplot(forggall,aes(DMSO_ratio,RMD_ratio))+geom_point(col=rgb(.5,.5,.5,.5)) + 
    theme_classic(base_size = global_size) + #rgb(.5,.5,.5,.5) ##eps cannot transparency ##col="grey70",shape=1
    #geom_text_repel(data = topkmers, aes(label=kmer, col="black"),max.overlaps = Inf,force = 4,force_pull=0.05,direction = "both", nudge_x = -0.03) +
    geom_text_repel(data = head(forggall,20), aes(label=kmer, col="black"),max.overlaps = Inf,force = 4,force_pull=0.05,direction = "both", nudge_x = -0.03) +
    geom_text_repel(data = subset(forggall,kmer%in%CPEs), aes(label=kmer, col="darkred"),max.overlaps = Inf,force = 4,force_pull=0.05,direction = "both", nudge_y = -0.01) +
    #geom_point(data = topkmers, aes(col="black")) +
    geom_point(data = subset(forggall,kmer%in%CPEs), aes(col="darkred")) +
    #geom_point(data = subset(forggall,canonical_CPE), col="darkred")+
    labs(x="normalized 6-mer count, DMSO", y="normalized 6-mer count, RMD") + scale_color_manual(values = c("black","darkred")) +
    theme(legend.position = "None") +
    #draw_label(label = "canonical CPEs", x = 7000, y = 500, color = "darkred",hjust = 0,size = 10) +
    theme(aspect.ratio=1) + ##square plot
    NULL
  print(g)
  return(g)
}

ggkmer <- plotkmer(forggall)

  # ggplot(forggall,aes(DMSO_ratio,RMD_ratio))+geom_point(col=rgb(.5,.5,.5,.5)) + 
  #   theme_classic(base_size = global_size) + #rgb(.5,.5,.5,.5) ##eps cannot transparency ##col="grey70",shape=1
  #   #geom_text_repel(data = topkmers, aes(label=kmer, col="black"),max.overlaps = Inf,force = 4,force_pull=0.05,direction = "both", nudge_x = -0.03) +
  #   geom_text_repel(data = head(forggall,20), aes(label=kmer, col="black"),max.overlaps = Inf,force = 4,force_pull=0.05,direction = "both", nudge_x = -0.03) +
  #   geom_text_repel(data = cpes, aes(label=kmer, col="darkred"),max.overlaps = Inf,force = 4,force_pull=0.05,direction = "both", nudge_y = -0.01) +
  #   #geom_point(data = topkmers, aes(col="black")) +
  #   geom_point(data = cpes, aes(col="darkred")) +
  #   #geom_point(data = subset(forggall,canonical_CPE), col="darkred")+
  #   labs(x="normalized 6-mer count, DMSO", y="normalized 6-mer count, RMD") + scale_color_manual(values = c("black","darkred")) +
  #   theme(legend.position = "None") +
  #   #draw_label(label = "canonical CPEs", x = 7000, y = 500, color = "darkred",hjust = 0,size = 10) +
  #   theme(aspect.ratio=1) + ##square plot
  #   NULL
  # 
write_csv(forggall,file="data/source_data/Fig_4D.csv")

save_plot("plots/Fig_4D.pdf",ggkmer)


## add pseudocounts (does not help)

# forggall %<>% 
#   mutate(UTRp=UTR+10, DMSOp=DMSO+10, RMDp=RMD+10) %>% 
#   mutate(DMSO_ratio_p=DMSOp/UTRp,RMD_ratio_p=RMDp/UTRp) 
# 
# ggplot(forggall,aes(DMSO_ratio_p,RMD_ratio_p))+geom_point(col=rgb(.5,.5,.5,.5)) + 
#   theme_classic(base_size = global_size) + 
#   #geom_text_repel(data = topkmers, aes(label=kmer, col="black"),max.overlaps = Inf,force = 4,force_pull=0.05,direction = "both", nudge_x = -0.03) +
#   #geom_text_repel(data = cpes, aes(label=kmer, col="darkred"),max.overlaps = Inf,force = 4,force_pull=0.05,direction = "both", nudge_y = -0.01) +
#   #geom_point(data = topkmers, aes(col="black")) +
#   #geom_point(data = cpes, aes(col="darkred")) +
#   geom_text_repel(data = head(forggall,20), aes(label=kmer, col="black"),max.overlaps = Inf,force = 4,force_pull=0.05,direction = "both", nudge_x = -0.03)  +
#   labs(x="normalized 6-mer count+10, DMSO", y="normalized 6-mer count+10, RMD") +
#   scale_color_manual(values = c("black","darkred")) +
#   theme(legend.position = "None") +
#   theme(aspect.ratio=1) + ##square plot
#   NULL
# 


## ask which sites are with GC (both intron and UTR sites are fine)


plotkmer(forggallint)

plotkmer(forggall3utr)

plotkmer(forggall5utr)


## where the GC comes from - compare ctrl and rmd
# myxl_cond <- plyranges::mutate(width = mywindow, anchor_center(myxlc)) ## make crosslink centered regions of 50nt (xl will be 25)
# myseq <- RNAStringSet(x = getSeq(myxl_cond,x=FaFile(hg19)))
# myfreq <- oligonucleotideFrequency(myseq, width=k, step=1)
# par(mar=c(10,3,3,3))
# myxl_cond[myfreq[,"CGGCGG"]>0] %>% .$region %>% factor %>% summary %>% barplot(las=2)
## Ctrl: intron and 5'UTR, RMD: intron mostly, so these come from the few 5'UTR sites


#################### Fig 4E:  heatmap for top 100 motifs ########################

## run code for 4D first

mytopN=100 ##how many top k-mers (out of 4096 possible 6 mers)
mywindow=50 ## 
## function to make a matrix of kmer position matches
make_kmer_matrix <- function(kmer_cnt, myseqs, mytopN=100, k=6, mywindow=50){
  mykmers <- head(kmer_cnt, mytopN) ##frequency sorted k-mers, take top
  mym <- matrix(0,
                nrow=(mywindow-k+1), ##nt position
                ncol=nrow(mykmers),
                dimnames = list(1:(mywindow-k+1), mykmers$kmer)) ## placeholder for counts per position
  closest_downstream=list() ## list of minimal distance to a downstream match
  for(i in 1:mytopN){ ##for each kmer find its positions across all windows
    myindex  <- Biostrings::vmatchPattern(
      pattern = RNAString(mykmers[i,"kmer"]),
      subject = myseqs) ##creates MIndex object
    ### this is asking for only positive values (only count kmers downstream) and then take shortest distance to XL (this is -26)
    ## like this we have, per cluster, either NA (no match or upstream match including XL itself so U-containting A-mers that XL will be lost or I can use position 21 as a cutoff to catch those) or closest downstream start of the k-mer
    ## position 21 start means : end of the 6-mer is the crosslink site
    ## generalize: center position is window/2+1 
    closest_downstream[[i]]=unlist(lapply(Biostrings::startIndex(myindex), function(x) ifelse(is.null(x),NA,
                                                                                              ifelse(any((x-mywindow/2+1)>0),
                                                                                                     min((x-mywindow/2+1)[x-mywindow/2+1>0]),NA)))
    )
    m <- table(unlist(Biostrings::startIndex(myindex)))%>%data.frame ## this will just count all matches so multiple match per cluster/window is allowed
    mym[m$Var1,mykmers[i,"kmer"]] <- m$Freq
  }
  names(closest_downstream)=mykmers$kmer
  sm.scaled = t(mym)%>%as("ScoreMatrix")%>%scaleScoreMatrix(.)#genomation
  return(list(closest_downstream=closest_downstream,km_matrix=mym, scaled_mat=sm.scaled))
}

## OK we also have A-mers starting in the middle because for clusters without crosslink we simply took this middle.
## we can keep it as it is (and then try to find a U-mer first, and then distance to A-mer in this cluster)
## alternatively we filter only sites that have a crosslink event

kmersc <- make_kmer_matrix(kmer_cnt = forggc$kmer_cnt, myseqs = forggc$sequences, mytopN = mytopN, k = k, mywindow = mywindow)
kmersr <- make_kmer_matrix(kmer_cnt = forggr$kmer_cnt, myseqs = forggr$sequences, mytopN=mytopN, k=k, mywindow = mywindow)


## heatmap and cluster motif logos (use 2 clusters)
myblue=hcl.colors(n=10,palette = "Blues",rev = T)  
set.seed(42) ##for clusters
plot.new()
heatMatrix(kmersc$scaled_mat, xcoords = (1:(mywindow-k+1))-25+2,
           group = as.list(rownames(kmersc$km_matrix)),
           group.col = my_palette[c(1,3)],
           clustfun = function(x) kmeans(x,2)$cluster,
           main=paste("DMSO"), col = myblue ,
           legend.name = "scaled k-mer count",
           xlab = "position around crosslink", grid=T)
heatc <- grid.grab()
plot.new()
set.seed(42) ##for clusters
heatMatrix(kmersr$scaled_mat, xcoords = c(-25, 20),
           group = as.list(rownames(kmersr$km_matrix)),
           group.col = my_palette[c(1,3)],
           clustfun = function(x) kmeans(x,2)$cluster,
           main=paste("RMD"), col = myblue,
           legend.name = "scaled k-mer count",
           xlab = "position around crosslink", grid=T)
heatr <- grid.grab()
plot.new()

## PWMs for clusters
set.seed(42) ##for clusters
kmc <- kmeans(kmersc$scaled_mat,2)
set.seed(42) ##for clusters
kmr <- kmeans(kmersr$scaled_mat,2)

pwmlistc <- lapply(1:2,function(x){
  c <- Biostrings::RNAStringSet(names(kmc$cluster%>%.[.==x]))
  rna <- Biostrings::consensusMatrix(c)[1:4,]
  ##generate an object of the class pcm directly from count matrix
  motif <- new("pcm", mat=as.matrix(rna), name=paste("DMSO cluster",x))
  return(motif)
})

pwmlistr <- lapply(1:2,function(x){
  c <- Biostrings::RNAStringSet(names(kmr$cluster%>%.[.==x]))
  rna <- Biostrings::consensusMatrix(c)[1:4,]
  ##generate an object of the class pcm directly from count matrix
  motif <- new("pcm", mat=as.matrix(rna), name=paste("RMD cluster",x))
  return(motif)
})

## ggplot pwm, provide coordinates in a data frame
df <- data.frame(xmin=c(.05, .05), ymin=c(.7, .2), xmax=c(.25, .25), ymax=c(.9, .4))
df$motif<-pwmlistc
pwmc=ggplot(df, aes(xmin=xmin, ymin=ymin, xmax=xmax, ymax=ymax, motif=motif)) + 
  geom_motif(ic.scale=F) + theme_void() + ylim(0, 1) + xlim(0, .3)

df$motif<-pwmlistr
pwmr=ggplot(df, aes(xmin=xmin, ymin=ymin, xmax=xmax, ymax=ymax, motif=motif)) + 
  geom_motif(ic.scale=F) + theme_void() + ylim(0, 1) + xlim(0, .3)

fullheatc=plot_grid(pwmc,heatc, rel_widths = c(1,3))

### source data provided as unscaled score matrix with cluster number
write.csv(cbind(as.data.frame(kmersc$km_matrix%>%t()),kmc$cluster), file="data/source_data/Fig_4E.csv")

save_plot("plots/Fig_4E.pdf",fullheatc)


################# 4C: site density metaplots over 3'UTR (RCAS package) ##################################

## input : bed and gtf (functions importBed, importGtf)
## create features list which is needed as input to the RCAS function
if(!file.exists("data/txdbFeatures.RData")){
  txdbFeatures <- RCAS::getTxdbFeaturesFromGRanges(gtf)
  save(txdbFeatures, file="data/txdbFeatures.RData")
}
load("data/txdbFeatures.RData")

## only 3UTR is informative so plot only this 
## restrict to 5,000 examples to speed up
# dcvgListc <- RCAS::calculateCoverageProfileList(queryRegions = ctrl2%>% .[grepl("3'UTR",.$region)] %>% .[order(-.$score)] %>% .[1:5000]
# , targetRegionsList = txdbFeatures[c("threeUTRs")], sampleN = 10000)
# dcvgListr <- RCAS::calculateCoverageProfileList(queryRegions = rmd2%>% .[grepl("3'UTR",.$region)] %>% .[order(-.$score)] %>% .[1:5000]
# , targetRegionsList = txdbFeatures[c("threeUTRs")], sampleN = 10000)

## plot windows of the same size around croslsink center instead of actual bed peaks
## the sequence is included in the objects above when I counted k-mers %>% 
## this still does not make it equal - maybe because it randomly selects from more sites?
# dcvgListc <- RCAS::calculateCoverageProfileList(queryRegions = plyranges::mutate(width = mywindow, anchor_center(myxlc)) %>% .[order(-.$score)] %>% .[1:10000]
#                                                 , targetRegionsList = txdbFeatures[c("threeUTRs")], sampleN = 10000)
# dcvgListr <- RCAS::calculateCoverageProfileList(queryRegions = plyranges::mutate(width = mywindow, anchor_center(myxlr)) %>% .[order(-.$score)] %>% .[1:10000]
#                                                 , targetRegionsList = txdbFeatures[c("threeUTRs")], sampleN = 10000)

## only crosslinks (also not equal)
#dcvgListc <- RCAS::calculateCoverageProfileList(queryRegions = myxlc, targetRegionsList = txdbFeatures[c("threeUTRs")], sampleN = 10000)
#dcvgListr <- RCAS::calculateCoverageProfileList(queryRegions = myxlr, targetRegionsList = txdbFeatures[c("threeUTRs")], sampleN = 10000)

## only 3'UTR sites
dcvgListc <- RCAS::calculateCoverageProfileList(queryRegions = plyranges::mutate(width = mywindow, anchor_center(myxlc)) %>% .[grepl("3'UTR",.$region)] %>% .[order(-.$score)] %>% .[1:5000]
                                                , targetRegionsList = txdbFeatures[c("threeUTRs")], sampleN = 10000)
dcvgListr <- RCAS::calculateCoverageProfileList(queryRegions = plyranges::mutate(width = mywindow, anchor_center(myxlr)) %>% .[grepl("3'UTR",.$region)] %>% .[order(-.$score)] %>% .[1:5000]
                                                , targetRegionsList = txdbFeatures[c("threeUTRs")], sampleN = 10000)

dcvgListc$condition="DMSO"
dcvgListr$condition="RMD"

## plot together
rcaspl=ggplot2::ggplot(rbind(dcvgListc,dcvgListr),aes(x = bins, y = meanCoverage)) + 
  geom_ribbon(aes(fill=condition,ymin = meanCoverage - standardError * 1.96, ymax = meanCoverage + standardError * 1.96)) + 
  geom_line(aes(color = condition)) + theme_classic2() +
  labs(y="CPEB4 sites \n mean coverage in 3'UTR",x=" relative 3'UTR length (%)",title="")+
  scale_fill_manual(values = c(alpha("blue",.1),alpha("red",.1)))+
  scale_color_manual(values = c("darkblue","darkred"))+
  theme(aspect.ratio = 1)+
  NULL

write_csv(rbind(dcvgListc,dcvgListr), file="data/source_data/Fig_4C.csv")

save_plot("plots/Fig_4C.pdf",rcaspl)




#################### Fig 4F: distance to other RBPs ##################################

## just use ctrl (DMSO) sample
cpeb4 <- ctrl2
cpeb4$name <- "CPEB4"
##import rest of RBP
helas <- rtracklayer::import("data/beds/all_HeLa_RBPs.bed") ## pre-merged HeLa CLIP peaks
helas$name <- sub("_.+","",helas$name)
## important: remove EWSR because it is still not published
helas%<>%.[.$name!="EWS"]
allrbps <- c(helas[,c("name","thick")], cpeb4[,c("name","thick")])
## import ncbi utrs (gencode is too full of crappy annotations) (done above)
## get center of RBP peaks 
RBPcenters <- GRanges(seqnames = seqnames(allrbps),
                      ranges=allrbps$thick,
                      strand = strand(allrbps),
                      name=allrbps$name)
names(RBPcenters) <- RBPcenters$name
RBPcenters <-unique(RBPcenters) 
## map RBP site centers to UTR
RBPCentersOnUtrs <- mapToTranscripts(x = RBPcenters, transcripts = red_utrs) 
myRBPs=levels(factor(names(RBPCentersOnUtrs)))
## inside RBP centers mapped to 3UTRs, find distance to nearest other side (NA if no other site in the same UTR)
dist2RBP <- function(RBP1,RBP2,RBPonTx){
  return(distanceToNearest(x = GenomicRanges::reduce(subset(RBPonTx, names(RBPonTx)%in%RBP1)),
                           subject = GenomicRanges::reduce(subset(RBPonTx, names(RBPonTx)%in%RBP2)))@elementMetadata$distance)}
# filter sites on last 500nt of 3'UTR
RBPCentersOnUtrs$distance_from_utr_end=width(red_utrs[RBPCentersOnUtrs$transcriptsHits])-RBPCentersOnUtrs@ranges@start
RBPCentFilt <- RBPCentersOnUtrs%>%.[.$distance_from_utr_end<=500] 

## take distances only upstream of PABP (but 3'UTRs are reduced so of course we have internal APASs)
pabp = GenomicRanges::reduce(subset(RBPCentFilt, names(RBPCentFilt)%in%"PABP"))
omnidistUp <- list()
for(i in 1:length(myRBPs)) {
  if( myRBPs[i]=="PABP" ) next  #remove self self distances
  #print(myRBPs[i]) #monitor progress
  myrbp = GenomicRanges::reduce(subset(RBPCentFilt, names(RBPCentFilt)%in%myRBPs[i]))
  ### need to be on the same UTR
  c=subset(myrbp,seqnames(myrbp)%in%seqnames(pabp))
  p=subset(pabp,seqnames(pabp)%in%seqnames(myrbp))
  idx <- precede(c,p) #these are idx of pabp which correcponding myrbp precedes
  #omnidistUp[myRBPs[i]] = distance(c[!is.na(idx)],p[na.omit(idx)])
  if( length(idx)<10 ) next  
  #print(distance(c[!is.na(idx)],p[na.omit(idx)]))
  omnidistUp[[myRBPs[i]]] <- distance(c[!is.na(idx)],p[na.omit(idx)])
}  
omnidistUp %<>% .[lengths(.)>100] ##filter for RBPs with at least n=(100) distances
## distances to PABP
mypabp <- data.frame(Distance=unlist(omnidistUp, use.names = F),
                     RBP=paste(rep(names(omnidistUp),lengths(omnidistUp))))
## add median distance
mypabp1 <- dplyr::inner_join(mypabp,mypabp%>%group_by(RBP)%>%dplyr::summarize(med=median(Distance)))
## plot violins
pabpdistPlot <- 
  ggplot(mypabp1,aes(RBP,Distance))+
  geom_violin(aes(fill=med),trim=T)+
  geom_boxplot(width=.1, outlier.shape = NA)+
  geom_point(aes(RBP,med), shape=21,fill="white",color="black",stroke=.5) +
  theme_classic(base_size=global_size)+
  theme(axis.text.x = element_text(angle=45, hjust = 1),
        legend.title = element_text(size=10,hjust = 0))+
  labs(x="",y="Distance to PABP\n on last 500nt of 3'UTR") +
  scale_fill_gradient(high = "white", low = "coral4",
                      name="median\ndistance (nt)")+
  NULL


## source data
mypabp1%>%rename(MedianDistanceToPABP=med)%>%write_csv(file="data/source_data/Fig_4F.csv")

save_plot("plots/Fig_4F.pdf",pabpdistPlot)

########### Fig 5 ##########

### import raw half lives and connect to target gene names and groups

## import raw half lives, remove ERCCs
hlraw <- read.csv("data/half_lives/hl_norm_to_ERCC92.csv"
                  , stringsAsFactors = F)%>%
  dplyr::select(X,ends_with("_hl"),ends_with("_R2"))%>%
  dplyr::filter(!grepl("^ERCC-",.$X))
## add diagnostic events and target groups
hlraw$DEs <- target_genes$sumDEPerGene[match(hlraw$X,target_genes$gene_names)]
hlraw$target_group <- Hmisc::cut2(x = hlraw$DEs,cuts = mycuts)%>%gsub("\\[| ","",.)%>%gsub("10\\)","9",.)%>%gsub("100\\)","99",.)%>%gsub("1000\\)","999",.)%>%gsub("5000\\)","4999",.)%>%gsub("46027\\]","more",.)%>%gsub("\\,","\\-",.)%>%forcats::fct_explicit_na(., "nontarget")
## add ARE
hlraw$ARE <- ifelse(hlraw$X%in%AREdb$GeneName, "ARE", "non-ARE")
## add IEG
hlraw$IEG <- ifelse(hlraw$X%in%fantom_ieg$Hs_symbol, "IEG", "non-IEG")

### calculate average log FC as average between 2 siRNAs and 2 controls (as Fabian did for thesis)
fab_half <- hlraw%>%
  filter_at(vars(ends_with("_R2")), all_vars(.>=0.5))%>%
  filter_at(vars(ends_with("_hl")), all_vars(.>=0))%>%
  mutate(mean_siCPEB4_DMSO=rowMeans(dplyr::select(., condHeLa_S140_DMSO_rep1_hl, condHeLa_S181_DMSO_rep1_hl)),
         mean_siCPEB4_RMD=rowMeans(dplyr::select(., condHeLa_S140_RMD_rep1_hl, condHeLa_S181_RMD_rep1_hl)),
         mean_Ctrl_DMSO=rowMeans(dplyr::select(., condHeLa_C2_DMSO_rep1_hl, condHeLa_DMSO_rep1_hl)),
         mean_Ctrl_RMD=rowMeans(dplyr::select(., condHeLa_C2_RMD_rep1_hl, condHeLa_RMD_rep1_hl)),
         meanlogFC_DMSO=log2(mean_siCPEB4_DMSO/mean_Ctrl_DMSO),
         meanlogFC_RMD=log2(mean_siCPEB4_RMD/mean_Ctrl_RMD))%>%
  mutate(across(.cols = starts_with("mean_"),.fns = list(log2=log2)))%>% ##for scatterplot need log2
  mutate(stabilized_DMSO=ifelse(IEG=="IEG" & meanlogFC_DMSO>log2(1.8),"IEG > 1.8x stabilized",
                                ifelse(IEG=="IEG" & meanlogFC_DMSO<=log2(1.8), "other IEG","non-IEG")))%>% ##need for highlighting IEG
  mutate(stabilized_RMD=ifelse(IEG=="IEG" & meanlogFC_RMD>log2(1.5),"IEG > 1.5x stabilized",
                               ifelse(IEG=="IEG" & meanlogFC_RMD<=log2(1.5), "other IEG","non-IEG"))) ##need for highlighting IEG

### source data

write_csv(fab_half, file="data/source_data/Fig_5DE_S7C.csv")


##################### 5A: violin plots for half-lives ###################

## split violin plot with only DMSO and RMD 
## select only DMSO and RMD raw condition and filter for positive hl/ R2>0.5  
hl2 <- hlraw%>%dplyr::select(X,target_group,IEG, ARE,
                             condHeLa_DMSO_rep1_hl,condHeLa_RMD_rep1_hl,
                             condHeLa_DMSO_rep1_R2,condHeLa_RMD_rep1_R2)%>%
  filter_at(vars(ends_with("_R2")), all_vars(.>0.5))%>%
  filter_at(vars(ends_with("_hl")), all_vars(.>0))%>%
  dplyr::select(X,ends_with("_hl"),target_group,ARE,IEG)%>%
  reshape2::melt(.)
hl2$target_group %<>% relevel(., ref = "nontarget") 
phl <- 
  ggplot(hl2, aes(target_group, log2(value),fill=variable)) +
  geom_split_violin(scale = "width", trim=F) + 
  scale_fill_brewer(type = "div", direction = -1, labels=c("DMSO","RMD")) +
  labs(y=expression("RNA "*tau[1/2]*" (log"[2]*")"), x="diagnostic events") +
  theme_classic() + 
  theme(axis.text.x = element_text(angle=90),
        legend.title = element_blank()) +
  geom_boxplot(width=0.15, outlier.shape = NA, coef=0) +
  #geom_point(shape=21, aes(fill=variable), color="black", stroke=.5, position=position_jitterdodge(dodge.width = 1.7, jitter.width = .1), alpha=1) + ##https://stackoverflow.com/a/24019668 ## use point because geom_jitter doesn't know fill aes
  NULL


write_csv(hl2, file = "data/source_data/Fig_5A.csv")

save_plot("plots/Fig_5A.pdf",phl)

##################### 5D: ecdf for IEGs #####################

ggsiIEG <- 
  ggplot(fab_half, aes(meanlogFC_DMSO, col=IEG)) +
  stat_ecdf(geom = "step", size=1) +
  labs(title="", y="",x=expression('RNA '*tau[1/2]*' (log'[2]*' fold change siCPEB4/ctrl)') ) +
  theme_classic() +
  theme(legend.position = c(0,1), legend.justification = c(0,1),
        legend.background = element_rect(fill=alpha("grey",0))) +
  scale_color_manual(values = my_palette[c(5,7)],
                     labels=summary(factor(fab_half$IEG))%>%paste(names(.),.,sep=": n=")) +
  NULL
ggsiIEG

save_plot("plots/Fig_5D.pdf",ggsiIEG)

################### 5E: ecdf for target groups #################

ggsimean <- 
  ggplot(fab_half, aes(meanlogFC_DMSO, col=target_group)) +
  stat_ecdf(geom = "step", size=1) +
  labs(title="", y="",x=expression('RNA '*tau[1/2]*' (log'[2]*' fold change siCPEB4/ctrl)') ) +
  theme_classic() +
  theme(legend.position = c(0,1), legend.justification = c(0,1),
        legend.background = element_rect(fill=alpha("grey",0))) +
  scale_color_manual(values = my_palette, name="Target group",
                     labels=summary(factor(fab_half$target_group))%>%paste(names(.),.,sep=": n=")) +
  NULL
ggsimean

save_plot("plots/Fig_5E.pdf",ggsimean)

#################### 5F: GSEA plot ################################

## exclude non ambiguous gene names (when 2 genes share a peak):
csv_uniq <- target_genes%>%.[!grepl(";",.$gene_names),] 
forGsea <- csv_uniq%>%dplyr::select(gene_names,sumDEPerGene)%>%deframe
## first you need to download the GSEA own database for hallmark genes using gene symbols from here https://www.gsea-msigdb.org/gsea/msigdb/collections.jsp#H (login with email)
pathways <- gmtPathways("data/GSEA/h.all.v7.4.symbols.gmt") ## all pathways
## add IEGs from Fantom paper as custom pathway
pathways$FANTOM_IEG <- fantom_ieg$Hs_symbol
set.seed(42)
fgseaRes <- fgsea(pathways, forGsea, maxSize=500, nperm=10000)
leadEdge <- fgseaRes[pathway=="FANTOM_IEG"]$leadingEdge
## modified the original fgsea funciton (change colors)
plotEnrichment_ <- function(pathway, stats, gseaParam = 1, ticksSize = 0.2) 
{
  rnk <- rank(-stats)
  ord <- order(rnk)
  statsAdj <- stats[ord]
  statsAdj <- sign(statsAdj) * (abs(statsAdj)^gseaParam)
  statsAdj <- statsAdj/max(abs(statsAdj))
  pathway <- unname(as.vector(na.omit(match(pathway, names(statsAdj)))))
  pathway <- sort(pathway)
  gseaRes <- calcGseaStat(statsAdj, selectedStats = pathway, 
                          returnAllExtremes = TRUE)
  bottoms <- gseaRes$bottoms
  tops <- gseaRes$tops
  n <- length(statsAdj)
  xs <- as.vector(rbind(pathway - 1, pathway))
  ys <- as.vector(rbind(bottoms, tops))
  toPlot <- data.frame(x = c(0, xs, n + 1), y = c(0, ys, 0))
  diff <- (max(tops) - min(bottoms))/8
  x = y = NULL
  g <- ggplot(toPlot, aes(x = x, y = y)) + geom_point(color = RColorBrewer::brewer.pal(11,"BrBG")[10], 
                                                      size = 0.1) + geom_hline(yintercept = max(tops),
                                                                               colour = RColorBrewer::brewer.pal(11,"BrBG")[2], 
                                                                               linetype = "dashed") + geom_hline(yintercept = min(bottoms), 
                                                                                                                 colour = RColorBrewer::brewer.pal(11,"BrBG")[2], linetype = "dashed") + geom_hline(yintercept = 0, 
                                                                                                                                                                                                    colour = "black") + geom_line(colour = RColorBrewer::brewer.pal(11,"BrBG")[10]) + theme_bw() + 
    geom_segment(data = data.frame(x = pathway), mapping = aes(x = x, 
                                                               y = -diff/2, xend = x, yend = diff/2), size = ticksSize) + 
    theme(panel.border = element_blank(), panel.grid.minor = element_blank()) + 
    labs(x = "rank", y = "enrichment score")
  g
}
ggsea <- 
  plotEnrichment_(pathway = pathways$FANTOM_IEG, stats = forGsea) +
  labs(title="FANTOM IEG",
       subtitle = paste("padj =", round(fgseaRes[pathway=="FANTOM_IEG"]$padj,3), "\ngenes:", unlist(leadEdge)[1:6]%>%paste(collapse = ";"),"..."),
       x="target rank by diagnostic events")


### write source data as GSEA result for all pathways
gsea_df=as.data.frame(fgseaRes[,-8])
gsea_df$leadingEdge=unlist(lapply(fgseaRes$leadingEdge,function(x) paste(x,collapse=",")))

write_csv(gsea_df%>%.[order(.$padj),], file="data/source_data/Fig_5F.csv")

save_plot("plots/Fig_5F.pdf",ggsea)

################################ 5G: example genome browser shots #############################################


#library(Gviz)
options(ucscChromosomeNames=FALSE)

## sequence track - will see T-to-C conversions
#library(BSgenome.Hsapiens.UCSC.hg19) ## to get sequence
sTrack <- SequenceTrack(Hsapiens,fontfamily="Arial",fontfamily.title="Arial")

## clip tracks
cliptrack1 <- AlignmentsTrack(allbams[1], isPaired = F, name = "CLIP r1 p32"
                              #,stacking = "hide"
                              , fontfamily="Arial",fontfamily.title="Arial"
                              , referenceSequence=sTrack
)
cliptrack2 <- AlignmentsTrack(allbams[3], isPaired = F, name = "CLIP r3 ir"
                              , fontfamily="Arial",fontfamily.title="Arial"
                              , referenceSequence=sTrack
)

#rename tracks easily
#cliptrack1@name <- "CLIP r1"
#cliptrack2@name <- "CLIP r2"

## gene region track
## gencode has too many transcripts for plotting, use meta-transcript feature
mytxdb <- loadDb("data/mytxdb.RData")
grtrack <- GeneRegionTrack(
  mytxdb
  , fill="darkblue"
  , collapseTranscripts="longest" ## show "meta" transcript to save space
  ,fontfamily = myfont 
  ,fontfamily.group = myfont 
  ,fontfamily.title = myfont 
)


## need JUNB,FOS,EGR1
# just look up coordinates by hand
mygenes <- data.frame(row.names = c("EGR1","FOS","JUNB"),
                      gene=factor(c("EGR1","FOS","JUNB"),levels=c("EGR1","FOS","JUNB")),
                      chr=c("chr5","chr14","chr19"),
                      afrom=c(137800454,75744708,12901667),
                      ato=c(137805547,75749523,12905204))
igvs <- 
  xyplot(1 ~ gene | gene, data = 
           data.frame(row.names = c("EGR1","FOS","JUNB"),
                      gene=factor(c("EGR1","FOS","JUNB"),levels=c("EGR1","FOS","JUNB")),
                      chr=c("chr5","chr14","chr19"),
                      afrom=c(137800454,75744708,12901667),
                      ato=c(137805547,75749523,12905204)), #par.settings=list(axis.text=list(fontfamily="Comic Sans MS")),
         panel = function(x) {
           afrom=mygenes$afrom[x]
           ato=mygenes$ato[x]
           chr=mygenes$chr[x]
           plotTracks(c(cliptrack1 ,cliptrack2, grtrack),from = afrom,to = ato,chromosome = chr,type="coverage",showId = F, geneSymbol=F,cex = 0.5, sizes = c(1,1,.1),
                      main = paste0(chr,":",afrom,"-",ato), cex.main=1, col.main = "grey30",background.title="transparent",col.title="grey30",col.axis="grey30",add=T)
         }
         , layout=c(3,1)
         , scales = list(draw = FALSE)
         , xlab = NULL, ylab = NULL)
## be careful! if not full screen there wont be enough place and will get weird viewport error


## need font otherwise error saving figure
loadfonts(device = "postscript")
par(family="Arial")
ps.options(family="Arial")
myfont="Arial"

fig5g=plot_grid(igvs)
save_plot("plots/Fig_5G.pdf",fig5g)

########### Fig S5 ##########

################# S5B,C: CLIP gel images #################

# ## add CLIP gel images
# p32 <- "data/gels/p32_clip_fabian_annotated.png"
# p32_panel <- ggdraw() + draw_image(p32, scale = 0.8)
# ir <- "data/gels/20200820-100605_irCLIP_crop.png"
# ir_panel <- ggdraw() + draw_image(ir, scale = .5, x = -.1)
# #annotate irCLIP figure
# ir_anno <-
#   ir_panel +
#   draw_line(x = c(.4,.6),y=.75) +
#   draw_label("3xF-CPEB4", size = 18, x = .5, y = .8) +
#   ## NIR ladder 70, 95, 130, 250 Thermo #26635
#   draw_label("NIR", size = 17, x = .23, y = .75) +
#   draw_label("70", size = 14, x = .1, y = .3) +
#   draw_label("95", size = 14, x = .1, y = .41) +
#   draw_label("130", size = 14, x = .1, y = .52) +
#   draw_label("250", size = 14, x = .1, y = .63) +
#   ## annotate bands with adapters
#   draw_label("CPEB4:RNA:1 adapter", size = 10, x = .67, y = .51, hjust = 0) +
#   draw_line(x = .66, y = c(.54,.48)) +
#   draw_label("CPEB4:RNA:2 adapters", size = 10, x = .67, y = .57, hjust = 0) +
#   draw_line(x = .66, y = c(.55,.59)) +
#   NULL
# 
# options(digits = 3, scipen = -2) #for consistent notation
# blots=plot_grid(p32_panel,ir_anno,labels = "AUTO", scale = .7)
# 
# save_plot("plots/Fig_S5BC.pdf",blots)

#################### S5D: pairs plot for CLIP replicates ###################

## reproducibility of replicates 
## use counts within peaks
## do keep in mind that genes with 0 count in at least one replicate are excluded from each pair
## function to display cor from pairs help:
panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...)
{ usr <- par("usr"); on.exit(par(usr))
par(usr = c(0, 1, 0, 1))
r <- abs(cor(x, y))
txt <- format(c(r, 0.123456789), digits = digits)[1]
txt <- paste0(prefix, txt)
if(missing(cex.cor)) cex.cor <- 0.5/strwidth(txt)
text(0.5, 0.5, txt, cex = cex.cor * r)
}

## log counts and add pseudo
logCnt=log10(counts(ddscl)+1)
colnames(logCnt)%<>%paste0(.,"_log10_p1")

## record plots as objects for ctrl and RMD
plot.new()
pairs(logCnt%>%.[,1:3], lower.panel = panel.smooth, upper.panel = panel.cor,
      gap=0, row1attop=FALSE, labels = c("ir","ir","p32"))
p_repr_c <- recordPlot()
plot.new()
pairs(logCnt%>%.[,4:6], lower.panel = panel.smooth, upper.panel = panel.cor,
      gap=0, row1attop=FALSE, labels = c("ir","ir","p32"))
p_repr_r <- recordPlot() 

write.csv(as.data.frame(logCnt), file = "data/source_data/Fig_S5D.csv")

save_plot("plots/Fig_S5D.pdf",plot_grid(p_repr_c,p_repr_r,ncol=1,scale=.7),base_height = 10, base_width = 6)

########### S5E: TC conversions and mutations barplots ################

## sum diagnostic events: how many T.C and T.deletion compared to other mutations (use txt file from omniCLIP)
plot.new()
barplot(colSums(ctrl[,8:27]),main="DMSO",las=2)
pTC1 <- recordPlot() ## to save plot in an object
plot.new()
barplot(colSums(rmd[,8:27]),main="RMD",las=2)
pTC2 <- recordPlot()
plot.new()

TCpl=plot_grid(pTC1,pTC2,labels = "AUTO", scale=.7)

write.csv(data.frame(DMSO=colSums(ctrl[,8:27]),RMD=colSums(rmd[,8:27])), file="data/source_data/Fig_S5E.csv")

save_plot("plots/Fig_S5E.pdf",TCpl, base_height = 4, base_width = 8)

############# S5F: full transcript categories barplot ###########

##### please run the code for main fig 4B first

## bar plot for the figure
mybarc <- mcols(ctrl2)[,c("tx_region","target_group")]%>%as.data.frame%>%group_by(target_group,tx_region)%>%tally()%>%mutate(sum=sum(n), perc=n/sum, condition="DMSO")
mybarr <- mcols(rmd2)[,c("tx_region","target_group")]%>%as.data.frame%>%group_by(target_group,tx_region)%>%tally()%>%mutate(sum=sum(n), perc=n/sum, condition="RMD")
mybar <- rbind(mybarc,mybarr)
pbar <-  
  ggplot(mybar, aes(x=target_group, y=perc, fill=tx_region))+
  geom_bar(stat="identity", color="grey") +
  theme_classic()+
  scale_fill_manual(values = myCol, name = "")+
  labs(x= "diagnostic events" , y="percentage") +
  theme(axis.text.x = element_text(angle=90))  +
  facet_wrap(~condition) +
  NULL #trick to comment out things


write_csv(mybar, file="data/source_data/Fig_S5F.csv")

save_plot("plots/Fig_S5F.pdf",pbar)

########### Fig S6 ##########

######################## S6A: DESeq foldhcange RNA-seq vs CLIP ##############################

## load previously calculated DEseq objects for RNA and CLIP 
load("data/DEseq.RData") 
## do not use log fold change shrinkage, merge RNA and CLIP fold changes and plot against each other
unshrunk1 <- data.frame(results(dds, name="condition_R_vs_C"))
unshrunkcl1 <- data.frame(results(ddscl, name="condition_R_vs_C"))
m <- merge(data.frame(unshrunk1), data.frame(unshrunkcl1), by=0, all=F)%>%dplyr::rename(gene_name=Row.names)
colnames(m) <- gsub(".y",".CLIP", gsub(".x",".RNAseq", colnames(m)))
mythreshold <- 0.01 # p adj cut off for color highlighting (cut off for CLIP)
m1 <- m%>%dplyr::select(log2FoldChange.RNAseq, gene_name, log2FoldChange.CLIP, padj.CLIP)
m1$sign.CLIP <- !(m1$padj.CLIP>mythreshold|is.na(m1$padj.CLIP)) ##trick to avoid NA
m1$padj.CLIP <- NULL
m1 <- na.omit(m1)
ggscat <- 
  ggplot(m1, aes(log2FoldChange.RNAseq,log2FoldChange.CLIP)) +
  geom_point(aes(color=sign.CLIP)) +
  scale_color_manual(values = c(rgb(.5,.5,.5,.3),"darkred")) +
  labs(title="Genes significantly changing in CLIP"
       ,subtitle = paste("padj (CLIP) <=", mythreshold,"\n","n=",nrow(m1))
       , x=expression("log"[2]*"FC RNAseq (RMD/DMSO)") 
       ,  y=expression("log"[2]*"FC CLIP (RMD/DMSO)") 
  ) +
  geom_text_repel(data=subset(m1, sign.CLIP &  abs(log2FoldChange.CLIP)>4 ),
                  aes(log2FoldChange.RNAseq,log2FoldChange.CLIP, label=gene_name)) +
  theme_classic() + theme(legend.position = "none") +
  geom_abline(lty=2, col="darkred") +
  NULL

write_csv(m1, file="data/source_data/Fig_S6A.csv")

save_plot("plots/Fig_S6A.pdf",ggscat)

##################### S6B: motif heatmap for RMD #######################

## run the code for main Fig 4E first
fullheatr=plot_grid(pwmr,heatr, rel_widths = c(1,3))
write.csv(cbind(as.data.frame(kmersr$km_matrix%>%t()),kmr$cluster), file="data/source_data/Fig_S6B.csv")
save_plot("plots/Fig_S6B.pdf",fullheatr)

########### Fig S7 ##########

############## S7B: boxplot for half-lives DMSO and RMD ####################

## boxplot from Fabian 
## select necessary columns and melt
hl4 <- fab_half%>%
  dplyr::select(X,IEG,starts_with("mean_"))%>%
  dplyr::select(!ends_with("_log2"))%>%
  reshape2::melt(.)
hl4$condition <- hl4$variable%>%sub(".+_","",.)
hl4$siCPEB4 <- hl4$variable%>%sub("mean_","",.)%>%sub("_.+","",.)
pd = position_dodge(width = .9) ## needed to align captions 
## do groups to get the right order 
hl4$mygroup <- paste(hl4$condition,hl4$IEG,hl4$siCPEB4,sep="_")
group_lengths <- summary(factor(hl4$mygroup))
noSam=length(group_lengths) ##number of samples

## add significance stars 
## select only groups that you want to compare
my_comparisons <- list( c("DMSO_IEG_siCPEB4", "DMSO_IEG_Ctrl")
                        ,c("DMSO_non-IEG_siCPEB4","DMSO_non-IEG_Ctrl")
                        ,c("RMD_IEG_siCPEB4", "RMD_IEG_Ctrl")
                        ,c("RMD_non-IEG_siCPEB4","RMD_non-IEG_Ctrl")
)
## shift and y coordinateto position labels and stars:
myshift=.6 
myy=log2(max(hl4$value))+myshift

bpstar <- 
  ggplot(hl4, aes(mygroup, log2(value), fill=IEG))+
  stat_boxplot(geom = 'errorbar', width=.3, position = pd) +
  geom_boxplot(outlier.colour="black", outlier.shape=16, 
               outlier.size=2, notch=FALSE, position = pd) +
  scale_fill_brewer(type = "div") +
  theme_classic(base_size = global_size) +
  labs(y=expression('RNA '*tau[1/2]*' (h, log'[2]*')'), x="siCPEB4", fill="") +
  ## make Fabian-like line annotations over the plot ## depends on group lengths
  annotate(geom = "text", x=c(.25*noSam+.5,.75*noSam+.5), y = myy+6*myshift, label = c("DMSO","RMD"), size = 4) +
  annotate(geom = "line", x = c(1,noSam/2), y = myy+5*myshift) +
  annotate(geom = "line", x = c(noSam/2+1,noSam), y = myy+5*myshift) +
  annotate(geom = "text", x=c(.125*noSam+.5,.375*noSam+.5, .625*noSam+.5, .875*noSam+.5), y = myy+4*myshift, label = rep(c("IEG","non IEG"),2),size = 4) + 
  annotate(geom = "line", x = c(1,noSam/4), y = myy+3*myshift) +
  annotate(geom = "line", x = c(noSam/4+1,noSam/2), y = myy+3*myshift) +
  annotate(geom = "line", x = c(noSam/2+1,0.75*noSam), y = myy+3*myshift) +
  annotate(geom = "line", x = c(0.75*noSam+1,noSam), y = myy+3*myshift) +
  annotate(geom = "text", x=c(.125*noSam+.5,.375*noSam+.5, .625*noSam+.5, .875*noSam+.5), y = myy+2*myshift, label = rep(unique(group_lengths),2), size = 3) +
  scale_x_discrete(labels=names(group_lengths)%>%sub(".+_","",.)%>%gsub("S140","+",.)%>%gsub("C2","-",.)%>%sub("Ctrl","-",.)%>%sub("siCPEB4","+",.)) +
  theme(legend.position = "none") +
  stat_compare_means(comparisons = my_comparisons, label = "p.signif", label.y = myy) +
  NULL

write_csv(hl4, file = "data/source_data/Fig_S7B.csv")

save_plot("plots/Fig_S7B.pdf",bpstar)

################### S7C: ecdf for ARE #########################

## run code for Fig 5 (preface) first
ggsiARE <- 
  ggplot(fab_half, aes(meanlogFC_DMSO, col=ARE)) +
  stat_ecdf(geom = "step", size=1) +
  labs(title="", y="",x=expression('RNA '*tau[1/2]*' (log'[2]*' fold change siCPEB4/ctrl)') ) +
  theme_classic() +
  theme(legend.position = c(0,1), legend.justification = c(0,1),
        legend.background = element_rect(fill=alpha("grey",0))) +
  scale_color_manual(values = my_palette[c(3,7)],
                     labels=summary(factor(fab_half$ARE))%>%paste(names(.),.,sep=": n=")) +
  NULL

save_plot("plots/Fig_S7C.pdf",ggsiARE)





########################reproducible environment ####################################

## print all packages and versions
options(max.print = 1000) ## so that all packages fit
sink("R_sessionInfo.txt")
sessionInfo()
sink()
options(max.print = 50)



