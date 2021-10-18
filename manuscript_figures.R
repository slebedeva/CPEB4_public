#### make figures for CPEB4 manuscript

## working dir is where this source file is
basedir <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(basedir)
## put plots separately in a folder
plotdir <- file.path(basedir,"plots")
if(! dir.exists(plotdir)){dir.create(plotdir)}

## additional funcitons
source("CLIP.Rprofile")

## library path
# for containers
.libPaths("CLIP_Rlibs") ## local path depends on mountpoint inside container!

###################################################### load environment ########################################################################################

## important! otherwise Hs.org.db fails ##https://support.bioconductor.org/p/9136239/#9136319
options(connectionObserver = NULL)

mypackages <- c("plyranges"
                ,"nVennR"
                ,"tidyverse"
                ,"reshape2"
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
                #,"TxDb.Hsapiens.UCSC.hg19.knownGene"
                #,"Homo.sapiens"
                ,"genomation"
                #,"RCAS"
                ,"SRAdb"
                ,"fgsea"
                ,"cowplot"
                ,"gridGraphics"
                ,"magick"
                ,"rio"
                )
suppressPackageStartupMessages(lapply(mypackages, require, character.only=T))

############## load data ##############

## define cut offs for target groups by diagnostic events
mycuts=c(0,1,10,100,1000,5000)

## get genome sequence (only once)
if(!file.exists("GRCh37.p13.genome.fa")|file.exists("GRCh37.p13.genome.fa.gz")|file.exists("GRCh37.p13.genome.fa.bgz")){
  system("wget http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/GRCh37.p13.genome.fa.gz")
  system("gunzip -c GRCh37.p13.genome.fa.gz | bgzip  > GRCh37.p13.genome.fa.bgz")
  }
hg19=("GRCh37.p13.genome.fa.bgz")

## annotation
load("data/gtf.RData")

## transcript reigons
load("data/myregions.RData")

## merged CLIP peaks
load("data/bed_reduced.RData") ## these peaks are reduced (overlapping peak calls are merged)
load("data/reduced_beds_for_igv.RData") ## these peaks have transcript regions assigned

## differential expresison
load("DEseq.RData") 


## bam file list
allbams <- list.files("data/bams")

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
system("rm wget-log")


########### Figures ################

########### Fig 4 ##########

## B: pie charts transcript categories for sites > 10 DEs



################# 4C: site density metaplots over 3'UTR (RCAS package) ##################################

## input : bed and gtf (functions importBed, importGtf)
#txdbFeatures <- RCAS::getTxdbFeaturesFromGRanges(gtf)
#save(txdbFeatures, file="data/txdbFeatures.RData")
load("data/txdbFeatures.RData")

## only 3UTR is informative so plot only this ## restrict to 10,000 examples to speed up
dcvgListc <- RCAS::calculateCoverageProfileList(queryRegions = ctrl2, targetRegionsList = txdbFeatures[c("threeUTRs")], sampleN = 10000)
dcvgListr <- RCAS::calculateCoverageProfileList(queryRegions = rmd2, targetRegionsList = txdbFeatures[c("threeUTRs")], sampleN = 10000)

dcvgListc$condition="DMSO"
dcvgListr$condition="RMD"

# ## different plots separately
# ggplot2::ggplot(dcvgListc, aes(x = bins, y = meanCoverage)) + 
#   geom_ribbon(fill=alpha("blue",.1),aes(ymin = meanCoverage - standardError * 1.96, ymax = meanCoverage + standardError * 1.96)) + 
#   geom_line(color = 'darkblue') + theme_classic2() +   labs(y="Mean coverage",x="3'UTR",title="DMSO")
# ggplot2::ggplot(dcvgListr, aes(x = bins, y = meanCoverage)) + 
#   geom_ribbon(fill=alpha("red",.1),aes(ymin = meanCoverage - standardError * 1.96, ymax = meanCoverage + standardError * 1.96)) + 
#   geom_line(color = 'darkred') + theme_classic2() +   labs(y="Mean coverage",x="3'UTR",title="RMD")
## together
ggplot2::ggplot(rbind(dcvgListc,dcvgListr),aes(x = bins, y = meanCoverage)) + 
  geom_ribbon(aes(fill=condition,ymin = meanCoverage - standardError * 1.96, ymax = meanCoverage + standardError * 1.96)) + 
  geom_line(aes(color = condition)) + theme_classic2() +
  labs(y="CPEB4 sites \n mean coverage in 3'UTR",x=" relative 3'UTR length (%)",title="")+
  scale_fill_manual(values = c(alpha("blue",.1),alpha("red",.1)))+
  scale_color_manual(values = c("darkblue","darkred"))+
  theme(aspect.ratio = 1)+
  NULL

write_csv(rbind(dcvgListc,dcvgListr), file="data/source_data/Fig_4C.csv")

################# 4D: scatterplot of motif counts #################################

### calculate background: what is kmer count in all UTRs?
## use reduced UTRs to avoid double counting same sequence
## make red utrs 
ncbi_utrs <- rtracklayer::import("ncbi_refseq_utrs.bed") ## downloaded from ucsc table browser
red_utrs <- GenomicRanges::reduce(ncbi_utrs,with.revmap=T) ##names needed for mapToTx
canonical_chr <- paste0("chr",c(1:22,"X","Y","M"))
red_utrs <- red_utrs%>%subset(., .@seqnames%in%canonical_chr)
red_utrs$name <- ncbi_utrs$name[unlist(lapply(red_utrs$revmap, function(x) return(x[[1]])))]
names(red_utrs) <- red_utrs$name

##
myseq <- RNAStringSet(x = getSeq(red_utrs,x=FaFile(hg19)))
myfreq <- oligonucleotideFrequency(myseq, width=6, step=1) ##6mers
mykmercount <- colSums(myfreq)%>%sort(decreasing=T)
forgg <- mykmercount%>%as.data.frame%>%rownames_to_column(var="kmer")

## can also count in random 50 nt windows in UTRs

###make granges for crosslink centers ("thick" position= highest DE per cluster, not filtered out sites with 0 DEs) ##
myxlc <- GRanges(seqnames = ctrl2@seqnames,
                 ranges = ctrl2$thick,
                 strand = ctrl2@strand)
myxlr <- GRanges(seqnames = rmd2@seqnames,
                 ranges = rmd2$thick,
                 strand = rmd2@strand)


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
## unite ctrl, rmd and background
forggall <- forggc$kmer_cnt%>%
  dplyr::full_join(forggr$kmer_cnt,by="kmer")%>%
  dplyr::full_join(forgg,.,by="kmer")%>%set_colnames(c("kmer","UTR","DMSO","RMD"))%>%
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
#ggkmer <- 
ggplot(forggall,aes(DMSO_ratio,RMD_ratio))+geom_point(col=rgb(.5,.5,.5,.5)) + theme_classic(base_size = global_size) + #rgb(.5,.5,.5,.5) ##eps cannot transparency ##col="grey70",shape=1
  geom_text_repel(data = topkmers, aes(label=kmer, col="black"),max.overlaps = Inf,force = 4,force_pull=0.05,direction = "both", nudge_x = -0.03) +
  geom_text_repel(data = cpes, aes(label=kmer, col="darkred"),max.overlaps = Inf,force = 4,force_pull=0.05,direction = "both", nudge_y = -0.01) +
  geom_point(data = topkmers, aes(col="black")) +
  geom_point(data = cpes, aes(col="darkred")) +
  #geom_point(data = subset(forggall,canonical_CPE), col="darkred")+
  labs(x="normalized 6-mer count, DMSO", y="normalized 6-mer count, RMD") + scale_color_manual(values = c("black","darkred")) +
  theme(legend.position = "None") +
  draw_label(label = "canonical CPEs", x = 7000, y = 500, color = "darkred",hjust = 0,size = 10) +
  theme(aspect.ratio=1) + ##square plot
  NULL

write_csv(forggall,file="data/source_data/Fig_4D.csv")

#################### Fig 4E:  heatmap for top 100 motifs ########################

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

plot_grid(pwmc,heatc, rel_widths = c(1,3))

### source data provided as unscaled score matrix with cluster number
write_csv(cbind(as.data.frame(kmersc$km_matrix%>%t()),kmc$cluster), file="data/source_data/Fig_4E.csv")


#################### Fig 4F: distance to other RBPs ##################################

## just use ctrl (DMSO) sample
cpeb4 <- ctrl2
cpeb4$name <- "CPEB4"
##import rest of RBP
helas <- rtracklayer::import("data/allrbps.bed") ## pre-merged HeLa CLIP peaks
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
  print(myRBPs[i]) #monitor progress
  myrbp = GenomicRanges::reduce(subset(RBPCentFilt, names(RBPCentFilt)%in%myRBPs[i]))
  ### need to be on the same UTR
  c=subset(myrbp,seqnames(myrbp)%in%seqnames(pabp))
  p=subset(pabp,seqnames(pabp)%in%seqnames(myrbp))
  idx <- precede(c,p) #these are idx of pabp which correcponding myrbp precedes
  #omnidistUp[myRBPs[i]] = distance(c[!is.na(idx)],p[na.omit(idx)])
  if( length(idx)<10 ) next  
  print(distance(c[!is.na(idx)],p[na.omit(idx)]))
  omnidistUp[[myRBPs[i]]] <- distance(c[!is.na(idx)],p[na.omit(idx)])
}  
omnidistUp %<>% .[lengths(.)>100] ##filter for RBPs with at least n=(100) distances
## distances to PABP
mypabp <- data.frame(Distance=unlist(omnidistUp, use.names = F),
                     RBP=paste(rep(names(omnidistUp),lengths(omnidistUp))))
## add median distance
mypabp1 <- inner_join(mypabp,mypabp%>%group_by(RBP)%>%summarize(med=median(Distance)))
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
pabpdistPlot

## source data
mypabp1%>%rename(MedianDistanceToPABP=med)%>%write_csv(file="data/source_data/Fig_4F.csv")

########### Fig 5 ##########

## A: violin plots for half-lives

## D: ecdf for IEGs

## E: ecdf for target groups

## F: GSEA plot

## G: example genome browser shots


########### Fig S5 ##########

################# S5B,C: CLIP gel images #################

## add CLIP gel images
p32 <- "data/gels/p32_clip_fabian_annotated.png"
p32_panel <- ggdraw() + draw_image(p32, scale = 0.8)
ir <- "data/gels/20200820-100605_irCLIP_crop.png"
ir_panel <- ggdraw() + draw_image(ir, scale = .5, x = -.1)
#annotate irCLIP figure
ir_anno <- 
  ir_panel +
  draw_line(x = c(.4,.6),y=.75) +
  draw_label("3xF-CPEB4", size = 18, x = .5, y = .8) +
  ## NIR ladder 70, 95, 130, 250 Thermo #26635
  draw_label("NIR", size = 17, x = .23, y = .75) +
  draw_label("70", size = 14, x = .1, y = .3) +
  draw_label("95", size = 14, x = .1, y = .41) +
  draw_label("130", size = 14, x = .1, y = .52) +
  draw_label("250", size = 14, x = .1, y = .63) +
  ## annotate bands with adapters
  draw_label("CPEB4:RNA:1 adapter", size = 10, x = .67, y = .51, hjust = 0) +
  draw_line(x = .66, y = c(.54,.48)) +
  draw_label("CPEB4:RNA:2 adapters", size = 10, x = .67, y = .57, hjust = 0) +
  draw_line(x = .66, y = c(.55,.59)) +
  NULL

options(digits = 3, scipen = -2) #for consistent notation
plot_grid(p32_panel,ir_anno,labels = "AUTO", scale = .7)

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

write_csv(as.data.frame(logCnt), file = "data/source_data/Fig_S5D.csv")


########### S5E: TC conversions and mutations barplots ################

## sum diagnostic events: how many T.C and T.deletion compared to other mutations (use txt file from omniCLIP)
plot.new()
barplot(colSums(ctrl[,8:27]),main="DMSO",las=2)
pTC1 <- recordPlot() ## to save plot in an object
plot.new()
barplot(colSums(rmd[,8:27]),main="RMD",las=2)
pTC2 <- recordPlot()
plot.new()

plot_grid(pTC1,pTC2,labels = "AUTO", scale=.7)
write_csv(data.frame(DMSO=colSums(ctrl[,8:27]),RMD=colSums(rmd[,8:27])), file="data/source_data/Fig_S5E.csv")

############# S5F: full transcript categories barplot ###########

## same color scheme everywhere - put all shared into "multiple" category otherwise too many colors!
ctrl2$tx_region <- ifelse(!grepl("_",ctrl2$region),ctrl2$region,
                          ifelse(grepl("3'UTR",ctrl2$region),"3UTR or other","multiple non-3'UTR"))
rmd2$tx_region <- ifelse(!grepl("_",rmd2$region),rmd2$region,
                         ifelse(grepl("3'UTR",rmd2$region),"3UTR or other","multiple non-3'UTR"))
allreg <- levels(factor(c(ctrl2$tx_region,rmd2$tx_region)))
myCol <- rev(RColorBrewer::brewer.pal(n = length(allreg), "BrBG"))
names(myCol) <- allreg

## bar plot for the figure
ctrl2$target_group %<>% gsub("5000-38127","5000-more", .)
rmd2$target_group %<>% gsub("5000-9399","5000-more",.)
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
pbar

write_csv(mybar, file="data/source_data/Fig_S5F.csv")


########### Fig S6 ##########

######################## S6A: DESeq foldhcange RNA-seq vs CLIP ##############################

## load previously calculated DEseq objects for RNA and CLIP 
load("DEseq.RData") 
## do not use log fold change shrinkage, merge RNA and CLIP fold changes and plot against each other
unshrunk1 <- data.frame(results(dds, name="condition_R_vs_C"))
unshrunkcl1 <- data.frame(results(ddscl, name="condition_RMD_CLIP_vs_Ctrl_CLIP"))
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
ggscat

write_csv(m1, file="data/source_data/Fig_S6A.csv")


##################### S6B: motif heatmap for RMD #######################

## run code for main Fig 4E first
plot_grid(pwmr,heatr, rel_widths = c(1,3))
write_csv(cbind(as.data.frame(kmersr$km_matrix%>%t()),kmr$cluster), file="data/source_data/Fig_S6B.csv")

########### Fig S7 ##########

## B: boxplot for half-lives DMSO and RMD

## C: ecdf for ARE










########################reproducible environment ####################################

sink(paste0(Sys.Date(),"_manuscript_figures_sessionInfo.txt"))
sessionInfo()
sink()



