### Figures for Revision

### March 2022

setwd("/ohler/Stoecklin/manuscript_draft/CPEB4_public/")
## where packages are
.libPaths("/ohler/containers/CLIP_Rlibs/")

plotdir <- file.path("plots","revision")
if(! dir.exists(plotdir)){dir.create(plotdir)}

source("load_data.R")

### I am responsible for points raised for Figures 4,5 

## reviewer's (Jack Keene? :-) ) comments and my notes from meeting with Georg

'
Figure 4:
  -Panel C: Suggest making it explicit on the figure that these are crosslink sites
  as opposed to the 6-mers in panel D. Why would there be more sites in DMSO? 
  or is that a reflection of the number of mRNA?
  
  plot top X 1000 each? check how many sites

-in which case add the numbers. Eitherway, as you see, I am not confident I understood this layout.
 
 Maybe it is because there is less RNA in which to detect the sites?

-Panel D: Maybe the GC rich 6-mers deserve a mention

check what those are! just say it is interesting to note (georg want that consensus is different) -
add to discussion



-Panel E: It is slightly unsatisfying to have the plots un-tethered from each other 
so that it is not clear if the same genes that have the CPEB4 motif have the AAAAAA motif.
An additional panel (F) that visually summarises the positioning of the motifs 
with respect to the end of the transcript would be good. Based on panel C,
not all the sites are right at/overlapping the canonical AAUAAA motif
which is normally within ~25 bases of the 3UTR-polyA junction. 
This would be good to see, -given the ideas around the position of these elements 
and functionality in translational activation [ref 96].

can we make this plot for genes? those where it is close vs far

maybe do it for a set of k-mers which are well enriched

how many genes have it plus what is the distance

(have a line per gene or something)

include canonical AAUAAA and distance to 3end

show if the group is close to 3end has better consensus?
  
  the real question is whther they really occur on the same transcript? so we need gene by gene basis

maybe have a color representation of this?
grren and red dots which get close and those which can be further apart

Figure 5:
  Panel G: Suggest putting in a zoomed-in section to highlight the (UUUUAAAA) codes in these transcripts. 
  A quick look in IGV failed to provide a satisfying picture 
  that fits the premise of Fig 4 E for any given transcript and
  that should be discussed so that we are not left 
  with ambiguous literature again regarding the precise positioning of these elements.

zoom in to show motifs

if those do not have we can add another 3 genes which have this

highlight the AAUAAA too

separate the crosslink track
'


## run first functions from manuscript_figures.R to load packages and check that data is there
## to do: (need to put into a separate scirpt)



##### Figure 4 panel C: density of the sites on 3UTR - output how many sites are plotted ####
## OK I already plotted 10k out of 15k RMD and 37k ctrl
## what if I plot all sites or only sites with >10 DEs like for pie charts?
#dcvgListc <- RCAS::calculateCoverageProfileList(queryRegions = ctrl2%>% .[!.$target_group%in%c("0","1-9")] , targetRegionsList = txdbFeatures[c("threeUTRs")])
#dcvgListr <- RCAS::calculateCoverageProfileList(queryRegions = rmd2%>% .[!.$target_group%in%c("0","1-9")] , targetRegionsList = txdbFeatures[c("threeUTRs")])


## I can try to do RBP centers so that less sites fall off the end
## I can also plot top X 1000 sorted by reverse XL events

###make granges for crosslink centers ("thick" position= highest DE per cluster, not filtered out sites with 0 DEs) ##
## try later: make thick2 = sumTC/TD (before it was max between 3 replicates, but I noticed on the igv that many of those are off)
ctrl2$thick2 <- IRanges(start=ifelse(!is.na(ctrl2$thickTCSum),ctrl2$thickTCSum,
                                    ifelse(!is.na(ctrl2$thickTDSum),ctrl2$thickTDSum,
                                           (start(ranges(ctrl2))+round(width(ranges(ctrl2))/2)))),width=1)

rmd2$thick2 <- IRanges(start=ifelse(!is.na(rmd2$thickTCSum),rmd2$thickTCSum,
                                     ifelse(!is.na(rmd2$thickTDSum),rmd2$thickTDSum,
                                            (start(ranges(rmd2))+round(width(ranges(rmd2))/2)))),width=1)

myxlc <- GRanges(seqnames = ctrl2@seqnames,
                 ranges = ctrl2$thick2,
                 strand = ctrl2@strand,
                 score=ctrl2$sumDE,
                 name=ctrl2$name)

myxlr <- GRanges(seqnames = rmd2@seqnames,
                 ranges = rmd2$thick2,
                 strand = rmd2@strand,
                 score=rmd2$sumDE,
                 name=rmd2$name)


#top2k
#dcvgListc <- RCAS::calculateCoverageProfile(queryRegions =  myxlc[order(myxlc$score, decreasing = T)][1:2000], targetRegions = txdbFeatures[[c("threeUTRs")]])
#dcvgListr <- RCAS::calculateCoverageProfileList(queryRegions = myxlr[order(myxlr$score, decreasing = T)][1:2000], targetRegionsList = txdbFeatures[c("threeUTRs")])
## weird! in theory 1811/2000 xl centers are in 3UTR! but this says 1608 windows fall off the target?
## it relies on genomation function genomation::ScoreMatrixBin
## because of error need to filter out UTRs which are shorter than 10bp
## it still removed most of the sites! function constrainRanges (but documentation very poor)

## maybe first filter for those that do overlap 3UTR 

#ctrl2 %>% .[grepl("3'UTR",.$region)] %>% .[order(-.$score)]%>% .[!.$target_group%in%c("0")]
##8k 

#rmd2 %>% .[grepl("3'UTR",.$region)] %>% .[order(-.$score)]%>% .[!.$target_group%in%c("0")]
##5k

## more strict would be to simply subset by overlaps to make sure no sites are outside 3UTR

#ctrl2 %>% subsetByOverlaps(ranges = txdbFeatures$threeUTRs, type = "within") # .[!.$target_group%in%c("0")]


## however genomation scorematrix is originally 
## thought to be used with coverage files such as bam,bw
## so the size of the regions probably matters too
## and I should rather use bam and no bed in this plot

# all sites or sites with cutoff
# because sampleN is random I would rather select top5k regions before this

## try crosslink centers (instead of ctrl2)
## still there is 1000 windows falling off the target even though I require them to be "within" ???

dcvgListc <- calculateCoverageProfileList(queryRegions =  myxlc  %>% subsetByOverlaps(ranges = txdbFeatures$threeUTRs)%>% .[order(-.$score)] %>% .[1:5000]# %>% .[grepl("3'UTR",.$region)] %>% .[order(-.$score)]%>% .[!.$target_group%in%c("0")]#.[!.$target_group%in%c("0"#,"1-9")] 
                                                , targetRegions = txdbFeatures[c("threeUTRs")])
dcvgListr <- calculateCoverageProfileList(queryRegions = myxlr  %>% subsetByOverlaps(ranges = txdbFeatures$threeUTRs)%>% .[order(-.$score)]%>% .[1:5000] #%>% .[grepl("3'UTR",.$region)] %>% .[order(-.$score)]%>% .[!.$target_group%in%c("0")] 
                                                , targetRegionsList = txdbFeatures[c("threeUTRs")])

dcvgListc$condition="DMSO"
dcvgListr$condition="RMD"

## plot together
# mylbl=c(paste("DMSO",length(ctrl2%>% .[!.$target_group%in%c("0"#,"1-9"
#                                                             )] )), ##%>% .[!.$target_group%in%c("0","1-9")] 
#       paste("RMD",length(rmd2 %>% .[!.$target_group%in%c("0"#,"1-9"
#                                                          )] )))

#cutof=1
rcaspl=ggplot2::ggplot(rbind(dcvgListc,dcvgListr),aes(x = bins, y = meanCoverage)) + 
  geom_ribbon(aes(fill=condition,ymin = meanCoverage - standardError * 1.96, ymax = meanCoverage + standardError * 1.96)) + 
  geom_line(aes(color = condition)) + theme_classic2() +
  labs(y="CPEB4 crosslink centers \n mean coverage in 3'UTR",x=" relative 3'UTR length (%)",title="")+
  scale_fill_manual(values = c(alpha("blue",.1),alpha("red",.1))
                    #,labels=mylbl
                    )+
  scale_color_manual(values = c("darkblue","darkred")
                     #,labels=mylbl  
                     )+
  theme(aspect.ratio = 1)+
  NULL
rcaspl



#save_plot(filename = file.path(plotdir, "Fig_5D_allsites.pdf"),rcaspl) ## all sites
#save_plot(filename = file.path(plotdir, "Fig_5D_10DE_sites.pdf"),rcaspl)## only > 10DE sites
#save_plot(filename = file.path(plotdir, "Fig_5D_top5k_sites.pdf"),rcaspl)## top5k select myself
save_plot(filename = file.path(plotdir, "Fig_5D_top5k_c.pdf"),rcaspl)## centers


#### to do: rewrite the funciton by hand?
#### at least make sure to write down how many windows fell off
#### I still do not get why it removes so many windows

### aparently it looks like that DMSO clusters are just "longer" 
## which is probably still related to larger depth of this library
## so that HMM of omniCLIP just had a longer call with those "tails"

width(ctrl2 %>% subsetByOverlaps(ranges = txdbFeatures$threeUTRs)%>% .[order(-.$score)] %>% .[1:5000]) %>% sum
#[1] 278583
width(rmd2 %>% subsetByOverlaps(ranges = txdbFeatures$threeUTRs)%>% .[order(-.$score)] %>% .[1:5000]) %>% sum
#[1] 190183

####### Fig. 4 Panel D: Maybe the GC rich 6-mers deserve a mention  ######

## run the code for panel in manuscript_figures.R

## the thing is that those k-mers are really rare in both UTRs and CLIP
## but the ratio comes up high 

#subset(forggall,DMSO_ratio>0.1)
# kmer   UTR  DMSO  RMD DMSO_ratio  RMD_ratio canonical_CPE
# 1 AAAAAA 90608 11565 5164  0.1276377 0.05699276         FALSE
# 2 AUAAAA 33058  3709 1642  0.1121967 0.04967028         FALSE
# 3 AAUAAA 42818  4548 2092  0.1062170 0.04885796         FALSE
# 4 UAAUAA 19358  2045  885  0.1056411 0.04571753         FALSE
# 5 CGGCGG  2903   305   46  0.1050637 0.01584568         FALSE
# 6 GCGGCG  2925   303   42  0.1035897 0.01435897         FALSE
# 7 UUUGUA 25141  2545 1297  0.1012291 0.05158904         FALSE
# 8 UUUAUA 25103  2512 1229  0.1000677 0.04895829         FALSE

## maybe I can show that they occur in noise sites (<1 DEs)
## this is to check 1 by 1 which sites those are:
#myxl_cond[vcountPattern("CGGCGG", myseq) !=0]

# GRanges object with 142 ranges and 0 metadata columns:
#   seqnames              ranges strand
# <Rle>           <IRanges>  <Rle>
#   [1]     chr1     6845508-6845557      +
#   [2]     chr1   36273779-36273828      +
#   [3]     chr1   36273935-36273984      +
#   [4]     chr1 224622380-224622429      +
#   [5]     chr1   21503040-21503089      -
#   ...      ...                 ...    ...
# [138]    chr13   51484121-51484170      +width(ctrl2 %>% subsetByOverlaps(ranges = txdbFeatures$threeUTRs)%>% .[order(-.$score)] %>% .[1:5000]) %>% sum
# [1] 278583
# > width(rmd2 %>% subsetByOverlaps(ranges = txdbFeatures$threeUTRs)%>% .[order(-.$score)] %>% .[1:5000]) %>% sum
# [1] 190183width(ctrl2 %>% subsetByOverlaps(ranges = txdbFeatures$threeUTRs)%>% .[order(-.$score)] %>% .[1:5000]) %>% sum
# [1] 278583
# > width(rmd2 %>% subsetByOverlaps(ranges = txdbFeatures$threeUTRs)%>% .[order(-.$score)] %>% .[1:5000]) %>% sum
# [1] 190183width(ctrl2 %>% subsetByOverlaps(ranges = txdbFeatures$threeUTRs)%>% .[order(-.$score)] %>% .[1:5000]) %>% sum
# [1] 278583
# > width(rmd2 %>% subsetByOverlaps(ranges = txdbFeatures$threeUTRs)%>% .[order(-.$score)] %>% .[1:5000]) %>% sum
# [1] 190183width(ctrl2 %>% subsetByOverlaps(ranges = txdbFeatures$threeUTRs)%>% .[order(-.$score)] %>% .[1:5000]) %>% sum
# [1] 278583
# > width(rmd2 %>% subsetByOverlaps(ranges = txdbFeatures$threeUTRs)%>% .[order(-.$score)] %>% .[1:5000]) %>% sum
# [1] 190183
#   [139]    chr18       657711-657760      -
#   [140]    chr18   21594070-21594119      -
#   [141]    chr18   45663327-45663376      -
#   [142]    chr21   46913644-46913693      +

## make sure to load 32p CLIP in this thing, seems to be the source of those
## in any case, there should be 0 DEs in CG-rich things

myseq <- RNAStringSet(x = getSeq(red_utrs,x=FaFile(hg19)))
myfreq <- oligonucleotideFrequency(myseq, width=6, step=1) ##6mers
mykmercount <- colSums(myfreq)%>%sort(decreasing=T)
forgg <- mykmercount%>%as.data.frame%>%rownames_to_column(var="kmer")

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

## count k-mers for only >1/>10 DEs
cutof=1 ## this is already enough to remove the GC rich sites
forggc <- count_kmers(myxlc%>% .[.$score>=cutof],k,mywindow)#Ctrl
forggr <- count_kmers(myxlr%>% .[.$score>=cutof],k,mywindow)#RMD
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
#topkmers=subset(forggall,!canonical_CPE & (DMSO_ratio > quantile(forggall$DMSO_ratio,.99) &
#                                             RMD_ratio > forggall$RMD_ratio %>% quantile(.99)))
## just put top 10 or whatever
topk=20
km1=forggall %>% .[!.$canonical_CPE,] %>% .[order(-.$DMSO_ratio),] %>% head(20) %>% .$kmer
km2=forggall %>% .[!.$canonical_CPE,] %>% .[order(-.$RMD_ratio),] %>% head(20) %>% .$kmer
topkmers=subset(forggall,kmer%in%c(km1,km2))
## color A-rich differently
Arich=topkmers[which( (topkmers$kmer %>% str_count(.,"A") ) >=4 ),]
Urich=topkmers[which( (topkmers$kmer %>% str_count(.,"U") ) >=4 ),]

ggkmer <- 
  ggplot(forggall,aes(DMSO_ratio,RMD_ratio))+
  geom_point(col=rgb(.5,.5,.5,.5)) + 
  theme_classic(base_size = global_size) + #rgb(.5,.5,.5,.5) ##eps cannot transparency ##col="grey70",shape=1
  geom_text_repel(data = cpes, aes(label=kmer, col="darkgrey"),max.overlaps = Inf,force = 40,force_pull=0,direction = "both", nudge_y = -0.01) +
  geom_text_repel(data = Urich, aes(label=kmer, col="darkred"),max.overlaps = Inf,force = 40,force_pull=0,direction = "both", nudge_x = -0.03) +
  geom_text_repel(data = Arich, aes(label=kmer, col="darkblue"),max.overlaps = Inf,force = 40,force_pull=0,direction = "both", nudge_x = -0.03) +
  geom_point(data = cpes, aes(col="darkgrey")) +
  geom_point(data = Urich, aes(col="darkred")) +
  geom_point(data = Arich, aes(col="darkblue")) +
  #geom_point(data = subset(forggall,canonical_CPE), col="darkred")+
  labs(title = "k-mer count in all sites >= 1DE"
       ,x="normalized 6-mer count, DMSO"
       , y="normalized 6-mer count, RMD") + 
  scale_color_manual(values = c("darkred","darkgrey","darkblue","black")) +
  theme(legend.position = "None") +
  draw_label(label = "canonical CPEs", x = 0.06, y = 0.005, color = "darkgrey",hjust = 0,size = 10) +
  theme(aspect.ratio=1) + ##square plot
  NULL

ggkmer

save_plot("plots/revision/Fig_4D_1DE_cutoff_thickSum.pdf",ggkmer)

##### 4E: the question is if polyU and polyA are in the same genes or not

## first we need to see if it still holds with our >1DE/>10DE filtered dataset.
## just run the respective portion of manuscript_figures.R keeping in mind that our kmers are now with cutoff
fullheatc=plot_grid(pwmc,heatc, rel_widths = c(1,3))
fullheatr=plot_grid(pwmr,heatr, rel_widths = c(1,3))

save_plot("plots/revision/Fig_4E_DMSO_1DE_thickSum.pdf",fullheatc)
save_plot("plots/revision/Fig_4E_RMD_1DE_thickSum.pdf",fullheatr)

## now we need a heatmap of those kmers in individual genes.. thinking

## OK maybe I can start with plotting sequences of top100 targets directly but including the 3´end

## 0. Also include the 3 genes EGR1, FOS, JUNB, into this plot (OK)
mynames=c("EGR1","FOS") ##"JUNB"is already there

top100c=c(myxlc %>% .[grepl("3'UTR",.$region)] %>% .[order(-.$score)] %>% head(100)
          ,myxlc %>% .[.$name%in%paste0(mynames,"@ctrl")])


#get seq - length too variable, rather get window around center again
#mywindow=50
## OK also I am quite stupid in taking a window of 50 - it would mean the xl site is no.26 on plus strand and no.25 on minus ...
## so need a window of 51 ideally
mywindow=51
top100cseq=RNAStringSet(x = getSeq(plyranges::mutate(width = mywindow, anchor_center(top100c))
  ,x=FaFile(hg19)))
## check that centers are all U
RNAStringSet(x = getSeq(top100c,x=FaFile(hg19))) %>% unique ## should be only U
## now the center is (hopefully) always 26
subseq(top100cseq,26,26) %>% unique()
# always check if it is all U! why this cluster is on + although it is on -?? chr12   9092998
# apparently it is a wrong strand assignment by omniCLIP itself - but how is it possible if all reads map to the minus strand???
# so there is 1 single read on pplus strand and then it got assigned to the overlapping gene 
# cbed[cbed$name=="ENSG00000111752_25202"]
# here are the two genes which overlap, there is cluster both in + and - strand. But: the quesiton is why the + strand one got such a high score. Original score is low which is what it should be.
# so the error must be at the pileup stage

#ctrl_red %>% .[.$gene_names%in%c("M6PR","PHC1")]

# GRanges object with 6 ranges and 5 metadata columns:
#   seqnames          ranges strand |        revmap            orig_names       orig_scores  gene_names        name
# <Rle>       <IRanges>  <Rle> | <IntegerList>           <character>       <character> <character> <character>
#   [1]    chr12 9070724-9070748      + |         26652 ENSG00000111752_25198 0.626907152726154        PHC1        PHC1
# [2]    chr12 9074051-9074095      + |         13837 ENSG00000111752_25199  1.20895457784025        PHC1      PHC1.1
# [3]    chr12 9075541-9075557      + |         34614 ENSG00000111752_25200 0.418340603687857        PHC1      PHC1.2
# [4]    chr12 9076625-9076667      + |         16008 ENSG00000111752_25201  1.04480625622971        PHC1      PHC1.3
# [5]    chr12 9092989-9093009      + |         28014 ENSG00000111752_25202 0.578399822814859        PHC1      PHC1.4
# [6]    chr12 9092960-9093032      - |           298 ENSG00000003056_25203  103.783956926282        M6PR        M6PR

## OK the problem is, the pileup counts reads on all strands, and it does not include the strand into which_label, so the reads from the - strand of the true cluster get also counted for the cluster on the + strand
## that means that overlapping clusters are sharing reads, but the crosslink in the wrong cluster will be on "A"
## however it will still be counted as TC because they are valid AG on the - strand
## I need to "remember" what the strand was for original cluster or try to get unique names in the pileup and only assign TC if cluster strand corresponds to mapped read strand

## maybe I can try karyotype-like plot for top 100
## just use geom_rect (or maybe tile like genomation)
## let us split top kmers into A-rich (>=4A) and U-rich (>=4U)
## I can always expand to top 20 kmers to include polyA or even more kmers


## more convenient to have it as strings

#stringr::str_locate(as.character(top100cseq),"AAUAAA|AAAAAA")
#stringr::str_locate_all(as.character(top100cseq),paste(Arich,collapse="|"))

#colnames(df)=paste0( c(rep("Arich_",2),rep("Urich_",2)),colnames(df))


topn=length(top100c)

df=expand_grid(x=1:topn,y=1:51) ## x gene, y nucleotide
#df$xl=df$x==26 #xl center
df$gene=expand_grid(gene=top100c$name,y=1:51)$gene
df$nt=as.character(top100cseq) %>% strsplit("") %>% unlist #nucleotide

## include the PAS
pas=stringr::str_locate_all(as.character(top100cseq),"AAUAAA")
names(pas)=top100c$name
## OK some have multiple pas
#pas[lengths(pas)>2]
#do.call(rbind, pas) #info is lost here what is going to where

## bad for loop but whatever this is small
#x=pas[lengths(pas)>0][1]
#lapply(pas[lengths(pas)>0],function(x){
#})
#unlist(lapply(pas, function(x) return(x[,"start"])))

pasdf=data.frame(start=unlist(lapply(pas, function(x) return(x[,"start"]))),end=unlist(lapply(pas, function(x) return(x[,"end"]))))
##handle duplicates separately
##then merge to the big df
id1=match(sub(".start|1$|2$","",rownames(pasdf)), df$gene)
#id2=which(df$gene%in%sub(".start|1$","",rownames(pasdf)))#match(df$gene,sub(".start|1$","",rownames(pasdf))) 
#match does not work, only 1st pas returned
#View(subset(df,grepl("RPL34",df$gene)))
pasdf$x=df[id1,]$x

#mypal=RColorBrewer::brewer.pal(11,"BrBG")
mypal=RColorBrewer::brewer.pal(9,"Blues")
#g1=
  ggplot(df)+
  geom_rect(aes(xmin=y-1,xmax=y
                #,ymin=x-1,ymax=x
                ,ymin=topn-x,ymax=topn-x+1
                ,fill=nt)
            ,col="darkgrey")+
  geom_rect(data=pasdf,
            aes(xmin=start-1,xmax=end,ymin=topn-x,ymax=topn-x+1)
            ,fill=rgb(0,0,0,0)
            ,col="darkred")+
  scale_y_discrete(limits=levels(factor(1:topn))
                   ,labels = rev(top100c$name) %>% sub("@ctrl","",.)
                   )+ ## see how to reverse this
  theme(axis.text.y = element_text(size=6))+
  geom_vline(xintercept = c(25,26))+
  labs(x="nucleotide",y="gene", title="sequences of top 100 genes")+
  scale_fill_brewer()+
  # scale_fill_manual(values = c(mypal[2]
  #                              ,mypal[5]
  #                              ,mypal[5]
  #                              ,mypal[8]))+
  NULL
#g1
save_plot(filename = "plots/revision/heatmap_per_gene_sum.pdf",g1)


## 1. add other k-mers (old or new plot)
Arich=topkmers[which( (topkmers$kmer %>% str_count(.,"A") ) >=4 ),]
Urich=topkmers[which( (topkmers$kmer %>% str_count(.,"U") ) >=4 ),]


## include the kmers
## function returns dataframe to be input into ggplot
locate_kmers=function(kmer_regex="AAUAAA",top100cseq,top100c){
  pas=stringr::str_locate_all(as.character(top100cseq),kmer_regex)
  names(pas)=top100c$name
  pasdf=data.frame(start=unlist(lapply(pas, function(x) return(x[,"start"])))
                   ,end=unlist(lapply(pas, function(x) return(x[,"end"]))))
  ## anywhere with more than 1 match there is numbering just after the name
  id1=match(sub(".start|[0-9]+$","",rownames(pasdf)), df$gene)
  pasdf$x=df[id1,]$x
  return(pasdf)
}
##locate A-rich and U-rich kmers
adf=locate_kmers(kmer_regex=paste0(Arich$kmer,collapse = "|")
                                   ,top100c=top100c
                                   ,top100cseq=top100cseq)
udf=locate_kmers(kmer_regex=paste0(Urich$kmer,collapse = "|")
                 ,top100c=top100c
                 ,top100cseq=top100cseq)


##Plot

ggplot(df)+
  geom_rect(aes(xmin=y-1,xmax=y
                ,ymin=topn-x,ymax=topn-x+1)
            ,fill="darkgrey"
            #,col="darkgrey"
            )+
  ## A rich
  geom_rect(data=adf,
            aes(xmin=start-1,xmax=end,ymin=topn-x,ymax=topn-x+1
            ,fill="A-rich"#rgb(0,0,0,0)
            #,col="A-rich"
            ))+
  ## U rich
  geom_rect(data=udf,
            aes(xmin=start-1,xmax=end,ymin=topn-x,ymax=topn-x+1
            ,fill="U-rich"#rgb(0,0,0,0)
            #,col="U-rich"
            ))+
  scale_y_discrete(limits=levels(factor(1:topn))
                   ,labels = rev(top100c$name) %>% sub("@ctrl","",.)
                   )+ 
  theme(axis.text.y = element_text(size=6))+
  geom_vline(xintercept = c(25,26))+
  labs(x="nucleotide",y="gene"
       , title="kmer positions of top 100 genes in 3'UTRs"
       , fill="6-mers",colour="6-mers")+
  scale_fill_manual(values = c("darkred","darkcyan")
                    ,labels=c("A-rich","U-rich"))+
  # scale_color_manual(values = c("darkred","darkcyan")
  #                 ,labels=c("A-rich","U-rich"))+
  NULL

## 2. prolong sequence to the annotated 3'end
## for this, we actually need 3'UTR sites (some in top100 are 5'UTR) ## done
## let us take any closest, same strand, annotated 3' end
## this file ncbi_utrs is not good, need new UTRs
names(ncbi_utrs)=ncbi_utrs$name
## OK there are apparently huge problems with this ones, too (many exons on diff. chromosomes)
#refseq=makeTxDbFromGFF(file = "annotation/table_browser_Genes_hg19_NCBI_RefSeq.gtf.gz")
#saveDb(refseq,file="data/RefSeq_txdb.RData")
#refseq=loadDb("data/RefSeq_txdb.RData")

## use gencode because "SIK1.2@ctrl" could not be mapped to NCBI UTRs
## but those UTRs are super long!
mytxdb=loadDb("data/mytxdb.RData")
utr3s=threeUTRsByTranscript(mytxdb, use.names=T)


#ends_of_3utr=plyranges::mutate(width = 1, anchor_3p(utr3s))
#xl2utr=plyranges::join_precede_downstream(x=top100c,y=ends_of_3utr) ## meh that returns >1 end per site
#top100ends=precede(x=top100c,subject = ends_of_3utr,select="first",ignore.strand=F)
#myutrs=ncbi_utrs[top100ends] #ncbi_utrs[xl2utr$name.y]
##take 25 nt upstream and the 3UTR end downstream
#top100utr=plyranges::flank_upstream(top100c,25)
## replace end by the UTR end by hand (did not find function to do this)
# for(i in 1:length(top100utr)){
#   if(strand(top100utr[i])@values=="+"){
#     end(ranges(top100utr[i]))=end(ranges(myutrs[i]))
#   }
#   if(strand(top100utr[i])@values=="-"){
#     start(ranges(top100utr[i]))=start(ranges(myutrs[i]))
#   }
# }
# 
# ## check on browser
# export.bed(top100utr, "xl2utr.bed")
mynames=c("EGR1", "FOS")
## can try to cut off on distance to the UTR end even though it is quite a cheating 
## (and even internal sites have pA so then at least split the two)
## or plot only IEGs - they are guaranteed to be short (no they are not! )

#IEG=myxlc %>% .[(.$name %>% sub("@.+","",.))%in%fantom_ieg$Hs_symbol]

## use red_utrs or non-reduced ncbi_utrs or gencode
## take sites with at least 10 DEs and plot separately those close vs those far away
top100c=c(
  (myxlc %>% subsetByOverlaps(utr3s) %>% .[order(-.$score)] %>% .[.$score>=10])
          ,myxlc %>% .[.$name%in%paste0(mynames,"@ctrl")])
## use UTRs to which I mapped distance (reduced_utrs)
#top100map=mapToTranscripts(top100c, transcripts = ncbi_utrs)
## map crosslink sites to gencode UTRs
top100map=mapToTranscripts(top100c, transcripts = utr3s)

utrseq=extractTranscriptSeqs(x= FaFile(hg19), transcripts = utr3s[top100map$transcriptsHits])
#utrseq=getSeq(x=FaFile(hg19),ncbi_utrs[top100map$transcriptsHits])

mapped=data.frame(
  crosslink=top100c[top100map$xHits]$name
  ,score=top100c[top100map$xHits]$score
  ,txid=names(utrseq)
  ,utrlen=width(utrseq)
  ,utrseq=as.character(utrseq)
  ,map_pos=ranges(top100map)@start
  ,stringsAsFactors = F
) %>% dplyr::mutate(dis_to_end=utrlen-map_pos) #%>% 

## define universal offset (all of the sequences start this amount of nt before XL)
offset=min(25,mapped$map_pos)

## this selects closest 3'UTR end
seq2plot=mapped %>%   dplyr::group_by(crosslink) %>% 
  dplyr::summarize(min_dis=min(dis_to_end)
                   ,map_pos=map_pos[which.min(dis_to_end)]
                   ,utr=txid[which.min(dis_to_end)]
                   #,utrseq=utrseq[which.min(dis_to_end)]
                   ,score=score[which.min(dis_to_end)]
                   ,substr=substr(utrseq[which.min(dis_to_end)]
                                  ,start=map_pos-offset
                                  ,stop=utrlen[which.min(dis_to_end)]) %>% gsub("T","U",.))


## take subset of those which are close (1st quantile < ~50nt to end) and far (> median=241 < 1000 3rd q.)
close=subset(seq2plot, min_dis<=quantile(seq2plot$min_dis)[2]) %>% .[order(-.$score),] %>% head(100)

farther=subset(seq2plot, min_dis>quantile(seq2plot$min_dis)[2] & 
                 min_dis<=quantile(seq2plot$min_dis)[3]) %>% .[order(-.$score),] %>% head(100)

evenfarther=subset(seq2plot, min_dis>quantile(seq2plot$min_dis)[3] & 
                 min_dis<=quantile(seq2plot$min_dis)[4]) %>% .[order(-.$score),] %>% head(100)

far=subset(seq2plot, min_dis>quantile(seq2plot$min_dis)[4] ) %>% .[order(-.$score),] %>% head(100)

##
# top100cseq=RNAStringSet(x = getSeq(top100utr,x=FaFile(hg19)))
# ## check that centers are all U
# RNAStringSet(x = getSeq(top100c,x=FaFile(hg19))) %>% unique ## should be only U
# ## now the center is (hopefully) always 26
# subseq(top100cseq,26,26) %>% unique()

plotUTR=function(seq2plot,offset){
  
  ## restore order
  #seq2plot%<>%.[match(top100c$name,.$crosslink),] ##by score
  seq2plot %<>% .[order(.$min_dis),] ## by distance
  
  ## one could not be mapped (exclude by filtering 3UTR only sites before)
  #seq2plot%<>%na.omit(.)
  
  ##Plot
  ## make new df with individual lengths
  ## x is ordinal to count genes, y to count nucleotides in each sequence
  mylen=stringr::str_length(seq2plot$substr)
  df=data.frame()
  for(i in 1:nrow(seq2plot)){
    df=rbind(df,expand_grid(x=i,y=1:mylen[i])) ## x gene, y nucleotide
  }
  mygenes=data.frame()
  for(i in 1:nrow(seq2plot)){
    mygenes=rbind(mygenes,expand_grid(gene=seq2plot$crosslink[i],y=1:mylen[i]))
  }
  df$gene=mygenes$gene %>% sub("@.+","",.)
  #df$nt=as.character(top100cseq) %>% strsplit("") %>% unlist #nucleotide
  
  ##locate A-rich and U-rich kmers
  locate_kmers=function(kmer_regex="AAUAAA",sequences,names,df){
    pas=stringr::str_locate_all(sequences,kmer_regex)
    names(pas)=names
    pasdf=data.frame(start=unlist(lapply(pas, function(x) return(x[,"start"])))
                     ,end=unlist(lapply(pas, function(x) return(x[,"end"]))))
    ## anywhere with more than 1 match there is numbering just after the name
    id1=match(sub("@.+","",rownames(pasdf)), df$gene)
    pasdf$x=df[id1,]$x
    return(pasdf)
  }
  adf=locate_kmers(kmer_regex=paste0(Arich$kmer,collapse = "|")
                   ,sequences = seq2plot$substr
                   ,names=seq2plot$crosslink, df=df)
  udf=locate_kmers(kmer_regex=paste0(Urich$kmer,collapse = "|")
                   ,sequences = seq2plot$substr
                   ,names=seq2plot$crosslink,df=df)
  
  topn=df$gene %>% unique() %>% length()
  g2=ggplot(df)+
    geom_rect(aes(xmin=y-1,xmax=y
                  ,ymin=topn-x,ymax=topn-x+1)
              ,fill="darkgrey"
              #,col="darkgrey"
    )+
    ## A rich
    geom_rect(data=adf,
              aes(xmin=start-1,xmax=end,ymin=topn-x,ymax=topn-x+1
                  ,fill="A-rich"#rgb(0,0,0,0)
                  #,col="A-rich"
              ))+
    ## U rich
    geom_rect(data=udf,
              aes(xmin=start-1,xmax=end,ymin=topn-x,ymax=topn-x+1
                  ,fill="U-rich"#rgb(0,0,0,0)
                  #,col="U-rich"
              ))+
    scale_y_discrete(limits=levels(factor(1:topn))
                     ,labels = rev(unique(df$gene))
    )+ 
    theme(axis.text.y = element_text(size=6))+
    geom_vline(xintercept = c(offset-1,offset))+
    labs(x="nucleotide",y="gene"
         , title="kmer positions of top 100 genes in 3'UTRs"
         , fill="6-mers",colour="6-mers")+
    scale_fill_manual(values = c("darkred","darkcyan")
                      ,labels=c("A-rich","U-rich"))+
    # scale_color_manual(values = c("darkred","darkcyan")
    #                    ,labels=c("A-rich","U-rich"))+
    NULL
  print(g2)
  return(g2)
}

g1=plotUTR(seq2plot = close, offset)
g2=plotUTR(seq2plot = farther, offset)
g2+xlim(c(0,50))
g3=plotUTR(seq2plot = evenfarther, offset)
g3+xlim(c(0,50))
g4=plotUTR(seq2plot = far, offset)
g4+xlim(c(0,50))



## IEG UTRs are not that shorter after all! 
