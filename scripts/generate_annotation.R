### Analyze omniCLIP peaks for CPEB4

## This script will generate annotation objects necessary for the analysis

## Only use it within manuscript_figures.R otherwise all paths will be broken!
if(!exists("basedir")){message("Please source me within manuscript_figures.R")}
message("Generating annotation... This will take some time.")

###################################################### load environment ########################################################################################

## important! otherwise Hs.org.db fails ##https://support.bioconductor.org/p/9136239/#9136319
options(connectionObserver = NULL)
## where packages are
.libPaths("/ohler/containers/CLIP_Rlibs/")

mypackages <- c("config"
                ,"plyranges"
                ,"tidyverse"
                ,"reshape2"
                ,"magrittr"
                ,"Biostrings"
                ,"GenomicRanges"
                ,"GenomicFeatures"
                ,"rtracklayer"
                ,"Rsamtools"
                ,"rio"
                ,"data.table"
)

pl=suppressWarnings(suppressPackageStartupMessages(lapply(mypackages, require, character.only=T)))
if(sum(!unlist(pl))){message("failed to load: \n",paste(mypackages[!unlist(pl)],collapse = "\n"), "\nanyway proceeding...")}


############## load data ##############

config=config::get()

## get genome sequence (only once)
hg19=config$genome

if(!(file.exists(hg19))){
  message("Downloading hg19 genome....")
  system("wget http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/GRCh37.p13.genome.fa.gz")
  system("gunzip -c GRCh37.p13.genome.fa.gz | bgzip  > GRCh37.p13.genome.fa.bgz")
}


####### generate gencode v19 annotation data #################

## get gtf (only once)
if(!file.exists(file.path(ann_dir,"gencode.v19.annotation.gtf.gz"))){
  message("Downloading Gencode 19 annotation...")
  system("wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz")
  system(paste("mv gencode.v19.annotation.gtf.gz", ann_dir))}
gtf <- rtracklayer::import(con = "gencode.v19.annotation.gtf.gz")
canonical_chr <- paste0("chr",c(1:22,"X","Y","M"))
## gene2name
gene2name <- unique(data.frame(gene_id = remove_dots(gtf$gene_id), gene_name=gtf$gene_name, gene_status=gtf$gene_status, gene_type=gtf$gene_type, stringsAsFactors = F))
tx2name <- with(subset(gtf, type=="transcript"), unique(data.frame(transcript_id, gene_name, stringsAsFactors = F)))

save(gtf,canonical_chr,gene2name,tx2name,file=file.path(ann_dir,"gtf.RData"))

########## generate transcript regions ###################

if(!file.exists(file.path(ann_dir,"mytxdb.RData"))){
  mytxdb <- makeTxDbFromGFF(file=file.path(ann_dir,"gencode.v19.annotation.gtf.gz"))
  saveDb(mytxdb,file = file.path(ann_dir,"mytxdb.RData"))
}
mytxdb=loadDb(file.path(ann_dir,"mytxdb.RData"))

## create annotated , reduced transcript regions for the annotation of peaks
message("Generating transcript regions...")
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

## reduce 3'UTRs to avoid counting same sequence twice
red_3utrs=GenomicRanges::reduce(threutrs,ignore.strand=FALSE, with.revmap=T)
canonical_chr <- paste0("chr",c(1:22,"X","Y","M"))
red_3utrs <- red_3utrs%>%subset(., .@seqnames%in%canonical_chr)
## just name by the first one but keep unique names
red_3utrs$name <- names(threutrs)[unlist(lapply(red_3utrs$revmap, function(x) return(x[[1]])))] %>% make.names(.,unique = T)
names(red_3utrs) <- red_3utrs$name


save(myregions, threeutrs, red_3utrs, file=file.path(ann_dir,"myregions.RData"))


## create features list which is needed as input to the RCAS function
if(!file.exists(file.path(ann_dir,"txdbFeatures.RData"))){
  txdbFeatures <- RCAS::getTxdbFeaturesFromGRanges(gtf)
  save(txdbFeatures, file=file.path(ann_dir,"txdbFeatures.RData"))
}
