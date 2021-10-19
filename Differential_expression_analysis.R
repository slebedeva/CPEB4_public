### DESeq analysis for CPEB4 CLIP

## We want to know if there is any differential CPEB4 binding sites in HDACi treatment
## fold change is RMD over DMSO (Ctrl)
## we use DEseq with interaction term to see if there was any changes in CLIP not explained by underlying changes in RNA-seq

############# load packages ################

basedir <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(basedir)
.libPaths("CLIP_Rlibs") ## local path depends on mountpoint inside container!

source("CLIP.Rprofile")

mypackages <- c("tidyverse"
                ,"GenomicRanges"
                ,"GenomicAlignments"
                ,"reshape2"
                ,"magrittr"
                ,"DESeq2")

suppressPackageStartupMessages(lapply(mypackages, require, character.only=T))


################################## DESeq analysis #####################################

## import counts for RNA-seq
fabian_rnaseq <- list.files("data/RNAseq/", pattern = "_ReadsPerGene.out.tab", full.names = T, recursive = F)
rnaseq <- lapply(fabian_rnaseq, function(x) {
  y <- read.table(x,skip=4, header=F)
  rownames(y) <- y[,1]
  y <- y[,4, drop=F]
  return(y)})
rnaseqcnt <- do.call(cbind, rnaseq)
colnames(rnaseqcnt) <- paste0("rna",c(paste0(rep("C",3),1:3), paste0(rep("R",3), 1:3)))

## summarize counts per gene name 
rnaseqcnt$gene_name <- gene2name$gene_name[match(sub("\\..+","",rownames(rnaseqcnt)),gene2name$gene_id)]
rnacnt <- rnaseqcnt%>%group_by(gene_name)%>%summarize_all(sum)%>%column_to_rownames("gene_name")

## DESeq for RNA
dds <- DESeqDataSetFromMatrix(countData = rnacnt,
                              colData = data.frame(row.names=colnames(rnacnt),
                                                   replicate=rep(1:3,2),
                                                   condition=c(rep("C",3), rep("R",3))),
                              design=~condition)
dds$condition <- relevel(dds$condition, ref = "C")
dds <- DESeq(dds)

### count reads in CLIP target genes

allbams <- list.files("data/bams", pattern = ".bam$", full.names = T)
clipbams <- BamFileList(allbams)

## only include counts that overlap merged set of Ctrl and RMD target sites
## make reduced bed for all CLIP sites in both conditions
load("data/bed_reduced.RData")
all_red <- GenomicRanges::reduce(c(rmd2,ctrl2), drop.empty.ranges = T, with.revmap=T) 

## add biological gene names
## first collapse gene names which are the same
myid <- lengths(all_red$revmap)==1
all_red$gene_names <- NA
all_red$gene_names[myid] <- c(rmd2,ctrl2)$gene_names[unlist(all_red$revmap[myid])]
## collapse genes which are different - although I will not count them (we lose only 1k out of 45k sites)
aux_genen <- c(rmd2,ctrl2)$gene_names ## much faster than nesting GRanges of course :-P also use system.time()
all_red$gene_names[!myid] <- unlist(lapply(all_red$revmap[!myid], function(x){ return(
  paste(aux_genen[x],collapse=";")%>%stringr::str_split(.,";")%>%unlist()%>%sort%>%unique%>%paste(.,collapse = ";")
) })) ## 

## add universal unique name for peaks
all_red$name <- make.names(all_red$gene_names, unique = T)

clipcnt <- GenomicAlignments::summarizeOverlaps(features = all_red,
                             reads = clipbams,
                             mode = "Union",
                             ignore.strand = F)

## aggregate clip counts per gene 
clcnt <- assay(clipcnt)%>%data.frame()
clcnt$gene_name <- all_red$gene_names
clcnt$region <- all_red$regions
## ambiguous gene assigments will be excluded with this index
id <- grepl(";",clcnt$gene_name)

## DE analysis
ccnt <- clcnt[!id,]%>%group_by(gene_name)%>%summarize_all(sum)%>%column_to_rownames("gene_name")
ddscl <- DESeqDataSetFromMatrix(countData = ccnt,
                                colData = data.frame(
                                  row.names=colnames(ccnt),
                                  condition=ifelse(grepl("RMD",colnames(ccnt)),"R","C")),
                                  design=~condition)
ddscl$condition <- relevel(ddscl$condition, ref = "C")
ddscl <- DESeq(ddscl)

### Introduce interaction term ### 
## replace row names with gene name
cnt <- merge(counts(ddscl),counts(dds),by=0,all=F)%>%column_to_rownames("Row.names") #clip counts on all merged reduced sites
colnames(cnt)[1:6]<-paste0(ifelse(grepl("RMD",colnames(ccnt)),"clipR","clipC"),1:6)
ddsIA <- DESeqDataSetFromMatrix(
  countData = cnt,
  colData=data.frame(row.names = colnames(cnt),
                     group=c(rep("clip",6),rep("rna",6)),
                     condition=ifelse(grepl("C",colnames(cnt)),"Ctrl","RMD")),
  design = ~ group + condition + group:condition) 
## choose reference levels: Ctrl and RNAseq
ddsIA$condition <- relevel(ddsIA$condition, ref = "Ctrl")
ddsIA$group <- relevel(ddsIA$group, ref = "rna")
ddsIA <- DESeq(ddsIA, test="LRT", reduced=~group + condition) ## these should be DE genes which have difference logFC between CLIP and rnaseq
## shrinking logFC (type "normal" does not work for desing with interaction term) 
resIAshrunk <- as.data.frame(lfcShrink(ddsIA, "groupclip.conditionRMD", type='apeglm')[order(results(ddsIA)$padj),])%>%rownames_to_column("gene_name")

## save data
save(all_red, dds, ddscl, ddsIA, file="data/DEseq.RData")

