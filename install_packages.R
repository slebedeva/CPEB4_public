mypackages <- c("devtools"
		,"config"
		,"plyranges"
                ,"tidyverse"
                ,"reshape2"
                ,"Hmisc"
                ,"magrittr"
                ,"ggplot2"
                ,"ggrepel"
                ,"ggpubr"
                ,"lattice"
                ,"Biostrings"
                ,"GenomicRanges"
                ,"GenomicFeatures"
                ,"rtracklayer"
                ,"Rsamtools"
                ,"motifStack"
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
                ,"ggunchained"
                ,"DESeq2"
                ,"Gviz"
                )

#### try to install missing packages if any (can easily fail) #####


to_install=mypackages[!mypackages%in%installed.packages()]

if(length(to_install)>0){
  message("missing packages: ", paste(to_install,collapse=", "), "\ntrying to install...")
  BiocManager::install(pkgs = to_install,update = FALSE, ask = FALSE)
  failed=mypackages[!mypackages%in%installed.packages()]
  ## ggunchained separately (not in Bioc) github
  if("ggunchained"%in%failed){devtools::install_github("JanCoUnchained/ggunchained")}
}

failed=mypackages[!mypackages%in%installed.packages()]
if(length(failed)>0){message("failed to install: ", paste(failed,collapse=", "), "\nanyway proceeding...")}

