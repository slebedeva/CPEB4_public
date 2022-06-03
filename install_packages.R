## propmt user to choose directory because config may be not there yet
libdir=tcltk::tk_choose.dir(caption="Choose directory for R libraries (double click)" )
.libPaths(libdir)

tcltk::tkmessageBox(message=paste("Do you want to proceed installing into ",.libPaths()[1],"?"),type="okcancel")

message('your librariy path: ',.libPaths()[1])

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
		,"apeglm"
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
if(length(failed)>0){message("Failed to install: ", paste(failed,collapse=", "), "\nAnyway proceeding...")}

if(length(failed)==0){message('All packages successfully installed!')}

config=config::get()

if(is.null(config$Rlibdir)){
	message('Adding your library path to config...')
	command=paste0('echo ','\'', '  Rlibdir: ', libdir, '\'', ' >> config.yml')
	system(command)
	# if(!file.exists(".Rprofile")){
	#   file.create(".Rprofile")
	#   writeLines(text=paste0(".libPaths(\'",libdir,"\')"), con=".Rprofile")
	#   }
}
