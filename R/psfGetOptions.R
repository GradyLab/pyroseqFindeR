psfGetOptions <- function(){
  # dependancies
  suppressMessages(require('BSgenome.Hsapiens.UCSC.hg19'))
  suppressMessages(require("IlluminaHumanMethylation450kprobe"))
  
  x <- readline(message("Do you want to include annotations with your query? (enter y,n or default)"))
  if(x=="y")
  {
    data(IlluminaHumanMethylation450kprobe)
    anno <- IlluminaHumanMethylation450kprobe

    message(paste0("Select what HM450 Manifest annotations to include. Column options are:"))
    for(i in 1:length(colnames(anno)))
    {
      message(paste0(i,". ",colnames(anno)[i]))
    }
    
    anno.cols <- readline(message("Enter column numbers to include (separated by ;): "))
    anno.cols <- as.numeric(gsub("[^0-9]","",unique(unlist(strsplit(anno.cols,";")))))
    psfGetFlankseq(anno.cols=anno.cols)
  } else{
    if(x=="default"){
      psfGetFlankseq(anno.cols=c(1,2,3,7,8,18,19,20,24,26))
    } else{
      psfGetFlankseq()
    }
  }
}
