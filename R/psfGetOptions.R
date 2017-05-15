psfGetOptions <- function(){
  # dependancies
  message("#################################\n### Welcome to PyroseqFindeR! ###\n#################################")
  message("loading human genome, bld hg19..")
  suppressMessages(require('BSgenome.Hsapiens.UCSC.hg19'))
  message("loading HM450 annotation..")
  suppressMessages(require("IlluminaHumanMethylation450kanno.ilmn12.hg19"))
  suppressMessages(require("minfi"))
  
  x <- readline(message("Do you want to include annotations with your query? (enter y,n or default)"))
  if(x=="y")
  {
    message("Bringing up annotation options..")
    data(IlluminaHumanMethylation450kanno.ilmn12.hg19)
    usranno <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
 
    message(paste0("Select what HM450 Manifest annotations to include. Column options are:"))
    for(i in 1:length(colnames(usranno)))
    {
      message(paste0(i,". ",colnames(usranno)[i]))
    }
    
    anno.cols <- readline(message("Enter column numbers to include (separated by ;): "))
    anno.cols <- as.numeric(gsub("[^0-9]","",unique(unlist(strsplit(anno.cols,";")))))
    psfGetFlankseq(anno=usranno,anno.cols)
  } else{
    if(x=="default"){
      psfGetFlankseq(anno=usranno,anno.cols=c(1,2,3,7,8,18,19,20,24,26))
    } else{
      psfGetFlankseq(anno=NULL,anno.cols=NULL)
    }
  }
}
