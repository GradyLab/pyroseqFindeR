# purpose: functions for pyroseqFinder R package
#

# dependancies
require('BSgenome.Hsapiens.UCSC.hg19')
require("pyroseqFinder") # load HM450 annotation df 'anno', same as manifest from minfi package


psf.getOptions <- function(){
  x <- readline(message("Do you want to include annotations with your query? (enter y,n or default)"))
  if(x=="y"){
    message(paste0("Select what HM450 Manifest annotations to include. Column options are:"))
    for(i in 1:length(colnames(anno))){message(paste0(i,". ",colnames(anno)[i]))}
    anno.cols <- readline(message("Enter column numbers to include (separated by ;): "))
    anno.cols <- as.numeric(gsub("[^0-9]","",unique(unlist(strsplit(anno.cols,";")))))
    getFlankseq(anno.cols=anno.cols)
  } else{
    if(x=="default"){
      getFlankseq(anno.cols=c(1,2,3,7,8,18,19,20,24,26))
    } else{
      psf.getFlankseq()
    }
  }
}


psf.getFlankseq <- function(anno.cols=NULL){
  suppressMessages(suppressWarnings(require('BSgenome.Hsapiens.UCSC.hg19')))
  
  x <- readline(message("Provide one or more probe IDs to search (separated by ';'):"))
  probes <- unique(unlist(strsplit(x,";")))
  if(nrow(anno[rownames(anno) %in% probes,])<length(probes)){
    message("Not all of the entered probe IDs are valid!")
    message("...done!")
  } else{
    
    y <- as.numeric(readline(message("How many bp up/downstream do you want to look? ")))
    if(!is.numeric(y)){
      message("Invalid basepair range entered!")
      message("...done!")
    } else{
      
      return.matrix <- as.data.frame(matrix(nrow=length(probes),ncol=3))
      colnames(return.matrix) <- c("probename",paste0(y,"bp_upstream"),paste0(y,"bp_downstream"))
      
      for(i in 1:length(probes)){
        return.matrix[i,1] <- probes[i]
        return.matrix[i,2] <- as.character(getSeq(Hsapiens,anno[rownames(anno)==return.matrix[i,1],]$chr,
                                                  anno[rownames(anno)==return.matrix[i,1],]$pos-(y-1),
                                                  anno[rownames(anno)==return.matrix[i,1],]$pos))
        return.matrix[i,3] <- as.character(getSeq(Hsapiens,anno[rownames(anno)==return.matrix[i,1],]$chr,
                                                  anno[rownames(anno)==return.matrix[i,1],]$pos,
                                                  anno[rownames(anno)==return.matrix[i,1],]$pos+(y-1)))
      }
      if(!is.null(anno.cols)){
        anno.select <- anno[rownames(anno) %in% probes,]
        anno.select <- anno.select[order(match(anno.select$Name,probes)),]
        for(a in 1:length(anno.cols)){
          which.anno <- as.numeric(anno.cols[a])
          return.matrix <- cbind(return.matrix,
                                 anno.select[,which.anno])
          colnames(return.matrix)[ncol(return.matrix)] <- colnames(anno)[which.anno]
        }
      }
      
      message("Here is a preview of your probe info:")
      print(head(return.matrix))
      
      z <- readline(message("Save this matrix as a csv? (enter y/n) "))
      if(z=="y"){
        filename <- readline(message("Enter your filename:"))
        write.csv(return.matrix,file=paste0(filename,".csv"))
        return(message(paste0("Probe info saved to: '",getwd(),"/",filename,"'")))
        
      } else{
        message("Here is your final probe info:")
        return(return.matrix)
        
      }
    }
  }
}
