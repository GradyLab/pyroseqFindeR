psfGetFlankseq <- function(anno,anno.cols=NULL){
  
  x <- readline(message("Provide one or more probe IDs to search, separated by a semicolon(s) (eg. 'cg00035673;cg00040027;..'):"))
  probes <- unique(unlist(strsplit(x,";")))
  
  if(nrow(anno[rownames(anno) %in% probes,])<length(probes)){
    message("Error: Not all of the entered CpG probe IDs are valid.")
    message("Returning..")
  } else{

    y <- as.numeric(readline(message("How many bp up/downstream do you want to look? (reading 5' => 3' on pos/+ strand relative to CpG position + 1 nt; consult chromosome coordinates if unsure of results)")))
    
    if(!is.numeric(y)){
      message("Invalid basepair range entered!")
      message("...done!")
    } else{

      return.matrix <- as.data.frame(matrix(nrow=length(probes),ncol=3))
      colnames(return.matrix) <- c("probename",paste0(y,"nt_upstream"),paste0(y,"nt_downstream"))

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
