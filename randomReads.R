randomReads <- function(data,seqLength,coverage,meanWidht,repeatable,seed,redraw){
  
  numberOfReads = as.integer((sum(seqLength) *coverage)/meanWidht)
  
  if(repeatable){
    set.seed(seed)
  }
  
  if(numberOfReads > length(data$pos)){
    if(is.na(redraw)){
      print("the reads are not enough to produce the requested coverage would you like to draw again? [y|n]")
      answer = readline()
      redraw = FALSE
    }
    else{
      answer = "z"
    }
    if(answer == "y" || answer == "yes" || redraw){
      help = data.frame()
      while(numberOfReads > length(data$pos)){
        help = rbind(help,data)
        numberOfReads = numberOfReads - length(data$pos)
      }
      data = rbind(help,data[sample(c(1:length(data$pos)),numberOfReads),])
      rm(help)
    }
    else{
      numberOfReads = length(data$pos)
      newCov = as.integer(numberOfReads *meanWidht / sum(seqLength))
      print(paste("the new coverage value is:",newCov,sep = " "))
    }
  }
  else{
    print(paste((length(data$pos) -numberOfReads)," reads were thrown out"))
    data = data[sample(c(1:length(data$pos)),numberOfReads),]
  }
  return(data)
}