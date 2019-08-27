IRangesTranslator <- function(seq1,seq2){
  y = Sys.time()
  seq1_2 = reduce(seq1)
  seq2_2 = reduce(seq2)

  names(seq1_2) <- c(1:length(seq1_2))
  names(seq2_2) <- c(1:length(seq2_2))
  con1 = list()
  con2 = list()

  for(i in 1:length(seq1_2)){
    tmp = subsetByOverlaps(seq1_2[i],seq1)
    tmp1 = subsetByOverlaps(seq1,tmp)
    tmp2 = seq2[which(names(seq2) %in% names(tmp1))]
    tmp2 = subsetByOverlaps(seq2_2,tmp2)
    con1[[i]] = which(names(seq2_2) %in% names(tmp2))
  }

  x = TRUE

  for(i in 1:length(names(seq2_2))){
    for(j in 1:length(con1)){
      if(i %in% con1[[j]]){
        if(!x){
          con2[[i]] = append(con2[[i]],j)
        }
        else{
          con2[[i]] = j
        }
        x = FALSE
      }
    }
    if(x ){
      con2[[i]] = 0
    }
    x = TRUE
  }

  res = list(seq1_2,con1,seq2_2,con2)
  print(Sys.time() -y)
  return(res)
}
