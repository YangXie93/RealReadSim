getLotsOfN50s <- function(filename,coverages){
  x = list()
  res = c()
  for(i in 1:length(coverages)){
    print(paste("Durchlauf:",i,sep = " "))
    x = binningTestData(filename,coverages[i],readAsBams = TRUE)
    for(j in seq(3,length(x),by = 3)){
      print(x[[j]])
      res[length(res) +1] = x[[j]]
    }
  }
  return(res)
}

getMoreN50s <- function(filename,x){
  res = data.frame(Nr = c(1:20))
  for(i in 1:x){
    res[paste("N50s",i)] = getLotsOfN50s(filename,c(1:20))
  }
  return(res)
}