sequenceAndLength <- function(fastaName){
  
  fastaSequences = readDNAStringSet(fastaName)
  seq = c()
  seqLength = c()
  DNAString = c()
  
  for(j in 1:length(fastaSequences)){
    
    seq[j] = toString(fastaSequences[j])
    seqLength[j] = width(fastaSequences[j])
    DNAString[j] = names(fastaSequences[j])
  }
  res = data.frame(seq = seq,length = seqLength,DNAString = DNAString,stringsAsFactors = FALSE)
  return(res)
}