
binningTestData <- function(filenames_csv,coverage,prebuild = TRUE,bowtieDefaultOptions = TRUE,humanReadable = FALSE,asTempoaryDatabase = TRUE,readAsBams = TRUE, minMapq = 40, redraw = NA, repeatable = FALSE,seed = 0){

  starttime = Sys.time()
  filenames = read.table(filenames_csv,header = TRUE,sep = ",",stringsAsFactors = FALSE)
  temp = list()
  cov = list()
  params = ScanBamParam(what = c("pos","qwidth","rname"),mapqFilter = minMapq)  #(?)
  covAt = 1
  
  sequenceData = seqsAsTempoaryData(filenames[,1])
  print(typeof(sequenceData$seq))
  if(readAsBams && asTempoaryDatabase){
    data = bamsAsTempoaryData(filenames[,2],params,coverage,sequenceData$lengths,redraw)
    times = 1
  }
  else{
    times = length(filenames[,1])
  }
  
  for(i in 1:times){
    if(!asTempoaryDatabase){
      data = bamOrBowtie(filenames,readAsBams,i,params,bowtieDefaultOptions)
      data = randomReads(data,sequenceData$lengths[i],coverage[covAt],mean(data$width),repeatable,seed,redraw)
    }
    
    print(Sys.time() -starttime)
    
    for(j in 1:length(sequenceData$seqs)){
      partialData = subset(data, DNAString == unique(data$DNAString)[j])
      name = toString(unique(partialData$DNAString))
      partialSeqData = subset(sequenceData,substr(DNAString,1,nchar(name)) == toString(name))
      print(length(partialData$pos))
      print(unique(partialData$DNAString))
      if(length(partialData$pos) > 1){
        cov[[i]] = getCoverage(partialData$pos,partialData$width,partialSeqData$lengths[1])
        contigs = evalCoverage(partialData$pos, partialData$width, partialSeqData$lengths[1],partialSeqData$seqs[1])
        temp[[length(temp) +1]] = contigs
        temp[[length(temp) +1]] = meanCovToRange(contigs,cov[[i]])
        temp[[length(temp) +1]] = unique(data$DNAString)[j]
      }
      
    }
    
    if(covAt < length(coverage)){
      covAt = covAt +1
    }
    
    print(Sys.time() -starttime)
  }
  
  res = formatResults(temp)
  
  if(humanReadable){
    makeHumanReadable(cov,res)
  }
  return(res)
}


makeHumanReadable <- function(cov,res){
  par(mfrow = c(2,1))
  x = seq(1,length(res),by = 3)
  for(i in 1:length(cov)){
    cov[[i]] = slidingWindowMaker(cov[[i]])
    pdf(paste(names(res[[x[i]]][1]),".pdf",sep = ""))
    plot.default(x = cov[[i]],y = c(1:length(cov[[i]])),type = "l",xlab = "Coverage",ylab = "Position")
    hist(width(res[[x[i]]]),border = "black",col = "red",breaks = 50,xlab = "Contig Länge",main = "Contig Längenhistogramm")
    abline(v = res[[x[i] +2]],col = "blue")
    dev.off()
  }
  par(mfrow = c(1,1))
}


formatResults<-function(temp){
  res = list()
  i = 1
  for(n in seq(1,length(temp),by = 3) ){
    if( !(length(temp[[n]]) == 0)){
      res[[n]] = IRanges(start = temp[[n]][seq(1,length(temp[[n]]) -1,by = 2)],end = temp[[n]][seq(2,length(temp[[n]]),by = 2)],names = rep(temp[[n +2]],times = (length(temp[[n]])/2)) )
      res[[n +1]] = temp[[n+1]]
      res[[n +2]] = N50(width(res[[n]]))
    }
  }
  return(res)
}


bamsAsTempoaryData <- function(filenames,params,coverage,seqLenghts,redraw){
  temp = list()
  for(i in 1:length(filenames)){
    temp[[i]] = scanBam(filenames[i],param = params)
  }
  pos = c()
  width = c()
  DNAString = c()
  
  j = 1
  covAt = 1
  for(i in 1:length(temp)){
    current = temp[[i]][[1]]
    chromeNr = length(levels(current$rnames))
    if(chromeNr > 1){
      seqLeng = sum(seqLengths[c(j:(j+chromeNr -1))])
      j = j + chromeNr -1
    }
    else{
      seqLeng = seqLenghts[j]
      j = j+1
    }
    notUsable = which(is.na(current$pos)) *-1
    if(length(notUsable) > 0){
      current$pos = current$pos[notUsable]
      current$qwidth = current$qwidth[notUsable]
      current$rname = current$rname[notUsable]
    }

    
    sampleSize = as.integer((seqLeng *coverage[covAt])/mean(current$qwidth))
    set.seed(sampleSize)
    elected = sample(c(1:length(current$pos)),sampleSize)
    
    pos = append(pos,current$pos[elected])
    width = append(width,current$qwidth[elected])
    DNAString = append(DNAString,as.character(current$rname[elected]))
    
    if(covAt < length(coverage)){
      covAt = covAt +1
    }
  }
  
  return(data.frame(pos = pos,width = width,DNAString = DNAString))
}

seqsAsTempoaryData <- function(filenames){
  temp = list()
  for(i in 1:length(filenames)){
    temp[[i]] = sequenceAndLength(filenames[i]) 
  }
  seqs = c()
  lengths = c()
  DNAString = c()
  
  for(i in 1:length(temp)){
    seqs = append(seqs,temp[[i]]$seq)
    lengths = append(lengths,temp[[i]]$length)
    DNAString = append(DNAString,temp[[i]]$DNAString)
  }
  
  return(data.frame(seqs = seqs,lengths = lengths,DNAString = DNAString,stringsAsFactors = FALSE))
}