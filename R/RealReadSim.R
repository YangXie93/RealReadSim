#' realReadSim
#' @description This function simulates metagenomic data using real experimental reads
#'
#' @param filenames_csv This parameter takes the name of a file in which the names of the fasta files and the name of the bam or fastq files and or the bowtie prebuild name are contained
#' @param coverage A integer vector containing the wished coverages of the individual genomes
#' @param prebuild A boolean value that determines wether the prebuild function is to be used
#' @param bowtieDefaultOptions A boolean value determening wether the set bowtie2 options schould be used or user defined options
#' @param humanReadable A boolean value determening wether there will be graphical output or not
#' @param asTempoaryDatabase A boolean value determening wether the data will be read in all at once or one by one
#' @param readAsBams A boolean value determening wether filenames_csv should be read as containing bam file names or not
#' @param minMapq A integer value stating the minimal mapq score of reads to use
#' @param redraw A boolean value determening, if the asked for coverage is to much for the data set, reads from the same data set should be drawn again
#' @param repeatable A boolean value determening wether a seed will be set
#' @param seed A integer value used as seed when reapeatable = TRUE
#' @export

realReadSim <- function(filenames_csv,coverage,bowtieOptions = "--no-unal",humanReadable = FALSE,readAsBams = TRUE, minMapq = 40, redraw = FALSE, repeatable = FALSE,seed = 0){

  starttime = Sys.time()
  filenames = read.table(filenames_csv,header = TRUE,sep = ",",stringsAsFactors = FALSE)

  covAt = 1
  temp = list()
  cov = list()
  params = ScanBamParam(what = c("pos","qwidth","rname"),mapqFilter = minMapq)  #(?)
  seqNames = c()

  for(i in 1:length(filenames[,1])){
    if(readAsBams){
      print(filenames)
      tempNames = scanBamHeader(filenames[i,2])[[1]]$targets
      seqNames[i] = names(tempNames[which.max(tempNames)])
    }
    else{
      tempNames = readDNAStringSet(filenames[i,1])
      seqNames[i] = names(tempNames)[which.max(widht(tempNames))]
    }
  }


  if(readAsBams){
    catalogue = addToDataSystem(seqNames,bam = filenames[,2],fasta = filenames[,1],minMapq = minMapq,bowtieOptions)
  }
  else{
    if(length(filenames[,2][grep(" ",filenames[,2])]) > 0){
      catalogue = addToDataSystem(seqNames,fasta = filenames[,1],fastq1 = gsub(" .*","",filenames[,2]),fastq1 = gsub(".* ","",filenames[,2]),minMapq = minMapq,bowtieOptions)
     }
    else{
      catalogue = addToDataSystem(seqNames,fasta = filenames[,1],fastq1 = filenames[,2],minMapq = minMapq,bowtieOptions)
    }
  }

  print("assembly")

  #-------------------------------- assembly --------------------------------------------------------
  RRSDS = "~/RealReadSimDS"
  print(length(seqNames))
  for(i in 1:length(seqNames)){

    print("readiying the data")

    #--------------------------- readiyng the data -----------------------------------------------------
    data = readRDS(paste(RRSDS,"/",seqNames[i],"/",seqNames[i],".Rds",sep =""))
    data = data.table(pos = start(data),width = width(data),DNAString = names(data))
    data = randomReads(data,catalogue$totalLength[i],coverage = coverage[covAt],meanWidht = mean(data$width),repeatable,seed,redraw)
    seqData = readDNAStringSet(as.character(catalogue$fastaName[i]))

    print(length(seqData))
    print("get and save contigs")

    #-------------------------- get and save the contigs -----------------------------------------------
    partialSeqData = seqData[substr(names(seqData),1,nchar(seqNames[i])) == seqNames[i]]

    print(1)

    if(length(data$pos) > 0){

      print(length(partialSeqData))

      cov[[i]] = RealReadSim::getCoverage(data$pos,data$width,partialSeqData@ranges@width[1])
      contigs = RealReadSim::evalCoverage(data$pos, data$width, partialSeqData@ranges@width[1],toString(partialSeqData[[1]]))
      temptemp = RealReadSim::meanCovToRange(contigs,cov[[i]])
      temp[[length(temp) +1]] = IRanges(start = contigs,end = contigs,names = seqNames[i])
      temp[[length(temp)]]@metadata = temptemp
    }

    if(covAt < length(coverage)){
      covAt = covAt +1
    }
    print(Sys.time() -starttime)
  }

  print("co-assembly")

  #----------------------------------- co-assembly ----------------------------------------------------

  #res = coAssembleRRSDS(temp)

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


coAssembleRRSDS <- function(contigs){
  RRSDS = "~/RealReadSimDS/"
  seqNames = c()
  n = 0
  switch = TRUE
  sets = list()
  which = list()
  for(i in seq(3,length(contigs),step = 2)){
    if(names(contigs[[i]])[1] != names(contigs[[i-2]])[1]){
      switch = TRUE
      sets = append(sets,contigs[[i]])
      which = append(which,i)
    }
    if(switch){
      n = n +1
      seqNames[n] = names(contigs[[i]])
      sets[[n]] = append(sets[[n]],contigs[[i]])
      which[[n]] = append(which[[n]],i)
      switch = FALSE
    }
  }

  filter = paste(seqNames,collapse = TRUE,sep = "/")
  filter = paste("[",filter,"]","_X_","[",filter,"].*",sep = "")
  crossFiles = dir(paste(RRSDS,"Crossmaps",sep = ""))
  crossFiles = crossFiles[grep(filter,crossFiles)]
  res = list()
  switch = TRUE

  for(i in seq(1,length(contigs),2)){
    xCross = crossFiles[grep(names(contigs[[i]])[1],crossFiles)]
    n = which(seqNames == names(contigs[[i]])[1])
    if(length(xCross) > 0){
      filesI = xCross[grep(paste(seqNames[n],"_X_",sep = ""))]
      for(j in 1:n){
        fileIJ = filesI[grep(seqNames[j],filesI)]
        if(length(fileIJ) == 1){
          tmp = readRDS(paste(RRSDS,"/Crossmaps/",fileIJ,sep = ""))

          tmp1 = subsetByOverlaps(tmp[[1]],contigs[[i]])
          toCutStart = subset(tmp1,start < start(contigs[[i]]))
          toCutEnd = subset(tmp1,end > end(contigs[[i]]))
          tmp2 = subset(tmp[[2]],names %in% names(tmp1))
          start(tmp2)[names(tmp2) %in% names(toCutStart)] = start(tmp2)[names(tmp2) %in% names(toCutStart)] + (start(toCutStart) - start(contigs[[i]]))
          end(tmp2)[names(tmp2) %in% nmaes(toCutEnd)] = end(tmp2)[names(tmp2) %in% names(toCutStart)] + (end(contigs[[i]]) - end(toCutEnd))

          tmp3 = subsetByOverlaps(sets[[j]],tmp2)
          for(n in 1:length(tmp3)){
            contigToReads = data.frame(startDiff = start(tmp2)-start(tmp3[n]),endDiff = end(tmp3[n]) -end(tmp2),names = names(tmp2))
            contigToReads = subset(contigToReads,!(endDiff > width(tmp3[n])) || !(startDiff > width(tmp3[n])))
            readsToSeq = subset(tmp1,names %in% contigToReads$names)
            partsOfSeq = reduce(readsToSeq)
            guess = IRanges(start = subset(tmp2,names == contigToReads$name)$start -contigToReads$startDiff, end = subset(tmp2,names == contigToReads$name)$end +contigToReads$endDiff)


          }

          sets[[j]] = subset(sets[[j]],sets[[j]] != tmp3)
          swtich = FALSE
        }
      }
      if(switch){
        res = append(res,contigs[[i]])
        res = append(res,contigs[[i+1]])
      }
      switch = TRUE
    }
  }
  return(res)
}
