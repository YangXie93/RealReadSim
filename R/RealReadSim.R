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

realReadSim <- function(filenames_csv,coverage,bowtieOptions = "--no-unal",humanReadable = FALSE,readAsBams = TRUE, minMapq = 40, redraw = FALSE, repeatable = FALSE,seed = 0,minContigLength = 10000){

    starttime = Sys.time()
    filenames = read.table(filenames_csv,header = TRUE,sep = ",",stringsAsFactors = FALSE)

    covAt = 1
    temp = list()

    rangs = list()
    mCov = list()
    nms = c()

    cov = list()
    params = ScanBamParam(what = c("pos","qwidth","rname"),mapqFilter = minMapq)  #(?) mapqFilter
    seqNames = c()
    tempNames = list()

    for(i in 1:length(filenames[,1])){          # reading in the names of all sequences involved
        if(readAsBams){                         # reading from bam files
            tempNames[[i]] = scanBamHeader(filenames[i,2])[[1]]$targets
            seqNames[i] = names(tempNames[[i]][which.max(tempNames[[i]])])
            tempNames[[i]] = names(tempNames[[i]])
        }
        else{                                   # reading from fasta files
            seqData = readDNAStringSet(filenames[i,1])
            tempNames[[i]] = gsub(" .*","",names(seqData))
            seqNames[i] = tempNames[[i]][which.max(widht(seqData))]
        }
    }

    #--------------------------- checking if already in Datasystem and if not adding to datasystem -------------

    if(readAsBams){
        catalogue = RealReadSim::addToDataSystem(seqNames,bam = filenames[,2],fasta = filenames[,1],minMapq = minMapq,bowtieOptions)
    }
    else{
        if(length(filenames[,2][grep(" ",filenames[,2])]) > 0){
            catalogue = RealReadSim::addToDataSystem(seqNames,fasta = filenames[,1],fastq1 = gsub(".*","",filenames[,2]),fastq2 = gsub(".* ","",filenames[,2]),minMapq = minMapq,bowtieOptions)
         }
        else{
            catalogue = RealReadSim::addToDataSystem(seqNames,fasta = filenames[,1],fastq1 = filenames[,2],minMapq = minMapq,bowtieOptions)
        }
    }

    #------------------------------------------------------------------------------------------------------------

    print(tempNames)############################################################################################

    #-------------------------------- assembly --------------------------------------------------------
    RRSDS = "~/RealReadSimDS"
    print(length(seqNames))############################################################################
    for(i in 1:length(tempNames)){

        print("readiying the data")##################################################################

        #--------------------------- readiyng the data -----------------------------------------------------
        seqData = readDNAStringSet(as.character(catalogue$fastaName[i]))

        print(length(seqData))#################################################################
        print("get and save contigs")###########################################################

        #-------------------------- get and save the contigs -----------------------------------------------
        for(j in 1:length(tempNames[[i]])){
            print(paste(RRSDS,"/",seqNames[i],"/",tempNames[[i]][j],".Rds",sep =""))

            ttt = Sys.time()###################################################################

            data = readRDS(paste(RRSDS,"/",seqNames[i],"/",tempNames[[i]][j],".Rds",sep =""))

            print(Sys.time() -ttt)#########################################################################
            ttt = Sys.time()########################################

            data = data.table(pos = start(data),width = width(data),DNAString = names(data))

            print(Sys.time() -ttt)#################################################

            data = RealReadSim::randomReads(data,catalogue$totalLength[i],coverage = coverage[covAt],meanWidht = mean(data$width),repeatable,seed,redraw)

            print(substr(names(seqData),1,nchar(tempNames[[i]][j])))
            print(tempNames[[i]][j])

            partialSeqData = seqData[substr(names(seqData),1,nchar(tempNames[[i]][j])) == tempNames[[i]][j]]

            if(length(data$pos) > 0){

                print(length(partialSeqData))###############################################################

                cov[[length(cov) +1]] = RealReadSim::getCoverage(data$pos,data$width,partialSeqData@ranges@width[1])
                contigs = RealReadSim::evalCoverage(data$pos, data$width, partialSeqData@ranges@width[1],toString(partialSeqData[[1]]),minContigLength)
                if(length(contigs) > 0){
                    rangs[[length(rangs) +1]] = IRanges(start = contigs[seq(1,length(contigs),2)],end = contigs[seq(2,length(contigs),2)])
                    mCov[[length(mCov) +1]] = RealReadSim::meanCovToRange(contigs,cov[[length(cov)]])
                    nms[[length(nms) +1]] = tempNames[[i]][j]
                }
            }
        }
        if(covAt < length(coverage)){
            covAt = covAt +1
        }
        print(Sys.time() -starttime)#################################################################
    }
  temp = list(rangs,mCov,nms)
  print("co-assembly")################################################################################

  #----------------------------------- co-assembly ----------------------------------------------------

  #res = coAssembleRRSDS(temp)

  if(humanReadable){
    #RealReadSim::makeHumanReadable(cov,res)
  }
  return(temp)
}


makeHumanReadable <- function(cov,res){
  par(mfrow = c(2,1))
  x = seq(1,length(res),by = 3)
  for(i in 1:length(cov)){
    cov[[i]] = RealReadSim::slidingWindowMaker(cov[[i]])
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
  crossmaps = dir(paste(RRSDS,"Crossmaps",sep = ""))

  seqNames = contigs[[3]]

  partners = as.data.frame(strsplit(crossmaps,"_X_|\\.Rds"))

  pairs = list()
  whichFile = list()
  for(i in 1:length(seqNames)){

    col1 = which(partners[1,] == seqNames[i])
    col2 = which(partners[2,] == seqNames[i])

    other1 = partners[2,col1]
    other2 = partners[1,col2]

    partners = partners[-c(col1,col2)]

    other = c(other1,other2)
    whichFile[[i]] = c(col1,col2)
    pairs[[i]] = which(seqNames %in% other)
  }

  res = list()
  conts = list()
  mCov = list()

  for(i in 1:length(pairs)){
    if(length(pairs[[i]]) > 0){
      for(j in 1:length(pairs[[i]])){
        file = crossmaps[whichFile[[i]][j]]
        crossmap = readRDS(paste(RRSDS,"Crossmaps/",file,sep = ""))

        reads1 = crossmap[ findOverlaps(crossmap[[seqNames[i]]],contigs[[1]][[i]])]
        reads2 = crossmap[findOverlaps(crossmap[[seqNames[pairs[[i]][j]]]],contigs[[1]][[pairs[[i]][j]]])]

        overlaps1 = reduce(overlapsRanges(contigs[[1]][[i]],reads1))
        overlaps2 = reduce(overlapsRanges(contigs[[1]][[pairs[[i]][j]]],reads2))

        if(length(reads1) > 0){


          matching1 = overlapsRanges(contigs[[1]][[i]],reads1)
          matching2 = overlapsRanges(contigs[[1]][[pairs[[i]][j]]],reads2)



        }
        else{
          conts = contigs[[1]][[i]]
          mCov = contigs[[2]][[i]]
        }
      }
    }
    else{
      conts = contigs[[1]][[i]]
      mCov = contigs[[2]][[i]]
    }
  }
}
