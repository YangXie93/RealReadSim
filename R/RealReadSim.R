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

realReadSim <- function(filenames_csv,coverage = 0,takeAll = TRUE,nrOfSamples = 1,bowtieOptions = "--no-unal",humanReadable = FALSE,readAsBams = TRUE, minMapq = 40, redraw = FALSE, repeatable = TRUE,seed = 0,minContigLength = 500,minDist = 0.8){

    library(Rsamtools)
    library(data.table)
    library(seqinr)

    starttime = Sys.time()
    filenames = read.table(filenames_csv,header = TRUE,sep = ",",stringsAsFactors = FALSE)

    covAt = 1

    params = ScanBamParam(what = c("pos","qwidth","rname"),mapqFilter = minMapq)
    seqNames = c()

    for(i in 1:length(filenames[,1])){          # reading in the names of all sequences involved
        if(readAsBams){                         # reading from bam files
            tempNames = scanBamHeader(filenames[i,2])[[1]]$targets
            seqNames[i] = names(tempNames[which.max(tempNames)])
        }
        else{                                   # reading from fasta files
            seqData = readDNAStringSet(filenames[i,1])
            seqNames[i] = gsub(" .*","",names(seqData))[which.max(width(seqData))]
        }
    }

    #--------------------------- checking if already in Datasystem and if not adding to datasystem -------------

    if(readAsBams){
        catalogue = RealReadSim::addToDataSystem(seqNames,bam = filenames[,2],fasta = filenames[,1],minMapq = minMapq,bowtieOptions = bowtieOptions)
    }
    else{
        if(length(filenames) > 2){
            catalogue = RealReadSim::addToDataSystem(seqNames,fasta = filenames[,1],fastq1 = filenames[,2],fastq2 = filenames[,3],minMapq = minMapq,bowtieOptions =  bowtieOptions )
        }
        else{
            catalogue = RealReadSim::addToDataSystem(seqNames,fasta = filenames[,1],fastq1 = filenames[,2],minMapq = minMapq,bowtieOptions =  bowtieOptions)
        }
    }

    #------------------------------------------------------------------------------------------------------------

    #-------------------------------- assembly --------------------------------------------------------
    RRSDS = "~/RealReadSimDS"

    starts = list()
    ends = list()
    readsPerSample = list()
    sdCov = list()
    seqs = list()
    cov = list()
    sequ = list()
    DS = unique(catalogue$dir)

    n = 1
    for(i in 1:length(DS)){

        partialCatalogue = subset(catalogue,dir == DS[i])
        fullData = readRDS(partialCatalogue$data[1])
        fullData = RealReadSim::randomReads(fullData,partialCatalogue$totalLength[1],coverage = coverage[covAt,],meanWidht = partialCatalogue$meanWidth[1],repeatable,seed,redraw,nrOfSamples,takeAll)
        for(j in 1:length(partialCatalogue$name)){

            data = subset(fullData,seq == partialCatalogue$name[j])
            if(length(data$pos) > 0){
                print(paste(length(data$pos)," reads have been selected"))###########################################
                contigs = RealReadSim::evalCoverage(data$pos, data$width,data$sampleNr, partialCatalogue$length[j],partialCatalogue$minOverlap[j],minContigLength,nrOfSamples)
                if(length(contigs[[1]]) > 0){
                    fst = readDNAStringSet(partialCatalogue$fasta[j])
                    starts[[n]] = contigs[[1]]
                    ends[[n]] = contigs[[2]]
                    sequ[[n]] = subSeqs(toString(fst),starts[[n]],ends[[n]])
                    cov = append(cov,contigs[[3]])
                    readsPerSample = append(readsPerSample, contigs[[4]])
                    seqs[[n]] = rep(partialCatalogue$name[j],length(contigs[[1]]))
                    n = n+1
                }
            }

        }
        if(covAt < length(coverage)){
            covAt = covAt +1
        }

    }
    temp = data.table(start = unlist(starts),end = unlist(ends),coverage = cov,seqName = unlist(seqs),seq = unlist(sequ),rps = readsPerSample)
    print("co-assembly")################################################################################

    #----------------------------------- co-assembly ----------------------------------------------------

    res = coAssembleRRSDS(temp,minDist)
    res$length = (res$end-res$start)+1
    res = subset(res,length >= minContigLength)
    res$covVec = calcCovVec(res$rps,res$length)
    #---------------------------------------------------------------------------------------------------
    if(humanReadable){
        #RealReadSim::makeHumanReadable(cov,res)
    }

    print(Sys.time() -starttime)#################################################################
    return(res)
}

# function to simulate coassembly using coAssembly.cpp
#
coAssembleRRSDS <- function(contigs,minDist){
    RRSDS = paste0(Sys.getenv("HOME"),"/RealReadSimDS/")

    involved = unique(contigs$seqName)

    crossmaps = dir(paste0(RRSDS,"Crossmaps"))
    tmp = unlist(strsplit(crossmaps,"_X_|.Rds"))
    col1 = tmp[seq(1,length(tmp),2)]
    col2 = tmp[seq(2,length(tmp),2)]

    table = data.table(name1 = col1,name2 = col2,name1IsIn = (col1 %in% involved),name2IsIn = (col2 %in% involved),path = paste0(RRSDS,"Crossmaps/",crossmaps))
    table = subset(table,name1IsIn && name2IsIn)

    for(i in 1:length(table$name1)){
        same = readRDS(table$path[i])
        same1 = same[table$name1[i]][[1]]
        same2 = same[table$name2[i]][[1]]
        conts1 = subset(contigs,seqName == table$name1[i])
        conts2 = subset(contigs,seqName == table$name2[i])

        newConts = mkChimeras(conts1$start,conts1$end,conts1$coverage,conts2$start,conts2$end,conts2$coverage,start(same1),end(same1),start(same2),end(same2),conts1$seq,conts2$seq,conts1$seqName,conts2$seqName,conts1$rps,conts2$rps,minDist)

        contigs = subset(contigs,seqName != table$name1[i] && seqName != table$name2[i])
        contigs = data.table(start = c(contigs$start,newConts[[1]],newConts[[2]]),end = c(contigs$end,newConts[[3]],newConts[[4]]),seq = c(contigs$seq,newConts[[5]],newConts[[6]]),coverage = c(contigs$coverage,newConts[[7]],newConts[[8]]),seqName= c(contigs$seqName,newConts[[9]],newConts[[10]]),rps = c(contigs$rps,newConts[[11]],newConts[[12]]))
    }
    return(contigs)
}

#
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

# function to draw reads the necessary reads randomly
#
randomReads <- function(data,seqLength,coverage,meanWidht,repeatable,seed,redraw,sampleNr = 1,takeAll){


    whch = c()
    sampleNrs = c()
    for(i in 1:sampleNr){
        if(!takeAll){
            numberOfReads = as.integer((sum(seqLength) *coverage[i])/meanWidht)
        }
        else{
            numberOfReads = length(data$pos)
        }
        if(repeatable){
            set.seed(seed)
        }
        if(numberOfReads > length(data$pos) && !redraw){
            whch =c(whch,1:length(data$pos))
            sampleNrs = c(sampleNrs,rep(i,length(data$pos)))
        }
        else{
            whch = c(whch,sample(1:length(data$pos),numberOfReads,replace = redraw))
            sampleNrs = c(sampleNrs,rep(i,numberOfReads))
        }
    }
    names(whch) = sampleNrs
    whch = sort(whch)
    res = data[whch,]
    res = cbind(res, sampleNr = as.integer(names(whch)))
    return(res)
}


slidingWindowMaker <- function(vec){
    res = c()
    size = 1
    n = 1
    while(size <=  length(vec)/10000 ){
        size = size *10
    }
    if(size > 1){
        for(i in seq(1,length(vec),by = size/2)){
            if(i + size -1 > length(vec)){
                x = length(vec)
            }
            else{
                x = (i +size -1)
            }
            res[n] = mean(vec[i:x ])
            n = n +1
        }
    }
    return(res)
}
