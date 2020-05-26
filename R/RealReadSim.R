assignInNamespace("cedta.override",
                  c(data.table:::cedta.override,"RealReadSim"),
                  "data.table")

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

realReadSim <- function(filenames_csv = "",coverage = 0,takeAll = TRUE,nrOfSamples = 1,bowtieOptions = c("--no-unal","-p 3"),humanReadable = FALSE,readAsBams = TRUE, minMapq = 40, redraw = FALSE, repeatable = TRUE,seed = 0,minContigLength = 500,minDist = 0.99,minIdenticalLength = 2000,metagenomeDir = "~/RealReadSimDS"){

    library(Rsamtools)
    library(data.table)
    library(seqinr)

    require(data.table)

    if(substr(metagenomeDir,nchar(metagenomeDir),nchar(metagenomeDir)) == "/"){
        metagenomeDir = substr(metagenomeDir,1,nchar(metagenomeDir)-1)
    }

    starttime = Sys.time()


    covAt = 1
    if(filenames_csv != ""){
        filenames = read.table(filenames_csv,header = TRUE,sep = ",",stringsAsFactors = FALSE)
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
            addToDataSystem(seqNames,bam = filenames[,2],fasta = filenames[,1],minMapq = minMapq,bowtieOptions = bowtieOptions,minIdL = minIdenticalLength,metagenomeDir = metagenomeDir)
        }
        else{
            if(length(filenames) > 2){
                addToDataSystem(seqNames,fasta = filenames[,1],fastq1 = filenames[,2],fastq2 = filenames[,3],minMapq = minMapq,bowtieOptions =  bowtieOptions,minIdL = minIdenticalLength,metagenomeDir = metagenomeDir )
            }
            else{
                addToDataSystem(seqNames,fasta = filenames[,1],fastq1 = filenames[,2],minMapq = minMapq,bowtieOptions =  bowtieOptions,minIdL = minIdenticalLength,metagenomeDir = metagenomeDir)
            }
        }

        #------------------------------------------------------------------------------------------------------------
    }
    catalogue = readRDS(paste0(metagenomeDir,"/DSTable.Rds"))

    #-------------------------------- assembly --------------------------------------------------------
    RRSDS = metagenomeDir

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
        fullData = randomReads(fullData,partialCatalogue$totalLength[1],coverage = coverage[covAt,],meanWidht = partialCatalogue$meanWidth[1],repeatable,seed,redraw,nrOfSamples,takeAll)
        for(j in 1:length(partialCatalogue$name)){

            data = subset(fullData,seq == partialCatalogue$name[j])
            if(length(data$pos) > 0){
                print(paste(length(data$pos)," reads have been selected"))###########################################
                contigs = evalCoverage(data$pos, data$width,data$sampleNr, partialCatalogue$len[j],partialCatalogue$minOverlap[j],minContigLength,nrOfSamples)
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

    res = coAssembleRRSDS(temp,minDist,metagenomeDir)
    res$length = (res$end-res$start)+1
    res = subset(res,length >= minContigLength)
    res$covVec = calcCovVec(res$rps,res$length)
    #---------------------------------------------------------------------------------------------------
    if(humanReadable){
        #makeHumanReadable(cov,res)
    }

    print(Sys.time() -starttime)#################################################################
    return(res)
}

# function to simulate coassembly using coAssembly.cpp
#
coAssembleRRSDS <- function(contigs,minDist,metagenomeDir){

    involved = unique(contigs$seqName)

    crossmaps = dir(paste0(metagenomeDir,"/Crossmaps"))
    if(length(crossmaps) != 0){
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

            contigs = contigs[seqName != table$name1[i] && seqName != table$name2[i],]
            contigs = data.table(start = c(contigs$start,newConts[[1]],newConts[[2]]),end = c(contigs$end,newConts[[3]],newConts[[4]]),seq = c(contigs$seq,newConts[[5]],newConts[[6]]),coverage = c(contigs$coverage,newConts[[7]],newConts[[8]]),seqName= c(contigs$seqName,newConts[[9]],newConts[[10]]),rps = c(contigs$rps,newConts[[11]],newConts[[12]]))
            print(contigs[c(1,2,5,6)])
        }
    }
    return(contigs)
}

#
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


######################################### insert genomes into the datasystem #########################################

addToDataSystem <- function(seqNames,bams = character(0),fasta,fastq1 = character(0),fastq2 = character(0),minMapq,bowtieOptions = "--no-unal",minIdL,metagenomeDir){

    #--------------------------------------- if not done already initialize datasystem -------------------------------
    initTable = FALSE
    if(!dir.exists(metagenomeDir)){
        dir.create(metagenomeDir)
        dir.create(paste0(metagenomeDir,"/Crossmaps"))
        initTable = TRUE
    }
    #-----------------------------------------------------------------------------------------------------------
    library(data.table)

    params = ScanBamParam(what = c("pos","qwidth","rname"), mapqFilter = minMapq)

    lngth = c()
    seqFastas = c()
    seqNms = c()
    seqLngths = c()
    seqMnWdths = c()
    seqDir = c()
    seqData = c()
    seqMinOv = c()
    n = 1
    hasNew = FALSE

    for(i in 1:length(fasta)){#------------------------- for all genomes --------------------------------------
        print(paste0("genome nr: ",i))
        if(!is.inRRSDS(seqNames[i],metagenomeDir)){
            hasOut = FALSE                    # mark wether there is output to be removed later
            hasNew = TRUE

            #---------------------------- if not bams -> make bams | else use bams ------------------------------------

            if(length(bams) == 0){
                system(paste0("bowtie2-build ",fasta[i]," ",metagenomeDir,"/index"),ignore.stdout = TRUE)

                if(length(fastq2) == 0){
                    system(paste0("bowtie2 ",paste(bowtieOptions,collapse = " ")," -x ",metagenomeDir,"/index -U ",fastq1[i]," -S ",metagenomeDir,"/out.sam"))#,ignore.stdout = TRUE)#,ignore.stderr = TRUE)
                }
                else{
                    system(paste0("bowtie2 ",paste(bowtieOptions,collapse = " ")," -x ",metagenomeDir,"/index -1 ",fastq1[i],"-2 ",fastq2[i]," -S ",metagenomeDir,"/out.sam"))#,ignore.stdout = TRUE)#,ignore.stderr = TRUE)
                }

                #asBam("~/RealReadSimDS/out.sam","~/RealReadSimDS/out")
                system(paste0("samtools view -bS ",metagenomeDir,"/out.sam > ",metagenomeDir,"/out.bam"))
                system(paste0("rm ",metagenomeDir,"/out.sam"))
                bam = paste0(metagenomeDir,"/out.bam")
                hasOut = TRUE
            }
            else{
                bam = bams[i]
            }

            #-----------------------------------------------------------------------------------------------------------

            seqs = scanBamHeader(bam)[[1]]$targets
            DS = dir(metagenomeDir)

            reads = scanBam(bam,param =  params)
            if(length(reads[[1]]$pos) > 0){
                reads = data.table(pos = reads[[1]]$pos,width = reads[[1]]$qwidth,seq = reads[[1]]$rname,key = c("seq","pos"))
                #system(paste("rm",bam))         users choice?

                this = paste0(metagenomeDir,"/",seqNames[i])

                dir.create(this)

                dataPath = paste0(this,"/",names(seqs[which.max(seqs)]),".Rds")
                saveRDS(reads,dataPath)

                new.name.fasta = strsplit(fasta[i],"/")
                new.name.fasta = new.name.fasta[[1]][length(new.name.fasta[[1]])]
                system(paste0("cp ",fasta[i]," ",this))
                file.rename(paste0(this,"/",new.name.fasta),paste0(this,"/bin.fasta"))

                if(length(seqs) > 1){
                    sep.fastas(fasta[i],this)
                }
                else{
                    system(paste0("cp ",fasta[i]," ",this))
                    file.rename(paste0(this,"/",new.name.fasta),paste0(this,"/",names(seqs[1]),".fasta"))
                }

                length = 0
                meanWidth = mean(reads$width)
                sequences = seqinr::read.fasta(fasta[i],as.string = TRUE)
                names(sequences) = gsub(" .*","",names(sequences))

                readPos = list()
                readsToBe = c()

                for(j in 1:length(seqs)){
                    seqLngths[n] = unname(seqs[j])
                    seqNms[n] = names(seqs[j])
                    seqMnWdths[n] = meanWidth

                    seqFastas[n] = paste0(this,"/",names(seqs[j]),".fasta")
                    seqDir[n] = this
                    lngth[n] = sum(seqs)

                    seqData[n] = dataPath
                    subSequence = sequences[names(sequences) == names(seqs[j])][[1]]
                    seqMinOv[n] = calcMinOverlap(subSequence,meanWidth)

                    readPos[[j]] = seq(0,unname(seqs[j])-1,round(meanWidth*0.5))
                    readsToBe[j] = subSequence

                    n = n +1
                }
                if(substr(this,1,1) == "~"){
                    this = paste0(Sys.getenv("HOME"),substr(this,2,nchar(this)))
                }
                newFasta = paste0(this,"/readsForMapX.fasta")
                sequenceToFastaReads(readPos,readsToBe,meanWidth,newFasta,names(seqs))

            }

            if(hasOut){
                system(paste0("rm ",metagenomeDir,"/out*"))
            }
        }

    }#----------------------------------------- end all genomes --------------------------------

    tmpTable = data.table(name = seqNms,len = seqLngths,meanWidth = seqMnWdths,fasta = seqFastas,data = seqData,dir = seqDir,totalLength = lngth,minOverlap = seqMinOv,isCrossmapped = rep(FALSE,length(seqNms)))
    if(hasNew){
        if(initTable){
            saveRDS(tmpTable,paste0(metagenomeDir,"/DSTable.Rds"))
            initTable = FALSE
        }
        else{
            table = readRDS(paste0(metagenomeDir,"/DSTable.Rds"))
            tmpTable = rbind(tmpTable,table)
            saveRDS(tmpTable,paste0(metagenomeDir,"/DSTable.Rds"))
        }
    }

    crossMapRRSDS(minIdL,metagenomeDir)

}


sep.fastas <- function(fasta,dir){
    seqs = readDNAStringSet(fasta)
    for(i in 1:length(seqs)){
        seqinr::write.fasta(sequences = toString(seqs[i]),names =gsub(" .*","",names(seqs[i])),file.out = paste0(dir,"/",gsub(" .*","",names(seqs[i])),".fasta") )
    }
}


is.inRRSDS <- function(name,metagenomeDir){
    all = dir(metagenomeDir)
    res = FALSE
    for(i in 1:length(all)){
        if(all[i] == name || paste(metagenomeDir,all[i],sep = "") == name || paste(metagenomeDir,all[i],"/",sep = "") == name){
            res = TRUE
        }
    }
    return(res)
}

# #-------------------------- Mapping the reads of the shorter genome against the refseq of the larger one ----------------------
# crossMapRRSDS <- function(minL){
#   RRSDS = "~/RealReadSimDS"
#   table = readRDS(paste0(RRSDS,"/DSTable.Rds"))
#
#   query = which(table$isCrossmapped == FALSE && table$length > minL)
#
#   params = ScanBamParam(what = c("pos","qwidth","qname"),mapqFilter = 40)
#
#   jj = 1:length(table$name)
#   ii = 0
#   #------------------------------------------------------ going through the datasystem------------------------------------------
#   for(p in query){
#       jj = jj[-(p -ii)]
#       print(p)
#       for(j in jj){
#         print(j)
#         #---- choosing from which file to take the fasta and from which to take the reads and from which to take the seq --------
#         if(table$length[j] > table$length[p]){ # take the reads from the smaller one and the seq of the bigger one
#           from = table[j,]
#           to = table[p,]
#         }
#         else{
#           from = table[p,]
#           to = table[j,]
#         }
#         readsToBe = readDNAStringSet(from$fasta)
#
#         #----------------------------------------------------------------------------------------------------------------
#
#         #---------------------------------------- preparing the reads ------------------------------------------------
#         newFasta = paste0(from$dir,"/",from$name,"readsForMapX.fasta")
#         readfile = dir(from$dir)
#         if(length(readfile[grep(paste0(from$name,"readsForMapX.fasta"),readfile)]) == 0){
#           readPos = seq(0,width(readsToBe)-1,round(from$meanWidth*0.5))
#           sequenceToFastaReads(readPos,toString(readsToBe),from$meanWidth,newFasta,from$name)
#         }
#         # cutting all seqs of the smaller genome in peaces and writing these down as reads with name "seqName_positionOfCut"
#         readInput = paste0("-f ",newFasta)
#         #------------------------------ mapping the reads onto the other seq with bowtie2 -----------------------------------------------------------------
#
#         system(paste0("bowtie2-build --threads 3 ",to$fasta," ",to$dir,"/index"),ignore.stdout = TRUE)
#         system(paste0("bowtie2 -p 3 --no-unal -x ",to$dir,"/index ",readInput," -S ",from$dir,"/out.sam"),ignore.stdout = TRUE)
#
#         system(paste0("samtools view -bS ",from$dir,"/out.sam > ",from$dir,"/out.bam"))
#         #------------------------------- readiying and saving the read data -------------------------------------------
#         mapped = scanBam(paste0(from$dir,"/out.bam"),param = params)
#
#         map = data.table(start1 = mapped[[1]]$pos,end1 = mapped[[1]]$pos + mapped[[1]]$qwidth -1,start2 = as.integer(gsub(".*_","",mapped[[1]]$qname)),end2 =as.integer(gsub(".*_","",mapped[[1]]$qname)) + mapped[[1]]$qwidth -1,key = "start1" )
#
#         if(length(map$start1) > 2){
#           sameConts = getIdenticalSeqs(map$start1,map$end1,map$start2,map$end2,minL)
#           if(length(sameConts[[1]]) > 0){
#             conts1 = IRanges(start = sameConts[[2]],end = sameConts[[3]],names = sameConts[[1]])
#             conts2 = IRanges(start = sameConts[[4]],end = sameConts[[5]],names = sameConts[[1]])
#
#             toSave = list(conts2,conts1)
#             names(toSave) = c(from$name,to$name)
#             saveRDS(toSave,paste0(RRSDS,"/Crossmaps/",from$name,"_X_",to$name,".Rds"))
#           }
#         }
#         #-------------------------------- deleting everything that is of no further use ------------------------------------------
#
#         system(paste0("rm ",to$dir,"/*index*"))
#         system(paste0("rm ",from$dir,"/out.sam"))
#         system(paste0("rm ",from$dir,"/out.bam*"))
#       }
#       ii = ii +1
#       table$isCrossmapped[p] = TRUE
#   }
#   saveRDS(table,paste0(Sys.getenv("HOME"),"/RealReadSimDS/DSTable.Rds"))
# }



#-------------------------- Mapping the reads of the shorter genome against the refseq of the larger one ----------------------
crossMapRRSDS <- function(minL,metagenomeDir){

    table = readRDS(paste0(metagenomeDir,"/DSTable.Rds"))
    bins = table[, .("length" = sum(len),"isCrossmapped" = all(isCrossmapped)),by = "dir"]

    params = ScanBamParam(what = c("pos","qwidth","qname","rname"),mapqFilter = 40)
    #------------------------------------------------------ going through the datasystem------------------------------------------
    for(p in 1:(length(bins$dir)-1)){
        print(p)#######################################
        if(bins[p,]$isCrossmapped == FALSE){
            for(j in (p+1):length(bins$dir)){

                #---- choosing from which file to take the fasta and from which to take the reads and from which to take the seq --------
                if(bins$length[j] > bins$length[p]){ # take the reads from the smaller one and the seq of the bigger one
                    from = bins[j,]
                    to = bins[p,]
                }
                else{
                    from = bins[p,]
                    to = bins[j,]
                }

                # cutting all seqs of the smaller genome in peaces and writing these down as reads with name "seqName_positionOfCut"
                readInput = paste0("-f ",from$dir,"/readsForMapX.fasta")
                #------------------------------ mapping the reads onto the other seq with bowtie2 -----------------------------------------------------------------

                system(paste0("bowtie2-build --threads 3 ",to$dir,"/bin.fasta ",to$dir,"/index"),ignore.stdout = TRUE)
                system(paste0("bowtie2 -p 3 --no-unal -x ",to$dir,"/index ",readInput," -S ",from$dir,"/out.sam"),ignore.stdout = TRUE)

                system(paste0("samtools view -bS ",from$dir,"/out.sam > ",from$dir,"/out.bam"))
                #------------------------------- readiying and saving the read data -------------------------------------------
                mapped = scanBam(paste0(from$dir,"/out.bam"),param = params)

                map = data.table(names1 = as.character(mapped[[1]]$rname),
                                 start1 = mapped[[1]]$pos,
                                 end1 = mapped[[1]]$pos + mapped[[1]]$qwidth -1,
                                 names2 = gsub(";;.*","",mapped[[1]]$qname),
                                 start2 = as.integer(gsub(".*;;","",mapped[[1]]$qname)),
                                 end2 =as.integer(gsub(".*;;","",mapped[[1]]$qname)) + mapped[[1]]$qwidth -1,
                                 key = c("names1","names2"),
                                 stringsAsFactors = FALSE )

                map = map[,.(start1 = list(start1),end1 = list(end1),start2 = list(start2),end2 = list(end2)),by  = list(names1,names2)]

                if(length(map$names1) > 0){
                    sameConts = getIdenticalSeqsList(map$names1,map$start1,map$end1,map$names2,map$start2,map$end2,minL)
                    for(i in 1:length(sameConts)){
                        if(length(sameConts[[i]][[1]]) > 0){
                            conts1 = IRanges(start = sameConts[[i]][[2]],end = sameConts[[i]][[3]],names = sameConts[[i]][[1]])
                            conts2 = IRanges(start = sameConts[[i]][[4]],end = sameConts[[i]][[5]],names = sameConts[[i]][[1]])

                            toSave = list(conts2,conts1)
                            names(toSave) = c(sameConts[[i]][[6]],sameConts[[i]][[7]])
                            saveRDS(toSave,paste0(metagenomeDir,"/Crossmaps/",sameConts[[i]][[6]],"_X_",sameConts[[i]][[7]],".Rds"))
                        }
                    }
                }
                #-------------------------------- deleting everything that is of no further use ------------------------------------------

                system(paste0("rm ",to$dir,"/*index*"))
                system(paste0("rm ",from$dir,"/out.sam"))
                system(paste0("rm ",from$dir,"/out.bam*"))
            }
            table[dir == bins[p,]$dir]$isCrossmapped = TRUE
        }
    }
    saveRDS(table,paste0(metagenomeDir,"/DSTable.Rds"))
}

####

