
######################################### insert genomes into the datasystem #########################################

addToDataSystem <- function(seqNames,bams = character(0),fasta,fastq1 = character(0),fastq2 = character(0),minMapq,bowtieOptions = "--no-unal"){

  #--------------------------------------- if not done already initialize datasystem -------------------------------
  if(!dir.exists("~/RealReadSimDS")){
    dir.create("~/RealReadSimDS")
    dir.create("~/RealReadSimDS/Crossmaps")
  }
  #-----------------------------------------------------------------------------------------------------------

  params = ScanBamParam(what = c("pos","qname","qwidth","rname"), mapqFilter = minMapq)
  lngth = c()
  mnwdth = c()
  for(i in 1:length(fasta)){#------------------------- for all genomes --------------------------------------

    print(seqNames[i]) #####################################################################################

    if(!is.inRRSDS(seqNames[i])){
      hasOut = FALSE                    # mark wether there is output to be removed later

      #---------------------------- if not bams -> make bams | else use bams ------------------------------------

      if(length(bams) == 0){
        system(paste("bowtie2-build ",fasta," ~/RealReadSimDS/index"))

        if(length(fastq2) == 0){
          system(paste("bowtie2 ",paste(bowtieOptions,collapse = " ")," -x ~/RealReadSimDS/index -U ",fastq1," -S ~/RealReadSimDS/out.sam"),ignore.stdout = TRUE)
        }
        else{
          system(paste("bowtie2 ",paste(bowtieOptions,collapse = " ")," -x ~/RealReadSimDS/index -1 ",fastq1,"-2 ",fastq2," -S ~/RealReadSimDS/out.sam"),ignore.stdout = TRUE)
        }

        asBam("~/RealReadSimDS/out.sam","~/RealReadSimDS/out")
        system("rm ~/RealReadSimDS/out.sam")
        bam = "~/RealReadSimDS/out.bam"
        hasOut = TRUE
      }
      else{
        bam = bams[i]
      }

      #-----------------------------------------------------------------------------------------------------------

      seqs = scanBamHeader(bam)[[1]]$targets
      DS = dir("~/RealReadSimDS")

      reads = scanBam(bam,param =  params)
      reads = data.table(pos = reads[[1]]$pos,width = reads[[1]]$qwidth,names = reads[[1]]$qname,seqName = reads[[1]]$rname)

      #system(paste("rm",bam))         users choice?

      print(1)#####################################################################################

      this = paste("~/RealReadSimDS/",seqNames[i],sep = "")

      dir.create(this)
      system(paste("cp",fasta[i],this))
      length = 0

      for(j in 1:length(seqs)){
        data = subset(reads, seqName == names(seqs[j]))
        map = IRanges(start = data$pos,width = data$width,names = data$names)

        saveRDS(map,paste(this,"/",names(seqs[j]),".Rds",sep = "" ))
        length = length + unname(seqs[j])
      }

      saveRDS(length,paste(this,"/Length.Rds",sep = ""))
      meanWidth = mean(reads$width)
      mnwdth[i] = meanWidth
      saveRDS(meanWidth,paste(this,"/MeanWidth.Rds",sep = ""))

      if(hasOut){
        system("rm ~/RealReadSimDS/out*")
      }
      if(length(dir("~/RealReadSimDS")) > 2){
        crossMapRRSDS(this)
      }
    }
    if(substr(fasta[i],nchar(fasta[i]),nchar(fasta[i])) == "/"){
      fasta[i] = substr(fasta[i],nchar(fasta[i])-1,nchar(fasta[i])-1)
    }

    print(2)###################################################################################

    fasta[i] = paste("~/RealReadSimDS/",seqNames[i],"/",gsub(".*/","",fasta[i]),sep = "")
    lngth[i] = readRDS(paste("~/RealReadSimDS/",seqNames[i],"/Length.Rds",sep = ""))
    mnwdth[i] = readRDS(paste("~/RealReadSimDS/",seqNames[i],"/MeanWidth.Rds",sep = ""))

  }#----------------------------------------- end all genomes --------------------------------

  return(data.table(names = seqNames,totalLength = lngth,meanWidth = mnwdth,fastaName = fasta))
}


is.inRRSDS <- function(name){
  all = dir("~/RealReadSimDS")
  res = FALSE
  for(i in 1:length(all)){
    if(all[i] == name || paste("~/RealReadSimDS/",all[i],sep = "") == name || paste("~/RealReadSimDS/",all[i],"/",sep = "") == name){
      res = TRUE
    }
  }
  return(res)
}

#-------------------------- Mapping the reads of the shorter genome against the refseq of the larger one ----------------------
crossMapRRSDS <- function(pov){

  print("crossmap")

  if(!is.inRRSDS(pov)){
    stop("The requested file is not in RRSDS") ############# alle Fehler möglichkeiten (besonders das Löschen)
  }

  RRSDS = "~/RealReadSimDS"

  if(substr(pov,1,16) != "~/RealReadSimDS/"){
    pov = paste(RRSDS,"/",pov,sep ="")
  }
  if(substr(pov,nchar(pov),nchar(pov)) == "/"){
    pov = gsub("/$","",pov)
  }
  povName = gsub(".*/","",pov)

  DS = dir(RRSDS)

  fastas = c()
  nms = c()
  seqLength = c()
  meanWidths = c()
  dirName = c()

  n = 1
  for(i in 1:length(DS)){
    if(DS[i] != "Crossmaps"){
      dirPath = paste(RRSDS,"/",DS[i],"/",sep = "")
      dirCont = dir(dirPath)
      dirFasta = paste(dirPath,dirCont[grep(".*\\.fasta",dirCont)],sep = "")
      mnwdth = readRDS(paste(dirPath,"/",dirCont[grep("MeanWidth.Rds",dirCont)],sep = ""))
      tmp = readDNAStringSet(dirFasta)

      for(j in 1:length(tmp)){

        fastas[n] = dirFasta
        nms[n] = gsub(" .*","",names(tmp)[j])
        seqLength[n] = width(tmp)[j]
        meanWidths[n] = mnwdth
        dirName = DS[i]
        n = n +1
      }

    }
  }

  table = data.table::data.table(name = nms,fasta = fastas,length = seqLength,meanWidth = meanWidths,dirName = dirName)
  print(table$name)
  print(povName)
  povTable = subset(table,name == povName)
  table = subset(table,name != povName)

  povLength = readRDS(paste(pov,"/Length.Rds",sep =""))
  params = ScanBamParam(what = c("pos","qwidth","rname","qname"),mapqFilter = 40)

  print(1)######################################

  #------------------------------------------------------ going through the datasystem------------------------------------------
  for(p in 1:length(povTable$name)){
      for(j in 1:length(table$name)){

        #---- choosing from which file to take the fasta and from which to take the reads and from which to take the seq --------
        print("")#####################################
        print(p)######################################
        print(povTable[p,])############################
        print(table[j,])#################
        print("")#######################
        if(povTable$length[p] > table$length[j]){ # take the reads from the smaller one and the seq of the bigger one
          fastafile = dir(pov)
          readfile = dir(paste(RRSDS,"/",table$dirName[j],sep = ""))

          fastaPath = pov
          readPath = paste(RRSDS,"/",table$dirName[j],sep = "")

          from = table$dirName[j]
          to = povName
        }
        else{
          fastafile = dir(paste(RRSDS,"/",table$dirName[j],sep = ""))
          readfile = dir(pov)

          fastaPath = paste(RRSDS,"/",table$dirName[j],sep = "")
          readPath = pov

          from = povName
          to = table$dirName[j]
        }
        length = readRDS(paste(fastaPath,"/Length.Rds",sep =""))
        meanWidth = readRDS(paste(readPath,"/MeanWidth.Rds",sep =""))
        readsToBe = readDNAStringSet(paste(readPath,"/",readfile[grep(".fasta",readfile)],sep = ""))

         #----------------------------------------------------------------------------------------------------------------

        #---------------------------------------- preparing the reads ------------------------------------------------
        newFasta = paste(Sys.getenv("HOME"),"/RealReadSimDS/",from,"/readsForMapX.fasta",sep = "")

        if(length(readfile[grep("readsForMapX.fasta",readfile)]) == 0){
          for(j in 1:length(readsToBe)){
            readPos = seq(0,width(readsToBe[j])-1,round(meanWidth*0.5))
            RealReadSim::sequenceToFastaReads(readPos,toString(readsToBe[j]),meanWidth,newFasta,gsub(" .*","",names(readsToBe[j])))
            # cutting all seqs of the smaller genome in peaces and writing these down as reads with name "seqName_positionOfCut"
          }
        }
        readInput = paste("-f ",newFasta,sep = "")
        #------------------------------ mapping the reads onto the other seq with bowtie2 -----------------------------------------------------------------
        if(length(fastafile[grep(".*[^X].fasta",fastafile)]) > 1){
          stop("One of the fasta files given ends with an 'X' which is not permitted")
        }
        system(paste("bowtie2-build ",fastaPath,"/",fastafile[grep(".*[^X].fasta",fastafile)]," ",fastaPath,"/index",sep =""),ignore.stdout = TRUE,ignore.stderr = TRUE)

        system(paste("bowtie2 --no-unal -x ",fastaPath,"/index ",readInput," -S ",readPath,"/out.sam",sep =""),ignore.stdout = TRUE,ignore.stderr = TRUE)

          #------------------------------- readiying and saving the read data -------------------------------------------

        asBam(paste(readPath,"/out.sam",sep =""),paste(readPath,"/out",sep =""))

        mapped = scanBam(paste(readPath,"/out.bam",sep =""),param = params)
        mapped = data.table(pos = mapped[[1]]$pos,width = mapped[[1]]$qwidth,name = mapped[[1]]$qname,seqs = mapped[[1]]$rname)
        mapped = subset(mapped, !is.na(pos))

        mapped2 = data.table(pos = as.integer(gsub(".*_","",mapped$name)),width = mapped$width,name = mapped$name,seqs = gsub("_.*","",mapped$name))

        print(2)####################################################################
        print(mapped)##################################################################
        print(mapped2)##################################################################
        print(3)################################################################
        for(n in 1:length(unique(mapped$seqs))){
          sub1 = subset(mapped,seqs == unique(mapped$seqs)[n])

          print(sub1)###############################################################################

          if(length(sub1$pos) > 0){
            sub1IR = IRanges(start = sub1$pos,width = sub1$width,names = sub1$name)
            for(q in 1:length(unique(mapped2$seqs))){
              sub2 = subset(mapped2,name %in% names(sub1IR))

              print(sub2)###########################################################################

              if(length(sub2$pos) > 0){
                sub2IR = IRanges(start = sub2$pos,width = sub2$width,names = sub2$name)

                tmp1 = reduce(sub1IR)
                tmp2 = reduce(sub2IR)
                tmp1 = subset(tmp1,width > 0)############### minimum
                tmp2 = subset(tmp2,width > 0)############### minimum

                sub1IR = sub1IR[findOverlaps(sub1IR,tmp1)]
                sub2IR = sub2IR[findOverlaps(sub2IR,tmp2)]
                sub1IR = subset(sub1IR,names %in% names(sub2IR))
                sub2IR = subset(sub2IR,names %in% names(sub1IR))

                if(length(width(sub2IR)) > 0){
                  toSave = list(sub1IR,sub2IR)
                  names(toSave) = c(sub1$seqs[1],sub2$seqs[1])
                  saveRDS(toSave,paste(RRSDS,"/Crossmaps/",sub2$seqs[1],"_X_",sub1$seqs[1],".Rds",sep = ""))
                }
              }
            }
          }
        }

          #-------------------------------- deleting everything that is of no further use ------------------------------------------

        system(paste("rm ",fastaPath,"/*index*",sep = ""))
        system(paste("rm ",readPath,"/out.sam",sep = ""))
        system(paste("rm ",readPath,"/out.bam*",sep = ""))
    }
  }
}


IRangesTranslator <- function(reads1,reads2){
  y = Sys.time()
  contigs1 = reduce(reads1)
  contigs2 = reduce(reads2)
  names(contigs1) <- c(1:length(contigs1))
  names(contigs2) <- c(1:length(contigs2))
  con1 = list()
  con2 = list()

  rang1 = IRanges()
  rang2 = IRanges()

  for(i in 1:length(contigs1)){
    reads1i = subsetByOverlaps(reads1,contigs1[i])
    reads2i = reads2[which(names(reads2) %in% names(reads1i))]
    contigs2i = subsetByOverlaps(contigs2,reads2i)
    contigs1i = reduce(reads1i)
    tmp1 = IRanges()
    tmp2 = IRanges()
    n = 1
    switch = FALSE
    while(length(contigs2i) > 1 || length(contigs1i) >1 || switch){
      if(length(contigs1i) >1){
        reads1n = subsetByOverlaps(reads1i,contigs1i[n])
        reads2n = subset(reads2i,names %in% names(reads1n))
        contigs2n = subsetByOverlaps(contigs2i,reads2n)

      }
      else{

      }
    }
  }
  res = list(rang1,rang2)
  print(Sys.time() -y)
  return(res)
}
####
