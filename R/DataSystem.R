
#------------------------------------------ insert genomes into the datasystem ----------------------------------------------

addToDataSystem <- function(seqNames,bams = character(0),fasta,fastq1 = character(0),fastq2 = character(0),minMapq,bowtieOptions = "--no-unal"){

  if(!dir.exists("~/RealReadSimDS")){
    dir.create("~/RealReadSimDS")
    dir.create("~/RealReadSimDS/Crossmaps")
  }
  params = ScanBamParam(what = c("pos","qname","qwidth","rname"), mapqFilter = minMapq)
  lngth = c()
  mnwdth = c()
  for(i in 1:length(fasta)){
    print(seqNames[i])
    if(!is.inRRSDS(seqNames[i])){
      hasOut = FALSE
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

      seqs = scanBamHeader(bam)[[1]]$targets
      DS = dir("~/RealReadSimDS")

      reads = scanBam(bam,param =  params)
      reads = data.frame(pos = reads[[1]]$pos,width = reads[[1]]$qwidth,names = reads[[1]]$qname,seqName = reads[[1]]$rname)

      #system(paste("rm",bam))------------------- users choice? ---------------------------------------------------------------

      this = paste("~/RealReadSimDS/",seqNames[i],sep = "")

      dir.create(this)
      system(paste("cp ",fasta[i],this))
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
    fasta[i] = paste("~/RealReadSimDS/",seqNames[i],"/",gsub(".*/","",fasta[i]),sep = "")
    lngth[i] = readRDS(paste("~/RealReadSimDS/",seqNames[i],"/Length.Rds",sep = ""))
    mnwdth[i] = readRDS(paste("~/RealReadSimDS/",seqNames[i],"/MeanWidth.Rds",sep = ""))
  }
  return(data.frame(names = seqNames,totalLength = lngth,meanWidth = mnwdth,fastaName = fasta))
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

  povLength = readRDS(paste(pov,"/Length.Rds",sep =""))
  params = ScanBamParam(what = c("pos","qwidth","rname","qname"))
  #------------------------------------------------------ going through the datasystem------------------------------------------
  for(i in 1:length(DS)){
    if(DS[i] != "Crossmaps"){
      if(DS[i] != povName){
        otherLength = readRDS(paste(RRSDS,"/",DS[i],"/Length.Rds",sep =""))

        if(povLength > otherLength){#----------- choosing from which file to take the fasta and from which to take the reads --------
          fastafile = dir(pov)
          readfile = dir(paste(RRSDS,"/",DS[i],sep = ""))

          fastaPath = pov
          readPath = paste(RRSDS,"/",DS[i],sep = "")

          from = DS[i]
          to = povName
        }
        else{
          fastafile = dir(paste(RRSDS,"/",DS[i],sep = ""))
          readfile = dir(pov)

          fastaPath = paste(RRSDS,"/",DS[i],sep = "")
          readPath = pov

          from = povName
          to = DS[i]
        }
        length = readRDS(paste(fastaPath,"/Length.Rds",sep =""))
        meanWidth = readRDS(paste(readPath,"/MeanWidth.Rds",sep =""))
        readsToBe = readDNAStringSet(paste(readPath,"/",readfile[grep(".fasta",readfile)],sep = ""))

        #---------------------------------------- preparing the reads ------------------------------------------------
        newFasta = paste(Sys.getenv("HOME"),"/RealReadSimDS/",from,"/readsForMapX.fasta",sep = "")

        if(length(readfile[grep("readsForMapX.fasta",readfile)]) == 0){
          for(j in 1:length(readsToBe)){
            readPos = seq(0,width(readsToBe[j])-1,round(meanWidth*0.5))
            sequenceToFastaReads(readPos,toString(readsToBe[j]),meanWidth,newFasta,gsub(" .*","",names(readsToBe[j])))
          }
        }
        readInput = paste("-f ",newFasta,sep = "")
        #------------------------------ Using bowtie2 -----------------------------------------------------------------
        if(length(fastafile[grep(".*[^X].fasta",fastafile)]) > 1){
          stop("One of the fasta files given ends with an 'X' which is not permitted")
        }
        system(paste("bowtie2-build ",fastaPath,"/",fastafile[grep(".*[^X].fasta",fastafile)]," ",fastaPath,"/index",sep =""),ignore.stdout = TRUE,ignore.stderr = TRUE)

        system(paste("bowtie2 --no-unal -x ",fastaPath,"/index ",readInput," -S ",readPath,"/out.sam",sep =""),ignore.stdout = TRUE,ignore.stderr = TRUE)

        #------------------------------- readiying and saving the read data -------------------------------------------

        asBam(paste(readPath,"/out.sam",sep =""),paste(readPath,"/out",sep =""))

        mapped = scanBam(paste(readPath,"/out.bam",sep =""),param = params)
        mapped = data.frame(pos = mapped[[1]]$pos,width = mapped[[1]]$qwidth,name = mapped[[1]]$qname,seqs = mapped[[1]]$rname,stringsAsFactors = FALSE)
        mapped = subset(mapped, !is.na(pos))

        mapped2 = data.frame(pos = as.integer(gsub(".*_","",mapped$name)),width = mapped$width,name = mapped$name,seqs = gsub("_.*","",mapped$name),stringsAsFactors = FALSE)
        for(i in 1:length(unique(mapped$seqs))){
          sub1 = subset(mapped,seqs == unique(mapped$seqs)[i])
          print(sub1)
          if(length(sub1$pos) > 0){
            sub1 = IRanges(start = sub1$pos,width = sub1$width,names = sub1$name)
            for(j in 1:length(unique(mapped2$seqs))){
              sub2 = subset(mapped2,name %in% names(sub1))
              print(sub2)
              if(length(sub2$pos) > 0){
                sub2 = IRanges(start = sub2$pos,width = sub2$width,names = sub2$name)
                saveRDS(list(sub1,sub2),paste(RRSDS,"/Crossmaps/",names(x[[1]])[1],"_X_",names(x[[3]])[1],".Rds",sep = ""))
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
