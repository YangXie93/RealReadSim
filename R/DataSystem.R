
######################################### insert genomes into the datasystem #########################################

addToDataSystem <- function(seqNames,bams = character(0),fasta,fastq1 = character(0),fastq2 = character(0),minMapq,bowtieOptions = "--no-unal"){

  #--------------------------------------- if not done already initialize datasystem -------------------------------
  initTable = FALSE
  if(!dir.exists("~/RealReadSimDS")){
    dir.create("~/RealReadSimDS")
    dir.create("~/RealReadSimDS/Crossmaps")
    initTable = TRUE
  }
  #-----------------------------------------------------------------------------------------------------------


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
    if(!is.inRRSDS(seqNames[i])){
      hasOut = FALSE                    # mark wether there is output to be removed later
      hasNew = TRUE

      #---------------------------- if not bams -> make bams | else use bams ------------------------------------

      if(length(bams) == 0){
        system(paste("bowtie2-build ",fasta," ~/RealReadSimDS/index"),ignore.stdout = TRUE,ignore.stderr = TRUE)

        if(length(fastq2) == 0){
          system(paste("bowtie2 ",paste(bowtieOptions,collapse = " ")," -x ~/RealReadSimDS/index -U ",fastq1[i]," -S ~/RealReadSimDS/out.sam"),ignore.stdout = TRUE)#,ignore.stderr = TRUE)
        }
        else{
          print(paste(fasta[i],fastq1[i],fastq2[i]))
          system(paste("bowtie2 ",paste(bowtieOptions,collapse = " ")," -x ~/RealReadSimDS/index -1 ",fastq1[i],"-2 ",fastq2[i]," -S ~/RealReadSimDS/out.sam"),ignore.stdout = TRUE)#,ignore.stderr = TRUE)
        }

        #asBam("~/RealReadSimDS/out.sam","~/RealReadSimDS/out")
        system("samtools view -bS ~/RealReadSimDS/out.sam > ~/RealReadSimDS/out.bam")
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
      reads = data.table(pos = reads[[1]]$pos,width = reads[[1]]$qwidth,seq = reads[[1]]$rname,key = c("seq","pos"))
      #system(paste("rm",bam))         users choice?

      this = paste0(Sys.getenv("HOME"),"/RealReadSimDS/",seqNames[i])

      dir.create(this)

      dataPath = paste0(this,"/",names(seqs[which.max(seqs)]),".Rds")
      saveRDS(reads,dataPath)

      if(length(seqs) > 1){
        sep.fastas(fasta[i],this)
      }
      else{
        new.name.fasta = strsplit(fasta[i],"/")
        new.name.fasta = new.name.fasta[[1]][length(new.name.fasta[[1]])]
        system(paste0("cp ",fasta[i]," ",this))
        file.rename(paste0(this,"/",new.name.fasta),paste0(this,"/",seqNames[i],".fasta"))
      }

      length = 0
      meanWidth = mean(reads$width)
      print(meanWidth)###########################
      sequences = readDNAStringSet(fasta[i])
      names(sequences) = gsub(" .*","",names(sequences))

      for(j in 1:length(seqs)){
        seqLngths[n] = unname(seqs[j])
        seqNms[n] = names(seqs[j])
        seqMnWdths[n] = meanWidth

        seqFastas[n] = paste0(this,"/",names(seqs[j]),".fasta")
        seqDir[n] = this
        lngth[n] = sum(seqs)

        seqData[n] = dataPath
        subSequences = sequences[names(sequences) == names(seqs[j])]
        seqMinOv[n] = calcMinOverlap(toString(subSequences),meanWidth)

        n = n +1
      }


      if(hasOut){
        system("rm ~/RealReadSimDS/out*")
      }
    }

  }#----------------------------------------- end all genomes --------------------------------

  tmpTable = data.table(name = seqNms,length = seqLngths,meanWidth = seqMnWdths,fasta = seqFastas,data = seqData,dir = seqDir,totalLength = lngth,minOverlap = seqMinOv,isCrossmapped = rep(FALSE,length(seqNms)))
  if(hasNew){
    if(initTable){
      saveRDS(tmpTable,"~/RealReadSimDS/DSTable.Rds")
      initTable = FALSE
    }
    else{
      table = readRDS("~/RealReadSimDS/DSTable.Rds")
      tmpTable = rbind(table,tmpTable)
      saveRDS(tmpTable,"~/RealReadSimDS//DSTable.Rds")
    }
  }
  else{
    tmpTable = readRDS(paste0(Sys.getenv("HOME"),"/RealReadSimDS/DSTable.Rds"))
  }

  crossMapRRSDS()

  return(tmpTable)
}


sep.fastas <- function(fasta,dir){
  seqs = readDNAStringSet(fasta)
  for(i in 1:length(seqs)){
    seqinr::write.fasta(sequences = toString(seqs[i]),names =gsub(" .*","",names(seqs[i])),file.out = paste0(dir,"/",gsub(" .*","",names(seqs[i])),".fasta") )
  }
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
crossMapRRSDS <- function(){
  RRSDS = "~/RealReadSimDS"
  table = readRDS(paste0(RRSDS,"/DSTable.Rds"))

  query = which(table$isCrossmapped == FALSE)

  params = ScanBamParam(what = c("pos","qwidth","qname"),mapqFilter = 40)

  jj = 1:length(table$name)
  ii = 0
  #------------------------------------------------------ going through the datasystem------------------------------------------
  for(p in query){
      jj = jj[-(p -ii)]
      for(j in jj){
        #---- choosing from which file to take the fasta and from which to take the reads and from which to take the seq --------
        if(table$length[j] > table$length[p]){ # take the reads from the smaller one and the seq of the bigger one
          from = table[j,]
          to = table[p,]
        }
        else{
          from = table[p,]
          to = table[j,]
        }
        readsToBe = readDNAStringSet(from$fasta)

        #----------------------------------------------------------------------------------------------------------------

        #---------------------------------------- preparing the reads ------------------------------------------------
        newFasta = paste0(from$dir,"/",from$name,"readsForMapX.fasta")
        readfile = dir(from$dir)
        if(length(readfile[grep(paste0(from$name,"readsForMapX.fasta"),readfile)]) == 0){
          readPos = seq(0,width(readsToBe)-1,round(from$meanWidth*0.5))
          RealReadSim::sequenceToFastaReads(readPos,toString(readsToBe),from$meanWidth,newFasta,from$name)
        }
        # cutting all seqs of the smaller genome in peaces and writing these down as reads with name "seqName_positionOfCut"
        readInput = paste0("-f ",newFasta)
        #------------------------------ mapping the reads onto the other seq with bowtie2 -----------------------------------------------------------------

        system(paste0("bowtie2-build ",to$fasta," ",to$dir,"/index"),ignore.stdout = TRUE,ignore.stderr = TRUE)
        system(paste0("bowtie2 --no-unal -x ",to$dir,"/index ",readInput," -S ",from$dir,"/out.sam"),ignore.stdout = TRUE,ignore.stderr = TRUE)

        system(paste0("samtools view -bS ",from$dir,"/out.sam > ",from$dir,"/out.bam"))
        #------------------------------- readiying and saving the read data -------------------------------------------
        mapped = scanBam(paste0(from$dir,"/out.bam"),param = params)

        map = data.table(start1 = mapped[[1]]$pos,end1 = mapped[[1]]$pos + mapped[[1]]$qwidth -1,start2 = as.integer(gsub(".*_","",mapped[[1]]$qname)),end2 =as.integer(gsub(".*_","",mapped[[1]]$qname)) + mapped[[1]]$qwidth -1,key = "start1" )

        if(length(map$start1) > 2){
          sameConts = getIdenticalSeqs(map$start1,map$end1,map$start2,map$end2,2000)
          if(length(sameConts[[1]]) > 0){
            conts1 = IRanges(start = sameConts[[2]],end = sameConts[[3]],names = sameConts[[1]])
            conts2 = IRanges(start = sameConts[[4]],end = sameConts[[5]],names = sameConts[[1]])

            toSave = list(conts2,conts1)
            names(toSave) = c(from$name,to$name)
            saveRDS(toSave,paste0(RRSDS,"/Crossmaps/",from$name,"_X_",to$name,".Rds"))
          }
        }
        #-------------------------------- deleting everything that is of no further use ------------------------------------------

        system(paste0("rm ",to$dir,"/*index*"))
        system(paste0("rm ",from$dir,"/out.sam"))
        system(paste0("rm ",from$dir,"/out.bam*"))
      }
      ii = ii +1
      table$isCrossmapped[p] = TRUE
  }
  saveRDS(table,paste0(Sys.getenv("HOME"),"/RealReadSimDS/DSTable.Rds"))
}

####
