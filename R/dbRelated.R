# auf eis erstmal

setUpRealReadSimDB <-function(overwrite = FALSE){
  if(file.exists("~/RealReadSimDB.sqlite")){
    if(overwrite == FALSE){
      return()
    }
    else{
      system("rm ~/RealReadSimDB.sqlite")
    }
  }
  system("touch ~/RealReadSim.sqlite")

  db = dbConnect(SQLite(),"~/RealReadSimDB.sqlite")
  dbCreateTable(db,"Reads",c(Pos = "integer",Width = "integer",SeqID = "String"))
  dbCreateTable(db,"Genomes",c(Seq = "String",Length = "integer",SeqID = "String",Fasta = "String",Bam = "String"))
  dbCreateTable(db,"Intersections",c(GenomeFrom = "String",ReadsFrom = "String",Start = "Integer",Stop = "integer"))
  dbDisconnect(db)
}

addToRealReadSimDB <- function(bam,fasta){

  param = ScanBamParam(what = c("pos","qwidth","rname"))
  Reads = scanBam(bam,param = param)
  typeof(Reads$rname)
  Reads = data.frame(Pos = Reads[[1]]$pos,Width = Reads[[1]]$qwidth,SeqID = as.character(Reads[[1]]$rname),stringsAsFactors = FALSE)

  print(1)

  Genomes = readDNAStringSet(fasta)
  Genomes = data.frame(Seq = as.character(Genomes),Length = width(Genomes),SeqID = gsub(" .*$","",names(Genomes)),Fasta = fasta,Bam = bam, stringsAsFactors = FALSE)

  print(2)
  db = dbConnect(SQLite(),"~/RealReadSimDB.sqlite")

  dbAppendTable(db,"Reads",Reads)
  dbAppendTable(db,"Genomes",Genomes)

  if( !(length(dbGetQuery(db,"SELECT SeqID FROM Genomes")) < 1) ){

  }

  dbDisconnect(db)
}

mapAllOnOneAndOneOnAll <- function(fastq1 , fastq2 = "", fasta,length,mWidth){



  DB = dbConnect(SQLite(),"~/RealReadSimDB.sqlite")
  lengths = dbGetQuery(DB,"SELECT Length FROM Genomes")
  names = dbGetQuery(DB,"SELECT SeqID FROM Genomes")
  widths = c()
  for(i in 1:length(names)){
    widths[i] = mean(dbGetQuery(DB,paste("SELECT Width FROM Reads WHERE SeqID ==",names[i],seq = " ")))
  }
  nmbrs1 = as.integer(2*(lenght / widhts))

  nmbrs2 = as.integer(2*(lengths/mWidth))

  x = 1

  fastas = dbGetQuery(DB,"SELECT Fasta FROM Genomes")
  for(i in 1:length(fastas)){

    system(paste("bowtie2-build",fatsas[i],"1731070519",sep = " "))
    if(fastq2 != ""){
      x = 2
      writeFastq(FastqSampler(fastq2,n = (nmbrs1[i]/x)),"Newfastq2.fastq")
      fastqToBowtie = "-1 Newfastq.fastq -2 Newfastq2.fastq"
    }
    else{
      fastqToBowtie = "-U Newfastq.fastq"
    }

    writeFastq(FastqSampler(fastq1,n = (nmbrs1[i]/x)),"Newfastq.fastq")


    system(paste("bowtie2 --no-unal -x 1731070519",fastqToBowtie,"-S NewSam.sam",sep = " "))

    asBam("NewSam.sam","NewBam.bam")


  }




}
