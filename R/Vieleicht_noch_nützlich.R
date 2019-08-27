
####################################################### Absolute Pfade !!!!!!!!!!!
mapOneToAll <- function(fastq1,fastq2,pov,meanWidth){

  params = ScanBamParam(what = c("pos","qwidth","rname","qname"))
  DS = dir("~/RealReadSimDS")
  for(i in 1:length(DS)){
    if(DS[i] != pov){

      if(dir.exists(paste("~/RealReadSimDS/",DS[i]))){

        folder = dir(paste("~/RealReadSimDS/",DS[i]))
        fasta = folder[grep(".*.fasta$",folder)]

        if(length(fasta) < 1){
          print(paste("keine fasta Datei in ",DS[i]))
          return(FALSE)
        }
        ############################################ Anzahl der Reads die gemapped werden sollen (ca. 2 fache Deckung)
        system(paste("bowtie2-build",fasta, "index",seq = " "))
        lengs = fasta.seqlengths(fasta)
        num = as.integer((sum(lengs)/meanWidth) * 2)
        ############################################

        ############################################ sam erzeugen  (Bowtie options ?)
        if(fastq2 == ""){
          writeFastq(FastqSampler(fastq1,n = num),"rand.fastq")
          system(paste("bowtie2 --no-unal -x index -U rand.fastq -S res.sam", sep = " "))
        }
        else{
          writeFastq(FastqSampler(fastq1,n = num/2),"rand.fastq")
          writeFastq(FastqSampler(fastq2,n = num/2),"rand2.fastq")
          system(paste("bowtie2 --no-unal -x index -1 rand.fastq -2 rand2.fastq -S res.sam"))
        }
        ##########################################

        ########################################## Bam zu IRanges und dann speichern   (Lösung für meherere Sequenzen in einer Fasta)
        as.bam("res.sam","res.bam")
        system("rm res.sam")

        mapped = scanBam("res.bam",params)
        mapped = data.frame(pos = mapped[[1]]$pos,width = mapped[[1]]$qwidth,name = mapped[[1]]$qname,seqs = mapped[[1]]$rname)
        seq = unique(mapped$seqs)

        for(i in 1:length(seq)){
          temp = subset(mapped, seqs == seq[i])
          map = IRanges(start = mapped$pos,width = mapped$width,names = qname)
          save(map,paste("~/RealReadSimDS/",pov,"/",seq[i],"_X_",pov,"CrossMap"))
        }
        ######################################### Die Gescriebenen Fastqs und Bams wieder löschen
        system("rm res.bam")
        system("rm rand.fastq")
        system("rm rand2.fastq")
        system("rm index*")
      }
    }

  }
  return(TRUE)
}


mapAllToOne <- function(fasta,pov){

  DS = dir("~/RealReadSimDS")

  for(i in 1:length(DS)){
    if(DS[i] != pov){
      folder = dir(paste("~/RealReadSImDS/",DS[i]))
      fastqs = folder[grep(".*.fastq$")]

      if(length(fastqs > 2 || length(fastqs) < 1)){
        print(paste("Die Anzahl der Fastq Dateien in ",DS[i]," stimmt nicht"))
        return(FALSE)
      }


    }
  }
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
