
#Diese Funktion macht den Pre-build fuer bowtie2 aus R herraus
prebuildBowtie2 <- function(sequence_fasta, build_name){
  print(paste("bowtie2-build",sequence_fasta,build_name,sep = " "))
  system(paste("bowtie2-build",sequence_fasta,build_name,sep = " "))
  
  return(build_name)
}

#Diese Funktion dient der Benutzung von Bowtie2 aus R herraus 
useBowtie2 <- function(build_name, reads_fastq, reads2_fastq, output, paired_end = FALSE, fasta_reads = FALSE,defaultOptions = TRUE){
  if(paired_end){
    reads = paste("-1",reads_fastq,"-2",reads2_fastq,sep = " ")
  }
  else{
    if(fasta_reads){
      reads = paste("-f",reads_fastq,sep = " ")
    }
    else{
      reads = paste("-U",reads_fastq,sep = " ")  
    }
    
  }
  if(defaultOptions){
    options = "--no-unal"
  }
  else{
    print("wished Bowtie2 options (only the options section without the command or input/output):")
    options = readline()
  }
  output = paste(build_name,".sam",sep = "")
  print(paste("bowtie2",options ,"-x",build_name,reads,"-S",output,sep = " "))
  system(paste("bowtie2",options ,"-x",build_name,reads,"-S",output,sep = " "))
  
  return(output)
}


#In dieser Funktion werden die oberen beiden funktional in Reihegeschaltet
mkSams <- function(filenames, prebuild = TRUE, i = 1,defaultOptions = TRUE){
    
  if(prebuild){
    if(filenames[i,4] == "" || is.na(filenames[i,4])){
      build_name = prebuildBowtie2(filenames[i,1], paste("Bowtie_build",i,sep = ""))
    }
    else{
      build_name = prebuildBowtie2(filenames[i,1], filenames[i,4])
    }
  }
    
  if(filenames[i,3] == "" || is.na(filenames[i,3])){
    sam = useBowtie2(build_name,filenames[i,2],defaultOptions = defaultOptions)
  }
  else{
    sam = useBowtie2(build_name,filenames[i,2],filenames[i,3],paired_end = TRUE ,defaultOptions = defaultOptions)
  }
    
  return(sam)
}


bamOrBowtie <- function(filenames,readAsBams,i,params,defaultOptions = TRUE){
  
  
  if(!readAsBams){
    sam = mkSams(filenames,prebuild,i,defaultOptions)   
    asBam(sam)                                
    if(i > 1){
      rm(temp)
    }
    data = scanBam(paste(substr(sam,1,nchar(sam) -4), ".bam",sep = ""),param = params)
  }
  else{
    print(filenames[i,2])
    data = scanBam(filenames[i,2],param = params)
  }
  data = data.frame(pos = data[[1]]$pos,width = data[[1]]$qwidth,DNAString = data[[1]]$rname)
  data = data[!is.na(data$pos),]
  return(data)
}
