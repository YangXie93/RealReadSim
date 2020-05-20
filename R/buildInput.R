buildInputCsv <- function(fastaPath,readPath,mode = "fastq",out = "~"){

    if(substr(fastaPath,nchar(fastaPath),nchar(fastaPath)) != "/"){
        fastaPath = paste0(fastaPath,"/")
    }

    if(substr(readPath,nchar(readPath),nchar(readPath)) != "/"){
        readPath = paste0(readPath,"/")
    }

    fastas = dir(fastaPath)
    reads = dir(readPath)
    if(mode == "fastq"){
        namesReads = sub(".fq|_[1|2]_.*","",reads)
    }
    if(mode == "bam"){
        namesReads = sub(".fq|_[1|2]_.*","",reads)
    }
    namesFasta = sub(".fna|.fasta|.final.*|.gt1kb.*|.scaffolds|_run.*|.*.fai","",fastas)

    print(namesFasta[!(namesFasta %in% namesReads)])

    fastas = fastas[namesFasta %in% namesReads]
    namesFasta = namesFasta[namesFasta %in% namesReads]


    rLink1 = c()
    rLink2 = c()
    for(i in 1:length(namesFasta)){
        rLink1[i] = reads[which(namesReads == namesFasta[i])[1]]
        rLink2[i] = reads[which(namesReads == namesFasta[i])[2]]
    }
    for(i in 1:length(rLink2)){
        if(!is.na(rLink2[i])){
            rLink2[i] = paste0(readPath,rLink2[i])
        }
        else{
            rLink2[i] = ""
        }
    }
    print(length(namesFasta))
    print(length(rLink1))
    print(length(rLink2))
    if(paste0(rLink2,collapse = "") == ""){
        table = data.frame(fasta = paste0(fastaPath,fastas),readsFwd = paste0(readPath,rLink1))
    }
    else{
        table = data.frame(fasta = paste0(fastaPath,fastas),readsFwd = paste0(readPath,rLink1),readsRev = rLink2)
    }
    write.csv(table,out,row.names = FALSE,quote = FALSE)
}
