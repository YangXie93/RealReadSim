buildInputCsv <- function(fastaPath,readPath,mode = "fastq",out = "~"){

    fastas = dir(fastaPath)
    reads = dir(readPath)
    if(mode == "fastq"){
        namesReads = sub(".fq|_[1|2]_.*","",reads)
    }
    if(mode == "bam"){
        namesReads = sub(".fq|_[1|2]_.*","",reads)
    }
    namesFasta = sub(".fna|.fasta|.final.*|.gt1kb.*|.scaffolds|_run.*|.*.fai","",fastas)

    link = c()
    j = 1
    for(i in 1:length(namesReads)){
        n = which(namesFasta == namesReads[i])
        print(namesReads[i])
        if(length(n) > 0){
            link[j] = n
            j = j +1
        }
    }
    namesFasta = namesFasta[link]
    fastas = fastas[link]

    rLink1 = c()
    rLink2 = c()
    print(namesFasta)
    for(i in 1:length(namesFasta)){
        rLink1[i] = reads[which(namesReads == namesFasta[i])][1]
        rLink2[i] = reads[which(namesReads == namesFasta[i])][2]
    }
    for(i in 1:length(rLink2)){
        if(!is.na(rLink2[i])){
            rLink2[i] = paste0(readPath,rLink2)
        }
        else{
            rLink2[i] = ""
        }
    }
    if(paste0(rLink2,collapse = "") == ""){
        table = data.frame(fasta = paste0(fastaPath,fastas),readsFwd = paste0(readPath,rLink1))
    }
    else{
        table = data.frame(fasta = paste0(fastaPath,fastas),readsFwd = paste0(readPath,rLink1),readsRev = rLink2)
    }
    write.csv(table,out,row.names = FALSE)
}
