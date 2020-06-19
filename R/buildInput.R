buildInputCsv <- function(fastaPath,readPath,mode = "fastq",out = "~"){
    print(length(fastaPath))
    print(length(readPath))
    if(length(fastaPath) != length(readPath)){
        print("fastaPath and readPath have to have the same length")
        return(0)
    }

    res = c()

    for(j in 1:length(fastaPath)){
        if(substr(fastaPath[j],nchar(fastaPath[j]),nchar(fastaPath[j])) != "/"){
            fastaPath[j] = paste0(fastaPath[j],"/")
        }

        if(substr(readPath[j],nchar(readPath[j]),nchar(readPath[j])) != "/"){
            readPath[j] = paste0(readPath[j],"/")
        }

        fastas = dir(fastaPath[j])
        reads = dir(readPath[j])
        if(mode == "fastq"){
            namesReads = sub(".fastq|.fq|_[1|2]_.*","",reads)
        }
        if(mode == "bam"){
            namesReads = sub(".bam","",reads)
        }
        namesFasta = sub(".fna|.fasta|.final.*|.gt1kb.*|.scaffolds|_run.*|.*.fai","",fastas)

        print(namesFasta[!(namesFasta %in% namesReads)])#######################

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
                rLink2[i] = paste0(readPath[j],rLink2[i])
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
        write.csv(table,paste0(out,j,".txt"),row.names = FALSE,quote = FALSE)
        res[j] = paste0(out,j,".txt")
    }
    return(res)
}
