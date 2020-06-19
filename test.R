path = "/home/yang/CAMI-High-DS/"
dir = dir(path)
dir = dir[which(dir != "Crossmaps")]
dir = dir[which(dir != "DSTable.Rds")]

paste0(path,dir,"/bin.fasta",collapse = ",")

for(i in 1:length(dir)){
    tmpDir = paste0(path,dir[which(dir != dir[i])],"/bin.fasta",collapse = ",")
    system(paste0("bowtie2-build --threads 3 ",tmpDir," ",path,dir[i],"/index"),ignore.stdout = TRUE)
}



dir = dir(path)
dir = dir[which(dir != "Crossmaps")]
dir = dir[which(dir != "DSTable.Rds")]

for(i in 1:length(dir)){
    system(paste0("rm ",path,dir[i],"/index*"))
}
