 # Rcpp::sourceCpp('src/coAss.cpp')
#
# translateOverlap(7,10,3,10,6,11,4,9) # Lösung: (0,3,2,5,1)
#
# translateOverlap(7,11,12,19,5,10,13,18) # Lösung: (0,3,3,6,2)
#
# translateOverlap(4,9,7,13,2,7,9,14) # Lösung: (0,2,4,6,1)
#
# translateOverlap(3,11,9,19,5,10,12,17) #Lösung: (2,7,3,8,3)
#
# translateOverlap(6,10,11,16,5,11,10,16) #Lösung: (0,4,0,4,0)

#
# s1 = c(12,25)
# e1 = c(23,34)
# s2 = c(22,38)
# e2 = c(34,49)
# a1s = c(12,21,30)
# a1e = c(17,27,35)
# a2s = c(24,33,42)
# a2e = c(29,39,47)
# cov1 = list(c(rep(10,6),rep(11,6)),rep(5,10))
# cov2 = list(rep(1,13),rep(2,12))
# sq1 = c(paste0(rep("A",12),collapse = ""),paste0(rep("A",10),collapse = ""))
# sq2 = c(paste0(rep("T",13),collapse = ""),paste0(rep("T",12),collapse = ""))
# nm1 = c("a","a")
# nm2 = c("b","b")

 # s1 = c(3,21)
 # e1 = c(19,31)
 # s2 = c(9,18,31)
 # e2 = c(14,28,38)
 # a1s = c(4,13,22)
 # a1e = c(11,18,33)
 # a2s = c(8,17,26)
 # a2e = c(15,22,37)
 # cov1 = list(rep(10,17),rep(1,11))
 # cov2 = list(rep(2,6),rep(2,11),rep(11,8))
 # sq1 = c(paste0(rep("A",17),collapse = ""),paste0(rep("A",11),collapse = ""))
 # sq2 = c(paste0(rep("T",6),collapse = ""),paste0(rep("T",11),collapse = ""),paste0(rep("T",8),collapse = ""))
 # nm1 = c("a","a")
 # nm2 = c("b","b","b")

#  s1 = c(7,16,27)
#  e1 = c(13,24,33)
#  s2 = c(4,13)
#  e2 = c(11,25)
#  a1s = c(5,20)
#  a1e = c(17,28)
#  a2s = c(1,16)
#  a2e = c(13,24)
#  cov1 = list(rep(10,7),rep(1,9),rep(13,7))
#  cov2 = list(rep(2,8),rep(12,13))
#  sq1 = c(paste0(rep("A",7),collapse = ""),paste0(rep("A",9),collapse = ""),paste0(rep("A",7),collapse = ""))
#  sq2 = c(paste0(rep("T",8),collapse = ""),paste0(rep("T",13),collapse = ""))
# nm1 = c("a","a","a")
# nm2 = c("b","b")
#
# mkChimeras(s1,e1,cov1,s2,e2,cov2,a1s,a1e,a2s,a2e,sq1,sq2,nm1,nm2)

# library(RealReadSim)
#x = realReadSim(c("/home/yang/uni/BA-Projekt-Data/testBinning.txt","/home/yang/uni/BA-Projekt-Data/testBinning.txt"),takeAll = TRUE,metagenomeDir = "~/testMultiSampleInput")
#x = realReadSim(takeAll = TRUE,metagenomeDir = "~/testMetagenome",fileOutput = TRUE,outputFile = "~/out")


# metagenomeDir = "~/Cami_low_DS"
# bowtiebuildOptions = ""
#
# dr = dir(metagenomeDir)
# dr = dr[which(dr != "Crossmaps")]
# dr = dr[which(dr != "DSTable.Rds")]
# dr = dr[!str_detect(dr,"index.*")]
#
#
# x = realReadSim("~/uni/CAMI-Daten/low/testRRS.txt",readAsBams = FALSE,metagenomeDir = "~/Cami_low_DS",outputFile = "~/camiLowOut")
#
# low_table = readRDS("~/Cami_low_DS/DSTable.Rds")
# low_table$isCrossmapped = rep(FALSE,length(low_table$isCrossmapped))
# saveRDS(low_table,"~/Cami_low_DS/DSTable.Rds")
# crossMapRRSDS(minL = 2000,threads = 3,metagenomeDir = "~/Cami_low_DS/")
#
#
# buildInputCsv("~/uni/CAMI-Daten/medium/source_genomes_medium/source_genomes/",c("~/uni/CAMI-Daten/medium/"))
#
# addToDataSystem("~/uni/CAMI-Daten/medium/mediumInFile.txt",bowtieOptions = "",minIdL = 2000,metagenomeDir = "~/Cami_medium_DS",bowtiebuildOptions = "",threads = 3,readAsBams = FALSE)


#x = realReadSim(takeAll = TRUE,readAsBams = FALSE,metagenomeDir = "~/CamiDataDS",fileOutput = TRUE,outputFile = "~/out")
# n = 0
# for(i in 1:length(x$start)){
#         if(x$end[i]-x$start[i]+1 != length(x$coverage[[i]])){
#                 n = n +1
#                 print(paste("1",i,x$end[i]-x$start[i]+1,length(x$coverage[[i]])))
#         }
#         if(length(x$coverage[[i]]) != nchar(x$seq[i])){
#                 n = n +1
#                 print(paste("2",i,length(x$coverage[[i]]),nchar(x$seq[i])))
#         }
#         if(nchar(x$seq[i]) != x$end[i]-x$start[i]+1){
#                 n = n +1
#                 print(paste("3",i,nchar(x$seq[i]),x$end[i]-x$start[i]+1))
#         }
# }
#
# for(i in 1:length(x$start)){
#         print(x$coverage[[i]][length(x$coverage[[i]])])
# }
#
# library("Rsamtools")
# library("data.table")
# library("seqinr")
# library("pryr")
# library("stringr")


# library(RealReadSim)
# #crossMapRRSDS(2000,"~/CAMI-High-DS/")
# x = realReadSim(filenames_csv = c("~/highIn1.txt","~/highIn2.txt","~/highIn3.txt","~/highIn4.txt","~/highIn5.txt"),
#                 readAsBams = FALSE,metagenomeDir = "~/CAMI-High-DS",fileOutput = TRUE,
#                 outputFile = "~/camiHighout")



# test getIdenticalSeqs
# Rcpp::sourceCpp('src/getIdenticalSeqs.cpp')
# starts1 = c(31,33,34,35,36)
# ends1 = c(34,35,37,38,38)
# starts2 = c(8,7,8,11,10)
# ends2 = c(11,9,11,14,12)
#
# print(getIdenticalSeqs(starts1,ends1,starts2,ends2,"1","2"))
#
# crossMapRRSDS(1000,"~/CamiDataDS/")
# crossMapRRSDS(1000,"~/CAMI-High-DS/")
#
# path1 = "~/test1/"
# path2 = "~/test2/"
#
# medium1Ins5000dir = dir("~/test1/")
# medium1Ins270dir = dir("~/test2/")
#
# path1 = "~/uni/CAMI-Daten/medium/prepFq1_5000/"
# path2 = "~/uni/CAMI-Daten/medium/prepFq1_270/"
#
# medium1Ins5000dir = dir("~/uni/CAMI-Daten/medium/prepFq1_5000/")
# medium1Ins270dir = dir("~/uni/CAMI-Daten/medium/prepFq1_270/")
#
# path1 = "~/uni/CAMI-Daten/medium/prepFq2_5000/"
# path2 = "~/uni/CAMI-Daten/medium/prepFq2_270/"
#
# medium1Ins5000dir = dir("~/uni/CAMI-Daten/medium/prepFq2_5000/")
# medium1Ins270dir = dir("~/uni/CAMI-Daten/medium/prepFq2_270/")
#
# for(i in 1:length(medium1Ins270dir)){
#         x = which(medium1Ins270dir == medium1Ins5000dir[i])
#         if(length(x) > 0)
#                 system(paste0("cat ",path1,medium1Ins5000dir[i]," >> ",path2,medium1Ins270dir[x]))
# }

#getLast8("/home/yang/uni/CAMI-Daten/low/gs_read_mapping.binning")


Rcpp::sourceCpp('/home/yang/prepareCamiDataR.cpp')

prepareCamiData("/home/yang/uni/CAMI-Daten/high/prepInFile.txt")

#prepare(c("/home/yang/uni/CAMI-Daten/low/RL_S001__insert_270.fq","/home/yang/uni/CAMI-Daten/low/gs_read_mapping.binning","/home/yang/uni/CAMI-Daten/low/preparedFqs/"))

#99829122

# prepareCamiData("/home/yang/uni/CAMI-Daten/medium/prepInFile.txt")

# medium1Ins5000dir = dir("~/uni/CAMI-Daten/medium/prepFq2_5000/")
# medium1Ins270dir = dir("~/uni/CAMI-Daten/medium/prepFq2_270/")
#
# for(i in 1:length(medium1Ins270dir)){
#         x = which(medium1Ins270dir == medium1Ins5000dir[i])
#         system(paste0("cat ",path1,medium1Ins5000dir[i]," >> ",path2,medium1Ins270dir[x]))
# }

#prepareCamiData("/home/yang/uni/CAMI-Daten/high/prepInFile.txt")

#buildInputCsv(readPath =  "~/uni/CAMI-Daten/low/preparedFqs/",fastaPath = "~/uni/CAMI-Daten/low/source_genomes_low/source_genomes/",out = "~/lowInCsv")

#buildInputCsv(fastaPath = c("~/uni/CAMI-Daten/medium/source_genomes_medium/source_genomes/",
              #               "~/uni/CAMI-Daten/medium/source_genomes_medium/source_genomes/"),
              # readPath = c("~/uni/CAMI-Daten/medium/prepFqs1","~/uni/CAMI-Daten/medium/prepFqs2"))

#buildInputCsv(fastaPath = c("~/uni/CAMI-Daten/high/source_genomes_high/","~/uni/CAMI-Daten/high/source_genomes_high/",
              #               "~/uni/CAMI-Daten/high/source_genomes_high/","~/uni/CAMI-Daten/high/source_genomes_high/",
              #               "~/uni/CAMI-Daten/high/source_genomes_high/","~/uni/CAMI-Daten/high/source_genomes_high/"),
              # readPath = c("~/uni/CAMI-Daten/high/prepFqs1","~/uni/CAMI-Daten/high/prepFqs2","~/uni/CAMI-Daten/high/prepFqs3",
              #              "~/uni/CAMI-Daten/high/prepFqs4","~/uni/CAMI-Daten/high/prepFqs5"),
              # out = "~/highIncsv")
#
# namesLow = dir("~/uni/CAMI-Daten/low/preparedFqs")
# namesLowx = gsub("_1_.fq|_2_.fq","",namesLow)
# path = "~/uni/CAMI-Daten/low/"
#
# dotBinning = c()
# dotFq = c()
#
# for(i in (1:length(namesLow))){
#         print(i)
#         dotBinning[i] = system(paste0("cat ",path,"gs_read_mapping.binning | grep -c ",namesLowx[i]),intern = TRUE)
#         dotFq[i] = system(paste0("wc -l ",path,"preparedFqs/",namesLow[i]),intern = TRUE)
# }
#
# lineNrs = data.frame(binning = dotBinning,fastq = dotFq)
# print(data.frame(binning = dotBinning,fastq = (dotFq/4)))
# readNum = system(paste0("wc -l ",path,"/RL_S001__insert_270.fq"),intern = TRUE)
# readNum = str_split_fixed(readNum," ",n =2)[,1]
# readNum = as.integer(readNum)
# readNumTrue = readNum/4
#
# readNumBin = system(paste0("wc -l ",path,"/gs_read_mapping.binning"),intern = TRUE)
# readNumBin = str_split_fixed(readNumBin," ",n =2)[,1]
# readNumBin = as.integer(readNumBin)

#testReadHash("/home/yang/uni/CAMI-Daten/low/gs_read_mapping.binning.l8","/home/yang/uni/CAMI-Daten/low/gs_read_mapping.binning")
#

assignInNamespace("cedta.override",
                  c(data.table:::cedta.override,"RealReadSim"),
                  "data.table")

