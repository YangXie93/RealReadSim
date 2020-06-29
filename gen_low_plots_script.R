install.packages("devtools")
install.packages("BiocManager")
BiocManager::install("Rsamtools")
devtools::install_github("https://github.com/YangXie93/RealReadSim.git")
install.packages("data.table")
install.packages("seqinr")
install.packages("pryr")

library("Rsamtools")
library("data.table")
library("seqinr")
library("pryr")
library("stringr")


library(RealReadSim)


#------------------------------ CAMI low data set ------------------------------------------------------------

path_to_low_fq = "/home/yang/uni/CAMI-Daten/low/RL_S001__insert_270.fq"
path_to_low_binning = "/home/yang/uni/CAMI-Daten/low/gs_read_mapping.binning"
path_to_low_prepped = "/home/yang/uni/CAMI-Daten/low/preparedFqs/"
path_to_low_source_genomes = "/home/yang/uni/CAMI-Daten/low/source_genomes_low/source_genomes/"

Rcpp::sourceCpp('~/prepareCamiDataR.cpp')

prepare(c(path_to_low_fq,path_to_low_binning,path_to_low_prepped))

low_rrs_in_file = buildInputCsv(path_to_low_source_genomes,path_to_low_prepped,mode = "fastq",out = "/home/yang/uni/CAMI-Daten/low/lowInFile")


low_cami_mapq40 = realReadSim("~/uni/CAMI-Daten/low/testRRS.txt",readAsBams = FALSE,metagenomeDir = "~/Cami_low_DS",outputFile = "~/camiLowOut",threads = 3,minMapq = 40)

low_cami_mapq30 = realReadSim("~/uni/CAMI-Daten/low/testRRS.txt",readAsBams = FALSE,metagenomeDir = "~/Cami_low_DS",outputFile = "~/camiLowOut",threads = 3,minMapq = 30)

low_cami_mapq20 = realReadSim("~/uni/CAMI-Daten/low/testRRS.txt",readAsBams = FALSE,metagenomeDir = "~/Cami_low_DS",outputFile = "~/camiLowOut",threads = 3,minMapq = 20)

low_cami_mapq10 = realReadSim("~/uni/CAMI-Daten/low/testRRS.txt",readAsBams = FALSE,metagenomeDir = "~/Cami_low_DS",outputFile = "~/camiLowOut",threads = 3,minMapq = 10)

low_cami_mapq0 = realReadSim("~/uni/CAMI-Daten/low/testRRS.txt",readAsBams = FALSE,metagenomeDir = "~/Cami_low_DS",outputFile = "~/camiLowOut",threads = 3,minMapq = 0)


mapq_cumu_length_low = list()
mapq_nr_length_low = list()

ordered_lengths_mapq40 = sort(low_cami_mapq40$length,decreasing = TRUE)
cumulative_length_mapq40 = c()
for(i in 1:length(ordered_lengths_mapq40)){
    cumulative_length_mapq40[i] = sum(ordered_lengths_mapq40[1:i])
}
nr_of_conts_low_mapq40 = c(1:length(cumulative_length_mapq40))

mapq_cumu_length_low[[1]] = cumulative_length_mapq40
mapq_nr_length_low[[1]] = nr_of_conts_low_mapq40

ordered_lengths_mapq30 = sort(low_cami_mapq30$length,decreasing = TRUE)
cumulative_length_mapq30 = c()
for(i in 1:length(ordered_lengths_mapq30)){
    cumulative_length_mapq30[i] = sum(ordered_lengths_mapq30[1:i])
}
nr_of_conts_low_mapq30 = c(1:length(cumulative_length_mapq30))

mapq_cumu_length_low[[2]] = cumulative_length_mapq30
mapq_nr_length_low[[2]] = nr_of_conts_low_mapq30

ordered_lengths_mapq20 = sort(low_cami_mapq20$length,decreasing = TRUE)
cumulative_length_mapq20 = c()
for(i in 1:length(ordered_lengths_mapq20)){
    cumulative_length_mapq20[i] = sum(ordered_lengths_mapq20[1:i])
}
nr_of_conts_low_mapq20 = c(1:length(cumulative_length_mapq20))

mapq_cumu_length_low[[3]] = cumulative_length_mapq20
mapq_nr_length_low[[3]] = nr_of_conts_low_mapq20

ordered_lengths_mapq10 = sort(low_cami_mapq10$length,decreasing = TRUE)
cumulative_length_mapq10 = c()
for(i in 1:length(ordered_lengths_mapq10)){
    cumulative_length_mapq10[i] = sum(ordered_lengths_mapq10[1:i])
}
nr_of_conts_low_mapq10 = c(1:length(cumulative_length_mapq10))

mapq_cumu_length_low[[4]] = cumulative_length_mapq10
mapq_nr_length_low[[4]] = nr_of_conts_low_mapq10

ordered_lengths_mapq0 = sort(low_cami_mapq0$length,decreasing = TRUE)
cumulative_length_mapq0 = c()
for(i in 1:length(ordered_lengths_mapq0)){
    cumulative_length_mapq0[i] = sum(ordered_lengths_mapq0[1:i])
}
nr_of_conts_low_mapq0 = c(1:length(cumulative_length_mapq0))

mapq_cumu_length_low[[5]] = cumulative_length_mapq0
mapq_nr_length_low[[5]] = nr_of_conts_low_mapq0


nms = c("minimum mapq 40","minimum mapq 30","minimum mapq 20","minimum mapq 10","minimum mapq 0")
clrs = c("red","blue","green","yellow","navy")

plot.default(mapq_nr_length_low,mapq_cumu_length_low,type = "l",col = clrs )
abline(h=max(unlist(mapq_cumu_length_low))/2,clrs)
legend("topright", border="black",
       legend=nms,
       col=clrs,
       cex = 0.65,
       pch=c(19,19,19,19))
