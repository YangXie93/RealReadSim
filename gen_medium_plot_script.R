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


#------------------------------ CAMI medium data set -----------------------------------------------------


pah_to_medium_prep_in = "/home/yang/uni/CAMI-Daten/medium/prepInFile.txt"
# the prepIn file is a csv file in which in the first column the the fastq file is noted, in the second column the .binning file is noted
# and in the last column the output directory is noted (it is important to have different output directorys for each row and that these
# directorys already exist and are empty)


path_to_medium_prepped = "/home/yang/uni/CAMI-Daten/medium/preparedFqs/"
path_to_medium_source_genomes = "/home/yang/uni/CAMI-Daten/medium/source_genomes_low/source_genomes/"

Rcpp::sourceCpp('~/prepareCamiDataR.cpp')

prepareCamiData(path_to_medium_prep_in)

# -------------------- combining the long and short read files --------------------------------------------
#
path1 = "~/uni/CAMI-Daten/medium/prepFq1_5000/"
path2 = "~/uni/CAMI-Daten/medium/prepFq1_270/"

medium1Ins5000dir = dir("~/uni/CAMI-Daten/medium/prepFq1_5000/")
medium1Ins270dir = dir("~/uni/CAMI-Daten/medium/prepFq1_270/")

path1 = "~/uni/CAMI-Daten/medium/prepFq2_5000/"
path2 = "~/uni/CAMI-Daten/medium/prepFq2_270/"

medium1Ins5000dir = dir("~/uni/CAMI-Daten/medium/prepFq2_5000/")
medium1Ins270dir = dir("~/uni/CAMI-Daten/medium/prepFq2_270/")

for(i in 1:length(medium1Ins270dir)){
    x = which(medium1Ins270dir == medium1Ins5000dir[i])
    if(length(x) > 0)
        system(paste0("cat ",path1,medium1Ins5000dir[i]," >> ",path2,medium1Ins270dir[x]))
}


medium_prepped_path_vec = c("~/uni/CAMI-Daten/medium/prepFq1_270/","~/uni/CAMI-Daten/medium/prepFq2_270/")
medium_soure_genome_path = c("~/uni/CAMI-Daten/medium/source_genomes_medium/source_genomes/","~/uni/CAMI-Daten/medium/source_genomes_medium/source_genomes/")

medium_rrs_in_file = buildInputCsv(medium_soure_genome_path,medium_prepped_path_vec,mode = "fastq",out = "~/uni/CAMI-Daten/medium/mediumInFIle")

medium_cami = realReadSim(medium_rrs_in_file,"~/Cami_medium_DS",readAsBams = FALSE,outputFile = "~/CamiMediumOut",threads = 3,minMapq = 0)
