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


#------------------------------ CAMI high data set -------------------------------------------------------


prepareCamiData("/home/yang/uni/CAMI-Daten/high/prepInFile2.txt")

prepareCamiData("/home/yang/uni/CAMI-Daten/high/prepInFile.txt")

path_to_high_fqs = c("/home/yang/uni/CAMI-Daten/high/prepFqs1/","/home/yang/uni/CAMI-Daten/high/prepFqs2/",
                     "/home/yang/uni/CAMI-Daten/high/prepFqs3/","/home/yang/uni/CAMI-Daten/high/prepFqs4/",
                     "/home/yang/uni/CAMI-Daten/high/prepFqs5/")
path_to_high_source_genomes = c("/home/yang/uni/CAMI-Daten/high/source_genomes_high/source_genomes/",
                                "/home/yang/uni/CAMI-Daten/high/source_genomes_high/source_genomes/",
                                "/home/yang/uni/CAMI-Daten/high/source_genomes_high/source_genomes/",
                                "/home/yang/uni/CAMI-Daten/high/source_genomes_high/source_genomes/",
                                "/home/yang/uni/CAMI-Daten/high/source_genomes_high/source_genomes/")

high_rrs_in_files = buildInputCsv(path_to_high_source_genomes,path_to_high_fqs,mode = "fastq",
                                  out = "/home/yang/uni/CAMI-Daten/high/highInFile")

addToDataSystem(filenames_csv = high_rrs_in_files,metagenomeDir = "~/Cami-High-DS",bowtieOptions = "",
                minIdL = 2000,readAsBams = FALSE,threads = 3,bowtiebuildOptions = "")


metagenomeDir = "~/Cami-High-DS"
bowtiebuildOptions = ""

crossMapRRSDS(2000,3,"~/Cami-High-DS/")
