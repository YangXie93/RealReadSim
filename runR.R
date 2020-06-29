library(RealReadSim)
medium_prepped_path_vec = c("~/uni/CAMI-Daten/medium/prepFq1_270/","~/uni/CAMI-Daten/medium/prepFq2_270/")
medium_soure_genome_path = c("~/uni/CAMI-Daten/medium/source_genomes_medium/source_genomes/","~/uni/CAMI-Daten/medium/source_genomes_medium/source_genomes/")

medium_rrs_in_file = buildInputCsv(medium_soure_genome_path,medium_prepped_path_vec,mode = "fastq",out = "~/uni/CAMI-Daten/medium/mediumInFIle")

medium_cami = realReadSim(medium_rrs_in_file,"~/Cami_medium_DS",readAsBams = FALSE,outputFile = "~/CamiMediumOut",threads = 3,minMapq = 0)
