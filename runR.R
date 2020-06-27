library(RealReadSim)
path_to_low_fq = "/home/yang/uni/CAMI-Daten/low/RL_S001__insert_270.fq"
path_to_low_binning = "/home/yang/uni/CAMI-Daten/low/gs_read_mapping.binning"
path_to_low_prepped = "/home/yang/uni/CAMI-Daten/low/preparedFqs/"
path_to_low_source_genomes = "/home/yang/uni/CAMI-Daten/low/source_genomes_low/source_genomes/"
low_rrs_in_file = buildInputCsv(path_to_low_source_genomes,path_to_low_prepped,mode = "fastq",out = "/home/yang/uni/CAMI-Daten/low/lowInFile")
low_cami_mapq0 = realReadSim("~/uni/CAMI-Daten/low/testRRS.txt",readAsBams = FALSE,metagenomeDir = "~/Cami_low_DS",outputFile = "~/camiLowOut",threads = 3,minMapq = 0)
