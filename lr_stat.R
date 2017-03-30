
library(seqinr)

# functions
#==========

fasta_median <- function(read_object){
  result <- list()
  for( i in c(1:length(read_object)))
  {
    result[[i]] <- length(read_object[[i]])
  } 
  return(result)
}


#==================================
#==================================

# reading uncorrected FASTA data

ecoli_fasta <- read.fasta("/data/test_data_from_server_maize/proovread/ecoli/ecoli_pacbio.fa_0001.fasta") 
trypanosoma_fasta <- read.fasta("/data/test_data_from_server_maize/pbcr/trypanosoma/trypanosoma11_pacbio.fa_0001.fasta")
yeast_fasta <- read.fasta("/data/test_data_from_server_maize/pbcr/yeast/4/yeast4_pacbio.fa_0001.fasta")
rice_fasta <- read.fasta("/data/test_data_from_server_maize/pbcr/rice/rice_pacbio.fa_0001.fasta")
human_fasta <- read.fasta("/data/test_data_from_server_maize/pbcr/human21/human21_pacbio.fa_0001.fasta")

#==================================

# read corrected fasta data
#==========================

# 1- full corrected reads
#=======================


# 1- lsc
#=======

lsc_ecoli_corrected_full_length <- read.fasta("/data/test_data_from_server_maize/LSC/ecoli/lsc_ecoli_2016-05-02/output/full_LR.fa")
lsc_trypanosoma_corrected_full_length <- read.fasta("/data/test_data_from_server_maize/LSC/trypanosoma/lsc_trypanosoma_2016-05-02/output/full_LR.fa")
lsc_yeast_corrected_full_length <- read.fasta("/data/test_data_from_server_maize/LSC/yeast/4/lsc_yeast4_2016-05-02/output/full_LR.fa")
lsc_rice_corrected_full_length <- read.fasta("/data/test_data_from_server_maize/LSC/rice/output/full_LR.fa")
lsc_human_corrected_full_length <- read.fasta("/data/test_data_from_server_maize/LSC/human/lsc_human_2016-04-29/output/full_LR.fa")

# 2- lordec
#========

lordec_ecoli_corrected_full_length <- read.fasta("/data/test_data_from_server_maize/LoRDEC/ecoli/ecoli_correcr_LR.fa")
lordec_trypanosoma_corrected_full_length <- read.fasta("/data/test_data_from_server_maize/LoRDEC/trypnosoma/trypanosoma_correcr_LR.fa")
lordec_yeast_corrected_full_length <- read.fasta("/data/test_data_from_server_maize/LoRDEC/yeast/4/yeast4_correcr_LR.fa")
lordec_rice_corrected_full_length <- read.fasta("/data/test_data_from_server_maize/LoRDEC/rice/rice_correcr_LR.fa")
lordec_human_corrected_full_length <- read.fasta("/data/test_data_from_server_maize/LoRDEC/human21/human_correcr_LR.fa")

# 3- proovread
#=============

proovread_ecoli_corrected_full_length <- read.fasta("/data/test_data_from_server_maize/proovread/ecoli/correct_2016-08-31/ecoli/ecoli.untrimmed.fa")
proovread_trypanosoma_corrected_full_length <- read.fasta("/data/test_data_from_server_maize/proovread/trypanosoma/trypanosoma/correct_2016-09-01/trypanosoma/trypanosoma.untrimmed.fa")
proovread_yeast_corrected_full_length <- read.fasta("/data/test_data_from_server_maize/proovread/yeast/4/yeast4/correct_2016-09-01/yeast/yeast.untrimmed.fa")
proovread_rice_corrected_full_length <- read.fasta("/data/test_data_from_server_maize/proovread/rice/correct_2016-09-06/rice/rice.untrimmed.fa")
proovread_human_corrected_full_length <- read.fasta("/data/test_data_from_server_maize/proovread/human21/correct_2016-09-10/human/human.untrimmed.fa")



# 2- split and trimmed corrected reads
#====================================


# 1- lordec
#==========
lordec_ecoli_corrected_split_length <- read.fasta("/data/test_data_from_server_maize/LoRDEC/ecoli/ecoli_correcr_split_trim_LR.fa")
lordec_trypanosoma_corrected_split_length <- read.fasta("/data/test_data_from_server_maize/LoRDEC/trypnosoma/trypanosoma_correcr__trim_split_LR.fa")
lordec_yeast_corrected_split_length <- read.fasta("/data/test_data_from_server_maize/LoRDEC/yeast/4/yeast4_correcr_trim_split_LR.fa")
lordec_rice_corrected_split_length <- read.fasta("/data/test_data_from_server_maize/LoRDEC/rice/rice_correcr_trim_split_LR.fa")
lordec_human_corrected_split_length <- read.fasta("/data/test_data_from_server_maize/LoRDEC/human21/human_correcr_LR_trim_split.fa")


# 2- proovread
#=============
proovread_ecoli_corrected_split_length <- read.fasta("/data/test_data_from_server_maize/proovread/ecoli/correct_2016-08-31/ecoli/ecoli.trimmed.fa")
proovread_trypanosoma_corrected_split_length <- read.fasta("/data/test_data_from_server_maize/proovread/trypanosoma/trypanosoma/correct_2016-09-01/trypanosoma/trypanosoma.trimmed.fa")
proovread_yeast_corrected_split_length <- read.fasta("/data/test_data_from_server_maize/proovread/yeast/4/yeast4/correct_2016-09-01/yeast/yeast.trimmed.fa")
proovread_rice_corrected_split_length <- read.fasta("//data/test_data_from_server_maize/proovread/rice/correct_2016-09-06/rice/rice.trimmed.fa")
proovread_human_corrected_split_length <- read.fasta("/data/test_data_from_server_maize/proovread/human21/correct_2016-09-10/human/human.trimmed.fa")


# 3- pbcr
#========




# results
#========


# full read stat before correction
#================================

ecoli_fasta_list <- fasta_median(ecoli_fasta)
paste("E.coli Number of reads is: ", length(ecoli_fasta))
print("Read length summary")
summary(as.numeric(fasta_median(ecoli_fasta)))


paste("Trypanosoma Number of reads is: ", length(trypanosoma_fasta))
print("Read length summary")
summary(as.numeric(fasta_median(trypanosoma_fasta)))


paste("Yeast Number of reads is: ", length(yeast_fasta))
print("Read length summary")
summary(as.numeric(fasta_median(yeast_fasta)))

paste("Rice Number of reads is: ", length(rice_fasta))
print("Read length summary")
summary(as.numeric(fasta_median(rice_fasta)))

paste("Human Number of reads is: ", length(human_fasta))
print("Read length summary")
summary(as.numeric(fasta_median(human_fasta)))

# full read stat after correction
#================================

# 1- lsc
#=======

paste("LSC E.coli Number of reads after correction is: ", length(lsc_ecoli_corrected_full_length))
print("Read length summary")
summary(as.numeric(fasta_median(lsc_ecoli_corrected_full_length)))

paste("LSC trypanosoma Number of reads after correction is: ", length(lsc_trypanosoma_corrected_full_length))
print("Read length summary")
summary(as.numeric(fasta_median(lsc_trypanosoma_corrected_full_length)))


paste("LSC yeast Number of reads after correction is: ", length(lsc_yeast_corrected_full_length))
print("Read length summary")
summary(as.numeric(fasta_median(lsc_yeast_corrected_full_length)))


paste("LSC Rice Number of reads after correction is: ", length(lsc_rice_corrected_full_length))
print("Read length summary")
summary(as.numeric(fasta_median(lsc_rice_corrected_full_length)))

paste("LSC Human Number of reads after correction is: ", length(lsc_human_corrected_full_length))
print("Read length summary")
summary(as.numeric(fasta_median(lsc_human_corrected_full_length)))


# 2- proovread
#=============

proovread_ecoli_list_read_length <- fasta_median(proovread_ecoli_corrected_full_length)
paste("proovread E.coli Number of reads after correction is: ", length(proovread_ecoli_list_read_length))
print("Read length summary")
summary(as.numeric(proovread_ecoli_list_read_length))
rm(proovread_ecoli_corrected_full_length)


proovread_trypanosoma_list_read_length <- fasta_median(proovread_trypanosoma_corrected_full_length)
paste("proovread trypanosoma Number of reads after correction is: ", length(proovread_trypanosoma_list_read_length))
print("Read length summary")
summary(as.numeric(proovread_trypanosoma_list_read_length))
rm(proovread_trypanosoma_corrected_full_length)



proovread_yeast_list_read_length <- fasta_median(proovread_yeast_corrected_full_length)
paste("proovread yeast Number of reads after correction is: ", length(proovread_yeast_list_read_length))
print("Read length summary")
summary(as.numeric(proovread_yeast_list_read_length))
rm(proovread_yeast_corrected_full_length)



proovread_rice_list_read_length <- fasta_median(proovread_rice_corrected_full_length)
paste("proovread rice Number of reads after correction is: ", length(proovread_rice_list_read_length))
print("Read length summary")
summary(as.numeric(proovread_rice_list_read_length))
rm(proovread_rice_corrected_full_length)


proovread_human_list_read_length <- fasta_median(proovread_human_corrected_full_length)
paste("proovread human Number of reads after correction is: ", length(proovread_human_list_read_length))
print("Read length summary")
summary(as.numeric(proovread_human_list_read_length))
rm(proovread_human_corrected_full_length)

# 3- lordec
#==========

lordec_ecoli_list_read_length <- fasta_median(lordec_ecoli_corrected_full_length)
paste("lordec E.coli Number of reads after correction is: ", length(lordec_ecoli_list_read_length))
print("Read length summary")
summary(as.numeric(lordec_ecoli_list_read_length))
rm(lordec_ecoli_corrected_full_length)


lordec_trypanosoma_list_read_length <- fasta_median(lordec_trypanosoma_corrected_full_length)
paste("lordec trypanosoma Number of reads after correction is: ", length(lordec_trypanosoma_list_read_length))
print("Read length summary")
summary(as.numeric(lordec_trypanosoma_list_read_length))
rm(lordec_trypanosoma_corrected_full_length)


lordec_yeast_list_read_length <- fasta_median(lordec_yeast_corrected_full_length)
paste("lordec Yeast Number of reads after correction is: ", length(lordec_yeast_list_read_length))
print("Read length summary")
summary(as.numeric(lordec_yeast_list_read_length))
rm(lordec_yeast_corrected_full_length)


lordec_rice_list_read_length <- fasta_median(lordec_rice_corrected_full_length)
paste("lordec rice Number of reads after correction is: ", length(lordec_rice_list_read_length))
print("Read length summary")
summary(as.numeric(lordec_rice_list_read_length))
rm(lordec_rice_corrected_full_length)


lordec_human_list_read_length <- fasta_median(lordec_human_corrected_full_length)
paste("lordec human Number of reads after correction is: ", length(lordec_human_list_read_length))
print("Read length summary")
summary(as.numeric(lordec_human_list_read_length))
rm(lordec_human_corrected_full_length)


# split trimmed reads stat after correction
#==========================================

# 1- lordec
#==========

lordec_ecoli_list_read_split_length <- fasta_median(lordec_ecoli_corrected_split_length)
paste("lordec E.coli Number of reads after correction is: ", length(lordec_ecoli_list_read_split_length))
print("Read length summary")
summary(as.numeric(lordec_ecoli_list_read_split_length))
rm(lordec_ecoli_corrected_split_length)


lordec_yeast_list_read_split_length <- fasta_median(lordec_yeast_corrected_split_length)
paste("lordec yeast Number of reads after correction is: ", length(lordec_yeast_list_read_split_length))
print("Read length summary")
summary(as.numeric(lordec_yeast_list_read_split_length))
rm(lordec_yeast_corrected_split_length)

lordec_trypanosoma_list_read_split_length <- fasta_median(lordec_trypanosoma_corrected_split_length)
paste("lordec trypanosoma Number of reads after correction is: ", length(lordec_trypanosoma_list_read_split_length))
print("Read length summary")/data/test_data_from_server_maize/proovread/ecoli/correct_2016-08-31/ecoli
summary(as.numeric(lordec_trypanosoma_list_read_split_length))
rm(lordec_trypanosoma_corrected_split_length)

lordec_rice_list_read_split_length <- fasta_median(lordec_rice_corrected_split_length)
paste("lordec rice Number of reads after correction is: ", length(lordec_rice_list_read_split_length))
print("Read length summary")
summary(as.numeric(lordec_rice_list_read_split_length))
rm(lordec_rice_corrected_split_length)


lordec_human_list_read_split_length <- fasta_median(lordec_human_corrected_split_length)
paste("lordec human Number of reads after correction is: ", length(lordec_human_list_read_split_length))
print("Read length summary")
summary(as.numeric(lordec_human_list_read_split_length))
rm(lordec_human_corrected_split_length)

# 1- proovread
#=============

proovread_ecoli_list_read_split_length <- fasta_median(proovread_ecoli_corrected_split_length)
paste("proovread E.coli Number of reads after correction is: ", length(proovread_ecoli_list_read_split_length))
print("Read length summary")
summary(as.numeric(proovread_ecoli_list_read_split_length))
rm(proovread_ecoli_corrected_split_length)


proovread_trypanosoma_list_read_split_length <- fasta_median(proovread_trypanosoma_corrected_split_length)
paste("proovread E.coli Number of reads after correction is: ", length(proovread_trypanosoma_list_read_split_length))
print("Read length summary")
summary(as.numeric(proovread_trypanosoma_list_read_split_length))
rm(proovread_trypanosoma_corrected_split_length)


proovread_yeast_list_read_split_length <- fasta_median(proovread_yeast_corrected_split_length)
paste("proovread E.coli Number of reads after correction is: ", length(proovread_yeast_list_read_split_length))
print("Read length summary")
summary(as.numeric(proovread_yeast_list_read_split_length))
rm(proovread_yeast_corrected_split_length)

proovread_rice_list_read_split_length <- fasta_median(proovread_rice_corrected_split_length)
paste("proovread E.coli Number of reads after correction is: ", length(proovread_rice_list_read_split_length))
print("Read length summary")
summary(as.numeric(proovread_rice_list_read_split_length))
rm(proovread_rice_corrected_split_length)

proovread_human_list_read_split_length <- fasta_median(proovread_human_corrected_split_length)
paste("proovread human Number of reads after correction is: ", length(proovread_human_list_read_split_length))
print("Read length summary")
summary(as.numeric(proovread_human_list_read_split_length))
rm(proovread_human_corrected_split_length)


# rest I got it from assess_stat.R as the dataframe there actuallly contains all the info  ex: summary(pbcr_ecoli_local$subLength)

