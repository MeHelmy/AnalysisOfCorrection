#==================================================
# global constant
column_names <- c("ch", "begin", "end", "name", "chr_ref", "begin_ref", "enf_ref", "repeat", "score", "strand", "repeat_family", "intersect_length")


human_NN_reads <- read.delim("/data/test_data_from_server_maize/assessment_2016-08-29/complementory_data/human21_pacbio_original_nucleotide_N_freq.bed")
rice_NN_reads <- read.delim("/data/test_data_from_server_maize/assessment_2016-08-29/complementory_data/rice_pacbio_original_nucleotide_N_freq.bed")
trypanosoma_NN_reads <- read.delim("/data/test_data_from_server_maize/assessment_2016-08-29/complementory_data/trypanosoma_pacbio_original_nucleotide_N_freq.bed")


#==================================================
# Functions

extract_rows <- function(whole_dataframe, small_dataframe, intersect_column_first, intersect_column_second, complement=TRUE){
  if(complement){
    return(subset(whole_dataframe, !(whole_dataframe[[intersect_column_first]] %in% small_dataframe[[intersect_column_second]])))
  }else{
    return(subset(whole_dataframe, (whole_dataframe[[intersect_column_first]] %in% small_dataframe[[intersect_column_second]])))
  }
  
}



analysis_repeat <- function(cor_repeat_data_frame, uncor_repeat_data_frame, columns_name, organism, NN_data){
  
  colnames(cor_repeat_data_frame) <- column_names
  colnames(uncor_repeat_data_frame) <- column_names
  
  if(!(missing(NN_data))){
  cor_repeat_data_frame <- separate(cor_repeat_data_frame, name, c('read', 'sigment'), "\\.", extra = "merge", remove = F, fill = "right")
  cor_repeat_data_frame <- extract_rows(whole_dataframe = cor_repeat_data_frame, small_dataframe = NN_data, intersect_column_first = "read", intersect_column_second = "name")
  uncor_repeat_data_frame <- separate(uncor_repeat_data_frame, name, c('read', 'sigment'), "\\.", extra = "merge", remove = F, fill = "right")
  uncor_repeat_data_frame <- extract_rows(whole_dataframe = uncor_repeat_data_frame, small_dataframe = NN_data, intersect_column_first = "read", intersect_column_second = "name")
  }
  
  cor_repeat_data_frame_sum_of_repeat <- sum(cor_repeat_data_frame$intersect_length)
  uncor_repeat_data_frame_sum_of_repeat  <- sum(uncor_repeat_data_frame$intersect_length)
  
  cor_repeat_data_frame_sum_of_length <- sum(cor_repeat_data_frame$end - cor_repeat_data_frame$begin)
  uncor_repeat_data_frame_sum_of_length  <- sum(uncor_repeat_data_frame$end - uncor_repeat_data_frame$begin)
  
  data_repeat_percentage <- data.frame(data_type=c('Corrected', 'Uncorrected'),
                                        SumOfRepeat=c(cor_repeat_data_frame_sum_of_repeat, uncor_repeat_data_frame_sum_of_repeat),
                                        SumOfLength=c(cor_repeat_data_frame_sum_of_length, cor_repeat_data_frame_sum_of_length)
  )
 
  return(mutate(data_repeat_percentage, percentage_of_repeat = (SumOfRepeat/SumOfLength)*100, Organism = organism))
}

#==================================================
# Analysis of repeats
library(ggplot2)
cor_repeat_yeast <- read.delim("/data/test_data_from_server_maize/proovread/yeast/4/yeast4/correct_2016-09-01/yeast/coordinate/corrected_repeat_intersect.bed", header = FALSE)
uncor_repeat_yeast <- read.delim("/data/test_data_from_server_maize/proovread/yeast/4/yeast4/correct_2016-09-01/yeast/coordinate/uncorrected_repeat_intersect.bed", header = FALSE)

cor_repeat_trypanosoma <- read.delim("/data/test_data_from_server_maize/proovread/trypanosoma/trypanosoma/correct_2016-09-01/trypanosoma/test_coordent/corrected_repeat_intersect.bed", header = FALSE)
uncor_repeat_trypanosoma <- read.delim("/data/test_data_from_server_maize/proovread/trypanosoma/trypanosoma/correct_2016-09-01/trypanosoma/test_coordent/uncorrected_repeat_intersect.bed", header = FALSE)

cor_repeat_rice <- read.delim("/data/test_data_from_server_maize/proovread/rice/correct_2016-09-06/rice/test_coordent/corrected_repeat_intersect.bed", sep="\t", header = F)
uncor_repeat_rice <- read.delim("/data/test_data_from_server_maize/proovread/rice/correct_2016-09-06/rice/test_coordent/uncorrected_repeat_intersect.bed", sep="\t", header = F)

cor_repeat_human <- read.delim("/data/test_data_from_server_maize/proovread/human21/correct_2016-09-10/human/test_coordent/corrected_repeat_intersect_21.bed", sep="\t", header = F)
uncor_repeat_human <- read.delim("/data/test_data_from_server_maize/proovread/human21/correct_2016-09-10/human/test_coordent/uncorrected_repeat_intersect_21.bed", sep="\t",  header = F)



yeast_repeat_percentage <- analysis_repeat(cor_repeat_data_frame = cor_repeat_yeast, uncor_repeat_data_frame = uncor_repeat_yeast, columns_name = column_names, organism = "Yeast")

trypanosoma_repeat_percentage <- analysis_repeat(cor_repeat_data_frame = cor_repeat_trypanosoma, uncor_repeat_data_frame = uncor_repeat_trypanosoma, columns_name = column_names,
                                                 organism = "Trypanosoma", NN_data = trypanosoma_NN_reads)

rice_repeat_percentage <- analysis_repeat(cor_repeat_data_frame = cor_repeat_rice, uncor_repeat_data_frame = uncor_repeat_rice, columns_name = column_names, organism = "Rice")
human_repeat_percentage <- analysis_repeat(cor_repeat_data_frame = cor_repeat_human, uncor_repeat_data_frame = uncor_repeat_human, columns_name = column_names, organism = "Human")


concat_data <- bind_rows(yeast_repeat_percentage, trypanosoma_repeat_percentage, rice_repeat_percentage, human_repeat_percentage)

gc_plot <- ggplot(concat_data, aes(x = data_type, y = percentage_of_repeat, fill = data_type))+
  geom_bar(stat="identity") + 
  facet_grid(~Organism) +
  labs(x = "Read Type", y = "Repeat percentage")+
  scale_fill_discrete(name  = NULL)

 ############################################
 
 # analysis of covergae
 
 corrected_coverage <- read.delim("/data/test_data_from_server_maize/proovread/yeast/4/yeast4/correct_2016-09-01/yeast/coordinate/corrected_coverage.txt", header = FALSE)
 uncorrected_covergae <-  read.delim("/data/test_data_from_server_maize/proovread/yeast/4/yeast4/correct_2016-09-01/yeast/coordinate/uncorrected_coverage.txt",  header = FALSE)
 
corrected_cov <- ggplot(data = corrected_coverage)
corrected_cov + geom_density(aes(x=V8))

uncorrected_cov <- ggplot(data = uncorrected_covergae)
uncorrected_cov + geom_density(aes(x=V8))

library(dplyr)
merged_data <- rbind(mutate(corrected_coverage, cor_type='corrected'), mutate(uncorrected_covergae, cor_type='uncorrected'))

ggplot(data = merged_data) + geom_density(aes(x = V8)) + facet_grid(cor_type~.) + labs(x = "Covergae", title = "coverage density")
ggplot(data = merged_data) + geom_density(aes(x = V8, colour = cor_type)) + labs(x = "Covergae", title = "coverage density") + scale_color_discrete(name = "corrected/uncorrected")



ggplot(data = merged_data) + geom_boxplot(aes(x = cor_type, y = V8)) + labs(x = "corrected/uncorrected", y = "covergae", title = "corrected/uncorrected comparsion ")

#####################################################

# analysis of average coverage

corrected_avg <- read.delim("/data/test_data_from_server_maize/proovread/yeast/4/yeast4/correct_2016-09-01/yeast/coordinate/corrected_coverage_average.txt", header = FALSE)
uncorrected_avg <- read.delim("/data/test_data_from_server_maize/proovread/yeast/4/yeast4/correct_2016-09-01/yeast/coordinate/uncorrected_coverage_average.txt", header = FALSE)
average_coverage_data <- rbind(mutate(corrected_avg, cor_type='corrected'), mutate(uncorrected_avg, cor_type='uncorrected'))

ggplot(data = average_coverage_data) + geom_boxplot(aes(x = cor_type, y = V2)) + labs(x = "corrected/uncorrected", y = "covergae", title = "corrected/uncorrected comparsion ")


#####################################################

# compoare mappabiity

corrected_mappability <- read.delim("/data/test_data_from_server_maize/proovread/yeast/4/yeast4/correct_2016-09-01/yeast/coordinate/corrected_read_coordinate_mappability_only_chr4.bed", header = F)
uncorrected_mappability  <- read.delim("/data/test_data_from_server_maize/proovread/yeast/4/yeast4/correct_2016-09-01/yeast/coordinate/uncorrected_read_coordinate_map_re_chr4_100.bed", header = F)
corrected_mappability <- mutate(corrected_mappability, kind = "Corrected") 
uncorrected_mappability <- mutate(uncorrected_mappability, kind = "Uncorrected")
correc_uncorrect_mappability <- rbind(corrected_mappability, uncorrected_mappability)
ggplot(data = correc_uncorrect_mappability, aes(kind, V5, fill = kind)) + geom_violin(draw_quantiles = c(0.25, 0.5, 0.75), scale = "area") +
  labs(x="Corrected/Uncorrected", y="Mappability")

# uncorrected read mappability length digram
uncor_mappability_length <- read.delim("/data/test_data_from_server_maize/proovread/yeast/4/yeast4/correct_2016-09-01/yeast/coordinate/uncorrected_read_coordinate_chr4_mappabilit.bed", header = FALSE)
head(uncor_mappability_length)
uncor_mappability_length <- mutate(uncor_mappability_length, leng=(V3-V2))
uncor_mappability_length_end_less_begin <-as.data.frame(uncor_mappability_length[uncor_mappability_length$V3 < uncor_mappability_length$V2,])
ggplot(data = uncor_mappability_length, aes(leng,V5)) + geom_jitter() + labs(x = "Length", y = "Mappability")

unc_minus_100 <- read.delim("/data/test_data_from_server_maize/proovread/yeast/4/yeast4/correct_2016-09-01/yeast/coordinate/uncorrected_read_coordinate_map_chr4_100.bed", header = FALSE)
unc_minus_100 <- mutate(unc_minus_100, diff = )


un_po <- read.delim("/data/test_data_from_server_maize/proovread/yeast/4/yeast4/correct_2016-09-01/yeast/coordinate/uncorrected_read_coordinate_map_re_chr4_100_positve_length.bed", header = F)
ggplot(data = un_po) + geom_violin(aes(x=V1, y=V5))


########################################################################

# compare error rate

cor_error_rate <- read.delim("/data/test_data_from_server_maize/proovread/yeast/4/yeast4/correct_2016-09-01/yeast/coordinate/corrected_error_rate.txt", header = F)
uncor_error_rate <- read.delim("/data/test_data_from_server_maize/proovread/yeast/4/yeast4/correct_2016-09-01/yeast/coordinate/uncorrected_error_rate.txt", header = F)

cor_error_rate <- mutate(cor_error_rate, all_change = (V3+V4),  percentage_of_mismatch = (V3/V2)*100, percentage_of_indel = (V4/V2)*100, all_percentage = (all_change/V2)*100, kind = "Corrected")
uncor_error_rate <- mutate(uncor_error_rate, all_change = (V3+V4),  percentage_of_mismatch = (V3/V2)*100, percentage_of_indel = (V4/V2)*100, all_percentage = (all_change/V2)*100,kind = "Uncorrected")

error_rate <- rbind(cor_error_rate, uncor_error_rate)

ggplot(error_rate) + geom_violin(aes(kind, percentage_of_indel, fill = kind))
ggplot(error_rate) + geom_violin(aes(kind, percentage_of_mismatch, fill = kind))
ggplot(error_rate) + geom_violin(aes(kind, all_percentage, fill = kind))
