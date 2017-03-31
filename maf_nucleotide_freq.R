# import lib

library(plyr)
library(dplyr)
library(ggplot2)
library(reshape2)

nucleotide_freq_data <- read.delim("/data/test_data_from_server_maize/assessment_2016-08-29/complementory_data/human21_pacbio_original_nucleotide_N_freq.bed")
nucleotide_freq_data <- mutate(nucleotide_freq_data, count_nucleotide_no_N = end - N)

histogram_N_vs_nonN <- ggplot(data = melt(nucleotide_freq_data[,c(3,8:9)]), mapping = aes(x = value)) + geom_histogram() + facet_wrap(~variable, scales = 'free_x')
