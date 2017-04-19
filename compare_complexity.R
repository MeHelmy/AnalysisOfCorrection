# python ~/Desktop/test_scripts/UncorrectedProovread.py -c yeast.trimmed.fa -f yeast.untrimmed.fa -o uncorrected_reads.fa
# python ~/Desktop/test_scripts/Assess/CalculateComplexity.py -i uncorrected_reads.fa -o uncorrected_reads_complexity.txt
# python ~/Desktop/test_scripts/Assess/CalculateComplexity.py -i yeast.trimmed.fa -o corrected_reads_complexity.txt

# import lib
library(ggplot2)
library(plyr)
library(grid)
library(gridExtra)
library(tidyr)

# ===============================================
# Global variable

data_location <- "/data/test_data_from_server_maize/assessment_2016-08-29/complexity/"


human_NN_reads <- read.delim("/data/test_data_from_server_maize/assessment_2016-08-29/complementory_data/human21_pacbio_original_nucleotide_N_freq.bed")
rice_NN_reads <- read.delim("/data/test_data_from_server_maize/assessment_2016-08-29/complementory_data/rice_pacbio_original_nucleotide_N_freq.bed")
trypanosoma_NN_reads <- read.delim("/data/test_data_from_server_maize/assessment_2016-08-29/complementory_data/trypanosoma_pacbio_original_nucleotide_N_freq.bed")


# ===============================================
# Functions

extract_rows <- function(whole_dataframe, small_dataframe, intersect_column_first, intersect_column_second, complement=TRUE){
  if(complement){
    return(subset(whole_dataframe, !(whole_dataframe[[intersect_column_first]] %in% small_dataframe[[intersect_column_second]])))
  }else{
    return(subset(whole_dataframe, (whole_dataframe[[intersect_column_first]] %in% small_dataframe[[intersect_column_second]])))
  }
  
}

# ===============================================
# Read Data

ecoli_corrected_complexity <- read.delim(paste(data_location, 'ecoli_corrected_reads_complexity.txt', sep = ''))
ecoli_uncorrected_complexity <- read.delim(paste(data_location, 'ecoli_uncorrected_reads_complexity.txt', sep = ''))

ecoli_corrected_complexity <- separate(ecoli_corrected_complexity, read_name, c('read', 'sigment'), "\\.", extra = "merge", remove = F, fill = "right")
ecoli_uncorrected_complexity <- separate(ecoli_uncorrected_complexity, read_name, c('read', 'sigment'), "\\.", extra = "merge", remove = F, fill = "right")


yeast_corrected_complexity <- read.delim(paste(data_location, 'yeast_corrected_reads_complexity.txt', sep = ''))
yeast_uncorrected_complexity <- read.delim(paste(data_location, 'yeast_uncorrected_reads_complexity.txt', sep = ''))

trypanosoma_corrected_complexity <- read.delim(paste(data_location, 'trypanosoma_corrected_reads_complexity.txt', sep = ''))
trypanosoma_uncorrected_complexity <- read.delim(paste(data_location, 'trypanosoma_uncorrected_reads_complexity.txt', sep = ''))

trypanosoma_corrected_complexity <- separate(trypanosoma_corrected_complexity, read_name, c('read', 'sigment'), "\\.", extra = "merge", remove = F, fill = "right")
trypanosoma_corrected_complexity <- extract_rows(whole_dataframe = trypanosoma_corrected_complexity, small_dataframe = trypanosoma_NN_reads, intersect_column_first = "read", intersect_column_second = "name", complement = T)
trypanosoma_uncorrected_complexity <- separate(trypanosoma_uncorrected_complexity, read_name, c('read', 'sigment'), "\\.", extra = "merge", remove = F, fill = "right")
trypanosoma_uncorrected_complexity <- extract_rows(whole_dataframe = trypanosoma_uncorrected_complexity, small_dataframe = trypanosoma_NN_reads, intersect_column_first = "read", intersect_column_second = "name", complement = T)


rice_corrected_complexity <- read.delim(paste(data_location, 'rice_corrected_reads_complexity.txt', sep = ''))
rice_uncorrected_complexity <- read.delim(paste(data_location, 'rice_uncorrected_reads_complexity.txt', sep = ''))

human_corrected_complexity <- read.delim(paste(data_location, 'human_corrected_reads_complexity.txt', sep = ''))
human_uncorrected_complexity <- read.delim(paste(data_location, 'human_uncorrected_reads_complexity.txt', sep = ''))

# ===============================================
# Calculations

ecoli_corrected_complexity <- mutate(ecoli_corrected_complexity, Type = 'Corrected', Organism = 'E.coli')
ecoli_uncorrected_complexity <- mutate(ecoli_uncorrected_complexity, Type = 'Uncorrected', Organism = 'E.coli')

yeast_corrected_complexity <- mutate(yeast_corrected_complexity, Type = 'Corrected', Organism = 'Yeast')
yeast_uncorrected_complexity <- mutate(yeast_uncorrected_complexity, Type = 'Uncorrected', Organism = 'Yeast')

trypanosoma_corrected_complexity <- mutate(trypanosoma_corrected_complexity, Type = 'Corrected', Organism = 'Trypanosoma')
trypanosoma_uncorrected_complexity <- mutate(trypanosoma_uncorrected_complexity, Type = 'Uncorrected', Organism = 'Trypanosoma')

rice_corrected_complexity <- mutate(rice_corrected_complexity, Type = 'Corrected', Organism = 'Rice')
rice_uncorrected_complexity <- mutate(rice_uncorrected_complexity, Type = 'Uncorrected', Organism = 'Rice')

human_corrected_complexity <- mutate(human_corrected_complexity, Type = 'Corrected', Organism = 'Human')
human_uncorrected_complexity <- mutate(human_uncorrected_complexity, Type = 'Uncorrected', Organism = 'Human')

complexity_data <- rbind(ecoli_corrected_complexity, ecoli_uncorrected_complexity, yeast_corrected_complexity, yeast_uncorrected_complexity,
                         trypanosoma_corrected_complexity, trypanosoma_uncorrected_complexity,
                         rice_corrected_complexity, rice_uncorrected_complexity, human_corrected_complexity, human_uncorrected_complexity)

gc_plot <- ggplot(complexity_data, aes(x = Type, y = GC, fill = Type))+
  geom_violin(scale = "area") + 
  facet_grid(~Organism) +
  theme(axis.title.x = element_blank()) 

complexity_plot <- ggplot(complexity_data, aes(x = Type, y = complexity, fill = Type))+
  geom_violin(scale = "area") + 
  facet_grid(~Organism) +
  theme(axis.title.x = element_blank())
  
gc_complexity_compare <-  grid.arrange(gc_plot, complexity_plot, ncol = 2)

gc_complexity_compare1 <-  grid.arrange(gc_plot, complexity_plot, ncol = 1)


# yeast_corrected <- read.delim('/data/test_data_from_server_maize/proovread/yeast/4/yeast4/correct_2016-09-01/yeast/corrected_reads_complexity.txt')
# yeast_uncorrected <- read.delim('/data/test_data_from_server_maize/proovread/yeast/4/yeast4/correct_2016-09-01/yeast/uncorrected_reads_complexity.txt')
# 
# yeast_corrected <- mutate(yeast_corrected, Type = 'correct')
# yeast_uncorrected <- mutate(yeast_uncorrected, Type = 'uncorrect')
# 
# yeast <- rbind(yeast_corrected, yeast_uncorrected)
# 
# comlexity_plot <- ggplot(yeast, aes(Type, GC))
# comlexity_plot + geom_violin()
# 
# ggplot(data = yeast, aes(Type, GC))  + geom_violin(scale = "area", aes(fill = Type))
# ggplot(data = yeast, aes(Type, complexity))  + geom_violin(scale = "area", aes(fill = Type))

# test
separate(ecoli_corrected_complexity, read_name, c('read', 'sigment'), "\\.", extra = "merge", remove = F, fill = "right")
