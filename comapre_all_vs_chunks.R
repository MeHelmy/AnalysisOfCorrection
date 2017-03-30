# Read data

# chunk data

yeast_local_chunk <- read.delim("/data/test_data_from_server_maize/assessment_2016-08-29/all_vs_chunks/yeast_local_chunk.txt", header = T)
yeast_local_all <- read.delim("/data/test_data_from_server_maize/assessment_2016-08-29/all_vs_chunks/yeast_local_all.txt", header = T)

yeast_global_chunk <- read.delim("/data/test_data_from_server_maize/assessment_2016-08-29/all_vs_chunks/yeast_global_chunk.txt", header = T)
yeast_global_all <- read.delim("/data/test_data_from_server_maize/assessment_2016-08-29/all_vs_chunks/yeast_global_all.txt", header = T)


yeast_local_chunk <- mutate(yeast_local_chunk, Type = "Divided Read")
yeast_local_all <- mutate(yeast_local_all , Type = "All Reads")
yeast_global_chunk <- mutate(yeast_global_chunk, Type = "Divided Read")
yeast_global_all <- mutate(yeast_global_all , Type = "All Reads")

all_raeds <- bind_rows(yeast_local_all, yeast_local_chunk)
all_raeds_global <- bind_rows(yeast_global_all, yeast_global_chunk)

ggplot(all_raeds) + geom_density(aes(x = subLength, fill = Type), alpha= 0.5) + theme(legend.title=element_blank()) + labs(x = "Read Length", y = "Density")
subread_length_local <- ggplot(all_raeds) + geom_violin(aes(x = Type, y = subLength, fill = Type)) + theme(axis.title.x=element_blank()) + scale_fill_discrete(name  = NULL) + labs(y="Subread Length", title="Distribution of read length in case of trimmed reads")
correction_length_local <- ggplot(all_raeds) + geom_violin(aes(x = Type, y = subsimilarity, fill = Type)) + theme(axis.title.x=element_blank()) + scale_fill_discrete(name  = NULL) + labs(y="Correction Efficiency", title = "Efficiency of correction in case of trimmed reads")

correction_length_global <- ggplot(all_raeds_global) + geom_violin(aes(x = Type, y = subsimilarity, fill = Type)) + theme(axis.title.x=element_blank()) + scale_fill_discrete(name  = NULL) + labs(y="Correction Efficiency", title = "Efficiency of correction in case of full reads")





number_of_reads <- data.frame(data_type = c("All reads", "Divided Reads"),
                              data_count = c(nrow(yeast_local_all), nrow(yeast_local_chunk))
)

comare_read_num <- ggplot(number_of_reads, aes(x = data_type, y = data_count, fill = data_type)) + geom_bar(stat="identity") + theme(legend.title=element_blank(), axis.title.x = element_blank())+ labs(y = "Number of reads", title = "Number of reads after correction")
grid.arrange(subread_length_local, correction_length_local, correction_length_global, comare_read_num, ncol=2)
