# libraryies
library(plyr)
library(ggplot2)
#loading data

sv_coverage_245 <- read.table("/data/correct_pacbio_2016-03-31/analysis_from_server/245_sv_accumilation_coverage.txt")
sv_coverage_79757 <- read.table("/data/correct_pacbio_2016-03-31/analysis_from_server/79757_sv_accumilation_coverage.txt")

colnames(sv_coverage_245) <- c("frequence","coverage")
colnames(sv_coverage_79757) <- c("frequence","coverage")

# mutate 245 table
sv_coverage_245 <- mutate(sv_coverage_245, cumulative=cumsum(frequence))
sv_coverage_245 <- mutate(sv_coverage_245, cumulative_per=cumsum(frequence)/sum(frequence))
sv_coverage_245 <- mutate(sv_coverage_245, cumulative_rev_per=(1-cumulative_per))
sv_coverage_245 <- mutate(sv_coverage_245, cumulative_rev=(sum(frequence)-cumulative))
sv_coverage_245 <- mutate(sv_coverage_245, response="Tolerant")

# mutate 79757 table
sv_coverage_79757 <- mutate(sv_coverage_79757, cumulative=cumsum(frequence))
sv_coverage_79757 <- mutate(sv_coverage_79757, cumulative_per=cumsum(frequence)/sum(frequence))
sv_coverage_79757 <- mutate(sv_coverage_79757, cumulative_rev_per=(1-cumulative_per))
sv_coverage_79757 <- mutate(sv_coverage_79757, cumulative_rev=(sum(frequence)-cumulative))
sv_coverage_79757 <- mutate(sv_coverage_79757, response="Sensitive")

# both data sets
sv_coverage_both <- rbind(sv_coverage_245, sv_coverage_79757)
sv_coverage_both$response <- as.factor(sv_coverage_both$response)
sv_coverage_both$response <- factor(sv_coverage_both$response, levels = rev(levels(sv_coverage_both$response)))

# plotting

# plot cumulative percentage
ggplot(data = sv_coverage_245, aes(x=coverage, y=cumulative_rev)) + 
  geom_line() + 
  scale_x_continuous(breaks=seq(2,10,by=1), limits=c(2,10)) + 
  scale_y_continuous(breaks=seq(0,0.4,by=0.1), limits=c(0,0.4)) + 
  theme_bw()

# plot count cumultive

ggplot(data = sv_coverage_245, aes(x=coverage, y=cumulative_rev)) + 
  geom_line()+ 
  scale_x_continuous(breaks=seq(2,10,by=1), limits=c(2,10)) + 
  scale_y_continuous(breaks=seq(0,5000,by=200), limits=c(0,5000)) + 
  theme_bw()


# plot cumulative count both

ggplot(data = sv_coverage_both, aes(x=coverage, y=cumulative_rev)) + 
  geom_line(aes(color=response))+ scale_colour_manual(values = c('#0000ff', '#008000') )+
  scale_x_continuous(breaks=seq(2,10,by=1), limits=c(2,10)) + 
  scale_y_continuous(breaks=seq(0,max(sv_coverage_both$cumulative_rev),by=200), limits=c(0,max(sv_coverage_both$cumulative_rev))) + 
  theme_bw()
