# import lib

library(plyr)
library(dplyr)
library(ggplot2)
library(magrittr)
require(reshape2)
library(RColorBrewer)
library(grid)
library(gridExtra)  # http://rstudio-pubs-static.s3.amazonaws.com/2852_379274d7c5734f979e106dcf019ec46c.html

#########################################################################################################################################################
# define global variable

data_location_correct <- "/data/test_data_from_server_maize/assessment_2016-08-29/"
data_location_original <- "/data/test_data_from_server_maize/assessment_2016-08-29/original_align/"

# reading NNNNN reads 

human_NN_reads <- read.delim("/data/test_data_from_server_maize/assessment_2016-08-29/complementory_data/human21_pacbio_original_nucleotide_N_freq.bed")
rice_NN_reads <- read.delim("/data/test_data_from_server_maize/assessment_2016-08-29/complementory_data/rice_pacbio_original_nucleotide_N_freq.bed")
trypanosoma_NN_reads <- read.delim("/data/test_data_from_server_maize/assessment_2016-08-29/complementory_data/trypanosoma_pacbio_original_nucleotide_N_freq.bed")

# number of reads from each organism

ecoli_number_of_reads <- 30364
ecoli_number_of_nucleotide <- 91173200
trypanosoma_number_of_reads <- 35025 - nrow(trypanosoma_NN_reads)
trypanosoma_number_of_nucleotide <- 105236020 - sum(trypanosoma_NN_reads$end)
yeast_number_of_reads <- 10198
yeast_number_of_nucleotide <- 30638660
rice_number_of_reads <- 243707 - nrow(rice_NN_reads)
rice_number_of_nucleotide <- 728276380 - sum(rice_NN_reads$end)
human_number_of_reads <- 311204 - nrow(human_NN_reads)
human_number_of_nucleotide <- 934199660 - sum(human_NN_reads$end)

# number of reads after correction full length

lordec_ecoli_number_of_all_corrected_reads <- 30364
lsc_ecoli_number_of_all_corrected_reads <- 30317
proovread_ecoli_number_of_all_corrected_reads <-  30360
  
lordec_trypanosoma_number_of_all_corrected_reads <- 35025
lsc_trypanosoma_number_of_all_corrected_reads <- 34965
proovread_trypanosoma_number_of_all_corrected_reads <- 35021
  
lordec_yeast_number_of_all_corrected_reads <- 10198
lsc_yeast_number_of_all_corrected_reads <- 10183
proovread_yeast_number_of_all_corrected_reads <- 10196
  
lordec_rice_number_of_all_corrected_reads <- 243707
lsc_rice_number_of_all_corrected_reads <- 243243
proovread_rice_number_of_all_corrected_reads <- 243678 
  
lordec_human_number_of_all_corrected_reads <- 311204
lsc_human_number_of_all_corrected_reads <- 267063
proovread_human_number_of_all_corrected_reads <-  311167  
########################################################################################################################################################

# define function

merge_corrected_original <- function(corrected_data_frame, original_data_frame){

  da <- plyr::rename(corrected_data_frame, c('sequence'='id'))
  db <- plyr::rename(original_data_frame, c('seq_name'='id'))
  return(plyr::rename(merge(da, db[,c(1,3:6)], by = "id", all.x=TRUE), c('id'='seq_name')))
}

merge_all_species_corrected_data <- function(lordec, proovread, pbcr, lsc, original_align, align='local'){
  if(align == "local"){
    lordec_data <- merge_corrected_original(corrected_data_frame = lordec, original_data_frame = original_align)
    proovread_data <- merge_corrected_original(corrected_data_frame = proovread, original_data_frame = original_align)
    pbcr_data <- merge_corrected_original(corrected_data_frame = pbcr, original_data_frame = original_align)
    return(rbind(mutate(lordec_data, software='LoRDEC'), mutate(proovread_data, software='proovread'), mutate(pbcr_data, software='PBcR')))
  } else {
    lordec_data <- merge_corrected_original(corrected_data_frame = lordec, original_data_frame = original_align)
    proovread_data <- merge_corrected_original(corrected_data_frame = proovread, original_data_frame = original_align)
    lsc_data <- merge_corrected_original(corrected_data_frame = lsc, original_data_frame = original_align)
    return(rbind(mutate(lordec_data, software='LoRDEC'), mutate(proovread_data, software='proovread'), mutate(lsc_data, software='lsc')))
  }
}

lost_data <- function(local_data_frame){
  lordec_lost_data <- mutate(subset(local_data_frame[local_data_frame$software=='LoRDEC',], !duplicated(seq_name)),lost_data=artificial_seq_length-sizeOfAllSub )
  proovread_lost_data <- mutate(subset(local_data_frame[local_data_frame$software=='proovread',], !duplicated(seq_name)) ,lost_data=artificial_seq_length-sizeOfAllSub)
  pbcr_lost_data <- mutate(subset(local_data_frame[local_data_frame$software=='PBcR',], !duplicated(seq_name)) ,lost_data=artificial_seq_length-sizeOfAllSub)
  
  return( bind_rows(lordec_lost_data, proovread_lost_data, pbcr_lost_data) )
}

plot_efficiency_distribtion <- function(dist_data_frame, organism, align){
  plot_title <- ""
  if(align == "local"){
    if(organism == "ecoli"){
      plot_title <- "E.coli local alignment efficiency"
    }else if(organism == "trypanosoma"){
      plot_title <- "Trypanosoma local alignment efficiency"
    }else if(organism == "yeast"){
      plot_title <- "Yeast local alignment efficiency"
    }else if(organism == "rice"){
      plot_title <- "Rice local alignment efficiency"
    }else{
      plot_title <- "Human local alignment efficiency"
    }
  }else
    {
    if(organism == "ecoli"){
      plot_title <- "E.coli global alignment efficiency"
    }else if(organism == "trypanosoma"){
      plot_title <- "Trypanosoma global alignment efficiency"
    }else if(organism == "yeast"){
      plot_title <- "Yeast global alignment efficiency"
    }else if(organism == "rice"){
      plot_title <- "Rice global alignment efficiency"
    }else{
      plot_title <- "Human global alignment efficiency"
    }
  }
     
  return( ggplot(data = dist_data_frame, aes(software, subsimilarity))  +
            geom_violin(scale = "area", aes(fill = software)) + my_scale_manual_fill(name = "Software", color_manual = myColors) +
            labs(title = plot_title, x = "Software", y = "Efficiency")+
           scale_y_continuous(breaks=seq(0,100,by=25), limits=c(0,100)) + theme_bw() )
}

plot_length_distribution <- function(organism_data_frame, organism){
  if(organism == "ecoli"){
    plot_title <- "E.coli distribution of corrected sequence length"
  }else if(organism == "trypanosoma"){
    plot_title <- "Trypanosoma distribution of corrected sequence length"
  }else if(organism == "yeast"){
    plot_title <- "Yeast distribution of corrected sequence length"
  }else if(organism == "rice"){
    plot_title <- "Rice distribution of corrected sequence length"
  }else{
    plot_title <- "Human distribution of corrected sequence length"
  }
    
  return(ggplot(organism_data_frame, aes(subLength)) + 
    geom_density(aes(fill = software, color = software), alpha = 0.1) + 
      my_scale_manual_fill(name = "Software", color_manual = myColors) + my_scale_manual_color(name = "Software", color_manual = myColors)+
    labs(title = plot_title , x = "sequence after correction"))
}

plot_corrected_length_vs_efficiency <- function(organism_data_frame, organism, align){
  if(align == "local"){
    if(organism == "ecoli"){
      plot_title <- "E.coli local length vs correction efficience"
    }else if(organism == "trypanosoma"){
      plot_title <- "Trypanosoma local length vs correction efficiency"
    }else if(organism == "yeast"){
      plot_title <- "Yeast local length vs correction efficiency"
    }else if(organism == "rice"){
      plot_title <- "Rice local length vs correction efficiency"
    }else{
      plot_title <- "Human local length vs correction efficiency"
    }
  }else
  {
    if(organism == "ecoli"){
      plot_title <- "E.coli global length vs correction efficiency"
    }else if(organism == "trypanosoma"){
      plot_title <- "Trypanosoma global length vs correction efficiency"
    }else if(organism == "yeast"){
      plot_title <- "Yeast global length vs correction efficiency"
    }else if(organism == "rice"){
      plot_title <- "Rice global length vs correction efficiency"
    }else{
      plot_title <- "Human global length vs correction efficiency"
    }
  }
  
  return( 
        ggplot(organism_data_frame, aes(subLength, subsimilarity)) +
        geom_point(aes(color=software), size = 1) + facet_grid(software ~ .) + my_scale_manual_color(name = "Software", color_manual = myColors) +
          
          # ggplot(ecoli_local, aes(subLength, subsimilarity)) +
          # geom_bin2d(binwidth = c(50,2)) + 
          # facet_grid(software ~ .) + 
          # scale_fill_gradientn(colours=rainbow(20)) +
         labs(title = plot_title, x = "Correction length", y = "Efficiency")+ scale_y_continuous(breaks=seq(0,100,by=25), limits=c(0,100)) + theme_bw() ) 
}

plot_read_percentage <- function(local_data_frame_ecoli, local_data_frame_trypanosoma,
                                 local_data_frame_yeast, local_data_frame_rice, local_data_frame_human , number_of_reads_ecoli,
                                 number_of_reads_trypanosoma, number_of_reads_yeast, number_of_reads_rice, number_of_reads_human){
  # if(organism == "ecoli"){
  #   plot_title <- "E.coli corrected read number percentage"
  # }else if(organism == "trypanosoma"){
  #   plot_title <- "Trypanosoma corrected read number percentage"
  # }else if(organism == "yeast"){
  #   plot_title <- "Yeast corrected read number percentage"
  # }else if(organism == "rice"){
  #   plot_title <- "Rice corrected read number percentage"
  # }else{
  #   plot_title <- "Human corrected read number percentage"
  # }
  
  # first extract reads that contains NN
  
  
  lordec_ecoli <- (nrow(subset(local_data_frame_ecoli[local_data_frame_ecoli$software=='LoRDEC',], !duplicated(seq_name)))/number_of_reads_ecoli)*100
  pbcr_ecoli <- (nrow(subset(local_data_frame_ecoli[local_data_frame_ecoli$software=='PBcR',], !duplicated(seq_name)))/number_of_reads_ecoli)*100
  proovread_ecoli <- (nrow(subset(local_data_frame_ecoli[local_data_frame_ecoli$software=='proovread',], !duplicated(seq_name)))/number_of_reads_ecoli)*100
  
  lordec_trypanosoma <- (nrow(subset(local_data_frame_trypanosoma[local_data_frame_trypanosoma$software=='LoRDEC',], !duplicated(seq_name)))/number_of_reads_trypanosoma)*100
  pbcr_trypanosoma <- (nrow(subset(local_data_frame_trypanosoma[local_data_frame_trypanosoma$software=='PBcR',], !duplicated(seq_name)))/number_of_reads_trypanosoma)*100
  proovread_trypanosoma <- (nrow(subset(local_data_frame_trypanosoma[local_data_frame_trypanosoma$software=='proovread',], !duplicated(seq_name)))/number_of_reads_trypanosoma)*100
  
  lordec_yeast <- (nrow(subset(local_data_frame_yeast[local_data_frame_yeast$software=='LoRDEC',], !duplicated(seq_name)))/number_of_reads_yeast)*100
  pbcr_yeast <- (nrow(subset(local_data_frame_yeast[local_data_frame_yeast$software=='PBcR',], !duplicated(seq_name)))/number_of_reads_yeast)*100
  proovread_yeast <- (nrow(subset(local_data_frame_yeast[local_data_frame_yeast$software=='proovread',], !duplicated(seq_name)))/number_of_reads_yeast)*100
  
  lordec_rice <- (nrow(subset(local_data_frame_rice[local_data_frame_rice$software=='LoRDEC',], !duplicated(seq_name)))/number_of_reads_rice)*100
  pbcr_rice <- (nrow(subset(local_data_frame_rice[local_data_frame_rice$software=='PBcR',], !duplicated(seq_name)))/number_of_reads_rice)*100
  proovread_rice <- (nrow(subset(local_data_frame_rice[local_data_frame_rice$software=='proovread',], !duplicated(seq_name)))/number_of_reads_rice)*100
  
  lordec_human <- (nrow(subset(local_data_frame_human[local_data_frame_human$software=='LoRDEC',], !duplicated(seq_name)))/number_of_reads_human)*100
  pbcr_human <- (nrow(subset(local_data_frame_human[local_data_frame_human$software=='PBcR',], !duplicated(seq_name)))/number_of_reads_human)*100
  proovread_human <- (nrow(subset(local_data_frame_human[local_data_frame_human$software=='proovread',], !duplicated(seq_name)))/number_of_reads_human)*100
  
  organism_correct_percentage <- data.frame(Tool=c('LoRDEC', 'PBcR', 'proovread'),
                                            E.coli=c(lordec_ecoli,pbcr_ecoli,proovread_ecoli),
                                            Trypanosma=c(lordec_trypanosoma,pbcr_trypanosoma,proovread_trypanosoma),
                                            Yeast=c(lordec_yeast,pbcr_yeast,proovread_yeast),
                                            Rice=c(lordec_rice,pbcr_rice,proovread_rice),
                                            Human=c(lordec_human,pbcr_human,proovread_human)
  )
  
  return( melt(organism_correct_percentage, id='Tool', variable.name='Organism' , value.name='Efficiency'))
}

get_full_read_percentage <- function(ecoli_lordec, ecoli_proovread, ecoli_lsc,
                                     yeast_lordec, yeast_proovread, yeast_lsc,
                                     trypanosoma_lordec, trypanosoma_proovread, trypanosoma_lsc,
                                     rice_lordec, rice_proovread, rice_lsc,
                                     human_lordec, human_proovread, human_lsc,
                                     ecoli_number_of_original_reads, yeast_number_of_original_reads, trypanosoma_number_of_original_reads,
                                     rice_number_of_original_reads, human_number_of_original_reads){
  
# remeber to remove the reads that contains N from trypanosoma, rics and human
  
  
organism_correct_percentage <- data.frame(Tool=c('LoRDEC', 'proovread', 'lsc'),
                                            E.coli=c((ecoli_lordec/ecoli_number_of_original_reads)*100, (ecoli_proovread/ecoli_number_of_original_reads)*100, (ecoli_lsc/ecoli_number_of_original_reads)*100),
                                            Trypanosma=c((trypanosoma_lordec/trypanosoma_number_of_original_reads)*100, (trypanosoma_proovread/trypanosoma_number_of_original_reads)*100, (trypanosoma_lsc/trypanosoma_number_of_original_reads)*100),
                                            Yeast=c((yeast_lordec/yeast_number_of_original_reads)*100, (yeast_proovread/yeast_number_of_original_reads)*100, (yeast_lsc/yeast_number_of_original_reads)*100),
                                            Rice=c((rice_lordec/rice_number_of_original_reads)*100, (rice_proovread/rice_number_of_original_reads)*100, (rice_lsc/rice_number_of_original_reads)*100),
                                            Human=c((human_lordec/human_number_of_original_reads)*100, (human_proovread/human_number_of_original_reads)*100, (human_lsc/human_number_of_original_reads)*100))
return( melt(organism_correct_percentage, id='Tool', variable.name='Organism' , value.name='Efficiency'))
}

plot_nucleotide_percentage <- function(local_data_frame_ecoli, local_data_frame_trypanosoma,
                                       local_data_frame_yeast, local_data_frame_rice, local_data_frame_human , number_of_nucleotide_ecoli,
                                       number_of_nucleotide_trypanosoma, number_of_nucleotide_yeast, number_of_nucleotide_rice, number_of_nucleotide_human){
  
  lordec_ecoli <- (sum(filter(local_data_frame_ecoli,software=='LoRDEC')$subLength)/number_of_nucleotide_ecoli)*100
  pbcr_ecoli <- (sum(filter(local_data_frame_ecoli,software=='PBcR')$subLength)/number_of_nucleotide_ecoli)*100
  proovread_ecoli <- (sum(filter(local_data_frame_ecoli,software=='proovread')$subLength)/number_of_nucleotide_ecoli)*100
  
  lordec_trypanosoma <- (sum(filter(local_data_frame_trypanosoma,software=='LoRDEC')$subLength)/number_of_nucleotide_trypanosoma)*100
  pbcr_trypanosoma <- (sum(filter(local_data_frame_trypanosoma,software=='PBcR')$subLength)/number_of_nucleotide_trypanosoma)*100
  proovread_trypanosoma <- (sum(filter(local_data_frame_trypanosoma,software=='proovread')$subLength)/number_of_nucleotide_trypanosoma)*100
  
  lordec_yeast <- (sum(filter(local_data_frame_yeast,software=='LoRDEC')$subLength)/number_of_nucleotide_yeast)*100
  pbcr_yeast <- (sum(filter(local_data_frame_yeast,software=='PBcR')$subLength)/number_of_nucleotide_yeast)*100
  proovread_yeast <- (sum(filter(local_data_frame_yeast,software=='proovread')$subLength)/number_of_nucleotide_yeast)*100
  
  lordec_rice <- (sum(filter(local_data_frame_rice,software=='LoRDEC')$subLength)/number_of_nucleotide_rice)*100
  pbcr_rice <- (sum(filter(local_data_frame_rice,software=='PBcR')$subLength)/number_of_nucleotide_rice)*100
  proovread_rice <- (sum(filter(local_data_frame_rice,software=='proovread')$subLength)/number_of_nucleotide_rice)*100
  
  lordec_human <- (sum(filter(local_data_frame_human,software=='LoRDEC')$subLength)/number_of_nucleotide_human)*100
  pbcr_human <- (sum(filter(local_data_frame_human,software=='PBcR')$subLength)/number_of_nucleotide_human)*100
  proovread_human <- (sum(filter(local_data_frame_human,software=='proovread')$subLength)/number_of_nucleotide_human)*100
  
  organism_correct_percentage <- data.frame(Tool=c('LoRDEC', 'PBcR', 'proovread'),
                                            E.coli=c(lordec_ecoli,pbcr_ecoli,proovread_ecoli),
                                            Trypanosma=c(lordec_trypanosoma,pbcr_trypanosoma,proovread_trypanosoma),
                                            Yeast=c(lordec_yeast,pbcr_yeast,proovread_yeast),
                                            Rice=c(lordec_rice,pbcr_rice,proovread_rice),
                                            Human=c(lordec_human,pbcr_human,proovread_human)
  )
  
  return( melt(organism_correct_percentage, id='Tool', variable.name='Organism' , value.name='Efficiency'))
  
}

compare_length <- function(organism_data_frame, organism_origina_align_data_frame){
  original <- select(organism_origina_align_data_frame, subLength = artificial_seq_length)
  original <- mutate(original, Type = "Original")
  return(bind_rows(select(organism_data_frame, subLength, Type = software), original))
}

# create global corrected reads vs original dataframe
get_global_alignment_comarison <- function(golbal_data_frame, original_data_frame){
  
  golbal_correction <- select(golbal_data_frame, subsimilarity, software)
  global_before_correction <- select(original_data_frame, global_align_value)
  global_before_correction["software"] <- c("Original reads")
  global_before_correction <- plyr::rename(global_before_correction, c("global_align_value" = "subsimilarity")  )
  
  return(rbind(golbal_correction, global_before_correction))
}

# draw plot for comparion the correction rate between the original uncorrected reads and corrected
compare_global_correction_plot <- function(compare_data_frame, title, x_axes, y_axes, bins = 30, plot_colour = myColors){
  
  #change order of legend

  compare_data_frame$software <- factor(compare_data_frame$software, levels=c("Original reads", "LoRDEC" ,  "lsc", "proovread"))

  return( ggplot(data = compare_data_frame, aes(subsimilarity, colour = software )) + geom_freqpoly(binwidth = bins) + theme_bw() + scale_x_continuous(breaks=seq(0,105,by=20), limits=c(0,105)) +
            labs(title = title, x = x_axes, y = y_axes , colour = "comnpare") +  my_scale_manual_color(name = "Software and original reads", color_manual = myColors_original) )
}

# function to draw the count fo the clipped reads

plot_clipped_length_distribution <- function(clipped_data_frame,  x_axes, legen_name, bins = 30, organism, manual_color, aggregate_data = TRUE , aggregate_value = 5000, clipped_filed){
  if(organism == "ecoli"){
    plot_title <- "E.coli clipped data"
  }else if(organism == "trypanosoma"){
    plot_title <- "Trypanosoma clipped data"
  }else if(organism == "yeast"){
    plot_title <- "Yeast clipped data"
  }else if(organism == "rice"){
    plot_title <- "Rice clipped data"
  }else{
    plot_title <- "Human clipped data"
  }
  if(aggregate_data){
    clipped_data_frame[clipped_data_frame[[clipped_filed]] > aggregate_value, clipped_filed] <- aggregate_value
  }
  
 return( ggplot(clipped_data_frame, aes(lost_data, color=software, fill=software)) +
    geom_freqpoly(binwidth = bins) + labs(title = plot_title, x = x_axes) +
    my_scale_manual_color(name = legen_name, color_manual = manual_color) + my_scale_manual_fill(name = legen_name, color_manual = manual_color)+ theme_bw() )
}

# function to plot the relarion between the lengt before and after correction
plot_corrected_length_vs_original <- function(clipped_data_frame,  x_axes = "Read length", legen_name = "Type", bins = 30, organism, manual_color = myColors, aggregate_data = TRUE , aggregate_value = 5000, clipped_filed)
{
  if(organism == "ecoli")
    {
    plot_title <- "E.coli Distribution of length after correction with the original length"
    }
    else if(organism == "trypanosoma")
    {
    plot_title <- "Trypanosoma Distribution of length after correction with the original length"
    }
    else if(organism == "yeast")
    {
    plot_title <- "Yeast Distribution of length after correction with the original length"
    }else if(organism == "rice")
    {
    plot_title <- "Rice Distribution of length after correction with the original length"
    } 
    else
    {
    plot_title <- "Human Distribution of length after correction with the original length"
    }
   
  #change order of legend
  clipped_data_frame$Type <- factor(clipped_data_frame$Type, levels=c("Original", "LoRDEC" ,  "PBcR", "proovread"))
  
  if(aggregate_data){
    clipped_data_frame[clipped_data_frame[[clipped_filed]] > aggregate_value, clipped_filed] <- aggregate_value
  }
  
  return(ggplot(clipped_data_frame, aes(subLength, color=Type, fill=Type)) +
    geom_freqpoly(binwidth = bins) + labs(title = plot_title, x = x_axes) +
    my_scale_manual_color(name = legen_name, color_manual = manual_color) + my_scale_manual_fill(name = legen_name, color_manual = manual_color)+ theme_bw())
}

# function to extrct complement rows from two dataframes

extract_rows <- function(whole_dataframe, small_dataframe, intersect_column_first, intersect_column_second, complement=TRUE){
  if(complement){
    return(subset(whole_dataframe, !(whole_dataframe[[intersect_column_first]] %in% small_dataframe[[intersect_column_second]])))
  }else{
    return(subset(whole_dataframe, (whole_dataframe[[intersect_column_first]] %in% small_dataframe[[intersect_column_second]])))
  }
  
}


# function to get the percentage of reads that was corrected with value 100%

similarity_great_than_99 <- function(my_data_frame, align_type = 'local', round_number=10){
  
  if(align_type == "local"){
    lordec <- (nrow(my_data_frame[my_data_frame$software == "LoRDEC" & my_data_frame$subsimilarity > 99,]) / nrow(my_data_frame[my_data_frame$software == "LoRDEC" ,]))*100
    proovread <- (nrow(my_data_frame[my_data_frame$software == "proovread" & my_data_frame$subsimilarity > 99,]) / nrow(my_data_frame[my_data_frame$software == "proovread" ,]))*100
    pbcr <- (nrow(my_data_frame[my_data_frame$software == "PBcR" & my_data_frame$subsimilarity > 99,]) / nrow(my_data_frame[my_data_frame$software == "PBcR" ,]))*100
    print(paste("Local percentage higher than 100% : lordec ->" , round(lordec, digits = round_number) , " proovread -> " , round(proovread, digits = round_number)  , "pbcr -> " , round(pbcr, digits = round_number) ))
  }
  else if(align_type == "global"){
    lordec <- (nrow(my_data_frame[my_data_frame$software == "LoRDEC" & my_data_frame$subsimilarity > 99,]) / nrow(my_data_frame[my_data_frame$software == "LoRDEC" ,]))*100
    proovread <- (nrow(my_data_frame[my_data_frame$software == "proovread" & my_data_frame$subsimilarity > 99,]) / nrow(my_data_frame[my_data_frame$software == "proovread" ,]))*100
    lsc <- (nrow(my_data_frame[my_data_frame$software == "lsc" & my_data_frame$subsimilarity > 99,]) / nrow(my_data_frame[my_data_frame$software == "lsc" ,]))*100
    print(paste("Global percentage higher than 100% : lordec ->" , round(lordec, digits = round_number) , " proovread -> " , round(proovread, digits = round_number)  , "lsc -> " , round(lsc, digits = round_number) ) )
  }
  
}

similarity_equal_great_than_number <- function(my_data_frame, align_type = 'local', round_number=10, compare_value=100){
  
  if(align_type == "local"){
    lordec <- (nrow(my_data_frame[my_data_frame$software == "LoRDEC" & my_data_frame$subsimilarity >= compare_value,]) / nrow(my_data_frame[my_data_frame$software == "LoRDEC" ,]))*100
    proovread <- (nrow(my_data_frame[my_data_frame$software == "proovread" & my_data_frame$subsimilarity >= compare_value,]) / nrow(my_data_frame[my_data_frame$software == "proovread" ,]))*100
    pbcr <- (nrow(my_data_frame[my_data_frame$software == "PBcR" & my_data_frame$subsimilarity >= compare_value,]) / nrow(my_data_frame[my_data_frame$software == "PBcR" ,]))*100
    print(paste("Local percentage higher than Or equal", compare_value, " : lordec ->" , round(lordec, digits = round_number) , " proovread -> " , round(proovread, digits = round_number)  , "pbcr -> " , round(pbcr, digits = round_number) ))
  }
  else if(align_type == "global"){
    lordec <- (nrow(my_data_frame[my_data_frame$software == "LoRDEC" & my_data_frame$subsimilarity >= compare_value,]) / nrow(my_data_frame[my_data_frame$software == "LoRDEC" ,]))*100
    proovread <- (nrow(my_data_frame[my_data_frame$software == "proovread" & my_data_frame$subsimilarity >= compare_value,]) / nrow(my_data_frame[my_data_frame$software == "proovread" ,]))*100
    lsc <- (nrow(my_data_frame[my_data_frame$software == "lsc" & my_data_frame$subsimilarity >= compare_value,]) / nrow(my_data_frame[my_data_frame$software == "lsc" ,]))*100
    print(paste("Global percentage higher than Or equal", compare_value, " : lordec ->" , round(lordec, digits = round_number) , " proovread -> " , round(proovread, digits = round_number)  , "lsc -> " , round(lsc, digits = round_number) ) )
  }
  
}

# color function


myColors <- c(brewer.pal(4,"Set1"), "#000000")
names(myColors) <- c("LoRDEC",  "lsc", "proovread", "PBcR", "Original")

myColors_original <- c(brewer.pal(4,"Set1"), "#000000")
names(myColors_original) <- c("LoRDEC",  "lsc", "proovread", "PBcR", "Original reads")

# myColors_original <- brewer.pal(5,"Set1")
# names(myColors_original) <- c("LoRDEC", "proovread", "PBcR", "lsc", "Original reads")

my_scale_manual_fill <- function(name = "Software", color_manual, ...) {
  scale_fill_manual(values=color_manual, name = name, ...)
}

my_scale_manual_color <- function(name = "Software", color_manual, ...) {
  scale_colour_manual(values=color_manual, name = name, ...)
}



#######################################################################################################################################################
# import data

####### LoRDEC ###########
#******* local **********#
lordec_ecoli_local <- read.delim(paste(data_location_correct, 'lordec_ecoli_local.txt', sep = ''))
lordec_trypanosoma_local <-   read.delim(paste(data_location_correct, 'lordec_trypanosoma_local.txt', sep = ''))
lordec_trypanosoma_local <- extract_rows(whole_dataframe = lordec_trypanosoma_local, small_dataframe = trypanosoma_NN_reads, intersect_column_first = "sequence", intersect_column_second = "name", complement = T)
lordec_yeast_local <-   read.delim(paste(data_location_correct, 'lordec_yeast_local.txt', sep = ''))
lordec_rice_local <- read.delim(paste(data_location_correct, 'lordec_rice_local.txt', sep = ''))
lordec_rice_local <- extract_rows(whole_dataframe = lordec_rice_local, small_dataframe = rice_NN_reads, intersect_column_first = "sequence", intersect_column_second = "name", complement = T)
lordec_human_local <-   read.delim(paste(data_location_correct, 'lordec_human_local.txt', sep = ''))
lordec_human_local <- extract_rows(whole_dataframe = lordec_human_local, small_dataframe = human_NN_reads, intersect_column_first = "sequence", intersect_column_second = "name", complement = T)

#******* global **********#
lordec_ecoli_global <- read.delim(paste(data_location_correct, 'lordec_ecoli_global.txt', sep = ''))
lordec_trypanosoma_global <-   read.delim(paste(data_location_correct, 'lordec_trypanosoma_global.txt', sep = ''))
lordec_trypanosoma_global <- extract_rows(whole_dataframe = lordec_trypanosoma_global, small_dataframe = trypanosoma_NN_reads, intersect_column_first = "sequence", intersect_column_second = "name", complement = T)
lordec_yeast_global <-   read.delim(paste(data_location_correct, 'lordec_yeast_global.txt', sep = ''))
lordec_rice_global <- read.delim(paste(data_location_correct, 'lordec_rice_global.txt', sep = ''))
lordec_rice_global <- extract_rows(whole_dataframe = lordec_rice_global, small_dataframe = rice_NN_reads, intersect_column_first = "sequence", intersect_column_second = "name", complement = T)
lordec_human_global <-   read.delim(paste(data_location_correct, 'lordec_human_global.txt', sep = ''))
lordec_human_global <- extract_rows(whole_dataframe = lordec_human_global, small_dataframe = human_NN_reads, intersect_column_first = "sequence", intersect_column_second = "name", complement = T)

####### proovread ########
#******* local **********#
proovread_ecoli_local <- read.delim(paste(data_location_correct, 'proovread_ecoli_local.txt', sep = ''))
proovread_trypanosoma_local <-   read.delim(paste(data_location_correct, 'proovread_trypanosoma_local.txt', sep = ''))
proovread_trypanosoma_local <- extract_rows(whole_dataframe = proovread_trypanosoma_local, small_dataframe = trypanosoma_NN_reads, intersect_column_first = "sequence", intersect_column_second = "name", complement = T)
proovread_yeast_local <-   read.delim(paste(data_location_correct, 'proovread_yeast_local.txt', sep = ''))
proovread_rice_local <- read.delim(paste(data_location_correct, 'proovread_rice_local.txt', sep = ''))
proovread_rice_local <- extract_rows(whole_dataframe = proovread_rice_local, small_dataframe = rice_NN_reads, intersect_column_first = "sequence", intersect_column_second = "name", complement = T)
proovread_human_local <-   read.delim(paste(data_location_correct, 'proovread_human_local.txt', sep = ''))
proovread_human_local <- extract_rows(whole_dataframe = proovread_human_local, small_dataframe = human_NN_reads, intersect_column_first = "sequence", intersect_column_second = "name", complement = T)

#******* global **********#
proovread_ecoli_global <- read.delim(paste(data_location_correct, 'proovread_ecoli_global.txt', sep = ''))
proovread_trypanosoma_global <-   read.delim(paste(data_location_correct, 'proovread_trypanosoma_global.txt', sep = ''))
proovread_trypanosoma_global <- extract_rows(whole_dataframe = proovread_trypanosoma_global, small_dataframe = trypanosoma_NN_reads, intersect_column_first = "sequence", intersect_column_second = "name", complement = T)
proovread_yeast_global <-   read.delim(paste(data_location_correct, 'proovread_yeast_global.txt', sep = ''))
proovread_rice_global <- read.delim(paste(data_location_correct, 'proovread_rice_global.txt', sep = ''))
proovread_rice_global <- extract_rows(whole_dataframe = proovread_rice_global, small_dataframe = rice_NN_reads, intersect_column_first = "sequence", intersect_column_second = "name", complement = T)
proovread_human_global <-   read.delim(paste(data_location_correct, 'proovread_human_global.txt', sep = ''))
proovread_human_global <- extract_rows(whole_dataframe = proovread_human_global, small_dataframe = human_NN_reads, intersect_column_first = "sequence", intersect_column_second = "name", complement = T)

####### PBcR #############
#******* local **********#
pbcr_ecoli_local <- read.delim(paste(data_location_correct, 'pbcr_ecoli_local.txt', sep = ''))
pbcr_trypanosoma_local <-   read.delim(paste(data_location_correct, 'pbcr_trypanosoma_local.txt', sep = ''))
pbcr_trypanosoma_local <- extract_rows(whole_dataframe = pbcr_trypanosoma_local, small_dataframe = trypanosoma_NN_reads, intersect_column_first = "sequence", intersect_column_second = "name", complement = T)
pbcr_yeast_local <-   read.delim(paste(data_location_correct, 'pbcr_yeast_local.txt', sep = ''))
pbcr_rice_local <- read.delim(paste(data_location_correct, 'pbcr_rice_local.txt', sep = ''))
pbcr_rice_local <- extract_rows(whole_dataframe = pbcr_rice_local, small_dataframe = rice_NN_reads, intersect_column_first = "sequence", intersect_column_second = "name", complement = T)
pbcr_human_local <-   read.delim(paste(data_location_correct, 'pbcr_human_local.txt', sep = ''))
pbcr_human_local <- extract_rows(whole_dataframe = pbcr_human_local, small_dataframe = human_NN_reads, intersect_column_first = "sequence", intersect_column_second = "name", complement = T)

######### LSC #############
#******* global **********#
lsc_ecoli_global <- read.delim(paste(data_location_correct, 'lsc_ecoli_global.txt', sep = ''))
lsc_trypanosoma_global <-   read.delim(paste(data_location_correct, 'lsc_trypanosoma_global.txt', sep = ''))
lsc_trypanosoma_global <- extract_rows(whole_dataframe = lsc_trypanosoma_global, small_dataframe = trypanosoma_NN_reads, intersect_column_first = "sequence", intersect_column_second = "name", complement = T)
lsc_yeast_global <-   read.delim(paste(data_location_correct, 'lsc_yeast_global.txt', sep = ''))
lsc_rice_global <- read.delim(paste(data_location_correct, 'lsc_rice_global.txt', sep = ''))
lsc_rice_global <- extract_rows(whole_dataframe = lsc_rice_global, small_dataframe = rice_NN_reads, intersect_column_first = "sequence", intersect_column_second = "name", complement = T)
lsc_human_global <-   read.delim(paste(data_location_correct, 'lsc_human_global.txt', sep = ''))
lsc_human_global <- extract_rows(whole_dataframe = lsc_human_global, small_dataframe = human_NN_reads, intersect_column_first = "sequence", intersect_column_second = "name", complement = T)

##### Read original alignment #######

ecoli_original_align <- read.delim(paste(data_location_original, 'ecoli_global.txt', sep = ''))
trypanosoma_original_align <- read.delim(paste(data_location_original, 'trypanosoma_global.txt', sep = ''))
trypanosoma_original_align <- extract_rows(whole_dataframe = trypanosoma_original_align, small_dataframe = trypanosoma_NN_reads, intersect_column_first = "seq_name", intersect_column_second = "name", complement = T)
yeast_original_align <- read.delim(paste(data_location_original, 'yeast_global.txt', sep = ''))
rice_original_align <- read.delim(paste(data_location_original, 'rice_global.txt', sep = ''))
rice_original_align <- extract_rows(whole_dataframe = rice_original_align, small_dataframe = rice_NN_reads, intersect_column_first = "seq_name", intersect_column_second = "name", complement = T)
human_original_align <- read.delim(paste(data_location_original, 'human_global.txt', sep = ''))
human_original_align <- extract_rows(whole_dataframe = human_original_align, small_dataframe = human_NN_reads, intersect_column_first = "seq_name", intersect_column_second = "name", complement = T)


############## extract Number of NN reads from the corrected read number ##########################

lordec_ecoli_number_of_all_corrected_reads <- 30364
lsc_ecoli_number_of_all_corrected_reads <- 30317
proovread_ecoli_number_of_all_corrected_reads <-  30360

lordec_trypanosoma_number_of_all_corrected_reads <- nrow(lordec_trypanosoma_global)
lsc_trypanosoma_number_of_all_corrected_reads <- nrow(lsc_trypanosoma_global)
proovread_trypanosoma_number_of_all_corrected_reads <- nrow(proovread_trypanosoma_global)

lordec_rice_number_of_all_corrected_reads <- nrow(lordec_rice_global)
lsc_rice_number_of_all_corrected_reads <- nrow(lsc_rice_global)
proovread_rice_number_of_all_corrected_reads <- nrow(proovread_rice_global)

lordec_human_number_of_all_corrected_reads <- nrow(lordec_human_global)
lsc_human_number_of_all_corrected_reads <- nrow(lsc_human_global)
proovread_human_number_of_all_corrected_reads <-  nrow(proovread_human_global)

#######################################################################################################################################################


# create data frame that gathers all local and global with the before correction alignment seprated 
# ecoli

ecoli_local <- merge_all_species_corrected_data(lordec = lordec_ecoli_local, proovread = proovread_ecoli_local, pbcr = pbcr_ecoli_local, original_align = ecoli_original_align, align = 'local')
ecoli_global <- merge_all_species_corrected_data(lordec = lordec_ecoli_global, proovread = proovread_ecoli_global, lsc = lsc_ecoli_global, original_align = ecoli_original_align, align = 'global')

# trypanosoma

trypanosoma_local <- merge_all_species_corrected_data(lordec = lordec_trypanosoma_local, proovread = proovread_trypanosoma_local, pbcr = pbcr_trypanosoma_local, original_align = trypanosoma_original_align, align = 'local')
trypanosoma_global <- merge_all_species_corrected_data(lordec = lordec_trypanosoma_global, proovread = proovread_trypanosoma_global, lsc = lsc_trypanosoma_global, original_align = trypanosoma_original_align, align = 'global')

# yeast

yeast_local <- merge_all_species_corrected_data(lordec = lordec_yeast_local, proovread = proovread_yeast_local, pbcr = pbcr_yeast_local, original_align = yeast_original_align, align = 'local')
yeast_global <- merge_all_species_corrected_data(lordec = lordec_yeast_global, proovread = proovread_yeast_global, lsc = lsc_yeast_global, original_align = yeast_original_align, align = 'global')

# rice

rice_local <- merge_all_species_corrected_data(lordec = lordec_rice_local, proovread = proovread_rice_local, pbcr = pbcr_rice_local, original_align = rice_original_align, align = 'local')
rice_global <- merge_all_species_corrected_data(lordec = lordec_rice_global, proovread = proovread_rice_global, lsc = lsc_rice_global, original_align = rice_original_align, align = 'global')

# human

human_local <- merge_all_species_corrected_data(lordec = lordec_human_local, proovread = proovread_human_local, pbcr = pbcr_human_local, original_align = human_original_align, align = 'local')
human_global <- merge_all_species_corrected_data(lordec = lordec_human_global, proovread = proovread_human_global, lsc = lsc_human_global, original_align = human_original_align, align = 'global')


# create data frame for locally clipped data

ecoli_clipped_data <- lost_data(ecoli_local)
trypanosoma_clipped_data <- lost_data(trypanosoma_local)
yeast_clipped_data <- lost_data(yeast_local)
rice_clipped_data <- lost_data(rice_local)
human_clipped_data <- lost_data(human_local)

# create data frame for percentage of data lost

percentage_of_reads <- plot_read_percentage(local_data_frame_ecoli = ecoli_local, local_data_frame_trypanosoma = trypanosoma_local,
                                            local_data_frame_yeast = yeast_local, local_data_frame_rice = rice_local, local_data_frame_human = human_local,
                                            number_of_reads_ecoli = ecoli_number_of_reads, number_of_reads_trypanosoma = trypanosoma_number_of_reads, number_of_reads_yeast = yeast_number_of_reads,
                                            number_of_reads_rice = rice_number_of_reads, number_of_reads_human = human_number_of_reads)

percentage_of_nucleotides <- plot_nucleotide_percentage(local_data_frame_ecoli = ecoli_local, local_data_frame_trypanosoma = trypanosoma_local,
                                            local_data_frame_yeast = yeast_local, local_data_frame_rice = rice_local, local_data_frame_human = human_local,
                                            number_of_nucleotide_ecoli = ecoli_number_of_nucleotide, number_of_nucleotide_trypanosoma = trypanosoma_number_of_nucleotide, number_of_nucleotide_yeast = yeast_number_of_nucleotide,
                                            number_of_nucleotide_rice = rice_number_of_nucleotide, number_of_nucleotide_human = human_number_of_nucleotide)

percentage_of_full_corrected_reads <- get_full_read_percentage(ecoli_lordec = lordec_ecoli_number_of_all_corrected_reads, ecoli_proovread = proovread_ecoli_number_of_all_corrected_reads,
                                                               ecoli_lsc = lsc_ecoli_number_of_all_corrected_reads,
                                                               yeast_lordec = lordec_yeast_number_of_all_corrected_reads, yeast_proovread = proovread_yeast_number_of_all_corrected_reads,
                                                               yeast_lsc = lsc_yeast_number_of_all_corrected_reads,
                                                               trypanosoma_lordec = lordec_trypanosoma_number_of_all_corrected_reads, trypanosoma_proovread = proovread_trypanosoma_number_of_all_corrected_reads,
                                                               trypanosoma_lsc = lsc_trypanosoma_number_of_all_corrected_reads,
                                                               rice_lordec = lordec_rice_number_of_all_corrected_reads, rice_proovread = proovread_rice_number_of_all_corrected_reads,
                                                               rice_lsc = lsc_rice_number_of_all_corrected_reads,
                                                               human_lordec = lordec_human_number_of_all_corrected_reads, human_proovread = proovread_human_number_of_all_corrected_reads,
                                                               human_lsc = lsc_human_number_of_all_corrected_reads,
                                                               ecoli_number_of_original_reads = ecoli_number_of_reads, yeast_number_of_original_reads = yeast_number_of_reads, trypanosoma_number_of_original_reads = trypanosoma_number_of_reads,
                                                               rice_number_of_original_reads = rice_number_of_reads,  human_number_of_original_reads = human_number_of_reads
                                                               )

# comare original with corrected read

ecoli_compare_original_with_corrected <- compare_length(organism_data_frame = ecoli_local, organism_origina_align_data_frame = ecoli_original_align)
trypanosoma_compare_original_with_corrected <- compare_length(organism_data_frame = trypanosoma_local, organism_origina_align_data_frame = trypanosoma_original_align)
yeast_compare_original_with_corrected <- compare_length(organism_data_frame = yeast_local, organism_origina_align_data_frame = yeast_original_align)
rice_compare_original_with_corrected <- compare_length(organism_data_frame = rice_local, organism_origina_align_data_frame = rice_original_align)
human_compare_original_with_corrected <- compare_length(organism_data_frame = human_local, organism_origina_align_data_frame = human_original_align)

######################################################################################################################################################

# draw relations

# 1- distribution of efficiency 

# ecoli
ecoli_local_efficiency_plot <- plot_efficiency_distribtion(dist_data_frame = ecoli_local, organism = "ecoli", align = "local")
ecoli_global_efficiency_plot <- plot_efficiency_distribtion(dist_data_frame = ecoli_global, organism = "ecoli", align = "global")

# trypanosoma
trypanosoma_local_efficiency_plot <- plot_efficiency_distribtion(dist_data_frame = trypanosoma_local, organism = "trypanosoma", align = "local")
trypanosoma_global_efficiency_plot <- plot_efficiency_distribtion(dist_data_frame = trypanosoma_global, organism = "trypanosoma", align = "global")

# yeast
yeast_local_efficiency_plot <- plot_efficiency_distribtion(dist_data_frame = yeast_local, organism = "yeast", align = "local")
yeast_global_efficiency_plot <- plot_efficiency_distribtion(dist_data_frame = yeast_global, organism = "yeast", align = "global")

# rice
rice_local_efficiency_plot <- plot_efficiency_distribtion(dist_data_frame = rice_local, organism = "rice", align = "local")
rice_global_efficiency_plot <- plot_efficiency_distribtion(dist_data_frame = rice_global, organism = "rice", align = "global")

# human
human_local_efficiency_plot <- plot_efficiency_distribtion(dist_data_frame = human_local, organism = "human", align = "local")
human_global_efficiency_plot <- plot_efficiency_distribtion(dist_data_frame = human_global, organism = "human", align = "global")


clipped_efficienct <- grid.arrange(ecoli_local_efficiency_plot, trypanosoma_local_efficiency_plot, yeast_local_efficiency_plot,
             rice_local_efficiency_plot, human_local_efficiency_plot, ncol = 2, top = "Clipped correction efficiency")

clipped_efficienct_global <- grid.arrange(ecoli_global_efficiency_plot, trypanosoma_global_efficiency_plot, yeast_global_efficiency_plot,
                                   rice_global_efficiency_plot, human_global_efficiency_plot, ncol = 2, top = "Clipped correction efficiency global")

                                         ##################################################

# 2- clipped data

# ecoli 
draw_ecoli_clipped_data <- plot_clipped_length_distribution(clipped_data_frame = ecoli_clipped_data, x_axes = "Clipped data", legen_name = "Software", bins = 150, organism = "ecoli", manual_color = myColors, aggregate_data = TRUE, aggregate_value = 5000, clipped_filed = "lost_data")

# trypanosoma
draw_trypanosoma_clipped_data <- plot_clipped_length_distribution(clipped_data_frame = trypanosoma_clipped_data, x_axes = "Clipped data", legen_name = "Software", bins = 150, organism = "trypanosoma", manual_color = myColors, aggregate_data = TRUE, aggregate_value = 5000, clipped_filed = "lost_data")

# yeast
draw_yeast_clipped_data <- plot_clipped_length_distribution(clipped_data_frame = yeast_clipped_data, x_axes = "Clipped data", legen_name = "Software", bins = 150, organism = "yeast", manual_color = myColors, aggregate_data = TRUE, aggregate_value = 5000, clipped_filed = "lost_data")

# rice
draw_rice_clipped_data <-plot_clipped_length_distribution(clipped_data_frame = rice_clipped_data, x_axes = "Clipped data", legen_name = "Software", bins = 150, organism = "rice", manual_color = myColors, aggregate_data = TRUE, aggregate_value = 5000, clipped_filed = "lost_data")

# human
draw_human_clipped_data <- plot_clipped_length_distribution(clipped_data_frame = human_clipped_data, x_axes = "Clipped data", legen_name = "Software", bins = 150, organism = "human", manual_color = myColors, aggregate_data = TRUE, aggregate_value = 5000, clipped_filed = "lost_data")

clipped_reads_dist <- grid.arrange(draw_ecoli_clipped_data, draw_trypanosoma_clipped_data, draw_yeast_clipped_data,
                                   draw_rice_clipped_data, draw_human_clipped_data, ncol = 2, top = "Clipped sequence distribution")

                                         ##################################################

# 3- distribution of length --> please see #6 first

ecoli_length_distribution <- plot_length_distribution(ecoli_local, organism = 'ecoli')
trypanosoma_length_distribution <- plot_length_distribution(trypanosoma_local, organism = 'trypanosoma')
yeast_length_distribution <- plot_length_distribution(yeast_local, organism = 'yeast')
rice_length_distribution <- plot_length_distribution(rice_local, organism = 'rice')
human_length_distribution <- plot_length_distribution(human_local, organism = 'human')

dist_of_read_length <- grid.arrange(ecoli_length_distribution, trypanosoma_length_distribution, yeast_length_distribution,
                                    rice_length_distribution, human_length_distribution, ncol = 2, top = "sequence length distripution after correction")

                                         ##################################################

# 4- length vs efficience

# ecoli
ecoli_correction_length_local_relation <- plot_corrected_length_vs_efficiency(organism_data_frame = ecoli_local, organism = 'ecoli', align = 'local')
ecoli_correction_length_global_relation <- plot_corrected_length_vs_efficiency(organism_data_frame = ecoli_global, organism = 'ecoli', align = 'global')

# trypanosoma
trypanosoma_correction_length_local_relation <- plot_corrected_length_vs_efficiency(organism_data_frame = trypanosoma_local, organism = 'trypanosoma', align = 'local')
trypanosoma_correction_length_global_relation <- plot_corrected_length_vs_efficiency(organism_data_frame = trypanosoma_global, organism = 'trypanosoma', align = 'global')

# yeast
yeast_correction_length_local_relation <- plot_corrected_length_vs_efficiency(organism_data_frame = yeast_local, organism = 'yeast', align = 'local')
yeast_correction_length_global_relation <- plot_corrected_length_vs_efficiency(organism_data_frame = yeast_global, organism = 'yeast', align = 'global')

# rice
rice_correction_length_local_relation <- plot_corrected_length_vs_efficiency(organism_data_frame = rice_local, organism = 'rice', align = 'local')
rice_correction_length_global_relation <- plot_corrected_length_vs_efficiency(organism_data_frame = rice_global, organism = 'rice', align = 'global')

# human
human_correction_length_local_relation <- plot_corrected_length_vs_efficiency(organism_data_frame = human_local, organism = 'human', align = 'local')
human_correction_length_global_relation <- plot_corrected_length_vs_efficiency(organism_data_frame = human_global, organism = 'human', align = 'global')

correction_vs_efficienct_local <- grid.arrange(ecoli_correction_length_local_relation, trypanosoma_correction_length_local_relation, yeast_correction_length_local_relation,
                                               rice_correction_length_local_relation, human_correction_length_local_relation, ncol = 2, top = "Correction vs Length local")

correction_vs_efficienct_global <- grid.arrange(ecoli_correction_length_global_relation, trypanosoma_correction_length_global_relation, yeast_correction_length_global_relation,
                                                rice_correction_length_global_relation, human_correction_length_global_relation, ncol = 2, top = "Correction vs Length global")

                                        ##################################################
# 5- percentage of corrected reads

corrected_reads_percentage_plot <- ggplot(percentage_of_reads,aes(Organism,Efficiency, fill=Tool)) + 
  geom_bar(stat="identity", position=position_dodge()) +   geom_text(aes(label = round(Efficiency,1)) , position=position_dodge(width=0.9), vjust=-0.25, size = 3) + 
  labs(title = "Percentage of corrected reads", x = "Organism") + my_scale_manual_fill(color_manual = myColors)

corrected_nucleotides_percentage_plot <- ggplot(percentage_of_nucleotides,aes(Organism,Efficiency, fill=Tool)) + 
  geom_bar(stat="identity", position=position_dodge()) +   geom_text(aes(label = round(Efficiency,1)) , position=position_dodge(width=0.9), vjust=-0.25, size = 3) + 
  labs(title = "Percentage of nucleotides in corrected reads", x = "Organism") + my_scale_manual_fill(color_manual = myColors) + scale_y_continuous(breaks=seq(0,101,by=25), limits = c(0,100))

# percentage_of_factory <- read.csv("/data/test_worker.csv")
# percentage_of_worker <- read.csv("/data/test_factory.csv")
# factory_plot <- ggplot(percentage_of_factory,aes(Factory,factory_efficiency, fill=Country)) + 
#   geom_bar(stat="identity", position=position_dodge()) +   geom_text(aes(label = round(factory_efficiency,1)) , position=position_dodge(width=0.9), vjust=-0.25, size = 3) + 
#   labs(title = "Factory performance percentage", x = "Factory") 
# 
# worker_plot <- ggplot(percentage_of_worker,aes(Factory,worker_efficiency, fill=Country)) + 
#   geom_bar(stat="identity", position=position_dodge()) +   geom_text(aes(label = round(worker_efficiency,1)) , position=position_dodge(width=0.9), vjust=-0.25, size = 3) + 
#   labs(title = "worker performance percentage", x = "Factory") 
# 
# c <- percentage_of_factory
# c$name <- "factory"
# colnames(c)[3] <- "Efficiency"
# p <- percentage_of_worker
# p$name <- "worker"
# colnames(p)[3] <- "Efficiency"
# l <- rbind(c, p)
# ggplot(l, aes(x=Factory, y=Efficiency , group=Country, fill=Country,  colour=name)) +
#   geom_bar(stat="identity", position="dodge") + scale_colour_brewer(palette = "Set1")
# ggplot(l)+geom_bar(aes(x=Factory,y=Efficiency,group=interaction(name,Country),fill=Country),position="dodge",stat="identity")
# ggplot(l)+geom_bar(aes(x=Factory,y=Efficiency,group=interaction(name,Country),fill=Country,alpha=name),position="dodge",colour="grey",stat="identity")


# ggplot(d, aes(Organism, Efficiency)) +   
#   geom_bar(aes(fill = name), position = "dodge", stat="identity")
#   
crp <- percentage_of_reads
crp$name <- "read"
cnp <- percentage_of_nucleotides
cnp$name <- "nucleotide"
read_nucelotide_comapre <- rbind(crp, cnp)




# p <- ggplot(d, aes(Organism,Efficiency, group=Tool, fill=Tool, alpha=name)) + geom_bar(position="dodge", stat="identity", colour="black") 
# 
# ggplot(d, aes(x=Organism, y=Efficiency , group=Tool, fill=Tool, alpha=name, colour=name)) +
#   geom_bar(stat="identity", position="dodge") +scale_colour_brewer(palette = "Set1")

# http://stackoverflow.com/questions/40102676/combine-two-data-frames-in-one-graph

corrected_reads_nucleotides_percentages <- ggplot(read_nucelotide_comapre)+geom_bar(aes(x=Organism,y=Efficiency,group=interaction(name,Tool),fill=Tool,alpha=name),position="dodge",colour="grey",stat="identity")+ scale_alpha_discrete(range = c(0.5, 1), name = "Nuceltoides/Reads") +
  labs(title = "Comparing tools performance in different organisms", y = "Percentage") + my_scale_manual_fill(color_manual = myColors) + theme_bw()

# ggplot(d)+geom_bar(aes(x=Organism,y=Efficiency,group=interaction(name,Tool),fill=Tool,alpha=name),position="dodge",colour="grey",stat="identity")+ scale_alpha_discrete(limits = c(0.2, 1))
# ggplot(d, aes(x=Organism, y=Efficiency , group=Tool, fill=Tool,  colour=name)) +
#   geom_bar(stat="identity", position="dodge") + scale_colour_brewer(palette = "Set1")


# ggplot(d, aes(x=Organism, y=Efficiency , group=Tool, fill=Tool,  colour=name)) +
#   geom_bar(stat="identity", position="dodge") + scale_colour_brewer(palette = "Set1")+scale_colour_discrete(name = "Nucleotide/Read", values=c("#999999", "#E69F00"))



ggplot(read_nucelotide_comapre, aes(x = name, y = Efficiency, group = Tool, fill = Tool))+
  geom_bar(position = "dodge", stat = "identity") + 
  facet_grid(~Organism) +
  theme(axis.title.x = element_blank())

ggplot(read_nucelotide_comapre, aes(x = Tool, y = Efficiency, fill = name))+
  geom_bar(position = "dodge", stat = "identity") + 
  facet_grid(~Organism) +
  theme(axis.title.x = element_blank())

ggplot(read_nucelotide_comapre, aes(x = Organism, y = Efficiency, fill = name))+
  geom_bar(position = "dodge", stat = "identity") + 
  facet_grid(~Tool) +
  theme(axis.title.x = element_blank())

corrected_full_reads_percentage_plot <- ggplot(percentage_of_full_corrected_reads,aes(Organism,Efficiency, fill=Tool)) + 
  geom_bar(stat="identity", position=position_dodge()) +   geom_text(aes(label = round(Efficiency,1)) , position=position_dodge(width=0.9), vjust=-0.25, size = 3) + 
  labs(title = "Percentage of full corrected reads ", x = "Organism") + my_scale_manual_fill(color_manual = myColors)

  
percentage_of_correction <- grid.arrange(  corrected_reads_percentage_plot,
                                           corrected_nucleotides_percentage_plot, ncol = 2, top = "Correction percentage")


# 6- compare lengthes with original  -- > related to number 3 this is takes original read in consideration
ecoli_compare_original_with_corrected_plot <- plot_corrected_length_vs_original(clipped_data_frame = ecoli_compare_original_with_corrected, organism = "ecoli", aggregate_data = TRUE, aggregate_value = 5000, clipped_filed = "subLength", bins = 150)

trypanosoma_compare_original_with_corrected_plot <- plot_corrected_length_vs_original(clipped_data_frame = trypanosoma_compare_original_with_corrected, organism = "trypanosoma", aggregate_data = TRUE, aggregate_value = 5000, clipped_filed = "subLength", bins = 150)

yeast_compare_original_with_corrected_plot <- plot_corrected_length_vs_original(clipped_data_frame = yeast_compare_original_with_corrected, organism = "yeast", aggregate_data = TRUE, aggregate_value = 5000, clipped_filed = "subLength", bins = 150)

rice_compare_original_with_corrected_plot <- plot_corrected_length_vs_original(clipped_data_frame = rice_compare_original_with_corrected, organism = "rice", aggregate_data = TRUE, aggregate_value = 5000, clipped_filed = "subLength", bins = 150)

human_compare_original_with_corrected_plot <- plot_corrected_length_vs_original(clipped_data_frame = human_compare_original_with_corrected, organism = "human", aggregate_data = TRUE, aggregate_value = 5000, clipped_filed = "subLength", bins = 150)

length_comparison <- grid.arrange(ecoli_compare_original_with_corrected_plot, trypanosoma_compare_original_with_corrected_plot, yeast_compare_original_with_corrected_plot,
                                  rice_compare_original_with_corrected_plot, human_compare_original_with_corrected_plot, ncol = 2, top = "Correction length with original")

#########################################################################################################################################################


# print datafame  https://cran.r-project.org/web/packages/gridExtra/vignettes/tableGrob.html
# library(gridExtra)
# library(grDevices)
# t <- ttheme_default(base_size = 5)
# pdf("data_output.pdf", height=11, width=8.5)
# grid.table(head(s, n = 20), theme=t)
# dev.off()

# c1 <-c %>% group_by(seq_name)
# subset(c, !duplicated(seq_name))

# 7 compare the global length correction efficience with the original one
ecoli_global_compare <- get_global_alignment_comarison(golbal_data_frame = ecoli_global, original_data_frame = ecoli_original_align)
ecoli_correction_rate_compare <- compare_global_correction_plot(ecoli_global_compare, "E.coli corrected vs not corrected reads similarity", "Similarity", "count", bins = 1)

yeast_global_compare <- get_global_alignment_comarison(golbal_data_frame = yeast_global, original_data_frame = yeast_original_align)
yeast_correction_rate_compare <- compare_global_correction_plot(yeast_global_compare, "Yeast corrected vs not corrected reads similarity", "Similarity", "count", bins = 1)

trypanosoma_global_compare <- get_global_alignment_comarison(golbal_data_frame = trypanosoma_global, original_data_frame = trypanosoma_original_align)
trypanosoma_correction_rate_compare <- compare_global_correction_plot(trypanosoma_global_compare, "Trypanosoma corrected vs not corrected reads similarity", "Similarity", "count", bins = 1)

rice_global_compare <- get_global_alignment_comarison(golbal_data_frame = rice_global, original_data_frame = rice_original_align)
rice_correction_rate_compare <- compare_global_correction_plot(rice_global_compare, "Rice corrected vs not corrected reads similarity", "Similarity", "count", bins = 1)

human_global_compare <- get_global_alignment_comarison(golbal_data_frame = human_global, original_data_frame = human_original_align)
human_correction_rate_compare <- compare_global_correction_plot(human_global_compare, "Human corrected vs not corrected reads similarity", "Similarity", "count", bins = 1)

organism_global_compare_length <- grid.arrange(ecoli_correction_rate_compare , trypanosoma_correction_rate_compare, yeast_correction_rate_compare,
                                               rice_correction_rate_compare, human_correction_rate_compare, ncol = 2, top = "Similarity before and after correction")


###########################  TEST REGION ##############################################################################

su <- function(df){
  print(nrow(df))
  return(summary(df$subLength))
}
#  extract_rows(whole_dataframe = proovread_human_local, small_dataframe = human_NN_reads, intersect_column_first = "sequence", intersect_column_second = "name", complement = T)

before_lor_human_local <- lordec_human_local
extract_lor_human_local <- extract_rows(whole_dataframe = before_lor_human_local, small_dataframe = human_NN_reads, intersect_column_first = "sequence", intersect_column_second = "name", complement = F)
after_lor_human_local <- extract_rows(whole_dataframe = before_lor_human_local, small_dataframe = human_NN_reads, intersect_column_first = "sequence", intersect_column_second = "name", complement = T)
before_pro_human_local <- proovread_human_local
extract_pro_human_local <- extract_rows(whole_dataframe = before_pro_human_local, small_dataframe = human_NN_reads, intersect_column_first = "sequence", intersect_column_second = "name", complement = F)
after_pro_human_local <- extract_rows(whole_dataframe = before_pro_human_local, small_dataframe = human_NN_reads, intersect_column_first = "sequence", intersect_column_second = "name", complement = T)



means <- aggregate(subLength ~  software, yeast_local, mean)

ggplot(data = yeast_local, aes(software,subLength, fill=software)) + geom_boxplot() + stat_summary(fun.y=mean, colour="darkred", geom="point", shape=18, size=3,show.legend = FALSE) + geom_text(data = means, aes(label = subLength, y = subLength + 0.08)) + labs(title = "Read length comparison", subtitle = "Yeast", x = "Software", y = "Length") 






