# This section extracts peak probabilities from the stromal cisTopic model, specifically the likelihood of a peak belonging to a stromal topic. These data are exported and used later on to help convert CAF epigenetic topics into develop gene signatures.

library(cisTopic)
library(dplyr)
set.seed(1)



########################################################################################################################
# FUNCTIONS
########################################################################################################################

read_topics <- function(topic_file_directory,
                        topic_file_count,
                        topic_file_prefix = "Topic_",
                        topic_file_suffix = ".bed"){
    topic_list <- list()
    for(topic_file in paste(topic_file_prefix,1:topic_file_count,topic_file_suffix, sep = "")){
        topic_list[[topic_file]] <- read.table(file = paste(topic_file_directory,topic_file,sep = ""), header = F, sep = '\t', quote = "", comment.char = "")
    }
    return(topic_list)
}

cat_topics <- function(topic_list = topic_list){
    for(topic in names(topic_list)){
        topic_df <- topic_list[[topic]]
        topic_df$Topic = topic
        topic_list[[topic]] <- topic_df
    }
    cat <- topic_list[[1]]
    for(topic in names(topic_list)){
        cat <- rbind(cat,topic_list[[topic]])
    }
    cat <- unique(cat)
    cols <- colnames(cat)
    cols <- cols[cols!="Topic"]
    cat <- cat[c("Topic",cols)]
    return(cat)
}



########################################################################################################################
# COLLECT DATA
########################################################################################################################


############################################################
# Stromal cisTopics

catTopics <- read_topics("~/cisTopic/Primary_TSSMixed_nFrags1k_FRIS0.4/peakset_MACS2ArchR_PrimaryTSSMixnFrags1kFRIS0.4/Fibroblast/cisTopic_50/Derivative_topics/topics27_bedfiles_binarized/",
                         topic_file_count = 27,
                         topic_file_prefix = "Topic_",
                         topic_file_suffix = ".bed")
catTopics <- cat_topics(catTopics)



# Make prettier
catTopics$Topic <- gsub("\\.bed","",catTopics$Topic)
catTopics$Topic <- gsub("_"," ",catTopics$Topic)

colnames(catTopics) <- gsub("^V1$","peakChromosome",colnames(catTopics))
colnames(catTopics) <- gsub("^V2$","peakStart",colnames(catTopics))
colnames(catTopics) <- gsub("^V3$","peakEnd",colnames(catTopics))
colnames(catTopics) <- gsub("^width$","peakWidth",colnames(catTopics))



############################################################
# Add peak-topic probability

cisTopicObject <- readRDS("~/cisTopic/Primary_TSSMixed_nFrags1k_FRIS0.4/peakset_MACS2ArchR_PrimaryTSSMixnFrags1kFRIS0.4/Fibroblast/cisTopic_50/Derivative_topics/cisTopic_alpha50_FibroblastTSSMixnFrags1kFRIS0.4_MACS2ArchR.cistopicObject_Binarized.rds")

# Get correct topic model
cisTopicObject <- selectModel(cisTopicObject, select = 27)
cisTopicObject <- getRegionsScores(cisTopicObject, method='NormTop', scale=TRUE)
cisTopicObject <- binarizecisTopics(cisTopicObject, thrP=0.975, plot=TRUE)

# Get region-topic probabilities
df <- cisTopicObject@region.data




output <- data.frame()
for(unique_topic in unique(catTopics$Topic)){
    output_slice <- catTopics[catTopics$Topic==unique_topic,]
    
    cistopic_colname <- paste("Scores_", gsub(" ","", unique_topic), sep = "")
    output_slice$Peak_Topic_Probability <- df[match(
        paste(output_slice$peakChromosome,":",output_slice$peakStart,"-",output_slice$peakEnd, sep = ""),
        paste(df$seqnames,":",df$start,"-",df$end, sep = "")), cistopic_colname]
    
    output <- rbind(output, output_slice)
}



output <- output %>%
    dplyr::relocate(Peak_Topic_Probability, .after = Topic)



# Write table
write.table(output, file = "~/cisTopic/Primary_TSSMixed_nFrags1k_FRIS0.4/peakset_MACS2ArchR_PrimaryTSSMixnFrags1kFRIS0.4/Fibroblast/cisTopic_50/Derivative_topics/Stromal_cisTopic_PeakTopicProbability.txt", row.names = F, col.names = T, sep = '\t', quote = F)


