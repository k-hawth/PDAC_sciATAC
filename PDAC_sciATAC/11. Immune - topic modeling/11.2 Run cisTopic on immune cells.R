library(cisTopic)
library(Matrix)
library(ggplot2)
set.seed(1)



########################################################################################################################
# SET VARIABLES
########################################################################################################################

# Where will outputs go?
ProjectFolder = "~/cisTopic/Primary_TSSMixed_nFrags1k_FRIS0.4/peakset_MACS2ArchR_PrimaryTSSMixnFrags1kFRIS0.4/Immune/cisTopic_50"

# Counts matrix ".values.gz" file (all 3 sparse matrix files should be in the same place)
CountsMatrixValues = "~/cisTopic/Primary_TSSMixed_nFrags1k_FRIS0.4/peakset_MACS2ArchR_PrimaryTSSMixnFrags1kFRIS0.4/Immune/ImmuneTSSMixnFrags1kFRIS0.4.filt.sparseMatrix.values.gz"

PlotPrefix = "cisTopic_alpha50_ImmuneTSSMixnFrags1kFRIS0.4_MACS2ArchR"

# cisTopic settings
cisTopic_topicCount = 10:50
cisTopic_cores = 10
cisTopic_thrP = 0.975
cisTopic_alpha = 50
cisTopic_beta = 0.1
cisTopic_iterations = 500



############################################################
# Read in sparse counts matrix

# Sparse matrix code obtained from scitools matrix-cistopic3 function
readSparseMatrix <- function(values,
                             cols = NULL,
                             rows = NULL){
    
    IN<-as.matrix(read.table(values))
    
    if(is.null(cols) & is.null(rows)){
        cols = gsub(".values.gz$",".cols.gz",values)
        rows = gsub(".values.gz$",".rows.gz",values)
    }
    
    COLS<-read.table(cols)
    ROWS<-read.table(rows)
    
    IN<-sparseMatrix(i=IN[,1],j=IN[,2],x=IN[,3])
    colnames(IN)<-COLS$V1
    row.names(IN)<-ROWS$V1
    row.names(IN)<-sub("_","-",sub("_",":",row.names(IN)))
    
    return(IN)
}

countsMatrix <- readSparseMatrix(CountsMatrixValues)

#Set up filtered binarized counts matrix
cisTopicObject <- createcisTopicObject(countsMatrix,min.cells=1,min.regions=1, keepCountsMatrix=FALSE, project.name=PlotPrefix)

# Remote intermediate objects
rm(countsMatrix)



############################################################
# Build cisTopic models

cisTopicObject <- runWarpLDAModels(cisTopicObject, topic= cisTopic_topicCount, seed=123, nCores = cisTopic_cores, alpha = cisTopic_alpha, beta = cisTopic_beta, iterations = cisTopic_iterations, addModels=FALSE, tmp = paste(PlotPrefix,".temp",sep = ""))

# Save cisTopic object
saveRDS(cisTopicObject, paste(PlotPrefix,".cistopicObject_NotBinarized.rds", sep = ""))

# Move temp to temp folder
system("mkdir cisTopic_temp_outputs")
system("mv *.temp* cisTopic_temp_outputs")

# Plot model selection
pdf(width = 7, height = 7, file = paste(PlotPrefix,"_ModelSelection.pdf",sep = ""))
par(mfrow=c(3,3))
cisTopicObject <- selectModel(cisTopicObject, type='maximum')
cisTopicObject <- selectModel(cisTopicObject, type='perplexity')
cisTopicObject <- selectModel(cisTopicObject, type='derivative')
dev.off()



############################################################
# Perplexity-based topic selection

system("mkdir Perplexity_topics")
setwd(paste(getwd(),"/Perplexity_topics/",sep=""))

# Select model and plot metrics
pdf(width = 7, height = 7, file = paste(PlotPrefix,"_ModelSelection_Perplexity.pdf",sep = ""))
cisTopicObject <- selectModel(cisTopicObject, type='perplexity')
dev.off()

topic_model = length(cisTopicObject@selected.model$topic_sums)



############################################################
# Get topic .bed files

# Binarize cisTopics
cisTopicObject <- getRegionsScores(cisTopicObject, method='NormTop', scale=TRUE)
cisTopicObject <- binarizecisTopics(cisTopicObject, thrP=cisTopic_thrP, plot=TRUE)

# Get metrics for topic-region selection
dev.off()

# Get binarized topics
getBedFiles(cisTopicObject, path = paste("topics",topic_model,"_bedfiles_binarized",sep=""))



############################################################
# Set names of cisTopic outputs

output_TopicMatrix = paste(PlotPrefix,"_Topics",topic_model,"_Matrix",".matrix",sep="")
output_TopicMatrixScaled = paste(PlotPrefix,"_Topics",topic_model,"_Matrix_Scaled",".matrix",sep="")



############################################################
# Get cell-topic score matrices

topic_weights <- cisTopicObject@selected.model$document_expects
topic_weights <- as.data.frame(topic_weights)
rownames(topic_weights) <- paste("Topic_",rownames(topic_weights),sep="")
write.table(topic_weights, file = output_TopicMatrix, row.names = T, col.names = T, sep = '\t', quote = F)

topic_weights <- cisTopicObject@selected.model$document_expects
topic_weights <- scale(topic_weights, center = TRUE, scale = TRUE)
topic_weights <- as.data.frame(topic_weights)
rownames(topic_weights) <- paste("Topic_",rownames(topic_weights),sep="")
write.table(topic_weights, file = output_TopicMatrixScaled,row.names = T, col.names = T, sep = '\t', quote = F)



############################################################
# Derivative-based topic selection

setwd(ProjectFolder)
system("mkdir Derivative_topics")
setwd(paste(getwd(),"/Derivative_topics/",sep=""))

# Plot model selection
pdf(width = 7, height = 7, file = paste(PlotPrefix,"_ModelSelection_Derivative.pdf",sep = ""))
cisTopicObject <- selectModel(cisTopicObject, type='derivative')
dev.off()

topic_model = length(cisTopicObject@selected.model$topic_sums)



############################################################
# Get topic .bed files

# Binarize cisTopics
cisTopicObject <- getRegionsScores(cisTopicObject, method='NormTop', scale=TRUE)
cisTopicObject <- binarizecisTopics(cisTopicObject, thrP=cisTopic_thrP, plot=TRUE)

# Get topic region selection metrics
dev.off()

# Get topic bed files - binarized
getBedFiles(cisTopicObject, path = paste("topics",topic_model,"_bedfiles_binarized",sep=""))



############################################################
# Set names of cisTopic outputs

output_TopicMatrix = paste(PlotPrefix,"_Topics",topic_model,"_Matrix",".matrix",sep="")
output_TopicMatrixScaled = paste(PlotPrefix,"_Topics",topic_model,"_Matrix_Scaled",".matrix",sep="")



############################################################
# Get cell-topic score matrices

topic_weights <- cisTopicObject@selected.model$document_expects
topic_weights <- as.data.frame(topic_weights)
rownames(topic_weights) <- paste("Topic_",rownames(topic_weights),sep="")
write.table(topic_weights, file = output_TopicMatrix, row.names = T, col.names = T, sep = '\t', quote = F)

topic_weights <- cisTopicObject@selected.model$document_expects
topic_weights <- scale(topic_weights, center = TRUE, scale = TRUE)
topic_weights <- as.data.frame(topic_weights)
rownames(topic_weights) <- paste("Topic_",rownames(topic_weights),sep="")
write.table(topic_weights, file = output_TopicMatrixScaled,row.names = T, col.names = T, sep = '\t', quote = F)



############################################################
# Save cisTopic R object

saveRDS(cisTopicObject, paste(PlotPrefix,".cistopicObject_Binarized.rds", sep = ""))