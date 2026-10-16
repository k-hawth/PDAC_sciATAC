# The goal of this analysis is to calculate the mean probability of a peak being associated with SIP, classical-NRP, and basal committed subtypes in primary PDAC tumors. This probability score is used to prioritize which subtype-specific peaks are used to create a proxy gene signature for transcriptomic analysis.

# To do this, we use the predictive distribution of the cisTopic epithelial/malignant topic model. This script extracts the mean probability of a peak in each epigenetic subtype. These data are exported to be used later on for gene signature development.



############################################################

library(cisTopic)
library(Matrix)
library(ggplot2)
set.seed(1)



############################################################
# Load cisTopic object

cisTopicObject <- readRDS("~/cisTopic/AllSamples_TSSMixed_nFrags1k_FRIS0.4/peakset_MACS3_AllSamplesTSSMixednFrags1kFRIS0.4/Epithelial/cisTopic100_Alpha50/Derivative_topics/cisTopic100Alpha50_EpithelialAllSamplesTSSMixnFrags1kFRIS0.4_AllSamplesMACS3.cistopicObject_Binarized.rds")

# Get correct topic model
cisTopicObject <- selectModel(cisTopicObject, select = 50)
cisTopicObject <- getRegionsScores(cisTopicObject, method='NormTop', scale=TRUE)
cisTopicObject <- binarizecisTopics(cisTopicObject, thrP=0.975, plot=TRUE)

# Get predictive distribution: probability of each peak in each cell
pred.matrix <- predictiveDistribution(cisTopicObject)

# region (col) cell (row)
pred.matrix <- t(pred.matrix)
pred.matrix <- as.data.frame(pred.matrix)



############################################################
# Load ArchRProject

library(ArchR)
set.seed(1)
addArchRGenome("hg38")
library(parallel)
library(dplyr)

proj <- loadArchRProject("~/ArchR/Epithelial_AllSamples_TSSMixed_nFrags1k_FRIS0.4_PeaksMACS3_cisTopic50")



############################################################
# Create celltype annotation where malignant subtypes from metastases are kept separate
# We do this because we are trying to find the mean probability of each peak in each primary tumor subtype

annot <- data.frame(cell = proj$cellNames,
                    sample = proj$Sample_ID,
                    celltype = proj$celltype)

# Make version where primary/met malignant subtypes are separated
annot$celltype_sep <- annot$celltype
annot$celltype_sep[ annot$celltype=="Basal Committed" & annot$sample %in% unique(annot$sample[grep("-M",annot$sample)]) ] = "Basal Committed (Met)"
annot$celltype_sep[ annot$celltype=="Classical-NRP" & annot$sample %in% unique(annot$sample[grep("-M",annot$sample)]) ] = "Classical-NRP (Met)"
annot$celltype_sep[ annot$celltype=="Subtype-Intermediate Progenitor" & annot$sample %in% unique(annot$sample[grep("-M",annot$sample)]) ] = "Subtype-Intermediate Progenitor (Met)"
annot$celltype_sep[ annot$celltype=="Malignant PAV" & annot$sample %in% unique(annot$sample[grep("-M",annot$sample)]) ] = "Malignant PAV (Met)"

# Add to ArchRProject
proj <- addCellColData(proj, cells = annot$cell, data = annot$celltype_sep, force = T, name = "celltype_sep")



############################################################
# Get mean peak probability by malignant subtype

pred.matrix$celltype <- proj$celltype_sep[match(rownames(pred.matrix), gsub("sears_tumors_hg38\\#","",proj$cellNames))]

test <- pred.matrix %>%
    group_by(celltype) %>%
    summarise_all("mean")

test <- as.data.frame(test)
rownames(test) <- test$celltype
test <- test[,colnames(test)!="celltype",drop=F]
test <- as.data.frame(t(test))
rownames(test) <- gsub(":","_",rownames(test))
rownames(test) <- gsub("-","_",rownames(test))



# Write table
write.table(test,
            file="~/cisTopic/AllSamples_TSSMixed_nFrags1k_FRIS0.4/peakset_MACS3_AllSamplesTSSMixednFrags1kFRIS0.4/Epithelial/cisTopic100_Alpha50/Perplexity_topics/cisTopic50_Epithelial_RegionCellScores_MeanByCelltypeSep.txt",
            row.names=T,col.names=T,sep='\t',quote=F)
