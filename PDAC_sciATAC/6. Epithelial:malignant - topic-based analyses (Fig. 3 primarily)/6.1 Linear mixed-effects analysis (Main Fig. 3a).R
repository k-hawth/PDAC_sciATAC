library(ArchR)
set.seed(1)
addArchRGenome("hg38")
library(parallel)
library(dplyr)


library(lme4)
library(emmeans)
library(pheatmap)



############################################################
# Load ArchRProject

proj <- loadArchRProject("~/ArchR/Epithelial_AllSamples_TSSMixed_nFrags1k_FRIS0.4_PeaksMACS3_cisTopic50")



########################################################################################################################
# Prepare cisTopic topic scores for linear mixed effects analysis
# In this section, each sample will have their epithelial/malignant cells separated by cell type.
# The mean topic score of each cell type subpopulation will then be calculated for the sample.
########################################################################################################################

# Load cisTopic data
cistopic <- as.data.frame(t(read.table("~/cisTopic/AllSamples_TSSMixed_nFrags1k_FRIS0.4/peakset_MACS3_AllSamplesTSSMixednFrags1kFRIS0.4/Epithelial/cisTopic100_Alpha50/Perplexity_topics/cisTopic100Alpha50_EpithelialAllSamplesTSSMixnFrags1kFRIS0.4_AllSamplesMACS3_Topics50_Matrix.matrix")))
rownames(cistopic) <- paste("sears_tumors_hg38",rownames(cistopic), sep = "#")



############################################################
# Function to collect mean topic value of each cell type within each sample

sampleGroup_Mean <- function(
    sample_group_cell_annot,
    cell_x_value_matrix
){
    
    final_output <- sample_group_cell_annot
    final_output <- final_output %>%
        group_by(sample) %>%
        mutate(count_total = n()) %>%
        group_by(sample,group) %>%
        mutate(count_group = n()) %>%
        select(-cell) %>%
        unique()
    
    for(feature in colnames(cell_x_value_matrix)){
        output <- sample_group_cell_annot
        output$feature <- cell_x_value_matrix[match(output$cell, rownames(cell_x_value_matrix)),feature]
        
        output <- output %>%
            group_by(sample,group) %>%
            mutate(mean_feature_score = mean(feature)) %>%
            select(-cell) %>%
            unique()
        
        final_output[[paste(feature,"_Mean",sep="")]] <- output$mean_feature_score[match(paste(final_output$sample,final_output$group),paste(output$sample,output$group))]
    }
    
    
    
    return(final_output)
    
}



############################################################
# Run function and filter for malignant subpopulations

df <- data.frame(sample = proj$Sample_ID,
                 group = proj$celltype,
                 cell = proj$cellNames)

df <- sampleGroup_Mean(df, cistopic)


# Discard tumor subpopulation if less than min cell threshold
df <- df[df$count_group>=10,]

# Filter for SIP, classical-NRP, and basal committed
df <- df[grep("Classical|Basal|Intermediate",df$group),]

# Discard malignant PAV
df <- df[df$sample %in% c("ST-00017310-T","ST-00019171-T")==F,]

# Discard metastases
df <- df[grep("-T",df$sample),]

# Convert groups to factors
df$group <- factor(df$group)



############################################################
# Run lmer() and emmeans() on each topic

emm_summary <- data.frame()
emm_emmeans_summary <- data.frame()

for(feature in grep("Topic",colnames(df),value = T)){
    equation <- paste(feature, " ~ group + (1 | sample)", sep = "")
    model <- lme4::lmer(equation, data = df)
    
    emm <- emmeans(model, pairwise ~ group, adjust = "none")
    
    # Get estimated marginal means
    emm_emmeans <- as.data.frame(emm$emmeans)
    emm_emmeans$feature <- gsub("_Mean$","",feature)
    emm_emmeans <- emm_emmeans[,c("feature",grep("^feature$",colnames(emm_emmeans),value = T, invert = T))]
    emm_emmeans_summary <- rbind(emm_emmeans_summary,emm_emmeans)
    
    # Get contrast results
    emm <- summary(emm$contrasts)
    emm$fdr <- p.adjust(emm$p.value, method = "fdr")
    emm_summary <- rbind(
        emm_summary,
        cbind(data.frame(feature = gsub("_Mean$","",feature)), emm)
    )
}



############################################################
# Identify subtype-associated topics
# A topic is subtype-associated if it is > each of the other subtypes, respectively

# Identify SIP topics
results <- emm_summary[
    (emm_summary$contrast=="Basal Committed - (Subtype-Intermediate Progenitor)" &
         emm_summary$fdr<=0.05 &
         emm_summary$t.ratio<0) |
        (emm_summary$contrast=="(Classical-NRP) - (Subtype-Intermediate Progenitor)" &
             emm_summary$fdr<=0.05 &
             emm_summary$t.ratio<0),]

results <- results %>%
    group_by(feature) %>%
    mutate(feature_count = n())
results <- results[results$feature_count==2,] # Both comparisons must be a hit

topics_SIP <- unique(results$feature)



# Identify classical-NRP topics
results <- emm_summary[
    (emm_summary$contrast=="Basal Committed - (Classical-NRP)" &
         emm_summary$fdr<=0.05 &
         emm_summary$t.ratio<0) |
        (emm_summary$contrast=="(Classical-NRP) - (Subtype-Intermediate Progenitor)" &
             emm_summary$fdr<=0.05 &
             emm_summary$t.ratio>0),]

results <- results %>%
    group_by(feature) %>%
    mutate(feature_count = n())
results <- results[results$feature_count==2,] # Both comparisons must be a hit

topics_classicalNRP <- unique(results$feature)



# Identify basal committed topics
results <- emm_summary[
    (emm_summary$contrast=="Basal Committed - (Classical-NRP)" &
         emm_summary$fdr<=0.05 &
         emm_summary$t.ratio>0) |
        (emm_summary$contrast=="Basal Committed - (Subtype-Intermediate Progenitor)" &
             emm_summary$fdr<=0.05 &
             emm_summary$t.ratio>0),]

results <- results %>%
    group_by(feature) %>%
    mutate(feature_count = n())
results <- results[results$feature_count==2,] # Both comparisons must be a hit

topics_basalCommitted <- unique(results$feature)



# Remove spurious topics
topics_basalCommitted <- topics_basalCommitted[topics_basalCommitted!="Topic_19"] # topic is specific to one sample only
topics_basalCommitted <- topics_basalCommitted[topics_basalCommitted!="Topic_24"] # topic is specific to pulmonary epithelial



############################################################
# Plot heatmap of estimated marginal means

# Collect emmean values into a matrix
heatmap <- data.frame()
for(unique_topic in unique(emm_emmeans_summary$feature)){
    for(unique_subtype in unique(emm_emmeans_summary$group)){
        heatmap[unique_subtype,unique_topic] <- emm_emmeans_summary$emmean[emm_emmeans_summary$feature==unique_topic & emm_emmeans_summary$group==unique_subtype]
    }
}

# Replace 
heatmap_signif <- heatmap
heatmap_signif[rownames(heatmap_signif),colnames(heatmap_signif)] <- ""
heatmap_signif["Basal Committed",topics_basalCommitted] <- "\u2217"
heatmap_signif["Classical-NRP",topics_classicalNRP] <- "\u2217"
heatmap_signif["Subtype-Intermediate Progenitor",topics_SIP] <- "\u2217"



p <- pheatmap(heatmap,
              clustering_distance_rows = "euclidean",
              clustering_distance_cols = "euclidean",
              clustering_method = "ward",
              
              fontsize = 6,
              angle_col = 90,
              
              cellwidth = 8,
              cellheight = 8,
              
              color = colorRampPalette(ArchR::ArchRPalettes$horizon)(1000),
              
              annotation_legend = T,
              
              treeheight_row = 15,
              treeheight_col = 15,
              
              display_numbers = heatmap_signif,
              fontsize_number = 8,
              number_color = "white",
              
              border_color = NA,
              
              labels_col = gsub("Topic_","T",colnames(heatmap))
)
dev.off()
dev.off()

cairo_pdf(width = 7.55, height = 1.01,
          file = "~/cisTopic50_EstimatedMarginalMeans_Heatmap.pdf",
          family = 'Arial Unicode MS')
p
dev.off()