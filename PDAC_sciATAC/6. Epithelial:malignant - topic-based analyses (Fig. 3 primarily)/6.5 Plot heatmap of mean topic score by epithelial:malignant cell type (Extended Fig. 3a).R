library(ArchR)
set.seed(1)
addArchRGenome("hg38")
library(parallel)
library(dplyr)

library(pheatmap)


proj <- loadArchRProject("~/ArchR/Epithelial_AllSamples_TSSMixed_nFrags1k_FRIS0.4_PeaksMACS3_cisTopic50/")



########################################################################################################################
# Plot heatmap of mean topic scores by epithelial/malignant cell type
########################################################################################################################

cistopic <- as.data.frame(read.table("~/cisTopic/AllSamples_TSSMixed_nFrags1k_FRIS0.4/peakset_MACS3_AllSamplesTSSMixednFrags1kFRIS0.4/Epithelial/cisTopic100_Alpha50/Perplexity_topics/cisTopic100Alpha50_EpithelialAllSamplesTSSMixnFrags1kFRIS0.4_AllSamplesMACS3_Topics50_Matrix.matrix"))

# Transpose: rows = cells, cols = topics
cistopic <- t(cistopic)

# Add ArchRProjct prefix
rownames(cistopic) <- paste("sears_tumors_hg38#",rownames(cistopic),sep="")

# Reorder cells (rows) to match in ArchRProject
cistopic <- cistopic[proj$cellNames,]

cistopic <- as.data.frame(cistopic)
cistopic$cluster <- proj$celltype[match(rownames(cistopic),proj$cellNames)]

# Calculate mean
result_df <- cistopic %>%
    group_by(cluster) %>%
    summarise(across(everything(), mean, na.rm = TRUE))
result_df <- as.data.frame(result_df)
rownames(result_df) <- result_df$cluster
result_df <- result_df[,colnames(result_df)!="cluster",drop=F]
result_df <- t(result_df)

# Remove "other" subpopulations
result_df <- result_df[,grep("^Other",colnames(result_df),value=T,invert=T)]




heatmap <- as.data.frame(t(result_df))

p <- pheatmap(heatmap,
              clustering_distance_rows = "euclidean",
              clustering_distance_cols = "euclidean",
              clustering_method = "ward",
              cluster_rows = F,
              
              fontsize = 6,
              angle_col = 90,
              
              cellwidth = 8,
              cellheight = 8,
              
              color = colorRampPalette(ArchR::ArchRPalettes$comet)(1000),
              
              annotation_legend = T,
              
              treeheight_row = 15,
              treeheight_col = 15,
              
              border_color = NA,
              
              labels_col = gsub("Topic_","T",colnames(heatmap))
              
)
p