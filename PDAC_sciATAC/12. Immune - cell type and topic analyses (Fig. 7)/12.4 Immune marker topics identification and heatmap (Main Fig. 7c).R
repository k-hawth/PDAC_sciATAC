library(ArchR)
set.seed(1)
addArchRGenome("hg38")
library(parallel)
library(dplyr)

library(tidyr)
library(dplyr)
library(pheatmap)

library(seriation)
library(dendextend)

proj <- loadArchRProject("~/ArchR/Primary_TSSMixed_nFrags1k_FRIS0.4/Immune_MACS2ArchR_cisTopic43")



############################################################
# ggplot theme

# Theme Settings
text_settings <- element_text(colour = "black", size = 6)
line_settings <- element_line(linewidth = 0.5)

# Default Theme
default_theme <- theme(
    axis.title.x = text_settings,
    axis.title.y = text_settings,
    
    axis.ticks.x = line_settings,
    axis.ticks.y = line_settings,
    
    axis.text.x = text_settings,
    axis.text.y = text_settings,
    
    #axis.line.x = element_line(linewidth = 0.5),
    #axis.line.y = element_line(linewidth = 0.5),
    
    panel.grid = element_blank(),
    panel.background = element_blank(),
    panel.border = element_rect(linewidth = 1, fill = NA),
    
    legend.text = text_settings,
    legend.title = text_settings
)



############################################################
# Load cisTopic matrix

cistopic <- read.table("~/cisTopic/Primary_TSSMixed_nFrags1k_FRIS0.4/peakset_MACS2ArchR_PrimaryTSSMixnFrags1kFRIS0.4/Immune/cisTopic_50/Perplexity_topics/cisTopic_alpha50_ImmuneTSSMixnFrags1kFRIS0.4_MACS2ArchR_Topics43_Matrix_Scaled.matrix")

# Transpose: rows = cells, cols = topics
cistopic <- as.data.frame(t(cistopic))

# Add ArchRProjct prefix
rownames(cistopic) <- paste("sears_tumors_hg38#",rownames(cistopic),sep="")

# Reorder cells (rows) to match in ArchRProject
cistopic <- cistopic[proj$cellNames,]



############################################################
# Perform marker feature analysis
# Topic scores in one cell type vs. topic scores in all other cells

tests <- cistopic
tests$celltype <- proj$celltype[match(rownames(tests),proj$cellNames)]

total_results <- data.frame()
for(unique_celltype in unique(tests$celltype)){
    celltype <- unique_celltype
    group1 <- celltype
    group2 <- unique(proj$celltype[proj$celltype %in% group1==F])
    
    group1 <- proj$cellNames[proj$celltype %in% group1]
    group2 <- proj$cellNames[proj$celltype %in% group2]
    
    results <- data.frame()
    for(unique_topic in colnames(tests)){
        if(unique_topic!="celltype"){
            
            results[length(rownames(results))+1,"Topic"] <- unique_topic
            group1_scores <- tests[group1,unique_topic]
            group2_scores <- tests[group2,unique_topic]
            
            results[length(rownames(results)),"group1_median"] <- median(group1_scores)
            results[length(rownames(results)),"group2_median"] <- median(group2_scores)
            
            results[length(rownames(results)),"t.test_p"] <- t.test(x = group1_scores, y = group2_scores)$p.val
            
        }
    }
    results$t.test_fdr <- p.adjust(results$t.test_p, method = "fdr")
    results$diff_median <- results$group1_median-results$group2_median
    results$log2FC_median <- log2(results$group1_median/results$group2_median)
    
    results$group1 <- celltype
    total_results <- rbind(total_results,results)
}
total_results$t.test_fdr <- p.adjust(total_results$t.test_p)



############################################################
# Convert values to a celltype x topic matrix

matrix_median <- total_results[,c("Topic","log2FC_median","group1")] %>%
    pivot_wider(
        names_from = Topic,
        values_from = log2FC_median
    )
matrix_median <- as.data.frame(matrix_median)
rownames(matrix_median) <- matrix_median$group1
matrix_median <- matrix_median[,-1,drop=F]



############################################################
# Identify which results in the matrix are significant

corr_signif <- matrix_median
corr_signif[corr_signif>0|corr_signif<=0] <- ""

# Apply corrected p-values to corr_signif matrix
for(unique_column in colnames(corr_signif)){
    for(unique_row in rownames(corr_signif)){
        if(total_results[total_results$Topic==unique_column & total_results$group1==unique_row,"t.test_fdr"]<=0.05){
            if(total_results[total_results$Topic==unique_column & total_results$group1==unique_row,"log2FC_median"]>=0.5){
                corr_signif[unique_row,unique_column] <- "\u2217"
            }
        }
    }
}



############################################################
# Plot heatmap

colnames(matrix_median) <- gsub("Topic_","T",colnames(matrix_median))

p <- pheatmap(matrix_median,
              clustering_distance_rows = "euclidean",
              clustering_distance_cols = "euclidean",
              clustering_method = "ward",
              scale = "none",
              color = colorRampPalette(c(ArchR::ArchRPalettes$greenBlue))(1200),
              breaks = seq(to = max(matrix_median),
                           from = min(matrix_median),
                           length.out = 1200),
              display_numbers = corr_signif,
              fontsize_number = 8,
              treeheight_row = 4,
              treeheight_col = 10,
              cellwidth = 7,
              cellheight = 7,
              fontsize = 6,
              number_color = "white",
              border_color = NA,
              angle = 90
)


col_dend <- p[[2]]
col_dend <- rotate(col_dend, order = paste(
    "T",
    c(2,38,36,6,40,41,31,7,25,39,4,32,42,21,33,1,10,8,30,12,17,20,11,35,14,22,
      13,9,24,
      16,19,34,
      28,5,43,26,
      15,23,29,37,
      3,27,18),
    sep = ""))


cairo_pdf(file = "~/immune_MarkerTopics.pdf", width = 5.5, height = 1.43,
          family = 'Arial Unicode MS')

pheatmap(matrix_median,
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "ward",
         scale = "none",
         color = colorRampPalette(c(ArchR::ArchRPalettes$greenBlue))(1200),
         breaks = seq(to = max(matrix_median),
                      from = min(matrix_median),
                      length.out = 1200),
         display_numbers = corr_signif,
         fontsize_number = 8,
         treeheight_row = 4,
         treeheight_col = 20,
         cellwidth = 7,
         cellheight = 7,
         fontsize = 6,
         number_color = "white",
         border_color = NA,
         angle = 90,
         cluster_cols = as.hclust(col_dend),
)
dev.off()