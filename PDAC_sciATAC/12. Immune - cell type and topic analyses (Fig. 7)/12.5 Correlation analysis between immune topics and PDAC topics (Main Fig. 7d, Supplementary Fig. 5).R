############################################################
# CISTOPIC CORRELATION FUNCTIONS
############################################################

# WORKFLOW:
#   1. Get cisTopic topic matrix -> rearrange as cols = topics, rows = barcodes
#   2. Get annot file
#   3. cisTopicSampleMetrics() -> correlateTopics() -> pheatmap



# Prepare topic medians/means (and ranked versions) by sample data frame for correlateTopics()
cisTopicSampleMetrics <- function(topic_matrix,
                                  annot_sample,  # col1 = cell, col2 = sample ID
                                  cellsToKeep   # list of cells to limit output to (optional)
){
    
    library(dplyr)
    
    # Rename annot_sample columns
    colnames(annot_sample) <- c("cell","sample")
    
    # Check that annot and topic matrix have barcode overlap
    if("TRUE" %in% names(table(annot_sample$cell %in% rownames(topic_matrix)))==F){
        stop("No barcodes in annot_sample exist in topic_matrix")
    }
    
    # Add sample info
    topic_matrix$sample <- annot_sample$sample[match(rownames(topic_matrix),annot_sample$cell)]
    
    # Subset topic_matrix and annot_sample by cellsToKeep (optional)
    if(missing("cellsToKeep")==F){
        if(!(class(cellsToKeep)=="character")){
            stop("cellsToKeep must be a list of barcodes")
        }
        if("TRUE" %in% names(table(cellsToKeep %in% rownames(topic_matrix)))==F){
            stop("No barcodes in cellsToKeep exist in the topic matrix")
        }
        annot_sample <- annot_sample[annot_sample$cell %in% cellsToKeep,,drop=F]
        topic_matrix <- topic_matrix[rownames(topic_matrix) %in% annot_sample$cell,,drop=F]
    }
    
    # Get counts per sample
    sample_counts <- topic_matrix %>%
        group_by(sample) %>%
        summarize(counts = n()) %>%
        as.data.frame()
    
    topic_median <- topic_matrix %>%
        group_by(sample) %>%
        summarize_all(median) %>%
        as.data.frame()
    colnames(topic_median) <- paste(colnames(topic_median),"Median",sep="_")
    colnames(topic_median) <- gsub("sample_Median","sample",colnames(topic_median))
    rownames(topic_median) <- topic_median$sample
    topic_median <- topic_median[-1]
    
    topic_mean <- topic_matrix %>%
        group_by(sample) %>%
        summarize_all(mean) %>%
        as.data.frame()
    colnames(topic_mean) <- paste(colnames(topic_mean),"Mean",sep="_")
    colnames(topic_mean) <- gsub("sample_Mean","sample",colnames(topic_mean))
    rownames(topic_mean) <- topic_mean$sample
    topic_mean <- topic_mean[-1]
    
    
    
    #############################
    # Topic rank mean/median
    topic_rank <- topic_matrix
    for(unique_col in colnames(topic_rank)){
        if(unique_col!="sample"){
            topic_rank[[unique_col]] <- rank(topic_rank[[unique_col]])
        }
        rm(unique_col)
    }
    
    topic_rank_mean <- topic_rank %>%
        group_by(sample) %>%
        summarize_all(mean) %>%
        as.data.frame()
    colnames(topic_rank_mean) <- paste(colnames(topic_rank_mean),"Rank_Mean",sep="_")
    colnames(topic_rank_mean) <- gsub("sample_Rank_Mean","sample",colnames(topic_rank_mean))
    
    rownames(topic_rank_mean) <- topic_rank_mean$sample
    topic_rank_mean <- topic_rank_mean[-1]
    
    topic_rank_median <- topic_rank %>%
        group_by(sample) %>%
        summarize_all(median) %>%
        as.data.frame()
    colnames(topic_rank_median) <- paste(colnames(topic_rank_median),"Rank_Median",sep="_")
    colnames(topic_rank_median) <- gsub("sample_Rank_Median","sample",colnames(topic_rank_median))
    
    rownames(topic_rank_median) <- topic_rank_median$sample
    topic_rank_median <- topic_rank_median[-1]
    
    
    
    return(cbind(sample_counts,topic_mean,topic_median,topic_rank_mean,topic_rank_median))
}





correlateTopics <- function(cistopic1,
                            cistopic2,
                            cistopic1_name = "Epithelial",
                            cistopic2_name = "Immune",
                            metric = "Median",
                            min_cells = 25,
                            returnCorrelationMatrix = T,
                            cor_method = "spearman"){
    
    # Filter matrices
    cistopic1 <- cistopic1[cistopic1$counts>=min_cells,]
    cistopic2 <- cistopic2[cistopic2$counts>=min_cells,]
    
    # Remove samples that are not in both matrices
    cistopic1 <- cistopic1[cistopic1$sample %in% cistopic2$sample,]
    cistopic2 <- cistopic2[cistopic2$sample %in% cistopic1$sample,]
    
    # Rename column names
    colnames(cistopic1)[-c(1,2)] = paste(colnames(cistopic1)[-c(1,2)], cistopic1_name, sep="_")
    colnames(cistopic2)[-c(1,2)] = paste(colnames(cistopic2)[-c(1,2)], cistopic2_name, sep="_")
    
    # Ensure sample order is identical
    cistopic1 <- cistopic1[order(cistopic1$sample),]
    cistopic2 <- cistopic2[order(cistopic2$sample),]
    
    if(!(identical(cistopic1$sample,cistopic2$sample))){
        stop("Sample order of matrices not identical")
    }
    
    if(metric %in% c("Median","median","mean","Mean")==F){
        stop("Acceptable metrics are mean or median")
    }
    
    if(metric %in% c("Median","median")){
        cistopic_corr <- cbind(cistopic1[grep("Median",colnames(cistopic1))],cistopic2[grep("Median",colnames(cistopic2))])
        cistopic_corr <- cistopic_corr[grep("Rank",rownames(cistopic_corr),invert = T),]
        cistopic_corr <- cistopic_corr[,grep("Rank",colnames(cistopic_corr),invert = T)]
    } else{
        if(metric %in% c("Mean","mean")){
            cistopic_corr <- cbind(cistopic1[grep("Mean",colnames(cistopic1))],cistopic2[grep("Mean",colnames(cistopic2))])
            cistopic_corr <- cistopic_corr[grep("Rank",rownames(cistopic_corr),invert = T),]
            cistopic_corr <- cistopic_corr[,grep("Rank",colnames(cistopic_corr),invert = T)]
        } else{
            stop("Acceptable metrics are mean or median")
        }
    }
    
    if(returnCorrelationMatrix){
        cistopic_corr <- as.data.frame(cor(cistopic_corr,method = cor_method))
        cistopic_corr <- cistopic_corr[grep(cistopic1_name,rownames(cistopic_corr)),grep(cistopic2_name,colnames(cistopic_corr))]
    }
    return(cistopic_corr)
}



########################################################################################################################
########################################################################################################################
########################################################################################################################
# CORRELATE TOPIC DATA
########################################################################################################################
########################################################################################################################
########################################################################################################################

############################################################
# Topic matrices

# Topic matrix
cistopic_imm <- as.data.frame(t(read.table("~/cisTopic/Primary_TSSMixed_nFrags1k_FRIS0.4/peakset_MACS2ArchR_PrimaryTSSMixnFrags1kFRIS0.4/Immune/cisTopic_50/Perplexity_topics/cisTopic_alpha50_ImmuneTSSMixnFrags1kFRIS0.4_MACS2ArchR_Topics43_Matrix_Scaled.matrix")))
rownames(cistopic_imm) <- paste("sears_tumors_hg38",rownames(cistopic_imm), sep = "#")

cistopic_epi <- as.data.frame(t(read.table("~/cisTopic/AllSamples_TSSMixed_nFrags1k_FRIS0.4/peakset_MACS3_AllSamplesTSSMixednFrags1kFRIS0.4/Epithelial/cisTopic100_Alpha50/Perplexity_topics/cisTopic100Alpha50_EpithelialAllSamplesTSSMixnFrags1kFRIS0.4_AllSamplesMACS3_Topics50_Matrix.matrix")))
rownames(cistopic_epi) <- paste("sears_tumors_hg38",rownames(cistopic_epi), sep = "#")



############################################################
library(ArchR)
set.seed(1)
addArchRGenome("hg38")
library(parallel)
library(dplyr)

library(pheatmap)



# Load Immune ArchRProject
proj_imm <- loadArchRProject("~/ArchR/Primary_TSSMixed_nFrags1k_FRIS0.4/Immune_MACS2ArchR_cisTopic43")

# Load Epithelial ArchRProject
proj_epi <- loadArchRProject("~/ArchR/Epithelial_AllSamples_TSSMixed_nFrags1k_FRIS0.4_PeaksMACS3_cisTopic50")




########################################################################################################################
# Select cells from primary PDAC samples (excluding ampulla of vater samples) with at least 10 T cells and 10 malignant cells

# Minimum cell count cutoff
min_cells_imm <- 10
min_cells_epi <- 10
celltype_imm <- "T cell"

count_imm <- as.data.frame(table(proj_imm$Sample_ID[proj_imm$celltype %in% celltype_imm]))
count_epi <- as.data.frame(table(proj_epi$Sample_ID[grep("Basal|Classical|Intermediate|Malignant PAV",proj_epi$celltype)]))
count_epi <- count_epi[grep("-T",count_epi$Var1),]

pass_imm <- proj_imm$cellNames[proj_imm$Sample_ID %in% count_imm$Var1[count_imm$Freq>=min_cells_imm]]
pass_imm <- pass_imm[pass_imm %in% proj_imm$cellNames[proj_imm$celltype %in% celltype_imm]]
pass_imm <- pass_imm[pass_imm %in% proj_imm$cellNames[proj_imm$Sample_ID %in% c("ST-00017310-T","ST-00019171-T")==F]]

pass_epi <- proj_epi$cellNames[proj_epi$Sample_ID %in% count_epi$Var1[count_epi$Freq>=min_cells_epi]]
pass_epi <- pass_epi[pass_epi %in% proj_epi$cellNames[grep("Basal|Classical|Intermediate|Malignant PAV",proj_epi$celltype)]]
pass_epi <- pass_epi[pass_epi %in% proj_epi$cellNames[proj_epi$Sample_ID %in% c("ST-00017310-T","ST-00019171-T")==F]]



immune_topics <- list(
    "T cell" = c(5,15,16,19,23,24,26,34,42,43),
    "B cell" = c(9,11,12,14,22,24,31),
    "Plasma B" = c(3,4,12,14,17,20,22,27,35),
    "Myeloid" = c(6,7,8,10,21,22,25,28,30,31,35,36,39),
    "Mast" = c(1,2,4,7,8,10,21,28,30,32,33,41,42),
    "Unknown" = c(1,2,6,9,13,21,30,33,36,38,40,41)
)
topics_to_keep_immune <- immune_topics[[celltype_imm]]
topics_to_keep_immune <- paste("Topic_",topics_to_keep_immune,"_Mean_Immune",sep="")


# Epigenetic subtype associated topics
topics_to_keep_epithelial <- c(9,22,23,47,
                               6,20,42,
                               3,10,18,29,30,34,37,40,41,44,45,46)
topics_to_keep_epithelial <- paste("Topic_",topics_to_keep_epithelial,"_Mean_Epithelial",sep="")





# Epithelial vs. Immune
corr_pheatmap <- correlateTopics(
    
    cistopic1_name = "Immune",
    cistopic1 = cisTopicSampleMetrics(topic_matrix = cistopic_imm[rownames(cistopic_imm) %in% pass_imm,],
                                      annot_sample = data.frame(cell = proj_imm$cellNames[proj_imm$cellNames %in% pass_imm],
                                                                sample = proj_imm$Sample_ID[proj_imm$cellNames %in% pass_imm],
                                                                cluster = proj_imm$celltype[proj_imm$cellNames %in% pass_imm])
    ),
    
    cistopic2_name = "Epithelial",
    cistopic2 = cisTopicSampleMetrics(topic_matrix = cistopic_epi[rownames(cistopic_epi) %in% pass_epi,],
                                      annot_sample = data.frame(cell = proj_epi$cellNames[proj_epi$cellNames %in% pass_epi],
                                                                sample = proj_epi$Sample_ID[proj_epi$cellNames %in% pass_epi],
                                                                cluster = proj_epi$celltype[proj_epi$cellNames %in% pass_epi])
    ),
    
    metric = "Mean",
    min_cells = 0,       # Set to zero since we manually set cutoffs
    returnCorrelationMatrix = TRUE,
    cor_method = "kendall"
)



# Epithelial vs. Immune, return DF for significance calculation
corr_signif <- correlateTopics(
    
    cistopic1_name = "Immune",
    cistopic1 = cisTopicSampleMetrics(topic_matrix = cistopic_imm[rownames(cistopic_imm) %in% pass_imm,],
                                      annot_sample = data.frame(cell = proj_imm$cellNames[proj_imm$cellNames %in% pass_imm],
                                                                sample = proj_imm$Sample_ID[proj_imm$cellNames %in% pass_imm],
                                                                cluster = proj_imm$celltype[proj_imm$cellNames %in% pass_imm])
    ),
    
    cistopic2_name = "Epithelial",
    cistopic2 = cisTopicSampleMetrics(topic_matrix = cistopic_epi[rownames(cistopic_epi) %in% pass_epi,],
                                      annot_sample = data.frame(cell = proj_epi$cellNames[proj_epi$cellNames %in% pass_epi],
                                                                sample = proj_epi$Sample_ID[proj_epi$cellNames %in% pass_epi],
                                                                cluster = proj_epi$celltype[proj_epi$cellNames %in% pass_epi])
    ),
    
    metric = "Mean",
    min_cells = 0,       # Set to zero since we manually set cutoffs
    returnCorrelationMatrix = FALSE,
    cor_method = "kendall"
)


corr_samples <- corr_signif
corr_signif <- corr_signif[,c(topics_to_keep_immune,topics_to_keep_epithelial),drop=F]
corr_signif <- corrplot::cor.mtest(corr_signif, conf.level = 0.95, method = "kendall")$p
corr_signif <- corr_signif[topics_to_keep_immune,topics_to_keep_epithelial]



# Extract all comparisons and their unique p-values into column format
corr_signif_padj <- data.frame()
i = 1
for(unique_row in rownames(corr_signif)){
    for(unique_col in colnames(corr_signif)){
        corr_signif_padj[i,"row"] <- unique_row
        corr_signif_padj[i,"col"] <- unique_col
        corr_signif_padj[i,"pval"] <- corr_signif[unique_row,unique_col]
        i = i+1
    }
}
rm(unique_row,unique_col,i)



# Calculate corrected p-values
corr_signif_padj$padj <- p.adjust(corr_signif_padj$pval, method = "fdr")



# Apply corrected p-values to corr_signif matrix
for(i in 1:length(rownames(corr_signif_padj))){
    temp_rowname <- corr_signif_padj[i,"row"]
    temp_colname <- corr_signif_padj[i,"col"]
    temp_padj <- corr_signif_padj[i,"padj"]
    corr_signif[temp_rowname, temp_colname] <- temp_padj
    rm(temp_rowname, temp_colname, temp_padj)
}


corr_signif_threshold <- 0.05
corr_signif[corr_signif>corr_signif_threshold] <- ""
corr_signif[corr_signif!=""] <- "\u2217" # Unicode centered asterisk





############################################################
# Annotate epi. topics by subtype

classicalNRP_topics <- c(3,10,18,29,30,34,37,40,41,44,45,46)
classicalNRP_topics <- paste("Topic_", classicalNRP_topics, "_Mean_Epithelial", sep = "")

basal_topics <- c(9,22,23,47)
basal_topics <- paste("Topic_", basal_topics, "_Mean_Epithelial", sep = "")

intermed_topics <- c(6,20,42)
intermed_topics <- paste("Topic_", intermed_topics, "_Mean_Epithelial", sep = "")


heatmap_annot <- data.frame(topic = colnames(corr_pheatmap[topics_to_keep_immune,topics_to_keep_epithelial]))
heatmap_annot$Subtype[heatmap_annot$topic %in% classicalNRP_topics] <- "Classical-NRP"
heatmap_annot$Subtype[heatmap_annot$topic %in% basal_topics] <- "Basal Committed"
heatmap_annot$Subtype[heatmap_annot$topic %in% intermed_topics] <- "Subtype-Intermediate Progenitor"
rownames(heatmap_annot) <- heatmap_annot$topic
heatmap_annot <- heatmap_annot[,-1,drop=F]







# Heatmap colors
heatmap_colors <- list(
    "Subtype" = c("Basal Committed" = "#bd3531",
                  "Classical-NRP" = "#f79747",
                  "Subtype-Intermediate Progenitor" = "#e1468a"))



library(pheatmap)
# For clade rotating
library(seriation)
library(dendextend)



p <- pheatmap(t(corr_pheatmap[topics_to_keep_immune,topics_to_keep_epithelial]),
              clustering_distance_rows = "euclidean",
              clustering_distance_cols = "euclidean",
              clustering_method = "ward",
              
              annotation_row = heatmap_annot,
              annotation_colors = heatmap_colors,
              
              # scale = "none",
              
              fontsize = 6,
              angle_col = 90,
              
              color = colorRampPalette(ArchR::ArchRPalettes$solarExtra)(100),
              breaks = seq(from = -max(abs(corr_pheatmap)),
                           to = max(abs(corr_pheatmap)),
                           length.out = 100),
              annotation_legend = T,
              
              treeheight_row = 15,
              treeheight_col = 15,
              
              display_numbers = t(corr_signif),
              fontsize_number = 8,
              number_color = "gray95",
              
              border_color = NA
)





col_dend <- p[[2]]
col_dend <- rotate(col_dend, order = paste(
    "Topic_",
    c(15,23,5,16,34,
      19,
      24,43,26,42
    ),
    "_Mean_Immune", sep = ""))






cairo_pdf(width = 6.5, height = 5,
          file = "~/Topics_Tcell_vs_PDAC.pdf",
          family = 'Arial Unicode MS')

ComplexHeatmap::pheatmap(t(corr_pheatmap[topics_to_keep_immune,topics_to_keep_epithelial]),
                         clustering_distance_rows = "euclidean",
                         clustering_distance_cols = "euclidean",
                         clustering_method = "ward",
                         
                         annotation_row = heatmap_annot,
                         annotation_colors = heatmap_colors,
                         
                         # scale = "none",
                         
                         cellwidth = 8,
                         cellheight = 8,
                         
                         fontsize = 6,
                         angle_col = "90",
                         
                         color = colorRampPalette(ArchR::ArchRPalettes$solarExtra)(1200),
                         breaks = seq(from = -max(abs(corr_pheatmap)),
                                      to = max(abs(corr_pheatmap)),
                                      length.out = 1200),
                         annotation_legend = T,
                         
                         labels_row = paste("T",gsub("^Topic_|_Mean_.*", "", rownames(t(corr_pheatmap[topics_to_keep_immune,topics_to_keep_epithelial]))),sep=""),
                         labels_col = paste("T",gsub("^Topic_|_Mean_.*", "", colnames(t(corr_pheatmap[topics_to_keep_immune,topics_to_keep_epithelial]))),sep=""),
                         
                         
                         treeheight_row = 15,
                         treeheight_col = 15,
                         
                         display_numbers = t(corr_signif),
                         fontsize_number = 8,
                         number_color = "gray95",
                         cluster_cols = as.hclust(col_dend),
                         
                         border_color = NA,
                         row_dend_side ="right",
                         row_names_side = "left",
                         heatmap_legend_param = list(
                             title = "Value",
                             at = c(-max(abs(corr_pheatmap)), max(abs(corr_pheatmap)))
                         )
)

dev.off()






########################################################################################################################
########################################################################################################################
########################################################################################################################
# CORRELATION HEATMAP: B-CELL
########################################################################################################################
########################################################################################################################
########################################################################################################################

celltype_imm <- "B cell"

count_imm <- as.data.frame(table(proj_imm$Sample_ID[proj_imm$celltype %in% celltype_imm]))
count_epi <- as.data.frame(table(proj_epi$Sample_ID[grep("Basal|Classical|Intermediate|Malignant PAV",proj_epi$celltype)]))
count_epi <- count_epi[grep("-T",count_epi$Var1),]

pass_imm <- proj_imm$cellNames[proj_imm$Sample_ID %in% count_imm$Var1[count_imm$Freq>=min_cells_imm]]
pass_imm <- pass_imm[pass_imm %in% proj_imm$cellNames[proj_imm$celltype %in% celltype_imm]]
pass_imm <- pass_imm[pass_imm %in% proj_imm$cellNames[proj_imm$Sample_ID %in% c("ST-00017310-T","ST-00019171-T")==F]]

pass_epi <- proj_epi$cellNames[proj_epi$Sample_ID %in% count_epi$Var1[count_epi$Freq>=min_cells_epi]]
pass_epi <- pass_epi[pass_epi %in% proj_epi$cellNames[grep("Basal|Classical|Intermediate|Malignant PAV",proj_epi$celltype)]]
pass_epi <- pass_epi[pass_epi %in% proj_epi$cellNames[proj_epi$Sample_ID %in% c("ST-00017310-T","ST-00019171-T")==F]]




topics_to_keep_immune <- immune_topics[[celltype_imm]]
topics_to_keep_immune <- paste("Topic_",topics_to_keep_immune,"_Mean_Immune",sep="")





# Epithelial vs. Immune
corr_pheatmap <- correlateTopics(
    
    cistopic1_name = "Immune",
    cistopic1 = cisTopicSampleMetrics(topic_matrix = cistopic_imm[rownames(cistopic_imm) %in% pass_imm,],
                                      annot_sample = data.frame(cell = proj_imm$cellNames[proj_imm$cellNames %in% pass_imm],
                                                                sample = proj_imm$Sample_ID[proj_imm$cellNames %in% pass_imm],
                                                                cluster = proj_imm$celltype[proj_imm$cellNames %in% pass_imm])
    ),
    
    cistopic2_name = "Epithelial",
    cistopic2 = cisTopicSampleMetrics(topic_matrix = cistopic_epi[rownames(cistopic_epi) %in% pass_epi,],
                                      annot_sample = data.frame(cell = proj_epi$cellNames[proj_epi$cellNames %in% pass_epi],
                                                                sample = proj_epi$Sample_ID[proj_epi$cellNames %in% pass_epi],
                                                                cluster = proj_epi$celltype[proj_epi$cellNames %in% pass_epi])
    ),
    
    metric = "Mean",
    min_cells = 0,       # Set to zero since we manually set cutoffs
    returnCorrelationMatrix = TRUE,
    cor_method = "kendall"
)



# Epithelial vs. Immune, return DF for significance calculation
corr_signif <- correlateTopics(
    
    cistopic1_name = "Immune",
    cistopic1 = cisTopicSampleMetrics(topic_matrix = cistopic_imm[rownames(cistopic_imm) %in% pass_imm,],
                                      annot_sample = data.frame(cell = proj_imm$cellNames[proj_imm$cellNames %in% pass_imm],
                                                                sample = proj_imm$Sample_ID[proj_imm$cellNames %in% pass_imm],
                                                                cluster = proj_imm$celltype[proj_imm$cellNames %in% pass_imm])
    ),
    
    cistopic2_name = "Epithelial",
    cistopic2 = cisTopicSampleMetrics(topic_matrix = cistopic_epi[rownames(cistopic_epi) %in% pass_epi,],
                                      annot_sample = data.frame(cell = proj_epi$cellNames[proj_epi$cellNames %in% pass_epi],
                                                                sample = proj_epi$Sample_ID[proj_epi$cellNames %in% pass_epi],
                                                                cluster = proj_epi$celltype[proj_epi$cellNames %in% pass_epi])
    ),
    
    metric = "Mean",
    min_cells = 0,       # Set to zero since we manually set cutoffs
    returnCorrelationMatrix = FALSE,
    cor_method = "kendall"
)


corr_samples <- corr_signif
corr_signif <- corr_signif[,c(topics_to_keep_immune,topics_to_keep_epithelial),drop=F]
corr_signif <- corrplot::cor.mtest(corr_signif, conf.level = 0.95, method = "kendall")$p
corr_signif <- corr_signif[topics_to_keep_immune,topics_to_keep_epithelial]



# Extract all comparisons and their unique p-values into column format
corr_signif_padj <- data.frame()
i = 1
for(unique_row in rownames(corr_signif)){
    for(unique_col in colnames(corr_signif)){
        corr_signif_padj[i,"row"] <- unique_row
        corr_signif_padj[i,"col"] <- unique_col
        corr_signif_padj[i,"pval"] <- corr_signif[unique_row,unique_col]
        i = i+1
    }
}
rm(unique_row,unique_col,i)



# Calculate corrected p-values
corr_signif_padj$padj <- p.adjust(corr_signif_padj$pval, method = "fdr")



# Apply corrected p-values to corr_signif matrix
for(i in 1:length(rownames(corr_signif_padj))){
    temp_rowname <- corr_signif_padj[i,"row"]
    temp_colname <- corr_signif_padj[i,"col"]
    temp_padj <- corr_signif_padj[i,"padj"]
    corr_signif[temp_rowname, temp_colname] <- temp_padj
    rm(temp_rowname, temp_colname, temp_padj)
}


corr_signif[corr_signif>corr_signif_threshold] <- ""
corr_signif[corr_signif!=""] <- "\u2217" # Unicode centered asterisk





############################################################
# Annotate epi. topics by subtype

heatmap_annot <- data.frame(topic = colnames(corr_pheatmap[topics_to_keep_immune,topics_to_keep_epithelial]))
heatmap_annot$Subtype[heatmap_annot$topic %in% classicalNRP_topics] <- "Classical-NRP"
heatmap_annot$Subtype[heatmap_annot$topic %in% basal_topics] <- "Basal Committed"
heatmap_annot$Subtype[heatmap_annot$topic %in% intermed_topics] <- "Subtype-Intermediate Progenitor"
rownames(heatmap_annot) <- heatmap_annot$topic
heatmap_annot <- heatmap_annot[,-1,drop=F]



cairo_pdf(width = 6.5, height = 5,
          file = "~/Topics_Bcell_vs_PDAC.pdf",
          family = 'Arial Unicode MS')

ComplexHeatmap::pheatmap(t(corr_pheatmap[topics_to_keep_immune,topics_to_keep_epithelial]),
                         clustering_distance_rows = "euclidean",
                         clustering_distance_cols = "euclidean",
                         clustering_method = "ward",
                         
                         annotation_row = heatmap_annot,
                         annotation_colors = heatmap_colors,
                         
                         # scale = "none",
                         
                         cellwidth = 8,
                         cellheight = 8,
                         
                         fontsize = 6,
                         angle_col = "90",
                         
                         color = colorRampPalette(ArchR::ArchRPalettes$solarExtra)(1200),
                         breaks = seq(from = -max(abs(corr_pheatmap)),
                                      to = max(abs(corr_pheatmap)),
                                      length.out = 1200),
                         annotation_legend = T,
                         
                         labels_row = paste("T",gsub("^Topic_|_Mean_.*", "", rownames(t(corr_pheatmap[topics_to_keep_immune,topics_to_keep_epithelial]))),sep=""),
                         labels_col = paste("T",gsub("^Topic_|_Mean_.*", "", colnames(t(corr_pheatmap[topics_to_keep_immune,topics_to_keep_epithelial]))),sep=""),
                         
                         
                         treeheight_row = 15,
                         treeheight_col = 15,
                         
                         display_numbers = t(corr_signif),
                         fontsize_number = 8,
                         number_color = "gray95",
                         
                         border_color = NA,
                         row_dend_side ="right",
                         row_names_side = "left",
                         heatmap_legend_param = list(
                             title = "Value",
                             at = c(-max(abs(corr_pheatmap)), max(abs(corr_pheatmap)))
                         )
)
dev.off()








########################################################################################################################
########################################################################################################################
########################################################################################################################
# CORRELATION HEATMAP: MYELOID
########################################################################################################################
########################################################################################################################
########################################################################################################################

celltype_imm <- "Myeloid"

count_imm <- as.data.frame(table(proj_imm$Sample_ID[proj_imm$celltype %in% celltype_imm]))
count_epi <- as.data.frame(table(proj_epi$Sample_ID[grep("Basal|Classical|Intermediate|Malignant PAV",proj_epi$celltype)]))
count_epi <- count_epi[grep("-T",count_epi$Var1),]

pass_imm <- proj_imm$cellNames[proj_imm$Sample_ID %in% count_imm$Var1[count_imm$Freq>=min_cells_imm]]
pass_imm <- pass_imm[pass_imm %in% proj_imm$cellNames[proj_imm$celltype %in% celltype_imm]]
pass_imm <- pass_imm[pass_imm %in% proj_imm$cellNames[proj_imm$Sample_ID %in% c("ST-00017310-T","ST-00019171-T")==F]]

pass_epi <- proj_epi$cellNames[proj_epi$Sample_ID %in% count_epi$Var1[count_epi$Freq>=min_cells_epi]]
pass_epi <- pass_epi[pass_epi %in% proj_epi$cellNames[grep("Basal|Classical|Intermediate|Malignant PAV",proj_epi$celltype)]]
pass_epi <- pass_epi[pass_epi %in% proj_epi$cellNames[proj_epi$Sample_ID %in% c("ST-00017310-T","ST-00019171-T")==F]]




topics_to_keep_immune <- immune_topics[[celltype_imm]]
topics_to_keep_immune <- paste("Topic_",topics_to_keep_immune,"_Mean_Immune",sep="")





# Epithelial vs. Immune
corr_pheatmap <- correlateTopics(
    
    cistopic1_name = "Immune",
    cistopic1 = cisTopicSampleMetrics(topic_matrix = cistopic_imm[rownames(cistopic_imm) %in% pass_imm,],
                                      annot_sample = data.frame(cell = proj_imm$cellNames[proj_imm$cellNames %in% pass_imm],
                                                                sample = proj_imm$Sample_ID[proj_imm$cellNames %in% pass_imm],
                                                                cluster = proj_imm$celltype[proj_imm$cellNames %in% pass_imm])
    ),
    
    cistopic2_name = "Epithelial",
    cistopic2 = cisTopicSampleMetrics(topic_matrix = cistopic_epi[rownames(cistopic_epi) %in% pass_epi,],
                                      annot_sample = data.frame(cell = proj_epi$cellNames[proj_epi$cellNames %in% pass_epi],
                                                                sample = proj_epi$Sample_ID[proj_epi$cellNames %in% pass_epi],
                                                                cluster = proj_epi$celltype[proj_epi$cellNames %in% pass_epi])
    ),
    
    metric = "Mean",
    min_cells = 0,       # Set to zero since we manually set cutoffs
    returnCorrelationMatrix = TRUE,
    cor_method = "kendall"
)



# Epithelial vs. Immune, return DF for significance calculation
corr_signif <- correlateTopics(
    
    cistopic1_name = "Immune",
    cistopic1 = cisTopicSampleMetrics(topic_matrix = cistopic_imm[rownames(cistopic_imm) %in% pass_imm,],
                                      annot_sample = data.frame(cell = proj_imm$cellNames[proj_imm$cellNames %in% pass_imm],
                                                                sample = proj_imm$Sample_ID[proj_imm$cellNames %in% pass_imm],
                                                                cluster = proj_imm$celltype[proj_imm$cellNames %in% pass_imm])
    ),
    
    cistopic2_name = "Epithelial",
    cistopic2 = cisTopicSampleMetrics(topic_matrix = cistopic_epi[rownames(cistopic_epi) %in% pass_epi,],
                                      annot_sample = data.frame(cell = proj_epi$cellNames[proj_epi$cellNames %in% pass_epi],
                                                                sample = proj_epi$Sample_ID[proj_epi$cellNames %in% pass_epi],
                                                                cluster = proj_epi$celltype[proj_epi$cellNames %in% pass_epi])
    ),
    
    metric = "Mean",
    min_cells = 0,       # Set to zero since we manually set cutoffs
    returnCorrelationMatrix = FALSE,
    cor_method = "kendall"
)


corr_samples <- corr_signif
corr_signif <- corr_signif[,c(topics_to_keep_immune,topics_to_keep_epithelial),drop=F]
corr_signif <- corrplot::cor.mtest(corr_signif, conf.level = 0.95, method = "kendall")$p
corr_signif <- corr_signif[topics_to_keep_immune,topics_to_keep_epithelial]



# Extract all comparisons and their unique p-values into column format
corr_signif_padj <- data.frame()
i = 1
for(unique_row in rownames(corr_signif)){
    for(unique_col in colnames(corr_signif)){
        corr_signif_padj[i,"row"] <- unique_row
        corr_signif_padj[i,"col"] <- unique_col
        corr_signif_padj[i,"pval"] <- corr_signif[unique_row,unique_col]
        i = i+1
    }
}
rm(unique_row,unique_col,i)



# Calculate corrected p-values
corr_signif_padj$padj <- p.adjust(corr_signif_padj$pval, method = "fdr")



# Apply corrected p-values to corr_signif matrix
for(i in 1:length(rownames(corr_signif_padj))){
    temp_rowname <- corr_signif_padj[i,"row"]
    temp_colname <- corr_signif_padj[i,"col"]
    temp_padj <- corr_signif_padj[i,"padj"]
    corr_signif[temp_rowname, temp_colname] <- temp_padj
    rm(temp_rowname, temp_colname, temp_padj)
}


corr_signif[corr_signif>corr_signif_threshold] <- ""
corr_signif[corr_signif!=""] <- "\u2217" # Unicode centered asterisk





############################################################
# Annotate epi. topics by subtype

heatmap_annot <- data.frame(topic = colnames(corr_pheatmap[topics_to_keep_immune,topics_to_keep_epithelial]))
heatmap_annot$Subtype[heatmap_annot$topic %in% classicalNRP_topics] <- "Classical-NRP"
heatmap_annot$Subtype[heatmap_annot$topic %in% basal_topics] <- "Basal Committed"
heatmap_annot$Subtype[heatmap_annot$topic %in% intermed_topics] <- "Subtype-Intermediate Progenitor"
rownames(heatmap_annot) <- heatmap_annot$topic
heatmap_annot <- heatmap_annot[,-1,drop=F]



cairo_pdf(width = 6.5, height = 5,
          file = "~/Topics_Myeloid_vs_PDAC.pdf",
          family = 'Arial Unicode MS')

ComplexHeatmap::pheatmap(t(corr_pheatmap[topics_to_keep_immune,topics_to_keep_epithelial]),
                         clustering_distance_rows = "euclidean",
                         clustering_distance_cols = "euclidean",
                         clustering_method = "ward",
                         
                         annotation_row = heatmap_annot,
                         annotation_colors = heatmap_colors,
                         
                         # scale = "none",
                         
                         cellwidth = 8,
                         cellheight = 8,
                         
                         fontsize = 6,
                         angle_col = "90",
                         
                         color = colorRampPalette(ArchR::ArchRPalettes$solarExtra)(1200),
                         breaks = seq(from = -max(abs(corr_pheatmap)),
                                      to = max(abs(corr_pheatmap)),
                                      length.out = 1200),
                         annotation_legend = T,
                         
                         labels_row = paste("T",gsub("^Topic_|_Mean_.*", "", rownames(t(corr_pheatmap[topics_to_keep_immune,topics_to_keep_epithelial]))),sep=""),
                         labels_col = paste("T",gsub("^Topic_|_Mean_.*", "", colnames(t(corr_pheatmap[topics_to_keep_immune,topics_to_keep_epithelial]))),sep=""),
                         
                         
                         treeheight_row = 15,
                         treeheight_col = 15,
                         
                         display_numbers = t(corr_signif),
                         fontsize_number = 8,
                         number_color = "gray95",
                         
                         border_color = NA,
                         row_dend_side ="right",
                         row_names_side = "left",
                         heatmap_legend_param = list(
                             title = "Value",
                             at = c(-max(abs(corr_pheatmap)), max(abs(corr_pheatmap)))
                         )
)
dev.off()









########################################################################################################################
########################################################################################################################
########################################################################################################################
# CORRELATION HEATMAP: MAST
########################################################################################################################
########################################################################################################################
########################################################################################################################

celltype_imm <- "Mast"

count_imm <- as.data.frame(table(proj_imm$Sample_ID[proj_imm$celltype %in% celltype_imm]))
count_epi <- as.data.frame(table(proj_epi$Sample_ID[grep("Basal|Classical|Intermediate|Malignant PAV",proj_epi$celltype)]))
count_epi <- count_epi[grep("-T",count_epi$Var1),]

pass_imm <- proj_imm$cellNames[proj_imm$Sample_ID %in% count_imm$Var1[count_imm$Freq>=min_cells_imm]]
pass_imm <- pass_imm[pass_imm %in% proj_imm$cellNames[proj_imm$celltype %in% celltype_imm]]
pass_imm <- pass_imm[pass_imm %in% proj_imm$cellNames[proj_imm$Sample_ID %in% c("ST-00017310-T","ST-00019171-T")==F]]

pass_epi <- proj_epi$cellNames[proj_epi$Sample_ID %in% count_epi$Var1[count_epi$Freq>=min_cells_epi]]
pass_epi <- pass_epi[pass_epi %in% proj_epi$cellNames[grep("Basal|Classical|Intermediate|Malignant PAV",proj_epi$celltype)]]
pass_epi <- pass_epi[pass_epi %in% proj_epi$cellNames[proj_epi$Sample_ID %in% c("ST-00017310-T","ST-00019171-T")==F]]




topics_to_keep_immune <- immune_topics[[celltype_imm]]
topics_to_keep_immune <- paste("Topic_",topics_to_keep_immune,"_Mean_Immune",sep="")





# Epithelial vs. Immune
corr_pheatmap <- correlateTopics(
    
    cistopic1_name = "Immune",
    cistopic1 = cisTopicSampleMetrics(topic_matrix = cistopic_imm[rownames(cistopic_imm) %in% pass_imm,],
                                      annot_sample = data.frame(cell = proj_imm$cellNames[proj_imm$cellNames %in% pass_imm],
                                                                sample = proj_imm$Sample_ID[proj_imm$cellNames %in% pass_imm],
                                                                cluster = proj_imm$celltype[proj_imm$cellNames %in% pass_imm])
    ),
    
    cistopic2_name = "Epithelial",
    cistopic2 = cisTopicSampleMetrics(topic_matrix = cistopic_epi[rownames(cistopic_epi) %in% pass_epi,],
                                      annot_sample = data.frame(cell = proj_epi$cellNames[proj_epi$cellNames %in% pass_epi],
                                                                sample = proj_epi$Sample_ID[proj_epi$cellNames %in% pass_epi],
                                                                cluster = proj_epi$celltype[proj_epi$cellNames %in% pass_epi])
    ),
    
    metric = "Mean",
    min_cells = 0,       # Set to zero since we manually set cutoffs
    returnCorrelationMatrix = TRUE,
    cor_method = "kendall"
)



# Epithelial vs. Immune, return DF for significance calculation
corr_signif <- correlateTopics(
    
    cistopic1_name = "Immune",
    cistopic1 = cisTopicSampleMetrics(topic_matrix = cistopic_imm[rownames(cistopic_imm) %in% pass_imm,],
                                      annot_sample = data.frame(cell = proj_imm$cellNames[proj_imm$cellNames %in% pass_imm],
                                                                sample = proj_imm$Sample_ID[proj_imm$cellNames %in% pass_imm],
                                                                cluster = proj_imm$celltype[proj_imm$cellNames %in% pass_imm])
    ),
    
    cistopic2_name = "Epithelial",
    cistopic2 = cisTopicSampleMetrics(topic_matrix = cistopic_epi[rownames(cistopic_epi) %in% pass_epi,],
                                      annot_sample = data.frame(cell = proj_epi$cellNames[proj_epi$cellNames %in% pass_epi],
                                                                sample = proj_epi$Sample_ID[proj_epi$cellNames %in% pass_epi],
                                                                cluster = proj_epi$celltype[proj_epi$cellNames %in% pass_epi])
    ),
    
    metric = "Mean",
    min_cells = 0,       # Set to zero since we manually set cutoffs
    returnCorrelationMatrix = FALSE,
    cor_method = "kendall"
)


corr_samples <- corr_signif
corr_signif <- corr_signif[,c(topics_to_keep_immune,topics_to_keep_epithelial),drop=F]
corr_signif <- corrplot::cor.mtest(corr_signif, conf.level = 0.95, method = "kendall")$p
corr_signif <- corr_signif[topics_to_keep_immune,topics_to_keep_epithelial]



# Extract all comparisons and their unique p-values into column format
corr_signif_padj <- data.frame()
i = 1
for(unique_row in rownames(corr_signif)){
    for(unique_col in colnames(corr_signif)){
        corr_signif_padj[i,"row"] <- unique_row
        corr_signif_padj[i,"col"] <- unique_col
        corr_signif_padj[i,"pval"] <- corr_signif[unique_row,unique_col]
        i = i+1
    }
}
rm(unique_row,unique_col,i)



# Calculate corrected p-values
corr_signif_padj$padj <- p.adjust(corr_signif_padj$pval, method = "fdr")



# Apply corrected p-values to corr_signif matrix
for(i in 1:length(rownames(corr_signif_padj))){
    temp_rowname <- corr_signif_padj[i,"row"]
    temp_colname <- corr_signif_padj[i,"col"]
    temp_padj <- corr_signif_padj[i,"padj"]
    corr_signif[temp_rowname, temp_colname] <- temp_padj
    rm(temp_rowname, temp_colname, temp_padj)
}


corr_signif[corr_signif>corr_signif_threshold] <- ""
corr_signif[corr_signif!=""] <- "\u2217" # Unicode centered asterisk





############################################################
# Annotate epi. topics by subtype

heatmap_annot <- data.frame(topic = colnames(corr_pheatmap[topics_to_keep_immune,topics_to_keep_epithelial]))
heatmap_annot$Subtype[heatmap_annot$topic %in% classicalNRP_topics] <- "Classical-NRP"
heatmap_annot$Subtype[heatmap_annot$topic %in% basal_topics] <- "Basal Committed"
heatmap_annot$Subtype[heatmap_annot$topic %in% intermed_topics] <- "Subtype-Intermediate Progenitor"
rownames(heatmap_annot) <- heatmap_annot$topic
heatmap_annot <- heatmap_annot[,-1,drop=F]



cairo_pdf(width = 6.5, height = 5,
          file = "~/Topics_Mast_vs_PDAC.pdf",
          family = 'Arial Unicode MS')

ComplexHeatmap::pheatmap(t(corr_pheatmap[topics_to_keep_immune,topics_to_keep_epithelial]),
                         clustering_distance_rows = "euclidean",
                         clustering_distance_cols = "euclidean",
                         clustering_method = "ward",
                         
                         annotation_row = heatmap_annot,
                         annotation_colors = heatmap_colors,
                         
                         # scale = "none",
                         
                         cellwidth = 8,
                         cellheight = 8,
                         
                         fontsize = 6,
                         angle_col = "90",
                         
                         color = colorRampPalette(ArchR::ArchRPalettes$solarExtra)(1200),
                         breaks = seq(from = -max(abs(corr_pheatmap)),
                                      to = max(abs(corr_pheatmap)),
                                      length.out = 1200),
                         annotation_legend = T,
                         
                         labels_row = paste("T",gsub("^Topic_|_Mean_.*", "", rownames(t(corr_pheatmap[topics_to_keep_immune,topics_to_keep_epithelial]))),sep=""),
                         labels_col = paste("T",gsub("^Topic_|_Mean_.*", "", colnames(t(corr_pheatmap[topics_to_keep_immune,topics_to_keep_epithelial]))),sep=""),
                         
                         
                         treeheight_row = 15,
                         treeheight_col = 15,
                         
                         display_numbers = t(corr_signif),
                         fontsize_number = 8,
                         number_color = "gray95",
                         
                         border_color = NA,
                         row_dend_side ="right",
                         row_names_side = "left",
                         heatmap_legend_param = list(
                             title = "Value",
                             at = c(-max(abs(corr_pheatmap)), max(abs(corr_pheatmap)))
                         )
)
dev.off()











########################################################################################################################
########################################################################################################################
########################################################################################################################
# CORRELATION HEATMAP: PLASMA B
########################################################################################################################
########################################################################################################################
########################################################################################################################


celltype_imm <- "Plasma B"

count_imm <- as.data.frame(table(proj_imm$Sample_ID[proj_imm$celltype %in% celltype_imm]))
count_epi <- as.data.frame(table(proj_epi$Sample_ID[grep("Basal|Classical|Intermediate|Malignant PAV",proj_epi$celltype)]))
count_epi <- count_epi[grep("-T",count_epi$Var1),]

pass_imm <- proj_imm$cellNames[proj_imm$Sample_ID %in% count_imm$Var1[count_imm$Freq>=min_cells_imm]]
pass_imm <- pass_imm[pass_imm %in% proj_imm$cellNames[proj_imm$celltype %in% celltype_imm]]
pass_imm <- pass_imm[pass_imm %in% proj_imm$cellNames[proj_imm$Sample_ID %in% c("ST-00017310-T","ST-00019171-T")==F]]

pass_epi <- proj_epi$cellNames[proj_epi$Sample_ID %in% count_epi$Var1[count_epi$Freq>=min_cells_epi]]
pass_epi <- pass_epi[pass_epi %in% proj_epi$cellNames[grep("Basal|Classical|Intermediate|Malignant PAV",proj_epi$celltype)]]
pass_epi <- pass_epi[pass_epi %in% proj_epi$cellNames[proj_epi$Sample_ID %in% c("ST-00017310-T","ST-00019171-T")==F]]




topics_to_keep_immune <- immune_topics[[celltype_imm]]
topics_to_keep_immune <- paste("Topic_",topics_to_keep_immune,"_Mean_Immune",sep="")





# Epithelial vs. Immune
corr_pheatmap <- correlateTopics(
    
    cistopic1_name = "Immune",
    cistopic1 = cisTopicSampleMetrics(topic_matrix = cistopic_imm[rownames(cistopic_imm) %in% pass_imm,],
                                      annot_sample = data.frame(cell = proj_imm$cellNames[proj_imm$cellNames %in% pass_imm],
                                                                sample = proj_imm$Sample_ID[proj_imm$cellNames %in% pass_imm],
                                                                cluster = proj_imm$celltype[proj_imm$cellNames %in% pass_imm])
    ),
    
    cistopic2_name = "Epithelial",
    cistopic2 = cisTopicSampleMetrics(topic_matrix = cistopic_epi[rownames(cistopic_epi) %in% pass_epi,],
                                      annot_sample = data.frame(cell = proj_epi$cellNames[proj_epi$cellNames %in% pass_epi],
                                                                sample = proj_epi$Sample_ID[proj_epi$cellNames %in% pass_epi],
                                                                cluster = proj_epi$celltype[proj_epi$cellNames %in% pass_epi])
    ),
    
    metric = "Mean",
    min_cells = 0,       # Set to zero since we manually set cutoffs
    returnCorrelationMatrix = TRUE,
    cor_method = "kendall"
)



# Epithelial vs. Immune, return DF for significance calculation
corr_signif <- correlateTopics(
    
    cistopic1_name = "Immune",
    cistopic1 = cisTopicSampleMetrics(topic_matrix = cistopic_imm[rownames(cistopic_imm) %in% pass_imm,],
                                      annot_sample = data.frame(cell = proj_imm$cellNames[proj_imm$cellNames %in% pass_imm],
                                                                sample = proj_imm$Sample_ID[proj_imm$cellNames %in% pass_imm],
                                                                cluster = proj_imm$celltype[proj_imm$cellNames %in% pass_imm])
    ),
    
    cistopic2_name = "Epithelial",
    cistopic2 = cisTopicSampleMetrics(topic_matrix = cistopic_epi[rownames(cistopic_epi) %in% pass_epi,],
                                      annot_sample = data.frame(cell = proj_epi$cellNames[proj_epi$cellNames %in% pass_epi],
                                                                sample = proj_epi$Sample_ID[proj_epi$cellNames %in% pass_epi],
                                                                cluster = proj_epi$celltype[proj_epi$cellNames %in% pass_epi])
    ),
    
    metric = "Mean",
    min_cells = 0,       # Set to zero since we manually set cutoffs
    returnCorrelationMatrix = FALSE,
    cor_method = "kendall"
)


corr_samples <- corr_signif
corr_signif <- corr_signif[,c(topics_to_keep_immune,topics_to_keep_epithelial),drop=F]
corr_signif <- corrplot::cor.mtest(corr_signif, conf.level = 0.95, method = "kendall")$p
corr_signif <- corr_signif[topics_to_keep_immune,topics_to_keep_epithelial]



# Extract all comparisons and their unique p-values into column format
corr_signif_padj <- data.frame()
i = 1
for(unique_row in rownames(corr_signif)){
    for(unique_col in colnames(corr_signif)){
        corr_signif_padj[i,"row"] <- unique_row
        corr_signif_padj[i,"col"] <- unique_col
        corr_signif_padj[i,"pval"] <- corr_signif[unique_row,unique_col]
        i = i+1
    }
}
rm(unique_row,unique_col,i)



# Calculate corrected p-values
corr_signif_padj$padj <- p.adjust(corr_signif_padj$pval, method = "fdr")



# Apply corrected p-values to corr_signif matrix
for(i in 1:length(rownames(corr_signif_padj))){
    temp_rowname <- corr_signif_padj[i,"row"]
    temp_colname <- corr_signif_padj[i,"col"]
    temp_padj <- corr_signif_padj[i,"padj"]
    corr_signif[temp_rowname, temp_colname] <- temp_padj
    rm(temp_rowname, temp_colname, temp_padj)
}


corr_signif[corr_signif>corr_signif_threshold] <- ""
corr_signif[corr_signif!=""] <- "\u2217" # Unicode centered asterisk





############################################################
# Annotate epi. topics by subtype

heatmap_annot <- data.frame(topic = colnames(corr_pheatmap[topics_to_keep_immune,topics_to_keep_epithelial]))
heatmap_annot$Subtype[heatmap_annot$topic %in% classicalNRP_topics] <- "Classical-NRP"
heatmap_annot$Subtype[heatmap_annot$topic %in% basal_topics] <- "Basal Committed"
heatmap_annot$Subtype[heatmap_annot$topic %in% intermed_topics] <- "Subtype-Intermediate Progenitor"
rownames(heatmap_annot) <- heatmap_annot$topic
heatmap_annot <- heatmap_annot[,-1,drop=F]



cairo_pdf(width = 6.5, height = 5,
          file = "~/Topics_PlasmaB_vs_PDAC.pdf",
          family = 'Arial Unicode MS')

ComplexHeatmap::pheatmap(t(corr_pheatmap[topics_to_keep_immune,topics_to_keep_epithelial]),
                         clustering_distance_rows = "euclidean",
                         clustering_distance_cols = "euclidean",
                         clustering_method = "ward",
                         
                         annotation_row = heatmap_annot,
                         annotation_colors = heatmap_colors,
                         
                         # scale = "none",
                         
                         cellwidth = 8,
                         cellheight = 8,
                         
                         fontsize = 6,
                         angle_col = "90",
                         
                         color = colorRampPalette(ArchR::ArchRPalettes$solarExtra)(1200),
                         breaks = seq(from = -max(abs(corr_pheatmap)),
                                      to = max(abs(corr_pheatmap)),
                                      length.out = 1200),
                         annotation_legend = T,
                         
                         labels_row = paste("T",gsub("^Topic_|_Mean_.*", "", rownames(t(corr_pheatmap[topics_to_keep_immune,topics_to_keep_epithelial]))),sep=""),
                         labels_col = paste("T",gsub("^Topic_|_Mean_.*", "", colnames(t(corr_pheatmap[topics_to_keep_immune,topics_to_keep_epithelial]))),sep=""),
                         
                         
                         treeheight_row = 15,
                         treeheight_col = 15,
                         
                         display_numbers = t(corr_signif),
                         fontsize_number = 8,
                         number_color = "gray95",
                         
                         border_color = NA,
                         row_dend_side ="right",
                         row_names_side = "left",
                         heatmap_legend_param = list(
                             title = "Value",
                             at = c(-max(abs(corr_pheatmap)), max(abs(corr_pheatmap)))
                         )
)
dev.off()


