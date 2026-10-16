library(ArchR)
set.seed(1)
addArchRGenome("hg38")
library(parallel)
library(dplyr)

library(pheatmap)



############################################################
# Load ArchRProject

proj <- loadArchRProject("~/ArchR/Epithelial_AllSamples_TSSMixed_nFrags1k_FRIS0.4_PeaksMACS3_cisTopic50/")



########################################################################################################################
# Euclidean distance analysis
########################################################################################################################

# For each cell, calculate the mean distance of the cell to acinar/ductal cells.

cisTopic_euclideanDistance <- function(
    cistopic_matrix,  # cols = cells, rows = topics
    query_cells,      # cells to get mean distances for
    target_cells      # distances are calculated to these cells
){
    
    distances <- cistopic_matrix[,colnames(cistopic_matrix) %in% c(query_cells,target_cells)]
    
    distances <- t(distances)
    
    distances <- as.matrix(dist(distances,method="euclidean"))
    
    
    distances <- as.data.frame(distances)
    distances <- distances[,colnames(distances) %in% query_cells]
    
    distances <- distances[rownames(distances) %in% target_cells,]
    
    distances <- tidyr::gather(distances)
    
    distances <- distances %>%
        group_by(key) %>%
        summarise(mean = mean(value))
    
    return(distances)
}



# Load cell-by-topic score matrix
cistopic <- read.table("~/cisTopic/AllSamples_TSSMixed_nFrags1k_FRIS0.4/peakset_MACS3_AllSamplesTSSMixednFrags1kFRIS0.4/Epithelial/cisTopic100_Alpha50/Perplexity_topics/cisTopic100Alpha50_EpithelialAllSamplesTSSMixnFrags1kFRIS0.4_AllSamplesMACS3_Topics50_Matrix.matrix")

# Transpose: rows = cells, cols = topics
cistopic <- t(cistopic)

# Add ArchRProjct prefix
rownames(cistopic) <- paste("sears_tumors_hg38#",rownames(cistopic),sep="")

# Reorder cells (rows) to match in ArchRProject
cistopic <- cistopic[proj$cellNames,]

# Transpose: rows = topics, cols = cells (for euclidean distance function)
cistopic <- t(cistopic)




target_cell_type <- c("Acinar/Ductal","Metaplastic NRP","Metaplastic Gastric","Intestinal Epithelial")

# Remove cells from metastases
proj <- proj[grep("-T",proj$Sample_ID)]

plot <- data.frame()
for(unique_celltype in target_cell_type){
    test <- cisTopic_euclideanDistance(
        cistopic_matrix = cistopic,
        query_cells = proj$cellNames[grep("Basal|Classical|Intermediate|Malignant PAV",proj$celltype)],
        target_cells = proj$cellNames[proj$celltype==unique_celltype])
    
    # Add metadata
    test$Subtype <- proj$celltype[match(test$key,proj$cellNames)]
    
    test <- test %>%
        group_by(Subtype) %>%
        mutate(mean_mean = mean(mean)) %>%
        select(Subtype, mean_mean) %>%
        unique()
    colnames(test) <- gsub("^mean_mean$",unique_celltype,colnames(test))
    
    if(length(colnames(plot))==0){
        plot[1:length(test$Subtype),"Subtype"] <- test$Subtype
    }
    
    plot <- dplyr::left_join(plot, test, by = "Subtype")
}
rownames(plot) <- plot$Subtype
plot <- plot[,colnames(plot)!="Subtype",drop=F]




# Prepare heatmap
pheatmap(plot,
         color = colorRampPalette(rev(ArchRPalettes$solarExtra))(1000),
         fontsize = 6,
         cellwidth = 40,
         cellheight = 12,
         angle_col = 0,
         treeheight_row = 8,
         treeheight_col = 8,
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "ward")
