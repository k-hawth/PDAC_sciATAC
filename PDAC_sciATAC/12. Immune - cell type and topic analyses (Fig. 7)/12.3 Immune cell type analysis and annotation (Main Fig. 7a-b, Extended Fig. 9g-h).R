library(ArchR)
library(dplyr)
library(parallel)
addArchRGenome("hg38")
addArchRThreads(10)
set.seed(1)



############################################################
# Load ArchRProject

proj <- loadArchRProject("~/ArchR/Primary_TSSMixed_nFrags1k_FRIS0.4/Immune_MACS2ArchR_cisTopic43")



############################################################
# Marker peak analysis between T cell clusters and unknown clusters

markersPeaks <- getMarkerFeatures(
    ArchRProj = proj[proj$cistopic43_pg_k15 %in% c("pg_11","pg_8","pg_4","pg_3","pg_5","pg_7","pg_13")], 
    useMatrix = "PeakMatrix", 
    groupBy = "cistopic43_pg_k15",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "binomial",
    binarize = T
)



###########################################################
# Heatmap of marker peak analysis

heatmapPeaks <- plotMarkerHeatmap(
    seMarker = markersPeaks,
    log2Norm = TRUE,
    scaleRows = TRUE,
    nLabel = 0,
    cutOff = "FDR <= 0.05 & Log2FC >= 0.5",
    transpose = TRUE,
)

ComplexHeatmap::draw(heatmapPeaks,
                     heatmap_legend_side = "bot",
                     annotation_legend_side = "bot")



########################################################################################################################
# Cell type enhancer enrichment in immune marker peaks
# This section counts the overlap between 222 cell type peak sets with marker peak sets and background peaks
# The 222 cell type peak sets come from the chromatin atlas by Zhang K and Hocker JD, et al. (2021) (https://doi.org/10.1016/j.cell.2021.10.024)
# Hypergeometric tests then determine significance

# In this specific section, we do not test all 222 cell type peak sets—only immune cell types:
# Mast, CD4 T, CD8 T, plasma B, memory B, macrophage, and NK cell
########################################################################################################################

############################################################
# Query:
# Turn all marker peak sets (for each cluster) into a GRanges object that goes into topic_list

markerPeaks <- as.data.frame(getMarkers(markersPeaks, cutOff = "FDR <= 0.05 & Log2FC >= 0.5"))

# Topic Peak Sets
topic_list <- list()
for(unique_topic in unique(markerPeaks$group_name)){
    peakmarkers <- markerPeaks[markerPeaks$group_name==unique_topic,]
    peakmarkers <- GRanges(seqnames = peakmarkers$seqnames, IRanges(start = peakmarkers$start, end = peakmarkers$end))
    topic_list[[unique_topic]] <- peakmarkers
    rm(peakmarkers)
}


############################################################
# Enrichment:
# Load 222 cell type peak sets in the chromatin atlas by Zhang K and Hocker JD, et al. (2021) (https://doi.org/10.1016/j.cell.2021.10.024)

peakset_file_path <- "~/peak_sets/FetalAdultAtlas_Zhang_2021/celltype_cCREs"
peakset_files <- list.files(path = peakset_file_path, full.names = F)
peakset_list <- list()
for(unique_peakset in peakset_files){
    peakmarkers <- read.table(paste(peakset_file_path, unique_peakset, sep = "/"), header = F, sep = '\t')
    peakmarkers <- GRanges(seqnames = peakmarkers$V1, IRanges(start = peakmarkers$V2, end = peakmarkers$V3))
    peakset_list[[unique_peakset]] <- peakmarkers
    rm(peakmarkers)
}


# Filter for only immune cell types
names(peakset_list) <- gsub("\\.txt","",names(peakset_list))
celltypes <- c("Mast","T_Lymphocyte_1_CD8","T_lymphocyte_2_CD4","Plasma_B","Memory_B", "Macrophage_General","Natural_Killer_T")
peakset_list <- peakset_list[names(peakset_list) %in% celltypes]



############################################################
# Background:
# Load background set of called peaks

peaks <- read.table("~/ArchR/Primary_TSSMixed_nFrags1k_FRIS0.4/AllCells_MACS2ArchR/Plots/Peaks_MACS2ArchR_PrimaryTSSMixnFrags1kFRIS0.4.bed")
peaks <- GRanges(seqnames = peaks$V1, ranges = IRanges(start = peaks$V2, end = peaks$V3))



#########################################################
# Function to perform hypergeometric tests

peakSetEnrichment <- function(topic_list,
                              peakset_list,
                              bgd){
    output <- data.frame()
    
    for(unique_topic in names(topic_list)){
        for(unique_peakset in names(peakset_list)){
            peakset_topic <- topic_list[[unique_topic]]
            peakset_cCRE <- peakset_list[[unique_peakset]]
            peakset_all_topics <- bgd
            topic_length = length(peakset_topic)
            
            total_balls_drawn = topic_length
            balls_drawn_white = length(peakset_topic[overlapsAny(query = peakset_topic, subject = peakset_cCRE)])
            total_white_balls = length(peakset_all_topics[overlapsAny(query = peakset_all_topics, subject = peakset_cCRE)])
            total_black_balls = length(peakset_all_topics) - total_white_balls
            
            p_ln <- phyper(q = balls_drawn_white - 1,
                           m = total_white_balls,
                           n = total_black_balls,
                           k = total_balls_drawn,
                           lower.tail = F,
                           log.p = T)
            p_log10 <- -(p_ln/log(10))
            p_archr <- p_log10
            enrichment <- (balls_drawn_white/total_balls_drawn)/(total_white_balls/length(peakset_all_topics))
            p_archr
            enrichment
            
            output[length(rownames(output))+1,"Test_Set"] <- unique_topic
            output[length(rownames(output)),"Test_Peaks"] <- total_balls_drawn
            output[length(rownames(output)),"Test_Peak_Hits"] <- balls_drawn_white
            output[length(rownames(output)),"Bgd_Peaks"] <- total_white_balls + total_black_balls
            output[length(rownames(output)),"Bgd_Peak_Hits"] <- total_white_balls
            
            output[length(rownames(output)),"Peak_Set"] <- unique_peakset
            output[length(rownames(output)),"Enrichment"] <- enrichment
            output[length(rownames(output)),"mlog10p"] <- p_archr
            output[length(rownames(output)),"p_val"] <- exp(p_ln)
        }
    }
    return(output)
}



#########################################################
# Perform enrichments

test <- peakSetEnrichment(topic_list = topic_list,
                          peakset_list = peakset_list,
                          bgd = peaks)



#########################################################
# Calculate p-values

results_all <- test

# Correct p-values using all results
balls_drawn_white = results_all$Test_Peak_Hits
total_white_balls = results_all$Bgd_Peak_Hits
total_black_balls = results_all$Bgd_Peaks - results_all$Bgd_Peak_Hits
total_balls_drawn = results_all$Test_Peaks

p_ln <- phyper(q = balls_drawn_white - 1,
               m = total_white_balls,
               n = total_black_balls,
               k = total_balls_drawn,
               lower.tail = F,
               log.p = T)
p_log10 <- -(p_ln/log(10))
p_log10_bonferroni <- p_log10 - log10( length(p_log10) )

# Convert negative values to 0
p_log10_bonferroni[p_log10_bonferroni<0] <- 0

results_all$mlog10p <- p_log10
results_all$mlog10Padj <- p_log10_bonferroni



#########################################################
# Tabulate results

table <- data.frame(Query_Peak_Set = results_all$Test_Set,
                    Celltype_cCRE_Peak_Set = results_all$Peak_Set,
                    pval_log10 = results_all$mlog10p,
                    pval_log10_Bonferroni = results_all$mlog10Padj,
                    Enrichment = results_all$Enrichment,
                    Background_PeakCount = results_all$Bgd_Peaks,
                    Background_Hits = results_all$Bgd_Peak_Hits,
                    Background_Proportion_Hits = results_all$Bgd_Peak_Hits/results_all$Bgd_Peaks,
                    Query_PeakCount = results_all$Test_Peaks,
                    Query_Hits = results_all$Test_Peak_Hits,
                    Query_Proportion_Hits = results_all$Test_Peak_Hits/results_all$Test_Peaks)

table$Celltype_cCRE_Peak_Set <- gsub("_"," ",table$Celltype_cCRE_Peak_Set)



#########################################################
# Convert results table to matrix

results <- table

matrix_mlog10Padj <- data.frame()
for(unique_topic in unique(results$Query_Peak_Set)){
    for(unique_celltype in unique(results$Celltype_cCRE_Peak_Set)){
        matrix_mlog10Padj[unique_topic,unique_celltype] <- results$pval_log10_Bonferroni[results$Query_Peak_Set==unique_topic & results$Celltype_cCRE_Peak_Set==unique_celltype]
    }
}



#########################################################
# Scale 0 - 100
mat <- matrix_mlog10Padj
mat <- as.data.frame(mat)
colnames(mat) <- gsub("\\.txt","",colnames(mat))

scale_0_100 <- function(x) {
    (x - min(x, na.rm = TRUE)) /
        (max(x, na.rm = TRUE) - min(x, na.rm = TRUE)) * 100
}

mat <- as.data.frame(t(mat))
mat[] <- lapply(mat, scale_0_100)



#########################################################
# Plot data

colnames(mat) <- gsub("pg_","C",colnames(mat))
rownames(mat) <- gsub("_"," ",rownames(mat))


ComplexHeatmap::pheatmap(mat[,c("C5","C3","C7","C13","C8","C4","C11"),drop=F],
                         color = colorRampPalette(ArchR::ArchRPalettes$solarExtra)(1000),
                         cellwidth = 8, cellheight = 8,
                         treeheight_row = 15,
                         treeheight_col = 15,
                         border_color = NA,
                         cluster_cols = F,
                         cluster_rows = T,
                         fontsize = 6
)




########################################################################################################################
# Immune cell type annotation
########################################################################################################################


############################################################
# Annotate cell types

annot <- data.frame(cell = proj$cellNames,
                    sample = proj$Sample_ID,
                    cluster = proj$cistopic43_pg_k15)
annot$cluster[annot$cluster %in% c("pg_1","pg_2")] = "Myeloid"
annot$cluster[annot$cluster %in% c("pg_6","pg_12")] = "Plasma B"
annot$cluster[annot$cluster %in% c("pg_9","pg_10")] = "B cell"
annot$cluster[annot$cluster %in% c("pg_3","pg_4","pg_5","pg_8","pg_13")] = "T cell"
annot$cluster[annot$cluster %in% c("pg_7")] = "Mast"
annot$cluster[annot$cluster %in% c("pg_11")] = "Unknown"

proj <- addCellColData(proj, cells = annot$cell, data = annot$cluster, name = "celltype", force = T)



############################################################
# ggplot2 settings

colorPalette <- c("Myeloid" = "#ff6945",
                  "Mast" = "#ff9d38",
                  "T cell" = "#06629e",
                  "B cell" = "#b5838d",
                  "Plasma B" = "#6d6875",
                  "Unknown" = "gray80"
)


# Theme Settings
line_width <- 0.75
point_size <- 0.65
text_size <- 6
title_size <- 6

line_settings <- element_line(linewidth = line_width, color = "black")
text_settings <- element_text(colour = "black", size = text_size)
title_settings <- element_text(colour = "black", size = title_size, vjust = -3)


# Default Theme
default_theme <- theme(
    axis.title.x = element_blank(),
    axis.title.y = text_settings,
    
    axis.ticks.x = element_line(linewidth = line_width*0.5, color = "black"),
    axis.ticks.y = element_line(linewidth = line_width*0.5, color = "black"),
    
    axis.text.x = text_settings,
    axis.text.y = text_settings,
    
    axis.line.x = element_blank(),
    axis.line.y = element_blank(),
    
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "white"),
    panel.border = element_rect(linewidth = line_width, fill = NA),
    
    plot.background = element_blank(),
    plot.title = title_settings,
    
    legend.text = text_settings,
    legend.title = title_settings
)



############################################################
# Plot UMAP of immune cell types

plotEmbedding(proj,
              name = "celltype",
              embedding = "cisTopic_43",
              pal = colorPalette,
              labelMeans = F)



############################################################
# Plot proportions of immune cell types by sample

plot <- data.frame(Sample = proj$Sample_ID,
                   Cluster = proj$celltype)

plot <- plot %>%
    group_by(Sample,Cluster) %>%
    mutate(Count = n()) %>%
    unique()

# Factors: samples
sample_count <- data.frame(Sample = proj$Sample_ID)
sample_count <- sample_count %>%
    group_by(Sample) %>%
    mutate(Count = n()) %>%
    unique()

sample_count <- sample_count[order(sample_count$Count,decreasing = T),]

plot$Sample <- factor(plot$Sample, levels = rev(sample_count$Sample))

pdf(height = 1.6, width = 4.84, file = "immune_cluster_sampleProportions.pdf")
ggplot(plot, aes(x = Sample, fill = Cluster, y = Count)) +
    geom_bar(stat = "identity", position = "stack") +
    default_theme +
    ylab("Count") +
    ggtitle("Immune Composition of Primary Tumor sciATAC Profiles") +
    theme(
        legend.text = element_text(size = 6),
        legend.key.height = unit(1, "null"),
        plot.title = element_text(size = 6, color = "black", vjust = -1.5),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)
    ) +
    scale_fill_manual(values = colorPalette)

dev.off()



############################################################
# Save ArchRProject

proj <- saveArchRProject(proj)