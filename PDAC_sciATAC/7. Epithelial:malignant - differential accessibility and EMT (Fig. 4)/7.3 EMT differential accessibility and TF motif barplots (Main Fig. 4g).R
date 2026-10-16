library(ArchR)
library(dplyr)
library(parallel)
addArchRGenome("hg38")
addArchRThreads(10)
set.seed(1)

PlotPrefix = "Epithelial_AllSamples_TSSMixed_nFrags1k_FRIS0.4_PeaksMACS3_cisTopic50"



############################################################
# Load ArchRProject

proj <- loadArchRProject("~/ArchR/Epithelial_AllSamples_TSSMixed_nFrags1k_FRIS0.4_PeaksMACS3_cisTopic50")



############################################################
# EMT groups

annot <- data.frame(cell = proj$cellNames,
                    group = proj$cistopic50_pg_k15)

# Epithelial
annot$group[annot$group %in% c("pg_1","pg_10","pg_12","pg_17")] = "Epithelial"

# Partial EMT
annot$group[annot$group %in% c("pg_5","pg_6","pg_22","pg_3","pg_14")] = "Partial EMT"

# Mesenchymal
annot$group[annot$group %in% c("pg_13","pg_19","pg_9")] = "Mesenchymal"

# Add to ArchRProject
proj <- addCellColData(proj, data = annot$group, cells = annot$cell, force = T, name = "EMT_Broad")



########################################################################################################################
# DIFFERENTIAL ACCESSIBILITY: Broad groups
########################################################################################################################

##################################################################
# Epithelial vs. Partial EMT

group1 = "Epithelial"
group2 = "Partial EMT"

set.seed(1)
markersPeaks <- getMarkerFeatures(
    ArchRProj = proj,
    useMatrix = "PeakMatrix", 
    groupBy = "EMT_Broad",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    useGroups = group1,
    bgdGroups = group2,
    testMethod = "binomial",
    binarize = T
)

markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.01 & abs(Log2FC) >= 1")
peakmarkers <- as.data.frame(markerList)
rownames(peakmarkers) <- 1:length(rownames(peakmarkers))
peakmarkers$group_name[peakmarkers$Log2FC>0] <- group1
peakmarkers$group_name[peakmarkers$Log2FC<0] <- group2
peakmarkers$peak <- paste(peakmarkers$seqnames,peakmarkers$start,peakmarkers$end, sep="_")

epi_vs_partial <- peakmarkers



##################################################################
# Epithelial vs. Mesenchymal

group1 = "Epithelial"
group2 = "Mesenchymal"

set.seed(1)
markersPeaks <- getMarkerFeatures(
    ArchRProj = proj,
    useMatrix = "PeakMatrix", 
    groupBy = "EMT_Broad",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    useGroups = group1,
    bgdGroups = group2,
    testMethod = "binomial",
    binarize = T
)

markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.01 & abs(Log2FC) >= 1")
peakmarkers <- as.data.frame(markerList)
rownames(peakmarkers) <- 1:length(rownames(peakmarkers))
peakmarkers$group_name[peakmarkers$Log2FC>0] <- group1
peakmarkers$group_name[peakmarkers$Log2FC<0] <- group2
peakmarkers$peak <- paste(peakmarkers$seqnames,peakmarkers$start,peakmarkers$end, sep="_")

epi_vs_mesenchymal <- peakmarkers



##################################################################
# Partial EMT vs. Mesenchymal

group1 = "Partial EMT"
group2 = "Mesenchymal"

set.seed(1)
markersPeaks <- getMarkerFeatures(
    ArchRProj = proj,
    useMatrix = "PeakMatrix", 
    groupBy = "EMT_Broad",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    useGroups = group1,
    bgdGroups = group2,
    testMethod = "binomial",
    binarize = T
)

markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.01 & abs(Log2FC) >= 1")
peakmarkers <- as.data.frame(markerList)
rownames(peakmarkers) <- 1:length(rownames(peakmarkers))
peakmarkers$group_name[peakmarkers$Log2FC>0] <- group1
peakmarkers$group_name[peakmarkers$Log2FC<0] <- group2
peakmarkers$peak <- paste(peakmarkers$seqnames,peakmarkers$start,peakmarkers$end, sep="_")

partial_vs_mesenchymal <- peakmarkers



########################################################################################################################
# TF Motif Enrichment
# This section uses ArchR functions to count the number of peaks (in a given peak set) that have a TF motif
# First, we have ArchR identify TF motifs present in the background peak set (all MACS v3 called peaks)
# Second, we use ArchR's peakAnnoEnrichment to count the number of peaks in each peak set with each TF motif
# Third, we perform a hypergeometric test for each TF motif examined and correct the p-values
########################################################################################################################


############################################################
# Function: use ArchR to count the number of query peaks and background peaks containing a specific motif

count_motif_hits <- function(peaks_to_test = peaks_to_test,
                             peakset_name = "Peak_Set",
                             name_annot = "Motif",
                             ArchRProject = proj,
                             SE_object = markersPeaks){
    
    # Flag all peaks in peaks_to_test, unflag all other background peaks
    SE_peaks <- as.data.frame(rowData(SE_object))
    keep <- assay(SE_object,"FDR")
    colnames(keep) <- "keep"
    keep$keep <- 0
    keep$peak <- paste(SE_peaks$seqnames,SE_peaks$start,SE_peaks$end,sep="_")
    keep$keep[keep$peak %in% peaks_to_test$peak] = 1
    keep <- keep[,"keep",drop=F]
    assay(SE_object,"keep", withDimnames = F) <- keep
    
    # Run ArchR peakAnnoEnrichment
    enrichEncode <- peakAnnoEnrichment(
        seMarker = SE_object,
        ArchRProj = ArchRProject,
        peakAnnotation = name_annot,
        cutOff = "keep == 1", # keep only peaks in peaks_to_test
        background = "all"
    )
    
    results <- data.frame(TF = rownames(enrichEncode))
    for(assay in names(enrichEncode@assays)){
        temp1 <- as.data.frame(assay(enrichEncode,assay))
        colnames(temp1) <- "assay"
        results[assay] <- temp1[match(results$TF,rownames(temp1)),"assay"]
    }
    rm(assay,temp1)
    results <- results[order(results$Enrichment,decreasing = T),]
    results <- data.frame(Query_Peak_Set = peakset_name,
                          TF_Motif = results$TF,
                          Enrichment = results$Enrichment,
                          Query_Hits = results$CompareFrequency,
                          Query_PeakCount = results$nCompare,
                          Query_Proportion_Hits = results$CompareProportion,
                          Background_Hits = results$BackgroundFrequency,
                          Background_PeakCount = results$nBackground)
    return(results)
}



############################################################
# Initiate random marker peaks test
# the results of this test do not matter; this is only to generate the SummarizedExperiment object that the function above will take to perform the custom enrichment

group1 = "pg_1"
group2 = "pg_2"

set.seed(1)
markersPeaks <- getMarkerFeatures(
    ArchRProj = proj,
    useMatrix = "PeakMatrix", 
    groupBy = "cistopic50_pg_k15",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    useGroups = group1,
    bgdGroups = group2,
    testMethod = "binomial",
    binarize = T
)



############################################################
# Count overlap of TF motifs in differentially accessible peaks

results_all <- data.frame()



# TF motif enrichment in peaks that are Epithelial > Partial EMT
peaks_to_test <- epi_vs_partial
peaks_to_test <- peaks_to_test[peaks_to_test$FDR<=0.01 & abs(peaks_to_test$Log2FC)>=1.0 & peaks_to_test$group_name=="Epithelial",]
results <- count_motif_hits(peaks_to_test = peaks_to_test,
                            peakset_name = "Epithelial > Partial EMT",
                            name_annot = "Motif",
                            SE_object = markersPeaks)
results$TF_Motif <- gsub("_.*","",results$TF_Motif)
results_all <- rbind(results_all, results)



# TF motif enrichment in peaks that are Partial EMT > Epithelial
peaks_to_test <- epi_vs_partial
peaks_to_test <- peaks_to_test[peaks_to_test$FDR<=0.01 & abs(peaks_to_test$Log2FC)>=1.0 & peaks_to_test$group_name=="Partial EMT",]
results <- count_motif_hits(peaks_to_test = peaks_to_test,
                            peakset_name = "Partial EMT > Epithelial",
                            name_annot = "Motif",
                            SE_object = markersPeaks)
results$TF_Motif <- gsub("_.*","",results$TF_Motif)
results_all <- rbind(results_all, results)



# TF motif enrichment in peaks that are Partial EMT > Mesenchymal
peaks_to_test <- partial_vs_mesenchymal
peaks_to_test <- peaks_to_test[peaks_to_test$FDR<=0.01 & abs(peaks_to_test$Log2FC)>=1.0 & peaks_to_test$group_name=="Partial EMT",]
results <- count_motif_hits(peaks_to_test = peaks_to_test,
                            peakset_name = "Partial EMT > Mesenchymal",
                            name_annot = "Motif",
                            SE_object = markersPeaks)
results$TF_Motif <- gsub("_.*","",results$TF_Motif)
results_all <- rbind(results_all, results)



# TF motif enrichment in peaks that are Mesenchymal > Partial EMT
peaks_to_test <- partial_vs_mesenchymal
peaks_to_test <- peaks_to_test[peaks_to_test$FDR<=0.01 & abs(peaks_to_test$Log2FC)>=1.0 & peaks_to_test$group_name=="Mesenchymal",]
results <- count_motif_hits(peaks_to_test = peaks_to_test,
                            peakset_name = "Mesenchymal > Partial EMT",
                            name_annot = "Motif",
                            SE_object = markersPeaks)
results$TF_Motif <- gsub("_.*","",results$TF_Motif)
results_all <- rbind(results_all, results)



############################################################
# Select EMT and pro-epithelial TF motifs

motifs_epithelial <- c("ONECUT2","HNF1B","GATA6","ELF3")
motifs_emt <- c("ZEB1","TWIST2","TWIST1","SNAI2","SNAI1")

results_all <- results_all[results_all$TF_Motif %in% c(motifs_epithelial,motifs_emt),]



############################################################
# Perform hypergeometric tests and correct p-values

balls_drawn_white = results_all$Query_Hits
total_white_balls = results_all$Background_Hits
total_black_balls = results_all$Background_PeakCount - results_all$Background_Hits
total_balls_drawn = results_all$Query_PeakCount

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

results_all$pval_log10 <- p_log10
results_all$pval_log10_Bonferroni <- p_log10_bonferroni






########################################################################################################################
# Plot barplot of TF motif enrichment results
########################################################################################################################

############################################################
# Prepare data for plotting

emt <- results_all

emt$Celltype <- gsub(" >.*","", emt$Query_Peak_Set)

# Set non-significant enrichment scores to 0
emt$Enrichment_log2 <- log2(emt$Enrichment)
emt$Enrichment_log2[emt$pval_log10_Bonferroni<(-log10(0.05))] <- 0

# Set left-hand bars to "negative" enrichment scores (for plotting purposes)
emt$Enrichment_log2[emt$Query_Peak_Set %in% c("Epithelial > Partial EMT", "Partial EMT > Mesenchymal")] <- -(emt$Enrichment_log2[emt$Query_Peak_Set %in% c("Epithelial > Partial EMT", "Partial EMT > Mesenchymal")])



############################################################
# ggplot2 themes and settings

# Color palette
pal <- c("Mesenchymal" = "#7e2bb7",
         "Partial EMT" = "#00a0e1",
         "Epithelial" = "#43963a")



# Themes
text_settings <- element_text(colour = "black", size = 6)
line_settings <- element_line(linewidth = 0.5)

default_theme <- theme(
    axis.title.x = text_settings,
    axis.title.y = text_settings,
    
    axis.ticks.x = line_settings,
    axis.ticks.y = line_settings,
    
    axis.text.x = text_settings,
    axis.text.y = text_settings,
    
    #axis.line.x = element_line(linewidth = 0.5),
    #axis.line.y = element_line(linewidth = 0.5),
    panel.border = element_rect(linewidth = 0.75),
    
    plot.margin = margin(t = 2,  # Top margin
                         r = 2,  # Right margin
                         b = 0,  # Bottom margin
                         l = 0)
    )



############################################################
# Plot data

plot <- emt[emt$Query_Peak_Set %in% c("Epithelial > Partial EMT", "Partial EMT > Epithelial"),]
ggplot(plot, aes(x = Enrichment_log2, y = TF, fill = Celltype)) +
    geom_bar(stat = "identity") +
    geom_vline(xintercept = 0, linewidth = 0.5, color = "black") +
    geom_hline(yintercept = 4.5, linewidth = 0.5, color = "gray30", linetype = 2) +
    theme_ArchR() +
    default_theme +
    theme(
        legend.position = "none",
        axis.title.y = element_blank()
    ) +
    scale_fill_manual(values = pal) +
    xlim(-max(abs(plot$Enrichment_log2)),
         max(abs(plot$Enrichment_log2))) +
    xlab("Enrichment (log2)")



plot <- emt[emt$Query_Peak_Set %in% c("Partial EMT > Mesenchymal", "Mesenchymal > Partial EMT"),]
ggplot(plot, aes(x = Enrichment_log2, y = TF, fill = Celltype)) +
    geom_bar(stat = "identity") +
    geom_vline(xintercept = 0, linewidth = 0.5, color = "black") +
    geom_hline(yintercept = 4.5, linewidth = 0.5, color = "gray30", linetype = 2) +
    theme_ArchR() +
    default_theme +
    theme(
        legend.position = "none",
        axis.title.y = element_blank()
    ) +
    scale_fill_manual(values = pal) +
    xlim(-max(abs(plot$Enrichment_log2)),
         max(abs(plot$Enrichment_log2))) +
    xlab("Enrichment (log2)")


