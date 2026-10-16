library(ArchR)
library(dplyr)
library(parallel)
addArchRGenome("hg38")
addArchRThreads(10)
set.seed(1)

PlotPrefix = "cisTopic27_FibroblastTSSMixnFrags1kFRIS0.4_MACS2ArchR"



############################################################
# Load ArchRProject

proj <- loadArchRProject("~/ArchR/Primary_TSSMixed_nFrags1k_FRIS0.4/Fibroblast_MACS2ArchR_cisTopic27/")





########################################################################################################################
# Marker peak analysis between clusters
########################################################################################################################


############################################################
# Marker peak analysis between clusters

set.seed(1)
markersPeaks <- getMarkerFeatures(
    ArchRProj = proj, 
    useMatrix = "PeakMatrix", 
    groupBy = "cistopic27_pg_k50_ClusterFormat",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "binomial",
    binarize = T
)

markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.05 & Log2FC >= 0.5")
peakmarkers <- as.data.frame(markerList)



############################################################
# Save output

PlotDestination = paste(getOutputDirectory(proj), "Plots", sep = "/")

saveRDS(markersPeaks, file = paste(PlotDestination, "/", PlotPrefix, "_MarkerPeaks_pgk50Clusters.rds",sep=""))
write.table(peakmarkers, file = paste(PlotDestination, "/", PlotPrefix, "_MarkerPeaks_pgk50Clusters.txt",sep=""), row.names=F, col.names = T, sep = '\t', quote = F)



############################################################
# Plot heatmap of marker peaks between clusters

markersPeaks <- readRDS("~/ArchR/Primary_TSSMixed_nFrags1k_FRIS0.4/Fibroblast_MACS2ArchR_cisTopic27/Plots/cisTopic27_FibroblastTSSMixnFrags1kFRIS0.4_MACS2ArchR_MarkerPeaks_pgk50Clusters.rds")

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
# Marker peak analysis between broad cell types
########################################################################################################################


############################################################
# Annotate cell types

annot <- data.frame(cell = proj$cellNames,
                    sample = proj$Sample_ID,
                    cluster = proj$cistopic27_pg_k50)
annot$cluster[annot$cluster %in% c("pg_5","pg_7")] = "Stellate-like"
annot$cluster[annot$cluster %in% c("pg_4")] = "Unknown"
annot$cluster[annot$cluster %in% c("pg_6")] = "Endothelial"
annot$cluster[annot$cluster %in% c("pg_1")] = "myCAF 1"
annot$cluster[annot$cluster %in% c("pg_3")] = "myCAF 2"
annot$cluster[annot$cluster %in% c("pg_2","pg_8")] = "iCAF"

proj <- addCellColData(proj, cells = annot$cell, data = annot$cluster, name = "celltype", force = T)



############################################################
# Exclude the unknown cluster from marker peak analysis

proj <- proj[proj$celltype!="Unknown"]



############################################################
# Marker peak analysis between clusters

markersPeaks <- getMarkerFeatures(
    ArchRProj = proj, 
    useMatrix = "PeakMatrix", 
    groupBy = "celltype",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "binomial",
    binarize = T
)

peakmarkers <- as.data.frame(getMarkers(markersPeaks, cutOff = "FDR <= 0.01 & Log2FC >= 1"))
rownames(peakmarkers) <- 1:length(rownames(peakmarkers))
peakmarkers$peak <- paste(peakmarkers$seqnames,peakmarkers$start,peakmarkers$end,sep="_")
peakmarkers <- chipseek(peakmarkers, "seqnames", "start", "end")



############################################################
# Save output

PlotDestination = paste(getOutputDirectory(proj), "Plots", sep = "/")
saveRDS(markersPeaks, file = paste(PlotDestination, "/", PlotPrefix, "_MarkerPeaks_CAFSubtype_NoUnknown_RObject.rds",sep=""))



############################################################
# Plot heatmap of marker peaks between cell types

markersPeaks <- readRDS("~/ArchR/Primary_TSSMixed_nFrags1k_FRIS0.4/Fibroblast_MACS2ArchR_cisTopic27/Plots/cisTopic27_FibroblastTSSMixnFrags1kFRIS0.4_MACS2ArchR_MarkerPeaks_CAFSubtype_NoUnknown_RObject.rds")

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