# This analysis identifies subtype-specific differentially accessible peaks, which are exported to be used later on for gene signature development.


library(ArchR)
library(dplyr)
library(parallel)
addArchRGenome("hg38")
addArchRThreads(10)
set.seed(1)



############################################################
# Save ArchRProject with only primary PDAC tumors, excluding PAV tumors

proj <- loadArchRProject("~/ArchR/Epithelial_AllSamples_TSSMixed_nFrags1k_FRIS0.4_PeaksMACS3_cisTopic50")

# Remove ampulla of vater samples
proj <- proj[proj$Sample_ID %in% c("ST-00017310-T","ST-00019171-T")==F]

# Remove metastatic samples
proj <- proj[grep("-T",proj$Sample_ID)]

proj <- saveArchRProject(proj, dropCells = T, outputDirectory = "~/ArchR/Epithelial_TSSMixed_nFrags1k_FRIS0.4_PeaksMACS3_cisTopic50_PrimaryPDACOnly")



########################################################################################################################
# DIFFERENTIAL ACCESSIBILITY: Primary tumor epigenetic subtypes
# In this section, differential accessibility tests will be performed between each epigenetic subtype 
########################################################################################################################

############################################################
# Differential accessibility: basal vs. classical (primaries only)

group1 = "Basal Committed"
group2 = "Classical-NRP"

set.seed(1)
markersPeaks <- getMarkerFeatures(
    ArchRProj = proj,
    useMatrix = "PeakMatrix", 
    groupBy = "celltype",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    useGroups = group1,
    bgdGroups = group2,
    testMethod = "binomial",
    binarize = T
)

markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.2 & abs(Log2FC) >= 0.5")
peakmarkers <- as.data.frame(markerList)
rownames(peakmarkers) <- 1:length(rownames(peakmarkers))
peakmarkers$group_name[peakmarkers$Log2FC>0] <- group1
peakmarkers$group_name[peakmarkers$Log2FC<0] <- group2
peakmarkers$peak <- paste(peakmarkers$seqnames,peakmarkers$start,peakmarkers$end, sep="_")

basal_vs_classicalNRP <- peakmarkers

write.table(peakmarkers, file = paste(getOutputDirectory(proj),"/Plots/","DifferentialAccessibility_PrimaryPDAC_Basal_vs_ClassicalNRP_FDR0.01_Log2FC1.0.txt",sep=""), row.names = F, col.names = T, sep = '\t', quote = F)



############################################################
# Differential accessibility: intermediate vs. basal (primaries only)

group1 = "Subtype-Intermediate Progenitor"
group2 = "Basal Committed"

set.seed(1)
markersPeaks <- getMarkerFeatures(
    ArchRProj = proj,
    useMatrix = "PeakMatrix", 
    groupBy = "celltype",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    useGroups = group1,
    bgdGroups = group2,
    testMethod = "binomial",
    binarize = T
)

markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.2 & abs(Log2FC) >= 0.5")
peakmarkers <- as.data.frame(markerList)
rownames(peakmarkers) <- 1:length(rownames(peakmarkers))
peakmarkers$group_name[peakmarkers$Log2FC>0] <- group1
peakmarkers$group_name[peakmarkers$Log2FC<0] <- group2
peakmarkers$peak <- paste(peakmarkers$seqnames,peakmarkers$start,peakmarkers$end, sep="_")

intermediate_vs_basal <- peakmarkers

write.table(peakmarkers, file = paste(getOutputDirectory(proj),"/Plots/","DifferentialAccessibility_PrimaryPDAC_Intermediate_vs_Basal_FDR0.01_Log2FC1.0.txt",sep=""), row.names = F, col.names = T, sep = '\t', quote = F)



############################################################
# Differential accessibility: classical vs. intermediate (primaries only)

group1 = "Classical-NRP"
group2 = "Subtype-Intermediate Progenitor"

set.seed(1)
markersPeaks <- getMarkerFeatures(
    ArchRProj = proj,
    useMatrix = "PeakMatrix", 
    groupBy = "celltype",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    useGroups = group1,
    bgdGroups = group2,
    testMethod = "binomial",
    binarize = T
)

markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.2 & abs(Log2FC) >= 0.5")
peakmarkers <- as.data.frame(markerList)
rownames(peakmarkers) <- 1:length(rownames(peakmarkers))
peakmarkers$group_name[peakmarkers$Log2FC>0] <- group1
peakmarkers$group_name[peakmarkers$Log2FC<0] <- group2
peakmarkers$peak <- paste(peakmarkers$seqnames,peakmarkers$start,peakmarkers$end, sep="_")

classicalNRP_vs_intermediate <- peakmarkers

write.table(peakmarkers, file = paste(getOutputDirectory(proj),"/Plots/","DifferentialAccessibility_PrimaryPDAC_Classical_vs_Intermediate_FDR0.01_Log2FC1.0.txt",sep=""), row.names = F, col.names = T, sep = '\t', quote = F)



########################################################################################################################
# PEAK INTERSECTIONS: Primary PDAC malignant subtypes
# This section will identify peaks that are more accessible in one subtype compared to each of the other two subtypes.
########################################################################################################################

# Peak intersections function
intersectPeaks <- function(peakObjList){
    output <- peakObjList[[1]]
    for(peakset in peakObjList[2:length(peakObjList)]){
        output <- output[output$peak %in% peakset$peak,]
    }
    return(output)
}



# Peak quality filters
filter_fdr = 0.01
filter_log2fc = 1



############################################################
# Intersect peaks: SIP

peakSetList <- list()

filter_group = "Subtype-Intermediate Progenitor"
peakmarkers <- intermediate_vs_basal
peakmarkers <- peakmarkers[peakmarkers$group_name==filter_group & peakmarkers$FDR<=filter_fdr & abs(peakmarkers$Log2FC)>=filter_log2fc,]
peakSetList[[length(peakSetList)+1]] <- peakmarkers

filter_group = "Subtype-Intermediate Progenitor"
peakmarkers <- classicalNRP_vs_intermediate
peakmarkers <- peakmarkers[peakmarkers$group_name==filter_group & peakmarkers$FDR<=filter_fdr & abs(peakmarkers$Log2FC)>=filter_log2fc,]
peakSetList[[length(peakSetList)+1]] <- peakmarkers

specific_intermediate <- intersectPeaks(peakSetList)
dim(intersectPeaks(peakSetList))

write.table(intersectPeaks(peakSetList), file = paste(getOutputDirectory(proj),"/Plots/","DifferentialAccessibility_PrimaryPDAC_Specific_Intermediate.txt",sep=""), row.names = F, col.names = T, sep = '\t', quote = F)



############################################################
# Intersect peaks: basal committed

peakSetList <- list()

filter_group = "Basal Committed"
peakmarkers <- intermediate_vs_basal
peakmarkers <- peakmarkers[peakmarkers$group_name==filter_group & peakmarkers$FDR<=filter_fdr & abs(peakmarkers$Log2FC)>=filter_log2fc,]
peakSetList[[length(peakSetList)+1]] <- peakmarkers

filter_group = "Basal Committed"
peakmarkers <- basal_vs_classicalNRP
peakmarkers <- peakmarkers[peakmarkers$group_name==filter_group & peakmarkers$FDR<=filter_fdr & abs(peakmarkers$Log2FC)>=filter_log2fc,]
peakSetList[[length(peakSetList)+1]] <- peakmarkers

specific_basal <- intersectPeaks(peakSetList)
dim(intersectPeaks(peakSetList))

write.table(intersectPeaks(peakSetList), file = paste(getOutputDirectory(proj),"/Plots/","DifferentialAccessibility_PrimaryPDAC_Specific_Basal.txt",sep=""), row.names = F, col.names = T, sep = '\t', quote = F)



############################################################
# Intersect peaks: classical-NRP

peakSetList <- list()

filter_group = "Classical-NRP"
peakmarkers <- classicalNRP_vs_intermediate
peakmarkers <- peakmarkers[peakmarkers$group_name==filter_group & peakmarkers$FDR<=filter_fdr & abs(peakmarkers$Log2FC)>=filter_log2fc,]
peakSetList[[length(peakSetList)+1]] <- peakmarkers

filter_group = "Classical-NRP"
peakmarkers <- basal_vs_classicalNRP
peakmarkers <- peakmarkers[peakmarkers$group_name==filter_group & peakmarkers$FDR<=filter_fdr & abs(peakmarkers$Log2FC)>=filter_log2fc,]
peakSetList[[length(peakSetList)+1]] <- peakmarkers

specific_classical <- intersectPeaks(peakSetList)
dim(intersectPeaks(peakSetList))

write.table(intersectPeaks(peakSetList), file = paste(getOutputDirectory(proj),"/Plots/","DifferentialAccessibility_PrimaryPDAC_Specific_ClassicalNRP.txt",sep=""), row.names = F, col.names = T, sep = '\t', quote = F)



############################################################
# Save ArchRProject

proj <- saveArchRProject(proj)