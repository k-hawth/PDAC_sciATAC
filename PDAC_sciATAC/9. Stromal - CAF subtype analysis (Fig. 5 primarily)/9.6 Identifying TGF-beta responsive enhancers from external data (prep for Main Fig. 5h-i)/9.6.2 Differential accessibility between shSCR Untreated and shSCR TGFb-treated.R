########################################################################################################################
# Analysis of external sci-ATAC-seq dataset from Bowman CL and Daniel CJ, et al. (2025) (https://doi.org/10.1158/0008-5472.can-24-3437)

# In these data, there is a shPIN1 condition (PIN1 knockout) and a shSCR condition (control).
# Here, I use ArchR to perform differential accessibility between only the shSCR control conditions (shSCR untreated PSCs vs. shSCR TGFb-treated PSCs)

########################################################################################################################

library(ArchR)
set.seed(1)
addArchRGenome("hg38")
library(parallel)

proj <- loadArchRProject("~/PIN1_ArchR/Project_PDAC_sciATAC/")




############################################################
# Differential accessibility

group1 <- "shSCR UN"
group2 <- "shSCR TGFβ"

markersPeaks <- getMarkerFeatures(
    ArchRProj = proj, 
    useMatrix = "PeakMatrix", 
    groupBy = "sample_ID",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "binomial",
    binarize = T,
    useGroups = group1,
    bgdGroups = group2,
)

peakmarkers <- as.data.frame(getMarkers(markersPeaks, cutOff = "FDR <= 0.01 & abs(Log2FC) >= 1"))
rownames(peakmarkers) <- 1:length(rownames(peakmarkers))
peakmarkers$group_name[peakmarkers$Log2FC>0] <- group1
peakmarkers$group_name[peakmarkers$Log2FC<0] <- group2



############################################################
# Save differentially accessible peak sets

peakset <- peakmarkers[peakmarkers$group_name=="shSCR UN",c("seqnames","start","end")]
write.table(peakset, file = "~/peak_sets/PIN1_TGFb_2026/shSCR_UN_vs_shSCR_TGFb_Binomial_FDR0.01_log2FC1.0.txt", row.names = F, col.names = F, sep = '\t', quote = F)

peakset <- peakmarkers[peakmarkers$group_name=="shSCR TGFβ",c("seqnames","start","end")]
write.table(peakset, file = "~/peak_sets/PIN1_TGFb_2026/shSCR_TGFb_vs_shSCR_UN_Binomial_FDR0.01_log2FC1.0.txt", row.names = F, col.names = F, sep = '\t', quote = F)
