library(ArchR)
library(dplyr)
library(parallel)
addArchRGenome("hg38")
addArchRThreads(20)
set.seed(1)



############################################################
# Prepare new ArchRProject

proj <- loadArchRProject("~/ArchR/Primary_TSSMixed_nFrags1k_FRIS0.4/AllCells")

proj <- saveArchRProject(proj, outputDirectory = "~/ArchR/Primary_TSSMixed_nFrags1k_FRIS0.4/AllCells_MACS2ArchR")



############################################################
# Call peaks

library(BSgenome.Hsapiens.UCSC.hg38)

# Add pseudobulk replicates for Clusters
proj <- addGroupCoverages(ArchRProj = proj, groupBy = "Clusters")

# Call peaks via MACS2 + ArchR iterative overlap method
proj <- addReproduciblePeakSet(
    ArchRProj = proj, 
    groupBy = "Clusters", 
    pathToMacs2 = "~/python_environments/MACS2/bin/macs2"
)

# Add called peak set as a PeakMatrix
proj <- addPeakMatrix(proj)

# Save peak set
peaks <- as.data.frame(getPeakSet(proj), row.names = 1:length(getPeakSet(proj)))
peaks <- peaks[c(1:3)]
write.table(peaks, file = paste(getOutputDirectory(proj), "/Plots/Peaks_MACS2ArchR_PrimaryTSSMixnFrags1kFRIS0.4.bed", sep = ""), row.names = F, col.names = F, sep = '\t', quote = F)

