library(ArchR)
library(dplyr)
library(parallel)
addArchRGenome("hg38")
addArchRThreads(10)
set.seed(1)



############################################################
# Load ArchRProject

proj <- loadArchRProject("~/ArchR/Primary_TSSMixed_nFrags1k_FRIS0.4/Fibroblast_MACS2ArchR_cisTopic27/")



############################################################
# Add final cell type annotation

annot <- data.frame(cell = proj$cellNames,
                    sample = proj$Sample_ID,
                    pg = proj$cistopic27_pg_k50_ClusterFormat)

annot$celltype <- NA
annot$celltype[annot$pg %in% c("C1")] <- "LRRC15+ myCAF"
annot$celltype[annot$pg %in% c("C2")] <- "NGFR+ myCAF"
annot$celltype[annot$pg %in% c("C3","C4")] <- "iCAF"
annot$celltype[annot$pg %in% c("C5","C6")] <- "PSC"
annot$celltype[annot$pg %in% c("C7")] <- "Endothelial"
annot$celltype[annot$pg %in% c("C8")] <- "Unknown"

# Add to ArchRProject
proj <- addCellColData(proj, cells = annot$cell, data = annot$celltype, force = T, name = "celltype")



############################################################
# Plot UMAP of cell types

colorPalette <- c(
    "LRRC15+ myCAF" = "#E23C3C",
    "NGFR+ myCAF" = "#b176ef",
    "iCAF" = "#FFB132",
    "PSC" = "#89DDFF",
    "Endothelial" = "#3ECE39",
    "Unknown" = "gray88"
)


plotEmbedding(proj,
              embedding = "cisTopic_27",
              name = "celltype",
              pal = colorPalette)



############################################################
# Save ArchRProject

proj <- saveArchRProject(proj)