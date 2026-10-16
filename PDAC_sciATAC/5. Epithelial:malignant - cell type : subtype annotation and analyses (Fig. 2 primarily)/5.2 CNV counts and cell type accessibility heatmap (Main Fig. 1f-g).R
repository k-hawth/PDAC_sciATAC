library(ArchR)
library(dplyr)
library(parallel)
addArchRGenome("hg38")
addArchRThreads(10)
set.seed(1)

PlotPrefix = "cisTopic50_AllSamples_TSSMixed_nFrags1k_FRIS0.4_FetalAdultAtlas_PeakAnno"



############################################################
# Load ArchRProject

proj <- loadArchRProject("~/ArchR/Epithelial_AllSamples_TSSMixed_nFrags1k_FRIS0.4_FetalAdultAtlas_cisTopic50_PeakAnno/")



############################################################
# Add background peak set by corresponding to all cell type/lineage cCREs collated by the chromatin atlas by Zhang K and Hocker JD, et al. (2021) (https://doi.org/10.1016/j.cell.2021.10.024)

peaks <- read.table("~/peak_sets/FetalAdultAtlas_Zhang_2021/cCREs_noAlt.bed")
peaks <- GRanges(seqnames = peaks$V1, ranges = IRanges(start = peaks$V2, end = peaks$V3))

proj <- addPeakSet(ArchRProj = proj, peakSet = peaks,force=T)
proj <- addPeakMatrix(proj,force=T)
proj <- addBgdPeaks(proj, method = "chromVAR", force = T)



#########################################################
# Annotate ArchRProject peaks using cell type/lineage peak sets in the chromatin atlas by Zhang K and Hocker JD, et al. (2021) (https://doi.org/10.1016/j.cell.2021.10.024)
# This section calculates per-cell chromVAR scores for module cCRE peak sets (150 total) by Zhang K and Hocker JD, et al. (2021)

peakAnnotName <- "module_cCRE"
peak_set_directory = "~/peak_sets/FetalAdultAtlas_Zhang_2021/module_cCREs/"

module_list <- data.frame(module_file = list.files(peak_set_directory))
module_list$module_name <- gsub(".txt","",module_list$module_file)
celltype_cCRE <- c(paste(peak_set_directory,module_list$module_file, sep = ""))
names(celltype_cCRE) <- module_list$module_name
rm(module_list)



# Annotate peaks and run chromVAR
proj <- addPeakAnnotations(ArchRProj = proj, regions = celltype_cCRE, name = peakAnnotName, force=T)

proj <- addDeviationsMatrix(
    ArchRProj = proj, 
    peakAnnotation = peakAnnotName,
    force = TRUE
)



#########################################################
# Export per-cell scores for module cCRE peak sets
 
matrix <- getMatrixFromProject(proj, useMatrix = paste(peakAnnotName,"Matrix",sep=""), logFile = NULL, binarize = F)
matrix <- assay(matrix, "z")
matrix <- t(matrix)

write.table(as.data.frame(as.matrix(matrix)), file = paste(PlotPrefix,"_chromVAR_Module_cCRE.txt",sep=""), row.names = T, sep = '\t', col.names = T, quote = F)



#########################################################
# Reload matrix and calculate mean scores per cluster

matrix <- read.table("~/cisTopic50_AllSamples_TSSMixed_nFrags1k_FRIS0.4_FetalAdultAtlas_PeakAnno_chromVAR_Module_cCRE.txt", header = T, sep = '\t', comment.char = "", quote = "", check.names = F)

# Add PhenoGraph clustering annotation
matrix$Cluster <- proj$cistopic50_pg_k15[match(rownames(matrix), proj$cellNames)]


# Calculate mean module cCRE scores per cluster
library(tidyr)
plot <- matrix %>%
    group_by(Cluster) %>%
    summarise(across(where(is.numeric), mean)) %>%
    pivot_longer(
        cols = -Cluster,
        names_to = "variable",
        values_to = "Score"
    )

heatmap <- matrix %>%
    group_by(Cluster) %>%
    summarise(across(where(is.numeric), mean))
heatmap <- as.data.frame(heatmap)
rownames(heatmap) <- heatmap$Cluster
heatmap <- heatmap %>% 
    select(`72_Pancreatic_Epithelial`,`93_Gastric_Epithelial`,`65_Pulm_Epithelial`,`74_GI_Epithelial`,`85_Neuroendocrine`)

# Order heatmap by clusters
heatmap <- heatmap[c("pg_11", # acinar/ductal
                     "pg_18","pg_20", # metaplastic NRP
                     "pg_15", # islet
                     "pg_21", # pulmonary epithelial
                     "pg_7", # metaplastic gastric
                     "pg_2","pg_8", # intestinal epithelial
                     "pg_23", # other PAV
                     "pg_4","pg_16", # malignant PAV
                     "pg_1","pg_3","pg_5","pg_6","pg_9","pg_10","pg_12","pg_13","pg_14","pg_17","pg_19","pg_22"),] # malignant


#########################################################
# Load CNV counts per cell and calculate mean total CNVs per cluster

cnv <- read.table("~/epiAneufinder_SummaryOfEvents_PerCell.txt", header = T, sep = '\t', comment.char = "", quote = "")

cnv$Cluster <- proj$cistopic50_pg_k15[match(cnv$Cell,proj$cellNames)]
cnv$Total <- cnv$Gain + cnv$Loss

cnv <- cnv[is.na(cnv$Cluster)==F,]

cnv <- cnv %>%
    group_by(Cluster) %>%
    mutate(mean_cnv = mean(Total)) %>%
    select(Cluster,mean_cnv) %>%
    unique()

cnv <- as.data.frame(cnv)
rownames(cnv) <- cnv$Cluster
cnv <- cnv[,"mean_cnv",drop=F]
cnv <- cnv[rownames(heatmap),,drop=F]
cnv$mean_cnv <- log10(cnv$mean_cnv)


#########################################################
# Plot heatmap and CNV annotation

library(pheatmap)
pheatmap(heatmap,
         cluster_rows = F,
         cluster_cols = F,
         annotation_row = cnv,
         annotation_colors = list(mean_cnv = colorRampPalette(c(ArchRPalettes$blueYellow))(4000)),
         cellwidth = 9,
         cellheight = 8,
         fontsize_row = 6,
         fontsize_col = 6,
         labels_row = gsub("pg_","",rownames(heatmap)),
         border_color = NA,
         color = colorRampPalette(c("white",ArchRPalettes$comet,"black"))(100)
)