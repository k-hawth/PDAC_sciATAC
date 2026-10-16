library(gplots)
library(grid)
library(ComplexHeatmap)
library(circlize)
library(cluster)

# Set output plot name
OutputFileName = "AllSamplesTSSMixednFrags1kFRIS0.4_cisTopic50"



############################################################
# Read in epiAneufinder results

matrix <- read.table("~/epiAneufinder/AllSamples_TSSMixed_nFrags1k_FRIS0.4/epiAneufinder/results_table.tsv")



############################################################
# Get all cells cluster annotation

annot_allcells <- read.table("~/ArchR/AllSamples_TSSMixed_nFrags1k_FRIS0.4/Plots/AllSamples_TSSMixed_nFrags1k_FRIS0.4_ClusterAnnotation.annot")

# Remove any cells that are not in matrix
annot_allcells <- annot_allcells[annot_allcells $V1 %in% gsub("cell.","", colnames(matrix)),]

# Remove extra columns
annot_allcells <- annot_allcells[c(1,2)]

# Order cells by cluster number
cluster_prefix = "C"
annot_allcells <- annot_allcells[order(as.numeric(gsub(cluster_prefix,"", annot_allcells$V2))),]



############################################################
# Get epithelial/malignant PhenoGraph clusters (determined from topic modeling)

annot_epi <- read.table("~/cisTopic/AllSamples_TSSMixed_nFrags1k_FRIS0.4/peakset_MACS3_AllSamplesTSSMixednFrags1kFRIS0.4/Epithelial/cisTopic100_Alpha50/Perplexity_topics/cisTopic100Alpha50_EpithelialAllSamplesTSSMixnFrags1kFRIS0.4_AllSamplesMACS3_Topics50_Scaled.k15.pg.annot")

# Add distinct cluster prefix, so we can combine annot_epi with annot_allcells
annot_epi$V2 <- paste("E", annot_epi$V2, sep = "")
annot_epi$V2 <- gsub("pg_","",annot_epi$V2)

# Remove any cells that are not in matrix
annot_epi <- annot_epi[annot_epi $V1 %in% gsub("cell.","", colnames(matrix)),]

# Remove extra columns
annot_epi <- annot_epi[c(1,2)]

# Order cells by cluster number
cluster_prefix = "E"
annot_epi <- annot_epi[order(as.numeric(gsub(cluster_prefix,"", annot_epi$V2))),]



############################################################
# Replace epithelial cell cluster in annot_allcells with subclusters from annot_epi

# Remove epithelial cells from annot_allcells
annot_allcells <- annot_allcells[!(annot_allcells$V1 %in% annot_epi$V1),]

# Combine annots
annot_cluster <- rbind(annot_epi, annot_allcells)

rm(annot_allcells, annot_epi)



############################################################
# Sample IDs

# Get sample ID annot
annot_sample <- read.table("~/sample_id.annot")

# Remove any cells that are not in matrix
annot_sample <- annot_sample[annot_sample$V1 %in% gsub("cell.","", colnames(matrix)),]

# Order cells
annot_sample <- annot_sample[order(annot_sample$V2),]



############################################################
# Convert epiAneufinder results table to compatible format

# Remove "cell." prefix from colnames (barcodes)
colnames(matrix) <- gsub("cell.","", colnames(matrix))

# Turn rownames from simple counts to genomic windows
rownames(matrix) <- paste(matrix$seq, matrix$start, matrix$end, sep = "-")

# Remove "seq", "start", "end" columns
matrix <- matrix[c(-1:-3)]

# Filter matrix to contain only cells within our cluster annotation file
matrix <- matrix[,colnames(matrix) %in% annot_cluster$V1]



############################################################
# Order matrix by cluster, then by sample

# Get dataframe of cells with cluster and sample info
matrix_order <- dplyr::left_join(annot_cluster, annot_sample, by = "V1")

# Reorder dataframe by sample
matrix_order <- matrix_order[order(matrix_order$V2.y),]

# Reorder dataframe by cluster (remove cluster prefixes first)
cluster_prefix = "C"
cluster_prefix2 = "E"
matrix_order <- matrix_order[order(as.numeric( gsub(cluster_prefix, "", gsub(cluster_prefix2,"",matrix_order$V2.x)) )),]

# Reorder dataframe by epithelial clusters
matrix_order$V3 <- "Other"
matrix_order$V3[grep(cluster_prefix2,matrix_order$V2.x)] = "Epi"
matrix_order <- matrix_order[order(matrix_order$V3),]

# Reorder epiAneufinder matrix
matrix <- matrix[c(matrix_order$V1)]



############################################################
# Prepare heatmap annotations

heatmapAnnot <- HeatmapAnnotation(Cluster = matrix_order$V2.x,
                                  Sample = matrix_order$V2.y)



############################################################
# Plot Heatmap

cluster_palette <- ArchR::paletteDiscrete(values = unique(annot_cluster$V2), set = "stallion", reverse = FALSE)
sample_palette <- ArchR::paletteDiscrete(values = unique(annot_sample$V2), set = "stallion", reverse = FALSE)

heatmapAnnot_top <- HeatmapAnnotation(Cluster = matrix_order$V2.x,
                                      col = list(Cluster = cluster_palette))
heatmapAnnot_bottom <- HeatmapAnnotation(Sample = matrix_order$V2.y,
                                         col = list(Sample = sample_palette))

H_unclustered <- Heatmap(matrix,
                         top_annotation = heatmapAnnot_top,
                         bottom_annotation = heatmapAnnot_bottom,
                         cluster_rows=FALSE,
                         cluster_columns = FALSE,
                         show_row_names=F,
                         show_column_names=F,
                         
                         row_split = factor(gsub("-.*","", rownames(matrix)),
                                            levels = unique(gsub("-.*","", rownames(matrix)))),
                         row_gap = unit(2.5, "mm"),
                         row_title_gp = gpar(fontsize = 6),
                         
                         cluster_row_slices = F,
                         column_split = factor(matrix_order$V2.x,
                                               levels = unique(matrix_order$V2.x)),
                         cluster_column_slices = F,
                         column_title_gp = gpar(fontsize = 6),
                         raster_quality = 10
)
pdf(paste(getwd(),"/",OutputFileName,"_Unclustered_CellTypes_ColorPalette.pdf", sep = ""),width=18,height=12)
draw(H_unclustered)
dev.off()