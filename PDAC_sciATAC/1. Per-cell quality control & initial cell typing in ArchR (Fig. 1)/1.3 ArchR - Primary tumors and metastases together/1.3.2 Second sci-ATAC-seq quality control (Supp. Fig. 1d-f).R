library(ArchR)
library(dplyr)
library(parallel)
addArchRGenome("hg38")
addArchRThreads(32)
set.seed(1)

PlotPrefix = "AllSamples_TSS2.5_nFrags1k_FRIS0.4"



############################################################
# Prepare new ArchRProject

proj <- loadArchRProject("~/ArchR/AllSamples_TSS2.5_nFrags1k")

proj <- saveArchRProject(proj, outputDirectory = "~/ArchR/AllSamples_TSS2.5_nFrags1k_FRIS0.4")



############################################################
# Apply filters

# Read in per-cell fraction of reads overlapping Meuleman W, et al. (2020) DHS sites
qualityMetrics <- read.table("~/ArchR/AllSamples_nFrags500_QualityMetrics.txt",header = T, comment.char="")

# Remove any cell with less than 40% fraction of reads overlapping lineage DHS sites
proj <- proj[proj$cellNames %in% qualityMetrics$cell[qualityMetrics$FRIS_IndexDHS>=0.4]]



############################################################
# Save filtered ArchRProject

proj <- saveArchRProject(proj, dropCells = T)



############################################################
# iLSI and UMAP clustering

proj <- addIterativeLSI(ArchRProj = proj, useMatrix = "TileMatrix", name = "IterativeLSI",force=TRUE,iterations=3, clusterParams = list(resolution=2),varFeatures = 100000)
proj <- addClusters(input = proj, reducedDims = "IterativeLSI",force=TRUE)
proj <- addUMAP(ArchRProj = proj, reducedDims = "IterativeLSI",force=TRUE,minDist = 0.0)
proj <- addImputeWeights(proj)



############################################################
# Plot UMAP and numbers

# UMAP theme
umap_theme <- theme(
    legend.position = "none",
    
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    
    plot.title = element_blank(),
    plot.background = element_blank(),

    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    
    plot.margin = margin(t = 0,  # Top margin
                         r = 0,  # Right margin
                         b = 0,  # Bottom margin
                         l = 0),
    aspect.ratio = 1
)



p_numbers <- plotEmbedding(proj,
                           embedding = "UMAP",
                           name = "Clusters",
                           labelSize = 8/.pt,
                           labelMeans = T,
                           size = 0,
                           alpha = 0) +
    umap_theme

p_dots <- plotEmbedding(proj,
                        embedding = "UMAP",
                        name = "Clusters",
                        labelSize = 8/.pt,
                        labelMeans = F,
                        size = 0,
                        alpha = 0) +
    geom_point(stroke = 0, size = point_size) +
    umap_theme

PlotPrefix <- "AllSamplesTSS2.5nFrags1kFRIS0.4"
PlotSuffix <- "_UMAP_Clusters.pdf"
pdf(width = 1.80, height = 1.80, file = paste(PlotPrefix,PlotSuffix,sep=""))
p_dots
p_numbers
dev.off()






######################################################################
# Save ArchRProject

proj <- saveArchRProject(proj)




############################################################
# Prepare new ArchRProject

proj <- loadArchRProject("~/ArchR/AllSamples_TSS2.5_nFrags1k_FRIS0.4")

proj <- saveArchRProject(proj, outputDirectory = "~/ArchR/AllSamples_TSS2.5_nFrags1k_FRIS0.4_IndexDHS")



############################################################
# Add Index DHS peak set to ArchRProject

# Add peak set - DHS index
# Peaks obtained from Meuleman W, et al. (2020) (https://doi.org/10.1038/s41586-020-2559-3))
peaks <- read.table("~/peak_sets/DHS_Index_Meuleman_2020/DHS_Index_and_Vocabulary_hg38_WM20190703.non_overlap.bed")
peaks <- GRanges(seqnames = peaks$V1, ranges = IRanges(start = peaks$V2, end = peaks$V3))

proj <- addPeakSet(ArchRProj = proj, peakSet = peaks,force=T)
proj <- addPeakMatrix(proj,force=T)

# Annotate pan-human peak set by component
Index_DHS_components <- c(
    Tissueinvariant = '~/peak_sets/DHS_Index_Meuleman_2020/Tissueinvariant.txt',
    Placentaltrophoblast = '~/peak_sets/DHS_Index_Meuleman_2020/Placentaltrophoblast.txt',
    Neural = '~/peak_sets/DHS_Index_Meuleman_2020/Neural.txt',
    Primitiveembryonic = '~/peak_sets/DHS_Index_Meuleman_2020/Primitiveembryonic.txt',
    Digestive = '~/peak_sets/DHS_Index_Meuleman_2020/Digestive.txt',
    Lymphoid = '~/peak_sets/DHS_Index_Meuleman_2020/Lymphoid.txt',
    Musculoskeletal = '~/peak_sets/DHS_Index_Meuleman_2020/Musculoskeletal.txt',
    StromalB = '~/peak_sets/DHS_Index_Meuleman_2020/StromalB.txt',
    Myeloiderythroid = '~/peak_sets/DHS_Index_Meuleman_2020/Myeloiderythroid.txt',
    Cancerepithelial = '~/peak_sets/DHS_Index_Meuleman_2020/Cancerepithelial.txt',
    Pulmonarydevel = '~/peak_sets/DHS_Index_Meuleman_2020/Pulmonarydevel.txt',
    Renalcancer = '~/peak_sets/DHS_Index_Meuleman_2020/Renalcancer.txt',
    StromalA = '~/peak_sets/DHS_Index_Meuleman_2020/StromalA.txt',
    Organdevel.renal = '~/peak_sets/DHS_Index_Meuleman_2020/Organdevel.renal.txt',
    Vascularendothelial = '~/peak_sets/DHS_Index_Meuleman_2020/Vascularendothelial.txt',
    Cardiac = '~/peak_sets/DHS_Index_Meuleman_2020/Cardiac.txt'
)

proj <- addPeakAnnotations(ArchRProj = proj, regions = Index_DHS_components, name = "Index_DHS",force=T)

proj <- addBgdPeaks(proj,method = "ArchR",force=T)

proj <- addDeviationsMatrix(
    ArchRProj = proj, 
    peakAnnotation = "Index_DHS",
    force = TRUE
)









############################################################
# Plot per-cluster TSS enrichment

p_tss <- plotGroups(proj, groupBy = "Clusters", name = "TSSEnrichment", maxCells = 10000000, plotAs = "violin", quantCut = c(0,1)) + geom_hline(yintercept = 4)

colorPalette <- p_tss[["plot_env"]][["pal"]]
names(colorPalette) <- gsub("C","",names(colorPalette))





############################################################
# Extract IndexDHS matrix

matrix <- getMatrixFromProject(proj, useMatrix = "Index_DHSMatrix", logFile = NULL, binarize = F)
matrix <- assay(matrix, "z")
matrix <- t(matrix)



############################################################
# ggplot theme

# Theme Settings
text_settings <- element_text(colour = "black", size = 6)
line_settings <- element_line(linewidth = 0.5)

# Default Theme
default_theme <- theme(
    axis.title.x = text_settings,
    axis.title.y = text_settings,
    
    axis.ticks.x = line_settings,
    axis.ticks.y = line_settings,
    
    axis.text.x = text_settings,
    axis.text.y = text_settings,
    
    #axis.line.x = element_line(linewidth = 0.5),
    #axis.line.y = element_line(linewidth = 0.5),
    
    panel.grid = element_blank(),
    panel.background = element_blank(),
    panel.border = element_rect(linewidth = 1, fill = NA),
    
    legend.text = text_settings,
    legend.title = text_settings
)



############################################################
# Plot cancer/epithelial accessibility

gene = "Cancerepithelial"
matrix_gene <- data.frame(row.names = proj$cellNames,
                          gene_imputed = matrix[proj$cellNames,gene])

matrix_gene <- as.matrix(matrix_gene)
matrix_imputed <- imputeMatrix(
    mat = t(matrix_gene),
    imputeWeights = getImputeWeights(proj),
    threads = getArchRThreads(),
    verbose = FALSE,
    logFile = createLogFile("imputeMatrix")
)
matrix_imputed <- as.data.frame(matrix_imputed)

annot <- data.frame(cell = rownames(matrix_imputed),
                    score = matrix_imputed[[1]])
proj <- addCellColData(proj, data = annot$score, cells = annot$cell, name = "imputed_feature", force = T)

plot <- data.frame(group = proj$Clusters,
                   feature = proj$imputed_feature)
plot$group <- gsub("C","",plot$group)

plot_order <- plot %>%
    group_by(group) %>%
    summarise(median = median(feature))
plot_order <- plot_order[order(plot_order$median, decreasing = T),]


plot$group <- factor(plot$group, levels = plot_order$group)

plot$celltype <- "Stromal & Immune Clusters"
plot$celltype[plot$group %in% plot_order$group[1:13]] <- "Epithelial Clusters"


p <- ggplot(plot, aes(x = group, y = feature)) +
    geom_violin(aes(color = group, fill = group),
                alpha = 1,
                linewidth = line_width/2) +
    geom_boxplot(aes(color = group),
                 outlier.alpha = 0,
                 lwd = line_width*0.5,
                 fatten = line_width*2,
                 width = 0.5) +
    scale_fill_manual(values = colorPalette) +
    scale_color_manual(values = colorPalette) +
    scale_x_discrete(labels = scales::label_wrap(10)) +
    default_theme +
    theme(
        legend.position = "none",
        axis.text.x = element_text(angle=90, hjust = 1, vjust = 0.5),
        axis.title.y = element_text(vjust = -0.3),
        
        strip.background = element_blank(),
        strip.text = element_text(size = 6, color = "black", vjust = -0.8),

        plot.title = element_blank()
    ) +
    ylab("Imputed\nDeviation Z-Score") +
    ggtitle(paste(gene," Accessibility")) +
    geom_hline(yintercept = 0, linetype = 2, color = "gray25", linewidth = 0.50) +
    facet_grid(cols = vars(celltype), scales = "free", space = "free")




############################################################
# Plot TSS enrichment

plot <- data.frame(cell = proj$cellNames,
                   group = proj$Clusters,
                   feature = proj$TSSEnrichment)
plot$group <- gsub("C","",plot$group)
plot$group <- factor(plot$group, levels = plot_order$group)

plot$celltype <- "Stromal & Immune Clusters"
plot$celltype[plot$group %in% plot_order$group[1:13]] <- "Epithelial Clusters"

plot <- plot[plot$feature<=quantile(proj$TSSEnrichment,0.998)[[1]],]


p <- ggplot(plot, aes(x = group, y = feature)) +
    geom_violin(aes(color = group, fill = group),
                alpha = 1,
                linewidth = line_width/2) +
    geom_boxplot(aes(color = group),
                 outlier.alpha = 0,
                 lwd = line_width*0.5,
                 fatten = line_width*2,
                 width = 0.5) +
    scale_fill_manual(values = colorPalette) +
    scale_color_manual(values = colorPalette) +
    scale_x_discrete(labels = scales::label_wrap(10)) +
    default_theme +
    theme(
        legend.position = "none",
        axis.text.x = element_text(angle=90, hjust = 1, vjust = 0.5),
        axis.title.y = element_text(vjust = -0.3),
        
        strip.background = element_blank(),
        strip.text = element_text(size = 6, color = "black", vjust = -0.8),
        
        plot.title = element_blank()
    ) +
    ylab("TSS Enrichment") +
    ggtitle(paste(gene," Accessibility")) +
    geom_hline(yintercept = 4, linetype = 2, color = "gray25", linewidth = 0.50) +
    facet_grid(cols = vars(celltype), scales = "free", space = "free")