library(ArchR)
library(dplyr)
library(parallel)
addArchRGenome("hg38")
addArchRThreads(32)
set.seed(1)



############################################################
# Prepare new ArchRProject

proj <- loadArchRProject("~/ArchR/AllSamples_TSSMixed_nFrags1k_FRIS0.4")

proj <- saveArchRProject(proj, outputDirectory = "~/ArchR/AllSamples_TSSMixed_nFrags1k_FRIS0.4_IndexDHS")



############################################################
# Add Index DHS peak set to ArchRProject

# Add peak set - DHS index
# Peaks obtained from Meuleman W, et al. (2020) (https://doi.org/10.1038/s41586-020-2559-3))
peaks <- read.table("~/peak_sets/DHS_Index_Mueleman_2020/DHS_Index_and_Vocabulary_hg38_WM20190703.non_overlap.bed")
peaks <- GRanges(seqnames = peaks$V1, ranges = IRanges(start = peaks$V2, end = peaks$V3))

proj <- addPeakSet(ArchRProj = proj, peakSet = peaks,force=T)
proj <- addPeakMatrix(proj,force=T)

# Annotate pan-human peak set by component
Index_DHS_components <- c(
    Tissueinvariant = '~/peak_sets/DHS_Index_Mueleman_2020/Tissueinvariant.txt',
    Placentaltrophoblast = '~/peak_sets/DHS_Index_Mueleman_2020/Placentaltrophoblast.txt',
    Neural = '~/peak_sets/DHS_Index_Mueleman_2020/Neural.txt',
    Primitiveembryonic = '~/peak_sets/DHS_Index_Mueleman_2020/Primitiveembryonic.txt',
    Digestive = '~/peak_sets/DHS_Index_Mueleman_2020/Digestive.txt',
    Lymphoid = '~/peak_sets/DHS_Index_Mueleman_2020/Lymphoid.txt',
    Musculoskeletal = '~/peak_sets/DHS_Index_Mueleman_2020/Musculoskeletal.txt',
    StromalB = '~/peak_sets/DHS_Index_Mueleman_2020/StromalB.txt',
    Myeloiderythroid = '~/peak_sets/DHS_Index_Mueleman_2020/Myeloiderythroid.txt',
    Cancerepithelial = '~/peak_sets/DHS_Index_Mueleman_2020/Cancerepithelial.txt',
    Pulmonarydevel = '~/peak_sets/DHS_Index_Mueleman_2020/Pulmonarydevel.txt',
    Renalcancer = '~/peak_sets/DHS_Index_Mueleman_2020/Renalcancer.txt',
    StromalA = '~/peak_sets/DHS_Index_Mueleman_2020/StromalA.txt',
    Organdevel.renal = '~/peak_sets/DHS_Index_Mueleman_2020/Organdevel.renal.txt',
    Vascularendothelial = '~/peak_sets/DHS_Index_Mueleman_2020/Vascularendothelial.txt',
    Cardiac = '~/peak_sets/DHS_Index_Mueleman_2020/Cardiac.txt'
)

proj <- addPeakAnnotations(ArchRProj = proj, regions = Index_DHS_components, name = "Index_DHS",force=T)

proj <- addBgdPeaks(proj,method = "ArchR",force=T)

proj <- addDeviationsMatrix(
    ArchRProj = proj, 
    peakAnnotation = "Index_DHS",
    force = TRUE
)



######################################################################
# Save ArchRProject

proj <- saveArchRProject(proj)




############################################################
# Plot each lineage DHS component

markerDHS <- getFeatures(proj, useMatrix = "Index_DHSMatrix")
markerDHS <- sort(grep("z:", markerDHS, value = TRUE))

for(feature_name in markerDHS){
    color_by = "Index_DHSMatrix"
    embedding_name = "UMAP"
    
    PlotPrefix = "AllSamplesTSSMixednFrags1kFRIS0.4_IndexDHS"
    PlotSuffix = paste("_",feature_name,".pdf",sep="")
    PlotDestination <- "~/ArchR/AllSamples_TSSMixed_nFrags1k_FRIS0.4_IndexDHS/Plots/"
    PlotName = paste(PlotPrefix, PlotSuffix, sep = "")
    PlotPath = paste(PlotDestination, PlotName, sep = "/")
    PlotPath <- gsub("//","/",PlotPath)
    PlotPath <- gsub("\\:","_",PlotPath)
    
    PlotPointSize = 0.3
    
    
    
    plot_list <- list()
    # Initialize plot DF, but don't show any data points (due to point size issue)
    p <- plotEmbedding(proj,
                       embedding = embedding_name,
                       colorBy = color_by,
                       name = feature_name,
                       imputeWeights = getImputeWeights(proj),
                       quantCut = c(0.05,0.95),
                       baseSize = 0,
                       alpha = 0,
                       labelMeans = F,
                       rastr = F,
                       plotAs = "points")
    
    
    
    # Set breakpoints for plot legend
    break_low <- min(p[["data"]][["color"]])
    break_high <- max(p[["data"]][["color"]])
    break_mid <- (break_high+break_low)/2
    
    break_low <- break_low + ((break_high-break_low)*0.25)
    break_high <- break_high*0.75
    
    break_low <- signif(break_low, digits = 2)
    break_high <- signif(break_high, digits = 2)
    break_mid <- signif(break_mid, digits = 2)
    
    
    
    # Add data points back manually (to circumvent point size issue)
    p <- p +
        geom_point(size = PlotPointSize, aes(stroke = 0)) +
        
        theme(
            #legend.position = "none",
            legend.text = element_text(vjust = 6.5,
                                       color = "white",
                                       #face = "bold",
                                       size = 6),
            
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank(),
            
            plot.title = element_blank(),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            
            plot.margin = margin(t = 0,  # Top margin
                                 r = 0,  # Right margin
                                 b = 0,  # Bottom margin
                                 l = 0)
        ) +
        scale_color_gradientn(colors = colorRampPalette(ArchR::ArchRPalettes$solarExtra)(1000000),
                              breaks = c(break_low,
                                         #break_mid,
                                         break_high), labels = as.character(c(break_low,break_high))) +
        
        guides(colour=guide_colourbar(barwidth=4*0.65, 
                                      barheight = 1*0.65,
                                      ticks.linewidth = 0.5))
    
    # Add to plot list
    plot_list[[length(plot_list)+1]] <- p
    
    # Add plot w/ legend to plot list
    plot_list[[length(plot_list)+1]] <- p +
        geom_point(size = PlotPointSize, aes(stroke = 0)) +
        
        theme(
            legend.position = "none",
            
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank(),
            
            plot.title = element_blank(),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            
            plot.margin = margin(t = 0,  # Top margin
                                 r = 0,  # Right margin
                                 b = 0,  # Bottom margin
                                 l = 0)
        ) +
        scale_color_gradientn(colors = colorRampPalette(ArchR::ArchRPalettes$solarExtra)(1000000))
    
    
    
    # Save plots
    pdf(file = PlotPath,
        width = 1.5, height = 1.5)
    print(plot_list)
    dev.off()
    dev.off()
}



############################################################
# Annotate per-cluster cell type

annot <- data.frame(cell = proj$cellNames,
                    sample = proj$Sample_ID,
                    cluster = proj$Clusters)

annot$celltype <- "Epithelial/Malignant"

annot$celltype[annot$cluster %in% c("C3")] <- "Plasma B"
annot$celltype[annot$cluster %in% c("C2")] <- "B-cell"
annot$celltype[annot$cluster %in% c("C6","C7","C11")] <- "T-cell"
annot$celltype[annot$cluster %in% c("C8","C9")] <- "Stromal"
annot$celltype[annot$cluster %in% c("C4","C12")] <- "Myeloid"

annot$celltype[annot$cluster %in% c("C5")] <- "Other"

proj <- addCellColData(proj, data = annot$celltype, cells = annot$cell, name = "celltype", force = T)



############################################################
# Calculate per-sample cell type proportions

proportions <- data.frame()

i = 0
for(unique_sample in unique(annot$sample)){
    
    for(unique_celltype in unique(annot$celltype)){
        i = i+1
        proportions[i,"Sample"] <- unique_sample
        proportions[i,"Celltype"] <- unique_celltype
        proportions[i,"Total_Celltype"] <- length(annot$sample[annot$sample==unique_sample & annot$celltype==unique_celltype])
        proportions[i,"Total_Cells"] <- length(annot$sample[annot$sample==unique_sample])
    }
}
rm(unique_sample,unique_celltype,i)
proportions$Percent_Celltype <- proportions$Total_Celltype/proportions$Total_Cells




############################################################
# Plot cell type proportions by sample

plot <- proportions

plot$Specimen_Type <- "Other"
plot$Specimen_Type[grep("-T",plot$Sample)] <- "Primary PDAC"
plot$Specimen_Type[plot$Sample %in% c("ST-00017310-T","ST-00019171-T")] <- "Primary PAV"
plot$Specimen_Type[plot$Sample %in% c("ST-00015839-M")] <- "Ovary Met"
plot$Specimen_Type[plot$Sample %in% c("ST-00017838-M","ST-00020190-M","ST-00020281-M","ST-00020824-M")] <- "Liver Met"
plot$Specimen_Type[plot$Sample %in% c("ST-00004898-M","ST-00017915-M")] <- "Lung Met"



# Factors
plot$Specimen_Type <- factor(plot$Specimen_Type, levels = c("Primary PDAC","Primary PAV","Lung Met","Liver Met","Ovary Met"))

plot$Sample <- factor(plot$Sample, levels = unique(plot$Sample[order(plot$Total_Cells,decreasing = T)]))

plot$Celltype <- factor(plot$Celltype, levels = rev(c("Epithelial/Malignant","Stromal","T-cell","B-cell","Plasma B","Myeloid","Other")))


pdf(width = 7.5, height = 1.75, file = "CellTypeProportions.pdf")
ggplot(plot, aes(x = Sample, y = Total_Celltype, fill = Celltype)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = as.character(ArchR::ArchRPalettes$stallion2[1:length(unique(plot$Celltype))])) +
    default_theme +
    theme(
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 60, vjust =1, hjust = 1),
        strip.text = element_text(size=6, color = "black", vjust = -1),
        strip.background = element_blank(),
        panel.spacing.x=unit(1, "lines")
    ) +
    ylab("Count") +
    facet_grid(cols = vars(plot$Specimen_Type),scales = "free_x",space = "free")
dev.off()