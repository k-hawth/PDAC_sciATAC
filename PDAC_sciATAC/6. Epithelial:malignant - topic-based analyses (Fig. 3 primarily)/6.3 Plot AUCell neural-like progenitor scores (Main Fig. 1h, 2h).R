library(ArchR)
set.seed(1)
addArchRGenome("hg38")
library(parallel)
library(dplyr)



############################################################
# Load ArchRProject

proj <- loadArchRProject("~/ArchR/Epithelial_AllSamples_TSSMixed_nFrags1k_FRIS0.4_PeaksMACS3_cisTopic50")



############################################################
# Load AUC per-cell results

auc <- read.table("~/cisTopic/AllSamples_TSSMixed_nFrags1k_FRIS0.4/peakset_MACS3_AllSamplesTSSMixednFrags1kFRIS0.4/Epithelial/cisTopic100_Alpha50/Perplexity_topics/AUC_Epithelial_Signatures/Epithelial_cisTopic50_cisTopicAUC_AllMalignantSignatures_AUCMaxRank0.03_AUCScores.txt", header = T, sep = '\t', comment.char = "")
rownames(auc) <- paste("sears_tumors_hg38#",rownames(auc),sep="")

# add to ArchRProject
proj <- addCellColData(proj, cells = rownames(auc), data = auc$HwangWL2022_MalignantLineage_NeuralLikeProgenitor, force = T, name = "auc")




feature_name = "auc"
PlotPrefix = "cisTopic50"
PlotSuffix = paste("_",feature_name,".pdf",sep="")
PlotDestination <- "~/"
PlotName = paste(PlotPrefix, PlotSuffix, sep = "")
PlotPath = paste(PlotDestination, PlotName, sep = "/")
PlotPath <- gsub("//","/",PlotPath)

PlotPointSize = 0.3

plot_list <- list()




############################################################
# Plot UMAP of neural-like progenitor scores

# Initialize plot DF, but don't show any data points (due to point size issue)
p <- plotEmbedding(proj,
                   embedding = "cisTopic_50",
                   name = feature_name,
                   imputeWeights = NULL,
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
break_high <- break_high*0.85

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
plot_list
dev.off()



############################################################
# Plot classical-basal scatterplot, but color by neural-like progenitor scores

# Load module scores
modules <- read.table("~/ArchR/Epithelial_AllSamples_TSSMixed_nFrags1k_FRIS0.4_PeaksMACS3_cisTopic50/Plots/cisTopic50_EpithelialPrimary_TSSMixed_nFrags1k_FRIS0.4_PeaksMACS3_CellIdentities_ModuleScores_Imputed.txt", header = T, sep = '\t', quote = "", comment.char = "")

# Filter for cell types of interest
modules$Cell_Type <- proj$celltype[match(modules$cellNames, proj$cellNames)]
modules <- modules[grep("Basal|Classical|Intermediate|Malignant PAV|Acinar",modules$Cell_Type),]

# Add neural-like progenitor scores
modules$HwangWL2022_MalignantLineage_NeuralLikeProgenitor <- auc$HwangWL2022_MalignantLineage_NeuralLikeProgenitor[match(modules$cellNames,rownames(auc))]



# Theme Settings
line_width <- 0.75
point_size <- 0.65
text_size <- 6
title_size <- 6

line_settings <- element_line(linewidth = line_width, color = "black")
text_settings <- element_text(colour = "black", size = text_size)
title_settings <- element_text(colour = "black", size = title_size, vjust = -3)

ggplot(modules, aes(x = Module.MoffittRA2015_PDAC_Classical, y = Module.MoffittRA2015_PDAC_Basal_Like, color = HwangWL2022_MalignantLineage_NeuralLikeProgenitor)) + geom_point(size = point_size, stroke = 0) +
    scale_color_gradientn(colors = colorRampPalette(ArchR::ArchRPalettes$solarExtra)(1000000)) +
    theme_bw() + 
    theme(aspect.ratio = 1,
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          
          axis.ticks.x = element_line(linewidth = line_width*0.25),
          axis.ticks.y = element_line(linewidth = line_width*0.25),
          
          axis.text.x = text_settings,
          axis.text.y = text_settings,
          
          panel.border = element_rect(linewidth = line_width/2, fill = NA),
          panel.grid.minor = element_line(size = 0.125*1.25), 
          panel.grid.major = element_line(size = 0.25*1.25),
          
          plot.background = element_blank(),
          plot.title = title_settings,
          
          legend.text = text_settings,
          legend.title = element_blank()
    )
