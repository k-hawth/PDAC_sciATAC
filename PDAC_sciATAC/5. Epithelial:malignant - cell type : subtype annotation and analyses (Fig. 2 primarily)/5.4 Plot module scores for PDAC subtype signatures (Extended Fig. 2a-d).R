library(ArchR)
library(dplyr)
library(parallel)
addArchRGenome("hg38")
addArchRThreads(10)
set.seed(1)



############################################################
# Load ArchRProject

proj <- loadArchRProject("~/ArchR/Epithelial_AllSamples_TSSMixed_nFrags1k_FRIS0.4_PeaksMACS3_cisTopic50")



############################################################
# Load module scores

modules <- read.table("~/ArchR/Epithelial_AllSamples_TSSMixed_nFrags1k_FRIS0.4_PeaksMACS3_cisTopic50/Plots/cisTopic50_EpithelialPrimary_TSSMixed_nFrags1k_FRIS0.4_PeaksMACS3_CellIdentities_ModuleScores_Imputed.txt", header = T, sep = '\t', quote = "", comment.char = "")

# Add module scores to ArchRProject
for(module in grep("^Module",colnames(modules), value = T)){
    proj <- addCellColData(proj, data = modules[[module]], cells = modules$cellNames, name = module, force = T)
    rm(module)
}



############################################################
# Plot module scores

# Set file name prefix and plot output location
PlotPrefix = "cisTopic50_EpithelialPrimary_TSSMixed_nFrags1k_FRIS0.4_PeaksMACS3"
PlotDestination = paste(getOutputDirectory(proj), "/Plots/", sep = "")


for(feature_name in grep("^Module\\.",colnames(proj),value = T)){
    PlotSuffix = paste("_",feature_name,".pdf",sep="")
    
    PlotName = paste(PlotPrefix, PlotSuffix, sep = "")
    PlotPath = paste(PlotDestination, PlotName, sep = "/")
    PlotPath <- gsub("//","/",PlotPath)
    
    PlotPointSize = 0.3
    
    plot_list <- list()
    
    # Initialize plot DF, but don't show any data points (due to point size issue)
    p <- plotEmbedding(proj,
                       embedding = "cisTopic_50",
                       name = feature_name,
                       imputeWeights = NULL,
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
        width = 1.14, height = 1.14)
    print(plot_list)
    dev.off()
    rm(feature_name)
}