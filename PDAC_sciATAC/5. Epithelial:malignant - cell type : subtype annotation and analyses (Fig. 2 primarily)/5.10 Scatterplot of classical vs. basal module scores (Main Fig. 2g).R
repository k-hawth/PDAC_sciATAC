library(ArchR)
set.seed(1)
addArchRGenome("hg38")
library(parallel)
library(dplyr)



############################################################
# Load ArchRProject

proj <- loadArchRProject("~/ArchR/Epithelial_AllSamples_TSSMixed_nFrags1k_FRIS0.4_PeaksMACS3_cisTopic50/")



########################################################################################################################
# Plot scatterplot of imputed classical vs. basal module scores
########################################################################################################################

# Load module scores
modules <- read.table("~/ArchR/Epithelial_AllSamples_TSSMixed_nFrags1k_FRIS0.4_PeaksMACS3_cisTopic50/Plots/cisTopic50_EpithelialPrimary_TSSMixed_nFrags1k_FRIS0.4_PeaksMACS3_CellIdentities_ModuleScores_Imputed.txt", header = T, sep = '\t', quote = "", comment.char = "")

# Filter for cell types of interest
modules$Cell_Type <- proj$celltype[match(modules$cellNames, proj$cellNames)]
modules <- modules[grep("Basal|Classical|Intermediate|Malignant PAV|Acinar",modules$Cell_Type),]



# Theme Settings
line_width <- 0.75
point_size <- 0.65
text_size <- 6
title_size <- 6

line_settings <- element_line(linewidth = line_width, color = "black")
text_settings <- element_text(colour = "black", size = text_size)
title_settings <- element_text(colour = "black", size = title_size, vjust = -3)



colorPalette <- c("Basal Committed" = "#bd3531",
                  "Classical-NRP" = "#f79747",
                  "Subtype-Intermediate Progenitor" = "#e1468a",
                  "Malignant PAV" = "#77a00a",
                  "Acinar/Ductal" = "#174a7e",
                  "Metaplastic NRP" = "#7c177e",
                  "Islet" = "#c48eed",
                  "Intestinal Epithelial" = "#177e53",
                  "Other PAV" = "#77a00a",
                  "Metaplastic Gastric" = "#ffd898",
                  "Pulmonary Epithelial" = "#729292",                 
                  "Non.Highlighted" = "#D3D3D3"
)



ggplot(modules, aes(x = Module.MoffittRA2015_PDAC_Classical, y = Module.MoffittRA2015_PDAC_Basal_Like, color = Cell_Type)) + geom_point(size = point_size, stroke = 0) +
    scale_color_manual(values = colorPalette) +
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
