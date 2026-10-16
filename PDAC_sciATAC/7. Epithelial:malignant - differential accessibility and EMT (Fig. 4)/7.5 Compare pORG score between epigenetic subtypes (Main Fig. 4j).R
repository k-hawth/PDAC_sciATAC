library(readxl)
library(ggsignif)
library(ggplot2)



########################################################################################################################
# Compare pORG score between primary tumor subtypes
########################################################################################################################


############################################################
# Load table of primary tumors and their epigenetic subtype

subtype <- read.table("~/ArchR/Epithelial_AllSamples_TSSMixed_nFrags1k_FRIS0.4_PeaksMACS3_cisTopic50/Plots/cisTopic50_EpithelialPrimary_TSSMixed_nFrags1k_FRIS0.4_PeaksMACS3_PrimaryTumorEpigeneticSubtype.txt", header = T, sep = '\t')



############################################################
# Prepare to plot data

plot <- subtype

# Remove metastatic samples
plot <- plot[grep("-T",plot$Sample),]

# Remove PAV samples
plot <- plot[plot$Subtype!="Malignant PAV",]

# Remove samples with no subtype
plot <- plot[plot$Subtype!="No Subtype",]



############################################################
# Load pORG scores from Link JM, et al (2025) (https://doi.org/10.1038/s43018-024-00881-3)

# Supplemental_Dataset_1 obtained from https://github.com/engjen/Liver_Lung_PDAC
pORG <- readxl::read_xlsx("~/Downloads/Supplemental_Dataset_1.xlsx", sheet = 2)

# Add pORG score to plot dataframe
plot$Score <- as.numeric(pORG$pORG_Up_55_Primaries[match(plot$Sample, pORG$Public_Specimen_ID)])



############################################################
# Calculate p-values

group1 <- "Subtype-Intermediate Progenitor"
group2 <- "Basal Committed"
group3 <- "Classical-NRP"



pval <- t.test(
    x = plot$Score[plot$Subtype==group1],
    y = plot$Score[plot$Subtype==group2]
)[["p.value"]]


pval2 <- t.test(
    x = plot$Score[plot$Subtype==group1],
    y = plot$Score[plot$Subtype==group3]
)[["p.value"]]


pval3 <- t.test(
    x = plot$Score[plot$Subtype==group2],
    y = plot$Score[plot$Subtype==group3]
)[["p.value"]]



padj <- p.adjust(p = c(pval,pval2,pval3), method = "fdr")
padj <- round(padj,digits = 3)
padj <- paste("FDR=",padj,sep="")



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

colorPalette <- c("Basal Committed" = "#bd3531",
                  "Classical-NRP" = "#f79747",
                  "Subtype-Intermediate Progenitor" = "#e1468a"
)



############################################################
# Plot data

ylim_adjustment <- 1.12
step_adjustment <- 0.18

ggplot(plot, aes(x = Subtype, y = Score, fill = Subtype)) +
    scale_fill_manual(values = colorPalette) +
    geom_boxplot(outlier.alpha = 0) +
    geom_point() +
    geom_signif(comparisons = list(c(group2, group1),
                                   c(group3,group1),
                                   c(group3,group2)),
                test = "t.test",
                annotations = c(padj[1],padj[2],padj[3]),
                vjust = -0.1,
                step_increase = step_adjustment,
                #family = "Arial",
                textsize = 6/.pt,
                margin_top = 0.15) +
    default_theme +
    theme(
        legend.position = "none"
    ) +
    xlab("") + ylab("GSVA Score") +
    ylim(min(plot$Score),abs(max(plot$Score)-min(plot$Score))*ylim_adjustment)