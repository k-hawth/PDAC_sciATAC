library(ArchR)
library(dplyr)
library(parallel)
addArchRGenome("hg38")
addArchRThreads(10)
set.seed(1)



############################################################
# Load ArchRProject

proj <- loadArchRProject("~/ArchR/Primary_TSSMixed_nFrags1k_FRIS0.4/Fibroblast_MACS2ArchR_cisTopic27")


# Save new ArchRProject
proj <- saveArchRProject(proj, dropCells = T, outputDirectory = "~/ArchR/Primary_TSSMixed_nFrags1k_FRIS0.4/Fibroblast_TGFbPeaks_cisTopic27_PIN1PeakAnno")




########################################################################################################################
# Calculate per-cell accessibility of TGFb-responsive and TGFb-silenced PSC enhancers
########################################################################################################################

############################################################
# Add background set of peaks from Bowman CL and Daniel CJ, et al. (2025) (https://doi.org/10.1158/0008-5472.can-24-3437)
# This is the set of called peaks used in the paper linked above

peaks <- read.table("~/peak_sets/LangerE_TGFB_PSC/background_called_peaks_BowmanCL_DanielCJ_2025.txt")
peaks <- GRanges(seqnames = peaks$V1, ranges = IRanges(start = peaks$V2, end = peaks$V3))
proj <- addPeakSet(ArchRProj = proj, peakSet = peaks,force=T)
proj <- addPeakMatrix(proj,force=T)



############################################################
# Annotate peak set by the peaks that were found differentially accessible between shSCR untreated vs. shSCR TGFb-treated

peak_set_components <- c(list.files("~/peak_sets/PIN1_TGFb_2026/",full.names=T))
names(peak_set_components) <- gsub(".*/|.txt$","",peak_set_components)

proj <- addPeakAnnotations(ArchRProj = proj, regions = peak_set_components, name = "PSC_TGFB",force=T)
proj <- addBgdPeaks(proj,method = "ArchR",force=T)



############################################################
# Run chromVAR to calculate per-cell accessibility of the shSCR untreated peaks and the shSCR TGFb-treated peaks

proj <- addDeviationsMatrix(
    ArchRProj = proj, 
    peakAnnotation = "PSC_TGFB",
    force = TRUE
)



############################################################
# Extract matrix of per-cell chromVAR scores

matrix <- getMatrixFromProject(proj, "PSC_TGFBMatrix", useSeqnames="z", binarize = F)
matrix <- as.data.frame(t(as.matrix(assay(matrix))))

write.table(matrix, file = "~/ArchR/Primary_TSSMixed_nFrags1k_FRIS0.4/Fibroblast_TGFbPeaks_cisTopic27_PIN1PeakAnno/chromVAR_scores_PIN1_shSCR.txt", row.names = T, col.names = T, sep = '\t', quote = F)





########################################################################################################################
# Plot results
########################################################################################################################

tgfb <- read.table("~/ArchR/Primary_TSSMixed_nFrags1k_FRIS0.4/Fibroblast_TGFbPeaks_cisTopic27_PIN1PeakAnno/chromVAR_scores_PIN1_shSCR.txt", header = T, sep='\t', quote = "", comment.char = "")



############################################################
# ggplot theme

# Theme Settings
line_width <- 0.75
point_size <- 0.65
text_size <- 6
title_size <- 6

line_settings <- element_line(linewidth = line_width, color = "black")
text_settings <- element_text(colour = "black", size = text_size)
title_settings <- element_text(colour = "black", size = title_size, vjust = -3)


# Default Theme
default_theme <- theme(
    axis.title.x = element_blank(),
    axis.title.y = text_settings,
    
    axis.ticks.x = element_line(linewidth = line_width*0.5, color = "black"),
    axis.ticks.y = element_line(linewidth = line_width*0.5, color = "black"),
    
    axis.text.x = text_settings,
    axis.text.y = text_settings,
    
    axis.line.x = element_blank(),
    axis.line.y = element_blank(),
    
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "white"),
    panel.border = element_rect(linewidth = line_width, fill = NA),
    
    plot.background = element_blank(),
    plot.title = title_settings,
    
    legend.text = text_settings,
    legend.title = title_settings
)


colorPalette <- c(
    "LRRC15+ myCAF" = "#E23C3C",
    "NGFR+ myCAF" = "#b176ef",
    "iCAF" = "#FFB132",
    "PSC" = "#89DDFF",
    "Endothelial" = "#3ECE39",
    "Unknown" = "gray88"
)



############################################################
# Plot accessibility scores for peaks silenced by TGFb

plot <- data.frame(cellNames = proj$cellNames,
                   celltype = proj$celltype)

# Filter for CAF types
plot <- plot[grep("myCAF|iCAF|PSC",plot$celltype),]

# Order of plot
plot$celltype <- factor(plot$celltype, levels = c("LRRC15+ myCAF","NGFR+ myCAF","PSC","iCAF"))



plot$score <- tgfb$shSCR_UN_vs_shSCR_TGFb_Binomial_FDR0.01_log2FC1.0[match(plot$cellNames, rownames(tgfb))]

pval <- wilcox.test(
    x = plot$score[plot$celltype=="LRRC15+ myCAF"],
    y = plot$score[plot$celltype=="NGFR+ myCAF"],
)[["p.value"]]
pval <- round(pval,digits = 8)
pval <- paste("P=",pval,sep="")

pval2 <- wilcox.test(
    x = plot$score[plot$celltype=="NGFR+ myCAF"],
    y = plot$score[plot$celltype=="PSC"],
)[["p.value"]]
pval2 <- round(pval2,digits = 8)
pval2 <- paste("P=",pval2,sep="")

pval3 <- wilcox.test(
    x = plot$score[plot$celltype=="LRRC15+ myCAF"],
    y = plot$score[plot$celltype=="PSC"],
)[["p.value"]]
pval3 <- round(pval3,digits = 8)
pval3 <- paste("P=",pval3,sep="")

ggplot(plot, 
       aes(x = celltype, y = score)) +
    geom_boxplot(aes(fill = celltype)) +
    default_theme +
    theme(
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 50, hjust = 1, vjust = 1),
        title = element_text(size = 6),
        legend.position = "none"
    ) +
    geom_signif(comparisons = list(c("LRRC15+ myCAF","PSC"),
                                   c("LRRC15+ myCAF", "NGFR+ myCAF"),
                                   c("NGFR+ myCAF","PSC")),
                test = "wilcox.test",
                step_increase = 0.16,
                annotations = c(pval3,pval,pval2),
                vjust = -0.1,
                family = "Arial",
                textsize = 6/.pt,
                margin_top = 0.05) +
    ylim(0,8.5) +
    scale_x_discrete(labels = scales::label_wrap(10)) +
    scale_fill_manual(values = colorPalette) +
    ylab("Deviation Z-score")



############################################################
# Plot accessibility scores for peaks activated by TGFb

plot <- data.frame(cellNames = proj$cellNames,
                   celltype = proj$celltype)

# Filter for CAF types
plot <- plot[grep("myCAF|iCAF|PSC",plot$celltype),]

# Order of plot
plot$celltype <- factor(plot$celltype, levels = c("LRRC15+ myCAF","NGFR+ myCAF","PSC","iCAF"))



plot$score <- tgfb$shSCR_TGFb_vs_shSCR_UN_Binomial_FDR0.01_log2FC1.0[match(plot$cellNames, rownames(tgfb))]

pval <- wilcox.test(
    x = plot$score[plot$celltype=="LRRC15+ myCAF"],
    y = plot$score[plot$celltype=="NGFR+ myCAF"],
)[["p.value"]]
pval <- round(pval,digits = 8)
pval <- paste("P=",pval,sep="")

pval2 <- wilcox.test(
    x = plot$score[plot$celltype=="NGFR+ myCAF"],
    y = plot$score[plot$celltype=="PSC"],
)[["p.value"]]
pval2 <- round(pval2,digits = 8)
pval2 <- paste("P=",pval2,sep="")

pval3 <- wilcox.test(
    x = plot$score[plot$celltype=="LRRC15+ myCAF"],
    y = plot$score[plot$celltype=="PSC"],
)[["p.value"]]
pval3 <- round(pval3,digits = 8)
pval3 <- paste("P=",pval3,sep="")

ggplot(plot, 
       aes(x = celltype, y = score)) +
    geom_boxplot(aes(fill = celltype)) +
    default_theme +
    theme(
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 50, hjust = 1, vjust = 1),
        title = element_text(size = 6),
        legend.position = "none"
    ) +
    geom_signif(comparisons = list(c("LRRC15+ myCAF","PSC"),
                                   c("LRRC15+ myCAF", "NGFR+ myCAF"),
                                   c("NGFR+ myCAF","PSC")),
                test = "wilcox.test",
                step_increase = 0.11,
                annotations = c(pval3,pval,pval2),
                vjust = -0.1,
                family = "Arial",
                textsize = 6/.pt,
                margin_top = 0.05) +
    ylim(0,9.5) +
    scale_x_discrete(labels = scales::label_wrap(10)) +
    scale_fill_manual(values = colorPalette) +
    ylab("Deviation Z-score")