library(ArchR)
set.seed(1)
addArchRGenome("hg38")
library(parallel)
library(dplyr)

PlotPrefix = "cisTopic50_EpithelialPrimary_TSSMixed_nFrags1k_FRIS0.4_PeaksMACS3"



############################################################
# Load ArchRProject

proj <- loadArchRProject("~/ArchR/Epithelial_AllSamples_TSSMixed_nFrags1k_FRIS0.4_PeaksMACS3_cisTopic50")



########################################################################################################################
# Plot malignant subtype proportions by sample
########################################################################################################################

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



############################################################
# Calculate proportion of epithelial/malignant cell types per sample

annot <- data.frame(cell = proj$cellNames,
                    sample = proj$Sample_ID,
                    celltype = proj$celltype)

# Count cell types by sample
annot <- annot %>%
    group_by(sample) %>%
    mutate(count_total = n()) %>%
    group_by(sample,celltype) %>%
    mutate(count_celltype = n())
annot <- annot[,colnames(annot) %in% c("cell")==F,drop=F]
annot <- unique(annot)

prop <- data.frame()
for(unique_celltype in unique(annot$celltype)){
    for(unique_sample in unique(annot$sample)){
        prop <- rbind(prop,
                      data.frame(Sample = unique_sample,
                                 Celltype = unique_celltype))
    }
}

prop$count_total <- annot$count_total[match(prop$Sample, annot$sample)]
prop$count_celltype <- annot$count_celltype[match(paste(prop$Sample,prop$Celltype),paste(annot$sample,annot$celltype))]
prop$count_celltype[is.na(prop$count_celltype)] <- 0
prop$percent_celltype <- prop$count_celltype/prop$count_total





############################################################
# Determine epigenetic subtype of each tumor sample

subtype <- prop

subtype <- subtype %>%
    group_by(Sample) %>%
    mutate(count_malignant = sum(count_celltype[Celltype %in% c("Basal Committed","Classical-NRP","Subtype-Intermediate Progenitor","Malignant PAV")])) %>%
    
    group_by(Sample,Celltype) %>%
    mutate(percent_malignant = sum(count_celltype[Celltype %in% c("Basal Committed","Classical-NRP","Subtype-Intermediate Progenitor","Malignant PAV")]/count_malignant))

subtype$percent_malignant[is.nan(subtype$percent_malignant)] <- 0

subtype <- subtype[subtype$Celltype %in% c("Basal Committed","Classical-NRP","Subtype-Intermediate Progenitor","Malignant PAV"),]
subtype <- subtype[order(subtype$percent_malignant,decreasing = T),]
subtype$Subtype <- NA

for(unique_sample in unique(subtype$Sample)){
    top_subtype <- subtype$Celltype[subtype$Sample==unique_sample][1]
    top_subtype_perc <- subtype$percent_malignant[subtype$Sample==unique_sample][1]
    
    if(length(subtype$Celltype[subtype$percent_malignant==top_subtype_perc & subtype$Sample==unique_sample])>1){
        if(top_subtype_perc==0){
            subtype$Subtype[subtype$Sample==unique_sample] <- "Zero Malignant"
        } else{
            subtype$Subtype[subtype$Sample==unique_sample] <- "Split Subtypes"
            split_subtypes <- subtype$Celltype[subtype$Sample==unique_sample & subtype$percent_malignant==top_subtype_perc]
            print(paste("Sample ",unique_sample, " has split subtypes: ", paste(split_subtypes,collapse = ", "), sep =""))
        }
    } else{
        subtype$Subtype[subtype$Sample==unique_sample] <- top_subtype
    }
    
    rm(unique_sample,top_subtype, top_subtype_perc)
}

# Split subtype samples: if split with subtype committed, assign subtype as subtype committed
subtype$Subtype[subtype$Sample=="ST-00017762-T"] <- "Basal Committed"
subtype$Subtype[subtype$Sample=="ST-00020291-T"] <- "Classical-NRP"

# Only samples that meet minimum cell count threshold are subtyped
subtype$Subtype[subtype$count_malignant<10] <- "No Subtype"




############################################################
# Plot malignant subtype proportions by sample

plot <- subtype

plot_count_order <- unique(plot[plot$Celltype=="Subtype-Intermediate Progenitor",])
plot_count_order <- plot_count_order[order(plot_count_order$percent_malignant,decreasing = T),]
plot_count_order <- plot_count_order$Sample[plot_count_order$Sample %in% c("ST-00019171-T","ST-00017310-T")==F]

plot$Sample <- factor(plot$Sample, levels = c("ST-00019171-T","ST-00017310-T",plot_count_order))
plot$Celltype <- factor(plot$Celltype, levels = c("Subtype-Intermediate Progenitor","Classical-NRP","Basal Committed","Malignant PAV"))

plot <- plot[plot$count_malignant>=10,]

plot$type <- "1 Primary PDAC"
plot$type[plot$Sample %in% c("ST-00017310-T","ST-00019171-T")] = "2 Primary PAV"
plot$type[grep("-M",plot$Sample)] = "3 Metastasis"


ggplot(plot, aes(fill=Celltype, y=percent_malignant, x=Sample)) + 
    geom_bar(position="fill", stat="identity") +
    scale_fill_manual(values = colorPalette) +
    ylab("Percent of Malignant Cells") +
    default_theme + 
    theme(
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        axis.title.x = element_blank(),
        legend.position = "none",
        strip.background = element_blank(),
        strip.text = element_blank(),
        panel.border = element_rect(linewidth = 1, fill = NA),
    ) +
    facet_grid(cols = vars(type), scales = "free",space = "free")



############################################################
# Output list of primary tumors and their epigenetic subtype classification

subtype <- unique(subtype[,c("Sample","count_malignant","Subtype"),drop=F])

write.table(subtype, file = paste(getOutputDirectory(proj),"/Plots/",PlotPrefix,"_PrimaryTumorEpigeneticSubtype",".txt",sep=""), row.names = F, col.names = T, sep = '\t', quote = F)