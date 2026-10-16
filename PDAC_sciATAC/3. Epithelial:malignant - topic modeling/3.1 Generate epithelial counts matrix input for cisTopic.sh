########################################################################################################################
# Get list of epithelial/malignant cells
########################################################################################################################

ln -s ArchR/AllSamples_TSSMixed_nFrags1k_FRIS0.4/Plots/AllSamples_TSSMixed_nFrags1k_FRIS0.4_ClusterAnnotation.annot

# Start R
R

df <- read.table("AllSamples_TSSMixed_nFrags1k_FRIS0.4_ClusterAnnotation.annot")
df <- df[df$V2 %in% c("C2","C3","C4","C5","C6","C7","C8","C9","C11","C12")==F,]
df <- df[,"V1",drop=F]
write.table(df, file = "EpithelialAllSamplesTSSMixednFrags1kFRIS0.4_CellList.txt", row.names = F, col.names = F, sep = '\t', quote = F)

q()
n



########################################################################################################################
# Generate epithelial/malignant counts matrix for cisTopic
########################################################################################################################

# Filter BAM file to contain only epithelial/malignant cells
scitools bam-filter -L EpithelialAllSamplesTSSMixednFrags1kFRIS0.4_CellList.txt sears_tumors_hg38na.bbrd.q30.bam

# Generate sparse counts matrix using peaks called by MACS v3 on all passing cells
scitools atac-count -O EpithelialAllSamplesTSSMixednFrags1kFRIS0.4_AllSamplesMACS3.filt sears_tumors_hg38na.bbrd.q30.filt.bam AllSamplesTSSMixednFrags1kFRIS0.4_MACS3.bed