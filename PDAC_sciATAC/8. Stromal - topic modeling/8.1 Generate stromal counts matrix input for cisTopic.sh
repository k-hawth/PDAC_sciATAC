########################################################################################################################
# Get list of stromal cells from primary tumors
########################################################################################################################

ln -s ~/ArchR/Primary_TSSMixed_nFrags1k_FRIS0.4/AllCells/Plots/Primary_TSSMixed_nFrags1k_FRIS0.4_ClusterAnnotation.annot

# Start R
R

df <- read.table("Primary_TSSMixed_nFrags1k_FRIS0.4_ClusterAnnotation.annot")
df <- df[df$V2 %in% c("C9","C10"),]
df <- df[,1,drop=F]
write.table(df, "FibroblastTSSMixnFrags1kFRIS0.4_CellList.txt", row.names=F,col.names=F,sep='\t',quote=F)

q()
n



########################################################################################################################
# Generate stromal counts matrix for cisTopic
########################################################################################################################

# Filter BAM file to contain only primary tumor stromal cells
scitools bam-filter -L FibroblastTSSMixnFrags1kFRIS0.4_CellList.txt sears_tumors_hg38na.bbrd.q30.bam

# Generate sparse counts matrix using called peaks identified from ArchR's iterative implementation of MACS v2
scitools atac-count -O FibroblastTSSMixnFrags1kFRIS0.4.filt sears_tumors_hg38na.bbrd.q30.filt.bam Peaks_MACS2ArchR_PrimaryTSSMixnFrags1kFRIS0.4.bed