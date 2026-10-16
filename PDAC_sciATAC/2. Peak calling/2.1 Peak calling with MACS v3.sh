########################################################################################################################
# Filter BAM file to contain QC-passing cells from primary tumors and metastases 
########################################################################################################################

ln -s ~/ArchR/AllSamples_TSSMixed_nFrags1k_FRIS0.4/Plots/AllSamples_TSSMixed_nFrags1k_FRIS0.4_ClusterAnnotation.annot

# Create newline-separated list of cells in ArchRProject (which all passed QC)
awk '{print $1}' AllSamples_TSSMixed_nFrags1k_FRIS0.4_ClusterAnnotation.annot > AllSamples_TSSMixed_nFrags1k_FRIS0.4_CellList.txt

# Filter BAM file to contain only cells in the list
scitools bam-filter -O AllSamplesTSSMixednFrags1kFRIS0.4 -L AllSamples_TSSMixed_nFrags1k_FRIS0.4_CellList.txt sears_tumors_hg38na.bbrd.q30.bam

# Index BAM file
samtools index -@ 20 AllSamplesTSSMixednFrags1kFRIS0.4.filt.bam



########################################################################################################################
# Run MACS (v3.0.0b3)
########################################################################################################################

source python_environments/MACS3_3.0.0b3/bin/activate

macs3 callpeak -f BAMPE -g hs --keep-dup all --min-length 100 -n AllSamplesTSSMixednFrags1kFRIS0.4_MACS3 -t AllSamplesTSSMixednFrags1kFRIS0.4.filt.bam

awk '{print $1,$2,$3}' OFS='\t' AllSamplesTSSMixednFrags1kFRIS0.4_MACS3_peaks.narrowPeak > AllSamplesTSSMixednFrags1kFRIS0.4_MACS3.bed

