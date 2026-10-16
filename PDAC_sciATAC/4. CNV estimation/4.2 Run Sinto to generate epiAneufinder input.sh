# Run Sinto script to convert BAM file to scATAC-seq fragments file
bash Sinto.sh AllSamplesTSSMixednFrags1kFRIS0.4.filt.bam

# zcat fragments file for epiAneufinder input
zcat AllSamplesTSSMixednFrags1kFRIS0.4.filt.fragments.bed.gz > AllSamplesTSSMixednFrags1kFRIS0.4.filt.fragments.bed