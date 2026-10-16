library(epiAneufinder)
set.seed(1)



############################################################
# Run epiAneufinder

epiAneufinder(
    input= "AllSamplesTSSMixednFrags1kFRIS0.4.filt.fragments.bed", 
    outdir="~/epiAneufinder/AllSamples_TSSMixed_nFrags1k_FRIS0.4/epiAneufinder", 
    blacklist="~/epiAneufinder/hg38-blacklist.v2.bed", 
    windowSize=100000, 
    genome="BSgenome.Hsapiens.UCSC.hg38",
    exclude=c('chrX','chrY','chrM'), 
    reuse.existing=TRUE,
    title_karyo="Karyogram of sample data", 
    ncores=20,
    minFrags=0)