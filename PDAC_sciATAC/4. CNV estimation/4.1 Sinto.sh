# Input: BAM file
# Output: Compressed Fragments File

# TO RUN:
# bash Sinto.sh [bam file (with a .bam file extension) ]

bamInput=$1
trimmedName=${bamInput%.bam}
trimmedName=${trimmedName##*/}

source python_environments/Sinto/bin/activate

sinto fragments -b $bamInput -p 20 -f $trimmedName.Sinto.unsorted.bed --barcode_regex "[^:]*" --min_mapq 30 --max_distance 5000 --min_distance 10 --collapse_within

tabix="~/cellranger-7.1.0/external/anaconda/bin/tabix"
bgzip="~/cellranger-7.1.0/external/anaconda/bin/bgzip"
sort -k1,1 -k2,2n $trimmedName.Sinto.unsorted.bed > $trimmedName.fragments.bed
rm $trimmedName.Sinto.unsorted.bed
$bgzip -@ 8 $trimmedName.fragments.bed
$tabix -p bed $trimmedName.fragments.bed.gz