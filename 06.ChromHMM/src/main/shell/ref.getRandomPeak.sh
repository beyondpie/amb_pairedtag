#!/bin/bash

#!/bin/bash

# rawPeakFile="mba.whole.union.peaks.bed"
# peakFile="mba.whole.union.peaks.nospm.bed"
# awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$4}' ${rawPeakFile} > ${peakFile}

peakFile="mba.whole.union.peaks.nospm.bed"

# * randomly shuffled regions
shuffleBed -seed 10 -i ${peakFile} -g /projects/ps-renlab/yangli/genome/mm10/mm10.chrom.sizes \
           > mba.whole.shuffle.bed

cp /projects/ps-renlab/yangli/projects/CEMBA/01.joint_dat/rs1cemba/ovlpPeaks/ovlprDHS/mm10-rDHS-Filtered.bed ./

bedtools intersect -v -wa -a mba.whole.shuffle.bed -b mba.whole.union.peaks.nospm.bed \
    | bedtools intersect -v -wa -a - -b /projects/ps-renlab/yangli/genome/mm10/mm10.blacklist.bed.gz \
    | bedtools intersect -v -wa -a - -b mm10-rDHS-Filtered.bed \
    | bedtools intersect -v -wa -a - -b mm10-cCREs.v3.bed \
    | grep -v chrUn \
    | grep -v random \
    | grep -v chrM \
    | sort -k1,1 -k2,2n | uniq > mba.whole.shuffle.removeovlp.bed
sed -i -e "s/peak/random/g" mba.whole.shuffle.removeovlp.bed

awk 'BEGIN{FS=OFS="\t"}{print $1":"$2"-"$3,$4}' \
    mba.whole.shuffle.removeovlp.bed > mba.whole.random.coor2peak.txt

