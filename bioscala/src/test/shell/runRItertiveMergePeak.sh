#!/bin/bash

R="/home/szu/mambaforge/envs/seurat/bin/Rscript"
proj="/projects/ps-renlab2/szu/projects/amb_pairedtag"
script="${proj}/extratools/singlecell/iterativeMergePeak.R"

rscd="${proj}/src/test/resource"
inputfnm="${rscd}/peakfnm.list"
blfnm="${proj}/meta/mm10-blacklist.v2.bed"
chromfnm="${proj}/meta/mm10.chrom.sizes.lite"

$R $script --input ${inputfnm} -g mm10 -e 0 \
   --chromSize ${chromfnm} -o ${rscd}/mergedPeaks.byR.tsv \
   --ncpu 4
