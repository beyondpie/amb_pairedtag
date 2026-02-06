#!/bin/bash

Rscript="/home/szu/mambaforge/envs/seurat/bin/Rscript"
python="/home/szu/mambaforge/bin/python"

projd="/projects/ps-renlab2/szu/projects/amb_pairedtag/"
workd="${projd}/09.epimem"
actd="${workd}/out/actRegionWithInactSignal"

#sc="326_OPC_NN"
sc="327_Oligo_NN"

# map bigwig signals
$python ${workd}/src/main/python/02.mappingBigWigOnBed.py $sc H3K4me1 H3K4me1
$python ${workd}/src/main/python/02.mappingBigWigOnBed.py $sc H3K27me3 H3K27me3



# filter regions
$Rscript ../R/04.getEpiMemRegion.R \
    $sc H3K27me3,H3K9me3 H3K27ac,H3K4me1,ATAC intersect_all_inact2act mutual

$Rscript ../R/04.getEpiMemRegion.R \
    $sc H3K27me3 H3K27ac,H3K4me1,ATAC intersect_H3K27me3_inact2act mutual

$Rscript ../R/04.getEpiMemRegion.R \
    $sc H3K27me3,H3K9me3 H3K4me1,ATAC intersect_all_inact2act mutual

$Rscript ../R/04.getEpiMemRegion.R \
    $sc H3K27me3 H3K4me1,ATAC intersect_H3K27me3_inact2act mutual

# resize region for plotting tracks
$Rscript ../R/05.extendPeak4visual.R \
         ${actd}/${sc}.H3K27me3_H3K4me1.intersect_H3K27me3_inact2act.peak \
         ${actd}/${sc}.H3K27me3_H3K4me1.intersect_H3K27me3_inact2act.bed
         
$Rscript ../R/05.extendPeak4visual.R \
         ${actd}/${sc}.H3K27me3_H3K4me1.intersect_all_inact2act.peak \
         ${actd}/${sc}.H3K27me3_H3K4me1.intersect_all_inact2act.bed

$Rscript ../R/05.extendPeak4visual.R \
         ${actd}/${sc}.H3K27me3_H3K27ac.intersect_all_inact2act.peak \
         ${actd}/${sc}.H3K27me3_H3K27ac.intersect_all_inact2act.bed

$Rscript ../R/05.extendPeak4visual.R \
         ${actd}/${sc}.H3K27me3_H3K27ac.intersect_H3K27me3_inact2act.peak \
         ${actd}/${sc}.H3K27me3_H3K27ac.intersect_H3K27me3_inact2act.bed
