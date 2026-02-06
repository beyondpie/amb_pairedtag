#!/bin/bash

# work location
projd="/projects/ps-renlab2/szu/projects/amb_pairedtag"
workd="${projd}/09.epimem"
peakd="${workd}/out/actRegionWithInactSignal"
figd="${workd}/figure"

# genomic data location
genegtf="${projd}/genome/gencode.vM25.annotation.gtf.gz"
ATACbwd="${projd}/data/snATAC/subclass_bigwig_bamCoverage"
Hbwd="${projd}/data/ptDNAbam/bigwig"
RNAbwd="${projd}/data/ptRNAbam/bigwig"

#sc="326_OPC_NN"
sc="327_Oligo_NN"
outd="${figd}/${sc}"
prefix="H3K27ac"
suffix="all"
tracksfnm="${outd}/${sc}_act-${prefix}_inact-${suffix}.ini"
bed="${peakd}/${sc}.H3K27me3_${prefix}.intersect_${suffix}_inact2act.bed"
fignm="${outd}/${sc}_act-${prefix}_inact-${suffix}.pdf"

# prepare track configuration
mkdir -p ${outd}
# make_tracks_file --trackFiles \
#                  ${Hbwd}/${sc}.H3K27ac.e100.bs100.sm300.bw \
#                  ${Hbwd}/${sc}.H3K4me1.e100.bs100.sm300.bw \
#                  ${ATACbwd}/${sc}.ATAC.e100.bs100.sm300.bw \
#                  ${Hbwd}/${sc}.H3K27me3.e100.bs100.sm300.bw \
#                  ${Hbwd}/${sc}.H3K9me3.e100.bs100.sm1000.bw \
#                  ${RNAbwd}/${sc}.RPKM.bw \
#                  ${genegtf} \
#                  -o $tracksfnm

# plot figure
# To Plot One Region
# pyGenomeTracks --tracks $tracksfnm --region chr1:4487000-4496395 -o $fignm
# To Plot Multiple Regions
pyGenomeTracks --tracks $tracksfnm --BED ${bed} -o $fignm


