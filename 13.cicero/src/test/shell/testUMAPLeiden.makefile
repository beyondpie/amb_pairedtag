Rbin := /home/szu/mambaforge/envs/seurat/bin/Rscript
workd := /projects/ps-renlab2/szu/projects/amb_pairedtag/13.cicero
leidenRscript := ${workd}/src/main/R/09.clusteringOfATACPeak.R
outd := ${workd}/out/spectral
sc := 327_Oligo_NN


.PHONY: runLeiden
runLeiden:
	${Rbin} ${leidenRscript} 10 30 0.1 weight ${sc} \
     ${outd}/${sc}/spectralEmbed.k20.csv
