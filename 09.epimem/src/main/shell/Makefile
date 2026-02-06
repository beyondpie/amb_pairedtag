projd := /projects/ps-renlab2/szu/projects/amb_pairedtag
workd := ${projd}/09.epimem
epiMemScript := ${workd}/src/main/R/07.getEpiMemFromGenomeBin.R
highEpiMemScript := ${workd}/src/main/R/08.analyzeEpiMemFromGenomeBin.R
scf := ${workd}/src/main/resource/sc4epimem.txt
gf := ${workd}/src/main/resource/all_groups.txt

r := /home/szu/mambaforge/envs/seurat/bin/Rscript

runEpiMemByGenomeBin : ${epiMemScript} ${scf}
	parallel -j 8 $r $< {} :::: $(word 2,$^)

getEpiMemFromGenomeBin : ${highEpiMemScript} ${scf}
	parallel -j 8 $r $< {} :::: ${word 2,$^}


getH3K9me3GenomeBin : ${epiMemScript} ${gf}
	parallel -j 15 --no-run-if-empty $r $< {} :::: $(word 2,$^)


getAllGenomeBin : ${epiMemScript} ${gf}
	parallel -j 15 --no-run-if-empty $r $< {} H3K27ac,H3K27me3,H3K4me1 :::: $(word 2,$^)
